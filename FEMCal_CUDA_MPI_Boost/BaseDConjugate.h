#pragma once

extern int g_MatRowNum;


class CBaseDConjugate
{
public:
	CBaseDConjugate( size_t _RowNum = g_MatRowNum );
	virtual ~CBaseDConjugate(void);

	virtual void	HostMalloc();
	virtual void	DeviceMalloc();	
	/************************************************************************/
	/*                数据获取模块                                          */
	/************************************************************************/
	//生成随机的三对角对称阵
	void	genTridiag( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "");
	//生成随机的等式右端向量
	void	genRhs( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "..\\data\\check" );

	//从Mtx格式的文件读取三元数的稀疏矩阵信息
	HBXDef::DataAloc_t		SetStiffMatFromMTX( const char* _mat_filename = "gr_900_900_crg.mtx", const boost::filesystem::path& _dir = "..\\data\\check" );


	//从外部获取三个数组指针及其长度得到刚度矩阵的CSR格式
	HBXDef::DataAloc_t	SetStiffMat( const void *const _srcVal, const void *const _srcCol, const void *const _srcRow,
		size_t _nnA, size_t _nA, bool _bsave = false );
	//从文件中读取CSR格式的刚度矩阵
	HBXDef::DataAloc_t		SetStiffMat(const char* _NoneZeroVal = "NoneZeroVal.data",
		const char*	_ColSort = "ColSort.data",
		const char*	_ValRowSort = "NoneZeroRowSort.data",
		const boost::filesystem::path& _dir = "..\\data\\source");
	
	//从文件中读取载荷向量
	virtual HBXDef::DataAloc_t	SetLoadVec(const char* _FileName = "LoadVec.data", const boost::filesystem::path& _dir = "..\\data\\source");
	//从类外部获得载荷数组指针
	HBXDef::DataAloc_t	SetLoadVec( const void* _LoadVec, size_t _RowNum );

	/************************************************************************/
	/*                CUDA计算模块                                          */
	/************************************************************************/
public:
	//询问是否CUDA可用
	bool		isDevicePossible();
	//获取GPU设备属性，返回所有的设备数并在函数内获得最优设备号
	virtual	size_t	GetGPUDeviceQuery();
	//初始化描述器等对象
	virtual	bool		InitialDescr( cudaStream_t _stream = 0, cusparseMatrixType_t _MatrixType = CUSPARSE_MATRIX_TYPE_GENERAL);
	//设备和主机端的内存拷贝
	virtual	HBXDef::DataAloc_t		MemCpy(HBXDef::CopyType_t _temp = HBXDef::CopyType_t::HostToDevice);
	//各种特征值算法求解,返回计算时间，ms为单位.纯虚函数，每个派生类必有此函数用以求解
	//@_tol:残差
	//@_iter:迭代次数
	virtual	double	ConjugateWithGPU( const double &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER ) = 0;
	size_t		GetIters()const{ return m_iters; };
	double		GetsecUsed()const{ return msecUsed; };
	double		Getqaerr()const{return m_qaerr1;};
	//导出x的数值保存至文件
	int			ExportEigToFile(const char* _EigvalFile = "EigVal.data", boost::filesystem::path _path = "\\data\\export" );
	//释放设备端资源,不仅限于cublas，cusparse等描述器
	virtual	void		FreeGPUResource();		//各类释放的资源不一样，故用虚函数

protected:
	typedef boost::shared_ptr< HBXDef::_CSRInput<double> > CSRInput_ptr;

	bool		m_bCudaFree;		//判断是否在释放CUDA资源
	bool		m_bGenTridiag;	//判断数组是否已经malloc过，如果malloc过，hostmalloc函数里当检测指针不为0时，用delete[]；如果没有用delete
	bool		m_bGPUFlag;		//是否GPU可用
	int			m_DevID;				//所选择的最好的GPU的ID号
	bool		m_bSave;				//是否存盘

	cudaError_t			_cuError_id;		//cuda内存拷贝等错误返回值
	cublasStatus_t		_cublasError_id;		//cublas返回值
	cusparseStatus_t	_cusparseError_id;		//cusparse得返回值
	HBXDef::DataAloc_t		m_DataAloc;

	float	msecUsed;//时间统计double型变量
	double	m_qaerr1;		//残差
	size_t	m_iters;	//迭代次数

	int		m_nnzA;	//CSR格式第一个和第二个数组的长度;
	int		m_nA;	//CSR格式第三个数组的长度，等于矩阵某一维度长度+1;
	int		m_RowNum;	//矩阵的行数，M，方阵则为矩阵的某一维的长度
	int		m_ColNum;	//矩阵的列数，N

	//共轭梯度法的相关变量
	cudaStream_t		m_stream;		//当前对象所使用的流
	cublasHandle_t		cublasHandle;	//CUBLAS对象
	cusparseHandle_t	cusparseHandle;	//CUSPARSE对象
	cusparseMatDescr_t	Matdescr;		//A矩阵的描述器
	cusparseIndexBase_t	m_CSRIndexBase;	//矩阵的索引起始位，为0或者1

	//CSR格式矩阵相关参数,除中间变量的裸指针均放在结构体中
	//	HBXDef::_CSRInput<double>	m_CSRInput;	//存在问题，结构体内存连续可能导致赋值错误
	boost::filesystem::path	m_SavePath;
	//主机端buffer
	int			*h_iColSort, *h_iNonZeroRowSort;
	double		*h_NoneZeroVal;
	double		*h_x;
	double		*h_rhs;
	//显存内裸指针
	double*		d_NonZeroVal;
	int*		d_iColSort;
	int*		d_iNonZeroRowSort;
	double*		d_x;	//待求特征值
	double*		d_r;		//等式右边向量

	std::vector<double>	h_vNoneZeroVal;
	std::vector<int>	h_viColSort;
	std::vector<int>	h_viNonZeroRowSort;
	std::vector<double>	h_vRhs;

};

