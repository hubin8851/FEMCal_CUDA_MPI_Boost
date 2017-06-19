#pragma once
//在Eigen中，默认是列索引
class EigenCal
{
public:
	EigenCal(size_t _Row);
	~EigenCal(void);
public:
	void	HostMalloc();	//在初始化了m_rownum之后对向量的内存分配
	//生成三对角矩阵
	void	genTridiag( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "");
	//生成随机的等式右端向量
	void	genRhs( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "..\\data\\check" );

	//从外部获取三个数组指针及其长度得到刚度矩阵的CSR格式
	HBXDef::DataAloc_t	SetStiffMat( const void *const _srcVal, const void *const _srcCol, const void *const _srcRow, size_t _nnA, size_t _nA, bool _bsave = false );
	//从文件中读取CSR格式的刚度矩阵
	HBXDef::DataAloc_t	SetStiffMat( const char* _NoneZeroVal = "NoneZeroVal.data",
										const char*	_ColSort = "ColSort.data",
										const char*	_ValRowSort = "NoneZeroRowSort.data",
										const boost::filesystem::path& _dir = "..\\data\\source");
	//从类外部获得载荷数组指针
	HBXDef::DataAloc_t	SetLoadVec( const void* _LoadVec, size_t _RowNum, bool _bsave = false );
	//从文件中读取载荷向量
	HBXDef::DataAloc_t	SetLoadVec(const char* _FileName = "LoadVec.data", const boost::filesystem::path& _dir = "..\\data\\source");
	//A*x = b 已知A和b求x
	void	EigenCalc( HBXDef::SolveMethod_t _Mtd = HBXDef::SolveMethod_t::CG_Mtd );
	void	EigenCalc_SSOR();

	//导出x的数值至文件
	int		ExportEigToFile(const char* _EigvalFile = "Eigen_EigVal.data", boost::filesystem::path _path = "..\\data\\Eigen");
protected:
	typedef boost::shared_ptr< HBXDef::_CSRInput<double> > CSRInput_ptr;

	bool	m_bGenTridiag;	//判断数组是否已经malloc过，如果malloc过，hostmalloc函数里当检测指针不为0时，用delete[]；如果没有用delete
	bool	m_bSave;				//是否存盘
	boost::filesystem::path	m_SavePath;

	int		m_nnzA;	//CSR格式第一个和第二个数组的长度;
	int		m_nA;	//CSR格式第三个数组的长度，等于矩阵某一维度长度+1;
	int		m_RowNum;	//矩阵的行数，方阵则为矩阵的某一维的长度

	cusparseIndexBase_t	m_CSRIndexBase;	//矩阵的索引起始位，为0或者1


	std::vector<double>	m_vNoneZeroVal;
	std::vector<int>	m_viColSort;
	std::vector<int>	m_viNonZeroRowSort;
	std::vector< Eigen::Triplet<double> > m_TripletList;

	std::vector<double>	m_vRhs;
	Eigen::VectorXd		m_x,m_b;

	Eigen::SparseMatrix<double>		m_SparseMat;	//用于计算的矩阵A
	


};

