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
	/*                ���ݻ�ȡģ��                                          */
	/************************************************************************/
	//������������ԽǶԳ���
	void	genTridiag( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "");
	//��������ĵ�ʽ�Ҷ�����
	void	genRhs( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "..\\data\\check" );

	//��Mtx��ʽ���ļ���ȡ��Ԫ����ϡ�������Ϣ
	HBXDef::DataAloc_t		SetStiffMatFromMTX( const char* _mat_filename = "gr_900_900_crg.mtx", const boost::filesystem::path& _dir = "..\\data\\check" );


	//���ⲿ��ȡ��������ָ�뼰�䳤�ȵõ��նȾ����CSR��ʽ
	HBXDef::DataAloc_t	SetStiffMat( const void *const _srcVal, const void *const _srcCol, const void *const _srcRow,
		size_t _nnA, size_t _nA, bool _bsave = false );
	//���ļ��ж�ȡCSR��ʽ�ĸնȾ���
	HBXDef::DataAloc_t		SetStiffMat(const char* _NoneZeroVal = "NoneZeroVal.data",
		const char*	_ColSort = "ColSort.data",
		const char*	_ValRowSort = "NoneZeroRowSort.data",
		const boost::filesystem::path& _dir = "..\\data\\source");
	
	//���ļ��ж�ȡ�غ�����
	virtual HBXDef::DataAloc_t	SetLoadVec(const char* _FileName = "LoadVec.data", const boost::filesystem::path& _dir = "..\\data\\source");
	//�����ⲿ����غ�����ָ��
	HBXDef::DataAloc_t	SetLoadVec( const void* _LoadVec, size_t _RowNum );

	/************************************************************************/
	/*                CUDA����ģ��                                          */
	/************************************************************************/
public:
	//ѯ���Ƿ�CUDA����
	bool		isDevicePossible();
	//��ȡGPU�豸���ԣ��������е��豸�����ں����ڻ�������豸��
	virtual	size_t	GetGPUDeviceQuery();
	//��ʼ���������ȶ���
	virtual	bool		InitialDescr( cudaStream_t _stream = 0, cusparseMatrixType_t _MatrixType = CUSPARSE_MATRIX_TYPE_GENERAL);
	//�豸�������˵��ڴ濽��
	virtual	HBXDef::DataAloc_t		MemCpy(HBXDef::CopyType_t _temp = HBXDef::CopyType_t::HostToDevice);
	//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ.���麯����ÿ����������д˺����������
	//@_tol:�в�
	//@_iter:��������
	virtual	double	ConjugateWithGPU( const double &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER ) = 0;
	size_t		GetIters()const{ return m_iters; };
	double		GetsecUsed()const{ return msecUsed; };
	double		Getqaerr()const{return m_qaerr1;};
	//����x����ֵ�������ļ�
	int			ExportEigToFile(const char* _EigvalFile = "EigVal.data", boost::filesystem::path _path = "\\data\\export" );
	//�ͷ��豸����Դ,��������cublas��cusparse��������
	virtual	void		FreeGPUResource();		//�����ͷŵ���Դ��һ���������麯��

protected:
	typedef boost::shared_ptr< HBXDef::_CSRInput<double> > CSRInput_ptr;

	bool		m_bCudaFree;		//�ж��Ƿ����ͷ�CUDA��Դ
	bool		m_bGenTridiag;	//�ж������Ƿ��Ѿ�malloc�������malloc����hostmalloc�����ﵱ���ָ�벻Ϊ0ʱ����delete[]�����û����delete
	bool		m_bGPUFlag;		//�Ƿ�GPU����
	int			m_DevID;				//��ѡ�����õ�GPU��ID��
	bool		m_bSave;				//�Ƿ����

	cudaError_t			_cuError_id;		//cuda�ڴ濽���ȴ��󷵻�ֵ
	cublasStatus_t		_cublasError_id;		//cublas����ֵ
	cusparseStatus_t	_cusparseError_id;		//cusparse�÷���ֵ
	HBXDef::DataAloc_t		m_DataAloc;

	float	msecUsed;//ʱ��ͳ��double�ͱ���
	double	m_qaerr1;		//�в�
	size_t	m_iters;	//��������

	int		m_nnzA;	//CSR��ʽ��һ���͵ڶ�������ĳ���;
	int		m_nA;	//CSR��ʽ����������ĳ��ȣ����ھ���ĳһά�ȳ���+1;
	int		m_RowNum;	//�����������M��������Ϊ�����ĳһά�ĳ���
	int		m_ColNum;	//�����������N

	//�����ݶȷ�����ر���
	cudaStream_t		m_stream;		//��ǰ������ʹ�õ���
	cublasHandle_t		cublasHandle;	//CUBLAS����
	cusparseHandle_t	cusparseHandle;	//CUSPARSE����
	cusparseMatDescr_t	Matdescr;		//A�����������
	cusparseIndexBase_t	m_CSRIndexBase;	//�����������ʼλ��Ϊ0����1

	//CSR��ʽ������ز���,���м��������ָ������ڽṹ����
	//	HBXDef::_CSRInput<double>	m_CSRInput;	//�������⣬�ṹ���ڴ��������ܵ��¸�ֵ����
	boost::filesystem::path	m_SavePath;
	//������buffer
	int			*h_iColSort, *h_iNonZeroRowSort;
	double		*h_NoneZeroVal;
	double		*h_x;
	double		*h_rhs;
	//�Դ�����ָ��
	double*		d_NonZeroVal;
	int*		d_iColSort;
	int*		d_iNonZeroRowSort;
	double*		d_x;	//��������ֵ
	double*		d_r;		//��ʽ�ұ�����

	std::vector<double>	h_vNoneZeroVal;
	std::vector<int>	h_viColSort;
	std::vector<int>	h_viNonZeroRowSort;
	std::vector<double>	h_vRhs;

};

