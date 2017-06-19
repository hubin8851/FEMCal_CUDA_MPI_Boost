#pragma once

extern int g_MatRowNum;


class CBaseSConjugate
{
public:
	//���첢��ʼ���������CSR��ʽ�����x��b�������ڴ����
	CBaseSConjugate( size_t _RowNum = g_MatRowNum );
	virtual	~CBaseSConjugate(void);
	virtual void	HostMalloc();
	virtual void	DeviceMalloc();	
	/************************************************************************/
	/*                ���ݻ�ȡģ��                                          */
	/************************************************************************/
	//������������ԽǶԳ���
	void	genTridiag( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "");
	//��������ĵ�ʽ�Ҷ�����
	void	genRhs( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "..\\data\\check" );
	//���ⲿ��ȡ��������ָ�뼰�䳤�ȵõ��նȾ����CSR��ʽ
	HBXDef::DataAloc_t	SetStiffMat( const void *const _srcVal, const void *const _srcCol, const void *const _srcRow,
									size_t _nnA, size_t _nA, bool _bsave = false );
	//���ļ��ж�ȡCSR��ʽ�ĸնȾ���
	HBXDef::DataAloc_t		SetStiffMat(const char* _NoneZeroVal = "NoneZeroVal.data",
															const char*	_ColSort = "ColSort.data",
															const char*	_ValRowSort = "NoneZeroRowSort.data",
															const boost::filesystem::path& _dir = "..\\data\\source");
	//���ⲿ��ȡ��������ָ��õ��նȾ����CSR��ʽ
	void		SetStiffMat(const void* _srcVal, const void* _srcCol, const void* _srcRow, size_t _nnA, size_t _nA);

	//�����ⲿ����غ�����ָ��
	HBXDef::DataAloc_t		SetLoadVec(const void* _LoadVec, size_t _RowNum, bool _bsave = false );
	//���ļ��ж�ȡ�غ�����
	HBXDef::DataAloc_t		SetLoadVec(const char* _FileName = "LoadVec.data", const boost::filesystem::path& _dir = "..\\data\\source");

	/************************************************************************/
	/*                CUDA����ģ��                                          */
	/************************************************************************/
public:
	//ѯ���Ƿ�CUDA����
	bool		isDevicePossible();
	//��ȡGPU�豸���ԣ��������е��豸�����ں����ڻ�������豸��
	virtual	size_t	GetGPUDeviceQuery();
	//��ʼ���������ȶ���
	virtual	bool		InitialDescr(cusparseMatrixType_t _MatrixType = CUSPARSE_MATRIX_TYPE_GENERAL);
	//�豸�������˵��ڴ濽��
	virtual	HBXDef::DataAloc_t		MemCpy(HBXDef::CopyType_t _temp = HBXDef::CopyType_t::HostToDevice);
	//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ.���麯����ÿ����������д˺����������
	//@_tol:�в�
	//@_iter:��������,������Ϊ�����rank����Ϊ������������
	virtual	float	ConjugateWithGPU( const float &_tol = ERRORTOLERATE, const int &_iter = g_MatRowNum ) = 0;
	size_t		GetIters()const{ return m_iters; };
	float		GetsecUsed()const{ return msecUsed; };
	float		Getqaerr()const{return m_qaerr1;};
	//����x����ֵ�������ļ�
	int			ExportEigToFile(const char* _EigvalFile = "EigVal.data", boost::filesystem::path _path = "..//data//export" );
	//�ͷ��豸����Դ,��������cublas��cusparse��������
	virtual	void		FreeGPUResource();		//�����ͷŵ���Դ��һ���������麯��

protected:
	typedef boost::shared_ptr< HBXDef::_CSRInput<float> > CSRInput_ptr;

	bool		m_bCudaFree;		//�ж��Ƿ����ͷ�CUDA��Դ
	bool		m_bGenTridiag;	//�ж������Ƿ��Ѿ�malloc�������malloc����hostmalloc�����ﵱ���ָ�벻Ϊ0ʱ����delete[]�����û����delete
	bool		m_bGPUFlag;		//�Ƿ�GPU����
	int			m_DevID;				//��ѡ�����õ�GPU��ID��
	bool		m_bSave;				//�Ƿ����

	cudaError_t			_cuError_id;		//cuda�ڴ濽���ȴ��󷵻�ֵ
	cublasStatus_t		_cublasError_id;		//cublas����ֵ
	cusparseStatus_t	_cusparseError_id;		//cusparse�÷���ֵ
	HBXDef::DataAloc_t		m_DataAloc;

	float	msecUsed;//ʱ��ͳ��float�ͱ���
	float	m_qaerr1;		//�в�
	size_t	m_iters;	//��������

	int		m_nnzA;	//CSR��ʽ��һ���͵ڶ�������ĳ���;
	int		m_nA;	//CSR��ʽ����������ĳ��ȣ����ھ���ĳһά�ȳ���+1;
	int		m_RowNum;	//�����������������Ϊ�����ĳһά�ĳ���

	//�����ݶȷ�����ر���
	cublasHandle_t		cublasHandle;	//CUBLAS����
	cusparseHandle_t	cusparseHandle;	//CUSPARSE����
	cusparseMatDescr_t	Matdescr;		//A�����������
	cusparseIndexBase_t	m_CSRIndexBase;	//�����������ʼλ��Ϊ0����1

//CSR��ʽ������ز���,���м��������ָ������ڽṹ����
//	HBXDef::_CSRInput<float>	m_CSRInput;	//�������⣬�ṹ���ڴ��������ܵ��¸�ֵ����
	boost::filesystem::path	m_SavePath;
	//������buffer
	int			*h_iColSort, *h_iNonZeroRowSort;
	float		*h_NoneZeroVal;
	float		*h_x;
	float		*h_rhs;
	//�Դ�����ָ��
	float*		d_NonZeroVal;
	int*		d_iColSort;
	int*		d_iNonZeroRowSort;
	float*		d_x;	//��������ֵ
	float*		d_r;		//��ʽ�ұ�����

	std::vector<float>	h_vNoneZeroVal;
	std::vector<int>	h_viColSort;
	std::vector<int>	h_viNonZeroRowSort;
	std::vector<float>	h_vRhs;
//	HBXDef::_CSRMatrix<float>	m_CSRInput;
//�豸��thrust������ָ�뼰����.PS:�Դ�����deviece_ptr������Ҳ���õ�������malloc��...=��=��2016-8-1,hbx
//2016-10-25����Ϲ���ѳ�����ȫ������...
// 	thrust::host_vector<float>		h_vNonZeroVal;
// 	thrust::host_vector<int>			h_viColSort;
// 	thrust::host_vector<int>			h_viNonZeroRowSort;
// 	thrust::host_vector<float>		h_vfx;
// 	thrust::host_vector<float>		h_vRhs;
// 
// 	thrust::device_ptr<float>		d_pNonZeroVal;
// 	thrust::device_ptr<int>			d_piColSort;
// 	thrust::device_ptr<int>			d_piNonZeroRowSort;
// 	thrust::device_ptr<float>		d_pfx;
// 	thrust::device_ptr<float>		d_pr;


public:
	void runtest();
};

