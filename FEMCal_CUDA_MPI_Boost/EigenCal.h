#pragma once
//��Eigen�У�Ĭ����������
class EigenCal
{
public:
	EigenCal(size_t _Row);
	~EigenCal(void);
public:
	void	HostMalloc();	//�ڳ�ʼ����m_rownum֮����������ڴ����
	//�������ԽǾ���
	void	genTridiag( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "");
	//��������ĵ�ʽ�Ҷ�����
	void	genRhs( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "..\\data\\check" );

	//���ⲿ��ȡ��������ָ�뼰�䳤�ȵõ��նȾ����CSR��ʽ
	HBXDef::DataAloc_t	SetStiffMat( const void *const _srcVal, const void *const _srcCol, const void *const _srcRow, size_t _nnA, size_t _nA, bool _bsave = false );
	//���ļ��ж�ȡCSR��ʽ�ĸնȾ���
	HBXDef::DataAloc_t	SetStiffMat( const char* _NoneZeroVal = "NoneZeroVal.data",
										const char*	_ColSort = "ColSort.data",
										const char*	_ValRowSort = "NoneZeroRowSort.data",
										const boost::filesystem::path& _dir = "..\\data\\source");
	//�����ⲿ����غ�����ָ��
	HBXDef::DataAloc_t	SetLoadVec( const void* _LoadVec, size_t _RowNum, bool _bsave = false );
	//���ļ��ж�ȡ�غ�����
	HBXDef::DataAloc_t	SetLoadVec(const char* _FileName = "LoadVec.data", const boost::filesystem::path& _dir = "..\\data\\source");
	//A*x = b ��֪A��b��x
	void	EigenCalc( HBXDef::SolveMethod_t _Mtd = HBXDef::SolveMethod_t::CG_Mtd );
	void	EigenCalc_SSOR();

	//����x����ֵ���ļ�
	int		ExportEigToFile(const char* _EigvalFile = "Eigen_EigVal.data", boost::filesystem::path _path = "..\\data\\Eigen");
protected:
	typedef boost::shared_ptr< HBXDef::_CSRInput<double> > CSRInput_ptr;

	bool	m_bGenTridiag;	//�ж������Ƿ��Ѿ�malloc�������malloc����hostmalloc�����ﵱ���ָ�벻Ϊ0ʱ����delete[]�����û����delete
	bool	m_bSave;				//�Ƿ����
	boost::filesystem::path	m_SavePath;

	int		m_nnzA;	//CSR��ʽ��һ���͵ڶ�������ĳ���;
	int		m_nA;	//CSR��ʽ����������ĳ��ȣ����ھ���ĳһά�ȳ���+1;
	int		m_RowNum;	//�����������������Ϊ�����ĳһά�ĳ���

	cusparseIndexBase_t	m_CSRIndexBase;	//�����������ʼλ��Ϊ0����1


	std::vector<double>	m_vNoneZeroVal;
	std::vector<int>	m_viColSort;
	std::vector<int>	m_viNonZeroRowSort;
	std::vector< Eigen::Triplet<double> > m_TripletList;

	std::vector<double>	m_vRhs;
	Eigen::VectorXd		m_x,m_b;

	Eigen::SparseMatrix<double>		m_SparseMat;	//���ڼ���ľ���A
	


};

