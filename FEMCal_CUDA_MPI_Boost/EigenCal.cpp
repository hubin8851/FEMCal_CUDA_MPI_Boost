#include "StdAfx.h"
#include "HbxGloFunc.h"
#include "EigenCal.h"


EigenCal::EigenCal(size_t _Row):m_RowNum(std::move(_Row)),m_nA(std::move(_Row+1))
{

}


EigenCal::~EigenCal(void)
{
}

//�ڳ�ʼ����m_rownum֮����������ڴ����
void	EigenCal::HostMalloc()
{

}

void EigenCal::genTridiag( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath )
{
	using namespace std;
	srand((unsigned)time(NULL));
	cout<<"�������ԽǾ���ԭʼ�������ڴ���Ĩȥ..."<<endl;
	//��Ϊǿ����Ϊ���ԽǾ��󣬹�m_nnA,m_nA,rownum����Ҫ���¸�ֵ
	m_RowNum = (int)_genRowNum;
	m_nnzA = (int)(m_RowNum-2)*3+4; 
	m_nA = (int)m_RowNum+1;

#pragma region ���·����ڴ棬�����������
	//�ڴ����
	m_viColSort.resize( m_nnzA );
	m_vNoneZeroVal.resize( m_nnzA );
	m_viNonZeroRowSort.resize( m_nA );
	m_x.resize( m_RowNum );
	//��ֵ����
	m_viNonZeroRowSort[0] = 0, m_viColSort[0] = 0, m_viColSort[1] = 1;
	m_vNoneZeroVal[0] = (float)rand()/RAND_MAX + 10.0f;
	m_vNoneZeroVal[1] = (float)rand()/RAND_MAX;
	int start = 0;
	for (size_t i = 1; i < m_RowNum; i++)
	{
		if (i > 1)
		{
			m_viNonZeroRowSort[i] = m_viNonZeroRowSort[i-1]+3;
		}
		else
		{
			m_viNonZeroRowSort[1] = 2;
		}

		start = ((int)i-1)*3 + 2;//startΪ�����з������е�i�еĵ�һ����
		m_viColSort[start] = (int)i - 1;
		m_viColSort[start+1] = (int)i;

		if (i < m_RowNum-1)
		{
			m_viColSort[start+2] = (int)i + 1;
		}

		m_vNoneZeroVal[start] = m_vNoneZeroVal[start-1];//��һ�еĶԳ�ֵ
		m_vNoneZeroVal[start+1] = (float)rand()/RAND_MAX + 10.0f;	//�Խ����ϵ�ֵ
		if (i < m_RowNum-1)
		{
			m_vNoneZeroVal[start+2] = (float)rand()/RAND_MAX;
		}
	}
	m_viNonZeroRowSort[m_RowNum] = (m_RowNum-2)*3+4;	//���һ����
#pragma endregion ���·����ڴ棬�����������

//CSR��ʽ�����Ƿ�洢����ش��룬20161214,hbx
	if ( _bsave )
	{
		using namespace boost::gregorian;
		using namespace boost::posix_time;

		//		ptime	_pathtime = second_clock::local_time();
		if ("" == _savepath)
		{
			m_SavePath = "..\\data\\export\\";
		}
		else m_SavePath = _savepath;
		//		m_SavePath.append( to_iso_string(_pathtime) );
		if (exists(m_SavePath))
		{
			if ( boost::filesystem::is_empty(m_SavePath) )
			{
				remove(m_SavePath);
			}
		}
		else
		{
			remove_all(m_SavePath);
		}
		assert(!exists(m_SavePath));
		create_directories(m_SavePath);
		HBXDef::_ExportCSRdataToFile( m_vNoneZeroVal.data(), m_viColSort.data(), m_viNonZeroRowSort.data(),
			m_nnzA, m_RowNum, "NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data",
			m_SavePath );
	}
	m_bGenTridiag = true;
}

//��������ĵ�ʽ�Ҷ�����
void	EigenCal::genRhs( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath )
{
	m_b = Eigen::VectorXd::Random( m_RowNum );
	if (_bsave)
	{
		using namespace boost::gregorian;
		using namespace boost::posix_time;

		if (exists(m_SavePath))
		{
			if ( boost::filesystem::is_empty(m_SavePath) )
			{
				remove(m_SavePath);
			}
		}
		else
		{
			remove_all(m_SavePath);
		}
		assert(!exists(m_SavePath));
		create_directories(m_SavePath);
		HBXDef::OutputPointData( m_b.data(), m_RowNum, "LoadVec.data", m_SavePath );
	}
}

//���ⲿ��ȡ��������ָ��õ��նȾ����CSR��ʽ
//@_srcVal��Դ��ֵ����
//@_srcCol:����������
//@_srcRow:����������������
//@_nnA:����ֵ��Ŀ
//@_nA:rowNum+1
HBXDef::DataAloc_t	EigenCal::SetStiffMat( const void *const _srcVal, const void *const _srcCol, const void *const _srcRow, size_t _nnA, size_t _nA, bool _bsave)
{
	m_vNoneZeroVal.clear();
	m_viColSort.clear();
	m_viNonZeroRowSort.clear();
	using namespace HBXDef;
	const double *_tmpVal = static_cast<const double*>(_srcVal);
	const size_t *_tmpCol = static_cast<const size_t*>(_srcCol);
	const size_t *_tmpRow = static_cast<const size_t*>(_srcRow);
	double *h_NoneZeroVal = const_cast<double*>(_tmpVal);
	size_t *h_iColSort = const_cast<size_t*>(_tmpCol);
	size_t *h_iNonZeroRowSort = const_cast<size_t*>(_tmpRow);
	m_nnzA = (int)_nnA;
	m_nA = (int)_nA;
	m_vNoneZeroVal.resize(m_nnzA);
	m_viColSort.resize(m_nnzA);
	m_viNonZeroRowSort.resize(m_nA);
	//�ڴ���ÿ��Ԫ�ظ�ֵ����Ϊvector���ڴ���ܲ�����
	for (int i = 0; i < m_nnzA; i++)
	{
		m_vNoneZeroVal[i] = h_NoneZeroVal[i];
		m_viColSort[i] = (int)h_iColSort[i];
	}
	for (int i = 0; i < m_nA; i++)
	{
		m_viNonZeroRowSort[i] = (int)h_iNonZeroRowSort[i];
	}
	//�����浵�����ڶ������ݵĳ���У��
	if ( _bsave )
	{
		using namespace boost::gregorian;
		using namespace boost::posix_time;

		m_SavePath = "..\\data\\Eigen\\";
		//		m_SavePath.append( to_iso_string(_pathtime) );
		if (exists(m_SavePath))
		{
			if ( boost::filesystem::is_empty(m_SavePath) )
			{
				remove(m_SavePath);
			}
		}
		else
		{
			remove_all(m_SavePath);
		}
		assert(!exists(m_SavePath));
		create_directories(m_SavePath);
		HBXDef::_ExportCSRdataToFile( h_NoneZeroVal, h_iColSort, h_iNonZeroRowSort,
			m_nnzA, m_RowNum, "NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data",
			m_SavePath );
	}
	return DataAloc_t::DATAINMEM;
}

//���ⲿ��ȡ��������ָ��õ��նȾ����CSR��ʽ
//@_srcVal��Դ��ֵ
//@_srcCol:����������
//@_srcRow:����������������
//@_dir:������·��
HBXDef::DataAloc_t	EigenCal::SetStiffMat(const char* _NoneZeroVal,const char* _ColSort,const char* _ValRowSort, const boost::filesystem::path& _dir)
{
	if (0 == HBXDef::_ImportCSRdataFromFile<double, int>(  m_vNoneZeroVal, m_viColSort, m_viNonZeroRowSort, _NoneZeroVal, _ColSort, _ValRowSort, _dir) )
	{
		std::cout<<"δ����ȷ��ȡCSR�ļ�����..."<<std::endl;
		return HBXDef::DataAloc_t::INVALID;
	}
	if ( m_nA == m_viNonZeroRowSort.size() )
	{
		std::cout<<"��ȷ��ȡCSR����,ά����ȷ..."<<std::endl;
		if (0 == m_viColSort[0] && 0 == m_viNonZeroRowSort[0])
		{
			m_CSRIndexBase = CUSPARSE_INDEX_BASE_ZERO;
			std::cout<<"����Ϊ������ʼ..."<<std::endl;
		}
		else if (1 == m_viColSort[0] && 1 == m_viNonZeroRowSort[0])
		{
			m_CSRIndexBase = CUSPARSE_INDEX_BASE_ONE;
			std::cout<<"��һΪ������ʼ..."<<std::endl;
		}
		else return HBXDef::DataAloc_t::INVALID;
		
		return	HBXDef::DataAloc_t::DATAINMEM;
	}
	else
	{
		std::cout<<"δ����ȷ��ȡCSR�ļ����ݣ�����ά�Ȳ���..."<<std::endl;
		return HBXDef::DataAloc_t::INVALID;
	}
}

//�����ⲿ����غ�����ָ��
HBXDef::DataAloc_t	EigenCal::SetLoadVec( const void* _LoadVec, size_t _RowNum, bool _bsave )
{
	using namespace HBXDef;
	const double* _tmpload = static_cast<const double*>(_LoadVec);
	m_vRhs.resize(_RowNum);
	m_b.resize(_RowNum);
	double *h_rhs = const_cast<double*>(_tmpload);
	m_nA = (int)_RowNum + 1;
	//�ڴ���ÿ��Ԫ�ظ�ֵ����Ϊvector���ڴ���ܲ�����
	for (int i = 0; i < _RowNum; i++)
	{
		m_vRhs[i] = h_rhs[i];
		m_b[i] = h_rhs[i];
	}
	if (_bsave)
	{
		using namespace boost::gregorian;
		using namespace boost::posix_time;

		if (exists(m_SavePath))
		{
			if ( boost::filesystem::is_empty(m_SavePath) )
			{
				remove(m_SavePath);
			}
		}
		else
		{
			remove_all(m_SavePath);
		}
		assert(!exists(m_SavePath));
		create_directories(m_SavePath);
		HBXDef::OutputPointData( h_rhs, m_RowNum, "LoadVec.data", m_SavePath );
	}
	return DataAloc_t::DATAINMEM;
}

//���ļ��ж�ȡ�غ�����
HBXDef::DataAloc_t	EigenCal::SetLoadVec(const char* _FileName, const boost::filesystem::path& _dir)
{
	size_t _RowNum = HBXDef::ReadVecByRowFromFile( m_vRhs,	_FileName, _dir);
	if (0 == _RowNum)
	{
		std::cout<<"δ��ȡ�غ�����..."<<std::endl;
		return HBXDef::DataAloc_t::INVALID;
	}
	else if ( m_RowNum == _RowNum)
	{
		std::cout<<"�غ�����ά����ȷ..."<<std::endl;
		m_b.resize( m_RowNum );
		for (int _i = 0; _i < m_RowNum; _i++)
		{
			m_b[_i] = m_vRhs[_i];
		}
	//	std::cout<<m_b<<"\n";	//��֤��
		return HBXDef::DataAloc_t::DATAINMEM;
	}
	std::cout<<"�غ�����ά�ȴ���..."<<std::endl;
	return HBXDef::DataAloc_t::INVALID;
}

//A*x = b ��֪A��b��x
void	EigenCal::EigenCalc( HBXDef::SolveMethod_t _Mtd )
{
	using namespace Eigen;
	//ǰ������
	//��m_TripletList��ÿ��Ԫ�ظ�ֵ
	//Ԥ����Ŀ
	m_TripletList.reserve(m_nnzA);
	m_TripletList.clear();
	for ( size_t _rowIndex=0; _rowIndex<m_viNonZeroRowSort.size()-1; _rowIndex++ )//��ȡÿһ���ж��ٸ��������֣���ΪiR
	{
		size_t _begin = m_viNonZeroRowSort[_rowIndex];
		size_t _end = m_viNonZeroRowSort[_rowIndex+1];

		for ( size_t _colIndex=_begin; _colIndex<_end; _colIndex++ )//��ȡiR��
		{
			m_TripletList.emplace_back( Eigen::Triplet<double>((int)_rowIndex, (int)m_viColSort[(int)_colIndex], (double)m_vNoneZeroVal[_colIndex] ) );
		}

	}
	m_SparseMat.resize(m_RowNum, m_RowNum);
	m_SparseMat.setFromTriplets( m_TripletList.begin(), m_TripletList.end() );
	if ( m_RowNum < 20 )
	{
		std::cout<<m_SparseMat;
	}
	//���㲿��
	if (_Mtd == HBXDef::SolveMethod_t::QR_Mtd)
	{
		SparseQR< SparseMatrix<double>, AMDOrdering<int> > qr;
		boost::timer::cpu_timer timer;	//boost�߾��ȼ�ʱ��
		timer.start();
		//����ֽ�
		qr.compute( m_SparseMat );
		//��Ax = b
		m_x = qr.solve( m_b );
		std::cout << "use move" << timer.format(6) << std::endl;
	}
	else if ( _Mtd == HBXDef::SolveMethod_t::BICGSTAB_Mtd )	//�ο����ף�http://eigen.tuxfamily.org/dox-devel/classEigen_1_1BiCGSTAB.html
	{
		using namespace Eigen;
		BiCGSTAB< SparseMatrix<double> > solver;
		boost::timer::cpu_timer timer;	//boost�߾��ȼ�ʱ��
		timer.start();
		solver.compute( m_SparseMat );
		m_x = solver.solve(m_b);
		std::cout << "#iterations:     " << solver.iterations() << std::endl;
		std::cout << "estimated error: " << solver.error()      << std::endl;
		std::cout << "use move" << timer.format(6) << std::endl;
		Eigen::VectorXd	 _tmpb = m_SparseMat*m_x;
	}
	else if ( _Mtd == HBXDef::SolveMethod_t::CG_Mtd )	//�ο����ף�http://eigen.tuxfamily.org/dox-devel/classEigen_1_1ConjugateGradient.html
	{
		using namespace Eigen;
		ConjugateGradient< SparseMatrix<double>, Lower|Upper > cg;
		boost::timer::cpu_timer timer;	//boost�߾��ȼ�ʱ��
		timer.start();
		cg.compute( m_SparseMat );
		m_x = cg.solve(m_b);
		std::cout << "#iterations:     " << cg.iterations() << std::endl;
		std::cout << "estimated error: " << cg.error()      << std::endl;
		std::cout << "use move" << timer.format(6) << std::endl;
		Eigen::VectorXd	 _tmpb = m_SparseMat*m_x;
	}
	else if (_Mtd == HBXDef::SolveMethod_t::ILUT_Mtd)	//�ο����ף�http://eigen.tuxfamily.org/dox-devel/classEigen_1_1IncompleteLUT.html
	{

	}
	else if (_Mtd == HBXDef::SolveMethod_t::LeastSquareCG_Mtd)		//�ο����ף�http://eigen.tuxfamily.org/dox-devel/classEigen_1_1LeastSquaresConjugateGradient.html
	{
		Eigen::LeastSquaresConjugateGradient< SparseMatrix<double> > lscg;
		boost::timer::cpu_timer timer;	//boost�߾��ȼ�ʱ��
		timer.start();
		lscg.compute( m_SparseMat );
		m_x = lscg.solve( m_b );
		std::cout << "#iterations:     " << lscg.iterations() << std::endl;
		std::cout << "estimated error: " << lscg.error()      << std::endl;
		std::cout << "use move" << timer.format(6) << std::endl;
		Eigen::VectorXd	 _tmpb = m_SparseMat*m_x;
	}


}

void	EigenCal::EigenCalc_SSOR()
{
	using namespace Eigen;
	//ǰ������
	//��m_TripletList��ÿ��Ԫ�ظ�ֵ
	for ( size_t _rowIndex=0; _rowIndex<m_viNonZeroRowSort.size()-1; _rowIndex++ )//��ȡÿһ���ж��ٸ��������֣���ΪiR
	{
		size_t _begin = m_viNonZeroRowSort[_rowIndex];
		size_t _end = m_viNonZeroRowSort[_rowIndex+1];
		for ( size_t _colIndex=_begin; _colIndex<_end; _colIndex++ )//��ȡiR��
		{
			m_TripletList.emplace_back( Eigen::Triplet<double>((int)_rowIndex, (int)m_viColSort[_colIndex], (double)m_vNoneZeroVal[_colIndex] ) );
		}
	}
	m_SparseMat.resize(m_RowNum, m_RowNum);
	m_SparseMat.setFromTriplets( m_TripletList.begin(), m_TripletList.end() );

	//���㲿��
}

//����x����ֵ���ļ�
int		EigenCal::ExportEigToFile( const char* _EigvalFile, boost::filesystem::path _path )
{
	HBXDef::CheckUserDefErrors (HBXDef::OutputPointData( m_x.data(), m_RowNum, _EigvalFile, _path ) );
	std::cout<<"����ֵ�������"
		<<_EigvalFile<<"�ļ���"<<std::endl;
	return true;
}