#include "StdAfx.h"
#include "SSOR_PCGMatCal.h"


SSOR_PCGMatCal::SSOR_PCGMatCal(size_t _RowNum, double _w):CBaseDConjugate( _RowNum ), m_w(_w)
{
	h_Ax = nullptr;
	h_p = nullptr;
}


SSOR_PCGMatCal::~SSOR_PCGMatCal(void)
{
}

//���ļ��ж�ȡ�غ�����
HBXDef::DataAloc_t	SSOR_PCGMatCal::SetLoadVec(const char* _FileName, const boost::filesystem::path& _dir )
{
	if ( HBXDef::DataAloc_t::DATAINMEM != CBaseDConjugate::SetLoadVec( _FileName, _dir ) )
	{
		return HBXDef::DataAloc_t::INVALID;
	}
	if ( 0 != m_RowNum )
	{
		m_b.resize( m_RowNum );
		for (int _i = 0; _i < m_RowNum; _i++)
		{
			m_b[_i] = h_vRhs[_i];
		}
	}
	return HBXDef::DataAloc_t::DATAINMEM;
}

void	SSOR_PCGMatCal::HostMalloc( bool _useGPU )
{
	using namespace HBXDef;
	using namespace Eigen;

	//ʹ��EigenʱӦע�⡣Eigen��ϡ�����洢��ʽ������Ԫ������������ʹ��ʱ�������ǶԳƣ���δ��������Ԫ���������Ҫת�ã�Ĭ������Ԫ���μ�
	//http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html optionһ��
	if (!_useGPU)
	{
		//ת��Eienϡ�������ʽ
		m_TripletList.reserve(m_nnzA);
		m_TripletList.clear();
		for ( size_t _rowIndex=0; _rowIndex<h_viNonZeroRowSort.size()-1; _rowIndex++ )//��ȡÿһ���ж��ٸ��������֣���ΪiR
		{
			size_t _begin = h_viNonZeroRowSort[_rowIndex];
			size_t _end = h_viNonZeroRowSort[_rowIndex+1];

			for ( size_t _colIndex=_begin; _colIndex<_end; _colIndex++ )//��ȡiR��
			{
				m_TripletList.emplace_back( Eigen::Triplet<double>((int)_rowIndex, (int)h_viColSort[(int)_colIndex], (double)h_vNoneZeroVal[_colIndex] ) );
			}

		}
		m_SparseMat.resize(m_RowNum, m_RowNum);
		m_SparseMat.setFromTriplets( m_TripletList.begin(), m_TripletList.end() );

		//����...
		//SparseLU< SparseMatrix<double>, COLAMDOrdering<int> >  LUSolver;
//		MatrixXd FULL_A = m_SparseMat;
//		LDLT< MatrixXd > _ldlt( FULL_A );
//		_ldlt.matrixLDLT();
//		std::cout<<_ldlt.isPositive();
//	//	std::cout<<ldlt.transpositionsP();
//		MatrixXd low_Mat = _ldlt.matrixL();
//		MatrixXd upper_Mat = _ldlt.matrixU();
//		MatrixXd _diagMat = _ldlt.vectorD().asDiagonal();
//		MatrixXd _tempMat = low_Mat*low_Mat.transpose();
		//Diagonal

		//����2...
		MatrixXd FULL_A = m_SparseMat;
		LDLT< MatrixXd > _ldlt;
		_ldlt.compute(FULL_A);
		VectorXd _x = _ldlt.solve( m_b );
	}
	

}

int		SSOR_PCGMatCal::solveSSOR_PCG( )
{
	return 0;
}