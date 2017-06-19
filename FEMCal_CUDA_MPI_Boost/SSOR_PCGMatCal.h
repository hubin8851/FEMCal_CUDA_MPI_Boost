#pragma once

#include "HBXFEMDefStruct.h"
#include "basedconjugate.h"

class SSOR_PCGMatCal :
	public CBaseDConjugate
{
public:
	SSOR_PCGMatCal( size_t _RowNum = g_MatRowNum, double _w = 1.0 );
	~SSOR_PCGMatCal(void);

	void	genLaplace( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath = "" );
	void	HostMalloc( bool _useGPU = false );//CPU��Eigen�汾

	//���ļ��ж�ȡ�غ�����
	HBXDef::DataAloc_t	SetLoadVec(const char* _FileName = "LoadVec.data", const boost::filesystem::path& _dir = "..\\data\\source");

	int		solveSSOR_PCG( );

	//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ,��ʱ����
	//@_tol:�в�
	//@_iter:��������
	double	ConjugateWithGPU( const double &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER ){ return -1.0f; };

private:
	//����������ʱ����
	double m_w;	//�ɳ�����
	double *h_Ax, *h_p;	


#pragma region CPU������ز���������Eigen

	std::vector< Eigen::Triplet<double> > m_TripletList;
	//ʹ��EigenʱӦע�⡣Eigen��ϡ�����洢��ʽ������Ԫ������������ʹ��ʱ�������ǶԳƣ���δ��������Ԫ���������Ҫת�ã�Ĭ������Ԫ���μ�
	//http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html optionһ��
	Eigen::SparseMatrix< double >		m_SparseMat;	//���ڼ���ľ���A

	Eigen::VectorXd		m_x,m_b;

#pragma endregion 


};

