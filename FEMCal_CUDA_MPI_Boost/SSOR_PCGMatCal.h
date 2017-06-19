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
	void	HostMalloc( bool _useGPU = false );//CPU端Eigen版本

	//从文件中读取载荷向量
	HBXDef::DataAloc_t	SetLoadVec(const char* _FileName = "LoadVec.data", const boost::filesystem::path& _dir = "..\\data\\source");

	int		solveSSOR_PCG( );

	//各种特征值算法求解,返回计算时间，ms为单位,暂时不用
	//@_tol:残差
	//@_iter:迭代次数
	double	ConjugateWithGPU( const double &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER ){ return -1.0f; };

private:
	//求解的两个临时数组
	double m_w;	//松弛因子
	double *h_Ax, *h_p;	


#pragma region CPU计算相关参数，利用Eigen

	std::vector< Eigen::Triplet<double> > m_TripletList;
	//使用Eigen时应注意。Eigen的稀疏矩阵存储格式是列主元索引，所以在使用时矩阵若非对称，在未设置行主元的情况下需要转置，默认列主元，参见
	//http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html option一项
	Eigen::SparseMatrix< double >		m_SparseMat;	//用于计算的矩阵A

	Eigen::VectorXd		m_x,m_b;

#pragma endregion 


};

