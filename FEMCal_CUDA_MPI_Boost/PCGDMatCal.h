#pragma once
#include "basedconjugate.h"
class PCGDMatCal :
	public CBaseDConjugate
{
public:
	PCGDMatCal( size_t _RowNum = g_MatRowNum );
	~PCGDMatCal(void);
	void	FreeGPUResource();

	void	genLaplace( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath = "" );
	void	DeviceMalloc();//重载devicemalloc，因为多了临时数组

	//初始化描述器等对象
	bool	InitialDescr( cudaStream_t _stream = 0, cusparseMatrixType_t _MatrixType = CUSPARSE_MATRIX_TYPE_GENERAL);
	//各种特征值算法求解,返回计算时间，ms为单位
	//@_tol:残差
	//@_iter:迭代次数
	double	ConjugateWithGPU( const double &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER );

private:
	//预处理共轭梯度法的相关变量
	cusparseSolveAnalysisInfo_t	infoA;
	cusparseSolveAnalysisInfo_t	info_u;
	cusparseMatDescr_t descrL;
	cusparseMatDescr_t descrU;

	double *h_ILUvals;
	double *d_ILUvals;
	double *d_p, *d_omega, *d_y;	//PCG法求解的设备端临时数组
	double *d_zm1, *d_zm2, *d_rm2;
};

