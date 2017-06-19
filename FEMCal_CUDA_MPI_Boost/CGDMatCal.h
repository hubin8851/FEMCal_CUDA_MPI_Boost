#pragma once
#include "basedconjugate.h"
class CGDMatCal :
	public CBaseDConjugate
{
public:
	CGDMatCal( size_t _RowNum = g_MatRowNum );
	~CGDMatCal(void);
	void	FreeGPUResource();

	void	DeviceMalloc();//重载devicemalloc，因为多了两个临时数组

	//各种特征值算法求解,返回计算时间，ms为单位.纯虚函数，每个派生类必有此函数用以求解
	//@_tol:残差
	//@_iter:迭代次数
	double	ConjugateWithGPU( const double &_tol = ERRORTOLERATE, const int &_iter = g_MatRowNum );
private:
	//CG法求解的设备端两个临时数组
	double *d_Ax, *d_p;	
};

