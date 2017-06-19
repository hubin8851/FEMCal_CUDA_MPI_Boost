#pragma once
#include "baseSconjugate.h"
class CGSMatCal :
	public CBaseSConjugate	//因为可能出现菱形继承，所以在此用虚继承
{
public:
	CGSMatCal( size_t _RowNum = g_MatRowNum );
	~CGSMatCal(void);
	void	FreeGPUResource();

	void	DeviceMalloc();//重载devicemalloc，因为多了两个临时数组

	//各种特征值算法求解,返回计算时间，ms为单位.纯虚函数，每个派生类必有此函数用以求解
	//@_tol:残差
	//@_iter:迭代次数
	float	ConjugateWithGPU( const float &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER );

private:
	float *d_Ax, *d_p;	//CG法求解的设备端两个临时数组

#pragma region	thrust版本
//	thrust::device_ptr<float>	d_vp;	//共轭梯度法中间变量p
//	thrust::device_ptr<float>	d_vAx;	//就是Ax
#pragma endregion
};

