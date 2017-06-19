#pragma once
#include "BaseSconjugate.h"
class CGSMatCal_UM :	//统一内存下用CG法
	public CBaseSConjugate	//因为可能出现菱形继承，所以在此用虚继承
{
public:
	CGSMatCal_UM( size_t _RowNum = g_MatRowNum );
	~CGSMatCal_UM(void);
	void	FreeGPUResource();

	//获取GPU设备属性，返回所有的设备数并在函数内获得最优设备号
	virtual	size_t	GetGPUDeviceQuery();
	void		HostMalloc();
	void		DeviceMalloc();//重载devicemalloc，使用统一内存
	HBXDef::DataAloc_t	MemCpy(HBXDef::CopyType_t _temp = HBXDef::CopyType_t::HostToDevice);	
	//各种特征值算法求解,返回计算时间，ms为单位.纯虚函数，每个派生类必有此函数用以求解
	//@_tol:残差
	//@_iter:迭代次数
	float	ConjugateWithGPU( const float &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER );

private:
	float *d_Ax, *d_p;	//CG法求解的设备端两个临时数组
};

