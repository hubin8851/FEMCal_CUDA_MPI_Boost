#pragma once
#include "BaseSconjugate.h"
class CGSMatCal_UM :	//ͳһ�ڴ�����CG��
	public CBaseSConjugate	//��Ϊ���ܳ������μ̳У������ڴ�����̳�
{
public:
	CGSMatCal_UM( size_t _RowNum = g_MatRowNum );
	~CGSMatCal_UM(void);
	void	FreeGPUResource();

	//��ȡGPU�豸���ԣ��������е��豸�����ں����ڻ�������豸��
	virtual	size_t	GetGPUDeviceQuery();
	void		HostMalloc();
	void		DeviceMalloc();//����devicemalloc��ʹ��ͳһ�ڴ�
	HBXDef::DataAloc_t	MemCpy(HBXDef::CopyType_t _temp = HBXDef::CopyType_t::HostToDevice);	
	//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ.���麯����ÿ����������д˺����������
	//@_tol:�в�
	//@_iter:��������
	float	ConjugateWithGPU( const float &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER );

private:
	float *d_Ax, *d_p;	//CG�������豸��������ʱ����
};

