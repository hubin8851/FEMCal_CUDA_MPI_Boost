#pragma once
#include "baseSconjugate.h"
class CGSMatCal :
	public CBaseSConjugate	//��Ϊ���ܳ������μ̳У������ڴ�����̳�
{
public:
	CGSMatCal( size_t _RowNum = g_MatRowNum );
	~CGSMatCal(void);
	void	FreeGPUResource();

	void	DeviceMalloc();//����devicemalloc����Ϊ����������ʱ����

	//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ.���麯����ÿ����������д˺����������
	//@_tol:�в�
	//@_iter:��������
	float	ConjugateWithGPU( const float &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER );

private:
	float *d_Ax, *d_p;	//CG�������豸��������ʱ����

#pragma region	thrust�汾
//	thrust::device_ptr<float>	d_vp;	//�����ݶȷ��м����p
//	thrust::device_ptr<float>	d_vAx;	//����Ax
#pragma endregion
};

