#pragma once
#include "basedconjugate.h"
class CGDMatCal :
	public CBaseDConjugate
{
public:
	CGDMatCal( size_t _RowNum = g_MatRowNum );
	~CGDMatCal(void);
	void	FreeGPUResource();

	void	DeviceMalloc();//����devicemalloc����Ϊ����������ʱ����

	//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ.���麯����ÿ����������д˺����������
	//@_tol:�в�
	//@_iter:��������
	double	ConjugateWithGPU( const double &_tol = ERRORTOLERATE, const int &_iter = g_MatRowNum );
private:
	//CG�������豸��������ʱ����
	double *d_Ax, *d_p;	
};

