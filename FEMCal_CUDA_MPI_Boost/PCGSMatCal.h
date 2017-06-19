#pragma once
#include "baseSconjugate.h"
#include "CGSMatCal.h"

class PCGSMatCal:
	public CBaseSConjugate
{
public:
	PCGSMatCal( size_t _RowNum = g_MatRowNum );
	~PCGSMatCal(void);
	void FreeGPUResource();

	void	genLaplace( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath = "" );
	void	DeviceMalloc();//����devicemalloc����Ϊ������ʱ����

	//��ʼ���������ȶ���
	bool		InitialDescr(cusparseMatrixType_t _MatrixType = CUSPARSE_MATRIX_TYPE_GENERAL);
	//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ
	//@_tol:�в�
	//@_iter:��������
	float	ConjugateWithGPU( const float &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER );

private:
	//Ԥ�������ݶȷ�����ر���
	cusparseSolveAnalysisInfo_t	infoA;
	cusparseSolveAnalysisInfo_t	info_u;
	cusparseMatDescr_t descrL;
	cusparseMatDescr_t descrU;

	float *h_ILUvals;
	float *d_ILUvals;
	float *d_p, *d_omega, *d_y;	//PCG�������豸����ʱ����
	float *d_zm1, *d_zm2, *d_rm2;
	
};

