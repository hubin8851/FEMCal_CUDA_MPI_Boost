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
	void	DeviceMalloc();//����devicemalloc����Ϊ������ʱ����

	//��ʼ���������ȶ���
	bool	InitialDescr( cudaStream_t _stream = 0, cusparseMatrixType_t _MatrixType = CUSPARSE_MATRIX_TYPE_GENERAL);
	//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ
	//@_tol:�в�
	//@_iter:��������
	double	ConjugateWithGPU( const double &_tol = ERRORTOLERATE, const int &_iter = MAX_GPUITER );

private:
	//Ԥ�������ݶȷ�����ر���
	cusparseSolveAnalysisInfo_t	infoA;
	cusparseSolveAnalysisInfo_t	info_u;
	cusparseMatDescr_t descrL;
	cusparseMatDescr_t descrU;

	double *h_ILUvals;
	double *d_ILUvals;
	double *d_p, *d_omega, *d_y;	//PCG�������豸����ʱ����
	double *d_zm1, *d_zm2, *d_rm2;
};

