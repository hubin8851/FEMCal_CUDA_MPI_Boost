#pragma once	//2017-03��hbx
#include "basedconjugate.h"
class BiCGStabD :
	public CBaseDConjugate
{
public:
	BiCGStabD( size_t _RowNum = g_MatRowNum );
	~BiCGStabD(void);
	void	FreeGPUResource(){};

	void	HostMalloc();//����HostMalloc����Ϊ������Ҫ��������Ĳ���
	void	DeviceMalloc();//����devicemalloc����Ϊ������ʱ����

	//�豸�������˵��ڴ濽��,�ڴ˴����أ���Ϊ��Ҫ��������Ĳ���
	HBXDef::DataAloc_t		MemCpy(HBXDef::CopyType_t _temp = HBXDef::CopyType_t::HostToDevice);

	//��ʼ���������ȶ���
	bool		InitialDescr( cudaStream_t _stream = 0, cusparseMatrixType_t _MatrixType = CUSPARSE_MATRIX_TYPE_GENERAL );

	//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ.
	//@_tol:�в�
	//@_iter:��������
	double	ConjugateWithGPU( const double &_tol = ERRORTOLERATE, const int &_iter = g_MatRowNum );
private:
	//@maxit:����������
	//@tol:��������
	//@ttt_sv:�����Ԥ��������ʱ��
	void gpu_pbicgstab( int maxit, double tol, double ttt_sv);

	//Ԥ�������ݶȷ�����ر���
	cusparseSolveAnalysisInfo_t	info_l;
	cusparseSolveAnalysisInfo_t	info_u;
	cusparseMatDescr_t	descra;
	cusparseMatDescr_t	descrm;
	cusparseOperation_t	_trans;



	int arraySizeX, arraySizeF, arraySizeR, arraySizeRW, arraySizeP,  arraySizePW, arraySizeS, arraySizeT, arraySizeV;
	//�����˼����е��м����
	double *h_ptrF;
	double *h_ptrR, *h_ptrRW;
	double *h_ptrP, *h_ptrPW;
	double *h_ptrS, *h_ptrT, *h_ptrV, *h_ptrTX; 
	double *h_ptrMval;
	//�豸�˼����е��м����
	double *d_ptrF;
	double *d_ptrR, *d_ptrRW;
	double *d_ptrP, *d_ptrPW;
	double *d_ptrS, *d_ptrT, *d_ptrV;
	double *d_ptrMval;
};

