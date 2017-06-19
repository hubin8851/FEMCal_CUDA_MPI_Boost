#include "StdAfx.h"
#include "HbxGloFunc.h"
#include "CGDMatCal.h"


CGDMatCal::CGDMatCal( size_t _RowNum ):CBaseDConjugate( _RowNum )
{

}


CGDMatCal::~CGDMatCal(void)
{
	if ( false == m_bCudaFree )
	{
		FreeGPUResource();
	}
}

void CGDMatCal::FreeGPUResource()
{
	if ( false == m_bCudaFree )
	{
		_cuError_id = cudaFree(d_Ax);
		_cuError_id = cudaFree(d_p);
		CBaseDConjugate::FreeGPUResource();
		m_bCudaFree = true;
	}
	else return;
}

void	CGDMatCal::DeviceMalloc()
{
	_cuError_id = cudaMalloc((void **)&d_Ax, m_RowNum*sizeof(double));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	//���̵�ʽ�ұ������Դ����
	_cuError_id = cudaMalloc((void **)&d_p, m_RowNum*sizeof(double));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	CBaseDConjugate::DeviceMalloc();
	cudaDeviceSynchronize();
}

//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ.���麯����ÿ����������д˺����������
//@_tol:�в�
//@_iter:��������
double	CGDMatCal::ConjugateWithGPU( const double &_tol, const int &_iter )
{
	double a, b, na;
	double r0, r1;	//�����еĵ�����������������仯

	double alpha, beta, dot;	//�����е�ϵ��������������仯
	double alpham1;

	alpha = 1.0;
	alpham1 = -1;
	beta = 0.0;
	r0 = 0.0;
	int k = 0;//��������

	m_qaerr1 = 0.0;//�����

	// ������ʱ��CUDA�¼�
	//	StopWatchInterface *timer = NULL;
	//	sdkCreateTimer(&timer);
	cudaEvent_t	Event_Start;
	cudaEvent_t Event_Stop;
	_cuError_id = cudaEventCreate(&Event_Start);
	_cuError_id = cudaEventCreate(&Event_Stop);

	//��¼�����ݶȷ��Ŀ�ʼʱ��
	std::cout<<"˫���ȹ����ݶȷ�����..."<<std::endl;
	//	sdkStartTimer(&timer);
	checkCudaErrors(cudaEventRecord(Event_Start, 0));
	//����cublas�ĵ������,�����������μ�cublas_library.pdf
	//�㷨�μ����㷽������-����Ԫ�ṹ�������м��㣨����������ά̩��,P179�㷨5.2	
//	double* _tmp;
//	_tmp = (double*)malloc( m_RowNum*sizeof(double));
	//_cuError_id = cudaMemcpy( _tmp, d_r,	m_RowNum*sizeof(double), cudaMemcpyDeviceToHost );//������
	_cusparseError_id = cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, &alpha, Matdescr, 
					d_NonZeroVal,d_iNonZeroRowSort, d_iColSort, 
					d_x, &beta, d_Ax);//y = alpha * op(A) * x  + beta * y,������ʵ�����d_omega = [A]*[r]
//	_cuError_id = cudaMemcpy( _tmp, d_Ax,	m_RowNum*sizeof(double), cudaMemcpyDeviceToHost );//������
	_cublasError_id = cublasDaxpy(cublasHandle, m_RowNum, &alpham1, d_Ax, 1, d_r, 1);//y[j]=alpha*x[j]+y[j]
	_cublasError_id = cublasDdot(cublasHandle, m_RowNum, d_r, 1, d_r, 1, &r1);//alpha�ķ��ӣ����
//	_cuError_id = cudaMemcpy( _tmp,	d_r,		m_RowNum*sizeof(double),		cudaMemcpyDeviceToHost);//������

	while (r1>_tol*_tol && k <= _iter)
	{
		k++;
		if (1 == k)
		{
			_cublasError_id = cublasDcopy(cublasHandle, m_RowNum, d_r, 1, d_p, 1);
//			_cuError_id = cudaMemcpy( _tmp,	d_p,		m_RowNum*sizeof(double),		cudaMemcpyDeviceToHost);//������
		}
		else
		{
			b = r1/r0;
			_cublasError_id = cublasDscal(cublasHandle, m_RowNum, &b, d_p, 1);//������ÿԪ�س��Գ���
			_cublasError_id = cublasDaxpy(cublasHandle, m_RowNum, &alpha, d_r, 1, d_p, 1) ;//y[j]=alpha*x[j]+y[j]
//			_cuError_id = cudaMemcpy( m_CSRInput.h_rhs,			d_p,				m_nnzA*sizeof(double),		cudaMemcpyDeviceToHost);	//������
		}
		//y = alpha * op(A) * x  + beta * y,������ʵ�����d_omega = [A]*[r]
		_cusparseError_id = cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, 
			&alpha, Matdescr, 
			d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, d_p, &beta, d_Ax);//y = alpha * op(A) * x  + beta * y
//		_cuError_id = cudaMemcpy( _tmp,	d_Ax,		m_RowNum*sizeof(double),		cudaMemcpyDeviceToHost);//������

		checkCudaErrors(cublasDdot(cublasHandle,	m_RowNum, d_p, 1, d_Ax, 1, &dot));//����alpha�ķ�ĸ
		a = r1/dot;
		checkCudaErrors(cublasDaxpy(cublasHandle,	m_RowNum, &a, d_p, 1, d_x, 1));//����xi

		na = -a;
		_cublasError_id = cublasDaxpy(cublasHandle,	m_RowNum, &na, d_Ax, 1, d_r, 1);//����r��i+1��
//		_cuError_id = cudaMemcpy( _tmp,	d_r,		m_RowNum*sizeof(double),		cudaMemcpyDeviceToHost);//������

		r0 = r1;
		_cublasError_id = cublasDdot(cublasHandle,	m_RowNum, d_r, 1, d_r, 1, &r1);
		//		cudaThreadSynchronize();
		//		printf("iteration = %3d, residual = %e\n", k, sqrt(r1));//�����ã�ʵ����ɾ��			
	}
	cudaMemcpy( h_x, d_x, m_RowNum*sizeof(double), cudaMemcpyDeviceToHost );
	//��¼�����ݶȷ��Ľ���ʱ��
	cudaEventRecord(Event_Stop, 0);
	checkCudaErrors(cudaDeviceSynchronize());
	//�������ʱ��
	//	sdkStopTimer(&timer);
	checkCudaErrors(cudaEventElapsedTime(&msecUsed, Event_Start, Event_Stop));

	printf("Time= %.3f msec", msecUsed);
	printf("  iteration = %3d, residual = %e \n", k, sqrt(r1));
	/*�������*/
	m_qaerr1 = CheckResult(m_nnzA, m_nA, 
		h_NoneZeroVal, h_iColSort, h_iNonZeroRowSort, 
		h_x, h_rhs);
	printf("  �в�: %f \n", m_qaerr1);
	printf("  Convergence Test: %s \n", (k <= MAX_GPUITER) ? "OK" : "FAIL");
	m_iters = k;
	return msecUsed;
}