#include "StdAfx.h"
#include "CGSMatCal.h"
#include "HbxGloFunc.h"


CGSMatCal::CGSMatCal( size_t _RowNum ):CBaseSConjugate( _RowNum )
{
//	int i=0;//�����öϵ�
#pragma region thrust�汾��������
// 	thrust::device_ptr<float> d_vp = thrust::device_malloc<float>(m_RowNum);
// 	thrust::device_ptr<float> d_vAx = thrust::device_malloc<float>(m_RowNum);
// 	//��ָ������ȷ��ַ
// 	d_p = thrust::raw_pointer_cast( d_vp );		//�����ݶȷ��м����p
// 	d_Ax = thrust::raw_pointer_cast( d_vAx );	//����Ax
#pragma endregion
}


CGSMatCal::~CGSMatCal(void)
{
	if ( false == m_bCudaFree )
	{
		FreeGPUResource();
	}
}

void CGSMatCal::FreeGPUResource()
{
	if ( false == m_bCudaFree )
	{
		_cuError_id = cudaFree(d_Ax);
		_cuError_id = cudaFree(d_p);
		CBaseSConjugate::FreeGPUResource();
		m_bCudaFree = true;
	}
	else return;
}

void		CGSMatCal::DeviceMalloc()
{
	_cuError_id = cudaMalloc((void **)&d_Ax, m_RowNum*sizeof(float));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	//���̵�ʽ�ұ������Դ����
	_cuError_id = cudaMalloc((void **)&d_p, m_RowNum*sizeof(float));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	CBaseSConjugate::DeviceMalloc();
	cudaDeviceSynchronize();
}


//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ.���麯����ÿ����������д˺����������
//@_tol:�в�
//@_iter:��������
float	CGSMatCal::ConjugateWithGPU( const float &_tol, const int &_iter )
{
	float a, b, na;
	float r0, r1;	//�����еĵ�����������������仯

	float alpha, beta, dot;	//�����е�ϵ��������������仯
	float alpham1;

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
	 std::cout<<"�����ȹ����ݶȷ�����..."<<std::endl;
//	sdkStartTimer(&timer);
	checkCudaErrors(cudaEventRecord(Event_Start, 0));
	//����cublas�ĵ������,�����������μ�cublas_library.pdf
	//�㷨�μ����㷽������-����Ԫ�ṹ�������м��㣨����������ά̩��,P179�㷨5.2	
	//log��20170322.ʹ�õ����㷨5.1...
//	 float* _tmp;
//	 _tmp = (float*)malloc( m_RowNum*sizeof(float));
//	 _cuError_id = cudaMemcpy( _tmp, d_r,	m_RowNum*sizeof(float), cudaMemcpyDeviceToHost );//������
	_cusparseError_id = cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, &alpha, Matdescr, 
										d_NonZeroVal,d_iNonZeroRowSort, d_iColSort, 
										d_x, &beta, d_Ax);//y = alpha * op(A) * x  + beta * y,������ʵ�����d_omega = [A]*[r]
//	_cuError_id = cudaMemcpy( _tmp, d_Ax,	m_RowNum*sizeof(float), cudaMemcpyDeviceToHost );//������
	_cublasError_id = cublasSaxpy(cublasHandle, m_RowNum, &alpham1, d_Ax, 1, d_r, 1);//y[j]=alpha*x[j]+y[j]
	_cublasError_id = cublasSdot(cublasHandle, m_RowNum, d_r, 1, d_r, 1, &r1);//alpha�ķ��ӣ����
//	_cuError_id = cudaMemcpy( _tmp,	d_r,		m_RowNum*sizeof(float),		cudaMemcpyDeviceToHost);//������

	while (r1>_tol*_tol && k <= _iter)
	{
		k++;
		if (1 == k)
		{
			_cublasError_id = cublasScopy(cublasHandle, m_RowNum, d_r, 1, d_p, 1);
		}
		else
		{
			b = r1/r0;			//����(4)
			_cublasError_id = cublasSscal(cublasHandle, m_RowNum, &b, d_p, 1);//������ÿԪ�س��Գ�������(5)
			_cublasError_id = cublasSaxpy(cublasHandle, m_RowNum, &alpha, d_r, 1, d_p, 1) ;//y[j]=alpha*x[j]+y[j],����(6)
//			_cuError_id = cudaMemcpy( m_CSRInput.h_rhs,			d_p,				m_nnzA*sizeof(float),		cudaMemcpyDeviceToHost);	//������
		}
		//y = alpha * op(A) * x  + beta * y,������ʵ�����d_omega = [A]*[r]
		_cusparseError_id = cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, 
																&alpha, Matdescr, 
																d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, d_p, &beta, d_Ax);
		//_cuError_id = cudaMemcpy( _tmp,	d_Ax,		m_RowNum*sizeof(float),		cudaMemcpyDeviceToHost);//������

		checkCudaErrors(cublasSdot(cublasHandle,	m_RowNum, d_p, 1, d_Ax, 1, &dot));//����alpha�ķ�ĸ
		a = r1/dot;
		checkCudaErrors(cublasSaxpy(cublasHandle,	m_RowNum, &a, d_p, 1, d_x, 1));//����xi

		na = -a;
		_cublasError_id = cublasSaxpy(cublasHandle,	m_RowNum, &na, d_Ax, 1, d_r, 1);//����r��i+1��
//		_cuError_id = cudaMemcpy( _tmp,	d_r,		m_RowNum*sizeof(float),		cudaMemcpyDeviceToHost);//������

		r0 = r1;
		_cublasError_id = cublasSdot(cublasHandle,	m_RowNum, d_r, 1, d_r, 1, &r1);
//		cudaThreadSynchronize();
//		printf("iteration = %3d, residual = %e\n", k, sqrt(r1));//�����ã�ʵ����ɾ��			
	}
	cudaMemcpy( h_x, d_x, m_RowNum*sizeof(float), cudaMemcpyDeviceToHost );
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