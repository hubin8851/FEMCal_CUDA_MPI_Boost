#include "StdAfx.h"
#include "CGSMatCal_UM.h"


CGSMatCal_UM::CGSMatCal_UM( size_t _RowNum ):CBaseSConjugate( _RowNum )
{

}


CGSMatCal_UM::~CGSMatCal_UM(void)
{

}

void	CGSMatCal_UM::FreeGPUResource()
{
	if ( false == m_bCudaFree)
	{
		cudaFree( h_NoneZeroVal );
		cudaFree( h_iColSort );
		cudaFree( h_iNonZeroRowSort );
		cudaFree( h_x );
		cudaFree( h_rhs );
		cudaFree(d_r);
		cudaFree(d_p);
		cudaFree(d_Ax);
		m_bCudaFree = true;
	}
	else return;
}

//��ȡGPU�豸���ԣ��������е��豸�����ں����ڻ�������豸��
size_t	CGSMatCal_UM::GetGPUDeviceQuery()
{
	// This will pick the best possible CUDA capable device
	cudaDeviceProp deviceProp;
	int deviceCount = 0;
	cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

	if (error_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
		printf("Result = FAIL\n");
		return 0;
	}
	if (!deviceProp.managedMemory) { //�ж��Ƿ���ͳһ�ڴ�
		fprintf(stderr, "Unified Memory not supported on this device\n");
		exit(EXIT_WAIVED);
	}
	if (deviceCount == 0)
	{
		printf("\nThere are no available device(s) that support CUDA\n");
		return 0;
	}
	else
	{
		printf("Detected %d CUDA Capable device(s)\n", deviceCount);
	}
	int dev, driverVersion = 0, runtimeVersion = 0;
	cudaDeviceProp	_BestDevice;//��õ�GPU�����Խṹ��
	for (dev = 0; dev < deviceCount; ++dev)
	{
		cudaSetDevice(dev);
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, dev);

		printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

		// Console log
		cudaDriverGetVersion(&driverVersion);
		cudaRuntimeGetVersion(&runtimeVersion);
		printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n", driverVersion/1000, (driverVersion%100)/10, runtimeVersion/1000, (runtimeVersion%100)/10);
		printf("  CUDA Capability Major/Minor version number:    %d.%d\n", deviceProp.major, deviceProp.minor);

		char msg[256];
		sprintf_s(msg, "  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n",
			(float)deviceProp.totalGlobalMem/1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);
		printf("%s", msg);

		printf("  (%2d) Multiprocessors, (%3d) CUDA Cores/MP:     %d CUDA Cores\n",
			deviceProp.multiProcessorCount,
			_ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
			_ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount);
		printf("  GPU Max Clock rate:                            %.0f MHz (%0.2f GHz)\n", deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);
		if ((_BestDevice.major<=deviceProp.major)&&(_BestDevice.minor<=deviceProp.minor))
		{
			m_DevID = dev;
			// ��SM���߳����Ƶ����ӣ��߳̿��ڵ��߳���������Ӧ����
			int m_BlockSize = (deviceProp.major < 2) ? 16 : 32;
		}
	}
	return deviceCount;
}

void	CGSMatCal_UM::HostMalloc()
{
	_cuError_id = cudaMallocManaged((void **)&h_NoneZeroVal,		sizeof(float)*(m_nnzA+1));
	_cuError_id = cudaMallocManaged((void **)&h_iColSort,			sizeof(int)*m_nnzA);
	_cuError_id = cudaMallocManaged((void **)&h_iNonZeroRowSort, sizeof(float)*m_nA);
	_cuError_id = cudaMallocManaged((void **)&h_x,						sizeof(float)*m_RowNum);
	_cuError_id = cudaMallocManaged((void **)&h_rhs,					sizeof(float)*m_RowNum);
}

//ʹ��ͳһ�ڴ�
void	CGSMatCal_UM::DeviceMalloc()
{
	// CG ���õ���ʱ����
	_cuError_id = cudaMallocManaged((void **)&d_r, m_RowNum*sizeof(float));
	_cuError_id = cudaMallocManaged((void **)&d_p, m_RowNum*sizeof(float));
	_cuError_id = cudaMallocManaged((void **)&d_Ax, m_RowNum*sizeof(float));
	cudaDeviceSynchronize();
	m_bCudaFree = false;
}

HBXDef::DataAloc_t	CGSMatCal_UM::MemCpy( HBXDef::CopyType_t _temp )
{
	if (HBXDef::HostToDevice == _temp)
	{
		for (int i=0; i < m_RowNum; i++)
		{
			d_r[i] = h_rhs[i];
		}
		return HBXDef::DataAloc_t::DATAINGPU;
	}
	else if (HBXDef::DeviceToHost == _temp)
	{
		for (int i=0; i < m_RowNum; i++)
		{
			h_rhs[i] = d_r[i];
		}
		return HBXDef::DataAloc_t::DATAINMEM;
	}
	else return HBXDef::DataAloc_t::INVALID;
}

//��������ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ.���麯����ÿ����������д˺����������
//@_tol:�в�
//@_iter:��������
float	CGSMatCal_UM::ConjugateWithGPU( const float &_tol, const int &_iter )
{
	float a, b, na;
	float r0, r1;	//�����еĵ�����������������仯
	float alpha, beta, dot;	//�����е�ϵ��������������仯
	float alpham1;

	alpha = 1.0;
	alpham1 = -1.0;
	beta = 0.0;
	r0 = 0.;
	int k = 1;//��������

	float	msecUsed = 0.0f;//ʱ��ͳ��float�ͱ���
	m_qaerr1 = 0.0;//�����

	// ������ʱ��CUDA�¼�
	StopWatchInterface *timer = NULL;
	sdkCreateTimer(&timer);
	cudaEvent_t	Event_Start;
	cudaEvent_t Event_Stop;
	_cuError_id = cudaEventCreate(&Event_Start);
	_cuError_id = cudaEventCreate(&Event_Stop);

	//��¼�����ݶȷ��Ŀ�ʼʱ��
	std::cout<<"������ͳһ�ڴ�Ĺ����ݶȷ�����..."<<std::endl;
	sdkStartTimer(&timer);
	checkCudaErrors(cudaEventRecord(Event_Start, 0));
	//����cublas�ĵ������,�����������μ�cublas_library.pdf
	//�㷨�μ����㷽������-����Ԫ�ṹ�������м��㣨����������ά̩��,P179�㷨5.2	
//	float* _tmp;//��ʱ���Ա���
//	_tmp = (float*)malloc( m_RowNum*sizeof(float));//��ʱ���Ա���
	cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, &alpha, Matdescr, 
		h_NoneZeroVal,h_iNonZeroRowSort, h_iColSort, 
		h_x, &beta, d_Ax);
//	_cuError_id = cudaMemcpy( _tmp, d_Ax, m_RowNum*sizeof(float), cudaMemcpyDeviceToHost );	//��ʱ���Ա���
	cublasSaxpy(cublasHandle, m_RowNum, &alpham1, d_Ax, 1, d_r, 1);
	cublasSdot(cublasHandle, m_RowNum, d_r, 1, d_r, 1, &r1);//alpha�ķ���

	while (r1>_tol*_tol && k <= _iter)
	{
		if (1 == k)
		{
			cublasScopy(cublasHandle, m_RowNum, d_r, 1, d_p, 1);
		}
		else
		{
			b = r1/r0;
			_cublasError_id = cublasSscal(cublasHandle, m_RowNum, &b, d_p, 1);//������ÿԪ�س��Գ���
			_cublasError_id = cublasSaxpy(cublasHandle, m_RowNum, &alpha, d_r, 1, d_p, 1) ;//y[j]=alpha*x[j]+y[j]
		}
		//y = alpha * op(A) * x  + beta * y,������ʵ�����d_omega = [A]*[r]
		_cusparseError_id = cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, 
			&alpha, Matdescr, 
			h_NoneZeroVal, h_iNonZeroRowSort, h_iColSort, d_p, &beta, d_Ax);

		_cublasError_id = cublasSdot(cublasHandle,	m_RowNum, d_p, 1, d_Ax, 1, &dot);//����alpha�ķ�ĸ
		a = r1/dot;
		_cublasError_id = cublasSaxpy(cublasHandle,	m_RowNum, &a, d_p, 1, h_x, 1);//����xi

		na = -a;
		_cublasError_id = cublasSaxpy(cublasHandle,	m_RowNum, &na, d_Ax, 1, d_r, 1);//����r��i+1��

		r0 = r1;
		_cublasError_id = cublasSdot(cublasHandle,	m_RowNum, d_r, 1, d_r, 1, &r1);
//		cudaThreadSynchronize();
//		printf("iteration = %3d, residual = %e\n", k, sqrt(r1));//�����ã�ʵ����ɾ��			
		k++;
	}
	//��¼�����ݶȷ��Ľ���ʱ��
	cudaEventRecord(Event_Stop, 0);
	checkCudaErrors(cudaDeviceSynchronize());
	//�������ʱ��
	sdkStopTimer(&timer);
	_cuError_id = cudaEventElapsedTime(&msecUsed, Event_Start, Event_Stop);

	printf("Time= %.3f msec", msecUsed);
	printf("  iteration = %3d, residual = %e \n", k, sqrt(r1));
	//	m_qaerr1 = CheckResult(m_nnzA, m_nA, 
	//										h_NoneZeroVal, h_iColSort, h_iNonZeroRowSort, 
	//										h_x, h_rhs);
	printf("  �в�: %f \n", m_qaerr1);
	printf("  Convergence Test: %s \n", (k <= MAX_GPUITER) ? "OK" : "FAIL");
	m_iters = k;
	return msecUsed;
}