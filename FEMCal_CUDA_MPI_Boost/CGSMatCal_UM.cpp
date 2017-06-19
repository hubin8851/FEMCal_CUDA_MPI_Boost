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

//获取GPU设备属性，返回所有的设备数并在函数内获得最优设备号
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
	if (!deviceProp.managedMemory) { //判断是否有统一内存
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
	cudaDeviceProp	_BestDevice;//最好的GPU的属性结构体
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
			// 因SM内线程限制的增加，线程块内的线程数量可相应增加
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

//使用统一内存
void	CGSMatCal_UM::DeviceMalloc()
{
	// CG 所用的临时变量
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

//各种特征值算法求解,返回计算时间，ms为单位.纯虚函数，每个派生类必有此函数用以求解
//@_tol:残差
//@_iter:迭代次数
float	CGSMatCal_UM::ConjugateWithGPU( const float &_tol, const int &_iter )
{
	float a, b, na;
	float r0, r1;	//过程中的点积结果，随迭代次数变化
	float alpha, beta, dot;	//过程中的系数，随迭代次数变化
	float alpham1;

	alpha = 1.0;
	alpham1 = -1.0;
	beta = 0.0;
	r0 = 0.;
	int k = 1;//迭代次数

	float	msecUsed = 0.0f;//时间统计float型变量
	m_qaerr1 = 0.0;//误差结果

	// 创建计时的CUDA事件
	StopWatchInterface *timer = NULL;
	sdkCreateTimer(&timer);
	cudaEvent_t	Event_Start;
	cudaEvent_t Event_Stop;
	_cuError_id = cudaEventCreate(&Event_Start);
	_cuError_id = cudaEventCreate(&Event_Stop);

	//记录共轭梯度法的开始时间
	std::cout<<"单精度统一内存的共轭梯度法解算..."<<std::endl;
	sdkStartTimer(&timer);
	checkCudaErrors(cudaEventRecord(Event_Start, 0));
	//调用cublas的点积函数,具体参数意义参见cublas_library.pdf
	//算法参见计算方法丛书-有限元结构分析并行计算（周树奎、梁维泰）,P179算法5.2	
//	float* _tmp;//临时测试变量
//	_tmp = (float*)malloc( m_RowNum*sizeof(float));//临时测试变量
	cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, &alpha, Matdescr, 
		h_NoneZeroVal,h_iNonZeroRowSort, h_iColSort, 
		h_x, &beta, d_Ax);
//	_cuError_id = cudaMemcpy( _tmp, d_Ax, m_RowNum*sizeof(float), cudaMemcpyDeviceToHost );	//临时测试变量
	cublasSaxpy(cublasHandle, m_RowNum, &alpham1, d_Ax, 1, d_r, 1);
	cublasSdot(cublasHandle, m_RowNum, d_r, 1, d_r, 1, &r1);//alpha的分子

	while (r1>_tol*_tol && k <= _iter)
	{
		if (1 == k)
		{
			cublasScopy(cublasHandle, m_RowNum, d_r, 1, d_p, 1);
		}
		else
		{
			b = r1/r0;
			_cublasError_id = cublasSscal(cublasHandle, m_RowNum, &b, d_p, 1);//向量内每元素乘以常数
			_cublasError_id = cublasSaxpy(cublasHandle, m_RowNum, &alpha, d_r, 1, d_p, 1) ;//y[j]=alpha*x[j]+y[j]
		}
		//y = alpha * op(A) * x  + beta * y,在这里实则计算d_omega = [A]*[r]
		_cusparseError_id = cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, 
			&alpha, Matdescr, 
			h_NoneZeroVal, h_iNonZeroRowSort, h_iColSort, d_p, &beta, d_Ax);

		_cublasError_id = cublasSdot(cublasHandle,	m_RowNum, d_p, 1, d_Ax, 1, &dot);//计算alpha的分母
		a = r1/dot;
		_cublasError_id = cublasSaxpy(cublasHandle,	m_RowNum, &a, d_p, 1, h_x, 1);//计算xi

		na = -a;
		_cublasError_id = cublasSaxpy(cublasHandle,	m_RowNum, &na, d_Ax, 1, d_r, 1);//计算r（i+1）

		r0 = r1;
		_cublasError_id = cublasSdot(cublasHandle,	m_RowNum, d_r, 1, d_r, 1, &r1);
//		cudaThreadSynchronize();
//		printf("iteration = %3d, residual = %e\n", k, sqrt(r1));//测试用，实例中删除			
		k++;
	}
	//记录共轭梯度法的结束时间
	cudaEventRecord(Event_Stop, 0);
	checkCudaErrors(cudaDeviceSynchronize());
	//获得所花时间
	sdkStopTimer(&timer);
	_cuError_id = cudaEventElapsedTime(&msecUsed, Event_Start, Event_Stop);

	printf("Time= %.3f msec", msecUsed);
	printf("  iteration = %3d, residual = %e \n", k, sqrt(r1));
	//	m_qaerr1 = CheckResult(m_nnzA, m_nA, 
	//										h_NoneZeroVal, h_iColSort, h_iNonZeroRowSort, 
	//										h_x, h_rhs);
	printf("  残差: %f \n", m_qaerr1);
	printf("  Convergence Test: %s \n", (k <= MAX_GPUITER) ? "OK" : "FAIL");
	m_iters = k;
	return msecUsed;
}