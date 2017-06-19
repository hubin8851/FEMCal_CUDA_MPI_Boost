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
	//方程等式右边向量显存分配
	_cuError_id = cudaMalloc((void **)&d_p, m_RowNum*sizeof(double));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	CBaseDConjugate::DeviceMalloc();
	cudaDeviceSynchronize();
}

//各种特征值算法求解,返回计算时间，ms为单位.纯虚函数，每个派生类必有此函数用以求解
//@_tol:残差
//@_iter:迭代次数
double	CGDMatCal::ConjugateWithGPU( const double &_tol, const int &_iter )
{
	double a, b, na;
	double r0, r1;	//过程中的点积结果，随迭代次数变化

	double alpha, beta, dot;	//过程中的系数，随迭代次数变化
	double alpham1;

	alpha = 1.0;
	alpham1 = -1;
	beta = 0.0;
	r0 = 0.0;
	int k = 0;//迭代次数

	m_qaerr1 = 0.0;//误差结果

	// 创建计时的CUDA事件
	//	StopWatchInterface *timer = NULL;
	//	sdkCreateTimer(&timer);
	cudaEvent_t	Event_Start;
	cudaEvent_t Event_Stop;
	_cuError_id = cudaEventCreate(&Event_Start);
	_cuError_id = cudaEventCreate(&Event_Stop);

	//记录共轭梯度法的开始时间
	std::cout<<"双精度共轭梯度法解算..."<<std::endl;
	//	sdkStartTimer(&timer);
	checkCudaErrors(cudaEventRecord(Event_Start, 0));
	//调用cublas的点积函数,具体参数意义参见cublas_library.pdf
	//算法参见计算方法丛书-有限元结构分析并行计算（周树奎、梁维泰）,P179算法5.2	
//	double* _tmp;
//	_tmp = (double*)malloc( m_RowNum*sizeof(double));
	//_cuError_id = cudaMemcpy( _tmp, d_r,	m_RowNum*sizeof(double), cudaMemcpyDeviceToHost );//测试用
	_cusparseError_id = cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, &alpha, Matdescr, 
					d_NonZeroVal,d_iNonZeroRowSort, d_iColSort, 
					d_x, &beta, d_Ax);//y = alpha * op(A) * x  + beta * y,在这里实则计算d_omega = [A]*[r]
//	_cuError_id = cudaMemcpy( _tmp, d_Ax,	m_RowNum*sizeof(double), cudaMemcpyDeviceToHost );//测试用
	_cublasError_id = cublasDaxpy(cublasHandle, m_RowNum, &alpham1, d_Ax, 1, d_r, 1);//y[j]=alpha*x[j]+y[j]
	_cublasError_id = cublasDdot(cublasHandle, m_RowNum, d_r, 1, d_r, 1, &r1);//alpha的分子，点乘
//	_cuError_id = cudaMemcpy( _tmp,	d_r,		m_RowNum*sizeof(double),		cudaMemcpyDeviceToHost);//测试用

	while (r1>_tol*_tol && k <= _iter)
	{
		k++;
		if (1 == k)
		{
			_cublasError_id = cublasDcopy(cublasHandle, m_RowNum, d_r, 1, d_p, 1);
//			_cuError_id = cudaMemcpy( _tmp,	d_p,		m_RowNum*sizeof(double),		cudaMemcpyDeviceToHost);//测试用
		}
		else
		{
			b = r1/r0;
			_cublasError_id = cublasDscal(cublasHandle, m_RowNum, &b, d_p, 1);//向量内每元素乘以常数
			_cublasError_id = cublasDaxpy(cublasHandle, m_RowNum, &alpha, d_r, 1, d_p, 1) ;//y[j]=alpha*x[j]+y[j]
//			_cuError_id = cudaMemcpy( m_CSRInput.h_rhs,			d_p,				m_nnzA*sizeof(double),		cudaMemcpyDeviceToHost);	//测试用
		}
		//y = alpha * op(A) * x  + beta * y,在这里实则计算d_omega = [A]*[r]
		_cusparseError_id = cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, 
			&alpha, Matdescr, 
			d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, d_p, &beta, d_Ax);//y = alpha * op(A) * x  + beta * y
//		_cuError_id = cudaMemcpy( _tmp,	d_Ax,		m_RowNum*sizeof(double),		cudaMemcpyDeviceToHost);//测试用

		checkCudaErrors(cublasDdot(cublasHandle,	m_RowNum, d_p, 1, d_Ax, 1, &dot));//计算alpha的分母
		a = r1/dot;
		checkCudaErrors(cublasDaxpy(cublasHandle,	m_RowNum, &a, d_p, 1, d_x, 1));//计算xi

		na = -a;
		_cublasError_id = cublasDaxpy(cublasHandle,	m_RowNum, &na, d_Ax, 1, d_r, 1);//计算r（i+1）
//		_cuError_id = cudaMemcpy( _tmp,	d_r,		m_RowNum*sizeof(double),		cudaMemcpyDeviceToHost);//测试用

		r0 = r1;
		_cublasError_id = cublasDdot(cublasHandle,	m_RowNum, d_r, 1, d_r, 1, &r1);
		//		cudaThreadSynchronize();
		//		printf("iteration = %3d, residual = %e\n", k, sqrt(r1));//测试用，实例中删除			
	}
	cudaMemcpy( h_x, d_x, m_RowNum*sizeof(double), cudaMemcpyDeviceToHost );
	//记录共轭梯度法的结束时间
	cudaEventRecord(Event_Stop, 0);
	checkCudaErrors(cudaDeviceSynchronize());
	//获得所花时间
	//	sdkStopTimer(&timer);
	checkCudaErrors(cudaEventElapsedTime(&msecUsed, Event_Start, Event_Stop));

	printf("Time= %.3f msec", msecUsed);
	printf("  iteration = %3d, residual = %e \n", k, sqrt(r1));
	/*结果检验*/
	m_qaerr1 = CheckResult(m_nnzA, m_nA, 
		h_NoneZeroVal, h_iColSort, h_iNonZeroRowSort, 
		h_x, h_rhs);
	printf("  残差: %f \n", m_qaerr1);
	printf("  Convergence Test: %s \n", (k <= MAX_GPUITER) ? "OK" : "FAIL");
	m_iters = k;
	return msecUsed;
}