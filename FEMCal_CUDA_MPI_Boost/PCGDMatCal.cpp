#include "StdAfx.h"
#include "HbxGloFunc.h"
#include "PCGDMatCal.h"

//PCG法双精度运算类
PCGDMatCal::PCGDMatCal( size_t _RowNum ):CBaseDConjugate( _RowNum )
{
	//预处理共轭梯度法部分参数初始化
	infoA = 0;
	info_u = 0;
	descrL = 0;
	descrU = 0;

	h_ILUvals = nullptr;
	d_ILUvals = nullptr;
	d_zm1 = nullptr;
	d_zm2 = nullptr;
	d_rm2 = nullptr;
}


PCGDMatCal::~PCGDMatCal(void)
{

}

void PCGDMatCal::FreeGPUResource()
{
	if ( false == m_bCudaFree )
	{
		/* 销毁LU分解先关描述信息 */
		checkCudaErrors(cusparseDestroySolveAnalysisInfo(infoA));
		checkCudaErrors(cusparseDestroySolveAnalysisInfo(info_u));
		/* 释放设备内存 */
		_cuError_id = cudaFree(d_y);
		_cuError_id = cudaFree(d_p);
		_cuError_id = cudaFree(d_omega);
		_cuError_id = cudaFree(d_ILUvals);
		_cuError_id = cudaFree(d_zm1);
		_cuError_id = cudaFree(d_zm2);
		_cuError_id = cudaFree(d_rm2);
		CBaseDConjugate::FreeGPUResource();
		m_bCudaFree = true;
	}
	else return;
}

void	PCGDMatCal::genLaplace( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath )
{
	using namespace std;
	cout<<"生成PCG法所用矩阵，原始数据在内存中抹去..."<<endl;
	//因为强制设为三对角矩阵，故m_nnA,m_nA,rownum均需要重新赋值
	m_RowNum = _genRowNum;
	m_nnzA = 5*m_RowNum-4*(int)sqrt((double)m_RowNum);
	m_nA = m_RowNum+1;
	//重新分配内存,并完成内存内h_x和h_rhs的置0
	HostMalloc();

	int n = (int)sqrt((double)m_RowNum);
	assert(n*n == m_RowNum);
	printf("laplace dimension = %d\n", n);
	int idx = 0;

	// loop over degrees of freedom
	for (int i=0; i<m_RowNum; i++)
	{
		int ix = i%n;
		int iy = i/n;

		h_iNonZeroRowSort[i] = idx;

		// up
		if (iy > 0)
		{
			h_NoneZeroVal[idx] = 1.0;
			h_iColSort[idx] = i-n;
			idx++;
		}
		else
		{
			h_rhs[i] -= 1.0;
		}

		// left
		if (ix > 0)
		{
			h_NoneZeroVal[idx] = 1.0;
			h_iColSort[idx] = i-1;
			idx++;
		}
		else
		{
			h_rhs[i] -= 0.0;
		}

		// center
		h_NoneZeroVal[idx] = -4.0;
		h_iColSort[idx]=i;
		idx++;

		//right
		if (ix  < n-1)
		{
			h_NoneZeroVal[idx] = 1.0;
			h_iColSort[idx] = i+1;
			idx++;
		}
		else
		{
			h_rhs[i] -= 0.0;
		}

		//down
		if (iy  < n-1)
		{
			h_NoneZeroVal[idx] = 1.0;
			h_iColSort[idx] = i+n;
			idx++;
		}
		else
		{
			h_rhs[i] -= 0.0;
		}
	}
	h_iNonZeroRowSort[m_RowNum] = idx;
	//CSR格式矩阵是否存储的相关代码，20161103
	if ( _bsave )
	{
		using namespace boost::gregorian;
		using namespace boost::posix_time;
		boost::filesystem::path _SavePath;
		if ("" == _savepath)
		{
			_SavePath = "..\\data\\PCGCheck\\";
		}
		else 
		{
			_SavePath = _savepath;
		}
		if (exists(_SavePath))
		{
			if ( boost::filesystem::is_empty(_SavePath) )
			{
				remove(_SavePath);
			}
		}
		else
		{
			remove_all(_SavePath);
		}
		assert(!exists(_SavePath));
		create_directories(_SavePath);
		HBXDef::_ExportCSRdataToFile( h_NoneZeroVal, h_iColSort, h_iNonZeroRowSort,
			m_nnzA, m_RowNum, "NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data",
			_SavePath );
	}
}

void	PCGDMatCal::DeviceMalloc()
{
	_cuError_id = cudaMalloc((void **)&d_y, m_RowNum*sizeof(double));
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
	_cuError_id = cudaMalloc((void **)&d_omega, m_RowNum*sizeof(double));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	CBaseDConjugate::DeviceMalloc();
	cudaDeviceSynchronize();
	m_bCudaFree = false;
}

//初始化描述器等对象
bool	PCGDMatCal::InitialDescr( cudaStream_t _stream, cusparseMatrixType_t _MatrixType )
{
	m_stream = _stream;
	h_ILUvals = (double *) malloc(m_nnzA*sizeof(double));
	_cuError_id = cudaMalloc((void **)&d_ILUvals, m_nnzA*sizeof(double));
	_cuError_id = cudaMalloc((void **)&d_zm1, m_RowNum*sizeof(double));
	_cuError_id = cudaMalloc((void **)&d_zm2, m_RowNum*sizeof(double));
	_cuError_id = cudaMalloc((void **)&d_rm2, m_RowNum*sizeof(double));

	CBaseDConjugate::InitialDescr( _stream, _MatrixType );
	/* 创建矩阵A的分析信息 */
	_cusparseError_id = cusparseCreateSolveAnalysisInfo(&infoA);
	if ( CUSPARSE_STATUS_SUCCESS != _cusparseError_id)
	{
		printf("cusparseCreateSolveAnalysisInfo returned %d\n-> %s\n", (int)_cusparseError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	// 创建对该非转置矩阵的描述器,相当于绑定
	_cusparseError_id = cusparseDcsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		m_RowNum, m_nnzA, Matdescr, 
		d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, infoA);
	if ( CUSPARSE_STATUS_SUCCESS != _cusparseError_id)
	{
		printf("cusparseDcsrsv_analysis returned %d\n-> %s\n", (int)_cusparseError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}

	/* 将A数据作为ILU0的值传入*/
	cudaMemcpy(d_ILUvals, d_NonZeroVal, m_nnzA*sizeof(double), cudaMemcpyDeviceToDevice);
	/* 描述器，以CSR格式存储的A矩阵的不完全LU分解以便于后期计算 */
	_cusparseError_id = cusparseDcsrilu0(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, 
		Matdescr, d_ILUvals, d_iNonZeroRowSort, d_iColSort, infoA);
	/* 创建对LU分解预处理的相关信息 */
	_cusparseError_id = cusparseCreateSolveAnalysisInfo(&info_u);

	_cusparseError_id = cusparseCreateMatDescr(&descrL);
	_cusparseError_id = cusparseSetMatType(descrL,CUSPARSE_MATRIX_TYPE_GENERAL);
	_cusparseError_id = cusparseSetMatIndexBase(descrL, m_CSRIndexBase );
	_cusparseError_id = cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
	_cusparseError_id = cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);

	_cusparseError_id = cusparseCreateMatDescr(&descrU);
	_cusparseError_id = cusparseSetMatType(descrU,CUSPARSE_MATRIX_TYPE_GENERAL);
	_cusparseError_id = cusparseSetMatIndexBase(descrU, m_CSRIndexBase );
	_cusparseError_id = cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
	_cusparseError_id = cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
	_cusparseError_id = cusparseDcsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 
												m_RowNum, m_nnzA, descrU, d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, info_u);

	return true;
}

//GPU加速的预处理共轭梯度法/* 使用ILU的预处理共轭梯度法.*/
//算法参见计算方法丛书-有限元结构分析并行计算（周树奎、梁维泰）,P183算法5.4
//PCG特征值算法求解,返回计算时间，ms为单位
//@_tol:残差
//@_iter:最大迭代次数
double	PCGDMatCal::ConjugateWithGPU( const double &_tol, const int &_iter )
{
	printf("\n使用LU分解的预处理共轭梯度法，starting...: \n");
	// 创建计时的CUDA事件
	cudaEvent_t	Event_Start, Event_PreStart, Event_MemStart;
	cudaEvent_t Event_Stop, Event_PreEnd, Event_MemEnd;
	checkCudaErrors(cudaEventCreate(&Event_Start));
	checkCudaErrors(cudaEventCreate(&Event_Stop));
	checkCudaErrors(cudaEventCreate(&Event_PreStart));
	checkCudaErrors(cudaEventCreate(&Event_PreEnd));
	checkCudaErrors(cudaEventCreate(&Event_MemStart));
	checkCudaErrors(cudaEventCreate(&Event_MemEnd));

	double r0, r1, alpha, beta;
	int k = 0;//迭代次数
	double numerator, denominator, nalpha;
	const double floatone = 1.0;
	const double floatzero = 0.0;
	float	msecUsed = 0.0f;	//时间统计float型变量
	float  CopysecUsed = 0.0f;	//内存拷贝时间统计量
	m_qaerr1 = 0.0;//误差结果
	int	nzILU0 = 2*m_RowNum - 1;//具体在哪赋值待定，先放在这

	checkCudaErrors(cudaEventRecord(Event_Start, NULL));	//开始计时


	_cuError_id = cudaMemcpy(d_r,		h_rhs,	m_RowNum*sizeof(double), cudaMemcpyHostToDevice);
	_cuError_id = cudaMemcpy(d_x,		h_x,		m_RowNum*sizeof(double), cudaMemcpyHostToDevice);

	//记录内存拷贝的结束时间
	checkCudaErrors(cudaEventRecord(Event_MemEnd, NULL));
	// 等待所有的停止事件完成
	checkCudaErrors(cudaEventSynchronize(Event_MemEnd));
	checkCudaErrors(cudaEventElapsedTime(&CopysecUsed, Event_Start, Event_MemEnd));

	cublasDdot(cublasHandle, m_RowNum, d_r, 1, d_r, 1, &r1);
	//	cudaMemcpy(m_CSRInput.h_x, d_r, m_RowNum*sizeof(double), cudaMemcpyDeviceToHost);//测试用的，540M转560发送错误，供测试用...

	while (r1>_tol*_tol && k <= _iter)
	{
		// 预处理，使用A的L部分
		//向前解算，因为infoA是矩阵A的L部分，可重使用
		_cusparseError_id = cusparseDcsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, &floatone, descrL,
			d_ILUvals, d_iNonZeroRowSort, d_iColSort, infoA, d_r, d_y);
		//回溯子步
		_cusparseError_id = cusparseDcsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, &floatone, descrU,
			d_ILUvals, d_iNonZeroRowSort, d_iColSort, info_u, d_y, d_zm1);
		k++;
		if ( 1 == k )
		{
			cublasDcopy(cublasHandle, m_RowNum, d_zm1, 1, d_p, 1);
		}
		else
		{//步骤(2)
			cublasDdot(cublasHandle, m_RowNum, d_r, 1, d_zm1, 1, &numerator);
			cublasDdot(cublasHandle, m_RowNum, d_rm2, 1, d_zm2, 1, &denominator);
			beta = numerator/denominator;
			cublasDscal(cublasHandle, m_RowNum, &beta, d_p, 1);
			cublasDaxpy(cublasHandle, m_RowNum, &floatone, d_zm1, 1, d_p, 1) ;
		}
		/*		y = floatone * op(U) * p  + floatzero * y	*/
		_cusparseError_id = cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, 
			nzILU0, &floatone, descrU, 
			d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, d_p, &floatzero, d_omega);
		cublasDdot(cublasHandle, m_RowNum, d_r,	1, d_zm1,		1, &numerator);//步骤(3)
		cublasDdot(cublasHandle, m_RowNum, d_p,1, d_omega,	1, &denominator);//步骤(4)
		alpha = numerator / denominator;//步骤(4)
		cublasDaxpy(cublasHandle, m_RowNum, &alpha, d_p, 1, d_x, 1);//步骤(5)
		cublasDcopy(cublasHandle, m_RowNum, d_r, 1, d_rm2, 1);
		cublasDcopy(cublasHandle, m_RowNum, d_zm1, 1, d_zm2, 1);
		nalpha = -alpha;
		cublasDaxpy(cublasHandle, m_RowNum, &nalpha, d_omega, 1, d_r, 1);//步骤(6)
		cublasDdot(cublasHandle,  m_RowNum, d_r, 1, d_r, 1, &r1);
	}
	cudaMemcpy(h_x, d_x, m_RowNum*sizeof(double), cudaMemcpyDeviceToHost);
	//记录共轭梯度法的结束时间
	checkCudaErrors(cudaEventRecord(Event_Stop, NULL));
	// 等待所有的停止事件完成
	checkCudaErrors(cudaEventSynchronize(Event_Stop));
	checkCudaErrors(cudaEventElapsedTime(&msecUsed, Event_Start, Event_Stop));

	printf("内存拷贝时长Time= %.3f msec \n", 2*CopysecUsed);
	printf("所用总时长Time= %.3f msec \n", msecUsed);
	printf("  iteration = %3d, residual = %e \n", k, sqrt(r1));
	/*结果检验*/
	m_qaerr1 = CheckResult(m_nnzA, m_nA, 
		h_NoneZeroVal, h_iColSort, h_iNonZeroRowSort, 
		h_x, h_rhs);
	printf("  Convergence Test: %s \n", (k <= MAX_GPUITER) ? "OK" : "FAIL");
	m_iters = k;
	return (double)msecUsed;
}