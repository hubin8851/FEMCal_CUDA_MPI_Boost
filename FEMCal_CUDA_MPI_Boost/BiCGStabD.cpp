#include "StdAfx.h"
#include "HbxGloFunc.h"
#include "BiCGStabD.h"

extern double HBXDef::GetTimeStamp();

BiCGStabD::BiCGStabD( size_t _RowNum )
{
	_trans = CUSPARSE_OPERATION_NON_TRANSPOSE;
}


BiCGStabD::~BiCGStabD(void)
{
}

void	BiCGStabD::HostMalloc()
{
	h_ptrF  = (double *)malloc ( m_RowNum* sizeof(h_ptrF[0]));
	h_ptrR  = (double *)malloc ( m_RowNum* sizeof(h_ptrR[0]));
	h_ptrRW = (double *)malloc ( m_RowNum* sizeof(h_ptrRW[0]));
	h_ptrP  = (double *)malloc ( m_RowNum* sizeof(h_ptrP[0]));
	h_ptrPW = (double *)malloc ( m_RowNum* sizeof(h_ptrPW[0]));
	h_ptrS  = (double *)malloc ( m_RowNum* sizeof(h_ptrS[0]));
	h_ptrT  = (double *)malloc ( m_RowNum* sizeof(h_ptrT[0]));
	h_ptrV  = (double *)malloc ( m_RowNum* sizeof(h_ptrV[0]));
	h_ptrTX = (double *)malloc ( m_RowNum* sizeof(h_ptrTX[0]));
	h_ptrMval=(double *)malloc ( m_nnzA	 * sizeof(h_ptrMval[0]));
	CBaseDConjugate::HostMalloc();

	memset(h_ptrF,  0, m_RowNum * sizeof(h_ptrF[0]));
	memset(h_ptrR,  0, m_RowNum * sizeof(h_ptrR[0]));
	memset(h_ptrRW, 0, m_RowNum * sizeof(h_ptrRW[0]));
	memset(h_ptrP,  1, m_RowNum * sizeof(h_ptrP[0]));
	memset(h_ptrPW, 0, m_RowNum * sizeof(h_ptrPW[0]));
	memset(h_ptrS,  0, m_RowNum * sizeof(h_ptrS[0]));
	memset(h_ptrT,  0, m_RowNum * sizeof(h_ptrT[0]));
	memset(h_ptrV,  0, m_RowNum * sizeof(h_ptrV[0]));
	memset(h_ptrTX, 0, m_RowNum * sizeof(h_ptrTX[0]));

}

//重载devicemalloc，因为多了临时数组
void	BiCGStabD::DeviceMalloc()
{
	/* allocate device memory for csr matrix and vectors */
	HBXDef::CheckUserDefErrors( cudaMalloc ((void**)&d_ptrF, sizeof(d_ptrF[0]) * m_RowNum) );
	HBXDef::CheckUserDefErrors( cudaMalloc ((void**)&d_ptrR, sizeof(d_ptrR[0]) * m_RowNum) );
	HBXDef::CheckUserDefErrors( cudaMalloc ((void**)&d_ptrRW,sizeof(d_ptrRW[0])* m_RowNum) );
	HBXDef::CheckUserDefErrors( cudaMalloc ((void**)&d_ptrP, sizeof(d_ptrP[0]) * m_RowNum) );
	HBXDef::CheckUserDefErrors( cudaMalloc ((void**)&d_ptrPW,sizeof(d_ptrPW[0])* m_RowNum) );
	HBXDef::CheckUserDefErrors( cudaMalloc ((void**)&d_ptrS, sizeof(d_ptrS[0]) * m_RowNum) );
	HBXDef::CheckUserDefErrors( cudaMalloc ((void**)&d_ptrT, sizeof(d_ptrT[0]) * m_RowNum) );
	HBXDef::CheckUserDefErrors( cudaMalloc ((void**)&d_ptrV, sizeof(d_ptrV[0]) * m_RowNum) );
	HBXDef::CheckUserDefErrors( cudaMalloc ((void**)&d_ptrMval, sizeof(d_ptrMval[0]) * m_nnzA) );

	CBaseDConjugate::DeviceMalloc();
	cudaDeviceSynchronize();
	m_bCudaFree = false;
}

HBXDef::DataAloc_t	BiCGStabD::MemCpy( HBXDef::CopyType_t _temp )
{
	using namespace HBXDef;
	double start_matrix_copy = GetTimeStamp();
	
	if ( HBXDef::HostToDevice == _temp )
	{
		HBXDef::CheckUserDefErrors(cudaMemcpy (d_ptrMval, d_NonZeroVal  , (size_t)(m_nnzA * sizeof(d_ptrMval[0])), cudaMemcpyDeviceToDevice));//此算法多传入一个长度为nnz的数组
		//传入额外的中间变量,以下步骤是否可以省略？？？
		HBXDef::CheckUserDefErrors(cudaMemcpy (d_ptrF, h_ptrF, (size_t)(m_RowNum * sizeof(d_ptrF[0])), cudaMemcpyHostToDevice));
		HBXDef::CheckUserDefErrors(cudaMemcpy (d_ptrR, h_ptrR, (size_t)(m_RowNum * sizeof(d_ptrR[0])), cudaMemcpyHostToDevice));
		HBXDef::CheckUserDefErrors(cudaMemcpy (d_ptrRW,h_ptrRW,(size_t)(m_RowNum * sizeof(d_ptrRW[0])),cudaMemcpyHostToDevice));
		HBXDef::CheckUserDefErrors(cudaMemcpy (d_ptrP, h_ptrP, (size_t)(m_RowNum * sizeof(d_ptrP[0])), cudaMemcpyHostToDevice));
		HBXDef::CheckUserDefErrors(cudaMemcpy (d_ptrPW,h_ptrPW,(size_t)(m_RowNum * sizeof(d_ptrPW[0])),cudaMemcpyHostToDevice));
		HBXDef::CheckUserDefErrors(cudaMemcpy (d_ptrS, h_ptrS, (size_t)(m_RowNum * sizeof(d_ptrS[0])), cudaMemcpyHostToDevice));
		HBXDef::CheckUserDefErrors(cudaMemcpy (d_ptrT, h_ptrT, (size_t)(m_RowNum * sizeof(d_ptrT[0])), cudaMemcpyHostToDevice));
		HBXDef::CheckUserDefErrors(cudaMemcpy (d_ptrV, h_ptrV, (size_t)(m_RowNum * sizeof(d_ptrV[0])), cudaMemcpyHostToDevice));
	}							

	double stop_matrix_copy = GetTimeStamp();
	std::cout<<"双共轭梯度法内存拷贝耗时共计:"<<stop_matrix_copy-start_matrix_copy<<std::endl;
	return CBaseDConjugate::MemCpy(_temp);
}

//初始化描述器等对象
bool	BiCGStabD::InitialDescr( cudaStream_t _stream, cusparseMatrixType_t _MatrixType )
{
	m_stream = _stream;
	CBaseDConjugate::InitialDescr( _stream, _MatrixType );

	/* create three matrix descriptors */
	cusparseStatus_t	status1 = cusparseCreateMatDescr(&descra); 
	cusparseStatus_t	status2 = cusparseCreateMatDescr(&descrm); 
	if ((status1 != CUSPARSE_STATUS_SUCCESS) || (status2 != CUSPARSE_STATUS_SUCCESS))
	{
			fprintf( stderr, "!!!! CUSPARSE cusparseCreateMatDescr (coefficient matrix or preconditioner) error\n" );
			return EXIT_FAILURE;
	}

	//设定矩阵A和M的相关属性，比如转置与否，及首索引号
	HBXDef::CheckUserDefErrors( cusparseSetMatType(descra,CUSPARSE_MATRIX_TYPE_GENERAL) );
	if ( m_CSRIndexBase == CUSPARSE_INDEX_BASE_ZERO )
	{
		HBXDef::CheckUserDefErrors( cusparseSetMatIndexBase(descra,CUSPARSE_INDEX_BASE_ZERO) );
	}
	else
	{
		HBXDef::CheckUserDefErrors( cusparseSetMatIndexBase(descra,CUSPARSE_INDEX_BASE_ONE) );
	}

	HBXDef::CheckUserDefErrors( cusparseSetMatType(descrm,CUSPARSE_MATRIX_TYPE_GENERAL) );
	if ( m_CSRIndexBase == CUSPARSE_INDEX_BASE_ZERO )
	{
		HBXDef::CheckUserDefErrors( cusparseSetMatIndexBase(descrm,CUSPARSE_INDEX_BASE_ZERO) );
	}
	else
	{
		HBXDef::CheckUserDefErrors( cusparseSetMatIndexBase(descrm,CUSPARSE_INDEX_BASE_ONE) );
	}


	return true;
}

//@maxit:最大迭代次数
//@tol:许用误差精度
//@ttt_sv:传入的预处理所用时间
void BiCGStabD::gpu_pbicgstab( int maxit, double tol, double ttt_sv )
{
	using namespace HBXDef;
	double rho, rhop, beta, alpha, negalpha, omega, negomega, temp, temp2;
	double nrmr, nrmr0;
	rho = 0.0;
	double zero = 0.0;
	double one  = 1.0;
	double mone = -1.0;
	int i=0;
	int j=0;
	double ttl,ttl2,ttu,ttu2,ttm,ttm2;
	double ttt_mv=0.0;

	//在给定初值条件下计算r0=b-Ax0的残差
	checkCudaErrors(cudaThreadSynchronize());
	ttm = GetTimeStamp();
	//dr = A*x0
	checkCudaErrors(cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, m_RowNum, m_nnzA, &one, 
									descra, d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, d_x, &zero, d_ptrR));

	cudaThreadSynchronize();
	ttm2= GetTimeStamp();
	ttt_mv += (ttm2-ttm);
	//printf("matvec %f (s)\n",ttm2-ttm);

	//d_r = -d_r
	checkCudaErrors(cublasDscal(cublasHandle, m_RowNum, &mone, d_ptrR, 1));
	//d_r= d_f + d_r
	checkCudaErrors(cublasDaxpy(cublasHandle, m_RowNum, &one, d_ptrF, 1, d_ptrR, 1));

	//将r拷贝至r^{\hat} 和 p并计算残差
	checkCudaErrors(cublasDcopy(cublasHandle, m_RowNum, d_ptrR, 1, d_ptrRW, 1));
	checkCudaErrors(cublasDcopy(cublasHandle, m_RowNum, d_ptrR, 1, d_ptrP, 1)); 
	checkCudaErrors(cublasDnrm2(cublasHandle, m_RowNum, d_ptrR, 1, &nrmr0));

	for (i=0; i<maxit; )
	{
		rhop = rho;
		checkCudaErrors(cublasDdot(cublasHandle, m_RowNum, d_ptrRW, 1, d_ptrR, 1, &rho));	//d_rw .d_r
		if (i>0)	//即初次迭代后
		{
			beta= (rho/rhop) * (alpha/omega);
			negomega = -omega;
			checkCudaErrors(cublasDaxpy(cublasHandle,m_RowNum, &negomega, d_ptrV, 1, d_ptrP, 1));	//d_p = - omega * d_v + d_p
			checkCudaErrors(cublasDscal(cublasHandle,m_RowNum, &beta, d_ptrP, 1));					//d_p *= beta;
			checkCudaErrors(cublasDaxpy(cublasHandle,m_RowNum, &one, d_ptrR, 1, d_ptrP, 1));		//d_p = d_r+ d_p
		}
		//预处理阶段 (lower and upper triangular solve)
		checkCudaErrors(cudaThreadSynchronize());
		ttl  = GetTimeStamp();
		checkCudaErrors(cusparseSetMatFillMode(descrm,CUSPARSE_FILL_MODE_LOWER));
		checkCudaErrors(cusparseSetMatDiagType(descrm,CUSPARSE_DIAG_TYPE_UNIT));
		checkCudaErrors(cusparseDcsrsv_solve( cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
												m_RowNum,&one,descrm,
												d_ptrMval, d_iNonZeroRowSort, d_iColSort, info_l, d_ptrP, d_ptrT));
		checkCudaErrors(cudaThreadSynchronize());
		ttl2 = GetTimeStamp();
		ttu  = GetTimeStamp();

		checkCudaErrors(cusparseSetMatFillMode(descrm,CUSPARSE_FILL_MODE_UPPER));
		checkCudaErrors(cusparseSetMatDiagType(descrm,CUSPARSE_DIAG_TYPE_NON_UNIT));
		checkCudaErrors(cusparseDcsrsv_solve(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
												m_RowNum,&one,descrm,
												d_ptrMval, d_iNonZeroRowSort, d_iColSort, info_u, d_ptrT, d_ptrPW));
		checkCudaErrors(cudaThreadSynchronize());
		ttu2 = GetTimeStamp();
		ttt_sv += (ttl2-ttl)+(ttu2-ttu);

		checkCudaErrors(cudaThreadSynchronize());
		ttm = GetTimeStamp();
		// d_v = A * d_pw;
		checkCudaErrors(cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
												m_RowNum, m_RowNum, m_nnzA, &one, descra,
												d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, d_ptrPW, &zero, d_ptrV));
		checkCudaErrors(cudaThreadSynchronize());
		ttm2= GetTimeStamp();
		ttt_mv += (ttm2-ttm);
		//printf("matvec %f (s)\n",ttm2-ttm);

		checkCudaErrors(cublasDdot(cublasHandle, m_RowNum, d_ptrRW, 1, d_ptrV, 1,&temp));  // d_v *= d_rw
		alpha= rho / temp;
		negalpha = -(alpha);
		checkCudaErrors(cublasDaxpy(cublasHandle, m_RowNum, &negalpha, d_ptrV,  1, d_ptrR,	1));
		checkCudaErrors(cublasDaxpy(cublasHandle, m_RowNum, &alpha,    d_ptrPW, 1, d_x,		1));
		checkCudaErrors(cublasDnrm2(cublasHandle, m_RowNum, d_ptrR, 1, &nrmr));

		if (nrmr < tol*nrmr0)
		{
			j=5;
			break;
		}
		//预处理步骤
		checkCudaErrors(cudaThreadSynchronize());
		ttl  = GetTimeStamp();
		checkCudaErrors(cusparseSetMatFillMode(descrm,CUSPARSE_FILL_MODE_LOWER));
		checkCudaErrors(cusparseSetMatDiagType(descrm,CUSPARSE_DIAG_TYPE_UNIT));
		//d_t = M * d_r
		checkCudaErrors(cusparseDcsrsv_solve(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
											m_RowNum, &one,descrm,
											d_ptrMval, d_iNonZeroRowSort, d_iColSort, info_l, d_ptrR, d_ptrT));
		checkCudaErrors(cudaThreadSynchronize());
		ttl2 = GetTimeStamp();
		ttu  = GetTimeStamp();

		checkCudaErrors(cusparseSetMatFillMode(descrm,CUSPARSE_FILL_MODE_UPPER));
		checkCudaErrors(cusparseSetMatDiagType(descrm,CUSPARSE_DIAG_TYPE_NON_UNIT));
		//d_s = M * d_t
		checkCudaErrors(cusparseDcsrsv_solve(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
											m_RowNum, &one,descrm,
											d_ptrMval, d_iNonZeroRowSort, d_iColSort, info_u, d_ptrT, d_ptrS));
		checkCudaErrors(cudaThreadSynchronize());
		ttu2 = GetTimeStamp();
		ttt_sv += (ttl2-ttl)+(ttu2-ttu);
		//printf("solve lower %f (s), upper %f (s) \n",ttl2-ttl,ttu2-ttu);
		//matrix-vector multiplication
		checkCudaErrors(cudaThreadSynchronize());
		ttm = GetTimeStamp();
		//d_t = A * d_s;
		checkCudaErrors(cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 
										m_RowNum, m_RowNum, m_nnzA, &one, descra, 
										d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, d_ptrS, &zero, d_ptrT));
		checkCudaErrors(cudaThreadSynchronize());
		ttm2= GetTimeStamp();
		ttt_mv += (ttm2-ttm);
		//printf("matvec %f (s)\n",ttm2-ttm);

		checkCudaErrors(cublasDdot(cublasHandle, m_RowNum, d_ptrT, 1, d_ptrR, 1,&temp));
		checkCudaErrors(cublasDdot(cublasHandle, m_RowNum, d_ptrT, 1, d_ptrT, 1,&temp2));
		omega= temp / temp2;
		negomega = -(omega);
		checkCudaErrors(cublasDaxpy(cublasHandle, m_RowNum, &omega,		d_ptrS, 1, d_x,		1));
		checkCudaErrors(cublasDaxpy(cublasHandle, m_RowNum, &negomega,	d_ptrT, 1, d_ptrR,	1));

		checkCudaErrors(cublasDnrm2(cublasHandle, m_RowNum, d_ptrR, 1,&nrmr));

		if (nrmr < tol*nrmr0)
		{
			i++;
			j=0;
			break;
		}
		i++;

	}//end for

	printf("gpu total solve time %f (s), matvec time %f (s)\n",ttt_sv,ttt_mv);
}


//各种特征值算法求解,返回计算时间，ms为单位.
//@_tol:许用残差
//@_iter:许用迭代次数
double	BiCGStabD::ConjugateWithGPU( const double &_tol, const int &_iter )
{
	using namespace HBXDef;
	/* 创建对LU分解预处理的相关信息 */
	HBXDef::CheckUserDefErrors(cusparseCreateSolveAnalysisInfo(&info_l));
	HBXDef::CheckUserDefErrors(cusparseCreateSolveAnalysisInfo(&info_u));

	double ttl = HBXDef::GetTimeStamp();
	//对矩阵M进行下三角和上三角的分解
	HBXDef::CheckUserDefErrors(cusparseSetMatFillMode(descrm, CUSPARSE_FILL_MODE_LOWER));
	HBXDef::CheckUserDefErrors(cusparseSetMatDiagType(descrm, CUSPARSE_DIAG_TYPE_UNIT));
	HBXDef::CheckUserDefErrors(cusparseDcsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		m_RowNum, m_nnzA, descrm,
		d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, info_l));

	HBXDef::CheckUserDefErrors(cudaThreadSynchronize());
	double ttl2 = HBXDef::GetTimeStamp();

	double ttu = HBXDef::GetTimeStamp();

	HBXDef::CheckUserDefErrors(cusparseSetMatFillMode(descrm, CUSPARSE_FILL_MODE_UPPER));
	HBXDef::CheckUserDefErrors(cusparseSetMatDiagType(descrm, CUSPARSE_DIAG_TYPE_NON_UNIT));
	HBXDef::CheckUserDefErrors(cusparseDcsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		m_RowNum, m_nnzA, descrm,
		d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, info_u));

	HBXDef::CheckUserDefErrors(cudaThreadSynchronize());
	double ttu2= HBXDef::GetTimeStamp();
	double ttt_sv = 0.0;
	ttt_sv += (ttl2-ttl)+(ttu2-ttu);

	printf("下三角分解耗时 %f (s), 上三角分解耗时 %f (s) \n",ttl2-ttl,ttu2-ttu);

	/* 在GPU上用不完全的LU分解完成矩阵M的下三角和上三角分解 */
	double start_ilu, stop_ilu;
	printf("CUSPARSE csrilu0 ");  

	start_ilu= HBXDef::GetTimeStamp();

	HBXDef::CheckUserDefErrors(cusparseDcsrilu0(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
		m_RowNum,descra,
		d_ptrMval,d_iNonZeroRowSort,d_iColSort,info_l));

#ifdef TIME_INDIVIDUAL_LIBRARY_CALLS
	HBXDef::CheckUserDefErrors(cudaThreadSynchronize());
	stop_ilu = HBXDef::GetTimeStamp();
#endif // TIME_INDIVIDUAL_LIBRARY_CALLS

	fprintf (stdout, "time(s) = %10.8f \n",stop_ilu-start_ilu);

	double calcstart = HBXDef::GetTimeStamp();
	gpu_pbicgstab( m_RowNum, _tol, _iter );
	HBXDef::CheckUserDefErrors(cudaThreadSynchronize());
	double calcstop = HBXDef::GetTimeStamp();

	/* copy the result into host memory */
	checkCudaErrors( cudaMemcpy (h_ptrTX, d_x, (size_t)(m_RowNum * sizeof(h_ptrTX[0])), cudaMemcpyDeviceToHost) );

	return calcstop - calcstart;
}