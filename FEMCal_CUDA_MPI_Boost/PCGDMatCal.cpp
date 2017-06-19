#include "StdAfx.h"
#include "HbxGloFunc.h"
#include "PCGDMatCal.h"

//PCG��˫����������
PCGDMatCal::PCGDMatCal( size_t _RowNum ):CBaseDConjugate( _RowNum )
{
	//Ԥ�������ݶȷ����ֲ�����ʼ��
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
		/* ����LU�ֽ��ȹ�������Ϣ */
		checkCudaErrors(cusparseDestroySolveAnalysisInfo(infoA));
		checkCudaErrors(cusparseDestroySolveAnalysisInfo(info_u));
		/* �ͷ��豸�ڴ� */
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
	cout<<"����PCG�����þ���ԭʼ�������ڴ���Ĩȥ..."<<endl;
	//��Ϊǿ����Ϊ���ԽǾ��󣬹�m_nnA,m_nA,rownum����Ҫ���¸�ֵ
	m_RowNum = _genRowNum;
	m_nnzA = 5*m_RowNum-4*(int)sqrt((double)m_RowNum);
	m_nA = m_RowNum+1;
	//���·����ڴ�,������ڴ���h_x��h_rhs����0
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
	//CSR��ʽ�����Ƿ�洢����ش��룬20161103
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
	//���̵�ʽ�ұ������Դ����
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

//��ʼ���������ȶ���
bool	PCGDMatCal::InitialDescr( cudaStream_t _stream, cusparseMatrixType_t _MatrixType )
{
	m_stream = _stream;
	h_ILUvals = (double *) malloc(m_nnzA*sizeof(double));
	_cuError_id = cudaMalloc((void **)&d_ILUvals, m_nnzA*sizeof(double));
	_cuError_id = cudaMalloc((void **)&d_zm1, m_RowNum*sizeof(double));
	_cuError_id = cudaMalloc((void **)&d_zm2, m_RowNum*sizeof(double));
	_cuError_id = cudaMalloc((void **)&d_rm2, m_RowNum*sizeof(double));

	CBaseDConjugate::InitialDescr( _stream, _MatrixType );
	/* ��������A�ķ�����Ϣ */
	_cusparseError_id = cusparseCreateSolveAnalysisInfo(&infoA);
	if ( CUSPARSE_STATUS_SUCCESS != _cusparseError_id)
	{
		printf("cusparseCreateSolveAnalysisInfo returned %d\n-> %s\n", (int)_cusparseError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	// �����Ը÷�ת�þ����������,�൱�ڰ�
	_cusparseError_id = cusparseDcsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		m_RowNum, m_nnzA, Matdescr, 
		d_NonZeroVal, d_iNonZeroRowSort, d_iColSort, infoA);
	if ( CUSPARSE_STATUS_SUCCESS != _cusparseError_id)
	{
		printf("cusparseDcsrsv_analysis returned %d\n-> %s\n", (int)_cusparseError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}

	/* ��A������ΪILU0��ֵ����*/
	cudaMemcpy(d_ILUvals, d_NonZeroVal, m_nnzA*sizeof(double), cudaMemcpyDeviceToDevice);
	/* ����������CSR��ʽ�洢��A����Ĳ���ȫLU�ֽ��Ա��ں��ڼ��� */
	_cusparseError_id = cusparseDcsrilu0(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, 
		Matdescr, d_ILUvals, d_iNonZeroRowSort, d_iColSort, infoA);
	/* ������LU�ֽ�Ԥ����������Ϣ */
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

//GPU���ٵ�Ԥ�������ݶȷ�/* ʹ��ILU��Ԥ�������ݶȷ�.*/
//�㷨�μ����㷽������-����Ԫ�ṹ�������м��㣨����������ά̩��,P183�㷨5.4
//PCG����ֵ�㷨���,���ؼ���ʱ�䣬msΪ��λ
//@_tol:�в�
//@_iter:����������
double	PCGDMatCal::ConjugateWithGPU( const double &_tol, const int &_iter )
{
	printf("\nʹ��LU�ֽ��Ԥ�������ݶȷ���starting...: \n");
	// ������ʱ��CUDA�¼�
	cudaEvent_t	Event_Start, Event_PreStart, Event_MemStart;
	cudaEvent_t Event_Stop, Event_PreEnd, Event_MemEnd;
	checkCudaErrors(cudaEventCreate(&Event_Start));
	checkCudaErrors(cudaEventCreate(&Event_Stop));
	checkCudaErrors(cudaEventCreate(&Event_PreStart));
	checkCudaErrors(cudaEventCreate(&Event_PreEnd));
	checkCudaErrors(cudaEventCreate(&Event_MemStart));
	checkCudaErrors(cudaEventCreate(&Event_MemEnd));

	double r0, r1, alpha, beta;
	int k = 0;//��������
	double numerator, denominator, nalpha;
	const double floatone = 1.0;
	const double floatzero = 0.0;
	float	msecUsed = 0.0f;	//ʱ��ͳ��float�ͱ���
	float  CopysecUsed = 0.0f;	//�ڴ濽��ʱ��ͳ����
	m_qaerr1 = 0.0;//�����
	int	nzILU0 = 2*m_RowNum - 1;//�������ĸ�ֵ�������ȷ�����

	checkCudaErrors(cudaEventRecord(Event_Start, NULL));	//��ʼ��ʱ


	_cuError_id = cudaMemcpy(d_r,		h_rhs,	m_RowNum*sizeof(double), cudaMemcpyHostToDevice);
	_cuError_id = cudaMemcpy(d_x,		h_x,		m_RowNum*sizeof(double), cudaMemcpyHostToDevice);

	//��¼�ڴ濽���Ľ���ʱ��
	checkCudaErrors(cudaEventRecord(Event_MemEnd, NULL));
	// �ȴ����е�ֹͣ�¼����
	checkCudaErrors(cudaEventSynchronize(Event_MemEnd));
	checkCudaErrors(cudaEventElapsedTime(&CopysecUsed, Event_Start, Event_MemEnd));

	cublasDdot(cublasHandle, m_RowNum, d_r, 1, d_r, 1, &r1);
	//	cudaMemcpy(m_CSRInput.h_x, d_r, m_RowNum*sizeof(double), cudaMemcpyDeviceToHost);//�����õģ�540Mת560���ʹ��󣬹�������...

	while (r1>_tol*_tol && k <= _iter)
	{
		// Ԥ����ʹ��A��L����
		//��ǰ���㣬��ΪinfoA�Ǿ���A��L���֣�����ʹ��
		_cusparseError_id = cusparseDcsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, &floatone, descrL,
			d_ILUvals, d_iNonZeroRowSort, d_iColSort, infoA, d_r, d_y);
		//�����Ӳ�
		_cusparseError_id = cusparseDcsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_RowNum, &floatone, descrU,
			d_ILUvals, d_iNonZeroRowSort, d_iColSort, info_u, d_y, d_zm1);
		k++;
		if ( 1 == k )
		{
			cublasDcopy(cublasHandle, m_RowNum, d_zm1, 1, d_p, 1);
		}
		else
		{//����(2)
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
		cublasDdot(cublasHandle, m_RowNum, d_r,	1, d_zm1,		1, &numerator);//����(3)
		cublasDdot(cublasHandle, m_RowNum, d_p,1, d_omega,	1, &denominator);//����(4)
		alpha = numerator / denominator;//����(4)
		cublasDaxpy(cublasHandle, m_RowNum, &alpha, d_p, 1, d_x, 1);//����(5)
		cublasDcopy(cublasHandle, m_RowNum, d_r, 1, d_rm2, 1);
		cublasDcopy(cublasHandle, m_RowNum, d_zm1, 1, d_zm2, 1);
		nalpha = -alpha;
		cublasDaxpy(cublasHandle, m_RowNum, &nalpha, d_omega, 1, d_r, 1);//����(6)
		cublasDdot(cublasHandle,  m_RowNum, d_r, 1, d_r, 1, &r1);
	}
	cudaMemcpy(h_x, d_x, m_RowNum*sizeof(double), cudaMemcpyDeviceToHost);
	//��¼�����ݶȷ��Ľ���ʱ��
	checkCudaErrors(cudaEventRecord(Event_Stop, NULL));
	// �ȴ����е�ֹͣ�¼����
	checkCudaErrors(cudaEventSynchronize(Event_Stop));
	checkCudaErrors(cudaEventElapsedTime(&msecUsed, Event_Start, Event_Stop));

	printf("�ڴ濽��ʱ��Time= %.3f msec \n", 2*CopysecUsed);
	printf("������ʱ��Time= %.3f msec \n", msecUsed);
	printf("  iteration = %3d, residual = %e \n", k, sqrt(r1));
	/*�������*/
	m_qaerr1 = CheckResult(m_nnzA, m_nA, 
		h_NoneZeroVal, h_iColSort, h_iNonZeroRowSort, 
		h_x, h_rhs);
	printf("  Convergence Test: %s \n", (k <= MAX_GPUITER) ? "OK" : "FAIL");
	m_iters = k;
	return (double)msecUsed;
}