#include "stdafx.h"
#include "HbxGloFunc.h"
#include "BaseSConjugate.h"



CBaseSConjugate::CBaseSConjugate( size_t _RowNum ):m_RowNum(_RowNum)
{
	m_nA = (int)m_RowNum+1;
	m_bGenTridiag = false;
	m_bCudaFree = true;
	m_bGPUFlag = false;
	m_DevID = -1;
	m_DataAloc = HBXDef::INVALID;
	m_bSave = false;				//�Ƿ����

	//blas,sparse��matdescr�ĳ�ʼ��
	_cuError_id = cudaSuccess;
	_cublasError_id = CUBLAS_STATUS_SUCCESS;
	_cusparseError_id = CUSPARSE_STATUS_SUCCESS;
	cublasHandle	= 0;
	cusparseHandle	= 0;
	Matdescr		= 0;
	m_CSRIndexBase = cusparseIndexBase_t::CUSPARSE_INDEX_BASE_ZERO;

	msecUsed = 0.0f;
	m_iters = 0;

	h_NoneZeroVal = nullptr;
	h_iColSort = nullptr;
	h_iNonZeroRowSort = nullptr;
	h_x = nullptr;
	h_rhs = nullptr;
}

CBaseSConjugate::~CBaseSConjugate()
{
	if (NULL == cusparseHandle)
	{
		return;
	}
	else FreeGPUResource();
	//�ͷ�CPU��Դ��todo something
}

void	CBaseSConjugate::HostMalloc()
{
	if (m_bGenTridiag) return;
	//�ڴ����
	if ( nullptr != h_NoneZeroVal )	//�жϽṹ���������Ƿ��Ѿ�����ֵ
	{
		delete []h_NoneZeroVal;
		h_NoneZeroVal = nullptr;
	}
	h_NoneZeroVal = (float *) malloc(sizeof(float) * m_nnzA);
	if ( nullptr != h_iColSort )	
	{
		delete []h_iColSort;
		h_iColSort = nullptr;
	}
	h_iColSort = (int *) malloc(sizeof(int) * m_nnzA);
	if ( nullptr != h_iNonZeroRowSort )	
	{
		delete []h_iNonZeroRowSort;
		h_iNonZeroRowSort = nullptr;
	}
	h_iNonZeroRowSort = (int *) malloc(sizeof(int) * m_nA);
	if ( nullptr != h_x )	
	{
		delete []h_x;
		h_x = nullptr;
	}
	_cuError_id = cudaHostAlloc((void**)&h_x, sizeof(float)*m_RowNum, 0);
	h_x = (float *) malloc(sizeof(float) * m_RowNum);
	memset( h_x, 0, m_RowNum );
// 	for (int i=0; i<m_RowNum; i++)
// 	{
// 		h_x[i] = 0;
// 	}
	if ( nullptr != h_rhs )	
	{
		delete []h_rhs;
		h_rhs = nullptr;
	}
	_cuError_id = cudaHostAlloc((void**)&h_rhs, sizeof(float)*m_RowNum, 0);
//	h_rhs = (float *) malloc(sizeof(float) * m_RowNum);
	memset( h_rhs, 0, m_RowNum );
}

void	CBaseSConjugate::DeviceMalloc()
{
#pragma region  ��ָ��ʱCSR��ʽ������Դ����
	//CSR����ֵ�Դ����
	_cuError_id = cudaMalloc((void **)&d_NonZeroVal, m_nnzA*sizeof(float));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	//CSR������ֵ�Դ����
	_cuError_id = cudaMalloc((void **)&d_iColSort, m_nnzA*sizeof(int));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	//CSRÿ�е�һ�������������Դ����
	_cuError_id = cudaMalloc((void **)&d_iNonZeroRowSort, m_nA*sizeof(int));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	//���̵�ʽ�ұ������Դ����
	_cuError_id = cudaMalloc((void **)&d_r, m_RowNum*sizeof(float));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	//��������ֵ�Դ����
	_cuError_id = cudaMalloc((void **)&d_x, m_RowNum*sizeof(float));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	m_bCudaFree = false;
#pragma  endregion

	//Ϊ����ָ��ṹ�帳�ڴ�
	//	CSRInput_ptr	m_CSRInput = boost::make_shared<CSRInput<float>>();
#pragma region ʹ��thrust
	// 	h_vfx.resize(m_RowNum);
	// //���Դ��ʼ��thrust������ָ��
	// 	d_pNonZeroVal = thrust::device_malloc<float>(m_nnzA);
	// 	d_piColSort = thrust::device_malloc<int>(m_nnzA);
	// 	d_piNonZeroRowSort = thrust::device_malloc<int>(m_nA);
	// 	d_pfx = thrust::device_malloc<float>(m_RowNum);
	//  	d_pr = thrust::device_malloc<float>(m_RowNum);
	// //��ʼ��ԭʼָ��, host�ϵ�ԭʼָ�뱣����δ��ʼ��
	// 	m_CSRInput.d_NoneZeroVal = thrust::raw_pointer_cast( d_pNonZeroVal );
	// 	m_CSRInput.d_iColSort = thrust::raw_pointer_cast( d_piColSort );
	// 	m_CSRInput.d_iNonZeroRowSort = thrust::raw_pointer_cast( d_piNonZeroRowSort );
	// 	d_fx = thrust::raw_pointer_cast( d_pfx );
	// 	d_r = thrust::raw_pointer_cast( d_pr );
#pragma endregion 
}


/* ������������ԽǶԳ��� */
void		CBaseSConjugate::genTridiag( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath )
{
	using namespace std;
	srand((unsigned)time(NULL));
	cout<<"�������ԽǾ���ԭʼ�������ڴ���Ĩȥ..."<<endl;
	//��Ϊǿ����Ϊ���ԽǾ��󣬹�m_nnA,m_nA,rownum����Ҫ���¸�ֵ
	m_RowNum = (int)_genRowNum;
	m_nnzA = (int)(m_RowNum-2)*3+4; 
	m_nA = (int)m_RowNum+1;

#pragma region ���·����ڴ棬�����������
	HostMalloc();
	//��ֵ����
	h_iNonZeroRowSort[0] = 0, h_iColSort[0] = 0, h_iColSort[1] = 1;
	h_NoneZeroVal[0] = (float)rand()/RAND_MAX + 10.0f;
	h_NoneZeroVal[1] = (float)rand()/RAND_MAX;
	int start = 0;
 	for (size_t i = 1; i < m_RowNum; i++)
 	{
 		if (i > 1)
 		{
 			h_iNonZeroRowSort[i] = h_iNonZeroRowSort[i-1]+3;
 		}
 		else
 		{
 			h_iNonZeroRowSort[1] = 2;
 		}
 
 		start = ((int)i-1)*3 + 2;//startΪ�����з������е�i�еĵ�һ����
 		h_iColSort[start] = i - 1;
 		h_iColSort[start+1] = i;
 
 		if (i < m_RowNum-1)
 		{
 			h_iColSort[start+2] = i + 1;
 		}
 
		h_NoneZeroVal[start] =h_NoneZeroVal[start-1];//��һ�еĶԳ�ֵ
		h_NoneZeroVal[start+1] = (float)rand()/RAND_MAX + 10.0f;	//�Խ����ϵ�ֵ
 		if (i < m_RowNum-1)
 		{
 			h_NoneZeroVal[start+2] = (float)rand()/RAND_MAX;
 		}
 	}
 	h_iNonZeroRowSort[m_RowNum] = (m_RowNum-2)*3+4;	//���һ����
//CSR��ʽ�����Ƿ�洢����ش��룬20161101��С�����=��=
	if ( _bsave )
	{
		using namespace boost::gregorian;
		using namespace boost::posix_time;

//		ptime	_pathtime = second_clock::local_time();
		if ("" == _savepath)
		{
			m_SavePath = "..\\data\\export\\";
		}
		else m_SavePath = _savepath;
//		m_SavePath.append( to_iso_string(_pathtime) );
		if (exists(m_SavePath))
		{
			if ( boost::filesystem::is_empty(m_SavePath) )
			{
				remove(m_SavePath);
			}
		}
		else
		{
			remove_all(m_SavePath);
		}
		assert(!exists(m_SavePath));
		create_directories(m_SavePath);
		HBXDef::_ExportCSRdataToFile( h_NoneZeroVal, h_iColSort, h_iNonZeroRowSort,
																	m_nnzA, m_RowNum, "NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data",
																	m_SavePath );
	}
	m_bGenTridiag = true;

#pragma endregion

#pragma region thrust�汾
// 	h_viNonZeroRowSort.resize(m_RowNum+1);
// 	h_vNonZeroVal.resize( 3*(m_RowNum-2)+4 );
// 	h_viColSort.resize( 3*(m_RowNum-2)+4 );
// 	h_viNonZeroRowSort[0] = 0, h_viColSort[0] = 0, h_viNonZeroRowSort[1] = 1;
// 	h_vNonZeroVal[0] = (float)rand()/RAND_MAX + 10.0f;
// 	h_vNonZeroVal[1] = (float)rand()/RAND_MAX;
// 	int start;
// 
// 	for (size_t i = 1; i < m_RowNum; i++)
// 	{
// 		if (i > 1)
// 		{
// 			h_viNonZeroRowSort[i] = h_viNonZeroRowSort[i-1]+3;
// 		}
// 		else
// 		{
// 			h_viNonZeroRowSort[1] = 2;
// 		}
// 
// 		start = (i-1)*3 + 2;
// 		h_viColSort[start] = i - 1;
// 		h_viColSort[start+1] = i;
// 
// 		if (i < m_RowNum-1)
// 		{
// 			h_viColSort[start+2] = i + 1;
// 		}
// 
// 		h_vNonZeroVal[start] = h_vNonZeroVal[start-1];
// 		h_vNonZeroVal[start+1] = (float)rand()/RAND_MAX + 10.0f;
// 
// 		if (i < m_RowNum-1)
// 		{
// 			h_vNonZeroVal[start+2] = (float)rand()/RAND_MAX;
// 		}
// 	}
// 	h_viNonZeroRowSort[m_RowNum] = (m_RowNum-2)*3+4;
#pragma endregion
}

//��������ĵ�ʽ�Ҷ�����
void CBaseSConjugate::genRhs( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath )
{
	using namespace std;
	srand((unsigned)time(NULL));
	cout<<"���ɵ�ʽ�ұߣ�ԭʼ�������ڴ���Ĩȥ..."<<endl;
	m_SavePath = _savepath;
	h_vRhs.resize(_genRowNum);
	for (int i = 0; i < _genRowNum; i++)
	{
		h_rhs[i] = h_vRhs[i] = 10*(float)rand()/RAND_MAX+1.0f;
		h_x[i] = 0.0;
	}
	if (_bsave)
	{
		using namespace boost::gregorian;
		using namespace boost::posix_time;

		if (exists(m_SavePath))
		{
			if ( boost::filesystem::is_empty(m_SavePath) )
			{
				remove(m_SavePath);
			}
		}
		else
		{
			remove_all(m_SavePath);
		}
		assert(!exists(m_SavePath));
		create_directories(m_SavePath);
		HBXDef::OutputPointData( h_rhs, m_RowNum, "LoadVec.data", m_SavePath );
	}
}

//���ⲿ��ȡ��������ָ�뼰�䳤�ȵõ��նȾ����CSR��ʽ
HBXDef::DataAloc_t	CBaseSConjugate::SetStiffMat( const void *const _srcVal, const void *const _srcCol, const void *const _srcRow,
								size_t _nnA, size_t _nA, bool _bsave )
{
	h_vNoneZeroVal.clear();
	h_viColSort.clear();
	h_viNonZeroRowSort.clear();
	using namespace HBXDef;
	const float *_tmpVal = static_cast<const float*>(_srcVal);
	const size_t *_tmpCol = static_cast<const size_t*>(_srcCol);
	const size_t *_tmpRow = static_cast<const size_t*>(_srcRow);
	float *p_NoneZeroVal = const_cast<float*>(_tmpVal);
	size_t *p_iColSort = const_cast<size_t*>(_tmpCol);
	size_t *p_iNonZeroRowSort = const_cast<size_t*>(_tmpRow);
	m_nnzA = (int)_nnA;
	m_nA = (int)_nA;
	m_RowNum = m_nA-1;
	h_vNoneZeroVal.resize(m_nnzA);
	h_viColSort.resize(m_nnzA);
	h_viNonZeroRowSort.resize(m_nA);
	HostMalloc();
	//�ڴ���ÿ��Ԫ�ظ�ֵ����Ϊvector���ڴ���ܲ�����
	for (int i = 0; i < m_nnzA; i++)
	{
		h_NoneZeroVal[i] = h_vNoneZeroVal[i] = p_NoneZeroVal[i];
		h_iColSort[i] = h_viColSort[i] = (int)p_iColSort[i];
	}
	for (int i = 0; i < m_nA; i++)
	{
		h_viNonZeroRowSort[i] = h_viNonZeroRowSort[i] = (int)p_iNonZeroRowSort[i];
	}
	//�����浵�����ڶ������ݵĳ���У��
	if ( _bsave )
	{
		using namespace boost::gregorian;
		using namespace boost::posix_time;
		//�жϴ洢·���Ƿ�Ϊ��
		if (m_SavePath.empty())
		{
			m_SavePath = "..\\data\\UFget\\";
		}
		//		m_SavePath.append( to_iso_string(_pathtime) );
		if (exists(m_SavePath))
		{
			if ( boost::filesystem::is_empty(m_SavePath) )
			{
				remove(m_SavePath);
			}
		}
		else
		{
			remove_all(m_SavePath);
		}
		assert(!exists(m_SavePath));
		create_directories(m_SavePath);
		HBXDef::_ExportCSRdataToFile( h_NoneZeroVal, h_iColSort, h_iNonZeroRowSort,
			m_nnzA, m_RowNum, "NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data",
			m_SavePath );
	}
	return DataAloc_t::DATAINMEM;
}

HBXDef::DataAloc_t	CBaseSConjugate::SetStiffMat(const char* _NoneZeroVal,const char* _ColSort,const char* _ValRowSort, const boost::filesystem::path& _dir)
{
	if (0 == HBXDef::_ImportCSRdataFromFile<float, int>(  h_vNoneZeroVal, h_viColSort, h_viNonZeroRowSort, _NoneZeroVal, _ColSort, _ValRowSort, _dir) )
	{
		std::cout<<"δ����ȷ��ȡCSR�ļ�����..."<<std::endl;
		return HBXDef::DataAloc_t::INVALID;
	}
	m_nnzA = h_vNoneZeroVal.size();
	if (m_nnzA == h_vNoneZeroVal.size() && m_nnzA == h_viColSort.size() && m_nA == h_viNonZeroRowSort.size())//m_nA = N+1���ڴ�����
	{
		std::cout<<"��ȷ��ȡCSR����,ά����ȷ..."<<std::endl;
		HostMalloc();
		if (0 == h_viColSort[0] && 0 == h_viNonZeroRowSort[0])
		{
			m_CSRIndexBase = CUSPARSE_INDEX_BASE_ZERO;
			std::cout<<"����Ϊ������ʼ..."<<std::endl;
		}
		else if (1 == h_viColSort[0] && 1 == h_viNonZeroRowSort[0])
		{
			m_CSRIndexBase = CUSPARSE_INDEX_BASE_ONE;
			std::cout<<"��һΪ������ʼ..."<<std::endl;
		}
		else return HBXDef::DataAloc_t::INVALID;
		//�ڴ���ÿ��Ԫ�ظ�ֵ����Ϊvector���ڴ���ܲ�����
		for (int i = 0; i < m_nnzA; i++)
		{
			h_NoneZeroVal[i] = h_vNoneZeroVal[i];
			h_iColSort[i] = h_viColSort[i];
		}
		for (int i = 0; i < m_nA; i++)
		{
			h_iNonZeroRowSort[i] = h_viNonZeroRowSort[i];
		}
		return	HBXDef::DataAloc_t::DATAINMEM;
	}
	else
	{
		std::cout<<"δ����ȷ��ȡCSR�ļ����ݣ�����ά�Ȳ���..."<<std::endl;
		return HBXDef::DataAloc_t::INVALID;
	}
}

HBXDef::DataAloc_t	CBaseSConjugate::SetLoadVec(const void* _LoadVec, size_t _RowNum, bool _bsave )
{
	using namespace HBXDef;
	const double* _tmpload = static_cast<const double*>(_LoadVec);
	h_vRhs.resize(_RowNum);
	h_vRhs.resize(_RowNum);
	double *p_rhs = const_cast<double*>(_tmpload);
	m_nA = (int)_RowNum + 1;
	//�ڴ���ÿ��Ԫ�ظ�ֵ����Ϊvector���ڴ���ܲ�����
	for (int i = 0; i < _RowNum; i++)
	{
		h_vRhs[i] = p_rhs[i];
	}
	if (_bsave)
	{
		using namespace boost::gregorian;
		using namespace boost::posix_time;

		if (exists(m_SavePath))
		{
			if ( boost::filesystem::is_empty(m_SavePath) )
			{
				remove(m_SavePath);
			}
		}
		else
		{
			remove_all(m_SavePath);
		}
		assert(!exists(m_SavePath));
		create_directories(m_SavePath);
		HBXDef::OutputPointData( h_rhs, m_RowNum, "LoadVec.data", m_SavePath );
	}
	return DataAloc_t::DATAINMEM;
}

HBXDef::DataAloc_t	CBaseSConjugate::SetLoadVec(const char* _FileName,
											const boost::filesystem::path& _dir)
{
	size_t _RowNum = HBXDef::ReadVecByRowFromFile( h_vRhs,	_FileName, _dir);
	if (0 == _RowNum)
	{
		std::cout<<"δ��ȡ�غ�����,�Ƿ���Ҫ����..."<<std::endl;
		char ch;
		there:
		std::cout<<"����Y/N��ʾ���ɻ������Ҷ�����..."<<std::endl;
		ch = getchar();
		if ( !isalpha(ch) )
		{
			std::cout<<"��������ĸ..."<<std::endl;
			goto there;
		}
		else
		{
			switch (ch)
			{
			case 'y':case 'Y':
				genRhs( g_MatRowNum, false );
				return HBXDef::DataAloc_t::DATAINMEM;
			case  'n':case 'N':
				goto end;
				break;
			default:
				std::cerr<<"�������,����"<<std::endl;
				goto end;
			}
		}
		end:
		return HBXDef::DataAloc_t::INVALID;
	}
	else if ( m_RowNum == _RowNum)
	{
		std::cout<<"�غ�����ά����ȷ..."<<std::endl;
		memcpy(h_rhs,	h_vRhs.data(),	m_RowNum*sizeof(float) );
		return HBXDef::DataAloc_t::DATAINMEM;
	}
	std::cout<<"�غ�����ά�ȴ���..."<<std::endl;
	return HBXDef::DataAloc_t::INVALID;
}


/************************************************************************/
/*                CUDA����ģ��                                          */
/************************************************************************/
//ѯ���Ƿ�CUDA����
bool	CBaseSConjugate::isDevicePossible()
{
	return m_bGPUFlag;
}
//��ȡGPU�豸���ԣ��������е��豸�����ں����ڻ�������豸��
size_t		CBaseSConjugate::GetGPUDeviceQuery()
{
	int deviceCount = 0;
	cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

	if (error_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
		printf("Result = FAIL\n");
		return 0;
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
	checkCudaErrors(cudaGetDeviceProperties(&_BestDevice, 0));

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
		printf("\n�����ڴ�:%u Byte", deviceProp.sharedMemPerBlock);
	}
	return deviceCount;
}

//��û�����ʲô����ֵ��Ҳ���������ڴ濽����ö�٣��Ժ���˵
HBXDef::DataAloc_t		CBaseSConjugate::MemCpy(HBXDef::CopyType_t _temp)
{
	if (HBXDef::DeviceToHost == _temp)
	{
		_cuError_id = cudaMemcpy( h_NoneZeroVal,			d_NonZeroVal,				m_nnzA*sizeof(float),		cudaMemcpyDeviceToHost);
		if (_cuError_id != cudaSuccess)//�����з�����
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( h_iColSort,					d_iColSort,					m_nnzA*sizeof(int),			cudaMemcpyDeviceToHost);
		if (_cuError_id != cudaSuccess)//�����з�����������
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( h_iNonZeroRowSort,	d_iNonZeroRowSort,	m_nA*sizeof(int),			cudaMemcpyDeviceToHost);
		if (_cuError_id != cudaSuccess)//ÿ���׸�������������
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( h_x,			d_x,			m_RowNum*sizeof(float), cudaMemcpyDeviceToHost);
		if (_cuError_id != cudaSuccess)//����������
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( h_rhs,			d_r,			m_RowNum*sizeof(float),cudaMemcpyDeviceToHost);
		if (_cuError_id != cudaSuccess)//��ʽ��ֵ
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		return HBXDef::DataAloc_t::DATAINMEM;
	}
	else if (HBXDef::HostToDevice == _temp)
	{
		_cuError_id = cudaMemcpy( d_NonZeroVal,			h_NoneZeroVal,				m_nnzA*sizeof(float),		cudaMemcpyHostToDevice);
		if (_cuError_id != cudaSuccess)//����ֵ
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( d_iColSort,					h_iColSort,					m_nnzA*sizeof(int),		cudaMemcpyHostToDevice);
		if (_cuError_id != cudaSuccess)//������
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( d_iNonZeroRowSort,	h_iNonZeroRowSort,		m_nA*sizeof(int),			cudaMemcpyHostToDevice);
		if (_cuError_id != cudaSuccess)//������
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( d_x,								h_x,			m_RowNum*sizeof(float),		cudaMemcpyHostToDevice);
//		_cuError_id = cudaMemcpy( m_CSRInput.h_x,			d_x,			m_RowNum*sizeof(float), cudaMemcpyDeviceToHost);	//�����ã��޴���
		if (_cuError_id != cudaSuccess)//����������
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
 		}
		_cuError_id = cudaMemcpy(d_r,								h_rhs,		m_RowNum*sizeof(float),		cudaMemcpyHostToDevice);
		if (_cuError_id != cudaSuccess)//��ʽ��ֵ
		{
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
 		return HBXDef::DataAloc_t::DATAINGPU;
	}
	else return HBXDef::DataAloc_t::INVALID;
}

//��ʼ���������ȶ���
bool		CBaseSConjugate::InitialDescr( cusparseMatrixType_t _MatrixType )
{
	/* ����CUBLAS���� */
	_cublasError_id = cublasCreate(&cublasHandle);
	if (_cublasError_id != CUBLAS_STATUS_SUCCESS)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	int version; 
	cublasGetVersion(cublasHandle, &version);
	std::cout<<"��ǰcublas�汾Ϊ:"<<version<<std::endl;
	/* ����CUSPARSE���� */
	_cusparseError_id = cusparseCreate(&cusparseHandle);
	if (_cusparseError_id != CUSPARSE_STATUS_SUCCESS)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	/* A�����������*/
	_cusparseError_id = cusparseCreateMatDescr(&Matdescr);
	if (_cusparseError_id != CUSPARSE_STATUS_SUCCESS)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	/* �����������*/
	_cusparseError_id = cusparseSetMatType(Matdescr, _MatrixType);
	if (_cusparseError_id != CUSPARSE_STATUS_SUCCESS)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	_cusparseError_id = cusparseSetMatIndexBase(Matdescr, m_CSRIndexBase);
	if (_cusparseError_id != CUSPARSE_STATUS_SUCCESS)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
 	return true;
}

//�ͷ�GPU��Դ
void	CBaseSConjugate::FreeGPUResource()
{
#pragma region ��ָ���
	if ( false == m_bCudaFree )
	{
		/* �ͷ��豸�ڴ� */
		_cusparseError_id = cusparseDestroy(cusparseHandle);
		_cublasError_id = cublasDestroy(cublasHandle);
		_cuError_id = cudaFree(d_NonZeroVal);
		_cuError_id = cudaFree(d_iColSort);
		_cuError_id = cudaFree(d_iNonZeroRowSort);
		_cuError_id = cudaFree(d_x);
		_cuError_id = cudaFree(d_r);

		_cuError_id = cudaDeviceReset();	//����GPU�豸
		m_bCudaFree = true;
	}
	else return;
#pragma endregion

#pragma region thrust��
	//���ͷ�GPU��Դ��ͬʱresize�ڴ�������
// 	h_vNonZeroVal.clear();
// 	h_viColSort.clear();
// 	h_viNonZeroRowSort.clear();
// 	h_vfx.clear();
// 	h_vRhs.clear();
// 
// 	cusparseDestroy(cusparseHandle);
// 	cublasDestroy(cublasHandle);
// 	// �ڴ���Դ�����ͷţ���Ϊ�õ���vector
// 	cudaFree(m_CSRInput.d_NoneZeroVal);
// 	cudaFree(m_CSRInput.d_iColSort);
// 	cudaFree(m_CSRInput.d_iNonZeroRowSort);
// 	cudaFree(d_fx);
// 	cudaFree(d_r);
// 	cudaDeviceReset();	//����GPU�豸
#pragma endregion
}

//����x����ֵ���ļ�
int		CBaseSConjugate::ExportEigToFile(const char* _EigvalFile, boost::filesystem::path _path )
{
	HBXDef::CheckUserDefErrors (HBXDef::OutputPointData( h_x, m_RowNum, _EigvalFile, _path ) );
	std::cout<<"����ֵ�������"
			<<_EigvalFile<<"�ļ���"<<std::endl;
		return true;
}

void CBaseSConjugate::runtest()
{
	int M = 0, N = 0, nz = 0, *I = NULL, *J = NULL;
	float *val = NULL;
	const float tol = 1e-5f;
	const int max_iter = 10000;
	float *x;
	float *rhs;
	float a, b, na, r0, r1;
	int *d_col, *d_row;
	float *d_val, *d_x, dot;
	float *d_r, *d_p, *d_Ax;
	int k;
	float alpha, beta, alpham1;

	// This will pick the best possible CUDA capable device
	cudaDeviceProp deviceProp;
	int devID = 0;

	checkCudaErrors(cudaGetDeviceProperties(&deviceProp, devID));

	// Statistics about the GPU device
	printf("> GPU device has %d Multi-Processors, SM %d.%d compute capabilities\n\n",
		deviceProp.multiProcessorCount, deviceProp.major, deviceProp.minor);

	int version = (deviceProp.major * 0x10 + deviceProp.minor);

	if (version < 0x11)
	{
		exit(EXIT_SUCCESS);
	}

	/* Generate a random tridiagonal symmetric matrix in CSR format */
	M = N = 4;
	nz = (N-2)*3 + 4;
	I = (int *)malloc(sizeof(int)*(N+1));
	J = (int *)malloc(sizeof(int)*nz);
	val = (float *)malloc(sizeof(float)*nz);
	HBXDef::genTridiag( val, I, J, nz, N);

	x = (float *)malloc(sizeof(float)*N);
	rhs = (float *)malloc(sizeof(float)*N);

	for (int i = 0; i < N; i++)
	{
		rhs[i] = 1.0;
		x[i] = 0.0;
	}

	/* Get handle to the CUBLAS context */
	cublasHandle_t cublasHandle = 0;
	cublasStatus_t cublasStatus;
	cublasStatus = cublasCreate(&cublasHandle);

	checkCudaErrors(cublasStatus);

	/* Get handle to the CUSPARSE context */
	cusparseHandle_t cusparseHandle = 0;
	cusparseStatus_t cusparseStatus;
	cusparseStatus = cusparseCreate(&cusparseHandle);

	checkCudaErrors(cusparseStatus);

	cusparseMatDescr_t descr = 0;
	cusparseStatus = cusparseCreateMatDescr(&descr);

	checkCudaErrors(cusparseStatus);

	cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase( descr, m_CSRIndexBase );

	checkCudaErrors(cudaMalloc((void **)&d_col, nz*sizeof(int)));
	checkCudaErrors(cudaMalloc((void **)&d_row, (N+1)*sizeof(int)));
	checkCudaErrors(cudaMalloc((void **)&d_val, nz*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **)&d_x, N*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **)&d_r, N*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **)&d_p, N*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **)&d_Ax, N*sizeof(float)));
//	float* _tmp;	//������
//	_tmp = (float*)malloc( N*sizeof(float));//������

	cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, rhs, N*sizeof(float), cudaMemcpyHostToDevice);

	alpha = 1.0;
	alpham1 = -1.0;
	beta = 0.0;
	r0 = 0.;

	cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_x, &beta, d_Ax);
//	cudaMemcpy( _tmp,	d_Ax,	N*sizeof(float),		cudaMemcpyDeviceToHost);//������

	cublasSaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
	cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
//	cudaMemcpy( _tmp,	d_r,	N*sizeof(float),		cudaMemcpyDeviceToHost);//������

	k = 1;

	// ������ʱ��CUDA�¼�
	float	msecUsed = 0.0f;//ʱ��ͳ��float�ͱ���
	cudaError_t	_cuError_id=cudaSuccess;
	cudaEvent_t	Event_Start;
	cudaEvent_t Event_Stop;
	_cuError_id = cudaEventCreate(&Event_Start);
	_cuError_id = cudaEventCreate(&Event_Stop);

	while (r1 > tol*tol && k <= max_iter)
	{
		if (k > 1)
		{
			b = r1 / r0;
			cublasStatus = cublasSscal(cublasHandle, N, &b, d_p, 1);
			cublasStatus = cublasSaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
		}
		else
		{
			cublasStatus = cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
		}

		cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_p, &beta, d_Ax);
//		cudaMemcpy( _tmp,	d_Ax,	N*sizeof(float),		cudaMemcpyDeviceToHost);//������
		cublasStatus = cublasSdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
		a = r1 / dot;

		cublasStatus = cublasSaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
		na = -a;
		cublasStatus = cublasSaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);
//		cudaMemcpy( _tmp,	d_r,	N*sizeof(float),		cudaMemcpyDeviceToHost);//������

		r0 = r1;
		cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
		cudaThreadSynchronize();
		printf("iteration = %3d, residual = %e\n", k, sqrt(r1));
		k++;
	}

	cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

	//��¼�����ݶȷ��Ľ���ʱ��
	cudaEventRecord(Event_Stop, 0);
	checkCudaErrors(cudaEventSynchronize(Event_Stop));// �ȴ����е�ֹͣ�¼����
	checkCudaErrors(cudaEventElapsedTime(&msecUsed, Event_Start, Event_Stop));

	float rsum, diff, err = 0.0;

	for (int i = 0; i < N; i++)
	{
		rsum = 0.0;

		for (int j = I[i]; j < I[i+1]; j++)
		{
			rsum += val[j]*x[J[j]];
		}

		diff = fabs(rsum - rhs[i]);

		if (diff > err)
		{
			err = diff;
		}
	}

	cusparseDestroy(cusparseHandle);
	cublasDestroy(cublasHandle);

	free(I);
	free(J);
	free(val);
	free(x);
	free(rhs);
	cudaFree(d_col);
	cudaFree(d_row);
	cudaFree(d_val);
	cudaFree(d_x);
	cudaFree(d_r);
	cudaFree(d_p);
	cudaFree(d_Ax);

	printf("Test Summary:  Error amount = %f\n", err);
	exit((k <= max_iter) ? 0 : 1);
};