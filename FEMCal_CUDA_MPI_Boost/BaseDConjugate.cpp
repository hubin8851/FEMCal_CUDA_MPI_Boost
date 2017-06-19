#include "StdAfx.h"
#include "HbxGloFunc.h"
#include "BaseDConjugate.h"


CBaseDConjugate::CBaseDConjugate( size_t _RowNum ):m_RowNum((int)_RowNum)
{
	m_nA = (int)m_RowNum+1;
	m_bGenTridiag = false;
	m_bCudaFree = true;
	m_bGPUFlag = false;
	m_DevID = -1;
	m_DataAloc = HBXDef::INVALID;
	m_bSave = false;				//�Ƿ����

	//��������ʼ��
	m_stream = 0;
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


CBaseDConjugate::~CBaseDConjugate(void)
{
	if (NULL == cusparseHandle)
	{
		return;
	}
	else FreeGPUResource();
	//�ͷ�CPU��Դ��todo something
}

void	CBaseDConjugate::HostMalloc()
{
	if (m_bGenTridiag) return;
	//�ڴ����
	if ( nullptr != h_NoneZeroVal )	//�жϽṹ���������Ƿ��Ѿ�����ֵ
	{
		delete []h_NoneZeroVal;
		h_NoneZeroVal = nullptr;
	}
	_cuError_id = cudaHostAlloc((void**)&h_NoneZeroVal, sizeof(double)*m_nnzA, 0);
	if ( nullptr != h_iColSort )	
	{
		delete []h_iColSort;
		h_iColSort = nullptr;
	}
	_cuError_id = cudaHostAlloc((void**)&h_iColSort, sizeof(int)*m_nnzA, 0);
	if ( nullptr != h_iNonZeroRowSort )	
	{
		delete []h_iNonZeroRowSort;
		h_iNonZeroRowSort = nullptr;
	}
	_cuError_id = cudaHostAlloc((void**)&h_iNonZeroRowSort, sizeof(int)*m_nA, 0);
	if ( nullptr != h_x )	
	{
		delete []h_x;
		h_x = nullptr;
	}
	_cuError_id = cudaHostAlloc((void**)&h_x, sizeof(double)*m_RowNum, 0);
	h_x = (double *) malloc(sizeof(double) * m_RowNum);
	memset( h_x, 0, m_RowNum );
	if ( nullptr != h_rhs )	
	{
		delete []h_rhs;
		h_rhs = nullptr;
	}
	_cuError_id = cudaHostAlloc((void**)&h_rhs, sizeof(double)*m_RowNum, 0);
	//	h_rhs = (double *) malloc(sizeof(double) * m_RowNum);	//��Ϊ������ҳ�ڴ�ĵ�ַ����
	memset( h_rhs, 0, m_RowNum );
}

void	CBaseDConjugate::DeviceMalloc()
{
#pragma region  ��ָ��ʱCSR��ʽ������Դ����
	//CSR����ֵ�Դ����
	_cuError_id = cudaMalloc((void **)&d_NonZeroVal, m_nnzA*sizeof(double));
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
	_cuError_id = cudaMalloc((void **)&d_r, m_RowNum*sizeof(double));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	//��������ֵ�Դ����
	_cuError_id = cudaMalloc((void **)&d_x, m_RowNum*sizeof(double));
	if (_cuError_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	m_bCudaFree = false;
#pragma  endregion
}

/* ������������ԽǶԳ��� */
void		CBaseDConjugate::genTridiag( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath )
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
	h_NoneZeroVal[0] = (double)rand()/RAND_MAX + 10.0f;
	h_NoneZeroVal[1] = (double)rand()/RAND_MAX;
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
		h_iColSort[start] = (int)i - 1;
		h_iColSort[start+1] = (int)i;

		if (i < m_RowNum-1)
		{
			h_iColSort[start+2] = (int)i + 1;
		}

		h_NoneZeroVal[start] =h_NoneZeroVal[start-1];//��һ�еĶԳ�ֵ
		h_NoneZeroVal[start+1] = (double)rand()/RAND_MAX + 10.0f;	//�Խ����ϵ�ֵ
		if (i < m_RowNum-1)
		{
			h_NoneZeroVal[start+2] = (double)rand()/RAND_MAX;
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
}


//��������ĵ�ʽ�Ҷ�����
void CBaseDConjugate::genRhs( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath )
{
	using namespace std;
	srand((unsigned)time(NULL));
	cout<<"���ɵ�ʽ�ұߣ�ԭʼ�������ڴ���Ĩȥ..."<<endl;
	m_SavePath = _savepath;
	h_vRhs.resize(_genRowNum);
	
	for (int i = 0; i < _genRowNum; i++)
	{
		h_rhs[i] = h_vRhs[i] = 10*(double)rand()/RAND_MAX+1.0f;
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

//��Mtx��ʽ���ļ���ȡ��Ԫ����ϡ�������Ϣ
HBXDef::DataAloc_t		CBaseDConjugate::SetStiffMatFromMTX( const char* _mat_filename, const boost::filesystem::path& _dir )
{
	boost::filesystem::path _tmppath(_dir);
	_tmppath.append(_mat_filename);
// 	char *TmpPath = new char[_tmppath.string().size()+1];
// 	memcpy( TmpPath, _tmppath.string().data(), _tmppath.string().size() );
// 	TmpPath[_tmppath.string().size()] = '\0';
	if ( HBXDef::loadMMSparseMat( _tmppath.string().c_str(), true, &m_RowNum, &m_ColNum, &m_nnzA, h_vNoneZeroVal, h_viNonZeroRowSort, h_viColSort) )
	{
		fprintf (stderr, "!!!! cusparseLoadMMSparseMatrix FAILED\n");
		return HBXDef::DataAloc_t::INVALID;
	}

	HostMalloc();
	if (0 == h_viColSort[0] && 0 == h_viNonZeroRowSort[0])
	{
		m_CSRIndexBase = CUSPARSE_INDEX_BASE_ZERO;
		std::cout<<"ѹ����ʽ�ľ�������Ϊ������ʼ..."<<std::endl;
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
	return HBXDef::DataAloc_t::DATAINMEM;
}

//���ⲿ��ȡ��������ָ�뼰�䳤�ȵõ��նȾ����CSR��ʽ
HBXDef::DataAloc_t	CBaseDConjugate::SetStiffMat( const void *const _srcVal, const void *const _srcCol, const void *const _srcRow,
												 size_t _nnA, size_t _nA, bool _bsave )
{
	h_vNoneZeroVal.clear();
	h_viColSort.clear();
	h_viNonZeroRowSort.clear();
	using namespace HBXDef;
	const double *_tmpVal = static_cast<const double*>(_srcVal);
	const size_t *_tmpCol = static_cast<const size_t*>(_srcCol);
	const size_t *_tmpRow = static_cast<const size_t*>(_srcRow);
	double *p_NoneZeroVal = const_cast<double*>(_tmpVal);
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

//���ļ��ж�ȡ�նȾ���
HBXDef::DataAloc_t	CBaseDConjugate::SetStiffMat(const char* _NoneZeroVal,const char* _ColSort,const char* _ValRowSort, const boost::filesystem::path& _dir)
{
	if (0 == HBXDef::_ImportCSRdataFromFile<double, int>(  h_vNoneZeroVal, h_viColSort, h_viNonZeroRowSort, _NoneZeroVal, _ColSort, _ValRowSort, _dir) )
	{
		std::cout<<"δ����ȷ��ȡCSR�ļ�����..."<<std::endl;
		return HBXDef::DataAloc_t::INVALID;
	}
	//�����������Բ飬�������ڲ������Ƿ����
	HBXDef::CSRFormatError_t _tmpErr =  HBXDef::verify_pattern( h_viNonZeroRowSort.size()-1, h_viColSort.size(), h_viNonZeroRowSort.data(), h_viColSort.data() );
	m_nnzA = (int)h_vNoneZeroVal.size();
	if (m_nnzA == h_vNoneZeroVal.size() && m_nnzA == h_viColSort.size() && m_nA == h_viNonZeroRowSort.size())//m_nA = N+1���ڴ�����
	{
		std::cout<<"��ȷ��ȡCSR����,ά����ȷ..."<<std::endl;
		HostMalloc();
		if (0 == h_viColSort[0] && 0 == h_viNonZeroRowSort[0])
		{
			m_CSRIndexBase = CUSPARSE_INDEX_BASE_ZERO;
			std::cout<<"ѹ����ʽ�ľ�������Ϊ������ʼ..."<<std::endl;
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


HBXDef::DataAloc_t	CBaseDConjugate::SetLoadVec(const char* _FileName,
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

//�����ⲿ����غ�����ָ��
HBXDef::DataAloc_t	CBaseDConjugate::SetLoadVec( const void* _LoadVec, size_t _RowNum )
{
	using namespace HBXDef;
	const double* _tmpload = static_cast<const double*>(_LoadVec);
	h_vRhs.resize(_RowNum);
	h_rhs = const_cast<double*>(_tmpload);
	m_nA = (int)_RowNum + 1;
	//�ڴ���ÿ��Ԫ�ظ�ֵ����Ϊvector���ڴ���ܲ�����
	for (int i = 0; i < _RowNum; i++)
	{
		h_vRhs[i] = h_rhs[i];
	}
	return DataAloc_t::DATAINMEM;
}

/************************************************************************/
/*                CUDA����ģ��                                          */
/************************************************************************/
//ѯ���Ƿ�CUDA����
bool	CBaseDConjugate::isDevicePossible()
{
	return m_bGPUFlag;
}

//��ȡGPU�豸���ԣ��������е��豸�����ں����ڻ�������豸��
size_t		CBaseDConjugate::GetGPUDeviceQuery()
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
			(double)deviceProp.totalGlobalMem/1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);
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
		printf("\n�����ڴ�:%u Byte\n", deviceProp.sharedMemPerBlock);
	}
	return deviceCount;
}


//��û�����ʲô����ֵ��Ҳ���������ڴ濽����ö�٣��Ժ���˵
HBXDef::DataAloc_t		CBaseDConjugate::MemCpy(HBXDef::CopyType_t _temp)
{
	if (HBXDef::DeviceToHost == _temp)
	{
		_cuError_id = cudaMemcpy( h_NoneZeroVal,			d_NonZeroVal,				m_nnzA*sizeof(double),		cudaMemcpyDeviceToHost);
		if (_cuError_id != cudaSuccess)//�����з�����
		{
			printf("cudaMemcpy returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( h_iColSort,				d_iColSort,					m_nnzA*sizeof(int),			cudaMemcpyDeviceToHost);
		if (_cuError_id != cudaSuccess)//�����з�����������
		{
			printf("cudaMemcpy returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( h_iNonZeroRowSort,	d_iNonZeroRowSort,	m_nA*sizeof(int),			cudaMemcpyDeviceToHost);
		if (_cuError_id != cudaSuccess)//ÿ���׸�������������
		{
			printf("cudaMemcpy returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( h_x,			d_x,			m_RowNum*sizeof(double), cudaMemcpyDeviceToHost);
		if (_cuError_id != cudaSuccess)//����������
		{
			printf("cudaMemcpy returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( h_rhs,			d_r,			m_RowNum*sizeof(double),cudaMemcpyDeviceToHost);
		if (_cuError_id != cudaSuccess)//��ʽ��ֵ
		{
			printf("cudaMemcpy returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		return HBXDef::DataAloc_t::DATAINMEM;
	}
	else if (HBXDef::HostToDevice == _temp)
	{
		double start_matrix_copy = HBXDef::GetTimeStamp();
		_cuError_id = cudaMemcpy( d_NonZeroVal,			h_NoneZeroVal,				m_nnzA*sizeof(double),		cudaMemcpyHostToDevice);
		if (_cuError_id != cudaSuccess)//����ֵ
		{
			printf("cudaMemcpy returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( d_iColSort,					h_iColSort,					m_nnzA*sizeof(int),		cudaMemcpyHostToDevice);
		if (_cuError_id != cudaSuccess)//������
		{
			printf("cudaMemcpy returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( d_iNonZeroRowSort,	h_iNonZeroRowSort,		m_nA*sizeof(int),			cudaMemcpyHostToDevice);
		if (_cuError_id != cudaSuccess)//������
		{
			printf("cudaMemcpy returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy( d_x,								h_x,			m_RowNum*sizeof(double),		cudaMemcpyHostToDevice);
		//		_cuError_id = cudaMemcpy( m_CSRInput.h_x,			d_x,			m_RowNum*sizeof(double), cudaMemcpyDeviceToHost);	//�����ã��޴���
		if (_cuError_id != cudaSuccess)//����������
		{
			printf("cudaMemcpy returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		_cuError_id = cudaMemcpy(d_r,								h_rhs,		m_RowNum*sizeof(double),		cudaMemcpyHostToDevice);
		if (_cuError_id != cudaSuccess)//��ʽ��ֵ
		{
			printf("cudaMemcpy returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
			exit(EXIT_FAILURE);
		}
		double end_matrix_copy = HBXDef::GetTimeStamp();
		std::cout<<"˫���Ȼ�����ϵ�������Ҷ������ڴ濽����ʱ����:"<<end_matrix_copy-start_matrix_copy<<std::endl;
		return HBXDef::DataAloc_t::DATAINGPU;
	}
	else return HBXDef::DataAloc_t::INVALID;
}

//��ʼ���������ȶ���
bool		CBaseDConjugate::InitialDescr( cudaStream_t _stream, cusparseMatrixType_t _MatrixType )
{
	m_stream = _stream;
	/* ����CUBLAS���� */
	_cublasError_id = cublasCreate(&cublasHandle);
	if (_cublasError_id != CUBLAS_STATUS_SUCCESS)
	{
		printf("cublasCreate returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	int version; 
	cublasGetVersion(cublasHandle, &version);
	std::cout<<"��ǰcublas�汾Ϊ:"<<version<<std::endl;
	/* ����CUSPARSE���� */
	_cusparseError_id = cusparseCreate(&cusparseHandle);
	if (_cusparseError_id != CUSPARSE_STATUS_SUCCESS)
	{
		printf("cusparseCreate returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	//���cusparse��cublas�����趨��
	if (m_stream)
	{
		checkCudaErrors( cudaStreamCreate(&m_stream) );
		if ( cublasSetStream(cublasHandle, m_stream) != CUBLAS_STATUS_SUCCESS )
		{
			fprintf (stderr, "!!!! cannot set CUBLAS stream\n");
			return false;
		}
		if ( cusparseSetStream(cusparseHandle, m_stream) != CUSPARSE_STATUS_SUCCESS )
		{
			fprintf (stderr, "!!!! cannot set CUSPARSE stream\n");
			return false;
		}	
	}

	/* A�����������*/
	_cusparseError_id = cusparseCreateMatDescr(&Matdescr);
	if (_cusparseError_id != CUSPARSE_STATUS_SUCCESS)
	{
		printf("cusparseCreateMatDescr returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	/* �����������*/
	_cusparseError_id = cusparseSetMatType(Matdescr, _MatrixType);
	if (_cusparseError_id != CUSPARSE_STATUS_SUCCESS)
	{
		printf("cusparseSetMatType returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	_cusparseError_id = cusparseSetMatIndexBase(Matdescr, m_CSRIndexBase);
	if (_cusparseError_id != CUSPARSE_STATUS_SUCCESS)
	{
		printf("cusparseSetMatIndexBase returned %d\n-> %s\n", (int)_cuError_id, cudaGetErrorString(_cuError_id));
		exit(EXIT_FAILURE);
	}
	return true;
}

//�ͷ�GPU��Դ
void	CBaseDConjugate::FreeGPUResource()
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
}

//����x����ֵ���ļ�
int		CBaseDConjugate::ExportEigToFile(const char* _EigvalFile, boost::filesystem::path _path )
{
	HBXDef::CheckUserDefErrors( HBXDef::OutputPointData( h_x, m_RowNum, _EigvalFile, _path ) );
		std::cout<<"����ֵ�������"
			<<_EigvalFile<<"�ļ���"<<std::endl;
		return true;
}