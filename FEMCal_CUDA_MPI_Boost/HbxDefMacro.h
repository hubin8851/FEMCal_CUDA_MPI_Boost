#ifndef _HBXDEFMACRO_H
#define	_HBXDEFMACRO_H

 /******
     ����NDEBUG��Ϊ�������(debug)�ڼ�ʹ�õģ�����������ϣ������
 ���а�(release)��NDEBUG��ʧȥ���á�Ϊ���ܹ�ʹassert()�����ڷ���
 ����Ҳ����ʹ�ã��������������������NDEBUG����˼Ϊ��
     ����Ѿ�����NDEBUG�������Assert(x)����xΪ��ʱ��ִ��throw��
 ����ϵͳ�Ѷ���ķ�ʽ��������Assert(x)����Ϊ����assert(x)����
 ������ ����x���в��ԣ�xΪ��ʱ��������ֹ���򲢴�ӡ������Ϣ��
 *****/
 #ifdef NDEBUG			
  #define Assert(x)	if (!x)	throw;	
 #else			//������ϵͳ����ĺ���assert()����
  #include <cassert>
  #define Assert(x)	assert(x);
 #endif			//NDEBUG

#ifndef CUDA_CALL(x)	//�ж�CUDA������к�������ֵ
#define CUDA_CALL(x) { const cudaError_t a=(x); if (a!=cudaSuccess){ printf("\nCUDA Error:%s(err_num=%d)\n", cudaGetErrorString(a), a); cudaDeviceReset(); Assert(0);} }
#endif


#ifndef DEVICE_RESET
#define	DEVICE_RESET cudaDeviceReset();
#endif

#ifndef FLOATERROR
#define FLOATERROR 1.0e-6f
#endif // !FLOATERROR

#ifndef DOUBLEERROR
#define DOUBLEERROR 1.0e-15
#endif // !DOUBLEERROR

#ifndef LONGDOUBLEERROR
#define LONGDOUBLEERROR 1.0e-30
#endif // !LONGDOUBLEERROR

#ifndef MAXPATH_LGT
#define MAXPATH_LGT 256
#endif

//profiling the code
#define TIME_INDIVIDUAL_LIBRARY_CALLS


namespace HBXDef
{
	typedef double UserDefFloat;

	#ifndef MAX_GPUITER
	#define MAX_GPUITER	10000
	#endif

	#ifndef RUNCLC
	#define RUNCLC	1000
	#endif

	#ifndef	ERRORTOLERATE
	#define ERRORTOLERATE	1e-6
	#endif

	#define GRIDSIZE 8
	#define BLOCKSIZE 512
	#define ALLTHREAD 4096

	#define MAX_SHARED 48*1024

	const int MAX_GPU_COUNT = 32;	//��������豸����ĿǰҲֻ֧�����32���豸

	//�������㾫��ö��
	typedef enum
	{
		_PRECS_FLOAT_ = 1,
		_PRECS_DOUBLE_ = 2
	} PrecsType_t;

	//���л���������ö��
	typedef enum 
	{
		ROWMAJOR = 1,
		COLMAJOR = 2
	} SortType_t;	

	//�����ʽ��ʼ����ö��
	typedef enum IndexBase
	{
		INDEX_BASE_ZERO = 0, 
		INDEX_BASE_ONE = 1,
		INDEX_INIVIAL
	} IndexBase_t;

	//���ݿ�������ö��
	typedef enum CopyType
	{
		CopyFailed = -1,
		HostToDevice = 0,
		DeviceToHost = 1,
		HostToHost = 2,
		DeviceToDevice = 3
	} CopyType_t;	

	//�����ڴ�/�Դ��־λö��
	typedef	enum ResetType
	{
		RESETCPUMEM = 0,
		RESETGPUMEM = 1
	} ResetType_t;	

	//�غ�����ö��
	typedef enum LoadType
	{
		TRACTION = 0,
		CONCENTRATE = 1
	} LoadType_t;	//

	typedef enum
	{
		NORMAL = 0,		//��ͨ�ڴ����
		UNIFIEDMEM = 1,	//ͳһ�ڴ�
		PAGELOCK,		//��ҳ�ڴ�
		ARRAY2D,		//2ά����
	} CudaMalloc_t;

	//�洢��ʽö��
	typedef	enum SaveType
	{
		CSR = 0,
		CSC,
		COO,
		DIA,	//�ԽǸ�ʽ
		ELL,	//���и�ʽ
		JAD		//�ݶԽǸ�ʽ
	} SaveType_t;

	typedef enum CSRFormatError
	{
		CHECKSUCCESS = 0,
		NNZCHECKFAILED = 1,
		BASEINDEXCHECKFAILED,
		CORRUPTEDROW,
		COLVSBASEINDEXCHECK,
		SORTINGCOLINDEXCHECK
	} CSRFormatError_t;	//CSR��ʽ�������

	//�������
	typedef enum MMError
	{
		MM_COULD_NOT_READ_FILE = 1,	
		MM_PREMATURE_EOF,		
		MM_NOT_MTX,				
		MM_NO_HEADER,			
		MM_UNSUPPORTED_TYPE,		
		MM_LINE_TOO_LONG,		
		MM_COULD_NOT_WRITE_FILE
	} MMError_t;

	typedef enum FactorError
	{
		FACTORSUCCESS = 0,
		RANKERROR = 1,
		NOTSYMETRY = 2,
		NOTREGULAR = 3,
		PRINCIPALELEMZERO = 4,
		SYMETRY	= 5,
		REGULAR = 6,
		SYMETRYbutNOREGULAR = 7,
		SYMETRYandREGULAR = 8
	} FactorError_t;	//��ʽ�ֽ��������

	typedef enum DataAloc
	{
		INVALID = -1,
		DATAINMEM = 0,
		DATAINGPU = 1
	} DataAloc_t;	//������������λ��

	typedef enum ControlMark
	{
		RESET = 0,
		NODE = 10,
		ELEMENT = 20,
		MATERIAL = 30,
		DENSITY = 32,
		ELASTIC = 34,
		SECTION = 40,
		ELSET = 51,
		NSET = 52,
		STIFFPARA = 101,
		MASSPARA = 102
	} ControlMark_t;


	typedef enum
	{
		PLANESTRESS = 0,
		PLANESTRAIN
	} SimpleType_t;

	//�����㷨ö��
	typedef enum SolveMethod
	{
		//ֱ�ӷ�
		QR_Mtd = 0,
		//������
		BICGSTAB_Mtd,
		CG_Mtd,
		DiagPre_Mtd,
		ILUT_Mtd,
		PCG_Mtd,
		LeastSquareDiagPre_Mtd,
		LeastSquareCG_Mtd

	} SolveMethod_t;

	typedef enum 
	{
		TITLE = 0,
		NAME,
		ID,
		SparseMat,
		AUX
	} MatField_t;

//��������ҳ�ڴ棬���ô˺�
#ifndef	CLEANPAGELOCK()
#define CLEANPAGELOCK()	\
	do {    \
		if ( h_NoneZeroVal )	cudaFreeHost(h_NoneZeroVal);		\
		if ( h_iColSort )		cudaFreeHost(h_iColSort);			\
		if ( h_iNonZeroRowSort) cudaFreeHost(h_iNonZeroRowSort);	\
		if ( h_x)				cudaFreeHost(h_x);					\
		if ( h_rhs)				cudaFreeHost(h_rhs);				\
	} while (0);
#endif

//�μ�http://blog.csdn.net/stpeace/article/details/9063293,�������buffer��������Ϣ
#ifndef CLEANDEVICE()
#define CLEANDEVICE()	\
	do {				\
		if (d_NoneZeroVal) cudaFree(d_NoneZeroVal);								\
		if (d_iColSort) cudaFree(d_iColSort);									\
		if (d_iNonZeroRowSort) cudaFree(d_iNonZeroRowSort);						\
		if (d_x) cudaFree(d_x);													\
		if (d_rhs) cudaFree(d_rhs);												\
		if (stream)           checkCudaErrors(cudaStreamDestroy(stream));       \
		if (cublasHandle)     checkCudaErrors(cublasDestroy(cublasHandle));     \
		if (cusparseHandle)   checkCudaErrors(cusparseDestroy(cusparseHandle)); \
		fflush (stdout);														\
	} while (0)

#endif	//CLEANDEVICE()

}

typedef enum InputFileResult
{
	IRRT_OK = 0,
	IRRT_NOTFOUND,
	IRRT_BAD_FORMAT
} InputFileResult_t;

typedef enum UserStatusError
{
	USER_STATUS_SUCCESS = 0,	//�ɹ�
	USER_STATUS_NOT_INITIAL,	//δ��ʼ��
	USER_STATUS_ALLOC_FAILED,	//����ռ�ʧ��
	USER_STATUS_INVALID_VALUE,	//����Чֵ
	USER_EXCUTION_FAILED,		//�˺�������ʧ��
	USER_EMPTY_POINT			//��ָ��
} UserStatusError_t;


typedef enum XmlResult
{
	XML_SUCCESS = 0,
	TABLE_ERROR,
	BLOCK_ERROR
}XmlResult_t;

#endif //!_HBXDEFMACRO_H