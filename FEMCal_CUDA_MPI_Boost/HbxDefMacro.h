#ifndef _HBXDEFMACRO_H
#define	_HBXDEFMACRO_H

 /******
     由于NDEBUG是为程序调试(debug)期间使用的，当调试期完毕，程序出
 发行版(release)后，NDEBUG将失去作用。为了能够使assert()函数在发行
 版中也可以使用，则定义了下面的条件编译NDEBUG，意思为：
     如果已经定义NDEBUG，则定义宏Assert(x)：当x为假时，执行throw，
 交由系统已定义的方式处理；否则将Assert(x)定义为函数assert(x)，该
 函数对 参数x进行测试：x为假时，函数终止程序并打印错误信息。
 *****/
 #ifdef NDEBUG			
  #define Assert(x)	if (!x)	throw;	
 #else			//否则用系统定义的函数assert()处理
  #include <cassert>
  #define Assert(x)	assert(x);
 #endif			//NDEBUG

#ifndef CUDA_CALL(x)	//判断CUDA相关运行函数返回值
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

	const int MAX_GPU_COUNT = 32;	//单机最大设备数，目前也只支持最大32个设备

	//数据运算精度枚举
	typedef enum
	{
		_PRECS_FLOAT_ = 1,
		_PRECS_DOUBLE_ = 2
	} PrecsType_t;

	//主行或主列索引枚举
	typedef enum 
	{
		ROWMAJOR = 1,
		COLMAJOR = 2
	} SortType_t;	

	//矩阵格式起始索引枚举
	typedef enum IndexBase
	{
		INDEX_BASE_ZERO = 0, 
		INDEX_BASE_ONE = 1,
		INDEX_INIVIAL
	} IndexBase_t;

	//数据拷贝类型枚举
	typedef enum CopyType
	{
		CopyFailed = -1,
		HostToDevice = 0,
		DeviceToHost = 1,
		HostToHost = 2,
		DeviceToDevice = 3
	} CopyType_t;	

	//重置内存/显存标志位枚举
	typedef	enum ResetType
	{
		RESETCPUMEM = 0,
		RESETGPUMEM = 1
	} ResetType_t;	

	//载荷类型枚举
	typedef enum LoadType
	{
		TRACTION = 0,
		CONCENTRATE = 1
	} LoadType_t;	//

	typedef enum
	{
		NORMAL = 0,		//普通内存分配
		UNIFIEDMEM = 1,	//统一内存
		PAGELOCK,		//锁页内存
		ARRAY2D,		//2维向量
	} CudaMalloc_t;

	//存储格式枚举
	typedef	enum SaveType
	{
		CSR = 0,
		CSC,
		COO,
		DIA,	//对角格式
		ELL,	//按行格式
		JAD		//齿对角格式
	} SaveType_t;

	typedef enum CSRFormatError
	{
		CHECKSUCCESS = 0,
		NNZCHECKFAILED = 1,
		BASEINDEXCHECKFAILED,
		CORRUPTEDROW,
		COLVSBASEINDEXCHECK,
		SORTINGCOLINDEXCHECK
	} CSRFormatError_t;	//CSR格式导入错误

	//矩阵错误
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
	} FactorError_t;	//因式分解错误类型

	typedef enum DataAloc
	{
		INVALID = -1,
		DATAINMEM = 0,
		DATAINGPU = 1
	} DataAloc_t;	//检验数据所在位置

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

	//迭代算法枚举
	typedef enum SolveMethod
	{
		//直接法
		QR_Mtd = 0,
		//迭代法
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

//若采用锁页内存，调用此宏
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

//参见http://blog.csdn.net/stpeace/article/details/9063293,清除流的buffer内所有信息
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
	USER_STATUS_SUCCESS = 0,	//成功
	USER_STATUS_NOT_INITIAL,	//未初始化
	USER_STATUS_ALLOC_FAILED,	//分配空间失败
	USER_STATUS_INVALID_VALUE,	//非有效值
	USER_EXCUTION_FAILED,		//核函数调用失败
	USER_EMPTY_POINT			//空指针
} UserStatusError_t;


typedef enum XmlResult
{
	XML_SUCCESS = 0,
	TABLE_ERROR,
	BLOCK_ERROR
}XmlResult_t;

#endif //!_HBXDEFMACRO_H