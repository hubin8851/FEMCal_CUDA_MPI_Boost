#pragma once
//注意，只有模板类卸载内联文件中，非模板类全部定义在.cpp中！否则编译不通过
#include "stdafx.h"

extern int g_MatRowNum;

/************************************************************************/
/* 控制台占位符相关操作，设置参数                                                                     */
/************************************************************************/
//@ mat_size：矩阵维度
//@ user_def：是否使用用户自定义矩阵维度
void	getMatrixSize(int argc, TCHAR* argv[], unsigned int &mat_size, unsigned int &user_def);

//@ _iter：迭代次数
//@ user_def：是否使用用户自定义迭代上限
void	getIterTiming( int argc, TCHAR* argv[], unsigned int &_iter );

//@ precision：迭代精度
//@ user_def：是否使用用户自定义迭代精度
void	getPrecision( int argc, TCHAR* argv[], float &precision, unsigned int &user_def);

//@ _filename：待写入结果的文件
void	getResultFilename( int argc, TCHAR* argv[], char *&_filename);

//@_Increment：thread的递增量
void	getIncrement( int argc, TCHAR* argv[], unsigned int &_Increment);

//判断long double型的两数字是否相等
bool	LongDoubleEqual(long double _lhs, long double _rhs);


/************************************************************************/
/*               内存管理相关子函数                                      */
/************************************************************************/
template<class T>
void FreeVal(T* _val){ free(_val) };


/************************************************************************/
/*               矩阵计算的全局函数                                     */
/************************************************************************/
//检验特征值的计算结果,CPU和GPU计算结果均可使用该函数，故置于全局函数头文件内
//@_nnz:CSR第一个数组和第二个数组的长度
//@_N：CSR第三个数组的长度
//@_ValIn：CSR第一个数组
//@_iColSort：CSR第二个数组
//@_iRowSort：CSR第三个数组
float	CheckResult(int _nnz,int _N, std::vector<float> _ValIn, std::vector<int> _iColSort, std::vector<int> _iRowSort, float *_x, float *_rhs);

//检验特征值的计算结果,CPU和GPU计算结果均可使用该函数，故置于全局函数头文件内
//@_nnz:CSR第一个数组和第二个数组的长度
//@_N：CSR第三个数组的长度
//@_ValIn：CSR第一个数组指针
//@_iColSort：CSR第二个数组指针
//@_iRowSort：CSR第三个数组指针
float	CheckResult(int _nnz,int _N, float *_ValIn, int *_iColSort, int *_iRowSort, float *_x, float *_rhs);

//检验特征值的计算结果,CPU和GPU计算结果均可使用该函数，故置于全局函数头文件内
//@_nnz:CSR第一个数组和第二个数组的长度
//@_N：CSR第三个数组的长度
//@_ValIn：CSR第一个数组指针
//@_iColSort：CSR第二个数组指针
//@_iRowSort：CSR第三个数组指针
double	CheckResult(int _nnz,int _N, double *_ValIn, int *_iColSort, int *_iRowSort, double *_x, double *_rhs);

//获取设备属性及个数
size_t		GetDevPropAndCount();

//boost索引目录并读取数据文件
boost::optional<boost::filesystem::path>	find_file(const boost::filesystem::path& _dir,const std::string& _filename);

//不区分大小写进行字符串比较
int myStrcmp(const char* _lhs, const char* _rhs);

namespace HBXDef
{
	/************************************************************************/
	/*                     函数指针相关,对比用                               */
	/************************************************************************/
	int cmp_cooFormat_csr( struct _CooFormat<int> *s, struct _CooFormat<int> *t);

	int cmp_cooFormat_csc( struct _CooFormat<int> *s, struct _CooFormat<int> *t);

	typedef int (*FUNPTR) (const void*, const void*);
	typedef int (*FUNPTR2) ( struct _CooFormat<int> *s, struct _CooFormat<int> *t);
	//函数指针列表
	static FUNPTR2  fptr_array[2] = {
		cmp_cooFormat_csr,
		cmp_cooFormat_csc,
	};

	/************************************************************************/
	/*                    数据读取相关                                      */
	/************************************************************************/

	//从文件中按行读取数据并放入数组指针
	template<class T>
	size_t	ReadVecByRowFromFile(std::vector<T> &_OutputPtr, const char *_FileName ,const boost::filesystem::path& _dir);
	//从文件中按行读取数据并放入数组指针
	template<class T>
	size_t	ReadVecByRowFromFile(thrust::host_vector<T> &_OutputPtr, const char *_FileName ,const boost::filesystem::path& _dir);

	template<class T>
	HBXDef::CMatrix<T>& GetMatFrmMtx( std::ifstream& _file, std::vector<std::string> &_vstrLine, int &_subMatDim );

	//参考NVIDIA的Samples下mmio.h文件
	//@filename：文件名
	//@csrFormat:是否为压缩格式
	//m,n,nnz:矩阵行、列、非零数
	//val,row,col:非零值向量，行索引向量，列索引向量
	int loadMMSparseMat( const char *filename, bool csrFormat,
						 int *m, int *n, int *nnz,
						 std::vector<double>& _vVal, std::vector<int>& _iRowSort, std::vector<int>& _iColSort );

	int mm_read_mtx_crd( const char *fname, int *_M, int *_N, int *_nz,
						std::vector<double>& _vVal, std::vector<int>& _iRowSort, std::vector<int>& _iColSort );
	//读取matrix market的尺寸
	int mm_read_mtx_crd_size( std::ifstream *fname, int *_M, int *_N, int *_nz );
	//读取相应三元组数据
	int mm_read_mtx_crd_data( std::ifstream *fname, int *_M, int *_N, int *_nz,
						std::vector<double>& _vVal, std::vector<int>& _RowSort, std::vector<int>& _iColSort );

	static void compress_index(	const int *Ind, int nnz, int m, int *Ptr, int base)
	{
		int i;
		/* initialize everything to zero */
		for(i=0; i<m+1; i++){
			Ptr[i]=0;
		} 
		/* count elements in every row */
		Ptr[0]=base;
		for(i=0; i<nnz; i++){
			Ptr[Ind[i]+(1-base)]++;
		} 
		/* add all the values */
		for(i=0; i<m; i++){
			Ptr[i+1]+=Ptr[i];
		} 
	}

	//检验CSR格式矩阵的行列号是否有问题
	static CSRFormatError_t verify_pattern( int m, int nnz, int *csrRowPtr, int *csrColInd )
	{
		int i, col, start, end, base_index;
		int error_found = 0;

		if (nnz != (csrRowPtr[m] - csrRowPtr[0])){
			fprintf(stderr, "Error (nnz check failed): (csrRowPtr[%d]=%d - csrRowPtr[%d]=%d) != (nnz=%d)\n", 0, csrRowPtr[0], m, csrRowPtr[m], nnz);
			error_found = NNZCHECKFAILED;
		}

		base_index = csrRowPtr[0];
		if ((0 != base_index) && (1 != base_index)){
			fprintf(stderr, "Error (base index check failed): base index = %d\n", base_index);
			error_found = BASEINDEXCHECKFAILED;
		}
		for (i=0; (!error_found) && (i<m); i++){
			start = csrRowPtr[i  ] - base_index;
			end   = csrRowPtr[i+1] - base_index;
			if (start > end){
				fprintf(stderr, "Error (corrupted row): csrRowPtr[%d] (=%d) > csrRowPtr[%d] (=%d)\n", i, start+base_index, i+1, end+base_index);
				error_found = CORRUPTEDROW;
			}
			for (col=start; col<end; col++){
				if (csrColInd[col] < base_index){
					fprintf(stderr, "Error (column vs. base index check failed): csrColInd[%d] < %d\n", col, base_index);
					error_found = COLVSBASEINDEXCHECK;
				}
				if ((col < (end-1)) && (csrColInd[col] >= csrColInd[col+1])){
					fprintf(stderr, "Error (sorting of the column indecis check failed): (csrColInd[%d]=%d) >= (csrColInd[%d]=%d)\n", col, csrColInd[col], col+1, csrColInd[col+1]);
					error_found = SORTINGCOLINDEXCHECK;
				}
			}
		}
		return CHECKSUCCESS;
	}

	//Vector至向量的内存拷贝
	template<class T>
	void		MemcpyArrayTransToVec( T* _SourceArray, size_t _length, std::vector<T>& _DestVec);
	//向量至Vector的内存拷贝
	template<class T>
	int		MemcpyVecTransToArray(std::vector<T> &_vecIn,  T  *_arrayOut);

	//生成三对角正定对称阵
	void genTridiag( float *val, int *I, int *J, int _nz, int _N );

	//输出容器的数据
	template<class T>
	UserStatusError_t		OutputVectorData(std::vector<T> _OutputVec, const char *_FileName);
	//输出数组的数据
	template<class T>
	UserStatusError_t		OutputPointData(T*	_OutputPoint, size_t _length, const char *_ExportFileName = "LoadVec.data", boost::filesystem::path _dir="..\\data\\export" );
	//将CSR格式的三个数组分别记录进三个文件
	//@_Val：非零向量
	//@_ColSort：非零向量列索引
	//@_ValRowSort：非零向量行索引
	//@_nnA：_Val长度
	//@_N：矩阵的size
	template<class T1, class T2>
	void		_ExportCSRdataToFile(T1 *_Val,			
													T2 *_ColSort,		
													T2 *_ValRowSort,
													int _nnA, int _N,
													const char*	_NoneZeroValFileName = "NoneZeroVal.data",
													const char*	_ColSortFileName = "ColSort.data",
													const char*	_NoneZeroRowSortFileName = "NoneZeroRowSort.data",
													const boost::filesystem::path& _dir = "..\\data\\export");
	//将CSR格式的三个数组分别记录进三个文件
	//@_vecVal：非零向量数组
	//@_vecColSort：非零向量列索引
	//@_vecValRowSort：非零向量行索引
	template<class T1, class T2>
	void		_ExportCSRdataToFile(std::vector<T1> &_vecVal,			
							  std::vector<T2> &_vecColSort,		
							  std::vector<T2> &_vecValRowSort,
							  const char*	_NoneZeroValFileName = "NoneZeroVal.data",
							  const char*	_ColSortFileName = "ColSort.data",
							  const char*	_NoneZeroRowSortFileName = "NoneZeroRowSort.data",
							  const boost::filesystem::path& _dir = "..\\data\\export");

	//导入CSR格式文本的数据并存入三个数组,thrust版
	//@_vecVal：非零向量数组
	//@_vecColSort：非零向量列索引
	//@_vecValRowSort：非零向量行索引
	template<class T1,class T2>
	size_t	_ImportCSRdataFromFile(	thrust::host_vector<T1> &_vecVal,			
									thrust::host_vector<T2> &_vecColSort,		
									thrust::host_vector<T2> &_vecValRowSort,
									const char*	_NoneZeroValFileName = "NoneZeroVal.data",
									const char*	_ColSortFileName = "ColSort.data",
									const char*	_NoneZeroRowSortFileName = "NoneZeroRowSort.data",
									const boost::filesystem::path& _dir = "..\\data\\source");
	//导入CSR格式文本的数据并存入三个数组
	//@_vecVal：非零向量数组
	//@_vecColSort：非零向量列索引
	//@_vecValRowSort：非零向量行索引
	template<class T1, class T2>
	size_t	_ImportCSRdataFromFile(std::vector<T1> &_vecVal,			
								   std::vector<T2> &_vecColSort,		
								   std::vector<T2> &_vecValRowSort,
								   const char*	_NoneZeroValFileName = "NoneZeroVal.data",
								   const char*	_ColSortFileName = "ColSort.data",
								   const char*	_NoneZeroRowSortFileName = "NoneZeroRowSort.data",
								   const boost::filesystem::path& _dir = "..\\data\\source");

	//将完整的矩阵转换成CSR形式
	//@_InCMatric	输入矩阵类
	//@_dNonZeroVal CSR模式第一个<T>容器
	//@_iColSort	CSR模式第二个int容器
	//@_iNonZeroRowSort	CSR模式第三个int容器
	//参见基于CUDA的大规模稀疏矩阵的PCG算法优化-郑经纬-清华大学学报
	template<class _iT, class _oT>
	bool	_FullMatricToCSR( HBXDef::CMatrix<_iT>& _InCMatric, 
		std::vector<_oT>& _dNonZeroVal, 
		std::vector<int>& _iColSort, 
		std::vector<int>& _iNonZeroRowSort);


	/************************************************************************/
	/*     GPU部分                                                          */
	/************************************************************************/
	double GetTimeStamp(void);
}

#include "HbxGloFunc.inl" //类及相关函数的定义头文件

