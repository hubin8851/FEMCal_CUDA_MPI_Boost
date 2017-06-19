/************************************************************************/
/*  HBX,2015-10-10                                                      */
/************************************************************************/
#include "stdafx.h"
#include "HbxGloFunc.h"


/************************************************************************/
/*                                  获取控制台参数占位符等相关操作                                   */
/************************************************************************/

//@ mat_size：矩阵维度
//@ user_def：是否适应用户自定义矩阵维度
void	getMatrixSize(int argc, TCHAR* argv[], unsigned int &mat_size, unsigned int &user_def)
{
	int tmp = -1;
	if ( checkCmdLineFlag(argc, (const char**)argv, "matrix-size") )
	{
		tmp = getCmdLineArgumentInt(argc, (const char **) argv, "matrix-size" );
	}
	if (tmp > 0)
	{
		mat_size = (unsigned int) tmp;
		//全局变量g_Rownum为unsigned int
		assert(mat_size < (1 << 64));
		//矩阵大小应大于2
		assert(mat_size >= 2);
		user_def = 1;
	}
	std::cout<<"矩阵维度为:"<<mat_size<<"x"<<mat_size<<std::endl;
}

//@ _iter：迭代次数
//@ user_def：是否使用用户自定义迭代上限
void getIterTiming( int argc, TCHAR* argv[], unsigned int &_iters )
{
	int tmp = -1;
	if ( checkCmdLineFlag(argc, (const char**)argv, "iters") )
	{
		tmp = getCmdLineArgumentInt(argc, (const char **) argv, "iters" );
	}
	if (tmp > 0)
	{
		_iters = tmp;
	}
	std::cout<<"迭代次数上限为:"<<_iters<<std::endl;
}

//@ precision：迭代精度
//@ user_def：是否使用用户自定义迭代精度
void	getPrecision(int argc, TCHAR* argv[], float &precision, unsigned int &user_def)
{
	float tmp = -1;
	if ( checkCmdLineFlag(argc, (const char**)argv, "precision") )
	{
		tmp = getCmdLineArgumentFloat(argc, (const char **) argv, "precision" );
	}
	if (tmp > 0)
	{
		precision = tmp;
	}
	std::cout<<"迭代精度为:"<<precision<<std::endl;
}

//@ _filename：待写入结果的文件
void	getResultFilename(int argc, TCHAR* argv[], char *&_filename)
{
	char *tmp = nullptr;
	if ( checkCmdLineFlag(argc, (const char**)argv, "ResultFile") )
	{
		getCmdLineArgumentString(argc, (const char **) argv, "ResultFile", &tmp );
	}
	if ( nullptr != tmp)
	{
		_filename = (char *) malloc(sizeof(char) * strlen(tmp));
		strcpy(_filename, tmp);
		free(tmp);
	}
	std::cout<<"结果存储文件名为:"<<_filename<<std::endl;
}

//@_Increment：thread的递增量
void getIncrement(int argc, TCHAR* argv[], unsigned int &_Increment)
{
	unsigned int tmp = -1;
	if ( checkCmdLineFlag(argc, (const char**)argv, "iters") )
	{
		tmp = getCmdLineArgumentInt(argc, (const char **) argv, "iters" );
	}
	if (tmp > 0)
	{
		_Increment = (unsigned int)tmp;
	}
	std::cout<<"每次thread递增量为:"<<_Increment<<std::endl;
}


//判断long double型的两数字是否相等
bool	LongDoubleEqual(long double _lhs, long double _rhs)
{
	if ( fabsl(_lhs - _rhs) < LONGDOUBLEERROR )
		return true;
	else
		return false;
}



//检验特征值的计算结果,使用thrust版本
//@_nnz:CSR第一个数组和第二个数组的长度
//@_N：CSR第三个数组的长度
//@_ValIn：CSR第一个数组
//@_iColSort：CSR第二个数组
//@_iRowSort：CSR第三个数组
float	CheckResult(int _nnz,int _N, thrust::host_vector<float> _ValIn, thrust::host_vector<int> _iColSort, thrust::host_vector<int> _iRowSort, thrust::host_vector<float> _x, thrust::host_vector<float> _rhs)
{
	float err = 0.0;

	for (int i = 0; i < g_MatRowNum; i++)
	{
		float rsum = 0.0;

		for (int j = _iRowSort[i]; j < _iRowSort[i+1]; j++)
		{
			rsum += _ValIn[j]*_x[_iColSort[j]];
		}

		float diff = fabs(rsum - _rhs[i]);

		if (diff > err)
		{
			err = diff;
		}
	}
	return err;
}

//检验特征值的计算结果，使用裸指针
//@_nnz:CSR第一个数组和第二个数组的长度
//@_N：CSR第三个数组的长度
//@_ValIn：CSR第一个数组指针
//@_iColSort：CSR第二个数组指针
//@_iRowSort：CSR第三个数组指针
float	CheckResult(int _nnz,int _N, float *_ValIn, int *_iColSort, int *_iRowSort, float *_x, float *_rhs)
{
	float err = 0.0;

	for (int i = 0; i < g_MatRowNum; i++)
	{
		float rsum = 0.0;

		for (int j = _iRowSort[i]; j < _iRowSort[i+1]; j++)
		{
			rsum += _ValIn[j]*_x[_iColSort[j]];
		}

		float diff = fabs(rsum - _rhs[i]);

		if (diff > err)
		{
			err = diff;
		}
	}
	return err;
}

//检验特征值的计算结果,CPU和GPU计算结果均可使用该函数，故置于全局函数头文件内
//@_nnz:CSR第一个数组和第二个数组的长度
//@_N：CSR第三个数组的长度
//@_ValIn：CSR第一个数组指针
//@_iColSort：CSR第二个数组指针
//@_iRowSort：CSR第三个数组指针
double	CheckResult(int _nnz,int _N, double *_ValIn, int *_iColSort, int *_iRowSort, double *_x, double *_rhs)
{
	double err = 0.0;

	for (int i = 0; i < g_MatRowNum; i++)
	{
		double rsum = 0.0;

		for (int j = _iRowSort[i]; j < _iRowSort[i+1]; j++)
		{
			rsum += _ValIn[j]*_x[_iColSort[j]];
		}

		double diff = fabs(rsum - _rhs[i]);

		if (diff > err)
		{
			err = diff;
		}
	}
	return err;
}

//boost索引目录并读取数据文件
boost::optional<boost::filesystem::path>	find_file(const boost::filesystem::path& _dir,const std::string& _filename)
{
	using namespace boost::filesystem;
	typedef boost::optional<boost::filesystem::path> result_type;
	typedef recursive_directory_iterator rd_iterator;
	if ( !exists(_dir)||!is_directory(_dir))
	{
		return result_type();
	}
	rd_iterator end;
	for (rd_iterator pos(_dir); pos!=end; ++pos)
	{
		if (!is_directory(*pos)&&pos->path().filename() == _filename)
		{
			return result_type(pos->path());
		}
	}
	return result_type();
}

namespace HBXDef
{
	int cmp_cooFormat_csr( struct _CooFormat<int> *s, struct _CooFormat<int> *t)
	{
		if ( s->_i < t->_i ){
			return -1 ;
		}
		else if ( s->_i > t->_i ){
			return 1 ;
		}
		else{
			return s->_j - t->_j ;
		}
	}

	int cmp_cooFormat_csc( struct _CooFormat<int> *s, struct _CooFormat<int> *t)
	{
		if ( s->_j < t->_j ){
			return -1 ;
		}
		else if ( s->_j > t->_j ){
			return 1 ;
		}
		else{
			return s->_i - t->_i ;
		}
	}


	//@filename：文件名
	//@csrFormat:是否为压缩格式
	//m,n,nnz:矩阵行、列、非零数
	//val,row,col:非零值向量，行索引向量，列索引向量
	int loadMMSparseMat( const char *filename, bool csrFormat,
						 int *m, int *n, int *nnz,
						 std::vector<double>& _vVal, std::vector<int>& _iRowSort, std::vector<int>& _iColSort )
 	{
		IndexBase_t  _Base;
		HBXDef::_CooFormat<int>	*_work;
		std::vector<double> _tmpVal;
		std::vector<int>	_tmpRow, _tmpCol;
 		//读取mtx文件
 		int error = mm_read_mtx_crd( filename, m, n, nnz, _tmpVal, _tmpRow, _tmpCol );
 		if (error) 
		{
 			fprintf(stderr, "!!!! can not open file: '%s'\n", filename);
 			return 1;       
 		}

		//检查读取错误
		 _work = (struct _CooFormat<int> *)malloc(sizeof(struct _CooFormat<int>)*(*nnz));
		 if (NULL == _work){
			 fprintf(stderr, "!!!! allocation error, malloc failed\n");
			 return 1;
		 }
		 for(int i=0; i<(*nnz); i++){
			 _work[i]._i = _tmpRow[i];
			 _work[i]._j = _tmpCol[i];
			 _work[i]._val = i; // permutation is identity
		 }
		 if (csrFormat){
			 /* create row-major ordering of indices (sorted by row and within each row by column) */
			 qsort(_work, *nnz, sizeof(struct _CooFormat<int>), (FUNPTR)fptr_array[0] );
		 }else{
			 /* create column-major ordering of indices (sorted by column and within each column by row) */
			 qsort(_work, *nnz, sizeof(struct _CooFormat<int>), (FUNPTR)fptr_array[1] );
		 }

		 for(int i=0; i<(*nnz); i++){
			 _tmpRow[i] = _work[i]._i;
			 _tmpCol[i] = _work[i]._j;
		 }

		 if ( 0 == _work[0]._i || 0 == _work[0]._j )
		 {
			 _Base = IndexBase_t::INDEX_BASE_ZERO;
		 }
		 else if ( 1 == _work[0]._i || 1 == _work[0]._j )
		 {
			  _Base = IndexBase_t::INDEX_BASE_ONE;
		 }
		 else
		 {
			 printf("传入矩阵首索引既非0也非1，请检查... \n");
			 return 1;
		 }
		 if (csrFormat)
		 {   //获得压缩格式的行索引信息
			 //对输出的向量重新调整大小
			 _vVal.resize(*nnz);
			 _iRowSort.resize( (*m)+1 );
			 compress_index( _tmpRow.data(), *nnz, *m, _iRowSort.data(), _Base );
			 _iColSort.resize( *nnz );
		 }
		 else
		 {
			 //CSC格式的相关计算
		 }
		 for (int i=0; i<(*nnz); i++) 
		 {
			 if (csrFormat)
			 {
				 _iColSort = _tmpCol;
				 _vVal[i] =  _tmpVal[ _work[i]._val ] ;
			 }
			 
		 }
		 //数据读取检验
		 int error_found;
		 if (csrFormat){
			 error_found = verify_pattern(*m, *nnz, _iRowSort.data(), _iColSort.data() );
		 }else{
			 error_found = verify_pattern(*n, *nnz, _iColSort.data(), _iRowSort.data());
		 }
		 if (error_found){
			 fprintf(stderr, "CSR或CSC检验完成，格式读取有误...\n");
			 return 1;
		 }
		 free(_work);

		 return 0;
 	}

	int mm_read_mtx_crd(const char *fname, int *_M, int *_N, int *_nz,
						std::vector<double>& _vVal, std::vector<int>& _iRowSort, std::vector<int>& _iColSort)
	{
		int ret_code;
		std::ifstream  inFile;
		std::string stringLine;
		std::vector<std::string> vstrLine;
		inFile.open( fname );
		if (!inFile.is_open())
			return MMError_t::MM_COULD_NOT_READ_FILE;

		if ( (ret_code = mm_read_mtx_crd_size(&inFile, _M, _N, _nz)) != 0 )
			return ret_code;
		_vVal.resize(*_nz);
		_iRowSort.resize(*_nz);
		_iColSort.resize(*_nz);
		if ( (ret_code = mm_read_mtx_crd_data(&inFile, _M, _N, _nz, _vVal, _iRowSort, _iColSort)) != 0 )
			return ret_code;
		return 0;
	}

	int mm_read_mtx_crd_size( std::ifstream *fname, int *_M, int *_N, int *_nz )
	{
		std::string	 stringLine;
		//先将各参数置0
		*_M = *_N = *_nz = 0;
		do 
		{
			if ( NULL == getline(*fname, stringLine) )
				return MM_PREMATURE_EOF;
		}while (stringLine[0] == '%');
		if (sscanf(stringLine.c_str(), "%d %d %d", _M, _N, _nz) == 3)
			return 0;
	}

	//读取相应三元组数据
	int mm_read_mtx_crd_data( std::ifstream *fname, int *_M, int *_N, int *_nz,
								std::vector<double>& _vVal, std::vector<int>& _RowSort, std::vector<int>& _iColSort )
	{
		int i;
		std::string stringLine;
		std::vector<std::string> vstrLine;
		for (i=0; i<*_nz; i++)
		{
			std::getline( *fname, stringLine );
			if (sscanf(stringLine.c_str(), "%d %d %lg\n", &_RowSort[i], &_iColSort[i], &_vVal[i]) != 3) 
				return MM_PREMATURE_EOF;

		}
		return 0;
	}


	double GetTimeStamp(void)
	{
		LARGE_INTEGER t;
		static double oofreq;
		static int checkedForHighResTimer;
		static BOOL hasHighResTimer;

		if (!checkedForHighResTimer) {
			hasHighResTimer = QueryPerformanceFrequency (&t);
			oofreq = 1.0 / (double)t.QuadPart;
			checkedForHighResTimer = 1;
		}
		if (hasHighResTimer) {
			QueryPerformanceCounter (&t);
			return (double)t.QuadPart * oofreq;
		} else {
			return (double)GetTickCount() / 1000.0;
		}
	}



	//生成三对角正定对称阵
	void genTridiag( float *val, int *I, int *J, int _nz, int _N )
	{
		I[0] = 0, J[0] = 0, J[1] = 1;
		val[0] = (float)rand()/RAND_MAX + 10.0f;
		val[1] = (float)rand()/RAND_MAX;
		int start;

		for (int i = 1; i < _N; i++)
		{
			if (i > 1)
			{
				I[i] = I[i-1]+3;
			}
			else
			{
				I[1] = 2;
			}

			start = (i-1)*3 + 2;
			J[start] = i - 1;
			J[start+1] = i;

			if (i < _N-1)
			{
				J[start+2] = i + 1;
			}

			val[start] = val[start-1];
			val[start+1] = (float)rand()/RAND_MAX + 10.0f;

			if (i < _N-1)
			{
				val[start+2] = (float)rand()/RAND_MAX;
			}
		}

		I[_N] = _nz;
	}
}


