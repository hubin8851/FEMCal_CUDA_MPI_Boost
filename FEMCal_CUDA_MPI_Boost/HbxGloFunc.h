#pragma once
//ע�⣬ֻ��ģ����ж�������ļ��У���ģ����ȫ��������.cpp�У�������벻ͨ��
#include "stdafx.h"

extern int g_MatRowNum;

/************************************************************************/
/* ����̨ռλ����ز��������ò���                                                                     */
/************************************************************************/
//@ mat_size������ά��
//@ user_def���Ƿ�ʹ���û��Զ������ά��
void	getMatrixSize(int argc, TCHAR* argv[], unsigned int &mat_size, unsigned int &user_def);

//@ _iter����������
//@ user_def���Ƿ�ʹ���û��Զ����������
void	getIterTiming( int argc, TCHAR* argv[], unsigned int &_iter );

//@ precision����������
//@ user_def���Ƿ�ʹ���û��Զ����������
void	getPrecision( int argc, TCHAR* argv[], float &precision, unsigned int &user_def);

//@ _filename����д�������ļ�
void	getResultFilename( int argc, TCHAR* argv[], char *&_filename);

//@_Increment��thread�ĵ�����
void	getIncrement( int argc, TCHAR* argv[], unsigned int &_Increment);

//�ж�long double�͵��������Ƿ����
bool	LongDoubleEqual(long double _lhs, long double _rhs);


/************************************************************************/
/*               �ڴ��������Ӻ���                                      */
/************************************************************************/
template<class T>
void FreeVal(T* _val){ free(_val) };


/************************************************************************/
/*               ��������ȫ�ֺ���                                     */
/************************************************************************/
//��������ֵ�ļ�����,CPU��GPU����������ʹ�øú�����������ȫ�ֺ���ͷ�ļ���
//@_nnz:CSR��һ������͵ڶ�������ĳ���
//@_N��CSR����������ĳ���
//@_ValIn��CSR��һ������
//@_iColSort��CSR�ڶ�������
//@_iRowSort��CSR����������
float	CheckResult(int _nnz,int _N, std::vector<float> _ValIn, std::vector<int> _iColSort, std::vector<int> _iRowSort, float *_x, float *_rhs);

//��������ֵ�ļ�����,CPU��GPU����������ʹ�øú�����������ȫ�ֺ���ͷ�ļ���
//@_nnz:CSR��һ������͵ڶ�������ĳ���
//@_N��CSR����������ĳ���
//@_ValIn��CSR��һ������ָ��
//@_iColSort��CSR�ڶ�������ָ��
//@_iRowSort��CSR����������ָ��
float	CheckResult(int _nnz,int _N, float *_ValIn, int *_iColSort, int *_iRowSort, float *_x, float *_rhs);

//��������ֵ�ļ�����,CPU��GPU����������ʹ�øú�����������ȫ�ֺ���ͷ�ļ���
//@_nnz:CSR��һ������͵ڶ�������ĳ���
//@_N��CSR����������ĳ���
//@_ValIn��CSR��һ������ָ��
//@_iColSort��CSR�ڶ�������ָ��
//@_iRowSort��CSR����������ָ��
double	CheckResult(int _nnz,int _N, double *_ValIn, int *_iColSort, int *_iRowSort, double *_x, double *_rhs);

//��ȡ�豸���Լ�����
size_t		GetDevPropAndCount();

//boost����Ŀ¼����ȡ�����ļ�
boost::optional<boost::filesystem::path>	find_file(const boost::filesystem::path& _dir,const std::string& _filename);

//�����ִ�Сд�����ַ����Ƚ�
int myStrcmp(const char* _lhs, const char* _rhs);

namespace HBXDef
{
	/************************************************************************/
	/*                     ����ָ�����,�Ա���                               */
	/************************************************************************/
	int cmp_cooFormat_csr( struct _CooFormat<int> *s, struct _CooFormat<int> *t);

	int cmp_cooFormat_csc( struct _CooFormat<int> *s, struct _CooFormat<int> *t);

	typedef int (*FUNPTR) (const void*, const void*);
	typedef int (*FUNPTR2) ( struct _CooFormat<int> *s, struct _CooFormat<int> *t);
	//����ָ���б�
	static FUNPTR2  fptr_array[2] = {
		cmp_cooFormat_csr,
		cmp_cooFormat_csc,
	};

	/************************************************************************/
	/*                    ���ݶ�ȡ���                                      */
	/************************************************************************/

	//���ļ��а��ж�ȡ���ݲ���������ָ��
	template<class T>
	size_t	ReadVecByRowFromFile(std::vector<T> &_OutputPtr, const char *_FileName ,const boost::filesystem::path& _dir);
	//���ļ��а��ж�ȡ���ݲ���������ָ��
	template<class T>
	size_t	ReadVecByRowFromFile(thrust::host_vector<T> &_OutputPtr, const char *_FileName ,const boost::filesystem::path& _dir);

	template<class T>
	HBXDef::CMatrix<T>& GetMatFrmMtx( std::ifstream& _file, std::vector<std::string> &_vstrLine, int &_subMatDim );

	//�ο�NVIDIA��Samples��mmio.h�ļ�
	//@filename���ļ���
	//@csrFormat:�Ƿ�Ϊѹ����ʽ
	//m,n,nnz:�����С��С�������
	//val,row,col:����ֵ����������������������������
	int loadMMSparseMat( const char *filename, bool csrFormat,
						 int *m, int *n, int *nnz,
						 std::vector<double>& _vVal, std::vector<int>& _iRowSort, std::vector<int>& _iColSort );

	int mm_read_mtx_crd( const char *fname, int *_M, int *_N, int *_nz,
						std::vector<double>& _vVal, std::vector<int>& _iRowSort, std::vector<int>& _iColSort );
	//��ȡmatrix market�ĳߴ�
	int mm_read_mtx_crd_size( std::ifstream *fname, int *_M, int *_N, int *_nz );
	//��ȡ��Ӧ��Ԫ������
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

	//����CSR��ʽ��������к��Ƿ�������
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

	//Vector���������ڴ濽��
	template<class T>
	void		MemcpyArrayTransToVec( T* _SourceArray, size_t _length, std::vector<T>& _DestVec);
	//������Vector���ڴ濽��
	template<class T>
	int		MemcpyVecTransToArray(std::vector<T> &_vecIn,  T  *_arrayOut);

	//�������Խ������Գ���
	void genTridiag( float *val, int *I, int *J, int _nz, int _N );

	//�������������
	template<class T>
	UserStatusError_t		OutputVectorData(std::vector<T> _OutputVec, const char *_FileName);
	//������������
	template<class T>
	UserStatusError_t		OutputPointData(T*	_OutputPoint, size_t _length, const char *_ExportFileName = "LoadVec.data", boost::filesystem::path _dir="..\\data\\export" );
	//��CSR��ʽ����������ֱ��¼�������ļ�
	//@_Val����������
	//@_ColSort����������������
	//@_ValRowSort����������������
	//@_nnA��_Val����
	//@_N�������size
	template<class T1, class T2>
	void		_ExportCSRdataToFile(T1 *_Val,			
													T2 *_ColSort,		
													T2 *_ValRowSort,
													int _nnA, int _N,
													const char*	_NoneZeroValFileName = "NoneZeroVal.data",
													const char*	_ColSortFileName = "ColSort.data",
													const char*	_NoneZeroRowSortFileName = "NoneZeroRowSort.data",
													const boost::filesystem::path& _dir = "..\\data\\export");
	//��CSR��ʽ����������ֱ��¼�������ļ�
	//@_vecVal��������������
	//@_vecColSort����������������
	//@_vecValRowSort����������������
	template<class T1, class T2>
	void		_ExportCSRdataToFile(std::vector<T1> &_vecVal,			
							  std::vector<T2> &_vecColSort,		
							  std::vector<T2> &_vecValRowSort,
							  const char*	_NoneZeroValFileName = "NoneZeroVal.data",
							  const char*	_ColSortFileName = "ColSort.data",
							  const char*	_NoneZeroRowSortFileName = "NoneZeroRowSort.data",
							  const boost::filesystem::path& _dir = "..\\data\\export");

	//����CSR��ʽ�ı������ݲ�������������,thrust��
	//@_vecVal��������������
	//@_vecColSort����������������
	//@_vecValRowSort����������������
	template<class T1,class T2>
	size_t	_ImportCSRdataFromFile(	thrust::host_vector<T1> &_vecVal,			
									thrust::host_vector<T2> &_vecColSort,		
									thrust::host_vector<T2> &_vecValRowSort,
									const char*	_NoneZeroValFileName = "NoneZeroVal.data",
									const char*	_ColSortFileName = "ColSort.data",
									const char*	_NoneZeroRowSortFileName = "NoneZeroRowSort.data",
									const boost::filesystem::path& _dir = "..\\data\\source");
	//����CSR��ʽ�ı������ݲ�������������
	//@_vecVal��������������
	//@_vecColSort����������������
	//@_vecValRowSort����������������
	template<class T1, class T2>
	size_t	_ImportCSRdataFromFile(std::vector<T1> &_vecVal,			
								   std::vector<T2> &_vecColSort,		
								   std::vector<T2> &_vecValRowSort,
								   const char*	_NoneZeroValFileName = "NoneZeroVal.data",
								   const char*	_ColSortFileName = "ColSort.data",
								   const char*	_NoneZeroRowSortFileName = "NoneZeroRowSort.data",
								   const boost::filesystem::path& _dir = "..\\data\\source");

	//�������ľ���ת����CSR��ʽ
	//@_InCMatric	���������
	//@_dNonZeroVal CSRģʽ��һ��<T>����
	//@_iColSort	CSRģʽ�ڶ���int����
	//@_iNonZeroRowSort	CSRģʽ������int����
	//�μ�����CUDA�Ĵ��ģϡ������PCG�㷨�Ż�-֣��γ-�廪��ѧѧ��
	template<class _iT, class _oT>
	bool	_FullMatricToCSR( HBXDef::CMatrix<_iT>& _InCMatric, 
		std::vector<_oT>& _dNonZeroVal, 
		std::vector<int>& _iColSort, 
		std::vector<int>& _iNonZeroRowSort);


	/************************************************************************/
	/*     GPU����                                                          */
	/************************************************************************/
	double GetTimeStamp(void);
}

#include "HbxGloFunc.inl" //�༰��غ����Ķ���ͷ�ļ�

