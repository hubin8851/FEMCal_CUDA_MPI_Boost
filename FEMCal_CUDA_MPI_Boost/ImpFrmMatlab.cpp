#include "StdAfx.h"
#include "HbxGloFunc.h"
#include "ImpFrmMatlab.h"


CImpFrmMatlab::CImpFrmMatlab(void):m_pmatFile(nullptr),m_pMxArray(nullptr),
								m_RowNum(0),m_ColNum(0),
								//_pNoneZeroVal(nullptr), _piColSort(nullptr), _piNonZeroRowSort(nullptr),
								m_nnzA(0), m_nA(0)
{

}


CImpFrmMatlab::~CImpFrmMatlab(void)
{
	m_pmatFile = nullptr;
	m_pMxArray = nullptr;
	//_pNoneZeroVal = nullptr;
	//_piColSort = nullptr;
	//_piNonZeroRowSort = nullptr;
}

//该函数读入某一路径下的函数并会在子目录下索引
//@_matname:矩阵名称
//@_searchpath：索引路径
ImpMatError_t CImpFrmMatlab::ImpFrmMatlab( const char* _matname,  boost::filesystem::path _searchpath)
{

	boost::optional<boost::filesystem::path> r = find_file(_searchpath,_matname);
	if (!r)
	{
		std::cerr<<"boost搜索当前路径"<<_searchpath<<"下无该文件"<<std::endl;
		return ImpMatError_t::MATOPENFAILD;
	}
	else
	{
		std::cout<<"当前文件所在路径"<<*r<<std::endl;
	}
	int ndir;	//根目录下mxarray数目
	int nfiled = -10;	//mxarray下file数目  
	m_pmatFile = matOpen( r->string().c_str(),"r" );
	if (nullptr == m_pmatFile)
	{
		std::cerr<<"当前路径索引不正确未找到mat文件"<<std::endl;
		return ImpMatError_t::MATGETDIRERROR;
	}
	const char** dir = (const char **)matGetDir(m_pmatFile, &ndir);
	if ( nullptr == dir )
	{
		if ( 0 == ndir )
		{
			std::cerr<<"当前路径"<<_searchpath<<"下无该文件"<<std::endl;
			return ImpMatError_t::MATFILENULL;
		}
		else if ( 0 > ndir)
		{
			std::cerr<<"读取mat文件时发生错误..."<<std::endl;
			return ImpMatError_t::MATGETVARERROR;
		}
	}
 	else if ( boost::iequals( "problem", dir[0] ))
 	{
		for (int i = 0; i < ndir; i++)
		{
			m_pMxArray = matGetVariable(m_pmatFile,dir[i]);
			ReadObject( m_pMxArray );
		}
	}
	else if ( boost::iequals( "UF_Index", dir[0] ))
	{
		for (int i = 0; i < ndir; i++)
		{
			m_pMxArray = matGetVariable(m_pmatFile,dir[i]);
			ReadMsg( m_pMxArray );
		}
	}


#if 0
	m_pMxArray = matGetVariable(m_pmatFile,dir[0]);				//会增加内存占用
	nfiled = mxGetNumberOfFields(m_pMxArray);
	if ( 0 == nfiled )
	{
		std::cerr<<"mxArray下无结构体..."<<std::endl;
		return ImpMatError_t::ARRAYNULL;
	}
	std::vector<mxArray*> dataptr;
	for (size_t _index = 0;_index<nfiled;++_index)
	{
		mxArray *innerArray = mxGetFieldByNumber_730(m_pMxArray,0,_index);
		if (mxIsSparse(innerArray))
		{
			std::cerr<<"有稀疏矩阵"<<std::endl;
			ReadSparseMat(innerArray);
		}

		dataptr.push_back(innerArray);
	}
#endif
	return ImpMatError_t::IMPSUCCESS;
}

//该函数用于读入Index.mat，读入所有矩阵的信息
//@_searchpath：索引路径
//@_matname：默认名称为UF_Index.mat
// ImpMatError_t	CImpFrmMatlab::ImpFrmMatlab( boost::filesystem::path _searchpath, const char* _matname  )
// {
// 	boost::optional<boost::filesystem::path> r = find_file(_searchpath, _matname);
// 	if (!r)
// 	{
// 		std::cerr<<"boost搜索当前路径"<<_searchpath<<"下无该文件"<<std::endl;
// 		return ImpMatError_t::MATOPENFAILD;
// 	}
// 	else
// 	{
// 		std::cout<<"当前文件所在路径"<<*r<<std::endl;
// 	}
// 	int ndir;	//根目录下mxarray数目
// 	int nfiled = -10;	//mxarray下file数目  
// 	m_pmatFile = matOpen( r->string().c_str(),"r" );
// 	if (nullptr == m_pmatFile)
// 	{
// 		std::cerr<<"当前路径索引不正确未找到mat文件"<<std::endl;
// 		return ImpMatError_t::MATGETDIRERROR;
// 	}
// 
// 
// }


//输出CSR形式的系数矩阵模板函数，默认参数double
HBXDef::_CSRInput<double>& CImpFrmMatlab::GetStiffMat( bool _bsave, boost::filesystem::path _SavePath )
{
	using namespace boost;
	//备案存档，用于读入数据的程序校验
	if ( _bsave )
	{
		using namespace boost::gregorian;
		using namespace boost::posix_time;

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
		HBXDef::_ExportCSRdataToFile( _vNoneZeroVal.data(), _viColSort.data(), _viNonZeroRowSort.data(),
			m_MatMsg._nnzA, m_RowNum, "NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data",
			m_SavePath );
	}
	return m_MatMsg;
}

//输出向量形式的右端向量,因为库中有些不存在右端向量，故分开导出
double*	CImpFrmMatlab::GetRhsVec( bool _bsave )
{
	using namespace boost;
	if (!m_MatMsg.h_rhs)
	{
		std::cout<<"mat中无等式右端向量，是否需要生成？..."<<std::endl;
		char ch;
		there:
		std::cout<<"输入Y/N表示生成或不生成右端向量..."<<std::endl;
		ch = getchar();
		if ( !isalpha(ch) )
		{
			std::cout<<"请输入字母..."<<std::endl;
			goto there;
		}
		else
		{
			switch (ch)
			{
			case 'y':case 'Y':
				genRhs( g_MatRowNum );
				break;
			case  'n':case 'N':
				goto end;
				break;
			default:
				std::cerr<<"输入错误，跳出"<<std::endl;
				goto end;
				break;
			}
		}
	}
	if (_bsave && m_MatMsg.h_rhs)
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
		HBXDef::OutputPointData( m_MatMsg.h_rhs, m_RowNum, "LoadVec.data", m_SavePath );
	}
	end:
	return m_MatMsg.h_rhs;
}

//输出nnA的长度
int& CImpFrmMatlab::GetNoneZeroLgt()
{
	return m_nnzA;
}

//输出rowsort的长度，即rownum+1
int& CImpFrmMatlab::GetnALgt()
{
	return m_nA;
}

//用于读取单个矩阵的相关信息
void CImpFrmMatlab::ReadObject(const mxArray *const _Array)
{
	int nobj = mxGetNumberOfFields(_Array);
	if (0 != nobj)
	{
		for (int i = 0; i < nobj; i++)
		{
			mxChar *_tmpchar;
			double *_tmpd;
			mxArray *innerArray = mxGetFieldByNumber_730(_Array,0,i);
			const char *_fieldName = mxGetFieldNameByNumber(_Array, i);
			mxClassID _tmpid = mxGetClassID(innerArray);
			
			switch (_tmpid)
			{
			case mxCHAR_CLASS:
				_tmpchar = mxGetChars(innerArray);
				break;
			case mxDOUBLE_CLASS:
				if (mxIsSparse(innerArray))		//获取稀疏矩阵A
				{
					std::cout<<"有稀疏矩阵"<<std::endl;
					_mxArray_Map.insert( std::make_pair(_fieldName, innerArray) );		
					ReadSparseMat(innerArray);
				}
				else if( 0 == std::strcmp(_fieldName, "b") )	//获取右端向量
				{
					std::cout<<"有等式右端向量"<<std::endl;
					_mxArray_Map.insert( std::make_pair(_fieldName, innerArray) );
					ReadRhs(innerArray);
				}
				else if (0 == std::strcmp(_fieldName, "id") )	//获取ID号
				{
					_mxArray_Map.insert( std::make_pair(_fieldName, innerArray) );
					m_id = *mxGetPr(innerArray);
				}
				_tmpd = mxGetPr(innerArray);
				break;
			case mxSTRUCT_CLASS:	//aux的情形
				//待定
				break;
			case mxOBJECT_CLASS:
				break;
			default:
				break;
			}

			ReadObject(innerArray);
		}
	}

}

//用于读取UF_Index内的信息
void CImpFrmMatlab::ReadMsg( const mxArray* const _Array )
{
	using namespace boost;
	int nobj = mxGetNumberOfFields(_Array);
	if (0 != nobj)
	{
		for (int i = 0; i < nobj; i++)
		{
			mxChar *_tmpchar;
			double* _tmpd, _tmpTimeStamp;
			mxArray *innerArray = mxGetFieldByNumber_730(_Array,0,i);
			const char *_fieldName = mxGetFieldNameByNumber(_Array, i);
			mxClassID _tmpid = mxGetClassID(innerArray);

			if ( iequals( "LastRevisionDate", _fieldName ) )
			{
				_tmpchar = mxGetChars(innerArray);
			}
			else if ( iequals("DownloadTimeStamp", _fieldName) )
			{
				_tmpTimeStamp = *mxGetPr(innerArray);
			}
			else if ( iequals("group", _fieldName) )
			{
				ReadCell(innerArray);
			}


		}
	}
}

void CImpFrmMatlab::ReadSparseMat(const mxArray *const _dataptr)
{
	m_RowNum = mxGetM( _dataptr );
	m_ColNum = mxGetN( _dataptr );
	if ( m_RowNum == m_ColNum )
	{
		std::cout<<"矩阵为square，对全局变量g_MatRowNum赋值..."<<std::endl;
		std::cout<<"当前矩阵维度为"<<m_RowNum<<std::endl;
		//以下确定维度
		g_MatRowNum = (int)m_RowNum;
		m_MatMsg._nA = m_nA = (int)(m_RowNum)+1;
		m_MatMsg._nnzA = m_nnzA = (int)mxGetNzmax( _dataptr )-1;
	}
	m_MatMsg.h_NoneZeroVal = mxGetPr( _dataptr );
	m_MatMsg.h_iColSort = mxGetIr( _dataptr ); //rowsort = nna_z
	m_MatMsg.h_iNonZeroRowSort = mxGetJc( _dataptr ); //colsort = N+1

	SetStiffMat( m_MatMsg.h_NoneZeroVal, m_MatMsg.h_iColSort, m_MatMsg.h_iNonZeroRowSort, m_nnzA, m_nA );
}

//从外部获取三个数组指针得到刚度矩阵的CSR格式
void CImpFrmMatlab::SetStiffMat( const double* _srcVal, const size_t* _srcCol, size_t* _srcRow, size_t _nnA, size_t _nA )
{
	_vNoneZeroVal.resize(_nnA);
	_viColSort.resize(_nnA);
	_viNonZeroRowSort.resize(_nA);

	//在此需每个元素赋值，因为vector的内存可能不连续
	for (size_t i = 0; i < m_nnzA; i++)
	{
		_vNoneZeroVal[i] = _srcVal[i];
		_viColSort[i] = (int)_srcCol[i];
	}
	for (size_t i = 0; i < m_nA; i++)
	{
		_viNonZeroRowSort[i] = (int)_srcRow[i];
	}
	if (_viNonZeroRowSort.back() > _viColSort.size())
	{
		_viNonZeroRowSort.back()--;
		_srcRow[m_RowNum]--;
		std::cout<<"输入的CSR格式的系数矩阵首行索引向量最后一项有误需检查"<<std::endl;
	}
}

void CImpFrmMatlab::ReadRhs( const mxArray *const _Array )
{
	m_MatMsg.h_rhs = mxGetPr( _Array );
	_vRhsY.resize(g_MatRowNum);
	//在此需每个元素赋值，因为vector的内存可能不连续
	for (int i = 0; i < g_MatRowNum; i++)
	{
		_vRhsY[i] = m_MatMsg.h_rhs[i];
	}
}

void CImpFrmMatlab::ReadCell( const mxArray* const _Array )
{
	int cellM, cellN;
	cellM = (int)mxGetM(_Array);
	cellN = (int)mxGetN(_Array);
	for (int idx=0; idx<cellM*cellN; idx++)
	{
		mxArray*	mat1 = mxGetCell(_Array, idx);
		std::string _tmpstr = mxArrayToString(mat1);//指向mat1的第一个数据


	}
}

//生成随机的等式右端向量
void CImpFrmMatlab::genRhs( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath )
{
	using namespace std;
	srand((unsigned)time(NULL));
	cout<<"生成等式右边，原始数据在内存中抹去..."<<endl;
	m_SavePath = _savepath;
	_vRhsY.resize(g_MatRowNum);
	for (int i = 0; i < _genRowNum; i++)
	{
		_vRhsY[i] = 10*(double)rand()/RAND_MAX+1.0f;
	}
	m_MatMsg.h_rhs = _vRhsY.data();
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
		HBXDef::OutputPointData( m_MatMsg.h_rhs, m_RowNum, "LoadVec.data", m_SavePath );
	}
}