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

//�ú�������ĳһ·���µĺ�����������Ŀ¼������
//@_matname:��������
//@_searchpath������·��
ImpMatError_t CImpFrmMatlab::ImpFrmMatlab( const char* _matname,  boost::filesystem::path _searchpath)
{

	boost::optional<boost::filesystem::path> r = find_file(_searchpath,_matname);
	if (!r)
	{
		std::cerr<<"boost������ǰ·��"<<_searchpath<<"���޸��ļ�"<<std::endl;
		return ImpMatError_t::MATOPENFAILD;
	}
	else
	{
		std::cout<<"��ǰ�ļ�����·��"<<*r<<std::endl;
	}
	int ndir;	//��Ŀ¼��mxarray��Ŀ
	int nfiled = -10;	//mxarray��file��Ŀ  
	m_pmatFile = matOpen( r->string().c_str(),"r" );
	if (nullptr == m_pmatFile)
	{
		std::cerr<<"��ǰ·����������ȷδ�ҵ�mat�ļ�"<<std::endl;
		return ImpMatError_t::MATGETDIRERROR;
	}
	const char** dir = (const char **)matGetDir(m_pmatFile, &ndir);
	if ( nullptr == dir )
	{
		if ( 0 == ndir )
		{
			std::cerr<<"��ǰ·��"<<_searchpath<<"���޸��ļ�"<<std::endl;
			return ImpMatError_t::MATFILENULL;
		}
		else if ( 0 > ndir)
		{
			std::cerr<<"��ȡmat�ļ�ʱ��������..."<<std::endl;
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
	m_pMxArray = matGetVariable(m_pmatFile,dir[0]);				//�������ڴ�ռ��
	nfiled = mxGetNumberOfFields(m_pMxArray);
	if ( 0 == nfiled )
	{
		std::cerr<<"mxArray���޽ṹ��..."<<std::endl;
		return ImpMatError_t::ARRAYNULL;
	}
	std::vector<mxArray*> dataptr;
	for (size_t _index = 0;_index<nfiled;++_index)
	{
		mxArray *innerArray = mxGetFieldByNumber_730(m_pMxArray,0,_index);
		if (mxIsSparse(innerArray))
		{
			std::cerr<<"��ϡ�����"<<std::endl;
			ReadSparseMat(innerArray);
		}

		dataptr.push_back(innerArray);
	}
#endif
	return ImpMatError_t::IMPSUCCESS;
}

//�ú������ڶ���Index.mat���������о������Ϣ
//@_searchpath������·��
//@_matname��Ĭ������ΪUF_Index.mat
// ImpMatError_t	CImpFrmMatlab::ImpFrmMatlab( boost::filesystem::path _searchpath, const char* _matname  )
// {
// 	boost::optional<boost::filesystem::path> r = find_file(_searchpath, _matname);
// 	if (!r)
// 	{
// 		std::cerr<<"boost������ǰ·��"<<_searchpath<<"���޸��ļ�"<<std::endl;
// 		return ImpMatError_t::MATOPENFAILD;
// 	}
// 	else
// 	{
// 		std::cout<<"��ǰ�ļ�����·��"<<*r<<std::endl;
// 	}
// 	int ndir;	//��Ŀ¼��mxarray��Ŀ
// 	int nfiled = -10;	//mxarray��file��Ŀ  
// 	m_pmatFile = matOpen( r->string().c_str(),"r" );
// 	if (nullptr == m_pmatFile)
// 	{
// 		std::cerr<<"��ǰ·����������ȷδ�ҵ�mat�ļ�"<<std::endl;
// 		return ImpMatError_t::MATGETDIRERROR;
// 	}
// 
// 
// }


//���CSR��ʽ��ϵ������ģ�庯����Ĭ�ϲ���double
HBXDef::_CSRInput<double>& CImpFrmMatlab::GetStiffMat( bool _bsave, boost::filesystem::path _SavePath )
{
	using namespace boost;
	//�����浵�����ڶ������ݵĳ���У��
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

//���������ʽ���Ҷ�����,��Ϊ������Щ�������Ҷ��������ʷֿ�����
double*	CImpFrmMatlab::GetRhsVec( bool _bsave )
{
	using namespace boost;
	if (!m_MatMsg.h_rhs)
	{
		std::cout<<"mat���޵�ʽ�Ҷ��������Ƿ���Ҫ���ɣ�..."<<std::endl;
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
				genRhs( g_MatRowNum );
				break;
			case  'n':case 'N':
				goto end;
				break;
			default:
				std::cerr<<"�����������"<<std::endl;
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

//���nnA�ĳ���
int& CImpFrmMatlab::GetNoneZeroLgt()
{
	return m_nnzA;
}

//���rowsort�ĳ��ȣ���rownum+1
int& CImpFrmMatlab::GetnALgt()
{
	return m_nA;
}

//���ڶ�ȡ��������������Ϣ
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
				if (mxIsSparse(innerArray))		//��ȡϡ�����A
				{
					std::cout<<"��ϡ�����"<<std::endl;
					_mxArray_Map.insert( std::make_pair(_fieldName, innerArray) );		
					ReadSparseMat(innerArray);
				}
				else if( 0 == std::strcmp(_fieldName, "b") )	//��ȡ�Ҷ�����
				{
					std::cout<<"�е�ʽ�Ҷ�����"<<std::endl;
					_mxArray_Map.insert( std::make_pair(_fieldName, innerArray) );
					ReadRhs(innerArray);
				}
				else if (0 == std::strcmp(_fieldName, "id") )	//��ȡID��
				{
					_mxArray_Map.insert( std::make_pair(_fieldName, innerArray) );
					m_id = *mxGetPr(innerArray);
				}
				_tmpd = mxGetPr(innerArray);
				break;
			case mxSTRUCT_CLASS:	//aux������
				//����
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

//���ڶ�ȡUF_Index�ڵ���Ϣ
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
		std::cout<<"����Ϊsquare����ȫ�ֱ���g_MatRowNum��ֵ..."<<std::endl;
		std::cout<<"��ǰ����ά��Ϊ"<<m_RowNum<<std::endl;
		//����ȷ��ά��
		g_MatRowNum = (int)m_RowNum;
		m_MatMsg._nA = m_nA = (int)(m_RowNum)+1;
		m_MatMsg._nnzA = m_nnzA = (int)mxGetNzmax( _dataptr )-1;
	}
	m_MatMsg.h_NoneZeroVal = mxGetPr( _dataptr );
	m_MatMsg.h_iColSort = mxGetIr( _dataptr ); //rowsort = nna_z
	m_MatMsg.h_iNonZeroRowSort = mxGetJc( _dataptr ); //colsort = N+1

	SetStiffMat( m_MatMsg.h_NoneZeroVal, m_MatMsg.h_iColSort, m_MatMsg.h_iNonZeroRowSort, m_nnzA, m_nA );
}

//���ⲿ��ȡ��������ָ��õ��նȾ����CSR��ʽ
void CImpFrmMatlab::SetStiffMat( const double* _srcVal, const size_t* _srcCol, size_t* _srcRow, size_t _nnA, size_t _nA )
{
	_vNoneZeroVal.resize(_nnA);
	_viColSort.resize(_nnA);
	_viNonZeroRowSort.resize(_nA);

	//�ڴ���ÿ��Ԫ�ظ�ֵ����Ϊvector���ڴ���ܲ�����
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
		std::cout<<"�����CSR��ʽ��ϵ���������������������һ����������"<<std::endl;
	}
}

void CImpFrmMatlab::ReadRhs( const mxArray *const _Array )
{
	m_MatMsg.h_rhs = mxGetPr( _Array );
	_vRhsY.resize(g_MatRowNum);
	//�ڴ���ÿ��Ԫ�ظ�ֵ����Ϊvector���ڴ���ܲ�����
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
		std::string _tmpstr = mxArrayToString(mat1);//ָ��mat1�ĵ�һ������


	}
}

//��������ĵ�ʽ�Ҷ�����
void CImpFrmMatlab::genRhs( size_t _genRowNum, bool _bsave, boost::filesystem::path _savepath )
{
	using namespace std;
	srand((unsigned)time(NULL));
	cout<<"���ɵ�ʽ�ұߣ�ԭʼ�������ڴ���Ĩȥ..."<<endl;
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