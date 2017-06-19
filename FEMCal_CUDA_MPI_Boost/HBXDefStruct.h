#ifndef _HBXDEFSTRUCT
#define _HBXDEFSTRUCT


extern int g_MatRowNum;

namespace HBXDef
{
	template< unsigned int _base, unsigned int _exponent>
	struct TPOW
	{//?С���ʣ���Ϊ�ݹ飬Ϊ��û����ʼ��ֵ����TPOW��1������val = 1  ��
		enum
		{
			VAL = TPOW<_base, _exponent - 1>::VAL * _base
		};
	};

	template< unsigned int _base>
	struct TPOW<_base, 0>
	{
		enum
		{
			VAL = 1
		};
	};


	template<class T>
	struct _CooFormat
	{
		T _val;
		int _i;
		int _j;
	};


	template<class T1>
	struct _CSRInput
	{
		_CSRInput():h_NoneZeroVal(nullptr),h_iColSort(nullptr),h_iNonZeroRowSort(nullptr),h_x(nullptr),h_rhs(nullptr)
		{

		}
		_CSRInput(const T1* _val, size_t* _col, size_t* _row, const T1 *_x = nullptr, const T1 *_b = nullptr):
 			h_NoneZeroVal(_val),h_iColSort(_col),h_iNonZeroRowSort(_row),
 			h_x(_x),h_rhs(_b)
		{
		
		}
		//������buffer
		T1*		h_NoneZeroVal;
		size_t*	h_iColSort;
		size_t*	h_iNonZeroRowSort;
		T1*		h_x;
		T1*		h_rhs;
		int		_nnzA;
		int		_nA;
	};

/************************************************************************/
/*                                    mysql��ؽṹ��                                                              */
/************************************************************************/
#pragma region
	typedef struct 
	{
		int _matsize;
		int _iters;
		float _tol;
		char *_eigfilename;
		bool _bSuccess;
	} CalRecord_t;
#pragma endregion

/************************************************************************/
/*                      matlab��mx���ݶ�ȡ���                           */
/************************************************************************/
#pragma region	matlab���ݶ�ȡ���
	typedef struct mtxFrmMatlab
	{
		std::string	LastRevisionDate;
		double	TimeStamp;
		std::string	Group;
		std::string	Name;
		size_t _nRows;
		size_t _nCols;
		size_t _nnz;
		size_t _nZero;
	} mtxFrmMatlab_t;
#pragma endregion

/************************************************************************/
/*                          XML�������ݱ����                             */
/************************************************************************/
#pragma region  XML������ϵ�����ݱ����
	typedef struct _AeroBlock
	{
		char _name[MAXPATH_LGT];	//��ǰ������
		unsigned int		_dim;				//��ֵά��
		unsigned int*		_numperdim;			//ÿ��ά�ȵ��������
		HBXDef::UserDefFloat**		_corddata;			//ÿ��ά�ȵ�������ֵ
		HBXDef::UserDefFloat*		_data;

		_AeroBlock()
		{
			_numperdim = nullptr;
			_corddata = nullptr;
		}

		~_AeroBlock()
		{
			if (nullptr != _numperdim)
			{
				delete []_numperdim;
				_numperdim = nullptr;
			}
			if (nullptr != _corddata)
			{
				for (size_t i=0; i<_dim; i++)
				{
					delete []_corddata[i];
				}
				delete []_corddata;
				_corddata = nullptr;
			}
			if (nullptr != _data)
			{
				delete []_data;
				_data = nullptr;
			}
		}
	} _AEROBLOCK;


	typedef struct _AeroTable
	{
		char _path[MAXPATH_LGT];
		char _name[MAXPATH_LGT];
		size_t	_blocknum;
		_AeroBlock*	_blocks;
		_AeroTable()
		{
			_blocks = nullptr;
		}
		~_AeroTable()
		{
			if (nullptr != _blocks)
			{
				delete[] _blocks;
				_blocks = nullptr;
			}
		}
	} _AEROTABLE;

//Ϊ�˳������չ�ԣ��ɴ�XML���ȡÿ�ֵ�Ԫ��Ӧ������
typedef std::map < std::string, boost::shared_ptr<HBXDef::_AeroTable> >	_AeroTableInMap;
typedef HBXDef::_AeroTableInMap::iterator								_AeroTableMapIter;

#pragma endregion


}



#endif