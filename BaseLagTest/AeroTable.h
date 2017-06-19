#pragma once
#include "CuComponent.cuh"

namespace HBXDef
{
	typedef struct _AeroBlock:public _CuBase_uni_
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

	typedef struct _AeroTable:public _CuBase_uni_
	{
	public:
		char _path[MAXPATH_LGT];
		char _name[MAXPATH_LGT];
		size_t	_blocknum;
		_AeroBlock*	_blocks;
	public:
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

}
