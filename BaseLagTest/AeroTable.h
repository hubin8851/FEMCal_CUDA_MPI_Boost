#pragma once
#include "CuComponent.cuh"

namespace HBXDef
{
	typedef struct _AeroBlock:public _CuBase_uni_
	{
		char _name[MAXPATH_LGT];	//当前块名称
		unsigned int		_dim;				//插值维度
		unsigned int*		_numperdim;			//每个维度的坐标个数
		HBXDef::UserDefFloat**		_corddata;			//每个维度的坐标数值
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
