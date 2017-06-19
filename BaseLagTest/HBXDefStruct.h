#ifndef _HBXDEFSTRUCT
#define _HBXDEFSTRUCT
#include "stdafx.h"
#include "CuComponent.cuh"


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



#pragma endregion


}



#endif