#pragma once

#include "StdAfx.h"
#include "CuComponent.cuh"

//避免命名空间的污染
namespace HBXDef
{
#ifndef  TPOW2(_VAL)
#define  TPOW2(_VAL) TPOW<2, _VAL>::VAL
#endif

	template< unsigned int T >
	class CBaseLag
	{
	private:
		typedef CBaseLag<T> _SameBaseLag;
		enum 
		{
			VALUE = T
		};
		bool	m_isInit;
	protected:
		size_t	m_slct_block_num;	//所选的block
		size_t	m_dim;
		HBXDef::_AEROTABLE* m_table_p;
	public:
		unsigned int	m_numperdim[T];	//栈上分配,一个block内每个维度上的数值点个数
		unsigned int	m_lag_cordinate[T];//?
		HBXDef::UserDefFloat	m_lag[T];//?
		//
		HBXDef::UserDefFloat*	m_cordials;	//线性化后的每个维度上的插值节点数组
		HBXDef::UserDefFloat*	m_data;
	public:
		__host__ CBaseLag( HBXDef::_AEROTABLE* _IptTable, size_t _blkId  );
		__host__ CBaseLag( const _SameBaseLag& _rhs );
		~CBaseLag();

		//@__blkId:待索引的block的ID,以0为起始索引
		//@UserStatusError_t:返回状态，可用CheckUserDefErrors宏判断
		UserStatusError_t SetBlkId( size_t _blkId );

		//主要是把结构体_AeroTable中的二维数组转化为一维
		bool Initial();

		HBXDef::UserDefFloat	ReadTableData( HBXDef::UserDefFloat* _ArrayIn );

		//寻找待插值点所围成空间所构成的点
		void	findcordi( unsigned int _tmp_position[TPOW2(T)], HBXDef::UserDefFloat _tmp_proper[T], HBXDef::UserDefFloat _vecIn[T] );

		//返回当前Table的总维度
		__host__ size_t GetDemtion(){ return T; };

	};

}



#include "BaseLag.inl"
// .inl文件