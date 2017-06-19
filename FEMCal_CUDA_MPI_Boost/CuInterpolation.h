#pragma once

#include "StdAfx.h"
#include "BaseLag.h"
#include "CuComponent.cuh"


namespace HBXDef
{
	struct _USERDEFT
	{
		HBXDef::UserDefFloat	_loc[4];
	};


	template< unsigned int T, HBXDef::CudaMalloc_t M >
	class CuInterpolation:
		public CuComponent<M>
	{
	private:
		bool	m_isInit;
		size_t	m_slct_block_num;	//所选的block
		size_t	m_dim;
	protected:
		UserDefFloat*	m_cordials;	//每个维度上的坐标刻度
		UserDefFloat*	m_data;		//待插值的数据
		unsigned int		m_numperdim[T];		//插值坐标维度

		HBXDef::_AEROTABLE* m_table_p;	//获取的_AeroTable的指针

		unsigned int		m_lag_cordinate[T];	//插值坐标标位置
		_USERDEFT*	p_UserDefT;	

	public:
		//该构造函数仅传递插值表的指针及该表内的表块
		__host__ CuInterpolation( HBXDef::_AEROTABLE* _IptTable, size_t _blkId  );

		//该构造函数仅传递插值表的相关数据，并不传输待插值点向量等信息。
		__host__ CuInterpolation( const HBXDef::CBaseLag<T>& _rhs );

		__host__ ~CuInterpolation(void);

		//@__blkId:待索引的block的ID,以0为起始索引
		//@UserStatusError_t:返回状态，可用CheckUserDefErrors宏判断
		__host__ UserStatusError_t SetBlkId( size_t _blkId );

		//主要是确定待插值点的数目，并完成若干待插值点坐标向量在显存上的分配
		//@_ArrayLgth：待插值点的数目
		__host__ bool Initial( unsigned int _ArrayLgth );

		//从BaseLag中获取数据并初始化
		__host__ bool Initial( const CBaseLag<T>& _rhs );

		//对同一插值表内若干个点进行插值计算
		//__global__ void  ReadTableArray( _USERDEFT* _ArrayIn );

		//该函数置于CUDA环境内的global函数内调用
		//@_ArrayIn:传入的带插值的坐标向量
		__device__ HBXDef::UserDefFloat	ReadTableData( HBXDef::UserDefFloat* _ArrayIn );

		//寻找待插值点所围成空间所构成的点
		__host__ __device__ HBXDef::UserDefFloat	findcordi( unsigned int _tmp_position[TPOW2(T)], HBXDef::UserDefFloat _tmp_proper[T], HBXDef::UserDefFloat _vecIn[T] );
	};

	template<unsigned int T, HBXDef::CudaMalloc_t M >
	__global__	void cutestlag1( HBXDef::CuInterpolation<T, M>* _ArrayIn, HBXDef::UserDefFloat* _outData ); 

}





#include "CuInterpolation.inl"