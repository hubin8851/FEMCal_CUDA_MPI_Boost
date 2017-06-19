#pragma once
#include "stdafx.h"
#include "CuComponent.cuh"
#include "AeroTable.h"

//避免命名空间的污染
namespace HBXDef
{
#ifndef  TPOW2(_VAL)
#define  TPOW2(_VAL) TPOW<2, _VAL>::VAL
#endif

	template< unsigned int T, HBXDef::CudaMalloc_t M >
	class CBaseLag:
		public CuComponent<M>
	{
	private:
		typedef CBaseLag<T, M> _SameBaseLag;
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
		__host__ CBaseLag( HBXDef::_AEROTABLE* _IptTable, size_t _blkId  )
		{
			{
				Assert( T>0 )
				IFNULL( _IptTable, "导入气动数据表错误...");
				m_isInit = false;
				m_table_p = _IptTable;
				this->SetBlkId( _blkId );

				memset( m_lag_cordinate, 0, T * sizeof(unsigned int) );
			}
		}


		__host__ CBaseLag( const _SameBaseLag& _rhs )
		{
			Assert( T>0 )
				memcpy( m_numperdim, _rhs.m_numperdim, T * sizeof(unsigned int) );

			memset( m_lag_cordinate, 0, T * sizeof(unsigned int) );

			m_cordials = _rhs.m_cordials;
			m_data = _rhs.m_data;

			m_isInit = true;
		}

		~CBaseLag()
		{
			if (m_table_p)
			{
				m_table_p = nullptr;
				delete	m_table_p;
			}

			if (m_cordials)
			{
				delete	[]m_cordials;
				m_cordials = nullptr;
			}
		}

		//@__blkId:待索引的block的ID,以0为起始索引
		//@UserStatusError_t:返回状态，可用CheckUserDefErrors宏判断
		UserStatusError_t SetBlkId( size_t _blkId )
		{
			if (_blkId >= m_table_p->_blocknum)
			{
				return UserStatusError_t::USER_STATUS_INVALID_VALUE;
			}
			m_slct_block_num = _blkId;
			m_dim = m_table_p->_blocks[_blkId]._dim;
			if ( m_dim <= 0 || T != m_dim )
			{
				return UserStatusError_t::USER_STATUS_INVALID_VALUE;
			}

			return UserStatusError_t::USER_STATUS_SUCCESS;
		}


		//主要是把结构体_AeroTable中的二维数组转化为一维
		bool Initial()
		{
			memcpy( m_numperdim, m_table_p->_blocks[m_slct_block_num]._numperdim, sizeof(unsigned int)*T );

			unsigned int _cordlgth = 0;	//插值节点数组的长度
			for (int i = 0; i < m_dim; i++)
			{
				_cordlgth += m_numperdim[i];
			}
			m_cordials = new HBXDef::UserDefFloat[_cordlgth];

			_cordlgth = 0;	//在此做二维转一维数组的临时变量
			for (int i = 0; i < m_dim; i++)
			{
				for (int j = 0; j < m_numperdim[i]; j++)
				{
					m_cordials[_cordlgth] = m_table_p->_blocks[m_slct_block_num]._corddata[i][j];
					_cordlgth++;
				}
			}

			m_data = m_table_p->_blocks[m_slct_block_num]._data;

			m_isInit = true;
			return UserStatusError_t::USER_STATUS_SUCCESS;
		}

		HBXDef::UserDefFloat	ReadTableData( HBXDef::UserDefFloat* _ArrayIn )
		{
			unsigned int _tmp_position[ TPOW2(T) ];	//递归临时数组
			HBXDef::UserDefFloat	_tmp_reducdata[ TPOW2(T) ];	//递归数据临时数组
			HBXDef::UserDefFloat	_tmp_proper[ T ];	//递归比例系数

			memset((void*)&_tmp_position, 0, TPOW2(T) * sizeof(unsigned int) );

			findcordi( _tmp_position, _tmp_proper, _ArrayIn );
			//迭代递归算法
			unsigned int _tmpidx = TPOW2(T)/2;
			for (unsigned int idx = 0; idx < TPOW2(T); idx++)//在临时数组中
			{
				_tmp_reducdata[idx] = m_data[ _tmp_position[idx] ];
			}
			for (unsigned int idx = T-1; idx > 0 ; idx--)
			{
				for (int idx2 = 0; idx2 < _tmpidx; idx2++)
				{//循环内表达式简而言之化为b（i） = a（i+offset）+coef*(a(i+offset)-a(i)),两点插值
					_tmp_reducdata[idx2] = _tmp_reducdata[idx2 + _tmpidx] - _tmp_proper[idx] * (_tmp_reducdata[idx2 + _tmpidx] - _tmp_reducdata[idx2]);
				}
				_tmpidx >> 1;//通过左移操作对另一维度的数据进行插值
			}//在此没有将最后一步放入循环中，便于未优化情况下提高速度
			return _tmp_reducdata[1] - _tmp_proper[0] * (_tmp_reducdata[1] - _tmp_reducdata[0]);
		}

		//寻找待插值点所围成空间所构成的点
		__host__ void	findcordi( unsigned int _tmp_position[TPOW2(T)], HBXDef::UserDefFloat _tmp_proper[T], HBXDef::UserDefFloat _vecIn[T] )
		{
			{
				unsigned int tempmul = 1;		//每个维度插值区间上点个数的乘积
				unsigned int temppow = 1;		//用于控制迭代过程中vec每个维度的偏移量
				unsigned int _tmpoffset = 0;	//各维度起始坐标在一维数组中的索引
				unsigned int		_length = 0;
				HBXDef::UserDefFloat*	_cordials = nullptr;
				unsigned int		_cordial = 0;

				for (unsigned int idx1 = 0; idx1 < T; idx1++)
				{
					_length = m_numperdim[idx1];	//寻找当前维度下每个维度的节点数目

					_cordials = &m_cordials[_tmpoffset];	
					_cordial = m_lag_cordinate[idx1];

					//判断当前插值点在上一步的左方亦或右方
					if (_vecIn[idx1] > _cordials[_cordial])
					{
						for (int idx2 = _cordial; idx2 < _length; idx2++)
						{
							if ((_length - 1) == idx2)	//判断是否在右边界
							{
								_tmp_proper[idx1] = 0;	//递归因子置0
								m_lag_cordinate[idx1] = _length - 1; 
							}
							else
							{
								if (_vecIn[idx1] < _cordials[idx2+1])
								{
									m_lag_cordinate[idx1] = idx2 + 1;
									_tmp_proper[idx1] = (_cordials[idx2 + 1] - _vecIn[idx1]) / (_cordials[idx2 + 1] - _cordials[idx2]);
									break;
								}
							}
						}
					}
					else
					{//左方
						for (int idx2 = _cordial; idx2 >= 0; idx2--)
						{
							if (idx2 == 0)//若该插值点靠近左端边界
							{
								_tmp_proper[idx1] = 1;				//递归因数
								m_lag_cordinate[idx1] = 1;
								break;
							}
							else
							{
								if ( _vecIn[idx1] > _cordials[idx2-1] )	//若在_cordials[idx2-1]和_cordials[idx2]之间
								{
									m_lag_cordinate[idx1] = idx2;
									_tmp_proper[idx1] = (_cordials[idx2] - _vecIn[idx1]) / (_cordials[idx2] - _cordials[idx2-1]);
									break;
								}
							}
						}
					}
					for (int idx2 = 0; idx2 < TPOW2(T); idx2++)
					{
						if (0 != (idx2 / temppow) % 2)
							_tmp_position[idx2] += m_lag_cordinate[idx1] * tempmul;
						else
							_tmp_position[idx2] += (m_lag_cordinate[idx1]-1) * tempmul;
					}
					tempmul *= _length;
					_tmpoffset += _length;
					temppow = temppow << 1;
				}
			}
		}

		//返回当前Table的总维度
		__host__ size_t GetDemtion(){ return T; };
	};


}


// .inl文件