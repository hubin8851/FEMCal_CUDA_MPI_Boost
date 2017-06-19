#pragma once
#include "stdafx.h"

#include "BaseLag.hpp"
#include "CuComponent.cuh"
#include "AeroTable.h"


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
		enum 
		{
			VALUE = T
		};
	private:
		bool	m_isInit;
		size_t	m_slct_block_num;	//所选的block
		size_t	m_dim;
	protected:
		UserDefFloat*	m_cordials;	//每个维度上的坐标刻度
		UserDefFloat*	m_data;		//待插值的数据
		unsigned int	m_numperdim[T];		//插值坐标维度
		unsigned int	m_lag_cordinate[T];	//插值坐标标位置
		_USERDEFT*	p_UserDefT;	

	public:
		void* operator new(size_t len) {
			void* ptr;
			cudaMallocManaged(&ptr, sizeof(CuInterpolation)*len);
			return ptr;
		};

		__host__ CuInterpolation()
		{
			m_isInit = false;
		}

		//该构造函数仅传递插值表的指针及该表内的表块
		__host__ CuInterpolation( HBXDef::_AEROTABLE* _IptTable, size_t _blkId  )
		{
			IFNULL( _IptTable, "导入气动数据表错误...");
			m_table_p = _IptTable;
			this->SetBlkId( _blkId );
		}

		//该构造函数仅传递插值表的相关数据，并不传输待插值点向量等信息。
		__host__ CuInterpolation( const HBXDef::CBaseLag<T, M>& _rhs )
		{
			switch (M)
			{
			case NORMAL: 
			case PAGELOCK: 
			case UNIFIEDMEM:
				checkCudaErrors( cudaMemcpy(m_numperdim, _rhs.m_numperdim, T * sizeof(unsigned int), cudaMemcpyHostToDevice ) );
				checkCudaErrors( cudaMemcpy(m_lag_cordinate, _rhs.m_lag_cordinate, T * sizeof(unsigned int), cudaMemcpyHostToDevice) );
				checkCudaErrors( cudaMemcpy(&m_cordials, &_rhs.m_cordials, sizeof(void*), cudaMemcpyDeviceToDevice));
				checkCudaErrors( cudaMemcpy(&m_data, &_rhs.m_data, sizeof(void*), cudaMemcpyDeviceToDevice));
				break;
			case ARRAY2D:
				break;
			default:
				break;
			}

		}

		__host__ ~CuInterpolation(void)
		{

		}

		//@__blkId:待索引的block的ID,以0为起始索引
		//@UserStatusError_t:返回状态，可用CheckUserDefErrors宏判断
		__host__ UserStatusError_t SetBlkId( size_t _blkId )
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

		//主要是确定待插值点的数目，并完成若干待插值点坐标向量在显存上的分配
		//@_ArrayLgth：待插值点的数目
// 		__host__ bool Initial( unsigned int _ArrayLgth )
// 		{
// 			m_isInit = false;
// 			std::cout<<sizeof(_USERDEFT)<<std::endl;
// 			checkCudaErrors( cudaMalloc(&p_UserDefT, _ArrayLgth * sizeof(_USERDEFT) ) );
// 
// 			m_isInit = true;
// 			return true;
// 		}

		//从BaseLag中获取数据并初始化
		__host__ bool Initial( const CBaseLag<T, M>& _rhs )
		{
// 			m_isInit = false;
// 			switch (M)
// 			{
// 			case NORMAL: 
// 			case PAGELOCK: 
// 			case UNIFIEDMEM:
				cuerror( cudaMemcpy(m_numperdim, _rhs.m_numperdim, T * sizeof(unsigned int), cudaMemcpyHostToDevice ) );
				cuerror( cudaMemcpy(m_lag_cordinate, _rhs.m_lag_cordinate, T * sizeof(unsigned int), cudaMemcpyHostToDevice) );
				int _mult, _add;
				for (int i = 0; i < this->VALUE; i++)
				{
					_add += m_numperdim[i];
					_mult *= m_numperdim[i];
				}
				cuerror( cudaMemcpy(&m_cordials, &_rhs.m_cordials, sizeof(void*), cudaMemcpyDeviceToDevice) );
				cuerror( cudaMemcpy(&m_data, &_rhs.m_data, sizeof(void*), cudaMemcpyDeviceToDevice) );
// 				break;
// 			case ARRAY2D:
// 				break;
// 			default:
//				break;
// 			}
//			m_isInit = true;
			return true;
		}

		//对同一插值表内若干个点进行插值计算
		//__global__ void  ReadTableArray( _USERDEFT* _ArrayIn );

		//该函数置于CUDA环境内的global函数内调用
		//@_ArrayIn:传入的带插值的坐标向量
		__device__ HBXDef::UserDefFloat	ReadTableData( HBXDef::UserDefFloat* _ArrayIn )
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
				_tmpidx = _tmpidx >> 1;//通过左移操作对另一维度的数据进行插值
			}//在此没有将最后一步放入循环中，便于未优化情况下提高速度
			return _tmp_reducdata[1] - _tmp_proper[0] * (_tmp_reducdata[1] - _tmp_reducdata[0]);

		}

		//寻找待插值点所围成空间所构成的点
		__device__ HBXDef::UserDefFloat	findcordi( unsigned int _tmp_position[TPOW2(T)], HBXDef::UserDefFloat _tmp_proper[T], HBXDef::UserDefFloat _vecIn[T] )
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
	};

	template<unsigned int T, HBXDef::CudaMalloc_t M >
	__global__	void cutestlag1( HBXDef::CuInterpolation<T, M>* _ArrayIn, HBXDef::UserDefFloat* _outData )
	{
		using namespace HBXDef;
		unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

		while ( idx < MAX_GPUITER )
		{
#ifdef USERDEFSHARED//为了控制共享内存大小

#else // USERDEFSHARED//为了控制共享内存大小
#endif  //WITHTEST
			UserDefFloat intercord[3];
			for (size_t i = 0; i < RUNCLC; i++)
			{
				intercord[0] = -3.0 + i * 20.0 / RUNCLC + 20.0 * i / MAX_GPUITER;
				intercord[1] = 0.2 + i * 2.0 / RUNCLC + 2.0 * i / MAX_GPUITER;
				intercord[2] = -30 + i * 30.0 / RUNCLC + 30 * i / MAX_GPUITER;
//				intercord[3] = -5.0 + 5*i / RUNCLC + i*5 / MAX_GPUITER;

				_ArrayIn[idx].ReadTableData(intercord);
			}
			idx += gridDim.x*blockDim.x;
		}
		__syncthreads();
	}


}



