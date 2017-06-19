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
		size_t	m_slct_block_num;	//��ѡ��block
		size_t	m_dim;
	protected:
		UserDefFloat*	m_cordials;	//ÿ��ά���ϵ�����̶�
		UserDefFloat*	m_data;		//����ֵ������
		unsigned int	m_numperdim[T];		//��ֵ����ά��
		unsigned int	m_lag_cordinate[T];	//��ֵ�����λ��
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

		//�ù��캯�������ݲ�ֵ���ָ�뼰�ñ��ڵı��
		__host__ CuInterpolation( HBXDef::_AEROTABLE* _IptTable, size_t _blkId  )
		{
			IFNULL( _IptTable, "�����������ݱ����...");
			m_table_p = _IptTable;
			this->SetBlkId( _blkId );
		}

		//�ù��캯�������ݲ�ֵ���������ݣ������������ֵ����������Ϣ��
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

		//@__blkId:��������block��ID,��0Ϊ��ʼ����
		//@UserStatusError_t:����״̬������CheckUserDefErrors���ж�
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

		//��Ҫ��ȷ������ֵ�����Ŀ����������ɴ���ֵ�������������Դ��ϵķ���
		//@_ArrayLgth������ֵ�����Ŀ
// 		__host__ bool Initial( unsigned int _ArrayLgth )
// 		{
// 			m_isInit = false;
// 			std::cout<<sizeof(_USERDEFT)<<std::endl;
// 			checkCudaErrors( cudaMalloc(&p_UserDefT, _ArrayLgth * sizeof(_USERDEFT) ) );
// 
// 			m_isInit = true;
// 			return true;
// 		}

		//��BaseLag�л�ȡ���ݲ���ʼ��
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

		//��ͬһ��ֵ�������ɸ�����в�ֵ����
		//__global__ void  ReadTableArray( _USERDEFT* _ArrayIn );

		//�ú�������CUDA�����ڵ�global�����ڵ���
		//@_ArrayIn:����Ĵ���ֵ����������
		__device__ HBXDef::UserDefFloat	ReadTableData( HBXDef::UserDefFloat* _ArrayIn )
		{
			unsigned int _tmp_position[ TPOW2(T) ];	//�ݹ���ʱ����
			HBXDef::UserDefFloat	_tmp_reducdata[ TPOW2(T) ];	//�ݹ�������ʱ����
			HBXDef::UserDefFloat	_tmp_proper[ T ];	//�ݹ����ϵ��

			memset((void*)&_tmp_position, 0, TPOW2(T) * sizeof(unsigned int) );

			findcordi( _tmp_position, _tmp_proper, _ArrayIn );
			//�����ݹ��㷨
			unsigned int _tmpidx = TPOW2(T)/2;
			for (unsigned int idx = 0; idx < TPOW2(T); idx++)//����ʱ������
			{
				_tmp_reducdata[idx] = m_data[ _tmp_position[idx] ];
			}
			for (unsigned int idx = T-1; idx > 0 ; idx--)
			{
				for (int idx2 = 0; idx2 < _tmpidx; idx2++)
				{//ѭ���ڱ��ʽ�����֮��Ϊb��i�� = a��i+offset��+coef*(a(i+offset)-a(i)),�����ֵ
					_tmp_reducdata[idx2] = _tmp_reducdata[idx2 + _tmpidx] - _tmp_proper[idx] * (_tmp_reducdata[idx2 + _tmpidx] - _tmp_reducdata[idx2]);
				}
				_tmpidx = _tmpidx >> 1;//ͨ�����Ʋ�������һά�ȵ����ݽ��в�ֵ
			}//�ڴ�û�н����һ������ѭ���У�����δ�Ż����������ٶ�
			return _tmp_reducdata[1] - _tmp_proper[0] * (_tmp_reducdata[1] - _tmp_reducdata[0]);

		}

		//Ѱ�Ҵ���ֵ����Χ�ɿռ������ɵĵ�
		__device__ HBXDef::UserDefFloat	findcordi( unsigned int _tmp_position[TPOW2(T)], HBXDef::UserDefFloat _tmp_proper[T], HBXDef::UserDefFloat _vecIn[T] )
		{
			{
				unsigned int tempmul = 1;		//ÿ��ά�Ȳ�ֵ�����ϵ�����ĳ˻�
				unsigned int temppow = 1;		//���ڿ��Ƶ���������vecÿ��ά�ȵ�ƫ����
				unsigned int _tmpoffset = 0;	//��ά����ʼ������һά�����е�����
				unsigned int		_length = 0;
				HBXDef::UserDefFloat*	_cordials = nullptr;
				unsigned int		_cordial = 0;

				for (unsigned int idx1 = 0; idx1 < T; idx1++)
				{
					_length = m_numperdim[idx1];	//Ѱ�ҵ�ǰά����ÿ��ά�ȵĽڵ���Ŀ

					_cordials = &m_cordials[_tmpoffset];	
					_cordial = m_lag_cordinate[idx1];

					//�жϵ�ǰ��ֵ������һ����������ҷ�
					if (_vecIn[idx1] > _cordials[_cordial])
					{
						for (int idx2 = _cordial; idx2 < _length; idx2++)
						{
							if ((_length - 1) == idx2)	//�ж��Ƿ����ұ߽�
							{
								_tmp_proper[idx1] = 0;	//�ݹ�������0
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
					{//��
						for (int idx2 = _cordial; idx2 >= 0; idx2--)
						{
							if (idx2 == 0)//���ò�ֵ�㿿����˱߽�
							{
								_tmp_proper[idx1] = 1;				//�ݹ�����
								m_lag_cordinate[idx1] = 1;
								break;
							}
							else
							{
								if ( _vecIn[idx1] > _cordials[idx2-1] )	//����_cordials[idx2-1]��_cordials[idx2]֮��
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
#ifdef USERDEFSHARED//Ϊ�˿��ƹ����ڴ��С

#else // USERDEFSHARED//Ϊ�˿��ƹ����ڴ��С
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



