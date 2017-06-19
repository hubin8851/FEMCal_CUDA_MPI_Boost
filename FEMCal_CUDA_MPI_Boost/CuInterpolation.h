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
		size_t	m_slct_block_num;	//��ѡ��block
		size_t	m_dim;
	protected:
		UserDefFloat*	m_cordials;	//ÿ��ά���ϵ�����̶�
		UserDefFloat*	m_data;		//����ֵ������
		unsigned int		m_numperdim[T];		//��ֵ����ά��

		HBXDef::_AEROTABLE* m_table_p;	//��ȡ��_AeroTable��ָ��

		unsigned int		m_lag_cordinate[T];	//��ֵ�����λ��
		_USERDEFT*	p_UserDefT;	

	public:
		//�ù��캯�������ݲ�ֵ���ָ�뼰�ñ��ڵı��
		__host__ CuInterpolation( HBXDef::_AEROTABLE* _IptTable, size_t _blkId  );

		//�ù��캯�������ݲ�ֵ���������ݣ������������ֵ����������Ϣ��
		__host__ CuInterpolation( const HBXDef::CBaseLag<T>& _rhs );

		__host__ ~CuInterpolation(void);

		//@__blkId:��������block��ID,��0Ϊ��ʼ����
		//@UserStatusError_t:����״̬������CheckUserDefErrors���ж�
		__host__ UserStatusError_t SetBlkId( size_t _blkId );

		//��Ҫ��ȷ������ֵ�����Ŀ����������ɴ���ֵ�������������Դ��ϵķ���
		//@_ArrayLgth������ֵ�����Ŀ
		__host__ bool Initial( unsigned int _ArrayLgth );

		//��BaseLag�л�ȡ���ݲ���ʼ��
		__host__ bool Initial( const CBaseLag<T>& _rhs );

		//��ͬһ��ֵ�������ɸ�����в�ֵ����
		//__global__ void  ReadTableArray( _USERDEFT* _ArrayIn );

		//�ú�������CUDA�����ڵ�global�����ڵ���
		//@_ArrayIn:����Ĵ���ֵ����������
		__device__ HBXDef::UserDefFloat	ReadTableData( HBXDef::UserDefFloat* _ArrayIn );

		//Ѱ�Ҵ���ֵ����Χ�ɿռ������ɵĵ�
		__host__ __device__ HBXDef::UserDefFloat	findcordi( unsigned int _tmp_position[TPOW2(T)], HBXDef::UserDefFloat _tmp_proper[T], HBXDef::UserDefFloat _vecIn[T] );
	};

	template<unsigned int T, HBXDef::CudaMalloc_t M >
	__global__	void cutestlag1( HBXDef::CuInterpolation<T, M>* _ArrayIn, HBXDef::UserDefFloat* _outData ); 

}





#include "CuInterpolation.inl"