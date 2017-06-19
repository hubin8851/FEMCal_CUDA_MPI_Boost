#include "CuComponent.cuh"
#include "CuInterpolation.h"

#include <stdio.h>

namespace HBXDef
{

	template<unsigned int T, HBXDef::CudaMalloc_t M >
	__global__	void cutestlag1( HBXDef::CuInterpolation<T, M>* _ArrayIn, HBXDef::UserDefFloat* _outData )
	{
		unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

		while ( idx < MAX_GPUITER )
		{
#ifdef USERDEFSHARED//Ϊ�˿��ƹ����ڴ��С

#endif // USERDEFSHARED//Ϊ�˿��ƹ����ڴ��С

#ifdef WITHTEST	//����������������ָ��Ϊ��������
//		std::cerr<<"����ָ��Ϊ��..."<<std::endl;
#endif  
			HBXDef::UserDefFloat intercord[4];
			for (size_t i = 0; i < RUNCLC; i++)
			{
				intercord[0] = -3.0 + i * 20.0 / RUNCLC + 20.0 * i / MAX_GPUITER;
				intercord[1] = 0.2 + i * 2.0 / RUNCLC + 2.0 * i / MAX_GPUITER;
				intercord[2] = -30 + i * 30.0 / RUNCLC + 30 * i / MAX_GPUITER;
				intercord[3] = -5.0 + 5*i / RUNCLC + i*5 / MAX_GPUITER;

				_ArrayIn[idx].ReadTableData(intercord);
				
			}
			idx += gridDim.x*blockDim.x;
		}
		__syncthreads();
	}




}