#pragma once

#include "stdafx.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace HBXDef
{
	template< HBXDef::CudaMalloc_t M >
	class CuComponent
	{
	private:
		enum 
		{
			_TYPE = M
		};
	protected:
		size_t m_DevID;

	public:
		CuComponent(void){};
		~CuComponent(void){};

		void *operator new( size_t _lgt )
		{
			using namespace HBXDef;
			void* ptr;
			switch (M)
			{
			case NORMAL:
				cudaMalloc( &ptr, _lgt );
				break;
			case PAGELOCK:
				cudaHostAlloc( &ptr, _lgt, cudaHostAllocMapped );
				break;
			case UNIFIEDMEM:
				cudaMallocManaged( &ptr, _lgt );
				break;
			}
			return ptr;
		}

		void operator delete( void* _data)
		{
			switch (M)
			{
			case HBXDef::CudaMalloc_t::NORMAL: 
			case HBXDef::CudaMalloc_t::UNIFIEDMEM: 
				cudaFree( _data );
				break;
			case HBXDef::CudaMalloc_t::PAGELOCK:
				cudaFreeHost( _data );
				break;
			default:
				break;
			}
			return;
		}

		size_t	GetGPUDeviceQuery()
		{
			int deviceCount = 0;
			cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

			if (error_id != cudaSuccess)
			{
				printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
				printf("Result = FAIL\n");
				return 0;
			}
			if (deviceCount == 0)
			{
				printf("\nThere are no available device(s) that support CUDA\n");
				return 0;
			}
			else
			{
				printf("Detected %d CUDA Capable device(s)\n", deviceCount);
			}
			int dev, driverVersion = 0, runtimeVersion = 0;
			cudaDeviceProp	_BestDevice;//最好的GPU的属性结构体
			checkCudaErrors(cudaGetDeviceProperties(&_BestDevice, 0));

			for (dev = 0; dev < deviceCount; ++dev)
			{
				cudaSetDevice(dev);
				cudaDeviceProp deviceProp;
				cudaGetDeviceProperties(&deviceProp, dev);

				printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

				// Console log
				cudaDriverGetVersion(&driverVersion);
				cudaRuntimeGetVersion(&runtimeVersion);
				printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n", driverVersion/1000, (driverVersion%100)/10, runtimeVersion/1000, (runtimeVersion%100)/10);
				printf("  CUDA Capability Major/Minor version number:    %d.%d\n", deviceProp.major, deviceProp.minor);

				char msg[256];
				sprintf_s(msg, "  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n",
					(float)deviceProp.totalGlobalMem/1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);
				printf("%s", msg);

				printf("  (%2d) Multiprocessors, (%3d) CUDA Cores/MP:     %d CUDA Cores\n",
					deviceProp.multiProcessorCount,
					_ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
					_ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount);
				printf("  GPU Max Clock rate:                            %.0f MHz (%0.2f GHz)\n", deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);
				if ((_BestDevice.major<=deviceProp.major)&&(_BestDevice.minor<=deviceProp.minor))
				{
					m_DevID = dev;
					// 因SM内线程限制的增加，线程块内的线程数量可相应增加
					int m_BlockSize = (deviceProp.major < 2) ? 16 : 32;
				}
				printf("\n共享内存:%u Byte", deviceProp.sharedMemPerBlock);
			}
			return deviceCount;
		}
	};

	typedef CuComponent<HBXDef::CudaMalloc_t::UNIFIEDMEM>	_CuBase_uni_;
}









