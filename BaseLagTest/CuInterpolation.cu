#include <stdio.h>

#include "XMLManage.h"
#include "HbxDefMacro.h"
#include "CuComponent.cuh"

#include "CuInterpolation.cuh"


namespace HBXDef
{

// 	template<unsigned int T, HBXDef::CudaMalloc_t M >
// 	__global__	void cutestlag1( HBXDef::CuInterpolation<T, M>* _ArrayIn, HBXDef::UserDefFloat* _outData )
// 	{
// 		using namespace HBXDef;
// 		unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
// 
// 		while ( idx < MAX_GPUITER )
// 		{
// #ifdef USERDEFSHARED//为了控制共享内存大小
// 
// #endif // USERDEFSHARED//为了控制共享内存大小
// 
// #ifdef WITHTEST	//测试用例，若传入指针为空则生成
// //		std::cerr<<"传入指针为空..."<<std::endl;
// 
// 			UserDefFloat intercord[4];
// 			for (size_t i = 0; i < RUNCLC; i++)
// 			{
// 				intercord[0] = -3.0 + i * 20.0 / RUNCLC + 20.0 * i / MAX_GPUITER;
// 				intercord[1] = 0.2 + i * 2.0 / RUNCLC + 2.0 * i / MAX_GPUITER;
// 				intercord[2] = -30 + i * 30.0 / RUNCLC + 30 * i / MAX_GPUITER;
// 				intercord[3] = -5.0 + 5*i / RUNCLC + i*5 / MAX_GPUITER;
// 
// 				_outData[i + RUNCLC*idx] = _ArrayIn[idx].ReadTableData(intercord);
// 				
// 			}
// 			idx += gridDim.x*blockDim.x;
// #endif  //WITHTEST
// 		}
// 		__syncthreads();
// 	}



}



extern "C"
void testlag()
{
	using namespace HBXDef;
	CXMLManage	g_Xml;		//用于读取单元属性的xml
	g_Xml.ReadAeroCoef("F:\\data from HBX_phd\\database\\x33气动数据\\CA_d0.xml");

	const unsigned int g_T = 2;
	typedef	HBXDef::CBaseLag<g_T, HBXDef::CudaMalloc_t::UNIFIEDMEM>	_samelag;
	typedef	HBXDef::CuInterpolation<g_T, HBXDef::CudaMalloc_t::UNIFIEDMEM>	_sameCulag;
	float elapsedTime;

	HBXDef::UserDefFloat* g_opData = nullptr;
#ifdef WITHTEST
	std::cout<<cudaGetErrorName( cudaMallocManaged(&g_opData, RUNCLC * MAX_GPUITER* sizeof(HBXDef::UserDefFloat)) )<<std::endl;
#endif // WITHTEST


	_samelag*	g_BaseLag = new _samelag( &g_Xml.GetAeroTable("CA_d0"), 0 );	//用以插值
	g_BaseLag->SetBlkId(0);
	g_BaseLag->Initial();

	_sameCulag*	g_culag = nullptr;
	std::cout<<sizeof(_sameCulag)<<std::endl;
	std::cout<<cudaGetErrorName( cudaMalloc(&g_culag, MAX_GPUITER* sizeof(_sameCulag)) )<<std::endl;
	//	g_culag = (_sameCulag*)malloc(MAX_GPUITER* sizeof(_sameCulag));
	if (!g_culag)
	{
		cout << "内存分配错误" << endl;
		return;
	}

	for (int i = 0; i < MAX_GPUITER; i++)
	{
//		g_culag[i].Initial(*g_BaseLag);
	}

	//串行算法
	clock_t _beg, _end;
	HBXDef::UserDefFloat _rlt;
	_beg = clock();
	HBXDef::UserDefFloat  g_dataIn[4];
	for (size_t i = 0; i < MAX_GPUITER; i++) {
		for (size_t j = 0; j < RUNCLC; j++) {
			g_dataIn[0] = -3.0 + j * 20.0 / RUNCLC + 20.0 * i / MAX_GPUITER;
			g_dataIn[1] = 0.2 + j * 2.0 / RUNCLC + 2.0 * i / MAX_GPUITER;
			g_dataIn[2] = -30 + j * 30.0 / RUNCLC + 30 * i / MAX_GPUITER;
			g_dataIn[3] = -5.0 + 5*j / RUNCLC + i*5 / MAX_GPUITER;
			_rlt = g_BaseLag->ReadTableData(g_dataIn);
		}
	}
	_end = clock();
	cout << "CPU time:" << _end - _beg << endl;

#ifdef WITHTEST
	std::cout<<_rlt<<std::endl;
#else
	std::cout<<_rlt<<std::endl;
#endif // WITHTEST
	//并行算法

	cuerror(cudaSetDevice(0));

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);//开始计时

#ifdef WITHTEST //测试用例
	cutestlag1<g_T, CudaMalloc_t::UNIFIEDMEM><<<GRIDSIZE, BLOCKSIZE>>>(g_culag, g_opData);
	cudaEventRecord(stop, 0);//计时结束
	std::cerr<<cudaGetErrorName(cudaGetLastError())<<endl;
	cudaThreadSynchronize();
	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);//获取用时
#else
	cutestlag1<g_T, CudaMalloc_t::UNIFIEDMEM><<<GRIDSIZE, BLOCKSIZE, MAX_SHARED, 0 >>>(g_culag, nullptr);
	cudaEventRecord(stop, 0);//计时结束
	std::cerr<<cudaGetErrorString(cudaGetLastError())<<endl;
	cudaThreadSynchronize();
	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);//获取用时
#endif // WITHTEST


	std::cout << "time eclaps:" << elapsedTime << std::endl;
	system("pause");
}