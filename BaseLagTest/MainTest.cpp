#include "stdafx.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <helper_functions.h> // CUDA SDK Helper functions
#include <helper_cuda.h>// CUDA device initialization helper functions

#include "XMLManage.h"
#include "BaseLag.hpp"
//#include "CuComponent.cuh"
#include "CuInterpolation.cuh"
//#include "cubicTexture_kernel.cuh"

extern "C" void testTexture1D( float2* d_input, size_t _wide, size_t _height, dim3 _blckSize, dim3 _gridSize, int filter_mode, HBXDef::UserDefFloat* _output);
extern "C" void testTexture2D( float2* d_input, size_t _wide, size_t _height, dim3 _blckSize, dim3 _gridSize, int filter_mode, HBXDef::UserDefFloat* _output);
extern "C" void testlag();


int *pArgc = NULL;
char **pArgv = NULL;
HBXDef::FilterMode_t g_FilterMode = HBXDef::FilterMode_t::MODE_NEAREST;


int main(int argc, char** argv) {
	pArgc = &argc;
	pArgv = argv;

	unsigned int g_width = 512, g_height = 512;
	dim3 blockSize(16, 16);
	dim3 gridSize(g_width / blockSize.x, g_height / blockSize.y);

	if (checkCmdLineFlag(argc, (const char**)argv, "mode"))
	{
		g_FilterMode = (HBXDef::FilterMode_t)getCmdLineArgumentInt(argc, (const char**)argv, "mode");
		if (g_FilterMode < 0 || g_FilterMode >= HBXDef::FilterMode_t::NUM_MODES)
		{
			printf("Invalid Mode setting %d\n", g_FilterMode);
			exit(EXIT_FAILURE);
		}
	}

	float2*	h_input;
	HBXDef::UserDefFloat*	d_output;	//输出的插值的数值
	HBXDef::UserDefFloat*	h_output;

	checkCudaErrors( cudaMallocManaged((void**)&d_output, sizeof(HBXDef::UserDefFloat)*g_width*g_height) );
#ifdef WITHTEST
	//对输入的数组赋值
	checkCudaErrors( cudaMallocManaged((void**)&h_input, sizeof(float2)*g_width*g_height) );
	for (int i = 0; i < g_width; i++)
	{
		for (int j = 0; j < g_height; j++)
		{
//			h_input[j*g_width+i].x = -3.0 + j * 20.0 / RUNCLC + 20.0 * i / MAX_GPUITER;
//			h_input[j*g_width+i].y = 0.2 + j * 2.0 / RUNCLC + 2.0 * i / MAX_GPUITER;
			h_input[j*g_width+i].x = 0;
			h_input[j*g_width+i].y = 0.6;
		}
	}
#endif // WITHTEST

	testTexture2D( h_input, g_width, g_height, gridSize, blockSize, g_FilterMode, d_output);
//	checkCudaErrors( cudaMemcpy( h_output, d_output, sizeof(HBXDef::UserDefFloat)*g_width*g_height, cudaMemcpyDeviceToHost ));
	
	//CPU测试
	testlag();
}


