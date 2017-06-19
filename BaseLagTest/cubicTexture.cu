#ifndef _CUBICTEXTURE_CU_
#define _CUBICTEXTURE_CU_

#include <stdio.h>
#include "HbxDefMacro.h"

#include "device_launch_parameters.h"
#include <usercuda/helper_gl.h>
//#include <usercuda/GL/freeglut.h>

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
// Helper functions
#include <helper_functions.h> // CUDA SDK Helper functions
#include <helper_cuda.h>// CUDA device initialization helper functions
//#include <helper_cuda_gl.h>   // CUDA device + OpenGL initialization functions
#include "cubicTexture_kernel.cuh"



cudaArray *d_Array = 0;




//以下为一维数组的相关初始化操作
extern "C"
void initTexture1D( size_t _wide, size_t _height, float1* h_data )
{
	cudaChannelFormatDesc _channelDesc = cudaCreateChannelDesc<float>();
	checkCudaErrors( cudaMallocArray(&d_Array, &_channelDesc, _wide, _height) );

	//设置纹理属性
	tex1d.addressMode[0] = cudaAddressModeClamp;//寻址模式，cudaAddressModeClamp表示超出区域钳位到寻址的极值，cudaAddressModeWrap则表示将超出范围插值，且仅支持归一化坐标
	tex1d.addressMode[1] = cudaAddressModeClamp;//对非归一化的坐标，如果寻址的坐标超过了范围[0，N]，大于N的坐标将被钳位，设为N-1。
	tex1d.filterMode = cudaFilterModeLinear;	//滤波模式，cudaFilterModePoint返回最接近的点，cudaFilterModeLinear则返回线性插值的数值
	tex1d.normalized = false;
	checkCudaErrors( cudaBindTextureToArray(tex1d, d_Array, _channelDesc) );
}

//@_wide:二维插值表的宽度
//@_height：二维插值表的高度
//@_blckSize：块内线程数
//@_gridSize：grid内线程数
//@filter_mode：插值方法
//@_output：输出数组
extern "C"
void testTexture1D(float1* d_input, size_t _width, size_t _height, dim3 _blckSize, dim3 _gridSize, int filter_mode, HBXDef::UserDefFloat* _output)
{
	initTexture1D(_width, _height, d_input);
	switch (filter_mode)
	{
	case interpolate_t::NEAREST:
		texref.filterMode = cudaFilterModePoint;//改用就附近的点的插值方法
		d_render1D<<<_gridSize, _blckSize>>>( d_input, _width, _height, _output );
		break;
	default:
		break;
	}

}


//按二维数组处理数据
//@_wide:数据宽度
//@_height：数据高度
//@h_data：数据首地址，要求地址连续
extern "C"
void initTexture2D( size_t _wide, size_t _height, float* h_data )
{
	//将内存内数据拷贝至显存
	//	cudaChannelFormatDesc _channelDesc = cudaCreateChannelDesc(sizeof(HBXDef::UserDefFloat), 0, 0, 0, cudaChannelFormatKindFloat);
	cudaChannelFormatDesc _channelDesc = cudaCreateChannelDesc<float>();


	checkCudaErrors( cudaMallocArray(&d_Array, &_channelDesc, _wide, _height) );
	unsigned int tmp_size = _wide *_height *sizeof(HBXDef::UserDefFloat);


	checkCudaErrors( cudaMemcpyToArray(d_Array, 0, 0, h_data, tmp_size, cudaMemcpyHostToDevice) );
	//free(h_data);

	//设置纹理属性
	texref.addressMode[0] = cudaTextureAddressMode::cudaAddressModeBorder;//寻址模式，cudaAddressModeClamp表示超出区域钳位到寻址的极值，cudaAddressModeWrap则表示将超出范围插值，且仅支持归一化坐标
	texref.addressMode[1] = cudaTextureAddressMode::cudaAddressModeBorder;//对非归一化的坐标，如果寻址的坐标超过了范围[0，N]，大于N的坐标将被钳位，设为N-1。
	texref.filterMode = cudaFilterModeLinear;	//滤波模式，cudaFilterModePoint返回最接近的点，cudaFilterModeLinear则返回线性插值的数值
	texref.normalized = false;

	checkCudaErrors( cudaBindTextureToArray(texref, d_Array, _channelDesc) );

}


//@_wide:二维插值表的宽度
//@_height：二维插值表的高度
//@_blckSize：块内线程数
//@_gridSize：grid内线程数
//@filter_mode：插值方法
//@_output：输出数组
extern "C"
void testTexture2D(float2* d_input, size_t _width, size_t _height, dim3 _blckSize, dim3 _gridSize, int filter_mode, HBXDef::UserDefFloat* _output)
{
	float* interpotdata = new float[_width*_height];
	for (int i = 0; i < _height; i++)
	{
		for (int j = 0; j < _width; j++)
		{
			interpotdata[i*_width+j] = i*_width+j;
		}
	}

	initTexture2D(_width, _height, interpotdata);
	switch (filter_mode)
	{
	case interpolate_t::NEAREST:
//		texref.filterMode = cudaFilterModePoint;//改用就附近的点的插值方法
		d_render2D<<<_gridSize, _blckSize>>>( d_input, _width, _height, _output );
		break;
	default:
		break;
	}
	cudaThreadSynchronize();
	std::cerr<<cudaGetErrorName(cudaGetLastError())<<std::endl;
}


#endif