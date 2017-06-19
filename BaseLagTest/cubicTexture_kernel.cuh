#ifndef _CUBICTEXTURE_CUH_
#define _CUBICTEXTURE_CUH_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "device_launch_parameters.h"
#include <usercuda/helper_gl.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <helper_math.h>
#include <helper_cuda.h>

#include "texture_fetch_functions.h"

typedef enum interpolate
{
	NEAREST = 0,
	BILINEAR,
	BICUBIC,
	FAST_BICUBIC,
	CATROM
} interpolate_t;

//纹理内存的相关应用参见自己翻译的using texture memory of cuda gpu.doxc


texture<float, 1> tex1d;//一维纹理参照系

__global__ void d_render1D(float1* d_input, unsigned int _wide, unsigned int _height, HBXDef::UserDefFloat* d_output)
{
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int i = y * _wide + x;
	float c = tex1D(tex1d, d_input[i].x);
	d_output[i] = c;
}



//@NO.1:数据类型
//@NO.2：cudaTextureType*D，*表示纹理内存的维度。cudaTextureType*DLayered相应对应层级纹理内存的维度
//@NO.3：定义读取模式，默认值为cudaReadModeElementType
texture<float, 2, cudaReadModeElementType> texref;	//纹理参照系
texture<float, 2, cudaReadModeElementType> tex2;    // need to use cudaReadModeElementType for tex2Dgather

__global__ void d_render2D(float2* d_input, unsigned int _wide, unsigned int _height, HBXDef::UserDefFloat* d_output)
{
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int i = y * _wide + x;

//	d_output[i] = tex2D(texref, d_input[i].x+0.02*i, d_input[i].y+0.02*i);
	d_output[i] = tex2D(texref, d_input[i].x, d_input[i].y);
//	d_output[i] = i;
//	tex2D()
//	d_output[i] = 1;
}


__global__ void d_render3D(float3* d_input, unsigned int _wide, unsigned int _height, HBXDef::UserDefFloat* d_output)
{

}

#endif _CUBICTEXTURE_CUH_