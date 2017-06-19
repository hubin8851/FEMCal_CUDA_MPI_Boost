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




//����Ϊһά�������س�ʼ������
extern "C"
void initTexture1D( size_t _wide, size_t _height, float1* h_data )
{
	cudaChannelFormatDesc _channelDesc = cudaCreateChannelDesc<float>();
	checkCudaErrors( cudaMallocArray(&d_Array, &_channelDesc, _wide, _height) );

	//������������
	tex1d.addressMode[0] = cudaAddressModeClamp;//Ѱַģʽ��cudaAddressModeClamp��ʾ��������ǯλ��Ѱַ�ļ�ֵ��cudaAddressModeWrap���ʾ��������Χ��ֵ���ҽ�֧�ֹ�һ������
	tex1d.addressMode[1] = cudaAddressModeClamp;//�Էǹ�һ�������꣬���Ѱַ�����곬���˷�Χ[0��N]������N�����꽫��ǯλ����ΪN-1��
	tex1d.filterMode = cudaFilterModeLinear;	//�˲�ģʽ��cudaFilterModePoint������ӽ��ĵ㣬cudaFilterModeLinear�򷵻����Բ�ֵ����ֵ
	tex1d.normalized = false;
	checkCudaErrors( cudaBindTextureToArray(tex1d, d_Array, _channelDesc) );
}

//@_wide:��ά��ֵ��Ŀ��
//@_height����ά��ֵ��ĸ߶�
//@_blckSize�������߳���
//@_gridSize��grid���߳���
//@filter_mode����ֵ����
//@_output���������
extern "C"
void testTexture1D(float1* d_input, size_t _width, size_t _height, dim3 _blckSize, dim3 _gridSize, int filter_mode, HBXDef::UserDefFloat* _output)
{
	initTexture1D(_width, _height, d_input);
	switch (filter_mode)
	{
	case interpolate_t::NEAREST:
		texref.filterMode = cudaFilterModePoint;//���þ͸����ĵ�Ĳ�ֵ����
		d_render1D<<<_gridSize, _blckSize>>>( d_input, _width, _height, _output );
		break;
	default:
		break;
	}

}


//����ά���鴦������
//@_wide:���ݿ��
//@_height�����ݸ߶�
//@h_data�������׵�ַ��Ҫ���ַ����
extern "C"
void initTexture2D( size_t _wide, size_t _height, float* h_data )
{
	//���ڴ������ݿ������Դ�
	//	cudaChannelFormatDesc _channelDesc = cudaCreateChannelDesc(sizeof(HBXDef::UserDefFloat), 0, 0, 0, cudaChannelFormatKindFloat);
	cudaChannelFormatDesc _channelDesc = cudaCreateChannelDesc<float>();


	checkCudaErrors( cudaMallocArray(&d_Array, &_channelDesc, _wide, _height) );
	unsigned int tmp_size = _wide *_height *sizeof(HBXDef::UserDefFloat);


	checkCudaErrors( cudaMemcpyToArray(d_Array, 0, 0, h_data, tmp_size, cudaMemcpyHostToDevice) );
	//free(h_data);

	//������������
	texref.addressMode[0] = cudaTextureAddressMode::cudaAddressModeBorder;//Ѱַģʽ��cudaAddressModeClamp��ʾ��������ǯλ��Ѱַ�ļ�ֵ��cudaAddressModeWrap���ʾ��������Χ��ֵ���ҽ�֧�ֹ�һ������
	texref.addressMode[1] = cudaTextureAddressMode::cudaAddressModeBorder;//�Էǹ�һ�������꣬���Ѱַ�����곬���˷�Χ[0��N]������N�����꽫��ǯλ����ΪN-1��
	texref.filterMode = cudaFilterModeLinear;	//�˲�ģʽ��cudaFilterModePoint������ӽ��ĵ㣬cudaFilterModeLinear�򷵻����Բ�ֵ����ֵ
	texref.normalized = false;

	checkCudaErrors( cudaBindTextureToArray(texref, d_Array, _channelDesc) );

}


//@_wide:��ά��ֵ��Ŀ��
//@_height����ά��ֵ��ĸ߶�
//@_blckSize�������߳���
//@_gridSize��grid���߳���
//@filter_mode����ֵ����
//@_output���������
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
//		texref.filterMode = cudaFilterModePoint;//���þ͸����ĵ�Ĳ�ֵ����
		d_render2D<<<_gridSize, _blckSize>>>( d_input, _width, _height, _output );
		break;
	default:
		break;
	}
	cudaThreadSynchronize();
	std::cerr<<cudaGetErrorName(cudaGetLastError())<<std::endl;
}


#endif