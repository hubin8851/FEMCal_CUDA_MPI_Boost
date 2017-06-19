#ifndef _HELPER_HBX_
#define _HELPER_HBX_

#pragma once 
#include "HbxDefMacro.h"
#include "MatlabPreDef.h"


//可能会多重定义，在此需小心...
static const char *_cudaGetErrorEnum(ImpMatError error)
{
	switch (error)
	{
	case IMPSUCCESS:
		return "matlab mat文件打开成功";
	case MATOPENFAILD:
		return "matlab mat文件打开失败";
	case MATFILENULL:
		return "matlab mat文件为空";
	case MATGETDIRERROR:
		return "读取matlab mat文件时出错";
	case MATGETVARERROR:
		return "获取Array变量时出错";
	case ARRAYNULL:
		return "Array非结构体";

	default:
		return "unkown";
	}
}

//可能会多重定义，在此需小心...
static const char *_cudaGetErrorEnum( UserStatusError_t _err )
{
	switch (_err)
	{
	case USER_STATUS_NOT_INITIAL:
		return "未正确初始化";
	case USER_STATUS_ALLOC_FAILED:
		return "内存分配失败";
	case USER_STATUS_INVALID_VALUE:
		return "非有效值";
	case USER_EXCUTION_FAILED:
		return "函数调用失败";
	case USER_EMPTY_POINT:
		return "函数内有空指针";
	default:
		return "<unknown>";
	}
}

//可能会多重定义，在此需小心...
static const char *_cudaGetErrorEnum( HBXDef::FactorError_t _err )
{
	using namespace HBXDef;
	switch (_err)
	{
	case RANKERROR:
		return "因式分解rank错误";
	case NOTSYMETRY:
		return "矩阵非对称";
	case NOTREGULAR:
		return "矩阵非正定";
	case SYMETRYbutNOREGULAR:
		return "正定但非对称";
	default:
		break;
	}
}

//可能会多重定义，在此需小心...
static const char *_cudaGetErrorEnum( HBXDef::CSRFormatError_t _err )
{
	using namespace HBXDef;
	switch (_err)
	{
	case NNZCHECKFAILED:
		return "CSR格式的NNZ检测出错";
	case BASEINDEXCHECKFAILED:
		return "行索引首元素检测出错";
	case CORRUPTEDROW:
		return "行索引中某前一元素大于后一元素";
	case COLVSBASEINDEXCHECK:
		return "列索引越界，小于首元素";
	case SORTINGCOLINDEXCHECK:
		return "行内列索引某元素大于等于后者";
	default:
		break;
	}
}

namespace HBXDef
{
	//放在命名空间内是为了避免命名的污染，仅在当前命名空间内使用CheckUserDefErrors宏命令调用扩展的函数对象
	template< typename T >
	void check_userdef( T result, char const *const func, const char *const file, int const line )
	{
		if (result)
		{
			fprintf( stderr, "程序错误，在 %s:%d code=%d(%s) \"%s\" \n",
				file, line, static_cast<unsigned int>(result), _cudaGetErrorEnum(result), func );
			DEVICE_RESET	// 确认退出前所有设备已重置
				getchar();
				exit(EXIT_FAILURE);
		}
	}
#ifndef CheckUserDefErrors(val)
	//通过boost的log完成流的重定向，尚未修改，2017-03-30，hbx
#define	CheckUserDefErrors(val)	check_userdef( (val), #val, __FILE__, __LINE__ )
#endif



	template< typename T>
	void check_alloc( T* result, const char* _tip, char const *const file, int const line )
	{
		using namespace boost;
		if (!result)	
		{	
 			fprintf( stderr, "内存分配错误，在 %s:%d code=\"%s\" \n",
 				file, line, _tip );	
			exit(EXIT_FAILURE);
		}
	}


//放在命名空间内是为了避免命名的污染，仅在当前命名空间内使用IFNULL宏命令调用扩展的函数对象
#ifndef IFNULL(rslt， _tip)
#define IFNULL(rslt, _tip) check_alloc( rslt, _tip, __FILE__, __LINE__ );
#endif


}




#endif