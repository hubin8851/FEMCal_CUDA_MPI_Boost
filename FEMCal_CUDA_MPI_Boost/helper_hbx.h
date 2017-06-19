#ifndef _HELPER_HBX_
#define _HELPER_HBX_

#pragma once 
#include "HbxDefMacro.h"
#include "MatlabPreDef.h"


//���ܻ���ض��壬�ڴ���С��...
static const char *_cudaGetErrorEnum(ImpMatError error)
{
	switch (error)
	{
	case IMPSUCCESS:
		return "matlab mat�ļ��򿪳ɹ�";
	case MATOPENFAILD:
		return "matlab mat�ļ���ʧ��";
	case MATFILENULL:
		return "matlab mat�ļ�Ϊ��";
	case MATGETDIRERROR:
		return "��ȡmatlab mat�ļ�ʱ����";
	case MATGETVARERROR:
		return "��ȡArray����ʱ����";
	case ARRAYNULL:
		return "Array�ǽṹ��";

	default:
		return "unkown";
	}
}

//���ܻ���ض��壬�ڴ���С��...
static const char *_cudaGetErrorEnum( UserStatusError_t _err )
{
	switch (_err)
	{
	case USER_STATUS_NOT_INITIAL:
		return "δ��ȷ��ʼ��";
	case USER_STATUS_ALLOC_FAILED:
		return "�ڴ����ʧ��";
	case USER_STATUS_INVALID_VALUE:
		return "����Чֵ";
	case USER_EXCUTION_FAILED:
		return "��������ʧ��";
	case USER_EMPTY_POINT:
		return "�������п�ָ��";
	default:
		return "<unknown>";
	}
}

//���ܻ���ض��壬�ڴ���С��...
static const char *_cudaGetErrorEnum( HBXDef::FactorError_t _err )
{
	using namespace HBXDef;
	switch (_err)
	{
	case RANKERROR:
		return "��ʽ�ֽ�rank����";
	case NOTSYMETRY:
		return "����ǶԳ�";
	case NOTREGULAR:
		return "���������";
	case SYMETRYbutNOREGULAR:
		return "�������ǶԳ�";
	default:
		break;
	}
}

//���ܻ���ض��壬�ڴ���С��...
static const char *_cudaGetErrorEnum( HBXDef::CSRFormatError_t _err )
{
	using namespace HBXDef;
	switch (_err)
	{
	case NNZCHECKFAILED:
		return "CSR��ʽ��NNZ������";
	case BASEINDEXCHECKFAILED:
		return "��������Ԫ�ؼ�����";
	case CORRUPTEDROW:
		return "��������ĳǰһԪ�ش��ں�һԪ��";
	case COLVSBASEINDEXCHECK:
		return "������Խ�磬С����Ԫ��";
	case SORTINGCOLINDEXCHECK:
		return "����������ĳԪ�ش��ڵ��ں���";
	default:
		break;
	}
}

namespace HBXDef
{
	//���������ռ�����Ϊ�˱�����������Ⱦ�����ڵ�ǰ�����ռ���ʹ��CheckUserDefErrors�����������չ�ĺ�������
	template< typename T >
	void check_userdef( T result, char const *const func, const char *const file, int const line )
	{
		if (result)
		{
			fprintf( stderr, "��������� %s:%d code=%d(%s) \"%s\" \n",
				file, line, static_cast<unsigned int>(result), _cudaGetErrorEnum(result), func );
			DEVICE_RESET	// ȷ���˳�ǰ�����豸������
				getchar();
				exit(EXIT_FAILURE);
		}
	}
#ifndef CheckUserDefErrors(val)
	//ͨ��boost��log��������ض�����δ�޸ģ�2017-03-30��hbx
#define	CheckUserDefErrors(val)	check_userdef( (val), #val, __FILE__, __LINE__ )
#endif



	template< typename T>
	void check_alloc( T* result, const char* _tip, char const *const file, int const line )
	{
		using namespace boost;
		if (!result)	
		{	
 			fprintf( stderr, "�ڴ��������� %s:%d code=\"%s\" \n",
 				file, line, _tip );	
			exit(EXIT_FAILURE);
		}
	}


//���������ռ�����Ϊ�˱�����������Ⱦ�����ڵ�ǰ�����ռ���ʹ��IFNULL�����������չ�ĺ�������
#ifndef IFNULL(rslt�� _tip)
#define IFNULL(rslt, _tip) check_alloc( rslt, _tip, __FILE__, __LINE__ );
#endif


}




#endif