#ifndef _CHBXFEMBASESEC_HPP_
#define _CHBXFEMBASESEC_HPP_
#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//ȫ�ֱ������������ݶ�ȡ�ͼ���ľ���(float or double)

//2017-04-28, HBX, ���嵥Ԫ���ͺ�//

namespace HBXFEMDef
{

	//���涨��ṹ��
	template <class T = double>
	struct _BaseSection 
	{
	public:
		_BaseSection( T _thick=0, const char *_SecName = nullptr, const char *_elset = nullptr, const char *_MatName = nullptr ):sThick(_thick)
		{
			strcpy(this->chSectionName, _SecName);
			strcpy(this->elsetName, _elset);
			strcpy(this->chMaterialName, _MatName);
		}
		//�˴��Ժ���Ըĳ��������Խ�Լ�ڴ棬2017-4-26��hbx����ע
		T sThick;//������
		T dA;	//��������
		//T dIz;	//�������Ծ�
		T dWide;	//�������
		T dH;		//������
		T InterPoint;	//���ֵ㣬��Ҳ��֪����ɶ�õ�...
		char chSectionName[80];
		char elsetName[32];
		char chMaterialName[80];//�������ƺͽ���ָ���Ĳ�������
	};
}

#endif	_CHBXFEMBASESEC_HPP_