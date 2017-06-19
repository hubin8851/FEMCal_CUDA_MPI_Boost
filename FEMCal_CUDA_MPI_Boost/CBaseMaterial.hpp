#ifndef _BASEMATERIAL_HPP_
#define _BASEMATERIAL_HPP_

#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//ȫ�ֱ������������ݶ�ȡ�ͼ���ľ���(float or double)

namespace HBXFEMDef
{
	//���϶���ṹ��
	template <class T>
	class _Material
	{
	public:
		_Material(T _E=0, T _mu=0, T _Gama=0, const char *ch = nullptr ):dE(_E), dMu(_mu), dMaGama(_Gama)
		{
			strcpy(this->chname, ch);
			dG = dE / (2 * (1+dMu));
		}

		T dE;	//����ģ��
		T dMu;	//���ɱ�
		T dG;	//����ģ��
		T dAlpha;	//������ϵ��
		T dMaGama; //�ض�
		char   chname[40];//���������ַ�����
	};
}


#endif // !_BASEMATERIAL_HPP_
