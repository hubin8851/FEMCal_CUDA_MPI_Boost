#ifndef _BASEMATERIAL_HPP_
#define _BASEMATERIAL_HPP_

#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//全局变量，设置数据读取和计算的精度(float or double)

namespace HBXFEMDef
{
	//材料定义结构体
	template <class T>
	class _Material
	{
	public:
		_Material(T _E=0, T _mu=0, T _Gama=0, const char *ch = nullptr ):dE(_E), dMu(_mu), dMaGama(_Gama)
		{
			strcpy(this->chname, ch);
			dG = dE / (2 * (1+dMu));
		}

		T dE;	//弹性模量
		T dMu;	//泊松比
		T dG;	//剪切模量
		T dAlpha;	//线膨胀系数
		T dMaGama; //重度
		char   chname[40];//材料名称字符数组
	};
}


#endif // !_BASEMATERIAL_HPP_
