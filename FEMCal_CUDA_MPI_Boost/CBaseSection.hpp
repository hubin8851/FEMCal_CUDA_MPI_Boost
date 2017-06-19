#ifndef _CHBXFEMBASESEC_HPP_
#define _CHBXFEMBASESEC_HPP_
#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//全局变量，设置数据读取和计算的精度(float or double)

//2017-04-28, HBX, 定义单元类型宏//

namespace HBXFEMDef
{

	//截面定义结构体
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
		//此处以后可以改成联合体以节约内存，2017-4-26，hbx，备注
		T sThick;//截面厚度
		T dA;	//横截面面积
		//T dIz;	//横截面惯性矩
		T dWide;	//横截面宽度
		T dH;		//横截面高
		T InterPoint;	//积分点，我也不知道干啥用的...
		char chSectionName[80];
		char elsetName[32];
		char chMaterialName[80];//截面名称和截面指定的材料名称
	};
}

#endif	_CHBXFEMBASESEC_HPP_