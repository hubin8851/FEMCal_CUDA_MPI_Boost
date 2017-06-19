#pragma once
#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//全局变量，设置数据读取和计算的精度(float or double)

namespace HBXFEMDef
{
	//框架单元结构体
	template <class T>
	class _FrameElement:
		public _BaseElem<T>
	{
	protected:
		short iSection;	//单元截面索引号
		T dLength;	//单元长度
		T dSin, dCos;	//单元局部坐标x轴与整体坐标x轴的夹角的正余弦
		T daEndInterForce[6];	//单元杆端力向量

	public:
		void StiffMatCal(){};
		T VolumeCalc(){return 0;};
		void MassMatCal( MassType_t _MassType = MassType_t::LUMPED ){};
		void SetMaterial( std::pair< boost::shared_ptr<_Material<T>>, boost::shared_ptr<_BaseSection<T>> > _MatAndSec ){ };

	};

	typedef boost::shared_ptr<_FrameElement<float>> _FrameElement_f_ptr;
	typedef boost::shared_ptr<_FrameElement<double>> _FrameElement_d_ptr;
}