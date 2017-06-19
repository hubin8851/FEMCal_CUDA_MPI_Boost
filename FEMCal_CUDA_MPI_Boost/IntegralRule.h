#pragma once

#include "HBXFEMDefStruct.h"

namespace HBXFEMDef
{
	/* 针对不同单元，有多种积分形式（完全/缩减、高斯/牛顿-柯兹/复合辛普森 等）。本类隐藏于
	算法内，仅设置正确的积分点和权值。鉴于动力学内材料的属性会发生变化，需内置材料的shared_ptr，且与
	每个单元相对应
	*/
	class CIntegralRule
	{
	private:
	protected:
		boost::shared_ptr<HBXFEMDef::_BaseElem<double>> _elem;
		bool	isDynamic;//动力学标志位，如果非动力学不记录之前的积分点；动力学中需记录高斯积分点在运算过程中的变化情况

		virtual int SetUpPointOnLine( int iPoint, MaterialType_t _material_t )=0;
	public:
		CIntegralRule(void);
		virtual ~CIntegralRule(void);

		//
		int SetUpIntegrationPoint( Integral_t _Int_type  );
	};
}


