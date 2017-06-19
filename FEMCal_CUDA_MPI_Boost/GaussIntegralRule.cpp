#include "StdAfx.h"
#include "GaussIntegralRule.h"

namespace HBXFEMDef
{
	CGaussIntegralRule::CGaussIntegralRule(void)
	{
	}


	CGaussIntegralRule::~CGaussIntegralRule(void)
	{
	}


	int CGaussIntegralRule::SetUpPointOnLine( int iPoint, MaterialType_t _material_t )
	{
		double _weight;
		switch (iPoint)
		{
		case 1:
			break;
		case 2:
			break;
		default:
			break;
		}
		return 0;
	}

}

