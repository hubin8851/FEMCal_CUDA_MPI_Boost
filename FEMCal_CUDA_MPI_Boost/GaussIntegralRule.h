#pragma once
#include "HBXFEMDefStruct.h"
#include "IntegralRule.h"


namespace HBXFEMDef
{

	class CGaussIntegralRule :
		public CIntegralRule
	{
	private:
	protected:
		//@iPoint£º»ý·Ö½×Êý
		int SetUpPointOnLine( int iPoint, MaterialType_t _material_t );
	public:
		CGaussIntegralRule(void);
		~CGaussIntegralRule(void);

		int SetUpIntegrationPoint( Integral_t _Int_type  ){};
	};


}



