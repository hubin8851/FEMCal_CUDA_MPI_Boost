#ifndef _HBXFEMDEFMACRO_H
#define	_HBXFEMDEFMACRO_H

namespace HBXFEMDef
{
	typedef enum
	{
		POINT,
		LINE,
		TRI,
		QUAD,
		TETRA,
		HEX
	} Integral_t;//被积分类使用，以便决定积分域

	typedef enum 
	{
		CONSISTENT = 0,	//一致质量法
		LUMPED = 1		//集中质量法	
	} MassType_t;

	typedef enum
	{
		D2BEAM,
		D3BEAM,
		D3SHELL
	} MaterialType_t;


}

#endif	_HBXFEMDEFMACRO_H