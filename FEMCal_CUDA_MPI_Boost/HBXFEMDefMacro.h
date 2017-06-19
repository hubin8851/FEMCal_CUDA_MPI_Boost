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
	} Integral_t;//��������ʹ�ã��Ա����������

	typedef enum 
	{
		CONSISTENT = 0,	//һ��������
		LUMPED = 1		//����������	
	} MassType_t;

	typedef enum
	{
		D2BEAM,
		D3BEAM,
		D3SHELL
	} MaterialType_t;


}

#endif	_HBXFEMDEFMACRO_H