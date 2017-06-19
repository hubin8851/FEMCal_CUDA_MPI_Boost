#pragma once

#include "StdAfx.h"

template< HBXDef::CudaMalloc_t M >
class CuComponent
{
protected:
	size_t m_DevID;
public:
	CuComponent(void);
	~CuComponent(void);

	void *operator new( size_t _lgt );
	void operator delete( void* );

	size_t	GetGPUDeviceQuery();
};


#include "CuComponent.inl"





