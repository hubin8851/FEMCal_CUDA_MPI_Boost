#ifndef _MATLABPREDEF_H_
#define _MATLABPREDEF_H_


#pragma once

#include "stdafx.h"

#include <matlab/mat.h>

enum ImpMatError
{
	IMPSUCCESS = 0,
	MATOPENFAILD = 1,
	MATGETDIRERROR = 2,
	MATFILENULL = 3,
	MATGETVARERROR = 4,
	ARRAYNULL = 5	//mxarrayÏÂÎª¿Õ
};

typedef ImpMatError ImpMatError_t;

#endif