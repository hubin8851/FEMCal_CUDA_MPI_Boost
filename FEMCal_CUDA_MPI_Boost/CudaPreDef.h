#pragma once

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cusparse.h>
#include <cublas.h>
#include <cublas_v2.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#ifndef _CUDA_PREDEF_
#define _CUDA_PREDEF_


#include <helper_cusolver.h>
#include <cusolver_common.h>
#include <cusolverSp.h>
#include <device_launch_parameters.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>

namespace HBXDef
{

}

#endif