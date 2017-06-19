// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include "targetver.h"
#include <math.h>
#include <stdio.h>
#include <tchar.h>



// TODO: 在此处引用程序需要的其他头文件
#include <string.h>
#include <tchar.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <valarray>
#include <complex>
#include <math.h>
#include <array>
#include <time.h>	
//mysql头文件
#include<windows.h>
#include <mysql/mysql.h>
//Eigen头文件
#include "EigenPreDef.h"
//matlab头文件
#include "MatlabPreDef.h"
//boost头文件
#include "BoostPreDef.h"
//CUDA头文件
#include "CudaPreDef.h"
//自定义头文件
#include "HbxDefMacro.h"//宏定义在先
#include "HBXDefStruct.h"
#include "CMatrix.hpp"

#include "HBXFEMDefMacro.h"
#include "HBXFEMDefStruct.h"
#include "CTri2D3NElemt.hpp"
#include "C3D8RElemt.hpp"
#include "CFrameElement.hpp"
#include "CPlaneElement.hpp"

#include "helper_hbx.h"