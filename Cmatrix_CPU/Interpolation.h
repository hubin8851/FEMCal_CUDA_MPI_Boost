// Interpolation.h		插值头文件
// Ver 1.0.0.0
// 版权所有(C)  何渝(HE Yu) 2002

// 最后修改: 2002.5.31.

#ifndef _INTERPOLATION_H	//避免多次编译
#define _INTERPOLATION_H

#include <valarray>			//数组模板类标准头文件
#include <Matrix.h>			//矩阵类头文件
#include <comm.h>			//公共头文件
using namespace std;		//名字空间

//一元全区间不等距插值
template <class _Ty>
_Ty Interpolation1VariableNotIsometry(valarray<_Ty>& x, 
											valarray<_Ty>& y, _Ty t);

//一元全区间等距插值
template <class _Ty>
_Ty Interpolation1VariableIsometry(_Ty x0, _Ty h, valarray<_Ty>& y, _Ty t);

//一元三点不等距插值
template <class _Ty>
_Ty Interpolation1Variable3PointsNotIsometry(valarray<_Ty>& x, 
											valarray<_Ty>& y, _Ty t);

//一元三点等距插值
template <class _Ty>
_Ty Interpolation1Variable3PointsIsometry(_Ty x0, _Ty h, 
											valarray<_Ty>& y, _Ty t);

//连分式不等距插值
template <class _Ty>
_Ty InterpolationFractionNotIsometry(valarray<_Ty>& x, 
											valarray<_Ty>& y, _Ty t);

//连分式等距插值
template <class _Ty>
_Ty InterpolationFractionIsometry(_Ty x0, _Ty h, valarray<_Ty>& y, _Ty t);

//埃尔米特不等距插值
template <class _Ty>
_Ty InterpolationHermiteNotIsometry(valarray<_Ty>& x, 
							valarray<_Ty>& y, valarray<_Ty>& dy, _Ty t);

//埃尔米特等距插值
template <class _Ty>
_Ty InterpolationHermiteIsometry(_Ty x0, _Ty h, 
						valarray<_Ty>& y, valarray<_Ty>& dy, _Ty t);

//埃特金不等距逐步插值
template <class _Ty>
_Ty InterpolationAitkenNotIsometry(valarray<_Ty>& x, 
									valarray<_Ty>& y, _Ty t, _Ty eps);

//埃特金等距逐步插值
template <class _Ty>
_Ty InterpolationAitkenIsometry(_Ty x, _Ty h,  
									valarray<_Ty>& y, _Ty t, _Ty eps);

//光滑不等距插值
template <class _Ty>
void InterpolationSmoothNotIsometry(valarray<_Ty>& x, 
					valarray<_Ty>& y, int k, _Ty t, valarray<_Ty>& s);


//光滑等距插值
template <class _Ty>
void InterpolationSmoothIsometry(_Ty x, _Ty h, 
					valarray<_Ty>& y, int k, _Ty t, valarray<_Ty>& s);

//第一种边界条件的三次样条函数插值、微商与积分
template <class _Ty>
_Ty Interpolation3Spooling1stBoundary(valarray<_Ty>& x, 
		valarray<_Ty>& y, valarray<_Ty>& dy, valarray<_Ty>& ddy, 
			valarray<_Ty>& t, valarray<_Ty>& z, valarray<_Ty>& dz, 
													valarray<_Ty>& ddz);

//第二种边界条件的三次样条函数插值、微商与积分
template <class _Ty>
_Ty Interpolation3Spooling2ndBoundary(valarray<_Ty>& x, 
		valarray<_Ty>& y, valarray<_Ty>& dy, valarray<_Ty>& ddy, 
			valarray<_Ty>& t, valarray<_Ty>& z, valarray<_Ty>& dz, 
													valarray<_Ty>& ddz);

//第三种边界条件的三次样条函数插值、微商与积分
template <class _Ty>
_Ty Interpolation3Spooling3thBoundary(valarray<_Ty>& x, 
		valarray<_Ty>& y, valarray<_Ty>& dy, valarray<_Ty>& ddy, 
			valarray<_Ty>& t, valarray<_Ty>& z, valarray<_Ty>& dz, 
													valarray<_Ty>& ddz);


//二元三点插值
template <class _Ty>
_Ty Interpolation2Variable3Points(valarray<_Ty>& x, valarray<_Ty>& y, 
											matrix<_Ty> z, _Ty u, _Ty v);

//二元全区间插值
template <class _Ty>
_Ty Interpolation2VariableWholeInterval(valarray<_Ty>& x, 
						valarray<_Ty>& y, matrix<_Ty> z, _Ty u, _Ty v);

#include "Interpolation.inl"		//类及相关函数的定义头文件

#endif		//_INTERPOLATION_H

