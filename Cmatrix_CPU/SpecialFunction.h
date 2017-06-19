//SpecialFunction.h		特殊函数头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _SPECIALFUNCTION_H		//避免多次编译
#define _SPECIALFUNCTION_H

#include <iostream>			//模板类输入输出流标准头文件
#include <comm.h>			//公共头文件

using namespace std;		//名字空间

//伽马函数
template <class _Ty>
_Ty GammaFunction(_Ty x);

//不完全伽马函数
template <class _Ty>
_Ty IncompleteGammaFunction(_Ty a, _Ty x);

//误差函数
template <class _Ty>
_Ty ErrorFunction(_Ty x);

//第一类整数阶贝塞尔函数
template <class _Ty>
_Ty IntegerBessel1stFunction(int n, _Ty x);

//第二类整数阶贝塞尔函数
template <class _Ty>
_Ty IntegerBessel2ndFunction(int n, _Ty x);

//变形第一类整数阶贝塞尔函数
template <class _Ty>
_Ty TransformativeIntegerBessel1stFunction(int n,_Ty x);

//变形第二类整数阶贝塞尔函数
template <class _Ty>
_Ty TransformativeIntegerBessel2ndFunction(int n, _Ty x);

//不完全贝塔函数
template <class _Ty>
_Ty IncompleteBetaFunction(_Ty a, _Ty b, _Ty x);

//正态分布函数
template <class _Ty>
_Ty NormalDistributionFunction(_Ty a, _Ty d, _Ty x);

//t-分布函数
template <class _Ty>
_Ty tDistributionFunction(_Ty t, int n);

//X^2-分布函数
template <class _Ty>
_Ty X2DistributionFunction(_Ty x, int n);

//F-分布函数
template <class _Ty>
_Ty FDistributionFunction(_Ty f, int n1, int n2);

//正弦积分
template <class _Ty>
_Ty SineIntegralFunction(_Ty x);

//余弦积分
template <class _Ty>
_Ty CosineIntegralFunction(_Ty x);

//指数积分
template <class _Ty>
_Ty ExponentIntegralFunction(_Ty x);

//第一类椭圆积分
template <class _Ty>
_Ty Ellipse1stIntegral(_Ty k, _Ty f);

//第二类椭圆积分
template <class _Ty>
_Ty Ellipse2ndIntegral(_Ty k, _Ty f);

#include "SpecialFunction.inl"		//类及相关函数的定义头文件

#endif		// _SPECIALFUNCTION_H

