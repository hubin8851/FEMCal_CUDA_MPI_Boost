// Integral.h		数值积分头文件
// Ver 1.0.0.0
// 版权所有(C) 2002
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _INTEGRAL_H		//避免多次编译
#define _INTEGRAL_H

#include <comm.h>		//公共头文件
#include <random.h>		//随机算法头文件
#include <valarray>		//模板类valarray的标准头文件

using namespace std;	//名字空间


//变步长梯形法求积
template <class _Ty>
_Ty IntegralTrapezia(_Ty a, _Ty b, _Ty eps);

//变步长辛卜生法求积
template <class _Ty>
_Ty IntegralSimpson1D(_Ty a, _Ty b, _Ty eps);

//自适应梯形法求积
template <class _Ty>
_Ty IntegralTrapeziaSelfAdapt(_Ty a, _Ty b, _Ty eps, _Ty d);

//龙贝格法求积
template <class _Ty>
_Ty IntegralRomberg(_Ty a, _Ty b, _Ty eps);

//一维连分式法求积
template <class _Ty>
_Ty IntegralFraction1D(_Ty a, _Ty b, _Ty eps);

//高振荡函数法求积
template <class _Ty>
void IntegralSurge(_Ty a, _Ty b, int m, valarray<_Ty>& fa, 
								valarray<_Ty>& fb, valarray<_Ty>& s);
//勒让德-高斯法求积
template <class _Ty>
_Ty IntegralLegendreGauss(_Ty a, _Ty b, _Ty eps);

//拉盖尔-高斯法求积
template <class _Ty>
void IntegralLaguerreGauss(_Ty& dValue);

//埃尔米特-高斯法求积
template <class _Ty>
void IntegralHermiteGauss(_Ty& dValue);

//切比雪夫法求积
template <class _Ty>
_Ty IntegralChebyshev(_Ty a, _Ty b, _Ty eps);

//蒙特卡洛法求积
template <class _Ty >
_Ty IntegralMonteCarlo1D(_Ty a, _Ty b);

//二重变步长辛卜生法求积
template <class _Ty>
_Ty IntegralSimpson2D(_Ty a, _Ty b, _Ty eps);

//多重高斯法求积
template <class _Ty >
_Ty IntegralGaussMD(valarray<_Ty>& js);

//二重连分式法求积
template <class _Ty >
_Ty IntegralFraction2D(_Ty a, _Ty b, _Ty eps);

//多重蒙特卡洛法求积
template <class _Ty >
_Ty IntegralMonteCarlo2D(valarray<_Ty>& a, valarray<_Ty>& b);

#include "Integral.inl"		//类及相关函数的定义头文件

#endif		// _INTEGRAL_H

