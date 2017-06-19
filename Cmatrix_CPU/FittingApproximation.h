//FittingApproximation.h	拟合与逼近头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _FITTINGAPPROXIMATION_H		//避免多次编译
#define _FITTINGAPPROXIMATION_H

#include <valarray>			//模板类数组的标准头文件
#include <Matrix.h>			//模板类矩阵头文件
#include <comm.h>			//公共头文件

using namespace std;		//名字空间

//template <class _Ty = float>

//最小二乘曲线拟合 
template <class _Ty>
void FitCurveLeastSquares(valarray<_Ty>& x, valarray<_Ty>& y, 
								valarray<_Ty>& a, valarray<_Ty>& dt);

//切比雪夫曲线拟合
template <class _Ty>
void FitCurveChebyshev(valarray<_Ty>& x, valarray<_Ty>& y, valarray<_Ty>& a);

//最佳一致逼近多项式里米兹法
template <class _Ty>
void ApproximationRemez(_Ty a, _Ty b, valarray<_Ty>& p,_Ty eps);

//矩形域的最小二乘曲面拟合
template <class _Ty>
void FitSurfaceLeastSquares(valarray<_Ty>& x, valarray<_Ty>& y, 
				matrix<_Ty>& z, matrix<_Ty>& a,	valarray<_Ty>& dt);

#include "FittingApproximation.inl"		//类及相关函数的定义头文件

#endif		// _FITTINGAPPROXIMATION_H

