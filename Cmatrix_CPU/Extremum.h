// Extremum.h			计算极值函数头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.9.6.

#ifndef _EXTREMUM_H		//避免多次编译
#define _EXTREMUM_H

#include <comm.h>		//公共头文件
#include <random.h>		//随机数头文件
#include <Matrix.h>		//矩阵类及相关函数头文件
#include <valarray>		//模板类valarray的标准头文件
using namespace std;	//名字空间

//一维极值连分式法
template <class _Ty>
void ExtremumFraction1D(valarray<_Ty>& x, _Ty eps, 
										int k, valarray<int>& js);

//n维极值连分式法
template <class _Ty>
void ExtremumFractionND(valarray<_Ty>& x, _Ty eps, 
										int k, valarray<int>& js);

//线性规划
template <class _Ty>
int ExtremumLinePrograming(matrix<_Ty>& a, valarray<_Ty>& b, 
								valarray<_Ty>& c, valarray<_Ty>& x);

//n维极值单形调优法
template <class _Ty>
int ExtremumSimplexND(_Ty d, _Ty u, _Ty v, valarray<_Ty>& x,  
				_Ty eps, int k, matrix<_Ty>& xx, valarray<_Ty>& f);

//n维极值复形调优法
template <class _Ty>
int ExtremumComplexND(int m, valarray<_Ty>& a, valarray<_Ty>& b,
		_Ty alpha, _Ty eps, valarray<_Ty>& x, matrix<_Ty>& xx, int k);

//确定1维函数极小值点所在区间
template <class _Ty>
void MinimizerInterval(_Ty& ax, _Ty& bx, _Ty& cx, _Ty& fa, _Ty& fb, _Ty& fc);

//确定1维函数极小值点所在区间(重载)
template <class _Ty>
int MinimizerInterval(_Ty& ax, _Ty& bx, _Ty& cx, int& sign, 
									_Ty b, _Ty& fa, _Ty& fb, _Ty& fc);
//1维函数极小值的黄金分割法
template <class _Ty>
_Ty ExtremumGold1D(_Ty ax, _Ty bx, _Ty cx, int sign, _Ty tol, _Ty& xmin);

//1维函数极小值的不用导数布伦特(Brent)法
template <class _Ty>
_Ty ExtremumBrentNonDerivative1D(_Ty ax, _Ty bx, _Ty cx, int sign, 
												_Ty tol, _Ty& xmin);

#include "Extremum.inl"		//类及相关函数的定义头文件

#endif	// _EXTREMUM_H

