//Statistic.h	数据处理与回归分析头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _STATISTIC_H		//避免多次编译
#define _STATISTIC_H

#include <iostream>			//模板类输入输出流iostream标准头文件
#include <valarray>			//模板类valarray的标准头文件
#include <Matrix.h>			//模板类Matrix头文件
#include <LinearEquation.h>	//线性方程组头文件
#include <comm.h>			//公共comm头文件

using namespace std;		//名字空间

//template <class _Ty = float>
//随机样本分析
template <class _Ty>
void StatRandomSample(valarray<_Ty>& x, _Ty x0, _Ty h, int l, 
			valarray<_Ty>& dt, valarray<int>& g, valarray<int>& q);

//一元线性回归分析
template <class _Ty>
void LinearRegression1D(valarray<_Ty>& x, valarray<_Ty>& y, 
								valarray<_Ty>& a, valarray<_Ty>& dt);

//n元线性回归分析
template <class _Ty>
void LinearRegressionND(matrix<_Ty>& x, valarray<_Ty>& y, 
			valarray<_Ty>& a, valarray<_Ty>& dt, valarray<_Ty>& v);

//逐步回归分析
template <class _Ty>
void StepwiseRegression(matrix<_Ty>& x, _Ty f1, _Ty f2, 
		_Ty eps, valarray<_Ty>& xx, valarray<_Ty>& b, valarray<_Ty>& v, 
				valarray<_Ty>& s, valarray<_Ty>& dt, valarray<_Ty>& ye, 
									valarray<_Ty>& yr, matrix<_Ty>& r);
//半对数数据相关
template <class _Ty>
void HalfLogarithmCorrelation(valarray<_Ty>& x, valarray<_Ty>& y, 
											_Ty t, valarray<_Ty>& a);

//对数数据相关
template <class _Ty>
void LogarithmCorrelation(valarray<_Ty>& x, 
								valarray<_Ty>& y, valarray<_Ty>& a);

#include "Statistic.inl"		//类及相关函数的定义头文件

#endif		// _STATISTIC_H