// OrdinaryDifferentialEguation.h	求解常微分方程(组)头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _ORDINARYDIFFERENTIALEGUATION_H		//避免多次编译
#define _ORDINARYDIFFERENTIALEGUATION_H

#include <LinearEquation.h>	//线性方程组头文件
#include <valarray>			//模板类数组的标准头文件
#include <Matrix.h>			//模板类矩阵头文件
#include <comm.h>			//公共头文件

using namespace std;		//名字空间

//定步长全区间积分的欧拉法
template <class _Ty>
void ODE_EulerContentStep(_Ty t, valarray<_Ty>& y,  
										_Ty h, int k, matrix<_Ty>& z);

//变步长积分欧拉法
template <class _Ty>
void ODE_EulerVariationalStep(_Ty t, _Ty h, valarray<_Ty>& y, _Ty eps);

//全区间定步长积分维梯法
template <class _Ty>
void ODE_WittyContentStep(_Ty t, valarray<_Ty>& y,  
										_Ty h, int k, matrix<_Ty>& z);

//全区间定步长积分龙格-库塔法
template <class _Ty>
void ODE_RungeKuttaContentStep(_Ty t, valarray<_Ty>& y,  
										_Ty h, int k, matrix<_Ty>& z);

//变步长积分龙格-库塔法
template <class _Ty>
void ODE_RungeKuttaVariationalStep(_Ty t, _Ty h, 
										valarray<_Ty>& y, _Ty eps);

//变步长积分一步基尔法
template <class _Ty>
void ODE_GillVariationalStep(_Ty t, _Ty h, valarray<_Ty>& y,
										_Ty eps, valarray<_Ty>& q);

//全区间变步长积分基尔法
template <class _Ty>
void ODE_GillVariationalStepInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z);

//变步长积分默森法
template <class _Ty>
void ODE_MersonVariationalStep(_Ty t, _Ty h, int n, valarray<_Ty>& y,
										_Ty eps, int k, matrix<_Ty>& z);

//积分一步连分式法
template <class _Ty>
void ODE_Fraction(_Ty t, _Ty h, valarray<_Ty>& y, _Ty eps);

//全区间积分连分式法
template <class _Ty>
void ODE_FractionInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z);

//全区间积分双边法
template <class _Ty>
void ODE_BothSidesInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z);

//全区间积分阿当姆斯预报-校正法
template <class _Ty>
void ODE_AdamsInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z);

//全区间积分哈明法
template <class _Ty>
void ODE_HammingInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z);

//积分一步特雷纳法
template <class _Ty>
void ODE_Treanor(_Ty t, _Ty h, valarray<_Ty>& y);

//全区间积分特雷纳法
template <class _Ty>
void ODE_TreanorInterval(_Ty t, _Ty h, valarray<_Ty>& y,
											int k, matrix<_Ty>& z);

//积分刚性方程组吉尔法
template <class _Ty>
int ODE_Gear(_Ty a, _Ty b, _Ty hmin, _Ty hmax, _Ty h, _Ty eps, 
		 valarray<_Ty>& y0, int k, valarray<_Ty>& t, matrix<_Ty>& z);

//二阶微分方程边值问题数值解法
template <class _Ty>
void ODE_LinearBoundaryValude(_Ty a, _Ty b, _Ty ya, _Ty yb, int n, 
													valarray<_Ty>& y);

#include "OrdinaryDifferentialEguation.inl"		//类及相关函数的定义头文件

#endif	// _ORDINARYDIFFERENTIALEGUATION_H

