//NonLinearEquation.h	非线性方程(组)求解函数(方法)声明
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _NONLINEAREQUATION_H
#define _NONLINEAREQUATION_H

#include <comm.h>				//公共头文件
#include <math.h>				//数学头文件
#include <random.h>				//随机数头文件
#include <matrix.h>				//模板类matrix标准头文件
#include <LinearEquation.h>		//线性方程(组)求解头文件
#include <EigenvalueVector.h>	//计算特征值特征向量头文件

//template <class _Ty = float>

//二分法搜索方程f(x)=0在区间[a,b]内的全部实根
template <class _Ty>
inline size_t 
RootHalves(_Ty a, _Ty b, _Ty step, _Ty eps, valarray<_Ty>& x, size_t m);

//牛顿(Newton)法求解非线性方程一个实根
template <class _Ty>
inline int 
RootNewton( _Ty& x, _Ty eps, size_t js);

//埃特金(Aitken)法求解非线性方程一个实根
template <class _Ty>
inline int 
RootAitken(_Ty& x, _Ty eps, size_t js);

//连分式(Fraction)法求解非线性方程一个实根
template <class _Ty>
inline int 
RootFraction(_Ty& x, _Ty eps);

//QR法求代数方程全部根
template <class _Ty, class _Tz>
inline int 
RootQR(valarray<_Ty>& a, valarray<_Tz>& x, _Ty eps, size_t jt);

//牛顿下山(NewtonHillDown)法求解实系数代数方程全部根(实根和复根)
template <class _Ty, class _Tz>
inline int 
RootNewtonHillDown(valarray<_Ty>& a, valarray<_Tz>& cx);

//牛顿下山(NewtonHillDown)法求解复系数代数方程全部根(实根和复根)
//重载RootNewtonHillDown()
template <class _Ty>
inline int 
RootNewtonHillDown(complex<_Ty> a[], valarray< complex<_Ty> >& cx);

//梯度(Gradient)法(最速下降)求解非线性方程组一组实根
template <class _Ty>
inline int 
RootGradient(_Ty eps, valarray<_Ty>& x, size_t js);

//拟牛顿(QuasiNewton)法求解非线性方程组一组实根
template <class _Ty>
inline int 
RootQuasiNewton(_Ty eps, _Ty t, _Ty h, valarray<_Ty>& x, int k);

//非线性方程组最小二乘解的广义逆法
template <class _Ty>
int RootLeastSquareGeneralizedInverse(int m, _Ty eps1, _Ty eps2, 
										valarray<_Ty>& x, int ka);

//蒙特卡洛(MonteCarlo)法求解非线性方程f(x)=0的一个实根
//f(x)的自变量为与系数都为实数
template <class _Ty>
inline void 
RootMonteCarloComplex(_Ty& x, _Ty b, int m, _Ty eps);

//蒙特卡洛(MonteCarlo)法求解实(复)函数方程f(x)=0的一个复根
//f(x)的自变量为复数，或自变量与系数都为复数(不能都为实数)
template <class _Tz, class _Ty>
inline void 
RootMonteCarloComplex(_Tz& cxy, _Ty b, int m, _Ty eps);

//蒙特卡洛(MonteCarlo)法求解非线性方程组F(x)=0的一组实根
//f(x)的自变量为与系数都为实数
template <class _Ty>
inline void 
RootMonteCarloGroupReal(valarray<_Ty>& x, _Ty b, int m, _Ty eps);

#include <NonLinearEquation.inl>

#endif //_NONLINEAREQUATION_H

