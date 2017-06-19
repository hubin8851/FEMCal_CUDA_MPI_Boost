// LinearEquation.h		线性方程(组)求解函数(方法)声明
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31

#ifndef _LINEAREQUATION_H
#define _LINEAREQUATION_H

#include <comm.h>		//公共头文件
#include <valarray>		//模板类valarray的标准头文件
#include <Matrix.h>		//矩阵类头文件

//全选主元高斯消元法
template <class _Ty>		//奇异返回0，否则正常
int LE_TotalChoiceGauss(matrix<_Ty>& a, valarray<_Ty>& b);

//全选主元高斯-约当消元法
template <class _Ty>		//奇异返回0，否则正常
int LE_TotalChoiceGaussJordan(matrix<_Ty>& a, matrix<_Ty>& b);

//三对角方程组的追赶法
template <class _Ty>
int LE_TridiagonalEquationGauss(valarray<_Ty>& b, valarray<_Ty>& d);

//一般带型方程组求解
template <class _Ty>
int LE_StrapEquationGauss(matrix<_Ty>& b, matrix<_Ty>& d, int l, int il);

//对称方程组的分解法
template <class _Ty>
int LE_SymmetryEquation(matrix<_Ty>& a, matrix<_Ty>& c);

//对称正定方程组的平方根法
template <class _Ty>
int LE_SymmetryRegularEuationSquareRoot(matrix<_Ty>& a, matrix<_Ty>& d);

//大型稀疏方程组全选主元高斯-约当(Gauss-Jordan)法
template <class _Ty>
int LE_SparseEuationTotalChoiceGaussJordan(matrix<_Ty>& a, valarray<_Ty>& b);

//托伯利兹(Toeplitz)方程组的列文逊(Levinson)法
template <class _Ty>
int LE_ToeplitzEuationLevinson(valarray<_Ty>& t, valarray<_Ty>& b, valarray<_Ty>& x);

//高斯-赛德尔迭代法
template <class _Ty>
int LE_GaussSeidelIteration(matrix<_Ty>& a, valarray<_Ty>& b, valarray<_Ty>& x, _Ty eps);

//对称正定方程组的共轭梯度法
template <class _Ty>
int LE_SymmetryRegularEuationConjugateGradient(matrix<_Ty>& a, valarray<_Ty>& b, _Ty eps, valarray<_Ty>& x);

//线性最小二乘的豪斯荷尔德变换法
template <class _Ty>
int LE_LinearLeastSquareHouseholder(matrix<_Ty>& a, valarray<_Ty>& b, matrix<_Ty>& q);

//线性最小二乘的广义逆法
template <class _Ty>
int LE_LinearLeastSquareGeneralizedInverse(matrix<_Ty>& a, valarray<_Ty>& b, matrix<_Ty>& x, matrix<_Ty>& aa, _Ty eps, matrix<_Ty>& u, matrix<_Ty>& v);

//病态方程组求解
template <class _Ty>
int LE_IllConditionedEquation(matrix<_Ty>& a, valarray<_Ty>& b, _Ty eps, matrix<_Ty>& x);

#include <LinearEquation.inl>

#endif //_LINEAREQUATION_H

