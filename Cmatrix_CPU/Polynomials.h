//polynomials.h		多项式与连分式头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _POLYNOMIALS_H
#define _POLYNOMIALS_H

#include <comm.h>		 //公共部分头文件
#include <math.h>		 //数学头文件
#include <matrix.h>		 //矩阵类头文件

//求一维实(复)多项式值
template<class T, class U>
inline U 
PolyValueOneDim(valarray<T>& dCoff, size_t stNo, U dX);

//求一维多项式组值
template<class T, class V, class U>
inline void 
PolyValueOneDimGroup(valarray<T>& dCoff, valarray<V>& dX, valarray<U>& dValue);

//求二维实(复)多项式值
template<class T, class U>
inline U	//内联函数
PolyValueTwoDim(matrix<T>& dCoff, U dX, U dY);

//两一维多项式相乘
template<class T>
inline void
PolyMultip(valarray<T>& dCoffP, valarray<T>& dCoffQ, valarray<T>& dCoffS);

//两一维多项式除法
template<class T>
inline int
PolyDiv(valarray<T>& dCoffP, valarray<T>& dCoffQ, valarray<T>& dCoffS, valarray<T>& dCoffR);

// 计算连分式函数值
template<class T>
inline T
FractionValue(valarray<T>& dXpara, valarray<T>& dCoff, T dX);

#include "polynomials.inl"	//多项式及连分式相关函数定义头文件

#endif //_POLYNOMIALS_H

