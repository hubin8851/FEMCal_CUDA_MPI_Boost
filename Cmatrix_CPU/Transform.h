//Transform.h	数学变换头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _TRANSFORM_H		//避免多次编译
#define _TRANSFORM_H

#include <iostream>			//模板类输入输出流标准头文件
#include <valarray>			//模板类数组的标准头文件
#include <Matrix.h>			//模板类矩阵头文件
#include <comm.h>			//公共头文件

using namespace std;		//名字空间

//傅里叶级数逼近
template <class _Ty>
void FourierSeriesApproach(valarray<_Ty>& f, 
								valarray<_Ty>& a, valarray<_Ty>& b);

//快速傅里叶变换
template <class _Ty>
void FourierTransform(valarray<complex<_Ty> >& pp, 
						valarray<complex<_Ty> >& ff, int l, int il);

//快速沃什变换
template <class _Ty>
void WalshTransform(valarray<_Ty>& p, valarray<_Ty>& x);

//五点三次平滑(曲线拟合)
template <class _Ty>
void Smooth5_3(valarray<_Ty>& y, valarray<_Ty>& yy);

//离散随机线性系统的卡尔曼滤波
template <class _Ty>
int SievingKalman(matrix<_Ty>& f, matrix<_Ty>& q, matrix<_Ty>& r, 
					matrix<_Ty>& h, matrix<_Ty>& y, matrix<_Ty>& x, 
									matrix<_Ty>& p, matrix<_Ty>& g);

//a-b-r滤波
template <class _Ty>
void SievingABR(valarray<_Ty>& x, _Ty t, _Ty a, _Ty b,
											_Ty r, valarray<_Ty>& y);


#include "Transform.inl"	//类及相关函数的定义头文件

#endif	// _TRANSFORM_H
