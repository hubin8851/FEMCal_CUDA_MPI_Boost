// EigenvalueVector.h		计算特征值特征向量头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002

// 最后修改: 2002.8.10

#ifndef _EIGENVALUEVECTOR_H		//避免多次编译
#define _EIGENVALUEVECTOR_H

#include "Matrix.h"		//矩阵类及相关函数等的定义
#include <comm.h>		//公共头文件
#include <math.h>		//数学头文件

using namespace std;	//名字空间

//约化对称阵为对称三对角阵的豪斯荷尔德变换法
template <class _Ty>
int HouseholderTransform(matrix<_Ty>& a, matrix<_Ty>& q, 
									valarray<_Ty>& b, valarray<_Ty>& c);

//实对称三角阵全部特征值及特征向量QR法
template <class _Ty>
int EigenvalueVectorRealTriangleQR(valarray<_Ty>& b, 
					valarray<_Ty>& c, matrix<_Ty>& q, _Ty eps, int l);

//约化一般实矩阵为赫申伯格阵的初等相似变换法
//矩阵类型应是浮点型
template <class _Ty>
int HessenbergTransform(matrix<_Ty>& a);

//求赫申伯格阵全部特征值QR法
template <class _Ty>
int EigenvalueVectorHessenbergQR(matrix<_Ty>& a,  
						valarray<complex<_Ty> >& uv, _Ty eps, int jt);

//实对称阵特征值及特征向量雅可比法
template <class _Ty>
int EigenvalueVectorRealSymmetryJacobi(matrix<_Ty>& a,  
										matrix<_Ty>& v, _Ty eps, int jt);

//实对称阵特征值及特征向量雅可比过关法
template <class _Ty>
int EigenvalueVectorRealSymmetryJacobiB(matrix<_Ty>& a,  
												matrix<_Ty>& v, _Ty eps);

#include "EigenvalueVector.inl"		//类及相关函数的定义头文件

#endif		// _EIGENVALUEVECTOR_H

