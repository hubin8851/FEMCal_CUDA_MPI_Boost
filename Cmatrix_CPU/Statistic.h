//Statistic.h	���ݴ�����ع����ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _STATISTIC_H		//�����α���
#define _STATISTIC_H

#include <iostream>			//ģ�������������iostream��׼ͷ�ļ�
#include <valarray>			//ģ����valarray�ı�׼ͷ�ļ�
#include <Matrix.h>			//ģ����Matrixͷ�ļ�
#include <LinearEquation.h>	//���Է�����ͷ�ļ�
#include <comm.h>			//����commͷ�ļ�

using namespace std;		//���ֿռ�

//template <class _Ty = float>
//�����������
template <class _Ty>
void StatRandomSample(valarray<_Ty>& x, _Ty x0, _Ty h, int l, 
			valarray<_Ty>& dt, valarray<int>& g, valarray<int>& q);

//һԪ���Իع����
template <class _Ty>
void LinearRegression1D(valarray<_Ty>& x, valarray<_Ty>& y, 
								valarray<_Ty>& a, valarray<_Ty>& dt);

//nԪ���Իع����
template <class _Ty>
void LinearRegressionND(matrix<_Ty>& x, valarray<_Ty>& y, 
			valarray<_Ty>& a, valarray<_Ty>& dt, valarray<_Ty>& v);

//�𲽻ع����
template <class _Ty>
void StepwiseRegression(matrix<_Ty>& x, _Ty f1, _Ty f2, 
		_Ty eps, valarray<_Ty>& xx, valarray<_Ty>& b, valarray<_Ty>& v, 
				valarray<_Ty>& s, valarray<_Ty>& dt, valarray<_Ty>& ye, 
									valarray<_Ty>& yr, matrix<_Ty>& r);
//������������
template <class _Ty>
void HalfLogarithmCorrelation(valarray<_Ty>& x, valarray<_Ty>& y, 
											_Ty t, valarray<_Ty>& a);

//�����������
template <class _Ty>
void LogarithmCorrelation(valarray<_Ty>& x, 
								valarray<_Ty>& y, valarray<_Ty>& a);

#include "Statistic.inl"		//�༰��غ����Ķ���ͷ�ļ�

#endif		// _STATISTIC_H