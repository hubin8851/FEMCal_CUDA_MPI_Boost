// Extremum.h			���㼫ֵ����ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.9.6.

#ifndef _EXTREMUM_H		//�����α���
#define _EXTREMUM_H

#include <comm.h>		//����ͷ�ļ�
#include <random.h>		//�����ͷ�ļ�
#include <Matrix.h>		//�����༰��غ���ͷ�ļ�
#include <valarray>		//ģ����valarray�ı�׼ͷ�ļ�
using namespace std;	//���ֿռ�

//һά��ֵ����ʽ��
template <class _Ty>
void ExtremumFraction1D(valarray<_Ty>& x, _Ty eps, 
										int k, valarray<int>& js);

//nά��ֵ����ʽ��
template <class _Ty>
void ExtremumFractionND(valarray<_Ty>& x, _Ty eps, 
										int k, valarray<int>& js);

//���Թ滮
template <class _Ty>
int ExtremumLinePrograming(matrix<_Ty>& a, valarray<_Ty>& b, 
								valarray<_Ty>& c, valarray<_Ty>& x);

//nά��ֵ���ε��ŷ�
template <class _Ty>
int ExtremumSimplexND(_Ty d, _Ty u, _Ty v, valarray<_Ty>& x,  
				_Ty eps, int k, matrix<_Ty>& xx, valarray<_Ty>& f);

//nά��ֵ���ε��ŷ�
template <class _Ty>
int ExtremumComplexND(int m, valarray<_Ty>& a, valarray<_Ty>& b,
		_Ty alpha, _Ty eps, valarray<_Ty>& x, matrix<_Ty>& xx, int k);

//ȷ��1ά������Сֵ����������
template <class _Ty>
void MinimizerInterval(_Ty& ax, _Ty& bx, _Ty& cx, _Ty& fa, _Ty& fb, _Ty& fc);

//ȷ��1ά������Сֵ����������(����)
template <class _Ty>
int MinimizerInterval(_Ty& ax, _Ty& bx, _Ty& cx, int& sign, 
									_Ty b, _Ty& fa, _Ty& fb, _Ty& fc);
//1ά������Сֵ�Ļƽ�ָ
template <class _Ty>
_Ty ExtremumGold1D(_Ty ax, _Ty bx, _Ty cx, int sign, _Ty tol, _Ty& xmin);

//1ά������Сֵ�Ĳ��õ���������(Brent)��
template <class _Ty>
_Ty ExtremumBrentNonDerivative1D(_Ty ax, _Ty bx, _Ty cx, int sign, 
												_Ty tol, _Ty& xmin);

#include "Extremum.inl"		//�༰��غ����Ķ���ͷ�ļ�

#endif	// _EXTREMUM_H

