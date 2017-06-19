// Integral.h		��ֵ����ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) 2002
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _INTEGRAL_H		//�����α���
#define _INTEGRAL_H

#include <comm.h>		//����ͷ�ļ�
#include <random.h>		//����㷨ͷ�ļ�
#include <valarray>		//ģ����valarray�ı�׼ͷ�ļ�

using namespace std;	//���ֿռ�


//�䲽�����η����
template <class _Ty>
_Ty IntegralTrapezia(_Ty a, _Ty b, _Ty eps);

//�䲽�������������
template <class _Ty>
_Ty IntegralSimpson1D(_Ty a, _Ty b, _Ty eps);

//����Ӧ���η����
template <class _Ty>
_Ty IntegralTrapeziaSelfAdapt(_Ty a, _Ty b, _Ty eps, _Ty d);

//���������
template <class _Ty>
_Ty IntegralRomberg(_Ty a, _Ty b, _Ty eps);

//һά����ʽ�����
template <class _Ty>
_Ty IntegralFraction1D(_Ty a, _Ty b, _Ty eps);

//���񵴺��������
template <class _Ty>
void IntegralSurge(_Ty a, _Ty b, int m, valarray<_Ty>& fa, 
								valarray<_Ty>& fb, valarray<_Ty>& s);
//���õ�-��˹�����
template <class _Ty>
_Ty IntegralLegendreGauss(_Ty a, _Ty b, _Ty eps);

//���Ƕ�-��˹�����
template <class _Ty>
void IntegralLaguerreGauss(_Ty& dValue);

//��������-��˹�����
template <class _Ty>
void IntegralHermiteGauss(_Ty& dValue);

//�б�ѩ�����
template <class _Ty>
_Ty IntegralChebyshev(_Ty a, _Ty b, _Ty eps);

//���ؿ��巨���
template <class _Ty >
_Ty IntegralMonteCarlo1D(_Ty a, _Ty b);

//���ر䲽�������������
template <class _Ty>
_Ty IntegralSimpson2D(_Ty a, _Ty b, _Ty eps);

//���ظ�˹�����
template <class _Ty >
_Ty IntegralGaussMD(valarray<_Ty>& js);

//��������ʽ�����
template <class _Ty >
_Ty IntegralFraction2D(_Ty a, _Ty b, _Ty eps);

//�������ؿ��巨���
template <class _Ty >
_Ty IntegralMonteCarlo2D(valarray<_Ty>& a, valarray<_Ty>& b);

#include "Integral.inl"		//�༰��غ����Ķ���ͷ�ļ�

#endif		// _INTEGRAL_H

