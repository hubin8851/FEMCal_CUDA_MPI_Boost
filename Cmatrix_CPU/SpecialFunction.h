//SpecialFunction.h		���⺯��ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _SPECIALFUNCTION_H		//�����α���
#define _SPECIALFUNCTION_H

#include <iostream>			//ģ���������������׼ͷ�ļ�
#include <comm.h>			//����ͷ�ļ�

using namespace std;		//���ֿռ�

//٤����
template <class _Ty>
_Ty GammaFunction(_Ty x);

//����ȫ٤����
template <class _Ty>
_Ty IncompleteGammaFunction(_Ty a, _Ty x);

//����
template <class _Ty>
_Ty ErrorFunction(_Ty x);

//��һ�������ױ���������
template <class _Ty>
_Ty IntegerBessel1stFunction(int n, _Ty x);

//�ڶ��������ױ���������
template <class _Ty>
_Ty IntegerBessel2ndFunction(int n, _Ty x);

//���ε�һ�������ױ���������
template <class _Ty>
_Ty TransformativeIntegerBessel1stFunction(int n,_Ty x);

//���εڶ��������ױ���������
template <class _Ty>
_Ty TransformativeIntegerBessel2ndFunction(int n, _Ty x);

//����ȫ��������
template <class _Ty>
_Ty IncompleteBetaFunction(_Ty a, _Ty b, _Ty x);

//��̬�ֲ�����
template <class _Ty>
_Ty NormalDistributionFunction(_Ty a, _Ty d, _Ty x);

//t-�ֲ�����
template <class _Ty>
_Ty tDistributionFunction(_Ty t, int n);

//X^2-�ֲ�����
template <class _Ty>
_Ty X2DistributionFunction(_Ty x, int n);

//F-�ֲ�����
template <class _Ty>
_Ty FDistributionFunction(_Ty f, int n1, int n2);

//���һ���
template <class _Ty>
_Ty SineIntegralFunction(_Ty x);

//���һ���
template <class _Ty>
_Ty CosineIntegralFunction(_Ty x);

//ָ������
template <class _Ty>
_Ty ExponentIntegralFunction(_Ty x);

//��һ����Բ����
template <class _Ty>
_Ty Ellipse1stIntegral(_Ty k, _Ty f);

//�ڶ�����Բ����
template <class _Ty>
_Ty Ellipse2ndIntegral(_Ty k, _Ty f);

#include "SpecialFunction.inl"		//�༰��غ����Ķ���ͷ�ļ�

#endif		// _SPECIALFUNCTION_H

