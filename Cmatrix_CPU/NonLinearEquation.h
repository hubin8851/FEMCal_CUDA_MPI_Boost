//NonLinearEquation.h	�����Է���(��)��⺯��(����)����
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _NONLINEAREQUATION_H
#define _NONLINEAREQUATION_H

#include <comm.h>				//����ͷ�ļ�
#include <math.h>				//��ѧͷ�ļ�
#include <random.h>				//�����ͷ�ļ�
#include <matrix.h>				//ģ����matrix��׼ͷ�ļ�
#include <LinearEquation.h>		//���Է���(��)���ͷ�ļ�
#include <EigenvalueVector.h>	//��������ֵ��������ͷ�ļ�

//template <class _Ty = float>

//���ַ���������f(x)=0������[a,b]�ڵ�ȫ��ʵ��
template <class _Ty>
inline size_t 
RootHalves(_Ty a, _Ty b, _Ty step, _Ty eps, valarray<_Ty>& x, size_t m);

//ţ��(Newton)���������Է���һ��ʵ��
template <class _Ty>
inline int 
RootNewton( _Ty& x, _Ty eps, size_t js);

//���ؽ�(Aitken)���������Է���һ��ʵ��
template <class _Ty>
inline int 
RootAitken(_Ty& x, _Ty eps, size_t js);

//����ʽ(Fraction)���������Է���һ��ʵ��
template <class _Ty>
inline int 
RootFraction(_Ty& x, _Ty eps);

//QR�����������ȫ����
template <class _Ty, class _Tz>
inline int 
RootQR(valarray<_Ty>& a, valarray<_Tz>& x, _Ty eps, size_t jt);

//ţ����ɽ(NewtonHillDown)�����ʵϵ����������ȫ����(ʵ���͸���)
template <class _Ty, class _Tz>
inline int 
RootNewtonHillDown(valarray<_Ty>& a, valarray<_Tz>& cx);

//ţ����ɽ(NewtonHillDown)����⸴ϵ����������ȫ����(ʵ���͸���)
//����RootNewtonHillDown()
template <class _Ty>
inline int 
RootNewtonHillDown(complex<_Ty> a[], valarray< complex<_Ty> >& cx);

//�ݶ�(Gradient)��(�����½�)�������Է�����һ��ʵ��
template <class _Ty>
inline int 
RootGradient(_Ty eps, valarray<_Ty>& x, size_t js);

//��ţ��(QuasiNewton)���������Է�����һ��ʵ��
template <class _Ty>
inline int 
RootQuasiNewton(_Ty eps, _Ty t, _Ty h, valarray<_Ty>& x, int k);

//�����Է�������С���˽�Ĺ����淨
template <class _Ty>
int RootLeastSquareGeneralizedInverse(int m, _Ty eps1, _Ty eps2, 
										valarray<_Ty>& x, int ka);

//���ؿ���(MonteCarlo)���������Է���f(x)=0��һ��ʵ��
//f(x)���Ա���Ϊ��ϵ����Ϊʵ��
template <class _Ty>
inline void 
RootMonteCarloComplex(_Ty& x, _Ty b, int m, _Ty eps);

//���ؿ���(MonteCarlo)�����ʵ(��)��������f(x)=0��һ������
//f(x)���Ա���Ϊ���������Ա�����ϵ����Ϊ����(���ܶ�Ϊʵ��)
template <class _Tz, class _Ty>
inline void 
RootMonteCarloComplex(_Tz& cxy, _Ty b, int m, _Ty eps);

//���ؿ���(MonteCarlo)���������Է�����F(x)=0��һ��ʵ��
//f(x)���Ա���Ϊ��ϵ����Ϊʵ��
template <class _Ty>
inline void 
RootMonteCarloGroupReal(valarray<_Ty>& x, _Ty b, int m, _Ty eps);

#include <NonLinearEquation.inl>

#endif //_NONLINEAREQUATION_H

