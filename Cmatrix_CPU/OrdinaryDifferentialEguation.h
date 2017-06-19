// OrdinaryDifferentialEguation.h	��ⳣ΢�ַ���(��)ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _ORDINARYDIFFERENTIALEGUATION_H		//�����α���
#define _ORDINARYDIFFERENTIALEGUATION_H

#include <LinearEquation.h>	//���Է�����ͷ�ļ�
#include <valarray>			//ģ��������ı�׼ͷ�ļ�
#include <Matrix.h>			//ģ�������ͷ�ļ�
#include <comm.h>			//����ͷ�ļ�

using namespace std;		//���ֿռ�

//������ȫ������ֵ�ŷ����
template <class _Ty>
void ODE_EulerContentStep(_Ty t, valarray<_Ty>& y,  
										_Ty h, int k, matrix<_Ty>& z);

//�䲽������ŷ����
template <class _Ty>
void ODE_EulerVariationalStep(_Ty t, _Ty h, valarray<_Ty>& y, _Ty eps);

//ȫ���䶨��������ά�ݷ�
template <class _Ty>
void ODE_WittyContentStep(_Ty t, valarray<_Ty>& y,  
										_Ty h, int k, matrix<_Ty>& z);

//ȫ���䶨������������-������
template <class _Ty>
void ODE_RungeKuttaContentStep(_Ty t, valarray<_Ty>& y,  
										_Ty h, int k, matrix<_Ty>& z);

//�䲽����������-������
template <class _Ty>
void ODE_RungeKuttaVariationalStep(_Ty t, _Ty h, 
										valarray<_Ty>& y, _Ty eps);

//�䲽������һ��������
template <class _Ty>
void ODE_GillVariationalStep(_Ty t, _Ty h, valarray<_Ty>& y,
										_Ty eps, valarray<_Ty>& q);

//ȫ����䲽�����ֻ�����
template <class _Ty>
void ODE_GillVariationalStepInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z);

//�䲽������Ĭɭ��
template <class _Ty>
void ODE_MersonVariationalStep(_Ty t, _Ty h, int n, valarray<_Ty>& y,
										_Ty eps, int k, matrix<_Ty>& z);

//����һ������ʽ��
template <class _Ty>
void ODE_Fraction(_Ty t, _Ty h, valarray<_Ty>& y, _Ty eps);

//ȫ�����������ʽ��
template <class _Ty>
void ODE_FractionInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z);

//ȫ�������˫�߷�
template <class _Ty>
void ODE_BothSidesInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z);

//ȫ������ְ���ķ˹Ԥ��-У����
template <class _Ty>
void ODE_AdamsInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z);

//ȫ������ֹ�����
template <class _Ty>
void ODE_HammingInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z);

//����һ�������ɷ�
template <class _Ty>
void ODE_Treanor(_Ty t, _Ty h, valarray<_Ty>& y);

//ȫ������������ɷ�
template <class _Ty>
void ODE_TreanorInterval(_Ty t, _Ty h, valarray<_Ty>& y,
											int k, matrix<_Ty>& z);

//���ָ��Է����鼪����
template <class _Ty>
int ODE_Gear(_Ty a, _Ty b, _Ty hmin, _Ty hmax, _Ty h, _Ty eps, 
		 valarray<_Ty>& y0, int k, valarray<_Ty>& t, matrix<_Ty>& z);

//����΢�ַ��̱�ֵ������ֵ�ⷨ
template <class _Ty>
void ODE_LinearBoundaryValude(_Ty a, _Ty b, _Ty ya, _Ty yb, int n, 
													valarray<_Ty>& y);

#include "OrdinaryDifferentialEguation.inl"		//�༰��غ����Ķ���ͷ�ļ�

#endif	// _ORDINARYDIFFERENTIALEGUATION_H

