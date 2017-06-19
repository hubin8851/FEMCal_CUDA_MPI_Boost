//FittingApproximation.h	�����ƽ�ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _FITTINGAPPROXIMATION_H		//�����α���
#define _FITTINGAPPROXIMATION_H

#include <valarray>			//ģ��������ı�׼ͷ�ļ�
#include <Matrix.h>			//ģ�������ͷ�ļ�
#include <comm.h>			//����ͷ�ļ�

using namespace std;		//���ֿռ�

//template <class _Ty = float>

//��С����������� 
template <class _Ty>
void FitCurveLeastSquares(valarray<_Ty>& x, valarray<_Ty>& y, 
								valarray<_Ty>& a, valarray<_Ty>& dt);

//�б�ѩ���������
template <class _Ty>
void FitCurveChebyshev(valarray<_Ty>& x, valarray<_Ty>& y, valarray<_Ty>& a);

//���һ�±ƽ�����ʽ�����ȷ�
template <class _Ty>
void ApproximationRemez(_Ty a, _Ty b, valarray<_Ty>& p,_Ty eps);

//���������С�����������
template <class _Ty>
void FitSurfaceLeastSquares(valarray<_Ty>& x, valarray<_Ty>& y, 
				matrix<_Ty>& z, matrix<_Ty>& a,	valarray<_Ty>& dt);

#include "FittingApproximation.inl"		//�༰��غ����Ķ���ͷ�ļ�

#endif		// _FITTINGAPPROXIMATION_H

