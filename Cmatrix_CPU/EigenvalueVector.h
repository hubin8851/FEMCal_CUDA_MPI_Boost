// EigenvalueVector.h		��������ֵ��������ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002

// ����޸�: 2002.8.10

#ifndef _EIGENVALUEVECTOR_H		//�����α���
#define _EIGENVALUEVECTOR_H

#include "Matrix.h"		//�����༰��غ����ȵĶ���
#include <comm.h>		//����ͷ�ļ�
#include <math.h>		//��ѧͷ�ļ�

using namespace std;	//���ֿռ�

//Լ���Գ���Ϊ�Գ����Խ���ĺ�˹�ɶ��±任��
template <class _Ty>
int HouseholderTransform(matrix<_Ty>& a, matrix<_Ty>& q, 
									valarray<_Ty>& b, valarray<_Ty>& c);

//ʵ�Գ�������ȫ������ֵ����������QR��
template <class _Ty>
int EigenvalueVectorRealTriangleQR(valarray<_Ty>& b, 
					valarray<_Ty>& c, matrix<_Ty>& q, _Ty eps, int l);

//Լ��һ��ʵ����Ϊ���겮����ĳ������Ʊ任��
//��������Ӧ�Ǹ�����
template <class _Ty>
int HessenbergTransform(matrix<_Ty>& a);

//����겮����ȫ������ֵQR��
template <class _Ty>
int EigenvalueVectorHessenbergQR(matrix<_Ty>& a,  
						valarray<complex<_Ty> >& uv, _Ty eps, int jt);

//ʵ�Գ�������ֵ�����������ſɱȷ�
template <class _Ty>
int EigenvalueVectorRealSymmetryJacobi(matrix<_Ty>& a,  
										matrix<_Ty>& v, _Ty eps, int jt);

//ʵ�Գ�������ֵ�����������ſɱȹ��ط�
template <class _Ty>
int EigenvalueVectorRealSymmetryJacobiB(matrix<_Ty>& a,  
												matrix<_Ty>& v, _Ty eps);

#include "EigenvalueVector.inl"		//�༰��غ����Ķ���ͷ�ļ�

#endif		// _EIGENVALUEVECTOR_H

