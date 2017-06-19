//Transform.h	��ѧ�任ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _TRANSFORM_H		//�����α���
#define _TRANSFORM_H

#include <iostream>			//ģ���������������׼ͷ�ļ�
#include <valarray>			//ģ��������ı�׼ͷ�ļ�
#include <Matrix.h>			//ģ�������ͷ�ļ�
#include <comm.h>			//����ͷ�ļ�

using namespace std;		//���ֿռ�

//����Ҷ�����ƽ�
template <class _Ty>
void FourierSeriesApproach(valarray<_Ty>& f, 
								valarray<_Ty>& a, valarray<_Ty>& b);

//���ٸ���Ҷ�任
template <class _Ty>
void FourierTransform(valarray<complex<_Ty> >& pp, 
						valarray<complex<_Ty> >& ff, int l, int il);

//������ʲ�任
template <class _Ty>
void WalshTransform(valarray<_Ty>& p, valarray<_Ty>& x);

//�������ƽ��(�������)
template <class _Ty>
void Smooth5_3(valarray<_Ty>& y, valarray<_Ty>& yy);

//��ɢ�������ϵͳ�Ŀ������˲�
template <class _Ty>
int SievingKalman(matrix<_Ty>& f, matrix<_Ty>& q, matrix<_Ty>& r, 
					matrix<_Ty>& h, matrix<_Ty>& y, matrix<_Ty>& x, 
									matrix<_Ty>& p, matrix<_Ty>& g);

//a-b-r�˲�
template <class _Ty>
void SievingABR(valarray<_Ty>& x, _Ty t, _Ty a, _Ty b,
											_Ty r, valarray<_Ty>& y);


#include "Transform.inl"	//�༰��غ����Ķ���ͷ�ļ�

#endif	// _TRANSFORM_H
