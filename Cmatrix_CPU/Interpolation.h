// Interpolation.h		��ֵͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C)  ����(HE Yu) 2002

// ����޸�: 2002.5.31.

#ifndef _INTERPOLATION_H	//�����α���
#define _INTERPOLATION_H

#include <valarray>			//����ģ�����׼ͷ�ļ�
#include <Matrix.h>			//������ͷ�ļ�
#include <comm.h>			//����ͷ�ļ�
using namespace std;		//���ֿռ�

//һԪȫ���䲻�Ⱦ��ֵ
template <class _Ty>
_Ty Interpolation1VariableNotIsometry(valarray<_Ty>& x, 
											valarray<_Ty>& y, _Ty t);

//һԪȫ����Ⱦ��ֵ
template <class _Ty>
_Ty Interpolation1VariableIsometry(_Ty x0, _Ty h, valarray<_Ty>& y, _Ty t);

//һԪ���㲻�Ⱦ��ֵ
template <class _Ty>
_Ty Interpolation1Variable3PointsNotIsometry(valarray<_Ty>& x, 
											valarray<_Ty>& y, _Ty t);

//һԪ����Ⱦ��ֵ
template <class _Ty>
_Ty Interpolation1Variable3PointsIsometry(_Ty x0, _Ty h, 
											valarray<_Ty>& y, _Ty t);

//����ʽ���Ⱦ��ֵ
template <class _Ty>
_Ty InterpolationFractionNotIsometry(valarray<_Ty>& x, 
											valarray<_Ty>& y, _Ty t);

//����ʽ�Ⱦ��ֵ
template <class _Ty>
_Ty InterpolationFractionIsometry(_Ty x0, _Ty h, valarray<_Ty>& y, _Ty t);

//�������ز��Ⱦ��ֵ
template <class _Ty>
_Ty InterpolationHermiteNotIsometry(valarray<_Ty>& x, 
							valarray<_Ty>& y, valarray<_Ty>& dy, _Ty t);

//�������صȾ��ֵ
template <class _Ty>
_Ty InterpolationHermiteIsometry(_Ty x0, _Ty h, 
						valarray<_Ty>& y, valarray<_Ty>& dy, _Ty t);

//���ؽ𲻵Ⱦ��𲽲�ֵ
template <class _Ty>
_Ty InterpolationAitkenNotIsometry(valarray<_Ty>& x, 
									valarray<_Ty>& y, _Ty t, _Ty eps);

//���ؽ�Ⱦ��𲽲�ֵ
template <class _Ty>
_Ty InterpolationAitkenIsometry(_Ty x, _Ty h,  
									valarray<_Ty>& y, _Ty t, _Ty eps);

//�⻬���Ⱦ��ֵ
template <class _Ty>
void InterpolationSmoothNotIsometry(valarray<_Ty>& x, 
					valarray<_Ty>& y, int k, _Ty t, valarray<_Ty>& s);


//�⻬�Ⱦ��ֵ
template <class _Ty>
void InterpolationSmoothIsometry(_Ty x, _Ty h, 
					valarray<_Ty>& y, int k, _Ty t, valarray<_Ty>& s);

//��һ�ֱ߽���������������������ֵ��΢�������
template <class _Ty>
_Ty Interpolation3Spooling1stBoundary(valarray<_Ty>& x, 
		valarray<_Ty>& y, valarray<_Ty>& dy, valarray<_Ty>& ddy, 
			valarray<_Ty>& t, valarray<_Ty>& z, valarray<_Ty>& dz, 
													valarray<_Ty>& ddz);

//�ڶ��ֱ߽���������������������ֵ��΢�������
template <class _Ty>
_Ty Interpolation3Spooling2ndBoundary(valarray<_Ty>& x, 
		valarray<_Ty>& y, valarray<_Ty>& dy, valarray<_Ty>& ddy, 
			valarray<_Ty>& t, valarray<_Ty>& z, valarray<_Ty>& dz, 
													valarray<_Ty>& ddz);

//�����ֱ߽���������������������ֵ��΢�������
template <class _Ty>
_Ty Interpolation3Spooling3thBoundary(valarray<_Ty>& x, 
		valarray<_Ty>& y, valarray<_Ty>& dy, valarray<_Ty>& ddy, 
			valarray<_Ty>& t, valarray<_Ty>& z, valarray<_Ty>& dz, 
													valarray<_Ty>& ddz);


//��Ԫ�����ֵ
template <class _Ty>
_Ty Interpolation2Variable3Points(valarray<_Ty>& x, valarray<_Ty>& y, 
											matrix<_Ty> z, _Ty u, _Ty v);

//��Ԫȫ�����ֵ
template <class _Ty>
_Ty Interpolation2VariableWholeInterval(valarray<_Ty>& x, 
						valarray<_Ty>& y, matrix<_Ty> z, _Ty u, _Ty v);

#include "Interpolation.inl"		//�༰��غ����Ķ���ͷ�ļ�

#endif		//_INTERPOLATION_H

