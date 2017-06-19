//polynomials.h		����ʽ������ʽͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _POLYNOMIALS_H
#define _POLYNOMIALS_H

#include <comm.h>		 //��������ͷ�ļ�
#include <math.h>		 //��ѧͷ�ļ�
#include <matrix.h>		 //������ͷ�ļ�

//��һάʵ(��)����ʽֵ
template<class T, class U>
inline U 
PolyValueOneDim(valarray<T>& dCoff, size_t stNo, U dX);

//��һά����ʽ��ֵ
template<class T, class V, class U>
inline void 
PolyValueOneDimGroup(valarray<T>& dCoff, valarray<V>& dX, valarray<U>& dValue);

//���άʵ(��)����ʽֵ
template<class T, class U>
inline U	//��������
PolyValueTwoDim(matrix<T>& dCoff, U dX, U dY);

//��һά����ʽ���
template<class T>
inline void
PolyMultip(valarray<T>& dCoffP, valarray<T>& dCoffQ, valarray<T>& dCoffS);

//��һά����ʽ����
template<class T>
inline int
PolyDiv(valarray<T>& dCoffP, valarray<T>& dCoffQ, valarray<T>& dCoffS, valarray<T>& dCoffR);

// ��������ʽ����ֵ
template<class T>
inline T
FractionValue(valarray<T>& dXpara, valarray<T>& dCoff, T dX);

#include "polynomials.inl"	//����ʽ������ʽ��غ�������ͷ�ļ�

#endif //_POLYNOMIALS_H

