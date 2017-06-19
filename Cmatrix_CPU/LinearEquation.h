// LinearEquation.h		���Է���(��)��⺯��(����)����
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31

#ifndef _LINEAREQUATION_H
#define _LINEAREQUATION_H

#include <comm.h>		//����ͷ�ļ�
#include <valarray>		//ģ����valarray�ı�׼ͷ�ļ�
#include <Matrix.h>		//������ͷ�ļ�

//ȫѡ��Ԫ��˹��Ԫ��
template <class _Ty>		//���췵��0����������
int LE_TotalChoiceGauss(matrix<_Ty>& a, valarray<_Ty>& b);

//ȫѡ��Ԫ��˹-Լ����Ԫ��
template <class _Ty>		//���췵��0����������
int LE_TotalChoiceGaussJordan(matrix<_Ty>& a, matrix<_Ty>& b);

//���ԽǷ������׷�Ϸ�
template <class _Ty>
int LE_TridiagonalEquationGauss(valarray<_Ty>& b, valarray<_Ty>& d);

//һ����ͷ��������
template <class _Ty>
int LE_StrapEquationGauss(matrix<_Ty>& b, matrix<_Ty>& d, int l, int il);

//�ԳƷ�����ķֽⷨ
template <class _Ty>
int LE_SymmetryEquation(matrix<_Ty>& a, matrix<_Ty>& c);

//�Գ������������ƽ������
template <class _Ty>
int LE_SymmetryRegularEuationSquareRoot(matrix<_Ty>& a, matrix<_Ty>& d);

//����ϡ�跽����ȫѡ��Ԫ��˹-Լ��(Gauss-Jordan)��
template <class _Ty>
int LE_SparseEuationTotalChoiceGaussJordan(matrix<_Ty>& a, valarray<_Ty>& b);

//�в�����(Toeplitz)�����������ѷ(Levinson)��
template <class _Ty>
int LE_ToeplitzEuationLevinson(valarray<_Ty>& t, valarray<_Ty>& b, valarray<_Ty>& x);

//��˹-���¶�������
template <class _Ty>
int LE_GaussSeidelIteration(matrix<_Ty>& a, valarray<_Ty>& b, valarray<_Ty>& x, _Ty eps);

//�Գ�����������Ĺ����ݶȷ�
template <class _Ty>
int LE_SymmetryRegularEuationConjugateGradient(matrix<_Ty>& a, valarray<_Ty>& b, _Ty eps, valarray<_Ty>& x);

//������С���˵ĺ�˹�ɶ��±任��
template <class _Ty>
int LE_LinearLeastSquareHouseholder(matrix<_Ty>& a, valarray<_Ty>& b, matrix<_Ty>& q);

//������С���˵Ĺ����淨
template <class _Ty>
int LE_LinearLeastSquareGeneralizedInverse(matrix<_Ty>& a, valarray<_Ty>& b, matrix<_Ty>& x, matrix<_Ty>& aa, _Ty eps, matrix<_Ty>& u, matrix<_Ty>& v);

//��̬���������
template <class _Ty>
int LE_IllConditionedEquation(matrix<_Ty>& a, valarray<_Ty>& b, _Ty eps, matrix<_Ty>& x);

#include <LinearEquation.inl>

#endif //_LINEAREQUATION_H

