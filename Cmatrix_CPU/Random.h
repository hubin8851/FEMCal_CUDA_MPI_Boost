//random.h			�궨�弰���������������ԭ��ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _RANDOM_H
#define _RANDOM_H

#define RandCoef2053 2053	
#define RandCoef13849 13849
#define MODUL65536 65536

#include <comm.h>		//��������ͷ�ļ�

//����һ��[0,1]�����ھ��ȷֲ�α�����
inline double rand_01_One(double& seed);

//�������[0,1]�����ھ��ȷֲ�α�����
inline void 
rand_01_Series(double& seed, valarray<double>& dp, const size_t stCount);

//��������[a,b]������һ�����ȷֲ�α�������
inline size_t 
rand_ab_One(size_t a, size_t b, size_t& seed);

//��������[a,b]�����ھ��ȷֲ�α�����������
inline void 
rand_ab_Series(size_t a, size_t b, size_t& seed, valarray<size_t>& sp, size_t stCount);

//����һ�������ֵ�뷽�����̬�ֲ������
inline double 
rand_NormalDistributing_One(double mu, double ro, double& seed);

//���������ֵ�뷽�����̬�ֲ����������
inline void 
rand_NormalDistributing_Series(double mu, double ro, double seed, valarray<double>& dp, size_t stCount);

#include <random.inl>

#endif //_RANDOM_H

