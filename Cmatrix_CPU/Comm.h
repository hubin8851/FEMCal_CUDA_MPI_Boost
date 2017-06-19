// Comm.h				����ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31

#ifndef _COMM_H		//�����α���
#define _COMM_H

 /******
     ����NDEBUG��Ϊ�������(debug)�ڼ�ʹ�õģ�����������ϣ������
 ���а�(release)��NDEBUG��ʧȥ���á�Ϊ���ܹ�ʹassert()�����ڷ���
 ����Ҳ����ʹ�ã��������������������NDEBUG����˼Ϊ��
     ����Ѿ�����NDEBUG�������Assert(x)����xΪ��ʱ��ִ��throw��
 ����ϵͳ�Ѷ���ķ�ʽ��������Assert(x)����Ϊ����assert(x)����
 ������ ����x���в��ԣ�xΪ��ʱ��������ֹ���򲢴�ӡ������Ϣ��
 *****/
 #ifdef NDEBUG			
  #define Assert(x)	if (!x)	throw;	
 #else			//������ϵͳ����ĺ���assert()����
  #include <cassert>
  #define Assert(x)	assert(x);
 #endif			//NDEBUG

 #include <complex>		//ģ����complex�ı�׼ͷ�ļ�
 #include <valarray>	//ģ����valarray�ı�׼ͷ�ļ�
 using namespace std;	//���ֿռ�

 const float  FLOATERROR = 1.0e-6F;
 const double DOUBLEERROR = 1.0e-15;
 const long double LONGDOUBLEERROR = 1.0e-30;

 const double GoldNo = 0.618033399;		//�ƽ�ָ��(1.0-0.381966)
 
 //ȡx����ֵ
 template <class T>		
 long double Abs(const T& x);				

 //ȡx���ţ�+-��0
 template <class T>		
 T Sgn(const T& x);
 
 //�Ƚ���float���������
 bool FloatEqual(float lhs, float rhs);		
 //�Ƚ���float�����������
 bool FloatNotEqual(float lhs, float rhs);		
 //�Ƚ���double���������
 bool FloatEqual(double lhs, double rhs);		
 //�Ƚ���double�����������
 bool FloatNotEqual(double lhs, double rhs);	
 //�Ƚ���long double���������
 bool FloatEqual(long double lhs, long double rhs);	
 //�Ƚ���long double�����������
 bool FloatNotEqual(long double lhs, long double rhs);	

 //��x��y����Сֵ������С��
 template <class T>
 T Min(const T& x, const T& y);

 //��x��y�����ֵ�����ش���
 template <class T>
 T Max(const T& x, const T& y);

 //��ӡ����(����)����Ԫ��ֵ
 template <class T>
 void ValarrayPrint(const valarray<T>& vOut);
 
 //��ӡĳ��ָ������(����)Ԫ��ֵ
 template <class T>
 void ValarrayPrint(const valarray<T>& vOut, size_t sPosition);
 
 #include "Comm.inl"	//ʵ��

#endif	// _COMM_H

