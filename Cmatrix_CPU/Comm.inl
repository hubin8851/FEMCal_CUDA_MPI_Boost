// Comm.inl		��������������������
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31

#ifndef _COMM_INL
#define _COMM_INL

template <class T>
inline long double		//�жϾ���ֵ
Abs(const T& x)
{
	complex<long double> cld(x);
	long double ldAbs = abs(cld);
	return(ldAbs);
}

template <class T>
inline T			//ȡx���ţ�+-0
Sgn(const T& x)
{
	return x < T(0) ? T(-1) : (x > T(0) ? T(1) : T(0));
}


inline bool			//�ж�float���������
FloatEqual(float lhs, float rhs)
{
	if (Abs(lhs - rhs) < FLOATERROR)
		return true;
	else
		return false;
}

inline bool			//�ж�float�����������
FloatNotEqual(float lhs, float rhs)
{
	if (Abs(lhs - rhs) >= FLOATERROR)
		return true;
	else
		return false;
}

inline bool			//�ж�double���������
FloatEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) < DOUBLEERROR)
		return true;
	else
		return false;
}

inline bool			//�ж�double�����������
FloatNotEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) >= DOUBLEERROR)
		return true;
	else
		return false;
}

inline bool				//�Ƚ���long double���������
FloatEqual(long double lhs, long double rhs)
{
	if (Abs(lhs - rhs) < LONGDOUBLEERROR)
		return true;
	else
		return false;
}

inline bool				//�Ƚ���long double�����������
FloatNotEqual(long double lhs, long double rhs)
{
	if (Abs(lhs - rhs) >= LONGDOUBLEERROR)
		return true;
	else
		return false;
}
 
template <class T>
inline T
Min(const T& x, const T& y)			//�Ƚ�x��y������С��
{
	if(x < y)
		return x;
	else
		return y;
}

template <class T>
T Max(const T& x, const T& y)		//��x��y�����ֵ�����ش���
{
	if(x > y)
		return x;
	else
		return y;
}

//��ӡ����(����)����Ԫ��ֵ
template <class T>
void ValarrayPrint(const valarray<T>& vOut)
{
	size_t vsize = vOut.size();		//ȡ����Ԫ�صĸ���
	for(size_t st=0; st<vsize; st++)
	{
		cout.width(15);				//Ԫ�ض��룬��ÿ��Ԫ��ռ15��
		cout << vOut[st] << ' ';
	}
	cout << endl;
}

//��ӡĳ��ָ������(����)Ԫ��ֵ
template <class T>
void ValarrayPrint(const valarray<T>& vOut, size_t sPosition)
{
	cout.width(15);					//Ԫ�ض��룬��ÿ��Ԫ��ռ15��
	cout << vOut[sPosition] << endl;
}

#endif			// _COMM_INL