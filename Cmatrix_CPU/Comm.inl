// Comm.inl		公共函数（方法）定义
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31

#ifndef _COMM_INL
#define _COMM_INL

template <class T>
inline long double		//判断绝对值
Abs(const T& x)
{
	complex<long double> cld(x);
	long double ldAbs = abs(cld);
	return(ldAbs);
}

template <class T>
inline T			//取x符号，+-0
Sgn(const T& x)
{
	return x < T(0) ? T(-1) : (x > T(0) ? T(1) : T(0));
}


inline bool			//判断float浮点数相等
FloatEqual(float lhs, float rhs)
{
	if (Abs(lhs - rhs) < FLOATERROR)
		return true;
	else
		return false;
}

inline bool			//判断float浮点数不相等
FloatNotEqual(float lhs, float rhs)
{
	if (Abs(lhs - rhs) >= FLOATERROR)
		return true;
	else
		return false;
}

inline bool			//判断double浮点数相等
FloatEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) < DOUBLEERROR)
		return true;
	else
		return false;
}

inline bool			//判断double浮点数不相等
FloatNotEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) >= DOUBLEERROR)
		return true;
	else
		return false;
}

inline bool				//比较两long double浮点数相等
FloatEqual(long double lhs, long double rhs)
{
	if (Abs(lhs - rhs) < LONGDOUBLEERROR)
		return true;
	else
		return false;
}

inline bool				//比较两long double浮点数不相等
FloatNotEqual(long double lhs, long double rhs)
{
	if (Abs(lhs - rhs) >= LONGDOUBLEERROR)
		return true;
	else
		return false;
}
 
template <class T>
inline T
Min(const T& x, const T& y)			//比较x与y，返回小者
{
	if(x < y)
		return x;
	else
		return y;
}

template <class T>
T Max(const T& x, const T& y)		//求x与y的最大值，返回大者
{
	if(x > y)
		return x;
	else
		return y;
}

//打印数组(向量)所有元素值
template <class T>
void ValarrayPrint(const valarray<T>& vOut)
{
	size_t vsize = vOut.size();		//取数组元素的个数
	for(size_t st=0; st<vsize; st++)
	{
		cout.width(15);				//元素对齐，让每个元素占15列
		cout << vOut[st] << ' ';
	}
	cout << endl;
}

//打印某个指定数组(向量)元素值
template <class T>
void ValarrayPrint(const valarray<T>& vOut, size_t sPosition)
{
	cout.width(15);					//元素对齐，让每个元素占15列
	cout << vOut[sPosition] << endl;
}

#endif			// _COMM_INL