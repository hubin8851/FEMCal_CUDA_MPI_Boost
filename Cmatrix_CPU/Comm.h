// Comm.h				公共头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31

#ifndef _COMM_H		//避免多次编译
#define _COMM_H

 /******
     由于NDEBUG是为程序调试(debug)期间使用的，当调试期完毕，程序出
 发行版(release)后，NDEBUG将失去作用。为了能够使assert()函数在发行
 版中也可以使用，则定义了下面的条件编译NDEBUG，意思为：
     如果已经定义NDEBUG，则定义宏Assert(x)：当x为假时，执行throw，
 交由系统已定义的方式处理；否则将Assert(x)定义为函数assert(x)，该
 函数对 参数x进行测试：x为假时，函数终止程序并打印错误信息。
 *****/
 #ifdef NDEBUG			
  #define Assert(x)	if (!x)	throw;	
 #else			//否则用系统定义的函数assert()处理
  #include <cassert>
  #define Assert(x)	assert(x);
 #endif			//NDEBUG

 #include <complex>		//模板类complex的标准头文件
 #include <valarray>	//模板类valarray的标准头文件
 using namespace std;	//名字空间

 const float  FLOATERROR = 1.0e-6F;
 const double DOUBLEERROR = 1.0e-15;
 const long double LONGDOUBLEERROR = 1.0e-30;

 const double GoldNo = 0.618033399;		//黄金分割常数(1.0-0.381966)
 
 //取x绝对值
 template <class T>		
 long double Abs(const T& x);				

 //取x符号，+-或0
 template <class T>		
 T Sgn(const T& x);
 
 //比较两float浮点数相等
 bool FloatEqual(float lhs, float rhs);		
 //比较两float浮点数不相等
 bool FloatNotEqual(float lhs, float rhs);		
 //比较两double浮点数相等
 bool FloatEqual(double lhs, double rhs);		
 //比较两double浮点数不相等
 bool FloatNotEqual(double lhs, double rhs);	
 //比较两long double浮点数相等
 bool FloatEqual(long double lhs, long double rhs);	
 //比较两long double浮点数不相等
 bool FloatNotEqual(long double lhs, long double rhs);	

 //求x与y的最小值，返回小者
 template <class T>
 T Min(const T& x, const T& y);

 //求x与y的最大值，返回大者
 template <class T>
 T Max(const T& x, const T& y);

 //打印数组(向量)所有元素值
 template <class T>
 void ValarrayPrint(const valarray<T>& vOut);
 
 //打印某个指定数组(向量)元素值
 template <class T>
 void ValarrayPrint(const valarray<T>& vOut, size_t sPosition);
 
 #include "Comm.inl"	//实现

#endif	// _COMM_H

