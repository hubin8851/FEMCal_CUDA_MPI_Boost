// Matrix.h		矩阵模板类头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _MATRIX_H		//避免多次编译
#define _MATRIX_H

#include <valarray>		//模板类valarray的标准头文件
#include "Comm.h"		//公共头文件
#include <math.h>		//数学头文件

template <class _Ty>
class matrix			//矩阵类matrix
{
	typedef matrix<_Ty> _Myt;

 private:
	std::valarray<_Ty> m_Datas;	//定义一维数组对象m_Datas

	size_t m_stRow;		//矩阵行数变量
	size_t m_stCol;		//矩阵列数变量

 public:
	typedef _Ty value_type;

	//构造函数一(参数1, 2分别为矩阵的行与列数)
/******
矩阵类matrix的构造函数一
	构造函数中出现的m_Datas为valarray类的对象，申请stRow * stCol个单元，
单元内没赋值。对数组对象m_Datas使用了valarray类的构造函数：
				explicit valarray(size_t n)
    对私有变量m_stRow和m_stCol分别赋初值stRow, stCol。
******/
	matrix(size_t stRow, size_t stCol)
		: m_Datas(stRow * stCol),
			m_stRow(stRow), m_stCol(stCol)
	{
		m_Datas.resize(GetRowNum() * GetColNum(), _Ty(0));
		/*
		构造函数中出现的m_Datas为valarray类的对象，若在函数体内调用的
		resize(size_t n, const T& c = T())函数有两个参数(参阅valarray中的
		定义)，第一个参数申请有矩阵行*列那么多元素个数的存储空间，第二个
		参数对这些申请的空间赋予具有模式T的初值0。如果不调用该函数，则缺省
		情况下m_Datas长度为0；另外，调用该函数但不使用第二个参数，则不会对
		m_Datas的任何一个元素赋初值。
		*/
	}

	//构造函数二(参数1为指向矩阵的指针，参数2, 3为矩阵行与列数)
/******
 矩阵类matrix的构造函数二
     对私有变量m_stRow和m_stCol分别赋初值stRow, stCol。
     对数组对象m_Datas使用了valarray类的构造函数：
				valarray(const _Ty *p, size_t n)
 m_Datas初始化的第一参数为矩阵rhs指针，第二个参数为rhs的元素总个数，
 即rhs行数*列数
******/
	matrix(const _Ty* rhs, size_t stRow, size_t stCol)
		: m_Datas(rhs, stRow * stCol),
			m_stRow(stRow), m_stCol(stCol)
	{
	}

	//构造函数三(拷贝构造函数，参数为对矩阵matrix的引用)
/******
 矩阵类matrix的构造函数三
     用引用矩阵rhs的数组对象m_Datas初始化matrix所定义对象的m_Datas，
 用引用矩阵rhs的行数rhs.GetRowNum()和列数rhs.GetColNum()分别初始化私
 有变量m_stRow和m_stCol
******/
	matrix(const _Myt& rhs)
		: m_Datas(rhs.m_Datas),
			m_stRow(rhs.GetRowNum()), m_stCol(rhs.GetColNum())
	{
	}

	size_t GetRowNum() const	//返回矩阵行数的函数
	{
		return m_stRow;
	}

	size_t GetColNum() const	//返回矩阵列数的函数
	{
		return m_stCol;
	}

	/*****
	    重载运算符()，获取二维矩阵某元素在一维数组m_Datas的位置
	(matrix类定义的矩阵实际上已经被转化为一维数组m_Datas)
	    当数组元素(i, j)在运算符左边(被写)时，使用_Ty& operator ()
	当数组元素(i, j)在运算符右边(被读)时，使用_Ty operator ()    
	*****/
	_Ty& operator () (size_t stRow, size_t stCol)
	{
		Assert(stRow < GetRowNum());	//断定stRow不超实际矩阵行值
		Assert(stCol < GetColNum());	//断定stCol不超实际矩阵列值

		return m_Datas[stRow * GetColNum() + stCol];
	}

/******
 重载()运算符，返回二维矩阵中元素(stRow, stCol)在一维数组中的位置
 stRow与stCol分别为矩阵某元素行与列的位置数
 当矩阵元素在运算符右边(被读)时使用
******/
	const _Ty operator () (size_t stRow, size_t stCol) const
	{
		Assert(stRow < GetRowNum());	//断定stRow不超实际矩阵行值
		Assert(stCol < GetColNum());	//断定stCol不超实际矩阵列值

		return m_Datas[stRow * GetColNum() + stCol];
	}

	// 赋值操作符
	//矩阵与矩阵的自反*, +, -运算符
	_Myt& operator += (const _Myt& rhs)		//矩阵与矩阵的自反+
	{
		Assert(GetRowNum() == rhs.GetRowNum());	//断定两矩阵的行数相等
		Assert(GetColNum() == rhs.GetColNum());	//断定两矩阵的列数相等
		//利用类valarray定义，使其对象m_Datas定义进行两矩阵相加
		m_Datas += rhs.m_Datas;
		//结果放在左边矩阵中，返回指向左边矩阵的指针
		return *this;
	}

	_Myt& operator -= (const _Myt& rhs)		//矩阵与矩阵的自反-
	{
		Assert(GetRowNum() == rhs.GetRowNum());
		Assert(GetColNum() == rhs.GetColNum());
		//利用类valarray定义，使其对象m_Datas定义进行两矩阵相减
		m_Datas -= rhs.m_Datas;
		//结果放在左边矩阵中，返回指向左边矩阵的指针
		return *this;
	}

	_Myt& operator *= (const _Myt& rhs)		//矩阵与矩阵的自反*
	{
		MatrixMultiply(*this, *this, rhs);	//调用矩阵乘法函数

		return *this;						//最后结果放在左边矩阵中
	}

	//矩阵自反加、减、乘以、除以数
	_Myt& operator += (const _Ty& rhs)		//矩阵自加数
	{
		m_Datas += rhs;			//利用数组对象对每个元素加数
	
		return *this;			//结果放在原矩阵(数组m_Datas)中
	}

	_Myt& operator -= (const _Ty& rhs)		//矩阵自减数
	{
		m_Datas -= rhs;

		return *this;
	}

	_Myt& operator *= (const _Ty& rhs)		//矩阵自乘数
	{
		m_Datas *= rhs;
	
		return *this;
	}

	_Myt& operator /= (const _Ty& rhs)		//矩阵自除以数
	{
		m_Datas /= rhs;

		return *this;
	}

	// 一元操作符  对矩阵(每个元素)赋予+或-号
	_Myt operator + () const	//赋+号
	{
		return *this;			//不用作任何处理，维持原状
	}

	_Myt operator - () const	//赋-号
	{
		_Myt mat(*this);
		mat.m_Datas = -mat.m_Datas;		//每个元素赋-号

		return mat;
	}

	// 二元操作符
	//矩阵加数	mat = lhs + rhs
	friend _Myt operator + (const _Myt& lhs, const _Ty& rhs)
	{
		_Myt mat(lhs);			//新生成一新矩阵对象mat
		mat.m_Datas += rhs;		//对新矩阵对象每个元素加数

		return mat;				//返回新矩阵对象
	}

	//矩阵减数	mat = lhs - rhs
	friend _Myt operator - (const _Myt& lhs, const _Ty& rhs)
	{
		_Myt mat(lhs);
		mat.m_Datas -= rhs;

		return mat;
	}

	//矩阵乘以数	mat = lhs * rhs
	friend _Myt operator * (const _Myt& lhs, const _Ty& rhs)
	{
		_Myt mat(lhs);			//新生成一新矩阵对象mat
		mat.m_Datas *= rhs;		//对新矩阵对象每个元素乘以数

		return mat;				//返回新矩阵对象
	}

	//矩阵除以数	mat = lhs / rhs
	friend _Myt operator / (const _Myt& lhs, const _Ty& rhs)
	{
		_Myt mat(lhs);
		mat.m_Datas /= rhs;

		return mat;
	}

	//数加矩阵	mat = lhs + rhs
	friend _Myt operator + (const _Ty& lhs, const _Myt& rhs)
	{
		_Myt mat(rhs);			//新生成一新矩阵对象mat
		mat.m_Datas += lhs;		//数加上新矩阵对象的每个元素

		return mat;
	}

	//数减矩阵	mat = lhs - rhs
	friend _Myt operator - (const _Ty& lhs, const _Myt& rhs)
	{
		_Myt mat(rhs);			//新生成一新矩阵对象mat
		mat.m_Datas -= lhs;		//数减新矩阵对象的每个元素

		return mat;
	}

	//数乘以矩阵	mat = lhs * rhs
	friend _Myt operator * (const _Ty& lhs, const _Myt& rhs)
	{
		_Myt mat(rhs);			//新生成一新矩阵对象mat
		mat.m_Datas *= lhs;		//对新矩阵对象每个元素乘以数

		return mat;
	}

	//矩阵加矩阵	mat = lhs + rhs
	friend _Myt operator + (const _Myt& lhs, const _Myt& rhs)
	{
		_Myt mat(lhs);					//新生成一新矩阵对象mat, 用左边阵初始化
		mat.m_Datas += rhs.m_Datas;		//加上右边矩阵对象每个相应元素

		return mat;						//返回新矩阵对象
	}

	//矩阵减矩阵	mat = lhs - rhs
	friend _Myt operator - (const _Myt& lhs, const _Myt& rhs)
	{
		_Myt mat(lhs);					//新生成一新矩阵对象mat, 用左边阵初始化
		mat.m_Datas -= rhs.m_Datas;		//减去右边矩阵对象每个相应元素

		return mat;						//返回新矩阵对象
	}

	//矩阵乘以矩阵	mTmp = lhs * rhs
	friend _Myt operator * (const _Myt& lhs, const _Myt& rhs)
	{	//生成一个矩阵新对象mTmp		
		_Myt mTmp(lhs.GetRowNum(), rhs.GetColNum());	//没初始化

		MatrixMultiply(mTmp, lhs, rhs);			//调用矩阵乘法函数

		return mTmp;		//用新矩阵对象返回
	}

	//矩阵比较
	//比较两矩阵是否不相等，若相等返回true，否则返回false
	friend bool operator != (const _Myt& lhs, const _Myt& rhs)
	{
		if (lhs.GetRowNum() != rhs.GetRowNum())
			return true;
		if (lhs.GetColNum() != rhs.GetColNum())
			return true;

		for (size_t i = 0; i < lhs.m_Datas.size(); i ++)
		{
			if (lhs.m_Datas[i] != rhs.m_Datas[i])
				return true;
		}

		return false;
	}

	//比较两矩阵是否相等，若相等返回true，否则返回false
	friend bool operator == (const _Myt& lhs, const _Myt& rhs)
	{
		if (lhs.GetRowNum() != rhs.GetRowNum())
			return false;
		if (lhs.GetColNum() != rhs.GetColNum())
			return false;

		for (size_t i = 0; i < lhs.m_Datas.size(); i ++)
		{
			if (lhs.m_Datas[i] != rhs.m_Datas[i])
				return false;
		}

		return true;
	}
};

typedef matrix<float> matrixf;
typedef matrix<double> matrixd;
typedef matrix<long double> matrixld;

//矩阵乘法函数
template <class _Tyout, class _Tylhs, class _Tyrhs>	//最后结果在mOut中
matrix<_Tyout>& MatrixMultiply(matrix<_Tyout>& mOut, const matrix<_Tylhs>& lhs, const matrix<_Tyrhs>& rhs);

//输出矩阵		一行一行输出全矩阵
template <class _Ty>
void MatrixLinePrint(const matrix<_Ty>& mOut);

//输出矩阵		输出指定的某行
template <class _Ty>
void MatrixLinePrint(const matrix<_Ty>& mOut, size_t LineNo);

//矩阵转置		原阵在mIn，转置后的矩阵在mOut
template <class _Ty>
void MatrixTranspose(matrix<_Ty>& mIn, matrix<_Ty>& mOut);

//判断矩阵对称
template <class _Ty>		//对称返回true，否则返回false
bool MatrixSymmetry(const matrix<_Ty>& rhs);

//判断矩阵对称正定
template <class _Ty>
int MatrixSymmetryRegular(const matrix<_Ty>& rhs, int sym);

//全选主元法求矩阵行列式
template <class _Ty>		//返回值为行列式值
long double MatrixDeterminant(const matrix<_Ty>& rhs);

//全选主元高斯(Gauss)消去法求一般矩阵的秩
template <class _Ty>		//返回值为秩数
size_t MatrixRank(const matrix<_Ty>& rhs);

//用全选主元高斯-约当(Gauss-Jordan)消去法计算实(复)矩阵的逆矩阵
//矩阵类型必须是浮点型
template <class _Ty>
int MatrixInversionGS(matrix<_Ty>& rhs);

//用“变量循环重新编号法”法求对称正定矩阵逆
//矩阵类型必须是浮点型
template <class _Ty>
int MatrixSymmetryRegularInversion(matrix<_Ty>& rhs);

//特兰持(Trench)法求托伯利兹(Toeplitz)矩阵逆
//矩阵类型必须是浮点型
template <class _Ty>
int MatrixToeplitzInversionTrench(const valarray<_Ty>& t, 
								  const valarray<_Ty>& tt, matrix<_Ty>& rhs);

//实矩阵LU分解
//矩阵类型必须是浮点型
template <class _Ty>
int MatrixLU(const matrix<_Ty>& rhs, matrix<_Ty>& lhs, matrix<_Ty>& uhs);

//用豪斯荷尔德(Householder)变换对一般m*n阶的实矩阵进行QR分解
//矩阵类型必须是浮点型
template <class _Ty>
int MatrixQR(matrix<_Ty>& rhs, matrix<_Ty>& rhq);

//对称正定阵的乔里斯基(Cholesky)分解及求其行列式值
//矩阵类型必须是浮点型
template <class _Ty>
long double MatrixSymmetryRegularCholesky(const matrix<_Ty>& rhs);

//一般实矩阵的奇异值分解
//矩阵类型必须是浮点型
template <class _Ty>
int MatrixSingularValue(matrix<_Ty>& a, matrix<_Ty>& u, 
												matrix<_Ty>& v, _Ty eps);

//广义逆的奇异值分解
//矩阵类型必须是浮点型
template <class _Ty>
int GeneralizedInversionSingularValue(matrix<_Ty>& a, matrix<_Ty>& aa, 
								_Ty eps, matrix<_Ty>& u, matrix<_Ty>& v);

#include "Matrix.inl"		//类及相关函数的定义头文件

#endif		// _MATRIX_H

