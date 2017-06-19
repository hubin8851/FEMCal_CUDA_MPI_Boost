#ifndef _CMATRIX_H
#define _CMATRIX_H

#include "stdafx.h"

namespace HBXDef
{


//矩阵模板类
template <class T=float>
class CMatrix
{
	typedef CMatrix<T>	_SameClass;
private:
	HBXDef::SortType_t	m_eMajorType;//是行为主序还是列为主序,4B
	bool		Isuse;//多线程中为了mutex可能会用到,1B
	size_t		m_row;//行数,4B
	size_t		m_col;//列数,4B
	std::valarray<T>	m_data;//数据向量指针,16B
public:
	//默认构造函数
	CMatrix()
	{
		m_eMajorType = HBXDef::ROWMAJOR;
		Isuse =	false;
		m_row = 0;
		m_col = 0;
	}
	//两参数的构造函数,
	//@_row:申请的行数
	//@_row:申请的列数
	CMatrix(size_t _row, size_t _col, HBXDef::SortType_t _Srtty =  HBXDef::ROWMAJOR ):m_row(_row),m_col(_col),m_eMajorType(_Srtty) 
	{ 
		m_data.resize(m_row * m_col, (T)0); 
	} 
	//三参数的构造函数
	//@rhs:所要填方数组的首地址指针,采用了valarray(const _Ty *p, size_t n)的构造形式
	//@_row:所要分配的行数
	//@_col：索要分配的列数
	CMatrix(const T *rhs, size_t _row, size_t _col):m_data(rhs, _row*_col), m_row(_row), m_col(_col)
	{
		m_eMajorType = HBXDef::ROWMAJOR;
	};
	//拷贝构造函数
	CMatrix(const _SameClass &rhs):m_data(rhs.m_data), m_row(rhs.GetRowNum()), m_col(rhs.GetColNum())
	{
		m_eMajorType = HBXDef::ROWMAJOR;
	}; 
	~CMatrix(){};
	//二元符号重载
	_SameClass& operator =(const _SameClass& rhs)//重载=
	{
		if(m_row != rhs.m_row || m_col != rhs.m_col) error("矩阵类两者行列不尽相同！");
		if (this!=&rhs)
		{
			m_data = rhs.m_data;
		}
		return *this;
	}
	//转移赋值，C11新特性
	_SameClass&	operator =(_SameClass&& _rhs)//移动语义赋值
	{
		m_row = _rhs.m_row;
		m_col = _rhs.m_col;
		m_data = _rhs.m_data;
		_rhs.m_data.free();
		return *this;
	}
	friend std::istream& operator>>(std::istream&, const CMatrix&);//重载输入运算符
	friend std::ostream& operator<<(std::ostream&, const CMatrix&);//重载输出运算符
	T	maxnorm() const;//获得绝对值最大值
	T	twonorm() const;//获得所有元素方差

	//从下三角矩阵变为完全形式
	void	Symmetrized()
	{
	#ifdef DEBUG
		Assert( m_row == m_col );
	#endif
		for ( size_t i = 1; i < m_row; i++ ) {
			for ( size_t j = 0; j < i; j++ ) {
				this->m_data(i, j) = this->m_data(j, i);
			}
		}
	}

	size_t	ResetSize(const size_t &nRow, const size_t &nCol, HBXDef::SortType_t nFlag = HBXDef::ROWMAJOR)
	{
		m_eMajorType = nFlag;
		m_row = nRow;
		m_col = nCol;
		m_data.resize(m_row * m_col, (T)0);
		return true;
	}
	size_t	GetSize(int &_nRow, int &_nCol)
	{
		_nRow = (int)this->m_row;
		_nCol = (int)this->m_col;
		return true;
	}
	//返回行数
	size_t	GetRowNum() const
	{
		return m_row;
	}
	//返回列数
	size_t	GetColNum() const
	{
		return m_col;
	}
	//设置所有m行n列的元素为常数
	int		SetData(const size_t &_iRow, const size_t &_iCol,const T _value)
	{
		Assert(_iRow <= this->m_row && _iCol <= this->m_col);
		this->m_data[_iRow * this->GetColNum() + _iCol] = _value;
		return true;
	}
	T		GetData(size_t _iRow, size_t _iCol)
	{
		return CMatrix(_iRow, _iCol);
	}

	//重载函数调用运算符(),返回二维矩阵中元素(_row, _col)在一维数组中的位置,作为左值被写时使用
	T& operator()(size_t _row, size_t _col)
	{
		if (_row < this->GetRowNum() && _col < this->GetColNum())
		{
			return m_data[_row * GetColNum() + _col];
		}
		else
		{
			std::cerr<<"索引越界"<<std::endl;
			//这里要注意，没有返回值，如果索引越界可能会报错。
		}
	}

	//重载函数调用运算符(),返回二维矩阵中元素(_row, _col)在一维数组中的位置,作为右值被读时使用
	const T operator()(size_t _row, size_t _col) const
	{
		if (_row > this->GetRowNum() || _col > this->GetColNum())
		{
			//这里要注意，没有返回值，如果索引越界可能会报错。
			std::cout<<"索引越界"<<std::endl;
		}
//		Assert(_row < this->GetRowNum());//断定_row不超实际矩阵行值,Assert已经被重定义
//		Assert(_col < this->GetColNum());//断定_col不超实际矩阵行值,Assert已经被重定义
		return m_data[_row * GetColNum() + _col];
	}
	//赋值操作符,矩阵加等一个同行数同列数的矩阵
	_SameClass& operator +=(_SameClass& rhs)
	{
		Assert(this->m_row == rhs.GetRowNum());//断定_row不超实际矩阵行值,Assert已经被定义
		Assert(this->m_col == rhs.GetColNum());//断定_col不超实际矩阵行值,Assert已经被定义
		//利用类valarray自身的重载+=
		m_data += rhs.m_data;
		return *this;
	}
	//赋值操作符,矩阵加等一个常数值
	_SameClass& operator +=(const T& rhsVal)
	{
		//利用类valarray自身的重载+=
		m_data += rhsVal;
		return *this;
	}
	//赋值操作符,矩阵减等一个同行数同列数的矩阵
	_SameClass& operator -=(_SameClass& rhs)
	{
		Assert(this->m_row == rhs.GetRowNum());//断定_row不超实际矩阵行值,Assert已经被定义
		Assert(this->m_col == rhs.GetColNum());//断定_col不超实际矩阵行值,Assert已经被定义
		//利用类valarray自身的重载-=
		m_data -= rhs.m_data;
		return *this;
	}
	//赋值操作符,矩阵减等一个常数值
	_SameClass& operator -=(const T& rhsVal)
	{
		//利用类valarray自身的重载-=
		m_data -= rhsVal;
		return *this;
	}
	//赋值操作符,矩阵乘等一个矩阵（向量）
	_SameClass& operator *=(_SameClass& rhs)
	{
		MatrixMultiply(*this, *this, rhs);	//调用矩阵乘法函数
		return *this;
	}
	//赋值操作符,矩阵乘等一个常数值
	_SameClass& operator *=(const T& rhsVal)
	{
		//利用类valarray自身的重载*=
		m_data *= rhsVal;
		return *this;
	}
	//赋值操作符,矩阵除等一个常数值
	_SameClass& operator /=(const T& rhsVal)
	{
		if ( 0 == rhsVal)
		{
			std::cout<<"右值不能为0"<<std::endl;
			return *this;
		}
		//利用类valarray自身的重载/=
		m_data /= rhsVal;
		return *this;
	}
	//一元+操作符，所有值取绝对值
	_SameClass	operator +() const
	{
		_SameClass _tempMat(*this);
		_tempMat.m_data = abs(this->m_data);
		return _tempMat;
	}
	//一元操作符，做CMatrix内部数据取负值
	_SameClass	operator -() const
	{
		_SameClass _tempMat(*this);
		_tempMat.m_data = -_tempMat.m_data;
		return _tempMat;
	}
	//二元操作符
	//矩阵加矩阵
	_SameClass operator +(const _SameClass &rhs)
	{
		Assert(this->GetRowNum() == rhs.GetRowNum());
		Assert(this->GetColNum() == rhs.GetColNum());
		_SameClass _tempMat(*this);
		_tempMat.m_data += rhs.m_data;
		return _tempMat;
	}
	//矩阵加常数
	_SameClass operator +(T &rhsVal)
	{
		_SameClass _tempMat(*this);
		_tempMat += rhsVal;
		return _tempMat;
	}
	//矩阵减矩阵
	_SameClass operator -(const _SameClass &rhs)
	{
		Assert(this->GetRowNum() == rhs.GetRowNum());
		Assert(this->GetColNum() == rhs.GetColNum());
		_SameClass _tempMat(*this);
		_tempMat -= rhs;
		return _tempMat;
	}
	//矩阵减常数
	_SameClass operator -(const T &rhsVal)
	{
		_SameClass _tempMat(*this);
		_tempMat -= rhsVal;
		return _tempMat;
	}
	//矩阵乘常数
	_SameClass operator *(const T &rhsVal)
	{
		_SameClass _tempMat(*this);
		_tempMat *= rhsVal;
		return _tempMat;
	}
	//矩阵乘矩阵
	_SameClass	operator *(const _SameClass &rhs)
	{
		if ( this->GetColNum() != rhs.GetRowNum() )
		{
			std::cerr<<"相乘矩阵维度不一致"<<std::endl;
		}
		_SameClass	_tempmat( this->GetRowNum(), rhs.GetColNum() );
		MatrixMultiply(_tempmat, *this, rhs);
		return _tempmat;
	}
	//矩阵除常数
	_SameClass operator /(const T &rhsVal)
	{
		CMatrix _tempMat(*this);
		_tempMat /= rhsVal;
		return _tempMat;
	}
	//常数加矩阵
	friend	_SameClass	operator + (const T &lhsVal, const _SameClass &rhs)
	{
		_SameClass _tempMat(rhs);
		_tempMat.m_data += lhsVal;
		return _tempMat;
	}
	//常数减矩阵
	friend	_SameClass	operator - (const T &lhsVal, const _SameClass &rhs)
	{
		_SameClass _tempMat(rhs);
		_tempMat.m_data = lhsVal - _tempMat.m_data;
		return _tempMat;
	}
	//常数乘矩阵
	friend	_SameClass	operator * (const T &lhsVal, const _SameClass &rhs)
	{
		_SameClass _tempMat(rhs);
		_tempMat.m_data *= lhsVal;
		return _tempMat;
	}
	//矩阵比较
	//比较两矩阵是否不相等，若相等返回true，否则返回false
	bool	operator != (const _SameClass &rhs)
	{
		if (this->GetRowNum() != rhs.GetRowNum()) return true;
		if (this->GetColNum() != rhs.GetColNum()) return true;
		for (size_t i=0; i<rhs.m_data.size(); i++)
		{
			if (this->m_data[i] != rhs.m_data[i])
			{
				return true;
			}
		}
		return false;
	}
	//比较两矩阵是否相等，若相等返回true，否则返回false
	bool	operator == (const _SameClass &rhs)
	{
		if (this->GetRowNum() != rhs.GetRowNum()) return false;
		if (this->GetColNum() != rhs.GetColNum()) return false;
		for (size_t i=0; i<rhs.m_data.size(); i++)
		{
			if (this->m_data[i] != rhs.m_data[i])
			{
				return false;
			}
		}
		return true;
	}
	//重载转置
	_SameClass operator!() 
	{ 
		_SameClass temp( m_col, m_row); 
		for(int i=0; i<m_col; i++) 
			for(int j=0; j<m_row; j++) 
				temp(i,j)=(*this)(j,i); 
		return temp; 
	} 
	//不要轻易使用，矩阵大的时候用printf输出很慢，用于小矩阵测试时验证
	int Print();
};//class定义结束

typedef CMatrix<float> CMatrixf, CMatrix_f;
typedef CMatrix<double> CMatrixd, CMatrix_d;
typedef CMatrix<long double> CMatrixld, CMatrix_ld;

inline void error(char* v)
{
	std::cout<<v<<".program exited\n";
	getchar();
	exit(1);
}

//重载输入运算符 
template<class T> 
std::istream &operator >> (std::istream &in,CMatrix<T> &x) 
{ 
	for(int i=0; i<x.m_row(); i++) 
	{ 
		for(int j=0; j<x.m_col(); j++) 
			in>>x(i,j); 
	} 
	return in; 
} 

//重载输出<< 
template <class T> 
std::ostream &operator<<(std::ostream &out,CMatrix<T> &x) 
{ 
	for(int i=0; i<x.m_row(); i++) 
	{ 
		for(int j=0; j<x.m_col(); j++) 
			out<<x(i,j)<<' '; 
		out<<endl; 
	} 
	out<<endl; 
	return out; 
} 

//绝对值最大值索引
template <class T>
T CMatrix<T>::maxnorm() const
{
	T norm=fabs(m_data[0]);
	for (int i=1;i<lenth;i++)
		norm=max(norm, fabs(m_data[i]));
	return norm;
}

//求元素的方差,二范数
template <class T>
T CMatrix<T>::twonorm() const
{
	T norm = m_data[0] * m_data[0];
	int temp = m_row * m_col;
	for (int i=1;i<temp;i++) norm += m_data[i] * m_data[i];
	return sqrt(norm);
}

template <class T>
int CMatrix<T>::Print()
{
	for (int _row= 0;_row<m_row ;++_row)
	{
		for (int _col= 0;_col<m_col ;++_col)
		{
			printf("%15.5f	",m_data[_row*m_col+_col]);
		}
		std::cout<<endl;
	}
	return 0;
};


//矩阵乘法函数
template <class _TOut, class _Tlhs, class _Trhs>	//最后结果在mOut中
CMatrix<_TOut>& MatrixMultiply(CMatrix<_TOut>& matOut, const CMatrix<_Tlhs>& matlhs, const CMatrix<_Trhs>& matrhs);

//全选主元法求矩阵行列式
template<class _Ty>
long double MatrixDetrminant( const CMatrix<_Ty>& _rhs );

//用全选主元高斯-约当(Gauss-Jordan)消去法计算实(复)矩阵的逆矩阵
//矩阵类型必须是浮点型
template <class _Ty>
int MatrixInversionGS( CMatrix<_Ty>& _rhs );

}//end namespace hbxdef

#include "MatrixCal.inl"



#endif