#ifndef _CMATRIX_H
#define _CMATRIX_H

#include "stdafx.h"

namespace HBXDef
{


//����ģ����
template <class T=float>
class CMatrix
{
	typedef CMatrix<T>	_SameClass;
private:
	HBXDef::SortType_t	m_eMajorType;//����Ϊ��������Ϊ����,4B
	bool		Isuse;//���߳���Ϊ��mutex���ܻ��õ�,1B
	size_t		m_row;//����,4B
	size_t		m_col;//����,4B
	std::valarray<T>	m_data;//��������ָ��,16B
public:
	//Ĭ�Ϲ��캯��
	CMatrix()
	{
		m_eMajorType = HBXDef::ROWMAJOR;
		Isuse =	false;
		m_row = 0;
		m_col = 0;
	}
	//�������Ĺ��캯��,
	//@_row:���������
	//@_row:���������
	CMatrix(size_t _row, size_t _col, HBXDef::SortType_t _Srtty =  HBXDef::ROWMAJOR ):m_row(_row),m_col(_col),m_eMajorType(_Srtty) 
	{ 
		m_data.resize(m_row * m_col, (T)0); 
	} 
	//�������Ĺ��캯��
	//@rhs:��Ҫ�������׵�ַָ��,������valarray(const _Ty *p, size_t n)�Ĺ�����ʽ
	//@_row:��Ҫ���������
	//@_col����Ҫ���������
	CMatrix(const T *rhs, size_t _row, size_t _col):m_data(rhs, _row*_col), m_row(_row), m_col(_col)
	{
		m_eMajorType = HBXDef::ROWMAJOR;
	};
	//�������캯��
	CMatrix(const _SameClass &rhs):m_data(rhs.m_data), m_row(rhs.GetRowNum()), m_col(rhs.GetColNum())
	{
		m_eMajorType = HBXDef::ROWMAJOR;
	}; 
	~CMatrix(){};
	//��Ԫ��������
	_SameClass& operator =(const _SameClass& rhs)//����=
	{
		if(m_row != rhs.m_row || m_col != rhs.m_col) error("�������������в�����ͬ��");
		if (this!=&rhs)
		{
			m_data = rhs.m_data;
		}
		return *this;
	}
	//ת�Ƹ�ֵ��C11������
	_SameClass&	operator =(_SameClass&& _rhs)//�ƶ����帳ֵ
	{
		m_row = _rhs.m_row;
		m_col = _rhs.m_col;
		m_data = _rhs.m_data;
		_rhs.m_data.free();
		return *this;
	}
	friend std::istream& operator>>(std::istream&, const CMatrix&);//�������������
	friend std::ostream& operator<<(std::ostream&, const CMatrix&);//������������
	T	maxnorm() const;//��þ���ֵ���ֵ
	T	twonorm() const;//�������Ԫ�ط���

	//�������Ǿ����Ϊ��ȫ��ʽ
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
	//��������
	size_t	GetRowNum() const
	{
		return m_row;
	}
	//��������
	size_t	GetColNum() const
	{
		return m_col;
	}
	//��������m��n�е�Ԫ��Ϊ����
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

	//���غ������������(),���ض�ά������Ԫ��(_row, _col)��һά�����е�λ��,��Ϊ��ֵ��дʱʹ��
	T& operator()(size_t _row, size_t _col)
	{
		if (_row < this->GetRowNum() && _col < this->GetColNum())
		{
			return m_data[_row * GetColNum() + _col];
		}
		else
		{
			std::cerr<<"����Խ��"<<std::endl;
			//����Ҫע�⣬û�з���ֵ���������Խ����ܻᱨ��
		}
	}

	//���غ������������(),���ض�ά������Ԫ��(_row, _col)��һά�����е�λ��,��Ϊ��ֵ����ʱʹ��
	const T operator()(size_t _row, size_t _col) const
	{
		if (_row > this->GetRowNum() || _col > this->GetColNum())
		{
			//����Ҫע�⣬û�з���ֵ���������Խ����ܻᱨ��
			std::cout<<"����Խ��"<<std::endl;
		}
//		Assert(_row < this->GetRowNum());//�϶�_row����ʵ�ʾ�����ֵ,Assert�Ѿ����ض���
//		Assert(_col < this->GetColNum());//�϶�_col����ʵ�ʾ�����ֵ,Assert�Ѿ����ض���
		return m_data[_row * GetColNum() + _col];
	}
	//��ֵ������,����ӵ�һ��ͬ����ͬ�����ľ���
	_SameClass& operator +=(_SameClass& rhs)
	{
		Assert(this->m_row == rhs.GetRowNum());//�϶�_row����ʵ�ʾ�����ֵ,Assert�Ѿ�������
		Assert(this->m_col == rhs.GetColNum());//�϶�_col����ʵ�ʾ�����ֵ,Assert�Ѿ�������
		//������valarray���������+=
		m_data += rhs.m_data;
		return *this;
	}
	//��ֵ������,����ӵ�һ������ֵ
	_SameClass& operator +=(const T& rhsVal)
	{
		//������valarray���������+=
		m_data += rhsVal;
		return *this;
	}
	//��ֵ������,�������һ��ͬ����ͬ�����ľ���
	_SameClass& operator -=(_SameClass& rhs)
	{
		Assert(this->m_row == rhs.GetRowNum());//�϶�_row����ʵ�ʾ�����ֵ,Assert�Ѿ�������
		Assert(this->m_col == rhs.GetColNum());//�϶�_col����ʵ�ʾ�����ֵ,Assert�Ѿ�������
		//������valarray���������-=
		m_data -= rhs.m_data;
		return *this;
	}
	//��ֵ������,�������һ������ֵ
	_SameClass& operator -=(const T& rhsVal)
	{
		//������valarray���������-=
		m_data -= rhsVal;
		return *this;
	}
	//��ֵ������,����˵�һ������������
	_SameClass& operator *=(_SameClass& rhs)
	{
		MatrixMultiply(*this, *this, rhs);	//���þ���˷�����
		return *this;
	}
	//��ֵ������,����˵�һ������ֵ
	_SameClass& operator *=(const T& rhsVal)
	{
		//������valarray���������*=
		m_data *= rhsVal;
		return *this;
	}
	//��ֵ������,�������һ������ֵ
	_SameClass& operator /=(const T& rhsVal)
	{
		if ( 0 == rhsVal)
		{
			std::cout<<"��ֵ����Ϊ0"<<std::endl;
			return *this;
		}
		//������valarray���������/=
		m_data /= rhsVal;
		return *this;
	}
	//һԪ+������������ֵȡ����ֵ
	_SameClass	operator +() const
	{
		_SameClass _tempMat(*this);
		_tempMat.m_data = abs(this->m_data);
		return _tempMat;
	}
	//һԪ����������CMatrix�ڲ�����ȡ��ֵ
	_SameClass	operator -() const
	{
		_SameClass _tempMat(*this);
		_tempMat.m_data = -_tempMat.m_data;
		return _tempMat;
	}
	//��Ԫ������
	//����Ӿ���
	_SameClass operator +(const _SameClass &rhs)
	{
		Assert(this->GetRowNum() == rhs.GetRowNum());
		Assert(this->GetColNum() == rhs.GetColNum());
		_SameClass _tempMat(*this);
		_tempMat.m_data += rhs.m_data;
		return _tempMat;
	}
	//����ӳ���
	_SameClass operator +(T &rhsVal)
	{
		_SameClass _tempMat(*this);
		_tempMat += rhsVal;
		return _tempMat;
	}
	//���������
	_SameClass operator -(const _SameClass &rhs)
	{
		Assert(this->GetRowNum() == rhs.GetRowNum());
		Assert(this->GetColNum() == rhs.GetColNum());
		_SameClass _tempMat(*this);
		_tempMat -= rhs;
		return _tempMat;
	}
	//���������
	_SameClass operator -(const T &rhsVal)
	{
		_SameClass _tempMat(*this);
		_tempMat -= rhsVal;
		return _tempMat;
	}
	//����˳���
	_SameClass operator *(const T &rhsVal)
	{
		_SameClass _tempMat(*this);
		_tempMat *= rhsVal;
		return _tempMat;
	}
	//����˾���
	_SameClass	operator *(const _SameClass &rhs)
	{
		if ( this->GetColNum() != rhs.GetRowNum() )
		{
			std::cerr<<"��˾���ά�Ȳ�һ��"<<std::endl;
		}
		_SameClass	_tempmat( this->GetRowNum(), rhs.GetColNum() );
		MatrixMultiply(_tempmat, *this, rhs);
		return _tempmat;
	}
	//���������
	_SameClass operator /(const T &rhsVal)
	{
		CMatrix _tempMat(*this);
		_tempMat /= rhsVal;
		return _tempMat;
	}
	//�����Ӿ���
	friend	_SameClass	operator + (const T &lhsVal, const _SameClass &rhs)
	{
		_SameClass _tempMat(rhs);
		_tempMat.m_data += lhsVal;
		return _tempMat;
	}
	//����������
	friend	_SameClass	operator - (const T &lhsVal, const _SameClass &rhs)
	{
		_SameClass _tempMat(rhs);
		_tempMat.m_data = lhsVal - _tempMat.m_data;
		return _tempMat;
	}
	//�����˾���
	friend	_SameClass	operator * (const T &lhsVal, const _SameClass &rhs)
	{
		_SameClass _tempMat(rhs);
		_tempMat.m_data *= lhsVal;
		return _tempMat;
	}
	//����Ƚ�
	//�Ƚ��������Ƿ���ȣ�����ȷ���true�����򷵻�false
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
	//�Ƚ��������Ƿ���ȣ�����ȷ���true�����򷵻�false
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
	//����ת��
	_SameClass operator!() 
	{ 
		_SameClass temp( m_col, m_row); 
		for(int i=0; i<m_col; i++) 
			for(int j=0; j<m_row; j++) 
				temp(i,j)=(*this)(j,i); 
		return temp; 
	} 
	//��Ҫ����ʹ�ã�������ʱ����printf�������������С�������ʱ��֤
	int Print();
};//class�������

typedef CMatrix<float> CMatrixf, CMatrix_f;
typedef CMatrix<double> CMatrixd, CMatrix_d;
typedef CMatrix<long double> CMatrixld, CMatrix_ld;

inline void error(char* v)
{
	std::cout<<v<<".program exited\n";
	getchar();
	exit(1);
}

//������������� 
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

//�������<< 
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

//����ֵ���ֵ����
template <class T>
T CMatrix<T>::maxnorm() const
{
	T norm=fabs(m_data[0]);
	for (int i=1;i<lenth;i++)
		norm=max(norm, fabs(m_data[i]));
	return norm;
}

//��Ԫ�صķ���,������
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


//����˷�����
template <class _TOut, class _Tlhs, class _Trhs>	//�������mOut��
CMatrix<_TOut>& MatrixMultiply(CMatrix<_TOut>& matOut, const CMatrix<_Tlhs>& matlhs, const CMatrix<_Trhs>& matrhs);

//ȫѡ��Ԫ�����������ʽ
template<class _Ty>
long double MatrixDetrminant( const CMatrix<_Ty>& _rhs );

//��ȫѡ��Ԫ��˹-Լ��(Gauss-Jordan)��ȥ������ʵ(��)����������
//�������ͱ����Ǹ�����
template <class _Ty>
int MatrixInversionGS( CMatrix<_Ty>& _rhs );

}//end namespace hbxdef

#include "MatrixCal.inl"



#endif