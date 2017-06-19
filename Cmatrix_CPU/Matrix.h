// Matrix.h		����ģ����ͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _MATRIX_H		//�����α���
#define _MATRIX_H

#include <valarray>		//ģ����valarray�ı�׼ͷ�ļ�
#include "Comm.h"		//����ͷ�ļ�
#include <math.h>		//��ѧͷ�ļ�

template <class _Ty>
class matrix			//������matrix
{
	typedef matrix<_Ty> _Myt;

 private:
	std::valarray<_Ty> m_Datas;	//����һά�������m_Datas

	size_t m_stRow;		//������������
	size_t m_stCol;		//������������

 public:
	typedef _Ty value_type;

	//���캯��һ(����1, 2�ֱ�Ϊ�������������)
/******
������matrix�Ĺ��캯��һ
	���캯���г��ֵ�m_DatasΪvalarray��Ķ�������stRow * stCol����Ԫ��
��Ԫ��û��ֵ�����������m_Datasʹ����valarray��Ĺ��캯����
				explicit valarray(size_t n)
    ��˽�б���m_stRow��m_stCol�ֱ𸳳�ֵstRow, stCol��
******/
	matrix(size_t stRow, size_t stCol)
		: m_Datas(stRow * stCol),
			m_stRow(stRow), m_stCol(stCol)
	{
		m_Datas.resize(GetRowNum() * GetColNum(), _Ty(0));
		/*
		���캯���г��ֵ�m_DatasΪvalarray��Ķ������ں������ڵ��õ�
		resize(size_t n, const T& c = T())��������������(����valarray�е�
		����)����һ�����������о�����*����ô��Ԫ�ظ����Ĵ洢�ռ䣬�ڶ���
		��������Щ����Ŀռ丳�����ģʽT�ĳ�ֵ0����������øú�������ȱʡ
		�����m_Datas����Ϊ0�����⣬���øú�������ʹ�õڶ����������򲻻��
		m_Datas���κ�һ��Ԫ�ظ���ֵ��
		*/
	}

	//���캯����(����1Ϊָ������ָ�룬����2, 3Ϊ������������)
/******
 ������matrix�Ĺ��캯����
     ��˽�б���m_stRow��m_stCol�ֱ𸳳�ֵstRow, stCol��
     ���������m_Datasʹ����valarray��Ĺ��캯����
				valarray(const _Ty *p, size_t n)
 m_Datas��ʼ���ĵ�һ����Ϊ����rhsָ�룬�ڶ�������Ϊrhs��Ԫ���ܸ�����
 ��rhs����*����
******/
	matrix(const _Ty* rhs, size_t stRow, size_t stCol)
		: m_Datas(rhs, stRow * stCol),
			m_stRow(stRow), m_stCol(stCol)
	{
	}

	//���캯����(�������캯��������Ϊ�Ծ���matrix������)
/******
 ������matrix�Ĺ��캯����
     �����þ���rhs���������m_Datas��ʼ��matrix����������m_Datas��
 �����þ���rhs������rhs.GetRowNum()������rhs.GetColNum()�ֱ��ʼ��˽
 �б���m_stRow��m_stCol
******/
	matrix(const _Myt& rhs)
		: m_Datas(rhs.m_Datas),
			m_stRow(rhs.GetRowNum()), m_stCol(rhs.GetColNum())
	{
	}

	size_t GetRowNum() const	//���ؾ��������ĺ���
	{
		return m_stRow;
	}

	size_t GetColNum() const	//���ؾ��������ĺ���
	{
		return m_stCol;
	}

	/*****
	    ���������()����ȡ��ά����ĳԪ����һά����m_Datas��λ��
	(matrix�ඨ��ľ���ʵ�����Ѿ���ת��Ϊһά����m_Datas)
	    ������Ԫ��(i, j)����������(��д)ʱ��ʹ��_Ty& operator ()
	������Ԫ��(i, j)��������ұ�(����)ʱ��ʹ��_Ty operator ()    
	*****/
	_Ty& operator () (size_t stRow, size_t stCol)
	{
		Assert(stRow < GetRowNum());	//�϶�stRow����ʵ�ʾ�����ֵ
		Assert(stCol < GetColNum());	//�϶�stCol����ʵ�ʾ�����ֵ

		return m_Datas[stRow * GetColNum() + stCol];
	}

/******
 ����()����������ض�ά������Ԫ��(stRow, stCol)��һά�����е�λ��
 stRow��stCol�ֱ�Ϊ����ĳԪ�������е�λ����
 ������Ԫ����������ұ�(����)ʱʹ��
******/
	const _Ty operator () (size_t stRow, size_t stCol) const
	{
		Assert(stRow < GetRowNum());	//�϶�stRow����ʵ�ʾ�����ֵ
		Assert(stCol < GetColNum());	//�϶�stCol����ʵ�ʾ�����ֵ

		return m_Datas[stRow * GetColNum() + stCol];
	}

	// ��ֵ������
	//�����������Է�*, +, -�����
	_Myt& operator += (const _Myt& rhs)		//�����������Է�+
	{
		Assert(GetRowNum() == rhs.GetRowNum());	//�϶���������������
		Assert(GetColNum() == rhs.GetColNum());	//�϶���������������
		//������valarray���壬ʹ�����m_Datas����������������
		m_Datas += rhs.m_Datas;
		//���������߾����У�����ָ����߾����ָ��
		return *this;
	}

	_Myt& operator -= (const _Myt& rhs)		//�����������Է�-
	{
		Assert(GetRowNum() == rhs.GetRowNum());
		Assert(GetColNum() == rhs.GetColNum());
		//������valarray���壬ʹ�����m_Datas����������������
		m_Datas -= rhs.m_Datas;
		//���������߾����У�����ָ����߾����ָ��
		return *this;
	}

	_Myt& operator *= (const _Myt& rhs)		//�����������Է�*
	{
		MatrixMultiply(*this, *this, rhs);	//���þ���˷�����

		return *this;						//�����������߾�����
	}

	//�����Է��ӡ��������ԡ�������
	_Myt& operator += (const _Ty& rhs)		//�����Լ���
	{
		m_Datas += rhs;			//������������ÿ��Ԫ�ؼ���
	
		return *this;			//�������ԭ����(����m_Datas)��
	}

	_Myt& operator -= (const _Ty& rhs)		//�����Լ���
	{
		m_Datas -= rhs;

		return *this;
	}

	_Myt& operator *= (const _Ty& rhs)		//�����Գ���
	{
		m_Datas *= rhs;
	
		return *this;
	}

	_Myt& operator /= (const _Ty& rhs)		//�����Գ�����
	{
		m_Datas /= rhs;

		return *this;
	}

	// һԪ������  �Ծ���(ÿ��Ԫ��)����+��-��
	_Myt operator + () const	//��+��
	{
		return *this;			//�������κδ���ά��ԭ״
	}

	_Myt operator - () const	//��-��
	{
		_Myt mat(*this);
		mat.m_Datas = -mat.m_Datas;		//ÿ��Ԫ�ظ�-��

		return mat;
	}

	// ��Ԫ������
	//�������	mat = lhs + rhs
	friend _Myt operator + (const _Myt& lhs, const _Ty& rhs)
	{
		_Myt mat(lhs);			//������һ�¾������mat
		mat.m_Datas += rhs;		//���¾������ÿ��Ԫ�ؼ���

		return mat;				//�����¾������
	}

	//�������	mat = lhs - rhs
	friend _Myt operator - (const _Myt& lhs, const _Ty& rhs)
	{
		_Myt mat(lhs);
		mat.m_Datas -= rhs;

		return mat;
	}

	//���������	mat = lhs * rhs
	friend _Myt operator * (const _Myt& lhs, const _Ty& rhs)
	{
		_Myt mat(lhs);			//������һ�¾������mat
		mat.m_Datas *= rhs;		//���¾������ÿ��Ԫ�س�����

		return mat;				//�����¾������
	}

	//���������	mat = lhs / rhs
	friend _Myt operator / (const _Myt& lhs, const _Ty& rhs)
	{
		_Myt mat(lhs);
		mat.m_Datas /= rhs;

		return mat;
	}

	//���Ӿ���	mat = lhs + rhs
	friend _Myt operator + (const _Ty& lhs, const _Myt& rhs)
	{
		_Myt mat(rhs);			//������һ�¾������mat
		mat.m_Datas += lhs;		//�������¾�������ÿ��Ԫ��

		return mat;
	}

	//��������	mat = lhs - rhs
	friend _Myt operator - (const _Ty& lhs, const _Myt& rhs)
	{
		_Myt mat(rhs);			//������һ�¾������mat
		mat.m_Datas -= lhs;		//�����¾�������ÿ��Ԫ��

		return mat;
	}

	//�����Ծ���	mat = lhs * rhs
	friend _Myt operator * (const _Ty& lhs, const _Myt& rhs)
	{
		_Myt mat(rhs);			//������һ�¾������mat
		mat.m_Datas *= lhs;		//���¾������ÿ��Ԫ�س�����

		return mat;
	}

	//����Ӿ���	mat = lhs + rhs
	friend _Myt operator + (const _Myt& lhs, const _Myt& rhs)
	{
		_Myt mat(lhs);					//������һ�¾������mat, ��������ʼ��
		mat.m_Datas += rhs.m_Datas;		//�����ұ߾������ÿ����ӦԪ��

		return mat;						//�����¾������
	}

	//���������	mat = lhs - rhs
	friend _Myt operator - (const _Myt& lhs, const _Myt& rhs)
	{
		_Myt mat(lhs);					//������һ�¾������mat, ��������ʼ��
		mat.m_Datas -= rhs.m_Datas;		//��ȥ�ұ߾������ÿ����ӦԪ��

		return mat;						//�����¾������
	}

	//������Ծ���	mTmp = lhs * rhs
	friend _Myt operator * (const _Myt& lhs, const _Myt& rhs)
	{	//����һ�������¶���mTmp		
		_Myt mTmp(lhs.GetRowNum(), rhs.GetColNum());	//û��ʼ��

		MatrixMultiply(mTmp, lhs, rhs);			//���þ���˷�����

		return mTmp;		//���¾�����󷵻�
	}

	//����Ƚ�
	//�Ƚ��������Ƿ���ȣ�����ȷ���true�����򷵻�false
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

	//�Ƚ��������Ƿ���ȣ�����ȷ���true�����򷵻�false
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

//����˷�����
template <class _Tyout, class _Tylhs, class _Tyrhs>	//�������mOut��
matrix<_Tyout>& MatrixMultiply(matrix<_Tyout>& mOut, const matrix<_Tylhs>& lhs, const matrix<_Tyrhs>& rhs);

//�������		һ��һ�����ȫ����
template <class _Ty>
void MatrixLinePrint(const matrix<_Ty>& mOut);

//�������		���ָ����ĳ��
template <class _Ty>
void MatrixLinePrint(const matrix<_Ty>& mOut, size_t LineNo);

//����ת��		ԭ����mIn��ת�ú�ľ�����mOut
template <class _Ty>
void MatrixTranspose(matrix<_Ty>& mIn, matrix<_Ty>& mOut);

//�жϾ���Գ�
template <class _Ty>		//�ԳƷ���true�����򷵻�false
bool MatrixSymmetry(const matrix<_Ty>& rhs);

//�жϾ���Գ�����
template <class _Ty>
int MatrixSymmetryRegular(const matrix<_Ty>& rhs, int sym);

//ȫѡ��Ԫ�����������ʽ
template <class _Ty>		//����ֵΪ����ʽֵ
long double MatrixDeterminant(const matrix<_Ty>& rhs);

//ȫѡ��Ԫ��˹(Gauss)��ȥ����һ��������
template <class _Ty>		//����ֵΪ����
size_t MatrixRank(const matrix<_Ty>& rhs);

//��ȫѡ��Ԫ��˹-Լ��(Gauss-Jordan)��ȥ������ʵ(��)����������
//�������ͱ����Ǹ�����
template <class _Ty>
int MatrixInversionGS(matrix<_Ty>& rhs);

//�á�����ѭ�����±�ŷ�������Գ�����������
//�������ͱ����Ǹ�����
template <class _Ty>
int MatrixSymmetryRegularInversion(matrix<_Ty>& rhs);

//������(Trench)�����в�����(Toeplitz)������
//�������ͱ����Ǹ�����
template <class _Ty>
int MatrixToeplitzInversionTrench(const valarray<_Ty>& t, 
								  const valarray<_Ty>& tt, matrix<_Ty>& rhs);

//ʵ����LU�ֽ�
//�������ͱ����Ǹ�����
template <class _Ty>
int MatrixLU(const matrix<_Ty>& rhs, matrix<_Ty>& lhs, matrix<_Ty>& uhs);

//�ú�˹�ɶ���(Householder)�任��һ��m*n�׵�ʵ�������QR�ֽ�
//�������ͱ����Ǹ�����
template <class _Ty>
int MatrixQR(matrix<_Ty>& rhs, matrix<_Ty>& rhq);

//�Գ������������˹��(Cholesky)�ֽ⼰��������ʽֵ
//�������ͱ����Ǹ�����
template <class _Ty>
long double MatrixSymmetryRegularCholesky(const matrix<_Ty>& rhs);

//һ��ʵ���������ֵ�ֽ�
//�������ͱ����Ǹ�����
template <class _Ty>
int MatrixSingularValue(matrix<_Ty>& a, matrix<_Ty>& u, 
												matrix<_Ty>& v, _Ty eps);

//�����������ֵ�ֽ�
//�������ͱ����Ǹ�����
template <class _Ty>
int GeneralizedInversionSingularValue(matrix<_Ty>& a, matrix<_Ty>& aa, 
								_Ty eps, matrix<_Ty>& u, matrix<_Ty>& v);

#include "Matrix.inl"		//�༰��غ����Ķ���ͷ�ļ�

#endif		// _MATRIX_H

