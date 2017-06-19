#ifndef _MATRIXCAL_INL
#define _MATRIXCAL_INL

namespace HBXDef
{
	//矩阵乘法函数, 最后结果在mOut中
	template <class _TOut, class _Tlhs, class _Trhs>
	CMatrix<_TOut>& MatrixMultiply(CMatrix<_TOut>& matOut, const CMatrix<_Tlhs>& matlhs, const CMatrix<_Trhs>& matrhs)
	{
		//判断左边矩阵的列数与右边矩阵的行数相等
		Assert(matlhs.GetColNum() == matrhs.GetRowNum());
		//生成矩阵新对象，用lhs的行作为新阵的行数，用rhs的列数作为新阵的列数
		CMatrix<_TOut> _matTemp(matlhs.GetRowNum(), matrhs.GetColNum());

		for(size_t i = 0; i < _matTemp.GetRowNum(); i ++)
		{
			for(size_t j = 0; j < _matTemp.GetColNum(); j ++)
			{
				_matTemp(i, j) = _TOut(0);		//赋初值0
				for(size_t k = 0; k < matlhs.GetColNum(); k ++)
				{
					_matTemp(i, j) += matlhs(i, k) * matrhs(k, j);
				}
			}
		}
		matOut = _matTemp;		//将最后结果转放入mOut矩阵中
		return matOut;		//返回结果矩阵mOut
	}


	//全选主元法求矩阵行列式
	template<class _Ty>
	long double MatrixDetrminant( const CMatrix<_Ty>& _rhs )
	{
		long double MaxValue, tmp;
		size_t stRank = _rhs.GetColNum();
		//若不为方阵则返回0
		if ( stRank != _rhs.GetRowNum() )
		return long double(0);
		
		CMatrix< long double > m( stRank, stRank );

		for(size_t i=0; i<stRank; i++)
			for(size_t j=0; j<stRank; j++)
				m(i, j) = (long double)_rhs(i, j);	//初始化

		size_t iSign, jSign;				// 主元的位置标志

		long double Det(1);					// 行列式的值
		int	nSgn = 1;						// 符号

		for ( size_t k=0; k<stRank-1; k++ )
		{
			MaxValue = 0.0;
			for(size_t i = k; i < stRank; i ++)
			{
				for(size_t j = k; j < stRank; j ++)//对角线以上
				{
					tmp = fabsl(m(i, j));		//求m(i,j)绝对值
					if(tmp > MaxValue)
					{
						MaxValue = tmp;
						iSign = i;			//记下主元位置
						jSign = j;
					}
				}
			}
			//判断主元是否为0
			if ( LongDoubleEqual(MaxValue, 0) )
				return long double(0);
			
			if ( iSign != k )//绝对值最大元素不在当前行
			{
				nSgn = -nSgn;
				for (size_t j = k; j < stRank; j++)
					swap( m(k,j), m(iSign,j) );
			}

			if(jSign != k)			//绝对值最大元素不在当前列
			{
				nSgn = -nSgn;		//变换行列式符号
				for(size_t i = k; i < stRank; i ++)		//交换列
					swap(m(i, jSign), m(i, k));
			}

			Det *= m(k, k);					//对角元相乘
			for(size_t i = k + 1; i < stRank; i ++)
			{
				long double d(m(i, k) / m(k, k));		//消元因子
				for(size_t j = k + 1; j < stRank; j ++)	//将主元下方元素消为0
					m(i, j) -= d * m(k, j);	//当前主元行下行其余元素作变换
			}
		}

		return Det * nSgn * m(stRank - 1, stRank - 1);	//返回行列式数值
	}


	//全选主元高斯-约当(Gauss-Jordan)法求矩阵逆
	template <class _Ty>
	int MatrixInversionGS( CMatrix<_Ty>& _rhs )
	{
		size_t stRank = _rhs.GetColNum();	// 矩阵阶数
		if(stRank != _rhs.GetRowNum())
			return int(-1);					//rhs不是方阵

		std::valarray<size_t> is(stRank);		//行交换信息
		std::valarray<size_t> js(stRank);		//列交换信息

		CMatrix<_Ty> m(_rhs);					//生成一matrix对象

		for(size_t k = 0; k < stRank; k++)
		{	// 全选主元
			long double MaxValue(0);
			for(size_t i = k; i < stRank; i ++)
			{
				for(size_t j = k; j < stRank; j ++)
				{
					long double tmp(fabsl(m(i, j)));	//求m(i,j)绝对值
					if(tmp > MaxValue)				//主元不在对角线上
					{
						MaxValue = tmp;
						is[k] = i;					//记录主元行数
						js[k] = j;					//记录主元列数
					}
				}
			}

			if( LongDoubleEqual(MaxValue, 0) ) 
				return int(0);						//主元为0，矩阵奇异

			if(is[k] != k)							//主元不在当前行
			{
				for(size_t j = 0; j < stRank; j ++)	//交换行元素
					swap(m(k, j), m(is[k], j));
			}

			if(js[k] != k)							//主元不在当前列
			{
				for(size_t i = 0; i < stRank; i ++)	//交换列元素
					swap(m(i, k), m(i, js[k]));
			}

			m(k, k) = 1.0 / m(k, k);				//主元倒数
			for(size_t j = 0; j < stRank; j ++)
				if(j != k)
					m(k, j) *= m(k, k);

			for(size_t i = 0; i < stRank; i ++)
				if(i != k)
					for(size_t j = 0; j < stRank; j ++)	
						if(j != k)
							m(i, j) = m(i, j) - m(i, k) * m(k, j);

			for(size_t i = 0; i < stRank; i ++)
				if(i != k)
					m(i, k) = -m(i, k) * m(k, k);
		}

		for(int r = stRank - 1; r >= 0; r--)
		{
			if(js[r] != r)
				for(size_t j = 0; j < stRank; j ++)
					swap(m(r, j), m(js[r], j));
			if(is[r] != r)
				for(size_t i = 0; i < stRank; i ++)
					swap(m(i, r), m(i, is[r]));
		}

		_rhs = m;

		return int(1);
	}

}//namespace hbxdef

#endif //结束_MATRIXCAL_INL