// Matrix.inl		矩阵模板类函数(方法)定义
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31

#ifndef _MATRIX_INL
#define _MATRIX_INL
typedef unsigned int size_t;
//矩阵乘法函数
template <class _Tyout, class _Tylhs, class _Tyrhs>	//最后结果在mOut中
matrix<_Tyout>& MatrixMultiply(matrix<_Tyout>& mOut, const matrix<_Tylhs>& lhs, const matrix<_Tyrhs>& rhs)
{	//断定左边矩阵的列数与右边矩阵的行数相等
	Assert(lhs.GetColNum() == rhs.GetRowNum());
	//生成矩阵新对象，用lhs的行作为新阵的行数，用rhs的列数作为新阵的列数
	matrix<_Tyout> mTmp(lhs.GetRowNum(), rhs.GetColNum());

	for(size_t i = 0; i < mTmp.GetRowNum(); i ++)
	{
		for(size_t j = 0; j < mTmp.GetColNum(); j ++)
		{
			mTmp(i, j) = _Tyout(0);		//赋初值0
			for(size_t k = 0; k < lhs.GetColNum(); k ++)
			{
				mTmp(i, j) += lhs(i, k) * rhs(k, j);
			}
		}
	}

	mOut = mTmp;		//将最后结果转放入mOut矩阵中
	return mOut;		//返回结果矩阵mOut
}

//输出矩阵函数		按一行一行进行输出
template <class _Ty>
void MatrixLinePrint(const matrix<_Ty>& mOut)
{	
	size_t sR, sC;
	sR=mOut.GetRowNum();
	sC=mOut.GetColNum();

	for(size_t stR=0; stR<mOut.GetRowNum(); stR++)
	{
		for(size_t stC=0; stC<mOut.GetColNum(); stC++)
		{
			cout.width(15);				//元素对齐，让每个元素占15列
			cout << mOut(stR, stC) << ' ';
		}
		cout << endl;
	}
}

//输出矩阵函数		按指定行进行输出
template <class _Ty>
void MatrixLinePrint(const matrix<_Ty>& mOut, size_t LineNo)
{	
	size_t sR, sC;

	sR=mOut.GetRowNum();
	sC=mOut.GetColNum();
	
	for(size_t stC=0; stC<mOut.GetColNum(); stC++)
	{
		cout.width(15);				//元素对齐，让每个元素占15列
		cout << mOut(LineNo, stC) << ' ';
	}
	cout << endl;
}

//矩阵转置	==	原阵在mIn，转置后的矩阵在mOut  ==
template <class _Ty>
void MatrixTranspose(matrix<_Ty>& mIn, matrix<_Ty>& mOut)
{
	size_t sR, sC;
	sR = mIn.GetRowNum();			//取原矩阵行数
	sC = mIn.GetColNum();			//取原矩阵列数

	matrix<_Ty> mTemp(sC, sR);		//生成一新阵，行数与列数与原阵互换
	
	for(size_t stC=0; stC<sC; stC++)
		for(size_t stR=0; stR<sR; stR++)
			mTemp(stC, stR) = mIn(stR, stC);	//对新阵赋值

	mOut = mTemp;					//返回新的转置阵
}

//判断矩阵对称
template <class _Ty>	
bool MatrixSymmetry(const matrix<_Ty>& rhs)
{
	bool bSy = true;
	size_t stRow = rhs.GetRowNum();	//取矩阵行数

	if(rhs.GetColNum() == stRow)	// 必须是方阵
	{
		for(size_t i = 1; i < stRow; i ++)			//判断是否对称
			for(size_t j = 0; j < i; j ++)
				if(FloatNotEqual((long double)rhs(i, j), (long double)rhs(j, i)))
				{
					bSy =  false;
					goto RET;
				}
	}
	else
		bSy = false;
	RET: return bSy;	//矩阵对称
}

//判断矩阵对称正定	
template <class _Ty>	
int MatrixSymmetryRegular(const matrix<_Ty>& rhs, int sym)
{							
	long double ldDet;
	size_t i, j, k;

	size_t sC = rhs.GetColNum();	//矩阵列数
	size_t sR = rhs.GetRowNum();	//矩阵行数

	size_t stRank = sR;				// 矩阵阶数
	if(stRank != rhs.GetRowNum())
		return int(-1);				// 不是方阵

	if(sym > 0)
		if(MatrixSymmetry(rhs)==false)
			return int(-2);			//rhs不是对称阵

	cout << " K = 1 \t Determinant = " << rhs(0, 0) <<endl;
	
	for(k = 0; k < stRank; k ++)	//若要判别半正定，负定，这句要修改
	{
		if(FloatEqual(rhs(k, k), 0)||rhs(k, k) < 0)
			return int(-3);	//对角元不大于0，矩阵不是正定阵
	}

	for(k = 2; k <= sR; k++)
	{
		matrix<long double> m(k, k);	//生成一matrix对象

		for(i=0; i<k; i++)
			for(j=0; j<k; j++)
				m(i, j) = (long double)rhs(i, j);	//初始化

		ldDet = MatrixDeterminant(m);				// 顺序主子式的值

		cout << " K = " << k << "\t Determinant = " << ldDet << endl; 
		
		if(FloatEqual(ldDet,0) || ldDet < 0.0)
			return (0);					//不是正定阵
	}
	if(sym == 1) return int(2);			//矩阵为正定对称阵
	else		 return int(1);			//矩阵为正定阵
}

//全选主元法求矩阵行列式函数
template <class _Ty>
long double MatrixDeterminant(const matrix<_Ty>& rhs)		
{	
	long double  MaxValue, tmp;
	size_t stRank = rhs.GetColNum();	// 矩阵阶数
	if(stRank != rhs.GetRowNum())
		return long double(0);			//rhs不是方阵

	matrix<long double> m(stRank, stRank);	//生成一matrix对象

	for(size_t i=0; i<stRank; i++)
		for(size_t j=0; j<stRank; j++)
			m(i, j) = (long double)rhs(i, j);	//初始化

	size_t iSign, jSign;				// 主元的位置标志

	long double Det(1);					// 行列式的值

	int	nSgn = 1;						// 符号

	for(size_t k = 0; k < stRank - 1; k ++)	// 全选主元
	{	
		MaxValue = 0.0;
		for(i = k; i < stRank; i ++)
		{
			for(size_t j = k; j < stRank; j ++)
			{
				tmp = Abs(m(i, j));		//求m(i,j)绝对值
				if(tmp > MaxValue)
				{
					MaxValue = tmp;
					iSign = i;			//记下主元位置
					jSign = j;
				}
			}
		}
		
		if(FloatEqual(MaxValue, 0))	//绝对值最大元素为0，行列式也为0
			return long double(0);

		if(iSign != k)			//绝对值最大元素不在当前行
		{
			nSgn = -nSgn;		//变换行列式符号
			for(size_t j = k; j < stRank; j ++)		//交换行
				swap(m(k, j), m(iSign, j));
		}

		if(jSign != k)			//绝对值最大元素不在当前列
		{
			nSgn = -nSgn;		//变换行列式符号
			for(size_t i = k; i < stRank; i ++)		//交换列
				swap(m(i, jSign), m(i, k));
		}

		Det *= m(k, k);					//对角元相乘
		for(i = k + 1; i < stRank; i ++)
		{
			long double d(m(i, k) / m(k, k));		//消元因子
			for(size_t j = k + 1; j < stRank; j ++)	//将主元下方元素消为0
				m(i, j) -= d * m(k, j);	//当前主元行下行其余元素作变换
		}
	}

	return Det * nSgn * m(stRank - 1, stRank - 1);	//返回行列式数值
}

//全选主元高斯(Gauss)消去法求一般矩阵的秩
template <class _Ty>						//返回值为秩数
size_t MatrixRank(const matrix<_Ty>& rhs)
{
	size_t iSign, jSign;					//主元的位置标志
	size_t mRank = 0;						//矩阵秩数
	size_t stRow = rhs.GetRowNum();			//取矩阵行数
	size_t stCol = rhs.GetColNum();			//取矩阵列数
	size_t ColRowMin = Min(stRow, stCol);	//取行或列最小值

	matrix<_Ty> m(rhs);				//生成一matrix对象，用rhs初始化

	for(size_t k = 0; k < ColRowMin; k ++)
	{	// 全选主元
		long double MaxValue(0);
		for(size_t i = k; i < stRow; i ++)
		{
			for(size_t j = k; j < stCol; j ++)
			{
				long double tmp(Abs(m(i, j)));	//求m(i,j)绝对值
				if(tmp > MaxValue)
				{
					MaxValue = tmp;
					iSign = i;					//记下主元位置
					jSign = j;
				}
			}
		}
		
		//子阵元素绝对值最大者为0，	注意浮点数与0相等的定义，见comm.h
		if(FloatEqual(MaxValue, 0)) 
			break;	//return mRank;
		else				
			mRank++;		//子阵元素绝对值最大者不为0，矩阵秩加1

		if(k ==(ColRowMin - 1))		//已到最后一行(列)
			break;	//return mRank;
		
		if(iSign != k)				//主元不在当前行
		{
			for(size_t j = k; j < stCol; j ++)	//交换行元素
				swap(m(k, j), m(iSign, j));
		}

		if(jSign != k)				//主元不在当前列
		{
			for(size_t i = k; i < stRow; i ++)	//交换列元素
				swap(m(i, jSign), m(i, k));
		}

		for(i = k + 1; i < stRow; i ++)
		{
			const _Ty d(m(i, k) / m(k, k));		//消元因子
			for(size_t j = k + 1; j < stCol; j ++)	
				m(i, j) -= d * m(k, j);		//当前主元右下阵元素作变换
		}
	}
	return mRank;
}

//全选主元高斯-约当(Gauss-Jordan)法求矩阵逆
template <class _Ty>
int MatrixInversionGS(matrix<_Ty >& rhs)
{
	size_t stRank = rhs.GetColNum();	// 矩阵阶数
	if(stRank != rhs.GetRowNum())
		return int(-1);					//rhs不是方阵

	valarray<size_t> is(stRank);		//行交换信息
	valarray<size_t> js(stRank);		//列交换信息
	
	matrix<_Ty> m(rhs);					//生成一matrix对象

	for(size_t k = 0; k < stRank; k++)
	{	// 全选主元
		long double MaxValue(0);
		for(size_t i = k; i < stRank; i ++)
		{
			for(size_t j = k; j < stRank; j ++)
			{
				long double tmp(Abs(m(i, j)));	//求m(i,j)绝对值
				if(tmp > MaxValue)				//主元不在对角线上
				{
					MaxValue = tmp;
					is[k] = i;					//记录主元行数
					js[k] = j;					//记录主元列数
				}
			}
		}
		
		if(FloatEqual(MaxValue, 0)) 
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

		for(i = 0; i < stRank; i ++)
			if(i != k)
				for(size_t j = 0; j < stRank; j ++)	
					if(j != k)
						m(i, j) = m(i, j) - m(i, k) * m(k, j);

		for(i = 0; i < stRank; i ++)
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

	rhs = m;

	return int(1);
}

//用“变量循环重新编号法”求对称正定矩阵逆
//矩阵类型必须是浮点型
template <class _Ty>
int MatrixSymmetryRegularInversion(matrix<_Ty>& rhs)
{
	int iSym = MatrixSymmetryRegular(rhs, 1);	//判别是否对称正定
	if(iSym < 2)
		return (iSym);			//rhs不是对称正定阵

	size_t stRank = rhs.GetColNum();	// 矩阵阶数
	
	matrix<_Ty> msr(rhs);			//生成一matrix对象，用rhs初始化

	valarray<_Ty> b(stRank);
	
	for(size_t k=0; k<stRank; k++)
	{
		_Ty w= msr(0, 0);
		size_t m = stRank - k -1;
		for(size_t i = 1; i < stRank; i++)
		{
			_Ty g = msr(i, 0);
			b[i] = g / w;
			if(i <= m) b[i] = -b[i];
			for(size_t j = 1; j <= i; j ++)
				msr((i-1),(j-1)) = msr(i, j) + g * b[j];
		}
		msr(stRank-1, stRank-1) = 1.0 / w;
		for(i = 1; i < stRank; i ++)
			msr(stRank-1,(i-1)) =  b[i];
	}
	for(size_t i = 0; i < stRank-1; i ++)
		for(size_t j = i+1; j < stRank; j ++)
			msr(i,j) = msr(j, i);

	rhs = msr;

	return (iSym);
}


//特兰持(Trench)法求托伯利兹(Toeplitz)矩阵逆
//矩阵类型必须是浮点型
template <class _Ty>
int MatrixToeplitzInversionTrench(const valarray<_Ty>& t, const valarray<_Ty>& tuo, matrix<_Ty>& rhs)
{
	size_t stRank = rhs.GetColNum();	// 矩阵阶数
	if(stRank != rhs.GetRowNum())
		return int(-1);					//rhs不是方阵

	if(FloatEqual(t[0], 0))
		return int(0);	

	valarray<_Ty> c(stRank);
	valarray<_Ty> r(stRank);
	valarray<_Ty> p(stRank);

	_Ty a=t[0]; 
	c[0]=tuo[1]/t[0];
	r[0]=t[1]/t[0];

	matrix<_Ty> b(rhs);
	for(size_t k=0; k<=stRank-3; k++)
	{
		_Ty s=0.0;
		for(size_t j=1; j<=k+1; j++)
			s=s+c[k+1-j]*tuo[j];
		s=(s-tuo[k+2])/a;
		for(size_t i=0; i<=k; i++)
			p[i]=c[i]+s*r[k-i];
        c[k+1]=-s;
        s=0.0;
        for(j=1; j<=k+1; j++)
			s=s+r[k+1-j]*t[j];
        s=(s-t[k+2])/a;
        for(i=0; i<=k; i++)
        {
			r[i]=r[i]+s*c[k-i];
            c[k-i]=p[k-i];
        }
        r[k+1]=-s;
		a=0.0;
        for(j=1; j<=k+2; j++)
			a=a+t[j]*c[j-1];
        a=t[0]-a;
		if(FloatEqual(a, 0))
			return int(0);	
    }
    b(0,0)=1.0/a;
    for(size_t i=0; i<stRank-1; i++)
    {	
		k=i+1;
        b(0, k)=-r[i]/a;
		b(i+1,0)=-c[i]/a;
    }
    for(i=0; i<stRank-1; i++)
    for(size_t j=0; j<stRank-1; j++)
        b(i+1, j+1)=b(i,j)-c[i]*b(0,j+1)+c[stRank-j-2]*b(0,stRank-i-1);

	rhs = b;
    
    return int(1);
}

//实矩阵LU分解
//矩阵必须是浮点型
template <class _Ty>
int MatrixLU(const matrix<_Ty>& rhs, matrix<_Ty>& lhs, matrix<_Ty>& uhs)
{
	size_t stRank = rhs.GetColNum();	// 矩阵阶数
	if(stRank != rhs.GetRowNum())
		return int(-1);			//rhs不是方阵
	
	matrix<_Ty> m(rhs);
	for(size_t k=0; k<stRank-1; k++)
	{
		if(FloatEqual(m(k,k),0))
			return int(0);		//主元为0
		for(size_t i=k+1; i<stRank; i++)
			m(i,k) /= m(k,k);
		for(i=k+1; i<stRank; i++)
			for(size_t j=k+1; j<stRank; j++)
				m(i,j)=m(i,j)-m(i,k)*m(k,j);
	}
	//给上、下三角阵赋值
	for(size_t i=0; i<stRank; i++)
	{
		for(size_t j=0; j<i; j++)
		{
			lhs(i,j)=m(i,j);
			uhs(i,j)=0.0;
		}
		lhs(i,i)=1.0;
		uhs(i,i)=m(i,i);
		for(j=i+1; j<stRank; j++)
		{
			lhs(i,j)=0.0;
			uhs(i,j)=m(i,j);
		}
	}
	
	return (1);		//分解成功
}

//用豪斯荷尔德(Householder)变换对一般m*n阶的实矩阵进行QR分解
template <class _Ty>
int MatrixQR(matrix<_Ty>& rhs, matrix<_Ty>& rhq)
{
	size_t stRow = rhs.GetRowNum();	// 矩阵行数
	size_t stCol = rhs.GetColNum();	// 矩阵列数
	
	if(stRow < stCol)
		return (0);					//行不能小于列

	for(size_t i=0; i<stRow; i++)
    for(size_t j=0; j<stRow; j++)
    { 
		rhq(i,j)=0.0;
        if(i==j) rhq(i,j)=1.0;
    }

	size_t nn=stCol;
    
	if(stRow == stCol) nn=stRow-1;

	for(size_t k=0; k<nn; k++)
    {
		_Ty u=0.0;

        for(size_t i = k; i < stRow; i++)
        { 
			_Ty w = Abs(rhs(i,k));
            if(w > u) u = w;
        }
        _Ty alpha=0.0;
        for(i = k; i < stRow; i++)
        {
			_Ty t=rhs(i,k)/u;
			alpha=alpha+t*t;
		}
        
		if(rhs(k,k)>0.0) u=-u;
        
		alpha=u*sqrt(alpha);
        
		if(FloatEqual(alpha,0.0)) return(0);

        u=sqrt(2.0*alpha*(alpha-rhs(k,k)));

        if(FloatNotEqual(u,0.0))
        { 
			rhs(k,k)=(rhs(k,k)-alpha)/u;
            
			for(i=k+1; i<stRow; i++)
            	rhs(i,k) /= u;

			for(size_t j=0; j<stRow; j++)
            {
				_Ty t=0.0;

                for(size_t jj=k; jj<stRow; jj++)
					t=t+rhs(jj,k)*rhq(jj,j);

                for(i=k; i<stRow; i++)
					rhq(i,j)=rhq(i,j)-2.0*t*rhs(i,k);
            }
            
			for(j=k+1; j<stCol; j++)
            { 
				_Ty t=0.0;
            
				for(size_t jj=k; jj<stRow; jj++)
					t=t+rhs(jj,k)*rhs(jj,j);

                for(i=k; i<stRow; i++)
            		rhs(i,j)=rhs(i,j)-2.0*t*rhs(i,k);
			}

            rhs(k,k)=alpha;

            for(i=k+1; i<stRow; i++)
              rhs(i,k)=0.0;
          }
    }

	for(i=0; i<stRow-1; i++)
		for(size_t j=i+1; j<stRow;j++)
				swap(rhq(i,j), rhq(j,i));
	
    return (1);

}

//对称正定阵的乔里斯基(Cholesky)分解及求其行列式值
//矩阵与返回值类型必须是浮点型
template <class _Ty>
long double MatrixSymmetryRegularCholesky(matrix<_Ty>& rhs)
{
	int iReValue= MatrixSymmetryRegular(rhs, 1);	//判别对称正定
	if(iReValue < 2)
		return long double(0);		//rhs不是对称正定阵

	size_t stRank = rhs.GetColNum();	// 矩阵阶数

	matrix<_Ty> m(rhs);				//生成一matrix对象，用rhs初始化

	long double Det = m(0,0);					//计算行列式值
	m(0,0) = sqrt(m(0,0)); 
	
	for(size_t i=1; i<stRank; i++)
		m(i, 0) /= m(0, 0);
	for(size_t j=1; j<stRank; j++)
	{
		for(size_t k=0; k<j; k++)
			m(j,j) = m(j,j) - m(j,k) * m(j,k);
		Det *= m(j,j);
		m(j,j) = sqrt(m(j,j));
		for(i=j+1; i<stRank; i++)
		{
			for(k=0; k<j; k++)
				m(i,j) = m(i,j) -m(i,k) * m(j,k);
			m(i,j) /= m(j,j);
		}
	}
	for(i=0; i<stRank-1; i++)
	    for(j=i+1; j<stRank; j++)
			m(i,j)=0;

	rhs = m;		//返回Cholesky阵，原矩阵将被复盖
	return Det;		//返回行列式值
}

//一般实矩阵的奇异值分解
template <class _Ty>
int MatrixSingularValue(matrix<_Ty>& a, matrix<_Ty>& u, 
											matrix<_Ty>& v, _Ty eps)
{
	int i, it(60), kk, mm, nn, m1, ks, ka,j;
    _Ty d,dd,t,sm,sm1,em1,sk,ek,b,c,shh;

	int m = a.GetRowNum();
	int n = a.GetColNum();

	for(int j=0; j<m; j++) u(j, m-1) = _Ty(0);

	if(m > n) ka = m + 1;
	else ka = n + 1;
	
	valarray<_Ty> s(ka), e(ka), w(ka), fg(2), cs(2);
    
	int k = n;
    if(m-1<n) k=m-1;
    int l = m;
    if(n-2<m) l=n-2;
    if(l<0) l=0;
    int ll=k;
    if(l>k) ll=l;
    if(ll>=1)
    { 
		for(kk=1; kk<=ll; kk++)
        {
			if(kk<=k)
            {
				d=0.0;
                for(i=kk; i<=m; i++)
					d=d+a(i-1,kk-1)*a(i-1,kk-1);
                s[kk-1]=sqrt(d);
                if(FloatNotEqual(s[kk-1],0.0))
                {
                    if(FloatNotEqual(a(kk-1,kk-1),0.0))
                    {
						s[kk-1]=Abs(s[kk-1]);
                        if(a(kk-1,kk-1)<0.0) s[kk-1]=-s[kk-1];
                    }
                    for(i=kk; i<=m; i++)	a(i-1,kk-1)=a(i-1,kk-1)/s[kk-1];
                    a(kk-1,kk-1)=1.0+a(kk-1,kk-1);
                }
                s[kk-1]=-s[kk-1];
            }
            if(n>=kk+1)
            { 
				for(j=kk+1; j<=n; j++)
                {
					if((kk<=k)&&(FloatNotEqual(s[kk-1],0.0)))
                    {
						d=0.0;
                        for(i=kk; i<=m; i++)	d=d+a(i-1,kk-1)*a(i-1,j-1);
                        d=-d/a(kk-1,kk-1);
                        for(i=kk; i<=m; i++)	a(i-1,j-1)=a(i-1,j-1)+d*a(i-1,kk-1);
                    }
                    e[j-1]=a(kk-1,j-1);
                }
            }
            if(kk<=k)
				for(i=kk; i<=m; i++)	u(i-1,kk-1)=a(i-1,kk-1);
            if(kk<=l)
            { 
				d=0.0;
                for(i=kk+1; i<=n; i++)
                  d=d+e[i-1]*e[i-1];
                e[kk-1]=sqrt(d);
                if(FloatNotEqual(e[kk-1],0.0))
                {
					if(FloatNotEqual(e[kk],0.0))
                    {
						e[kk-1]=Abs(e[kk-1]);
                        if(e[kk]<0.0) e[kk-1]=-e[kk-1];
                    }
                    for(i=kk+1; i<=n; i++)
                      e[i-1]=e[i-1]/e[kk-1];
                    e[kk]=1.0+e[kk];
                }
                e[kk-1]=-e[kk-1];
                if((kk+1<=m)&&(FloatNotEqual(e[kk-1],0.0)))
                { 
					for(i=kk+1; i<=m; i++) w[i-1]=0.0;
                    for(j=kk+1; j<=n; j++)
                      for(i=kk+1; i<=m; i++)
                        w[i-1]=w[i-1]+e[j-1]*a(i-1,j-1);
                    for(j=kk+1; j<=n; j++)
                      for(i=kk+1; i<=m; i++)
                          a(i-1,j-1)=a(i-1,j-1)-w[i-1]*e[j-1]/e[kk];
                }
                for(i=kk+1; i<=n; i++)
                  v(i-1,kk-1)=e[i-1];
            }
        }
    }
    mm=n;
    if(m+1<n) mm=m+1;
    if(k<n) s[k]=a(k,k);
    if(m<mm) s[mm-1]=0.0;
    if(l+1<mm) e[l]=a(l,mm-1);
    e[mm-1]=0.0;
    nn=m;
    if(m>n) nn=n;
    if(nn>=k+1)
    {
		for(j=k+1; j<=nn; j++)
        {
			for(i=1; i<=m; i++)	u(i-1,j-1)=0.0;
            u(j-1,j-1)=1.0;
        }
    }
    if(k>=1)
    { 
		for(ll=1; ll<=k; ll++)
        {
			kk=k-ll+1;
            if(s[kk-1]!=0.0)
            {
				if(nn>=kk+1)
                  for(j=kk+1; j<=nn; j++)
                  {
					  d=0.0;
                      for(i=kk; i<=m; i++)
						d=d+u(i-1,kk-1)*u(i-1,j-1)/u(kk-1,kk-1);
                      d=-d;
                      for(i=kk; i<=m; i++)
                          u(i-1,j-1)=u(i-1,j-1)+d*u(i-1,kk-1);
                  }
                  for(i=kk; i<=m; i++)
					  u(i-1,kk-1)=-u(i-1,kk-1);
                  u(kk-1,kk-1)=1.0+u(kk-1,kk-1);
                  if(kk-1>=1)
                    for(i=1; i<kk; i++)	u(i-1,kk-1)=0.0;
            }
            else
            { 
				for(i=1; i<=m; i++)	u(i-1,kk-1)=0.0;
                u(kk-1,kk-1)=1.0;
            }
        }
    }
    for(ll=1; ll<=n; ll++)
    { 
		kk=n-ll+1; 
        if((kk<=l)&&(e[kk-1]!=0.0))
        {
			for(j=kk+1; j<=n; j++)
            {
				d=0.0;
                for(i=kk+1; i<=n; i++)
                    d=d+v(i-1,kk-1)*v(i-1,j-1)/v(kk,kk-1);
                d=-d;
                for(i=kk+1; i<=n; i++)
                    v(i-1,j-1)=v(i-1,j-1)+d*v(i-1,kk-1);
            }
        }
        for(i=1; i<=n; i++)	v(i-1,kk-1)=0.0;
        v(kk-1,kk-1)=1.0;
    }
    for(i=1; i<=m; i++)
		for(j=1; j<=n; j++)	a(i-1,j-1)=0.0;
    m1=mm; 
	it=60;
    while(1)
    { 
		if(mm==0)
        {
			_MSV_1(a,e,s,v,m,n);
			return(1);
        }
        if(it==0)
        { 
			_MSV_1(a,e,s,v,m,n);
			return(0);
        }
        kk=mm-1;
	while((kk!=0)&&(FloatNotEqual(e[kk-1],0.0)))
    { 
		d=Abs(s[kk-1])+Abs(s[kk]);
        dd=Abs(e[kk-1]);
        if(dd>eps*d) kk=kk-1;
        else e[kk-1]=0.0;
    }
    if(kk==mm-1)
    {
		kk=kk+1;
        if(s[kk-1]<0.0)
        {
			s[kk-1]=-s[kk-1];
            for(i=1; i<=n; i++)		v(i-1,kk-1)=-v(i-1,kk-1);
        }
        while((kk!=m1)&&(s[kk-1]<s[kk]))
        { 
			d=s[kk-1]; 
			s[kk-1]=s[kk];
			s[kk]=d;
            if(kk<n)
			    for(i=1; i<=n; i++)
                {
                  d=v(i-1,kk-1);
				  v(i-1,kk-1)=v(i-1,kk); 
				  v(i-1,kk)=d;
                }
             if(kk<m)
                  for(i=1; i<=m; i++)
                  { 
                      d=u(i-1,kk-1);
					  u(i-1,kk-1)=u(i-1,kk);
					  u(i-1,kk)=d;
                  }
                kk=kk+1;
            }
            it=60;
            mm=mm-1;
        }
        else
        { 
			ks=mm;
            while((ks>kk)&&(Abs(s[ks-1])!=0.0))
            {
				d=0.0;
                if(ks!=mm) d=d+Abs(e[ks-1]);
                if(ks!=kk+1) d=d+Abs(e[ks-2]);
                dd=Abs(s[ks-1]);
                if(dd>eps*d) ks=ks-1;
                else s[ks-1]=0.0;
            }
            if(ks==kk)
            { 
				kk=kk+1;
                d=Abs(s[mm-1]);
                t=Abs(s[mm-2]);
                if(t>d) d=t;
                t=Abs(e[mm-2]);
                if(t>d) d=t;
                t=Abs(s[kk-1]);
                if(t>d) d=t;
                t=Abs(e[kk-1]);
                if(t>d) d=t;
                sm=s[mm-1]/d; sm1=s[mm-2]/d;
                em1=e[mm-2]/d;
                sk=s[kk-1]/d; ek=e[kk-1]/d;
                b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                c=sm*em1; c=c*c; shh=0.0;
                if((b!=0.0)||(c!=0.0))
                {
					shh=sqrt(b*b+c);
                    if(b<0.0) shh=-shh;
                    shh=c/(b+shh);
                }
                fg[0]=(sk+sm)*(sk-sm)-shh;
                fg[1]=sk*ek;
                for(i=kk; i<=mm-1; i++)
                { 
					_MSV_2(fg,cs);
                    if(i!=kk) e[i-2]=fg[0];
                    fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                    e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                    fg[1]=cs[1]*s[i];
                    s[i]=cs[0]*s[i];
                    if((cs[0]!=1.0)||(cs[1]!=0.0))
                      for(j=1; j<=n; j++)
                      {
                          d=cs[0]*v(j-1,i-1)+cs[1]*v(j-1,i);
                          v(j-1,i)=-cs[1]*v(j-1,i-1)+cs[0]*v(j-1,i);
                          v(j-1,i-1)=d;
                      }
                    _MSV_2(fg,cs);
                    s[i-1]=fg[0];
                    fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
                    s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
                    fg[1]=cs[1]*e[i];
                    e[i]=cs[0]*e[i];
                    if(i<m)
                      if((cs[0]!=1.0)||(cs[1]!=0.0))
                        for(j=1; j<=m; j++)
                        { 
                            d=cs[0]*u(j-1,i-1)+cs[1]*u(j-1,i);
                            u(j-1,i)=-cs[1]*u(j-1,i-1)+cs[0]*u(j-1,i);
                            u(j-1,i-1)=d;
                        }
                }
                e[mm-2]=fg[0];
                it=it-1;
            }
            else
            { 
				if(ks==mm)
                {
					kk=kk+1;
                    fg[1]=e[mm-2]; 
					e[mm-2]=0.0;
                    for(ll=kk; ll<=mm-1; ll++)
                    {
						i=mm+kk-ll-1;
                        fg[0]=s[i-1];
                        _MSV_2(fg,cs);
                        s[i-1]=fg[0];
                        if(i!=kk)
                        {
							fg[1]=-cs[1]*e[i-2];
                            e[i-2]=cs[0]*e[i-2];
                        }
                        if((cs[0]!=1.0)||(cs[1]!=0.0))
                          for(j=1; j<=n; j++)
                          { 
                              d=cs[0]*v(j-1,i-1)+cs[1]*v(j-1,mm-1);
                              v(j-1,mm-1)=-cs[1]*v(j-1,i-1)+cs[0]*v(j-1,mm-1);
                              v(j-1,i-1)=d;
                          }
                    }
                }
                else
                { 
					kk=ks+1;
                    fg[1]=e[kk-2];
                    e[kk-2]=0.0;
                    for(i=kk; i<=mm; i++)
                    {
						fg[0]=s[i-1];
                        _MSV_2(fg,cs);
                        s[i-1]=fg[0];
                        fg[1]=-cs[1]*e[i-1];
                        e[i-1]=cs[0]*e[i-1];
                        if((cs[0]!=1.0)||(cs[1]!=0.0))
                          for(j=1; j<=m; j++)
                          {
                              d=cs[0]*u(j-1,i-1)+cs[1]*u(j-1,kk-2);
                              u(j-1,kk-2)=-cs[1]*u(j-1,i-1)+cs[0]*u(j-1,kk-2);
                              u(j-1,i-1)=d;
                          }
                    }
                }
            }
        }
    }
    return(1);
}

//一般实矩阵的奇异值分解辅助函数
template <class _Ty>
void _MSV_1(matrix<_Ty>& a, valarray<_Ty>& e, valarray<_Ty>& s, 
											matrix<_Ty>& v, int m, int n)
{
	int i,j;
    _Ty d;
    if(m>=n) i=n;
    else i=m;
    for(j=1; j<i; j++)
    {
		a(j-1,j-1)=s[j-1];
        a(j-1,j)=e[j-1];
    }
    a(i-1,i-1)=s[i-1];
    if(m<n) a(i-1,i)=e[i-1];
    for(i=1; i<n; i++)
    for(j=i+1; j<=n; j++)
    { 
        d=v(i-1,j-1);
		v(i-1,j-1)=v(j-1,i-1); 
		v(j-1,i-1)=d;
    }
}

//一般实矩阵的奇异值分解辅助函数
template <class _Ty>
void _MSV_2(valarray<_Ty>& fg, valarray<_Ty>& cs)
{
	_Ty r,d;
	r = Abs(fg[0])+Abs(fg[1]);
    if(FloatEqual(r,0.0))
    {
		cs[0]=1.0;
		cs[1]=0.0;
		d=0.0;
	}
    else 
    {
		d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
        if(Abs(fg[0])>Abs(fg[1]))
        {
			d=Abs(d);
            if(fg[0]<0.0) d=-d;
        }
        if(Abs(fg[1])>Abs(fg[0])||FloatEqual(Abs(fg[1]),Abs(fg[0])))
        { 
			d=Abs(d);
            if(fg[1]<0.0) d=-d;
        }
        cs[0]=fg[0]/d; 
		cs[1]=fg[1]/d;
    }
    r=1.0;
    if(Abs(fg[0])>Abs(fg[1])) r=cs[1];
    else
      if(FloatNotEqual(cs[0],0.0)) r=1.0/cs[0];
    fg[0]=d;
	fg[1]=r;
}

//广义逆的奇异值分解
template <class _Ty>
int GeneralizedInversionSingularValue(matrix<_Ty>& a, matrix<_Ty>& aa, 
							_Ty eps, matrix<_Ty>& u, matrix<_Ty>& v)
{
	int k(0), l, ka;
	int m = a.GetRowNum();
	int n = a.GetColNum();

	if(m > n) ka = m + 1;
	else ka = n + 1;
  
	int i=MatrixSingularValue(a,u,v,eps);
    if (i<0) return(0);
    int j = n;
    if(m < n) j = m;
    j = j - 1;
    while((k<=j)&&(FloatNotEqual(a(k,k),0.0))) k = k + 1;
    k = k - 1;
    for(i=0; i<n; i++)
    for(j=0; j<m; j++)
    {
		aa(i,j)=0.0;
        for(l=0; l<=k; l++)
        {
            aa(i,j)=aa(i,j)+v(l,i)*u(j,l)/a(l,l);
        }
    }
    return(1);
}

#endif			// _MATRIX_INL

