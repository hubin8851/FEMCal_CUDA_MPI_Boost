//LinearEquation.inl		线性方程(组)求解函数(方法)定义
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31

#ifndef _LINEAREQUATION_INL
#define _LINEAREQUATION_INL

//全选主元高斯消去法
template <class _Ty>
int LE_TotalChoiceGauss(matrix<_Ty>& a, valarray<_Ty>& b)  
{ 
	long double MaxValue, tmp;		//记录主元绝对值
	int l(1), i, j, is;
    bool yn;
    
	int n = a.GetColNum();			//方程组阶数
	
	valarray<int> js(n);			//保存换列位置
    
	for(int k = 0; k < n - 1; k++)	//全选主元
	{	
		MaxValue = 0.0;				//给保存主元绝对值变量赋初值
		        
		for(i = k; i < n; i++)
			for(j = k; j < n; j++)
            {		
				tmp = Abs(a(i, j));	//求m(i,j)绝对值
				if(tmp > MaxValue)	//发现一个更大的主元
				{ 
					MaxValue = tmp;	//保存新主元绝对值
					js[k] = j;		//新主元所在列
					is = i;			//新主元所在行
				}
            }
			
		yn = FloatEqual(MaxValue, 0);
        if(yn) l = 0;				//主元为0
		else
		{
			if(js[k] != k)			//换列
				for(i = 0; i < n; i++) swap(a(i, k), a(i, js[k]));
								
			if(is != k)				//换行
			{ 
				for (j = k; j < n; j++)	swap(a(k, j), a(is, j));

				swap(b[k], b[is]);	//方程组右边第k元素与第is元素交换
			}
		}
        
		if(l == 0)					//矩阵奇异(主元为0)
		{
			printf("fail 1\n");
            return 0;				// 求解失败，返回0值
		}
        
		MaxValue =  Abs(a(k, k));

        for(j = k + 1; j < n; j++)	a(k, j) /= a(k, k); //MaxValue;
        
		b[k] /= a(k, k); //MaxValue;
        for(i = k + 1; i < n; i++)
		{
			for(j = k + 1; j < n; j++)
			{
                a(i, j) = a(i, j) - a(i, k) * a(k, j);
			}
            
			b[i] = b[i] - a(i, k) * b[k];
		}
	}
    
	MaxValue = Abs(a((n - 1), (n - 1)));	//主元

	yn = FloatEqual(MaxValue, 0);
    if(yn)							//主元为0
	{
		cout<<"fail 2"<<endl;
        return(0);					//求解失败，返回0值
	}

	b[n - 1] /= a((n - 1), (n - 1));//求解方程组右边X的解

    for(i = n - 2; i >= 0; i--)		//回代过程
	{
		_Ty t = 0.0;
        
		for(j = i + 1; j < n; j++)	t = t + a(i, j) * b[j];
        
		b[i] = b[i] - t;
	}
    
	js[n - 1] = n - 1;				//X最后一个元素不用换
    for(k = n - 2; k >= 0; k --)	//k可以从n-2开始
		if(js[k] != k)				//交换X的元素位置(由全选换列产生的)
			swap(b[k], b[js[k]]);
    
	return(1);						//方程组求解成功！
}

//全选主元高斯-约当消去法
template <class _Ty>
int LE_TotalChoiceGaussJordan(matrix<_Ty> & a, matrix<_Ty> & b)
{ 
	long double MaxValue, tmp;	//主元绝对值
	int l(1), k, i, j, is;
	bool yn;
	
    int n = a.GetColNum();		//方程组阶数
	int m = b.GetColNum();		//方程组右端常数向量的个数

	valarray<int> js(n);  
	for(k = 0; k < n; k++)
	{
		MaxValue = long double(0.0);
        
		for(i = k; i < n; i++)
			for(j = k; j < n; j++)
            {
				tmp = Abs( a(i, j) );
				if(tmp > MaxValue)
				{ 
					MaxValue = tmp;
					js[k] = j;
					is = i;
				}
            }
		
		yn = FloatEqual(MaxValue, 0);
		if(yn) l = 0; 		
		else
		{
			if(js[k] != k)
				for(i = 0; i < n; i++)
					swap(a(i, k), a(i, js[k]));
				
			if(is != k)
			{
				for(j = k; j < n; j++)
					swap(a(k, j), a(is, j));
				
				for(j = 0; j < m; j++)
					swap(b(k, j), b(is, j));

			}
		}

		if(l == 0)
		{
			cout<<"fail"<<endl;
			return 0;
		}
		
		for(j = k + 1; j < n; j++)
			a(k, j) /= a(k, k);

		for(j = 0; j < m; j++)
			b(k, j) /= a(k, k);				

		for(j = k + 1; j < n; j++)
			for(i = 0; i < n; i++)
				if(i != k)
					a(i, j) -= a(i, k) * a(k, j);

		for(j = 0; j < m; j++)
			for(i = 0; i< n; i++)
				if(i != k)
					b(i, j) -= a(i, k) * b(k, j);					
	}

    for(k = n - 1; k >= 0; k --)
	{
		if(js[k] != k)
			for(j = 0 ; j < m;  j++)
				swap(b(k, j), b(js[k], j));

	}
	return 1;				//成功返回
}

//求解三对角线方程组的追赶法
template <class _Ty>
int LE_TridiagonalEquationGauss(valarray<_Ty>& b, valarray<_Ty>& d)
{
	int j, k;
    _Ty s;
	bool yn;

    int n = d.size();	//方程组阶数
	//m为n阶三对角矩阵三条对角线上的元素个数，也是一维数组b的长度。
	//它的值应为m=3n-2, 在本函数中要对此值作检查
	int m = b.size();	

    if(m != (3 * n - 2))	//验证m是否等于3*n-2
	{
		cout << "Error!" << endl;
		return(-2);
	}
    
	for(k = 0; k < n - 1; k++)
	{
		j = 3 * k;
		s = b[j];

		yn = FloatEqual(s, 0.0);
		if(yn)						//主元为0
		{
			cout << "Fail!" << endl;
			return(0);
		}

        b[j + 1] = b[j + 1] / s;
        d[k] = d[k] / s;
        b[j + 3] = b[j + 3] - b[j + 2] * b[j + 1];
        d[k + 1] = d[k + 1] - b[j + 2] * d[k];
	}
    
	s = b[3 * n - 3];
    
	yn = FloatEqual(s, 0.0);
	if(yn)						//主元为0
	{
		cout << "Fail!" << endl;
		return(0);
	}
    
	d[n - 1] = d[n - 1] / s;
    
	for(k = n - 2; k >= 0; k--)
		d[k] = d[k] - b[3 * k + 1] * d[k + 1];
    
	return (2);	//运算成功，正常返回
}

//一般带型方程组的求解
template <class _Ty>
int LE_StrapEquationGauss(matrix<_Ty>& b, matrix<_Ty>& d, int l, int il)
{
	int ls, k, i, j, is, u, v, x, y, hh, xx, yy, xxx, yyy;
    _Ty p, t, tt;
	bool yn;
    
	int n = b.GetRowNum();		//矩阵b的行数
    int nn = b.GetColNum();
	int mm = d.GetColNum();		//方程组右边常数向量的个数
	
    if(il != (2 * l + 1))		//参数中半带宽l与带宽il的关系不对
	{
		cout<<"fail 1"<<endl;
		return(-2);
	}
    
	ls = l;
    
	for(k = 0; k < n - 1; k++)
	{
		p = 0.0;

        for(i = k ; i <= ls; i++)
		{
			u = i * il;
			y = u % nn;  x = (u - y) / nn;
			t = Abs(b(x, y));
            if(t > p) 
			{ 
				p = t;
				is = i;
			}
		}
        
		if( FloatEqual(p, 0.0) )		
		{
			cout << "fail 2" << endl;
			return(0);
		}
        
		for(j = 0; j < mm; j++)
		{
			u = k * mm + j;
			v = is * mm + j;
			y = u % mm;
			x = (u - y) / mm;
			yy = v % mm;
			xx = (v - yy) / mm;
    		tt = d(x, y);
			d(x, y) = d(xx, yy);
			d(xx, yy) = tt;
        }
        
		for(j = 0; j < il; j++)
		{
			u = k * il + j; 
			v = is * il + j;
            y = u % nn;  
			x = (u - y) / nn;
			yy = v % nn; 
			xx = (v - yy) / nn;
    		tt = b(x, y);
			b(x, y) = b(xx, yy); 
			b(xx, yy) = tt;
		}

        for(j = 0; j < mm; j++)
		{
			u = k * mm + j; v = k * il;
			y = u % mm;  x = (u - y) / mm;
            yy = v % nn;  xx = (v - yy) / nn;
	
			d(x, y) = d(x, y) / b(xx, yy);
		}

        for(j = 1; j < il; j++)
		{
			u = k * il + j;
			v = k * il;
			y = u % nn;  x = (u - y) / nn;
			yy = v % nn;  xx = (v - yy) / nn;
			b(x, y) = b(x, y) / b(xx, yy);
		}
        
		for(i = k + 1; i <= ls; i++)
		{
            hh = i * il;
		    y = hh % nn;
			x = (hh - y) / nn;
			t = b(x, y);
            for(j = 0; j < mm; j++)
			{
				u = i * mm + j;
				v = k * mm + j;
				y = u % mm;
				x = (u - y) / mm;
				yy = v % mm;
				xx = (v - yy) / mm;    
                d(x, y) = d(x, y) - t * d(xx, yy);
			}
            
			for(j = 1; j < il; j++)
			{
				u=i * il  +j;
				v = k * il + j;
				y = u % nn;  
				x = (u - y) / nn;
				yy = v % nn; 
				xx = (v - yy) / nn;
                yyy = (u - 1) % nn;  xxx = (u - yyy - 1) / nn;
                b(xxx, yyy) = b(x, y) - t * b(xx, yy);
			}
            
			u = i * il + il - 1; 
			y = u % nn;  x = (u - y) / nn;
			
			b(x, y) = 0.0;
		}
        
		if(ls != (n - 1))	ls++;
	}
    
	u = (n - 1) * il; 
	y = u % nn;  
	x = (u - y) / nn;
	p = b(x, y);
    
	yn = FloatEqual(p, 0.0);
	if(yn)
	{
		cout<<"fail 3"<<endl;
		return(0);
	}
    
	for(j = 0; j < mm; j++)
	{
		u = (n - 1) * mm + j; 
		y = u % mm;  
		x = (u - y) / mm;
		d(x, y) = d(x, y) / p;
	}

    ls = 1;
    for(i = n  - 2; i >= 0; i--)
	{
		for(k = 0; k < mm; k++)
		{ 
			u = i * mm + k;
            for(j = 1; j <= ls; j++)
			{
				v = i * il +j;
				is = (i + j) * mm + k;
                y = u % mm;  
				x = (u - y) / mm;
				yy = v % nn;  
				xx = (v - yy) / nn;
                yyy = is % mm;  
				xxx = (is - yyy) / mm;
                d(x, y) = d(x, y) - b(xx, yy) * d(xxx, yyy);
			}
		}
        if(ls != (il - 1)) ls++;
	}
    
	return 2;	//运算成功，正常返回
}

//求解对称方程组的分解法
template <class _Ty>
int LE_SymmetryEquation(matrix<_Ty>& a, matrix<_Ty>& c)
{
	int i, j, k, k1, k2, k3;
    _Ty p;
	bool yn;

	if(MatrixSymmetry(a)!=true)
		return (-1);			//方程组不对称
    
	int n = a.GetColNum();		//方程组阶数
	int m = c.GetColNum();		//方程组右边常数向量的个数
    
	yn = FloatEqual(a(0,0), 0.0);
	if(yn)						//主元为0
	{
		cout << "fail" << endl;
		return (-2);
	}

    for(i = 1; i < n; i++)
		a(i, 0) /= a(0, 0);

    for(i = 1; i < n - 1; i++)
	{
        for(j = 1; j <= i; j++)
            a(i, i) -= a(i, (j - 1)) * a(i, (j - 1)) * a((j - 1), (j - 1));
        
		p = a(i,i);
        yn = FloatEqual(p, 0.0);
		if(yn)						//主元为0
		{
			cout<<"fail"<<endl;
			return (-2);
		}
        
		for(k = i + 1; k < n; k++)
		{
            for(j = 1; j <= i; j++)
				a(k, i) -= a(k, (j - 1)) * a(i,(j - 1)) * a((j - 1), (j - 1));
            a(k, i) /= p;
		}
	}
    
	
    for(j = 1; j < n; j++)
	    a((n-1), (n-1)) -= a((n - 1), (j - 1)) * a((n - 1), (j - 1)) * a((j - 1), (j - 1));
    
	p=a((n - 1), (n - 1));
    yn = FloatEqual(p, 0.0);
	if(yn)						//主元为0
	{
		cout<<"fail"<<endl;
		return (-2);
	}
    
	for(j = 0; j < m; j++)
		for(i = 1; i < n; i++)
			for(k = 1; k <= i; k++)
				c(i, j) -= a(i, (k - 1)) * c((k - 1),j);
		
	for(i = 1; i < n; i++)
		for(j = i; j < n; j++)
            a((i - 1), j) = a((i - 1), (i - 1)) * a(j, (i - 1));
    
	for(j = 0; j < m; j++)  
    {
		c((n - 1), j) /= p;
		for(k = 1; k < n; k++)
		{
			k1 = n - k; 
			k3 = k1 - 1; 
            for(k2 = k1; k2 < n; k2++)
                c(k3, j) -= a(k3, k2) * c(k2, j);
			c(k3, j) /= a(k3, k3);
		}
	}
    
	return (1);	//运算成功，正常返回
}

//求解对称正定方程组的平方根法
template <class _Ty>
int LE_SymmetryRegularEuationSquareRoot(matrix<_Ty> & a, matrix<_Ty> & d)
{
	int i, j, k;

	if(MatrixSymmetryRegular(a,1) != 2)	//判别矩阵a是否对称正定
		return(-1);						//矩阵a不对称正定
    
	int n = a.GetColNum();				//方程组阶数
	int m = d.GetColNum();				//方程组右边常数向量的个数
	
	a(0,0) = sqrt( a(0, 0) );
	
    for(j = 1; j < n; j++) a(0, j) /= a(0, 0);
    
	for(i = 1; i < n; i++)
	{
        for(j = 1; j <= i; j++)
            a(i, i) = a(i, i) - a((j - 1), i) * a((j - 1), i);
    
		a(i, i) = sqrt(a(i, i));
        if(i != (n - 1))
		{
			for(j = i + 1; j < n; j++)
			{
                for(k = 1; k <= i; k++)
					a(i, j) -= a((k - 1), i) * a((k - 1), j);
                
				a(i, j) /= a(i, i);
			}
		}
	}
    
	for(j = 0; j < m; j++)
	{
		d(0, j) = d(0, j) / a(0, 0);
        for(i = 1; i < n; i++)
		{
            for(k = 1; k <= i; k++)
				d(i, j) -= a((k - 1), i) * d((k - 1), j);
            
			d(i, j) /= a(i, i);
		}
	}
    
	for(j = 0; j < m; j++)
	{
        d((n - 1), j) = d((n - 1), j) / a((n - 1), (n - 1));
		
		for(k = n - 1; k >= 1; k --)
		{
		    for(i = k; i < n; i++)
		        d((k - 1), j) -= a((k - 1), i) * d(i, j);
        
            d((k - 1), j) /= a((k - 1), (k - 1));			
		}
	}
    
	return (1);	//运行成功，正常返回
}

//求解大型稀疏方程组的全选主元高斯-约当消去法
template <class _Ty>
int LE_SparseEuationTotalChoiceGaussJordan(matrix<_Ty> & a, valarray<_Ty> & b)
{ 
	int i, j, k, is;
    _Ty d, t;

    int n = a.GetColNum();		//方程组阶数
    valarray<int> js(n);
	
    for(k = 0; k < n; k++)
	{ 
		d = 0.0;
    	for(i = k; i < n; i++)
			for(j = k; j < n; j++)
			{
				t = Abs(a(i, j));
		
				if(t > d)
				{ 
					d = t;
					js[k] = j; 
					is = i;
				}
			}
		
		if( FloatEqual(d, 0.0) )	//主元为0
		{					
			cout<<"fail 1"<<endl;
			return (0);
		}
		
		if(is != k)
		{ 
			for(j = k; j < n; j++)
			{
				t = a(k, j);
				a(k, j) = a(is, j);
				a(is, j) = t;
			}
			
			t = b[k];
			b[k] = b[is];
			b[is] = t;
		}
		
		if(js[k] != k)
			for(i = 0; i < n; i++)
			{
				t = a(i, k);
				a(i, k) = a(i, js[k]); 
				a(i, js[k]) = t;
			}
		
		t = a(k, k);
		
		for(j = k + 1; j < n; j++)
		{
			if(a(k, j) != 0.0)
				a(k, j) = a(k, j) / t;
		}
		
		b[k] = b[k] / t;
		for(j = k + 1; j < n; j++)
		{ 
			if(a(k, j) != 0.0)
			{
				for(i = 0; i < n; i++)
				{
					if((i != k) && (a(i, k) != 0.0))
					{ 
						a(i, j) = a(i, j) - a(i, k) * a(k, j);
					}
				}
			}
		}
		
		for(i = 0; i <  n; i++)
		{
			if((i != k) && (a(i, k) != 0.0))
				b[i] = b[i] - a(i, k) * b[k];
		}
	}
	
	for(k = n - 1; k >= 0; k--)
		if(k != js[k])
		{
			t = b[k];
			b[k] = b[js[k]]; 
			b[js[k]] = t;
		}
	   
	return (1);	//运行成功，正常返回
}

//求解托伯利兹方程组的列文逊方法
template <class _Ty>
int LE_ToeplitzEuationLevinson(valarray<_Ty>& t,  valarray<_Ty>& b, valarray<_Ty>& x)
{
	int i, j, k;
    _Ty a, beta, q, c, h;
	bool yn;

    int n = t.size();			//方程组阶数
	valarray<_Ty> y(n);
	valarray<_Ty> s(n);

	a = t[0];
	yn = FloatEqual(a, 0.0);
    if(yn)
	{
		cout << "fail 1" << endl;
		return(-1);
	}
    
	y[0] = 1.0;
	x[0] = b[0] / a;

    for(k = 1; k < n; k++)
	{
		beta = 0.0;
		q = 0.0;
        
		for(j = 0; j < k; j++)
		{
			beta = beta + y[j] * t[j + 1];
            q = q + x[j] * t[k - j];
		}
        
		yn = FloatEqual(a, 0.0);
		if(yn)		
		{
			cout << "fail 2" << endl;
			return(-1);
		}
        
		c = -beta / a; 
		s[0] = c * y[k - 1]; 
		y[k] = y[k - 1];
        
		if(k != 1)
			for(i = 1; i < k; i++)
				s[i] = y[i - 1] + c * y[k - i - 1];
			
		a = a + c * beta;
        
		yn = FloatEqual(a, 0.0);
		if(yn)
		{
			cout << "fail 3" << endl;
			return (-1);
		}

        h = (b[k] - q) / a;
        
		for(i = 0; i < k ; i++)
		{
			x[i] = x[i] + h * s[i];
			y[i] = s[i];
		}
        
		x[k] = h * y[k];
	}
  
	return (1);	//运行成功，正常返回
}

//高斯-赛德尔迭代
template <class _Ty>
int LE_GaussSeidelIteration(matrix<_Ty>& a, valarray<_Ty>& b, 
											valarray<_Ty>& x, _Ty eps)
{
	int i, j;
    _Ty p, t, s, q;

	int n = a.GetColNum();		//方程组阶数
    
	for(i = 0; i < n; i++)
	{	
		p = 0.0;
		x[i] = 0.0;
        
		for(j = 0; j < n; j++)
			if(i != j)
				p = p + Abs(a(i, j));
        
		if(p >= Abs(a(i, i)) )	//主对角线不占绝对优势
		{
			cout<<"fail 1"<<endl;
			return(-1);
		}
	}
    
	p = eps + 1.0;
    while(p > eps  || FloatEqual(p, eps))
	{
		p = 0.0;
        for(i = 0; i < n; i++)
		{
			t = x[i];
			s = 0.0;
            for(j = 0; j < n; j++)
				if(j != i) 
					s = s + a(i, j) * x[j];
			
			x[i] = (b[i] - s) / a(i, i);
            q = Abs(x[i] - t) / (1.0 + Abs(x[i]));
            if(q > p) p = q;
		}
	}

	return (1);	//运行成功，正常返回
}
    
//求解对称正定方程组的共轭梯度法
template <class _Ty>
int LE_SymmetryRegularEuationConjugateGradient(matrix<_Ty> & a, 
						valarray<_Ty> & b, _Ty eps, valarray<_Ty> & x)
{
	int i, k;
    _Ty  alpha, beta, d, e;
	
	if(MatrixSymmetryRegular(a,1) != 2)	//判别矩阵a是否对称正定
		return(-1);						//矩阵a不对称正定
    
	int n = a.GetColNum();				//方程组阶数
    valarray<_Ty> r(n);
    
	matrix<double> q(n, 1);
    matrix<double> p(n, 1);
	matrix<double> s(n, 1);
    matrix<double> xx(n, 1);
    
	for(i = 0; i < n; i++) 
		xx(i, 0) = x[i];

    for(i = 0; i < n; i++)
	{
		xx(i,0) = 0.0;
		p(i, 0) = b[i];
		r[i] = b[i]; 
	}
    
	i = 0;

    while(i < n)
	{		
		MatrixMultiply(s, a, p);
        d = 0.0;
		e = 0.0;
        
		for(k = 0; k < n; k++)
		{
			d = d + p(k, 0) * b[k];
			e = e + p(k, 0) * s(k,0); 
		}
        
		alpha = d / e;
        
		for(k = 0; k < n; k++)
			xx(k,0) = xx(k, 0) + alpha * p(k, 0);
        		
		MatrixMultiply(q, a, xx);
        
		d = 0.0;
        
		for(k = 0; k < n; k++)
		{
			r[k] = b[k] - q(k, 0);
			d = d + r[k] * s(k, 0);
		}
        
		beta = d / e;
		d = 0.0;
        
		for(k = 0; k < n; k++) d = d + r[k] * r[k];
        
		d = sqrt(d);
        
		if(d < eps) 
		{
			for(i = 0; i < n; i++) 
				x[i] = xx(i, 0);

		}

        for(k = 0; k < n; k++) 
 			p(k, 0) = r[k] - beta * p(k, 0);
        
		i++;
	}

	for(i = 0; i < n; i++) 
		x[i] = xx(i, 0);

	return (1);	//运行成功，正常返回
}

//求解线性最小二乘问题的豪斯荷尔德变换法
template <class _Ty>
int LE_LinearLeastSquareHouseholder(matrix<_Ty>& a, valarray<_Ty>& b, matrix<_Ty>& q)
{
	int i,j;
    _Ty d;

    int n = a.GetColNum();		//系数矩阵a的行数，m>=n
	int m = a.GetRowNum();		//系数矩阵a的列数，n<=m
    
	valarray<_Ty> c(n);

    i = MatrixQR(a, q);			//一般m×n阶的实矩阵QR分解函数

    if(i == 0)
	{
		return (-1);
	}

	for(i = 0; i < n; i++)
	{
		d = 0.0;

        for(j = 0; j < m; j ++)
				d = d + q(j ,i) * b[j];
        
		c[i] = d;
	}
    
	b[n - 1] = c[n - 1] / a((n - 1), (n - 1));
	
    for(i = n - 2; i >= 0; i--)
	{
		d = 0.0;

        for(j = i + 1; j < n; j++)
				d = d + a(i,j) * b[j];
        
		b[i] = (c[i] - d) / a(i, i);
	}
	return (1);	//运行成功，正常返回
}

//求解线性最小二乘问题的广义逆法
template <class _Ty>
int LE_LinearLeastSquareGeneralizedInverse(matrix<_Ty>& a, 
				valarray<_Ty>& b, valarray<_Ty>& x, matrix<_Ty>& aa, 
								_Ty eps, matrix<_Ty>& u, matrix<_Ty>& v)
{
	int i, j, ii;
       
	int n = a.GetColNum();	//系数矩阵a的列数
	int m = a.GetRowNum();	//系数矩阵a的行数
    
	int ka = (n > m) ? n : m;	
	
	ka++;

	ii = GeneralizedInversionSingularValue(a, aa, eps, u, v);
    
	if(ii < 0)  return(-1);
    
	for(i = 0; i < n; i ++)
	{
		x[i] = 0.0;
        for(j = 0; j < m; j ++)
		{		
			x[i] += aa(i, j) * b[j];			
		}
	}
    
	return (1);	//运行成功，正常返回
}

//病态方程组的求解
template <class _Ty>
int LE_IllConditionedEquation(matrix<_Ty> & a, valarray<_Ty> & b, _Ty eps, valarray<_Ty> & x)
{
	int i(60), j, kk;
    _Ty q, qq;

    int n = a.GetColNum();			//方程组阶数

	matrix<_Ty> p(n, n);
	matrix<_Ty> ee(n, 1);
	matrix<_Ty> xx(n, 1);
	valarray<double> rr(n);
	    
	for(int k = 0; k < n; k++)
		for(j = 0; j < n; j++)
			p(k, j) = a(k, j);
	
	for(k = 0; k < n; k++) x[k] = b[k];
    kk = LE_TotalChoiceGauss(p, x);
	
	for(int w = 0; w < n; w++) xx(w, 0) = x[w];  
	
	if(kk == 0)
	{ 
		for( w = 0; w < n; w++) x[w] = xx(w, 0);
	    
		return 0;
	}
    
	q = 1.0 + eps;
    while(q > eps || FloatEqual(q, eps))
	{
		if(i == 0)
		{ 
			for( w = 0; w < n; w++) x[w] = xx(w, 0);
			
			return 0;
		}
        i--;        
		
		MatrixMultiply(ee, a, xx);
		for( w = 0; w < n; w++) x[w] = xx(w, 0);        
			
		for(k = 0; k < n; k++) rr[k] = b[k] - ee(k, 0);

        for(k = 0; k < n; k++)
			for(j = 0; j < n; j++) p(k, j) = a(k, j);
			
		kk = LE_TotalChoiceGauss(p, rr);	

        if(kk == 0)
		{ 
			for( w = 0; w < n; w++) x[w] = xx(w, 0);
	
			return 0;
		}
        
		q = 0.0;
        for(k = 0; k < n; k++)
		{
			qq = Abs(rr[k]) / (1.0 + Abs(x[k] + rr[k]));
            if(qq > q) q = qq;
		}
        
		for(k = 0; k < n; k++) xx(k, 0) = xx(k, 0) + rr[k];
		for(int w = 0; w < n; w++)  x[w] = xx(w, 0);	    
	}
		
	for( w = 0; w < n; w++) x[w] = xx(w, 0);
	
	return (1);	//运行成功，正常返回
}

#endif //_LINEAREQUATION_INL

