//NonLinearEquation.inl 	非线性方程(组)求解函数(方法)定义
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _NONLINEAREQUATION_INL
#define _NONLINEAREQUATION_INL

//用二分法搜索方程f(x)=0在区间[a,b]内的全部实根
template <class _Ty>
inline size_t 
RootHalves(_Ty a, _Ty b, _Ty step, _Ty eps, valarray<_Ty>& x, size_t m)
{
    int n(0), js;
    _Ty z=a, z1, y, y1, z0, y0;
   
	y = FunctionValueRH(z);
    while((z<=b+step/2.0)&&(n!=m))
	{ 
	  if(Abs(y) < eps)
      {
		n = n + 1;
		x[n-1] = z;
		z = z + step / 2.0;
		y=FunctionValueRH(z);
      }
      else
      {
		z1=z+step;
		y1=FunctionValueRH(z1);
        if(Abs(y1)<eps)
        {
			n=n+1;
			x[n-1]=z1;
            z=z1+step/2.0;
			y=FunctionValueRH(z);
        }
        else
			if(y*y1>0.0)
            {
				y=y1;
				z=z1;
			}
            else
            {
				js=0;
                while(js==0)
                {
					if(Abs(z1-z)<eps)
                    {
						n=n+1;
						x[n-1]=(z1+z)/2.0;
                        z=z1+step/2.0;
						y=FunctionValueRH(z);
                        js=1;
                    }
                    else
                    {
						z0=(z1+z)/2.0;
						y0=FunctionValueRH(z0);
                        if(Abs(y0)<eps)
                        {
							x[n]=z0; n=n+1; js=1;
                            z=z0+step/2.0;
							y=FunctionValueRH(z);
                        }
                        else
							if((y*y0) < 0.0)
							{
								z1 = z0;
								y1 = y0;
							}
							else
							{
								z = z0;
								y = y0;
							}
                      }
                  }
              }
          }
	}
    return(n);	//返回搜索到的实根个数
}

//牛顿(Newton)法求解非线性方程一个实根
template <class _Ty>
inline int 
RootNewton( _Ty& x, _Ty eps, size_t js)
{ 
    int k = 1, r = js;
    _Ty  x1;
	_Ty x0 = x;
	valarray<_Ty> y(2);

    FunctionValueRN(x0, y);
    _Ty d=eps + 1.0;
    while((d > eps || FloatEqual(d, eps))&&(r != 0))
    {
		if(FloatEqual(y[1],0))	return(-1);
        x1 = x0 - y[0] / y[1];
        FunctionValueRN(x1,y);
        d=Abs(x1 - x0);
		_Ty p=Abs(y[0]);
        if(p > d) d = p;
        x0 = x1;
		r -= 1;
    }
    x = x1;
    k = js - r;
    return(k);
}

//埃特金(Aitken)法求解非线性方程一个实根
template <class _Ty>
inline int 
RootAitken( _Ty& x, _Ty eps, size_t js)
{
    _Ty u, v, x0;
    int r(0), flag(0);

	x0 = x;

    while((flag == 0) && (r != js))
    {	
		r++; 
        u = FunctionValueRA(x0);
		v = FunctionValueRA(u);
        if(Abs(u-v)<eps)
		{
			x0 = v;
			flag = 1;
		}
        else x0=v-(v-u)*(v-u)/(v-2.0*u+x0);
    }
    x = x0;
	r = js - r;
    return(r);
}

//连分式(Fraction)法求解非线性方程一个实根
template <class _Ty>
inline int 
RootFraction(_Ty& x, _Ty eps)
{
	int r(10), i, j, m, it;
    _Ty x0, q(1.0e+35), h(0), z;
	valarray<_Ty> a(10), y(10);
	
	x0 = x; 

    while(r!=0)
    { 
		r = r - 1; 
		j = 0; 
		it = r;
        while(j<=7)
        {
			if(j<=2)
				z = x0 + 0.1 * j;
            else
				z = h;
            y[j] = FunctionValueRF(z);
            h = z;
            if(j==0)
				a[0] = z;
            else
            {
				m = 0;
				i = 0;
                while((m==0)&&(i<j))
                {
					if(FloatEqual((h-a[i]),0))
						m = 1;
                    else
						h=(y[j]-y[i])/(h-a[i]);
                    i++;
                }
                a[j] = h;
                if(m!=0) a[j]=q;
                h=0.0;
                for(i=j-1; i>=0; i--)
                {
					if(FloatEqual((a[i+1]+h),0))
						h = q;
                    else
						h=-y[i]/(a[i+1]+h);
                }
                h = h + a[0];
            }
            if(Abs(y[j])>=eps) j++;
            else 
			{
				j = 10;
				r = 0;
			}
        }
        x0 = h;
    }
    x = h;
    return(10-it);
}

//QR法求代数方程全部根
template <class _Ty, class _Tz>
inline int 
RootQR(valarray<_Ty>& a, valarray<_Tz>& x, _Ty eps, size_t jt)
{
	int i, j, n;

    n = x.size();			//多项式方程的次数

	matrix<_Ty> q(n, n);

    for(j=0; j<n; j++)
		q(0,j) = -a[n-j-1] / a[n];

    for(j=1; j<n; j++)
		for(i=0; i<n; i++)
			q(j,i) = 0;

	for(i=0; i<n-1; i++)
      q(i+1, i) = 1;
	
	i = EigenvalueVectorHessenbergQR(q,x,eps,jt);	//求全部特征值的QR方法
    
	return(i);
}

//牛顿下山(NewtonHillDown)法求解实系数代数方程全部根(实根和复根)
template <class _Ty, class _Tz>
inline int 
RootNewtonHillDown(valarray<_Ty>& a, valarray<_Tz>& cx)
{
	int m, i, jt, k, is, it, n;
    _Ty t, x, y, x1, y1, dx, dy, dd, dc, c, g, g1;
	_Ty w, u, p, q, pq, v, u1, v1;
  
    n = cx.size();			//多项式方程的次数
    m = n;

    while((m > 0) && (FloatEqual(a[m], 0.0))) m--;

    if(m <= 0) return(-1);
    
	for(i = 0; i <= m; i++) a[i] /= a[m];

    for(i = 0; i <= m / 2; i++)
    {
		w = a[i];
		a[i] = a[m - i];
		a[m - i] = w;
	}

    k = m;
	is = 0;
	w = 1.0;
    jt = 1;
    while(jt == 1)
    {
		pq = Abs(a[k]);
		while(pq < LONGDOUBLEERROR)
        {
			cx[k - 1] = complex<_Ty>(0.0, 0.0);
			k = k - 1;
            if(k == 1)
            {
				cx[0] = complex<_Ty>(-a[1] * w / a[0], 0.0);
                return(1);
            }
            pq = Abs(a[k]);
        }
		q = log(pq);
		q = q / (1.0 * k);
		q = exp(q);
        p = q;
		w = w * p;
        for(i = 1; i <= k; i++)
        {
			a[i] /= q;
			q *= p;
		}
        x = 0.0001;
		x1 = x; 
		y = 0.2;
		y1 = y; 
		dx = 1.0;
        g = 1.0e+37; 
 
zjn:    u = a[0];
		v = 0.0;
        for(i = 1; i <= k; i++)
        { 
			p = u * x1;
			q = v * y1;
            pq = (u + v) * (x1 + y1);
            u = p - q + a[i];
			v = pq - p - q;
        }
        g1 = u * u + v * v;

        if(g1 >= g)
        { 
			if(is != 0)
			{ 
				it = 1;
				g65(x, y, x1, y1, dx, dy, dd, dc, c, k, is, it);
				if(it == 0)
					goto zjn;
			}
			else
			{
				g60(t, x, y, x1, y1, dx, dy, p, q, k, it);
				if(t >= 1.0e-03)  goto zjn;
				if(g > 1.0e-18)
				{
					it = 0;
					g65(x, y, x1, y1, dx, dy, dd, dc, c, k, is, it);
					if(it == 0) goto zjn;
				}

			}
			g90(cx , a, x, y, p, q, w, k);
        }
        else
        {
			g = g1;
			x = x1;
			y = y1;
			is = 0;
            if(g <= 1.0e-22)	
				g90(cx, a, x, y, p, q, w, k);
            else
            {
				u1 = k * a[0]; 
				v1 = 0.0;
                for(i = 2; i <= k; i++)
                {
					p = u1 * x;
					q = v1 * y;
					pq = (u1 + v1) * (x + y);
                    u1 = p - q + (k - i + 1) * a[i - 1];
                    v1 = pq - p - q;
                }
                p = u1 * u1 + v1 * v1;
                if(p <= 1.0e-20)
                {
					it = 0;
                    g65(x, y, x1, y1, dx, dy, dd, dc, c, k, is, it);
                    if(it == 0)	goto zjn;
                    g90(cx, a, x, y, p, q, w, k);
                }
                else
                {
					dx = (u * u1 + v * v1) / p;
                    dy = (u1 * v - v1 * u) / p;
                    t = 1.0 + 4.0 / k;
                    g60(t, x, y, x1, y1, dx, dy, p, q, k, it);
                    if(t >= 1.0e-03) goto zjn;
                    if(g > 1.0e-18)
                    {
						it = 0;
                        g65(x, y, x1, y1, dx, dy, dd, dc, c, k, is, it);
                        if(it == 0) goto zjn;
                    }
                    g90(cx, a, x, y, p, q, w, k);
                }
            }
        }

        if(k == 1) jt = 0;
        else jt = 1;
    }
    return(1);
}

//牛顿下山(NewtonHillDown)法求解复系数代数方程全部根(实根和复根)
//重载RootNewtonHillDown()
template <class _Ty>
inline int 
RootNewtonHillDown(complex<_Ty> a[], valarray<complex<_Ty> >& cx)
{
	int m, i, jt, k, is, it, n;
    _Ty t, x, y, x1, y1, dx, dy, dd, dc, c, g, g1, mo, sb, xb;
	_Ty w, u, p, q, pq, v, u1, v1;
    
	n = cx.size();			//多项式方程的次数
    m = n;
	
    while((m > 0) && (FloatEqual(Abs(a[m]), 0))) m--;

    if(m <= 0) return(-1);

    mo = Abs(a[m]);

	for(i = 0; i <= m; i++) a[i] /= mo;

    for(i = 0; i <= m / 2; i++)
    {
		swap(a[i], a[m-i]);		
	}

    k = m;
	is = 0;
	w = 1.0;
    jt = 1;
    while(jt == 1)
    {
		pq = Abs(a[k]);
		while(pq < LONGDOUBLEERROR)
        {
			cx[k - 1] = complex<_Ty>(0.0, 0.0);
			k--;
            if(k == 1)
            {
				mo = a[0].real() * a[0].real() + a[0].imag() * a[0].imag();
				sb = -w * (a[0].real() * a[1].real() + a[0].imag() * a[1].imag()) / mo;
				xb =  w * (a[1].real() * a[0].real() - a[0].imag() * a[1].imag()) / mo;
				cx[0] = complex<_Ty>(sb, xb);
                return(1); 
            }
			pq = Abs(a[k]);
        }

		q = log(pq);
		q = q / (1.0 * k);
		q = exp(q);
        p = q;
		w = w * p;
        
		for(i = 1; i <= k; i++)
        {
			a[i] /= q;
			q *= p;
		}
        
		x = 0.0001;
		x1 = x; 
		y = 0.2;
		y1 = y; 
		dx = 1.0;
        g = 1.0e+37; 
 
zjn:   
		u = a[0].real();
		v = a[0].imag();
        
		for(i = 1; i <= k; i ++)
        { 
			p = u * x1;
			q = v * y1;
            pq = (u + v) * (x1 + y1);
            u = p - q + a[i].real();
			v = pq - p - q + a[i].imag();
        }
        g1 = u * u + v * v;

        if(g1 >= g)
        { 
			if(is != 0)
			{ 
				it = 1;
				g65(x, y, x1, y1, dx, dy, dd, dc, c, k, is, it);
				if(it == 0)
					goto zjn;
			}
			else
			{
				gg60(t, x, y, x1, y1, dx, dy, p, q, k, it);
				if(t >= 1.0e-03)  goto zjn;
				if(g > 1.0e-18)
				{
					it = 0;
					g65(x, y, x1, y1, dx, dy, dd, dc, c, k, is, it);
					if(it == 0) goto zjn;
				}
			}

			gg90(cx, a, x, y, p, q, w, k);
        }
        else
        {
			g = g1;
			x = x1;
			y = y1;
			is = 0;

            if(g <= 1.0e-22)	
				gg90(cx, a, x, y, p, q, w, k);
            else
            {
				u1 = k * a[0].real(); 
				v1 = a[0].imag();
                for(i = 2; i <= k; i++)
                {
					p = u1 * x;
					q = v1 * y;
					pq = (u1 + v1) * (x + y);
                    u1 = p - q + (k - i + 1) * a[i - 1].real();
                    v1 = pq - p - q + (k - i + 1) * a[i - 1].imag();
                }
                p = u1 * u1 + v1 * v1;
                if(p <= 1.0e-20)
                {
					it = 0;
                    g65(x, y, x1, y1, dx, dy, dd, dc, c, k, is, it);
                    if(it == 0)	goto zjn;
                    gg90(cx, a, x, y, p, q, w, k);
                }
                else
                {
					dx = (u * u1 + v * v1) / p;
                    dy = (u1 * v - v1 * u) / p;
                    t = 1.0 + 4.0 / k;
                    gg60(t, x, y, x1, y1, dx, dy, p, q, k, it);
                    if(t >= 1.0e-03) goto zjn;
                    if(g > 1.0e-18)
                    {
						it = 0;
                        g65(x, y, x1, y1, dx, dy, dd, dc, c, k, is, it);
                        if(it == 0) goto zjn;
                    }
                    gg90(cx, a, x, y, p, q, w, k);
                }
            }
        }
        if(k == 1) jt = 0;
        else jt = 1;
    }
    return 1;
}

//牛顿下山法(NewtonHillDown)子函数
template <class _Ty>
inline void
g60(_Ty& t, _Ty& x, _Ty& y, _Ty& x1, _Ty& y1, _Ty& dx, _Ty& dy, 
	_Ty& p, _Ty& q, int& k, int& it)
{ 
	it = 1;
    while(it == 1)
    {
		t = t / 1.67;
		it = 0;
        x1 = x - t * dx;
        y1 = y - t * dy;
        if(k >= 50)
		{
			p = sqrt(x1 * x1 + y1 * y1);
            q = exp(85.0 / k);
            if(p >= q) it = 1;
        }
    }
}

//牛顿下山法(NewtonHillDown)子函数
template <class _Ty>
inline void
gg60(_Ty& t, _Ty& x, _Ty& y, _Ty& x1, _Ty& y1, _Ty& dx, _Ty& dy, 
	_Ty& p, _Ty& q, int& k, int& it)
{ 
	it = 1;
    while(it == 1)
    {
		t = t / 1.67;
		it = 0;
        x1 = x - t * dx;
        y1 = y - t * dy;
        if(k >= 30)
		{
			p = sqrt(x1 * x1 + y1 * y1);
            q = exp(75.0 / k);
            if(p >= q) it = 1;
        }
    }
	return;
}

//牛顿下山法(NewtonHillDown)子函数
template <class _Ty, class _Tz>
inline void
g90(valarray<_Tz> &cx, valarray<_Ty> &a, _Ty& x, _Ty& y, _Ty& p, _Ty& q, _Ty& w, int &k)
{
	int i;
    if(Abs(y) <= 1.0e-06)
    {
		p = -x;
		y = 0.0;
		q = 0.0;
	}
    else
    {
		p = -2.0 * x;
		q = x * x + y * y;
        cx[k - 1] = complex<_Ty>(x * w, -y * w);
        k = k-1;
    }
    for(i = 1; i <= k; i ++)
    {
		a[i] = a[i] - a[i - 1] * p;
        a[i + 1] = a[i + 1] - a[i - 1] * q;
    }
    cx[k - 1] = complex<_Ty>(x * w, y * w);
    k = k - 1;
    if(k == 1)
		cx[0] = complex<_Ty>(-a[1] * w / a[0], 0.0);
}

//牛顿下山法(NewtonHillDown)子函数
template <class _Ty, class _Tz>
inline void
gg90(valarray<_Tz> &cx, complex<_Ty> a[], _Ty& x, _Ty& y, _Ty& p, _Ty& q, _Ty& w, int &k)
{
	int i;
    double mo, sb, xb;    
	for (i = 1; i <= k; i ++)
	{ 
	  xb = a[i].real() + a[i-1].real() * x - a[i - 1].imag() * y;
	  sb = a[i].imag() + a[i-1].real() * y + a[i - 1].imag() * x;
	  a[i] = complex<_Ty>(xb, sb); 
	}
    
    cx[k - 1] = complex<_Ty>(x * w, y * w); 
	k --;
	if(k == 1)
	{
		mo = Abs(a[0])*Abs(a[0]);
		sb = -w * (a[0].real() * a[1].real() + a[0].imag() * a[1].imag()) / mo;
		xb =  w * (a[1].real() * a[0].imag() - a[0].real() * a[1].imag()) / mo;
		cx[0] = complex<_Ty>(sb, xb);		
	}
    return;
}

//牛顿下山法(NewtonHillDown)子函数
template <class _Ty>
inline void
g65(_Ty& x, _Ty& y, _Ty& x1, _Ty& y1, _Ty& dx, _Ty& dy, _Ty& dd, _Ty& dc, _Ty& c, 
	int& k, int& is, int& it)
{
	if(it == 0)
    {
		is = 1;
        dd = sqrt(dx * dx + dy * dy);
        if(dd > 1.0) dd = 1.0;
        dc = 6.28 / (4.5 * k);
		c = 0.0;
    }
    while(true)
    {
		c += dc;
        dx = dd * cos(c); 
		dy = dd * sin(c);
        x1 = x + dx; 
		y1 = y + dy;
        if(c <= 6.29)
        {
			it = 0;
			return;
		}
		dd /= 1.67;
		if(dd <= 1.0e-07)
		{
			it = 1;
			return;
		}
		c = 0.0;
	}
}

//梯度(Gradient)法(最速下降)求解非线性方程组一组实根
template <class _Ty>
inline int 
RootGradient(_Ty eps, valarray<_Ty>& x, size_t js)
{ 
    int r, j;
    _Ty f, d, s;
	
	int n = x.size();			//方程的个数，也是未知数的个数
	std::valarray<_Ty> y(n);	//定义一维数组对象y

    r = js;
    f = FunctionValueRG(x,y);
    while(f>=eps)
    {
		r = r - 1;
        if(r==0) return(js);

        d = 0;
        
		for(j=0; j<n; j++)	d=d+y[j]*y[j];

        if(FloatEqual(d,0)) return(-1);

        s = f / d;
        for(j=0; j<n; j++)	x[j]=x[j]-s*y[j];

        f=FunctionValueRG(x, y);
    }
    return(js-r);
}

//拟牛顿(QuasiNewton)法求解非线性方程组一组实根
template <class _Ty>
inline int 
RootQuasiNewton(_Ty eps, _Ty t, _Ty h, valarray<_Ty>& x, int k)
{
    int i, j, l;
    _Ty am, z, beta, d;
	
	int n = x.size();		//方程组中各方程个数，也是未知数个数
	
	matrix<_Ty> a(n, n);	//定义二维数组对象a
	valarray<_Ty> b(n);		//定义一维数组对象b
	valarray<_Ty> y(n);		//定义一维数组对象y
	    
	l = k;
	am = 1.0 + eps;
    while(am>=eps)
    {
		FunctionValueRSN(x, b);
        am = 0.0;
        for(i=0; i<n; i++)
        {
			z = Abs(b[i]);
            if(z>am) am = z;
        }
        if(am>=eps)
        {
			l--;
            if(l==0)
			{
				cout << "Fail!" << endl;
				return(0);
			}

            for(j=0; j<n; j++)
            {
				z = x[j];
				x[j] += h;
                FunctionValueRSN(x, y);
                for(i=0; i<n; i++) 
					a(i,j) = y[i];
                x[j] = z;
            }
			if(LE_TotalChoiceGauss(a, b) == 0)
			{
				cout << "Fail!" << endl;
				return(-1);					//矩阵a奇异
			}

            beta = 1.0;

            for(i=0; i<n; i++) beta -= b[i];
            
			if(FloatEqual(beta,0))
			{
				cout << "Fail!" << endl;
				return(-2);					//beta=0
			}

            d = h / beta;
            
			for(i=0; i<n; i++) x[i]=x[i]-d*b[i];

            h = t * h;
        }
    }
    return(k-l);					//正常结束，返回迭代次数
}

//非线性方程组最小二乘解的广义逆法
template <class _Ty>
int RootLeastSquareGeneralizedInverse(int m, _Ty eps1, _Ty eps2, valarray<_Ty>& x, int ka)
{
    int i, j, k, l, kk, jt;
	bool yn;
    _Ty y[10],b[10],alpha,z,h2,y1,y2,y3,y0,h1;
	
	int n = x.size();		//非线性方程组中未知数个数
    
	matrix<_Ty> p(m, n);
    matrix<_Ty> pp(n, m);
    matrix<_Ty> u(m, m);
    matrix<_Ty> v(n, n);

	valarray<_Ty> w(ka);
	valarray<_Ty> d(m);
	valarray<_Ty> dx(n);

    l = 60; 
	alpha = 1.0;

    while(l > 0)
	{
		FunctionValueRLSGI(x, d);	//计算非线性方程组左边函数值
        JacobiMatrix(x, p);			//计算雅可比矩阵
        
		//最小二乘的广义逆法求解线性方程组
		jt = LE_LinearLeastSquareGeneralizedInverse(p, d, dx, pp, eps2, u, v);

        if(jt < 0) return(jt);
        
		j = 0;
		jt = 1;
		h2 = 0.0;
        
		while(jt == 1)
		{
			jt = 0;
            
			if(j < 3) z = alpha + 0.01 * j;
            else z = h2;
            
			for(i = 0; i < n; i++) w[i] = x[i] - z * dx[i];
            
			FunctionValueRLSGI(w, d);
            y1 = 0.0;
            
			for(i = 0; i < m; i++) y1 = y1 + d[i] * d[i];
            for(i = 0; i < n; i++) w[i] = x[i] - (z + 0.00001) * dx[i];
            
			FunctionValueRLSGI(w, d);
            
			y2 = 0.0;
            for(i = 0; i < m; i++) y2 = y2 + d[i] * d[i];
            
			y0 =(y2 - y1) / 0.00001;
           
			if(Abs(y0) > 1.0e-10)
			{
				h1 = y0;
				h2 = z;
                
				if(j == 0)
				{ 
					y[0] = h1;
					b[0] = h2;
				}
                else
				{
					y[j] = h1;
					kk = k = 0;
					
                    while((kk == 0) && (k < j))
					{
						y3 = h2 - b[k];
						
						yn = FloatEqual(y3, 0);
						if(yn) kk = 1;
						else h2 = (h1 - y[k]) / y3;
						
						k++;
					}
                    
					b[j] = h2;
                    
					if(kk != 0) b[j] = 1.0e+35;
                    
					h2 = 0.0;
                    
					for(k = j - 1; k >= 0; k--)
						h2 = -y[k] / (b[k + 1] + h2);
                    h2 = h2 + b[0];
				}
                
				j++;
                
				if(j <= 7) jt = 1;
                else z = h2;
			}
		}
        
		alpha = z;
		y1 = y2 = 0.0;
       
		for(i = 0; i <= n-1; i ++)
		{
			dx[i] = -alpha * dx[i];
            x[i] = x[i] + dx[i];
            y1 = y1 + Abs(dx[i]);
            y2 = y2 + Abs(x[i]);
		}
        
		if(y1 < eps1 * y2)	return(1);

        l--;
	}

    return(0);
}

//蒙特卡洛(MonteCarlo)法求解f(x)=0的一个实根
//f(x)的自变量为与系数都为实数
template <class _Ty>
inline void 
RootMonteCarloReal(_Ty& x, _Ty b, int m, _Ty eps)
{
   // extern double mrnd1();
    //extern double dmtclf();
    //int k;
    _Ty x1, y1;
    _Ty a = b;
	size_t k = 1;
	double r = 1.0;
	_Ty xx = x;
	_Ty y = FunctionValueMCR(xx);	//计算函数值

    while(a>eps||FloatEqual(a,eps))
    {
		x1 = rand_01_One(r);		//取随机数
		x1 = -a + 2.0 * a * x1;
        x1 = xx + x1; 
		y1 = FunctionValueMCR(x1);	//计算函数值
        k++;
        if(Abs(y1)>Abs(y)||FloatEqual(Abs(y1),Abs(y)))
        {
			if(k>m)
			{
				k = 1; 
				a /= 2.0;
			}
		}
        else
        {
			k = 1;
			xx = x1;
			y = y1;
            if(Abs(y) < eps)
            {
				x=xx;
				exit(0);
			}
        }
    }
    x = xx;
}

//蒙特卡洛(MonteCarlo)法求解f(x)=0的一个复根
//f(x)的自变量为复数，或自变量与系数都为复数(不能都为实数)
template <class _Tz, class _Ty>
inline void 
RootMonteCarloComplex(_Tz& cxy, _Ty b, int m, _Ty eps)
{
    size_t k = 1;
	_Tz xxyy(cxy);
    _Ty a = b;
	double r(1);
	_Ty z, z1;

    z = FunctionModule(xxyy);

    while(a > eps || FloatEqual(a,eps))
    {
		_Ty tempx = -a + 2.0 *a * rand_01_One(r);
		_Ty tempy = -a + 2.0 *a * rand_01_One(r);

		_Tz x1y1(tempx,tempy);
		x1y1 += xxyy;

        z1 = FunctionModule(x1y1);
        
		k++;

        if(z1 > z || FloatEqual(z1,z))
        {
			if(k > m)
			{
				k = 1;
				a = a / 2.0;
			}
		}
        else
        {
			k = 1;
			xxyy = x1y1;
			z = z1;
            if(z < eps)
			{
				cxy = xxyy;
				exit(0);
			}
          }
      }
    cxy = xxyy;
}

//蒙特卡洛(MonteCarlo)法求解f(x)=0的一组实根
//f(x)的自变量为与系数都为实数
template <class _Ty>
inline void 
RootMonteCarloGroupReal(valarray<_Ty>& x, _Ty b, int m, _Ty eps)
{
    _Ty a=b;
	size_t k=1;
	double r=1.0;

	int n = x.size();		//方程个数，也是未知量的个数
	valarray<_Ty> y(n);

	_Ty z = FunctionModule(x);

    while(a>eps||FloatEqual(a,eps))
    {
		for(size_t i = 0; i < n; i++)
			y[i] = -a + 2.0 * a * rand_01_One(r) + x[i];

        _Ty z1 = FunctionModule(y);
        
		k++;

        if(z1 > z || FloatEqual(z1,z))
        {
			if(k > m)
			{
				k = 1;
				a = a / 2.0;
			}
		}
        else
        {
			k = 1; 

            for(i = 0; i < n; i++)	x[i] = y[i];

            z = z1;
            
			if(z < eps)	exit(0);
        }
    }
}

#endif //_NONLINEAREQUATION_INL

