// Integral.inl			数值积分头文件
// Ver 1.0.0.0
// 版权所有(C) 2002
// 最后修改: 2002.5.31.

#ifndef _INTEGRAL_INL	//避免多次编译
#define _INTEGRAL_INL

//变步长梯形法求积
template <class _Ty>
_Ty IntegralTrapezia(_Ty a, _Ty b, _Ty eps)
{
	_Ty t;
    _Ty fa = FunctionValueIT(a);
	_Ty fb = FunctionValueIT(b);
    int n(1);
	_Ty h = b - a;
    _Ty t1 = h * (fa + fb) / 2.0;
    _Ty p = eps + 1.0;
    while(p > eps || FloatEqual(p, eps))
    {
		_Ty s(0);
        for(int k=0; k<n; k++)
        {
			_Ty x = a + (k + 0.5) * h;
            s = s + FunctionValueIT(x);
          }
        t = (t1 + h * s) / 2.0;
        p = Abs(t1 - t);
        t1 = t;
		n = n + n;
		h = h / 2.0;
    }
    return(t);
}

//变步长辛卜生法求积
template <class _Ty>
_Ty IntegralSimpson1D(_Ty a, _Ty b, _Ty eps)
{
    _Ty s2;
    int n(1);
	_Ty h = b - a;			//步长
    _Ty t1 = h * (FunctionValueIS1D(a) + FunctionValueIS1D(b)) / 2.0;
    _Ty s1 = t1;
    _Ty ep = eps + 1.0;
    while(ep > eps || FloatEqual(ep, eps))
    {
		_Ty p(0);
        for(int k=0; k<n; k++)
        {
			_Ty x = a + (k + 0.5) * h;
            p = p + FunctionValueIS1D(x);
        }
        _Ty t2 = (t1 + h * p) / 2.0;
        s2 = (4.0 * t2 - t1) / 3.0;
        ep = Abs(s2 - s1);
        t1 = t2;
		s1 = s2;
		n = n + n;
		h = h / 2.0;
    }
    return(s2);
}

//自适应梯形法求积
template <class _Ty>
_Ty IntegralTrapeziaSelfAdapt(_Ty a, _Ty b, _Ty eps, _Ty d)
{
    valarray<_Ty> t(2);
    _Ty h = b - a;
	t[0]=0.0;
    _Ty f0 = FunctionValueITSA(a);
	_Ty f1 = FunctionValueITSA(b);
    _Ty t0 = h * (f0 + f1) / 2.0;
    ITSA(a,b,h,f0,f1,t0,eps,d,t);
    _Ty z = t[0];
	return(z);
}

//自适应梯形法求积辅助函数
template <class _Ty>
int ITSA(_Ty x0, _Ty x1, _Ty h, _Ty f0, _Ty f1, _Ty t0, 
								_Ty eps, _Ty d, valarray<_Ty>& t)
{ 
    _Ty x = x0 + h / 2.0;
	_Ty f = FunctionValueITSA(x);
    _Ty t1=h*(f0+f)/4.0;
	_Ty t2=h*(f+f1)/4.0;
    _Ty p=Abs(t0-(t1+t2));
    if((p<eps)||(h/2.0<d))
	{
		t[0]=t[0]+(t1+t2);
		return(0);
	}
    else
    {
		_Ty g=h/2.0; 
		_Ty eps1=eps/1.4;
        ITSA(x0,x,g,f0,f,t1,eps1,d,t);
        ITSA(x,x1,g,f,f1,t2,eps1,d,t);
        return(1);
    }
}

//龙贝格法求积
template <class _Ty>
_Ty IntegralRomberg(_Ty a, _Ty b, _Ty eps)
{
    _Ty y[10], q;
    _Ty h = b - a;			//步长
    y[0] = h * (FunctionValueR(a) + FunctionValueR(b)) / 2.0;
    int m(1), n(1);
	_Ty ep=eps+1.0;
    while((ep>eps || FloatEqual(ep,eps)) && (m <= 9))
    {
		_Ty p(0);
        for(int i=0;i<n;i++)
        {
			_Ty x = a + (i + 0.5) * h;
            p = p + FunctionValueR(x);
        }
        p = (y[0] + h * p) / 2.0;
        _Ty s(1);
        for(int k=1; k<=m; k++)
        {
			s = 4.0 * s;
            q = (s * p - y[k-1]) / (s - 1.0);
            y[k-1] = p;
			p = q;
        }
        ep = Abs(q - y[m-1]);
        m = m + 1;
		y[m-1] = q; 
		n = n + n;
		h = h / 2.0;
    }
    return(q);
}

//一维连分式法求积
template <class _Ty>
_Ty IntegralFraction1D(_Ty a, _Ty b, _Ty eps)
{
    valarray<_Ty> h(10), bb(10);
    int m(1), n(1), k, l, j;
    _Ty g, hh = b - a;
	h[0] = hh;
    _Ty t1 = hh * (FunctionValueF1D(a) + FunctionValueF1D(b)) / 2.0;
    _Ty s1 = t1;
	bb[0]=s1;
	_Ty ep = 1.0 + eps;
    while((ep>eps || FloatEqual(ep,eps)) && (m <= 9))
    { 
		_Ty s(0);
        for(k=0; k<n; k++)
        {
			_Ty x = a + (k + 0.5) * hh;
            s = s + FunctionValueF1D(x);
        }
        _Ty t2 = (t1 + hh * s) / 2.0;
        m++;
        h[m-1] = h[m-2] / 2.0;
        g = t2;
        l = 0;
		j = 2;
        while((l==0) && (j<=m))
        {
			s = g - bb[j-2];
            if(FloatEqual(s,0)) l = 1;
            else g = (h[m-1] - h[j-2]) / s;
            j++;
        }
        bb[m-1] = g;
        if(l!=0) bb[m-1] = 1.0e+35;
        g = bb[m-1];
        for(j=m; j>=2; j--)	g = bb[j-2] - h[j-2] / g;
        ep = Abs(g-s1);
        s1 = g;
		t1 = t2;
		hh = hh / 2.0;
		n = n + n;
    }
    return(g);
}

//高振荡函数法求积
template <class _Ty>
void IntegralSurge(_Ty a, _Ty b, int m, valarray<_Ty>& fa, 
								valarray<_Ty>& fb, valarray<_Ty>& s)
{
    _Ty sa[4],sb[4],ca[4],cb[4];

	int n = fa.size();		//给定积分区间两端点上的导数最高阶数+1
	int mm=1;

    _Ty sma=sin(m*a);
	_Ty smb=sin(m*b);
    _Ty cma=cos(m*a);
	_Ty cmb=cos(m*b);
    
	sa[0] = ca[3] = sma;
	sa[1] = ca[0] = cma;
	sa[2] = ca[1] = -sma; 
	sa[3] = ca[2] = -cma;
    sb[0] = cb[3] = smb;
	sb[1] = cb[0] = cmb;
	sb[2] = cb[1] = -smb;
	sb[3] = cb[2] = -cmb;

    s[0] = s[1] = 0.0; 

    for(int k=0; k<n; k++)
    {
		int j = k;
        while(j>=4) j = j - 4;
        mm = mm * m;
        s[0] = s[0]+(fb[k]*sb[j]-fa[k]*sa[j])/(1.0*mm);
        s[1] = s[1]+(fb[k]*cb[j]-fa[k]*ca[j])/(1.0*mm);
    }
    
	s[1] = -s[1];
}

//勒让德-高斯法求积
template <class _Ty>
_Ty IntegralLegendreGauss(_Ty a, _Ty b, _Ty eps)
{
    _Ty t[5] = 
	{
		-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459
	};
    
	_Ty c[5] = 
	{
		0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851
	};

    int m(1);
	_Ty g;
    _Ty h = b - a;
	_Ty s = Abs(0.001 * h);
    _Ty p = 1.0e+35;
	_Ty ep = eps + 1.0;

    while((ep > eps || FloatEqual(ep,eps)) && (Abs(h) > s))
    {
		g = 0;
        for(int i=1; i<=m; i++)
        {
			_Ty aa = a + (i - 1.0) * h;
			_Ty bb = a + i * h;
            _Ty w(0);
            for(int j=0; j<5; j++)
            {
				_Ty x = ((bb - aa) * t[j] + (bb + aa)) / 2.0;
                w = w + FunctionValueLegG(x) * c[j];
            }
            g += w;
        }
        g = g * h / 2.0;
        ep = Abs(g - p) / (1.0 + Abs(g));
        p = g;
		m++;
		h = (b - a) / m;
    }
    return(g);
}

//拉盖尔-高斯法求积
template <class _Ty>
void IntegralLaguerreGauss(_Ty& dValue)
{
	_Ty x;
	
	dValue = 0.0;

    long double t[5] = 
	{ 
		0.26355990, 1.41340290, 3.59642600, 7.08580990, 12.64080000
	};
    
	long double c[5] = 
	{
		0.6790941054, 1.638487956, 2.769426772, 4.315944000, 7.104896230
	};
    
	for(int i=0; i<5; i++)
    {
		x = t[i];
		dValue = dValue + c[i] * FunctionValueLagG(x); 
	}
}

//埃尔米特-高斯法求积
template <class _Ty >
void IntegralHermiteGauss(_Ty& dValue)
{
    _Ty t[5] = 
	{
		-2.02018200, -0.95857190,	0.0, 0.95857190, 2.02018200
	};
    
	_Ty c[5] = 
	{
		1.181469599, 0.9865791417,0.9453089237, 0.9865791417,1.181469599
	};
    
	dValue = 0;
    
	for(int i=0; i<5; i++)
    {
		_Ty x = t[i];
		dValue = dValue + c[i] * FunctionValueHG(x);
	}
}

//切比雪夫法求积
template <class _Ty >
_Ty IntegralChebyshev(_Ty a, _Ty b, _Ty eps)
{
	_Ty t[5]={-0.8324975, -0.3745414, 0.0, 0.3745414, 0.8324975};
	int m(1);
	_Ty g;
    _Ty h = b - a;
	_Ty d = Abs(0.001 * h);
    _Ty p = 1.0e+35;
	_Ty ep = eps + 1.0;
    while((ep > eps || FloatEqual(ep,eps)) && (Abs(h) > d))
    {
		g = 0.0;
        for(int i=1; i<=m; i++)
        {
			_Ty aa = a + (i - 1.0) * h;
			_Ty bb = a + i * h;
            _Ty s(0);
			for(int j=0; j<5; j++)
            {
				_Ty x = ((bb - aa) * t[j] + (bb + aa)) / 2.0;
                s = s + FunctionValueCb(x);
            }
            g += s;
        }
        g = g * h / 5.0;
        ep = Abs(g - p) / (1.0 + Abs(g));
        p = g;
		m++;
		h = (b - a) / m;
    }
    return(g);
}

//一维蒙特卡洛法求积
template <class _Ty >
_Ty IntegralMonteCarlo1D(_Ty a, _Ty b)
{
    _Ty r(1), s(0), d(10000);

    for(int m=0; m<10000; m++)
    {
		_Ty x = a + (b - a) * rand_01_One(r);	//取随机数

        s = s + FunctionValueMC1D(x) / d;		//调用被积函数值
    }
    
	s = s * (b - a);
    
	return(s);
}

//二重变步长辛卜生法求积
template <class _Ty>
_Ty IntegralSimpson2D(_Ty a, _Ty b, _Ty eps)
{
	int n(1);
	_Ty h = 0.5 * (b - a);			//步长
	_Ty d = Abs((b - a) * 1.0e-06);
    _Ty s1 = IntegralSimp2(a, eps);	//调用IntegralSimp2函数
	_Ty s2 = IntegralSimp2(b, eps);
    _Ty t1 = h * (s1 + s2);
    _Ty s0(1.0e+35), s;
	_Ty ep = eps + 1.0;
    while((ep>eps||FloatEqual(ep,eps))&&((Abs(h)>d)||(n<16)))
	{
		_Ty x = a - h;
		_Ty t2 = 0.5 * t1;
        for(int j=1; j<=n; j++)
        {
			x = x + 2.0 * h;
            t2 = t2 + h * IntegralSimp2(x,eps);
        }
        s = (4.0 * t2 - t1) / 3.0;
        ep = Abs(s - s0) / (1.0 + Abs(s));
        n = n + n;
		s0 = s;
		t1 = t2; 
		h = h * 0.5;
    }
    return(s);
}

//二重变步长辛卜生法求积辅助函数
template <class _Ty>
_Ty IntegralSimp2(_Ty x, _Ty eps)
{
    valarray<_Ty> y(2);
	int n(1);
    FunctionBoundaryIS2D(x, y);
    _Ty h = 0.5 *(y[1] - y[0]);
    _Ty d = Abs(h * 2.0e-06);
    _Ty t1 = h * (FunctionValueIS2D(x, y[0]) + FunctionValueIS2D(x, y[1]));
    _Ty ep = 1.0 + eps;
	_Ty g0(1.0e+35), g;
    while((ep>eps||FloatEqual(ep,eps))&&((Abs(h)>d)||(n<16)))
    {
		_Ty yy = y[0] - h;
        _Ty t2 = 0.5 * t1;
        for(int i=1; i<=n; i++)
        {
			yy = yy + 2.0 * h;
            t2 = t2 + h * FunctionValueIS2D(x, yy);
        }
        g = (4.0 * t2 - t1) / 3.0;
        ep = Abs(g - g0) / (1.0 + Abs(g));
        n = n + n;
		g0 = g;
		t1 = t2;
		h = 0.5 * h;
    }
    return(g);
}

//多重高斯法求积
template <class _Ty >
_Ty IntegralGaussMD(valarray<_Ty>& js)
{
    valarray<_Ty> y(2);
	_Ty p, s;
    _Ty t[5] = 
	{
		-0.9061798459,-0.5384693101,0.0, 0.5384693101,0.9061798459
	};
    
	_Ty c[5] = 
	{
		0.2369268851,0.4786286705,0.5688888889, 0.4786286705,0.2369268851
	};
	
	int n = js.size();		//积分重数

	valarray<int> is(2*(n+1));
    valarray<_Ty> x(n);
    valarray<_Ty> a(2*(n+1));
    valarray<_Ty> b(n+1);
	int m(1), l(1);
    
	a[n] = a[2*n+1] = 1.0;

    while(l==1)
    {
		for(int j=m; j<=n; j++)
        {
			FunctionBoundaryGMD(j-1,x,y);
            a[j-1]=0.5*(y[1]-y[0])/js[j-1];
            b[j-1]=a[j-1]+y[0];
            x[j-1]=a[j-1]*t[0]+b[j-1];
            a[n+j]=0.0;
            is[j-1]=1; is[n+j]=1;
        }
        j = n;
		int q = 1;
        while(q==1)
        {
			int k = is[j-1];
            if(j==n) p = FunctionValueGMD(x);
            else p = 1.0;
            a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];
            is[j-1] = is[j-1] + 1;
            if(is[j-1] > 5)
              if(is[n+j] >= js[j-1])
              {
				  j = j - 1;
				  q = 1;
                  if(j==0)
                  {
					  s = a[n+1] * a[0]; 
                      return(s);
                  }
              }
              else
              {
				  is[n+j]=is[n+j]+1;
                  b[j-1]=b[j-1]+a[j-1]*2.0;
                  is[j-1]=1; k=is[j-1];
                  x[j-1]=a[j-1]*t[k-1]+b[j-1];
                  if(j==n) q = 1;
                  else q = 0;
              }
            else
            {
				k = is[j-1];
                x[j-1] = a[j-1]*t[k-1]+b[j-1];
                if(j==n) q = 1;
                else q = 0;
            }
        }
        m = j + 1;
    }
}

//二重连分式法求积
template <class _Ty >
_Ty IntegralFraction2D(_Ty a, _Ty b, _Ty eps)
{
    _Ty bb[10],h[10],x,s;
	int m(1), n(1);
    _Ty hh = b - a;
	h[0] = hh;
    _Ty s1 = IntegralFraction2D2(a, eps);
	_Ty s2 = IntegralFraction2D2(b, eps);
    _Ty t1 = hh * (s1 + s2) / 2.0;
    _Ty s0 = t1;
	bb[0] = t1;
	_Ty ep = 1.0 + eps;
    while((ep > eps || FloatEqual(ep,eps)) && ( m <= 9))
    {
		_Ty t2 = 0.5 * t1;
        for(int k=0; k<n; k++)
        {
			x = a + (k + 0.5) * hh;
            s1 = IntegralFraction2D2(x,eps);
            t2 = t2 + 0.5 * s1 * hh;
        }
        m = m + 1;
        h[m-1] = h[m-2] / 2.0;
        _Ty g = t2;
		int l(0), j(2);
        while((l==0)&&(j<=m))
         {
			s = g - bb[j-2];
            if(FloatEqual(s,0)) l = 1;
            else g = (h[m-1] - h[j-2]) / s;
            j = j + 1;
        }
        bb[m-1] = g;
        if(l!=0) bb[m-1] = 1.0e+35;
        s = bb[m-1];
        for(j=m; j>=2; j--) s = bb[j-2] -h [j-2] / s;
        ep = Abs(s-s0) / (1.0 + Abs(s));
        n = n + n;
		t1 = t2; 
		s0 = s;
		hh = hh / 2.0;
    }
    return(s);
}

//二重连分式法求积辅助函数
template <class _Ty>
_Ty IntegralFraction2D2(_Ty x, _Ty eps)
{
    _Ty b[10],h[10],yy,s;
	valarray<_Ty> y(2);
    int m(1), n(1);
    FunctionBoundaryF2D(x, y);
    _Ty hh = y[1] - y[0];
	h[0] = hh;
    _Ty t1=0.5*hh*(FunctionValueF2D(x,y[0])+FunctionValueF2D(x,y[1]));
    _Ty s0 = t1;
	b[0] = t1;
	_Ty ep = 1.0 + eps;
    while((ep > eps || FloatEqual(ep, eps)) && (m <= 9))
    {
		_Ty t2 = 0.5 * t1;
        for(int k=0; k<n; k++)
        {
			yy = y[0] +(k + 0.5) * hh;
            t2 = t2 + 0.5 * hh * FunctionValueF2D(x, yy);
		}
        m = m + 1;
        h[m-1] = h[m-2] / 2.0;
        _Ty g = t2;
		int l(0), j(2);
        while((l==0)&&(j<=m))
        {
			s = g - b[j-2];
            if(FloatEqual(s,0)) l = 1;
            else g = (h[m-1] - h[j-2]) / s;
            j = j + 1;
        }
        b[m-1] = g;
        if(l!=0) b[m-1] = 1.0e+35;
        s = b[m-1];
        for(j=m; j>=2; j--) s = b[j-2] - h[j-2] / s;
        ep = Abs(s - s0) / (1.0 + Abs(s));
        n = n + n;
		t1 = t2;
		s0 = s;
		hh = 0.5 * hh;
    }
    return(s);
}

//多重蒙特卡洛法求积
template <class _Ty >
_Ty IntegralMonteCarlo2D(valarray<_Ty>& a, valarray<_Ty>& b)
{
	int i;
	int n = a.size();		//积分的重数

	valarray<_Ty> x(n);
    _Ty r(1), d(10000), s(0);

    for(int m=0; m<10000; m++)
    {
		for(i=0; i<n; i++)
			x[i] = a[i] + (b[i] - a[i]) * rand_01_One(r);

        s = s + FunctionValueMC2D(n, x) / d;
    }
    
	for(i=0; i<n; i++)	s = s * (b[i] - a[i]);
    
	return(s);
}

#endif		// _INTEGRAL_INL

