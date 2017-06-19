//SpecialFunction.inl		特殊函数定义头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.


#ifndef _SPECIALFUNCTION_INL		//避免多次编译
#define _SPECIALFUNCTION_INL

//template <class _Ty = float>

//伽马函数
template <class _Ty>
_Ty GammaFunction(_Ty x)
{
    _Ty t;
    _Ty a[11] = 
	{
		 0.0000677106, -0.0003442342,  0.0015397681,
		-0.0024467480,  0.0109736958, -0.0002109075,
		 0.0742379071,  0.0815782188,  0.4118402518,
		 0.4227843370,  1.0
	};

    if(x<0.0 || FloatEqual(x,0.0))
    {
		cout << "Error **x<=0!" << endl; 
		return(-1.0);
	}
    _Ty y = x;
    if(y<1.0 || FloatEqual(y,1.0))
    {
		t=1.0/(y*(y+1.0)); 
		y=y+2.0;
	}
    else
		if(y<2.0 || FloatEqual(y,2.0))
		{
			t=1.0/y; 
			y=y+1.0;
		}
		else
			if(y<3.0 || FloatEqual(y,3.0)) t=1.0;
			else
			{
				t=1.0;
				while(y>3.0)
				{
					y=y-1.0;
					t=t*y;
				}
			}
    _Ty s=a[0]; 
	_Ty u = y - 2.0;
    for(int i=1; i<11; i++)	s = s * u + a[i];
    s = s * t;
    return(s);
}

//不完全伽马函数
template <class _Ty>
_Ty IncompleteGammaFunction(_Ty a, _Ty x)
{
	int n;
    _Ty p,d,s,s1,p0,q0,p1,q1;
    if((a<0.0 || FloatEqual(a,0.0)) || (x<0.0))
    { 
		if(a<0.0 || FloatEqual(a,0.0)) cout <<"Error **a<=0!" <<endl;
        if(x<0.0) cout << "Error **x<0!" << endl;
        return(-1.0);
    }
    if(FloatEqual(x,0.0)) return(0.0);
    if(x>1.0e+35) return(1.0);
    _Ty q = log(x); 
	q = a * q; 
	_Ty qq = exp(q);
    if(x<1.0+a)
    {
		p=a; 
		d=1.0/a; 
		s=d;
        for(n=1; n<101; n++)
        {
			p=1.0+p;
			d=d*x/p;
			s=s+d;
			if(Abs(d)<Abs(s)*FLOATERROR)
            {
				s=s*exp(-x)*qq/GammaFunction(a);
                return(s);
            }
        }
    }
    else
    { 
		s=1.0/x; 
		p0=0.0; 
		p1=1.0;
		q0=1.0;
		q1=x;
        for(n=1; n<101; n++)
        { 
			p0=p1+(n-a)*p0;
			q0=q1+(n-a)*q0;
            p=x*p0+n*p1;
			q=x*q0+n*q1;
            if(FloatNotEqual(q,0.0))
            {
				s1=p/q; 
				p1=p; 
				q1=q;
                if(Abs((s1-s)/s1)<FLOATERROR)
                { 
					s=s1*exp(-x)*qq/GammaFunction(a);
                    return(1.0-s);
                }
                s=s1;
            }
            p1=p; 
			q1=q;
        }
    }
    cout << "a too large !" << endl;
    s=1.0-s*exp(-x)*qq/GammaFunction(a);
    return(s);
}

//误差函数
template <class _Ty>
_Ty ErrorFunction(_Ty x)
{
	_Ty y;

    if(x>0.0 || FloatEqual(x,0.0))
      y=IncompleteGammaFunction(0.5,x*x);
    else
      y=-IncompleteGammaFunction(0.5,x*x);
    return(y);
}

//第一类整数阶贝塞尔函数
template <class _Ty>
_Ty IntegerBessel1stFunction(int n, _Ty x)
{
	int i,m;
    _Ty y,z,p,q,s,b0,b1;
    _Ty a[6]={ 57568490574.0,-13362590354.0,
             651619640.7,-11214424.18,77392.33017,-184.9052456};
    _Ty b[6]={ 57568490411.0,1029532985.0,
             9494680.718,59272.64853,267.8532712,1.0};
    _Ty c[6]={ 72362614232.0,-7895059235.0,
             242396853.1,-2972611.439,15704.4826,-30.16036606};
    _Ty d[6]={ 144725228443.0,2300535178.0,
             18583304.74,99447.43394,376.9991397,1.0};
    _Ty e[5]={ 1.0,-0.1098628627e-02,
             0.2734510407e-04,-0.2073370639e-05,0.2093887211e-06};
    _Ty f[5]={ -0.1562499995e-01,
             0.1430488765e-03,-0.6911147651e-05,
             0.7621095161e-06,-0.934935152e-07};
    _Ty g[5]={ 1.0,0.183105e-02,
             -0.3516396496e-04,0.2457520174e-05, -0.240337019e-06};
    _Ty h[5]={ 0.4687499995e-01,
             -0.2002690873e-03,0.8449199096e-05,
             -0.88228987e-06,0.105787412e-06};
    _Ty t=Abs(x);
    if(n<0) n=-n;
    if(n!=1)
    { 
		if(t<8.0)
        {
			y=t*t;
			p=a[5];
			q=b[5];
			for(i=4; i>=0; i--)
            { 
				p=p*y+a[i];
				q=q*y+b[i];
			}
			p=p/q;
        }
        else
        { 
			z=8.0/t;
			y=z*z;
            p=e[4];
			q=f[4];
            for(i=3; i>=0; i--)
            {
				p=p*y+e[i];
				q=q*y+f[i];
			}
            s=t-0.785398164;
            p=p*cos(s)-z*q*sin(s);
            p=p*sqrt(0.636619772/t);
        }
    }
    if(n==0) return(p);
    b0=p;
    if(t<8.0)
    { 
		y=t*t; 
		p=c[5];
		q=d[5];
        for(i=4; i>=0; i--)
        {
			p=p*y+c[i];
			q=q*y+d[i];
		}
        p=x*p/q;
    }
    else
    {
		z=8.0/t; 
		y=z*z;
        p=g[4];
		q=h[4];
        for(i=3; i>=0; i--)
        {
			p=p*y+g[i];
			q=q*y+h[i];
		}
        s=t-2.356194491;
        p=p*cos(s)-z*q*sin(s);
        p=p*x*sqrt(0.636619772/t)/t;
    }
    if(n==1) return(p);
    b1=p;
    if(x==0.0) return(0.0);
    s=2.0/t;
    if(t>1.0*n)
    {
		if(x<0.0) b1=-b1;
        for(i=1; i<n; i++)
        { 
			p=s*i*b1-b0; 
			b0=b1; 
			b1=p;
		}
    }
    else
    {
		m=(n+(int)sqrt(40.0*n))/2;
        m=2*m;
        p = q = b1 = 0.0;
		b0=1.0;
        for(i=m-1; i>=0; i--)
        { 
			t=s*(i+1)*b0-b1;
            b1=b0;
			b0=t;
            if(Abs(b0)>1.0e+10)
            {
				b0=b0*1.0e-10; 
				b1=b1*1.0e-10;
                p=p*1.0e-10; 
				q=q*1.0e-10;
            }
            if((i+2)%2==0) q=q+b0;
            if((i+1)==n) p=b1;
        }
        q=2.0*q-b0;
		p=p/q;
    }
    if((x<0.0)&&(n%2==1)) p=-p;
    return(p);
}

//第二类整数阶贝塞尔函数
template <class _Ty>
_Ty IntegerBessel2ndFunction(int n, _Ty x)
{
	int i;
    _Ty y,z,p,q,s,b0,b1;

    _Ty a[6]={ -2.957821389e+9,7.062834065e+9,
             -5.123598036e+8,1.087988129e+7,-8.632792757e+4,
             2.284622733e+2};
    _Ty b[6]={ 4.0076544269e+10,7.452499648e+8,
           7.189466438e+6,4.74472647e+4,2.261030244e+2,1.0};
    _Ty c[6]={ -4.900604943e+12,1.27527439e+12,
            -5.153438139e+10,7.349264551e+8,-4.237922726e+6,
             8.511937935e+3};
    _Ty d[7]={ 2.49958057e+13,4.244419664e+11,
            3.733650367e+9,2.245904002e+7,1.02042605e+5,
            3.549632885e+2,1.0};
    _Ty e[5]={ 1.0,-0.1098628627e-02,
             0.2734510407e-04,-0.2073370639e-05,
             0.2093887211e-06};
    _Ty f[5]={ -0.1562499995e-01,
             0.1430488765e-03,-0.6911147651e-05,
             0.7621095161e-06,-0.934935152e-07};
    _Ty g[5]={ 1.0,0.183105e-02,
             -0.3516396496e-04,0.2457520174e-05,
             -0.240337019e-06};
    _Ty h[5]={ 0.4687499995e-01,
             -0.2002690873e-03,0.8449199096e-05,
             -0.88228987e-06,0.105787412e-06};
    if(n<0) n=-n;
    if(x<0.0) x=-x;
    if(FloatEqual(x,0.0)) return(-1.0e+70);
    if(n!=1)
    { 
		if(x<8.0)
        {
			y=x*x; 
			p=a[5]; 
			q=b[5];
			for(i=4; i>=0; i--)
            {
				p=p*y+a[i];
				q=q*y+b[i];}
				p=p/q+0.636619772*IntegerBessel1stFunction(0,x)*log(x);
        }
        else
        { 
			z=8.0/x;
			y=z*z;
            p=e[4];
			q=f[4];
            for(i=3; i>=0; i--)
            {
				p=p*y+e[i];
				q=q*y+f[i];
			}
			s=x-0.785398164;
            p=p*sin(s)+z*q*cos(s);
            p=p*sqrt(0.636619772/x);
        }
    }
    if(n==0) return(p);
    b0=p;
    if(x<8.0)
    {
		y=x*x; 
		p=c[5];
		q=d[6];
        for(i=4; i>=0; i--)
        {
			p=p*y+c[i];
			q=q*y+d[i+1];
		}
        q=q*y+d[0];
        p=x*p/q+0.636619772*(IntegerBessel1stFunction(1,x)*log(x)-1.0/x);;
    }
    else
    { 
		z=8.0/x;
		y=z*z;
        p=g[4]; 
		q=h[4];
        for(i=3; i>=0; i--)
        {
			p=p*y+g[i];
			q=q*y+h[i];
		}
        s=x-2.356194491;
        p=p*sin(s)+z*q*cos(s);
        p=p*sqrt(0.636619772/x);
    }
    if(n==1) return(p);
    b1=p;
    s=2.0/x;
    for(i=1; i<n; i++)
    { 
		p=s*i*b1-b0; 
		b0=b1; 
		b1=p;
	}
    return(p);
}

//变形第一类整数阶贝塞尔函数
template <class _Ty>
_Ty TransformativeIntegerBessel1stFunction(int n,_Ty x)
{
	int i,m;
    _Ty t,y,p,b0,b1,q;
    _Ty a[7]={ 1.0,3.5156229,3.0899424,1.2067492,
                         0.2659732,0.0360768,0.0045813};
    _Ty b[7]={ 0.5,0.87890594,0.51498869,
              0.15084934,0.02658773,0.00301532,0.00032411};
    _Ty c[9]={ 0.39894228,0.01328592,0.00225319,
                        -0.00157565,0.00916281,-0.02057706,
                         0.02635537,-0.01647633,0.00392377};
    _Ty d[9]={ 0.39894228,-0.03988024,-0.00362018,
                        0.00163801,-0.01031555,0.02282967,
                        -0.02895312,0.01787654,-0.00420059};
    if(n<0) n=-n;
    t=Abs(x);
    if(n!=1)
    { 
		if(t<3.75)
        {
			y=(x/3.75)*(x/3.75);
			p=a[6];
            for(i=5; i>=0; i--)	p=p*y+a[i];
        }
        else
        { 
			y=3.75/t;
			p=c[8];
            for(i=7; i>=0; i--)	p=p*y+c[i];
            p=p*exp(t)/sqrt(t);
        }
    }
    if(n==0) return(p);
    q=p;
    if(t<3.75)
    {
		y=(x/3.75)*(x/3.75); 
		p=b[6];
        for(i=5; i>=0; i--) p=p*y+b[i];
        p=p*t;
    }
    else
    { 
		y=3.75/t;
		p=d[8];
        for(i=7; i>=0; i--) p=p*y+d[i];
        p=p*exp(t)/sqrt(t);
    }
    if(x<0.0) p=-p;
    if(n==1) return(p);
    if(FloatEqual(x,0.0)) return(0.0);
    y=2.0/t;
	t = b0 = 0.0;
	b1=1.0;
    m=n+(int)sqrt(40.0*n);
    m=2*m;
    for(i=m; i>0; i--)
    { 
		p=b0+i*y*b1; 
		b0=b1;
		b1=p;
        if(Abs(b1)>1.0e+10)
        {
			t=t*1.0e-10; 
			b0=b0*1.0e-10;
            b1=b1*1.0e-10;
        }
        if(i==n) t=b0;
    }
    p=t*q/b1;
    if((x<0.0)&&(n%2==1)) p=-p;
    return(p);
}

//变形第二类整数阶贝塞尔函数
template <class _Ty>
_Ty TransformativeIntegerBessel2ndFunction(int n, _Ty x)
{
	int i;
    _Ty y,p,b0,b1;

    _Ty a[7]={ -0.57721566,0.4227842,0.23069756,
           0.0348859,0.00262698,0.0001075,0.0000074};
    _Ty b[7]={ 1.0,0.15443144,-0.67278579,
           -0.18156897,-0.01919402,-0.00110404,-0.00004686};
    _Ty c[7]={ 1.25331414,-0.07832358,0.02189568,
           -0.01062446,0.00587872,-0.0025154,0.00053208};
    _Ty d[7]={ 1.25331414,0.23498619,-0.0365562,
           0.01504268,-0.00780353,0.00325614,-0.00068245};
    if(n<0) n=-n;
    if(x<0.0) x=-x;
    if(FloatEqual(x,0.0)) return(1.0e+70);
    if(n!=1)
    { 
		if(x < 2.0 || FloatEqual(x, 2.0))
		{ 
			y=x*x/4.0; p=a[6];
            for(i=5; i>=0; i--) p=p*y+a[i];
            p=p-TransformativeIntegerBessel1stFunction(0,x)*log(x/2.0);
        }
        else
        { 
			y=2.0/x; p=c[6];
            for(i=5; i>=0; i--) p=p*y+c[i];
            p=p*exp(-x)/sqrt(x);
        }
    }
    if(n==0) return(p);
    b0=p;
    if(x < 2.0 || FloatEqual(x, 2.0))
    { 
		y=x*x/4.0;
		p=b[6];
        for(i=5; i>=0; i--) p=p*y+b[i];
        p=p/x+TransformativeIntegerBessel1stFunction(1,x)*log(x/2.0);
    }
    else
    { 
		y=2.0/x;
		p=d[6];
        for(i=5; i>=0; i--) p=p*y+d[i];
        p=p*exp(-x)/sqrt(x);
    }
    if(n==1) return(p);
    b1=p;
    y=2.0/x;
    for(i=1; i<n; i++)
    { 
		p=b0+i*y*b1;
		b0=b1; 
		b1=p;
	}
    return(p);
}

//不完全贝塔函数
template <class _Ty>
_Ty IncompleteBetaFunction(_Ty a, _Ty b, _Ty x)
{
	_Ty y;

    if(a<0.0 || FloatEqual(a,0.0))
    { 
		cout << "err**a<=0!";
		return(-1.0);
	}
    if(b<0.0 || FloatEqual(b,0.0))
    { 
		cout << "err**b<=0!";
		return(-1.0);
	}
    if((x<0.0)||(x>1.0))
    {
		cout << "err**x<0 or x>1 !";
        return(1.0e+70);
    }
    if(FloatEqual(x,0.0)||FloatEqual(x,1.0)) y=0.0;
    else
    { 
		y=a*log(x)+b*log(1.0-x);
        y=exp(y);
        y=y*GammaFunction(a+b)/(GammaFunction(a)*GammaFunction(b));
    }
    if(x<(a+1.0)/(a+b+2.0))
      y=y*beta(a,b,x)/a;
    else
      y=1.0-y*beta(b,a,1.0-x)/b;
    return(y);
}

//不完全贝塔函数的从属函数
template <class _Ty>
_Ty beta(_Ty a, _Ty b, _Ty x)
{
    _Ty p0(0.0), q0(1.0), p1(1.0), q1(1.0), s1;

    for(int k=1; k<101; k++)
    {
		_Ty d=(a+k)*(a+b+k)*x;
		d=-d/((a+k+k)*(a+k+k+1.0));
        p0=p1+d*p0;
		q0=q1+d*q0;
		_Ty s0=p0/q0;
        d=k*(b-k)*x;
        d=d/((a+k+k-1.0)*(a+k+k));
        p1=p0+d*p1;
		q1=q0+d*q1;
		s1=p1/q1;
        if(Abs(s1-s0)<Abs(s1)*1.0e-07)	return(s1);
    }
    cout << "a or b too big !";
    return(s1);
}

//正态分布函数
template <class _Ty>
_Ty NormalDistributionFunction(_Ty a, _Ty d, _Ty x)
{
	_Ty y;
    
    if(d<0.0 || FloatEqual(d,0.0)) d=1.0e-10;
    y=0.5+0.5*ErrorFunction((x-a)/(sqrt(2.0)*d));
    return(y);
}

//t-分布函数
template <class _Ty>
_Ty tDistributionFunction(_Ty t, int n)
{
	_Ty y;
    
    if(t<0.0) t = -t;
    y=1.0-IncompleteBetaFunction(n/2.0, 0.5, n/(n+t*t));
    return(y);
}

//X^2-分布函数
template <class _Ty>
_Ty X2DistributionFunction(_Ty x, int n)
{
	_Ty y;
    
    if(x<0.0) x=-x;
    y=IncompleteGammaFunction(n/2.0,x/2.0);
    return(y);
}

//F-分布函数
template <class _Ty>
_Ty FDistributionFunction(_Ty f, int n1, int n2)
{
	_Ty y;
    //extern double lbeta();
    if(f<0.0) f=-f;
    y=IncompleteBetaFunction(n2/2.0,n1/2.0,n2/(n2+n1*f));
    return(y);
}

//正弦积分
template <class _Ty>
_Ty SineIntegralFunction(_Ty x)
{
	int m(1), i, j;
    _Ty aa,bb,w,xx;
    _Ty t[5]={-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459};
    _Ty c[5]={0.2369268851,0.4786286705,0.5688888889,
											0.4786286705,0.2369268851};

    if(FloatEqual(x,0.0)) return(0.0);
    _Ty h=Abs(x); 
	_Ty s=Abs(0.0001*h);
    _Ty p(1.0e+35), ep(FLOATERROR), g(0.0);
    while((ep>0.0000001 || FloatEqual(ep,0.0000001)) && (Abs(h)>s))
    { 
		g=0.0;
        for(i=1;i<=m;i++)
        { 
			aa=(i-1.0)*h; 
			bb=i*h;
            w=0.0;
            for(j=0;j<5;j++)
            {
				xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                w=w+sin(xx)/xx*c[j];
            }
            g=g+w;
        }
        g=g*h/2.0;
        ep=Abs(g-p)/(1.0+Abs(g));
        p=g;
		m=m+1;
		h=Abs(x)/m;
    }
    return(g);
}

//余弦积分
template <class _Ty>
_Ty CosineIntegralFunction(_Ty x)
{
	int m(1), i, j;
    _Ty ep(FLOATERROR),aa,bb,w,xx,g(0.0);
    _Ty t[5]={-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459};
    _Ty c[5]={0.2369268851,0.4786286705,0.5688888889,
											0.4786286705,0.2369268851};
    
    if(FloatEqual(x,0.0)) x = 1.0e-35;
    if(x<0.0) x = -x;
    _Ty r = 0.57721566490153286060651;
    _Ty q = r + log(x);
    _Ty h = x; 
	_Ty s = Abs(0.0001*h);
    _Ty p=1.0e+35;
    while((ep>0.0000001||FloatEqual(ep,0.0000001))&&(Abs(h)>s))
    { 
		g=0.0;
        for(i=1;i<=m;i++)
        { 
			aa=(i-1.0)*h; 
			bb=i*h;
            w=0.0;
            for(j=0;j<5;j++)
            {
				xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                w=w+(1.0-cos(xx))/xx*c[j];
            }
            g=g+w;
        }
        g=g*h/2.0;
        ep=Abs(g-p)/(1.0+Abs(g));
        p=g; m=m+1; h=x/m;
    }
    g=q-g;
    return(g);
}

//指数积分
template <class _Ty>
_Ty ExponentIntegralFunction(_Ty x)
{
	int m(1), i, j;
    _Ty ep(FLOATERROR),aa,bb,w,xx,g(0.0);
    _Ty t[5]={-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459};
    _Ty c[5]={0.2369268851,0.4786286705,0.5688888889,
											0.4786286705,0.2369268851};
    
    if(FloatEqual(x,0.0)) x=1.0e-10;
    if(x<0.0) x=-x;
    _Ty r=0.57721566490153286060651;
    _Ty q=r+log(x);
    _Ty h=x;
	_Ty s=Abs(0.0001*h);
    _Ty p=1.0e+35;
    while((ep>0.0000001||FloatEqual(ep,0.0000001))&&(Abs(h)>s))
    { 
		g=0.0;
        for(i=1;i<=m;i++)
        {
			aa=(i-1.0)*h; bb=i*h;
            w=0.0;
            for(j=0;j<5;j++)
            { 
				xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                w=w+(exp(-xx)-1.0)/xx*c[j];
            }
            g=g+w;
        }
        g=g*h/2.0;
        ep=Abs(g-p)/(1.0+Abs(g));
        p=g; m=m+1; h=x/m;
    }
    g=q+g;
    return(g);
}

//第一类椭圆积分
template <class _Ty>
_Ty Ellipse1stIntegral(_Ty k, _Ty f)
{
	int n(0);
    _Ty pi(3.1415926), e(1.0), ff;
    //_Ty fk();
    if(k<0.0) k=-k;
    if(k>1.0) k=1.0/k;
    _Ty y = Abs(f);
    
    while(y>pi||FloatEqual(y,pi))
    {
		n=n+1;
		y=y-pi;
	}
    
    if(y>pi/2.0||FloatEqual(y,pi/2.0))
    {
		n=n+1;
		e=-e; 
		y=pi-y;
	}
    if(n==0)
		ff=fk(k,y);
    else
    { 
		ff=fk(k,pi/2.0);
        ff=2.0*n*ff+e*fk(k,y);
    }
    if(f<0.0) ff=-ff;
    return(ff);
}

//第一类椭圆积分辅助函数
template <class _Ty>
_Ty fk(_Ty k, _Ty f)
{ 
	int m(1), i, j;
    _Ty p(1.0e+35),ep(FLOATERROR),aa,bb,w,xx,g(0.0),q;
    _Ty t[5]={-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459};
    _Ty c[5]={0.2369268851,0.4786286705,0.5688888889,
											0.4786286705,0.2369268851};
    
    _Ty h=Abs(f); 
	_Ty s=Abs(0.0001*h);
    //p=1.0e+35; ep=0.000001;
    while((ep>0.0000001||FloatEqual(ep,0.0000001))&&(Abs(h)>s))
    { 
		g=0.0;
        for(i=1;i<=m;i++)
        {
			aa=(i-1.0)*h; 
			bb=i*h;
            w=0.0;
            for(j=0;j<5;j++)
            { 
				xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                q=sqrt(1.0-k*k*sin(xx)*sin(xx));
                w=w+c[j]/q;
            }
            g=g+w;
        }
        g=g*h/2.0;
        ep=Abs(g-p)/(1.0+Abs(g));
        p=g; 
		m=m+m; 
		h=0.5*h;
    }
    return(g);
}

//第二类椭圆积分
template <class _Ty>
_Ty Ellipse2ndIntegral(_Ty k, _Ty f)
{
	int n(0);
    _Ty pi(3.1415926), e(1.0), ff;
    //_Ty ek();
    if(k<0.0) k=-k;
    if(k>1.0) k=1.0/k;
    _Ty y = Abs(f);
    
    while(y>pi||FloatEqual(y,pi))
    {
		n=n+1; 
		y=y-pi;
	}
    
    if(y>pi/2.0||FloatEqual(y,pi/2.0))
    {
		n=n+1;
		e=-e; 
		y=pi-y;
	}
    if(n==0)
		ff=ek(k,y);
    else
    { 
		ff=ek(k,pi/2.0);
        ff=2.0*n*ff+e*ek(k,y);
    }
    if(f<0.0) ff=-ff;
    return(ff);
}

//第二类椭圆积分辅助函数
template <class _Ty>
_Ty ek(_Ty k, _Ty f)
{	int m(1), i, j;
    _Ty p(1.0e+35),ep(FLOATERROR),aa,bb,w,xx,g(0.0),q;
    _Ty t[5]={-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459};
    _Ty c[5]={0.2369268851,0.4786286705,0.5688888889,
											0.4786286705,0.2369268851};
    
    _Ty h=Abs(f);
	_Ty s=Abs(0.0001*h);
    
    while((ep>0.0000001||FloatEqual(ep,0.0000001))&&(Abs(h)>s))
    { 
		g=0.0;
        for(i=1;i<=m;i++)
		{
			aa=(i-1.0)*h;
			bb=i*h;
            w=0.0;
            for(j=0;j<5;j++)
            { 
				xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                q=sqrt(1.0-k*k*sin(xx)*sin(xx));
                w=w+q*c[j];
            }
            g=g+w;
        }
        g=g*h/2.0;
        ep=Abs(g-p)/(1.0+Abs(g));
        p=g; 
		m=m+m; 
		h=0.5*h;
    }
    return(g);
}

#endif		// _SPECIALFUNCTION_INL