//OrdinaryDifferentialEguation.inl		求解常微分方程(组)头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _ORDINARYDIFFERENTIALEGUATION_INL	//避免多次编译
#define _ORDINARYDIFFERENTIALEGUATION_INL

//定步长全区间积分的欧拉法
template <class _Ty>
void ODE_EulerContentStep(_Ty t, valarray<_Ty>& y,  
										_Ty h, int k, matrix<_Ty>& z)
{
	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> d(n);
    for(int i=0; i<n; i++)
		z(i,0)=y[i];
    for(int j=1; j<k; j++)
	{ 
		_Ty x=t+(j-1)*h;
        FunctionValueECS(x,y,d);
        for(i=0; i<n; i++)
          y[i]=z(i,j-1)+h*d[i];
        x=t+j*h;
        FunctionValueECS(x,y,d);
        for(i=0; i<n; i++)
          d[i]=z(i,j-1)+h*d[i];
        for(i=0; i<n; i++)
		{ 
			y[i]=(y[i]+d[i])/2.0;
            z(i,j)=y[i];
		}
	}
}

//变步长积分欧拉法
template <class _Ty>
void ODE_EulerVariationalStep(_Ty t, _Ty h, valarray<_Ty>& y, _Ty eps)
{
	int j, m(1);
    _Ty p(1.0+eps), x, q, hh=h;

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
    valarray<_Ty> a(n), b(n), c(n), d(n);

    for(int i=0; i<n; i++) a[i] = y[i];

    while(p>eps || FloatEqual(p,eps))
	{ 
		for(i=0; i<n; i++)
		{
			b[i]=y[i];
			y[i]=a[i];
		}
        for(j=0; j<m; j++)
		{
			for(i=0; i<n; i++) c[i]=y[i];
            x=t+j*hh;
            FunctionValueEVS(x,y,d);
            for(i=0; i<n; i++)
				y[i]=c[i]+hh*d[i];
            x=t+(j+1)*hh;
            FunctionValueEVS(x,y,d);
            for(i=0; i<n; i++)
				d[i]=c[i]+hh*d[i];
            for(i=0; i<n; i++)
				y[i]=(y[i]+d[i])/2.0;
		}
        p = 0.0;
        for(i=0; i<n; i++)
		{ 
			q=Abs(y[i]-b[i]);
            if(q>p) p=q;
        }
        hh /= 2.0;
		m = m + m;
    }
}

//全区间定步长积分维梯法
template <class _Ty>
void ODE_WittyContentStep(_Ty t, valarray<_Ty>& y,  
										_Ty h, int k, matrix<_Ty>& z)
{
    int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> a(n), d(n);

    for(int i=0; i<n; i++) z(i,0)=y[i];
    FunctionValueWCS(t,y,d);
    for(int j=1; j<k; j++)
    { 
		for(i=0; i<n; i++)
			a[i]=z(i,j-1)+h*d[i]/2.0;	
        _Ty x=t+(j-0.5)*h;
        FunctionValueWCS(x,a,y);
        for(i=0; i<n; i++)
        { 
			d[i]=2.0*y[i]-d[i];
            z(i,j)=z(i,j-1)+h*y[i];
        }
    }
}

//全区间定步长积分龙格-库塔法
template <class _Ty>
void ODE_RungeKuttaContentStep(_Ty t, valarray<_Ty>& y,  
										_Ty h, int k, matrix<_Ty>& z)
{
    _Ty a[4], tt;

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> b(n), d(n);

    a[0] = h / 2.0;
	a[1] = a[0];
    a[2] = a[3] = h;

    for(int i=0; i<n; i++) z(i,0)=y[i];

    for(int l=1; l<k; l++)
    { 
		FunctionValueRKCS(t, y, d);
        for(i=0; i<n; i++) b[i]=y[i];
        for(int j=0; j<=2; j++)
        {
			for(i=0; i<n; i++)
            {
				y[i]=z(i,l-1)+a[j]*d[i];
                b[i]=b[i]+a[j+1]*d[i]/3.0;
            }
            tt=t+a[j];
            FunctionValueRKCS(tt, y, d);
        }
        for(i=0; i<n; i++)
          y[i]=b[i]+h*d[i]/6.0;
        for(i=0; i<n; i++)
          z(i,l)=y[i];
        t=t+h;
    }
}

//变步长积分龙格-库塔法
template <class _Ty>
void ODE_RungeKuttaVariationalStep(_Ty t, _Ty h, 
										valarray<_Ty>& y, _Ty eps)
{
	int m(1), j, k;
    _Ty hh(h), dt,x(t),tt,q,a[4];

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> g(n), b(n), c(n), d(n), e(n);
    _Ty p=1.0+eps;

    for(int i=0; i<n; i++) c[i]=y[i];
    while(p>eps || FloatEqual(p,eps))
    {
		a[0]=hh/2.0; 
		a[1]=a[0];
		a[2] = a[3] = hh;
        for(i=0; i<n; i++)
        {
			g[i]=y[i];
			y[i]=c[i];
		}
        dt=h/m; 
		t=x;
        for(j=0; j<m; j++)
        {
			FunctionValueRKVS(t, y, d);
            for(i=0; i<n; i++) 
            {
				b[i]=y[i]; 
				e[i]=y[i];
			}
            for(k=0; k<=2; k++)
            {
				for(i=0; i<n; i++)
				{
					y[i]=e[i]+a[k]*d[i];
                    b[i]=b[i]+a[k+1]*d[i]/3.0;
                }
                tt=t+a[k];
                FunctionValueRKVS(tt, y, d);
            }
            for(i=0; i<n; i++)
              y[i]=b[i]+hh*d[i]/6.0;
            t += dt;
        }
        p = 0.0;
        for(i=0; i<n; i++)
        {
			q = Abs(y[i]-g[i]);
            if(q>p) p = q;
        }
        hh /= 2.0;
		m = m + m;
    }
}

//变步长积分基尔法
template <class _Ty>
void ODE_GillVariationalStep(_Ty t, _Ty h, valarray<_Ty>& y,
										_Ty eps, valarray<_Ty>& q)
{
    int i,k,m(1),ii;
    _Ty x(t),hh(h),r,s,t0,dt,qq;
    _Ty a[4]={0.5, 0.29289321881, 1.7071067812, 0.166666667};
    _Ty b[4]={2.0,1.0,1.0,2.0};
    _Ty c[4],e[4]={0.5,0.5,1.0,1.0};

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> g(n), u(n), v(n), d(n);

    for(i=0; i<3; i++) c[i]=a[i];
    c[3]=0.5;
    _Ty p=1.0+eps;
	
    for(int j=0; j<n; j++) u[j]=y[j];
    while(p>eps || FloatEqual(p,eps))
    {
		for(j=0; j<n; j++)
		{
			v[j]=y[j];
			y[j]=u[j];
			g[j]=q[j];
		}
        dt=h/m; 
		t=x;
        for(k=0; k<m; k++)
        { 
			FunctionValueGVS(t,y,d);
            for(ii=0; ii<=3; ii++)
            {
				for(j=0; j<n; j++)
                  d[j]=d[j]*hh;
                for(j=0; j<n; j++)
                {
					r=(a[ii]*(d[j]-b[ii]*g[j])+y[j])-y[j];
                    y[j]=y[j]+r;
                    s=g[j]+3.0*r;
                    g[j]=s-c[ii]*d[j];
                }
                t0=t+e[ii]*hh;
                FunctionValueGVS(t0,y,d);
            }
            t=t+dt;
        }
        p=0.0;
        for(j=0; j<n; j++)
        {
			qq=Abs(y[j]-v[j]);
            if(qq>p) p=qq;
        }
        hh=hh/2.0; 
		m=m+m;
    }
    for(j=0; j<n; j++) q[j]=g[j];
}

//全区间变步长积分基尔法
template <class _Ty>
void ODE_GillVariationalStepInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z)
{
    int j,m,kk,ii;
    _Ty aa(t),hh,x,p,dt,r,s,t0,qq;
    _Ty a[4]={0.5, 0.29289321881, 1.7071067812, 0.166666667};
    _Ty b[4]={2.0,1.0,1.0,2.0};
    _Ty c[4],e[4]={0.5,0.5,1.0,1.0};

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> g(n), q(n), u(n), v(n), d(n);
    
    for(int i=0; i<3; i++) c[i]=a[i];
    c[3]=0.5;
    
    for(i=0; i<n; i++) 
    {
		z(i,0) = y[i];
		q[i] = 0.0;
	}
    for(i=2; i<=k; i++)
    {
		x=aa+(i-2)*h;
		m=1;
		hh=h;
        p=1.0+eps;
        for(j=0; j<n; j++) u[j]=y[j];
        while(p>eps || FloatEqual(p,eps))
        {
			for(j=0; j<n; j++)
			{ 
				v[j]=y[j];
				y[j]=u[j];
				g[j]=q[j];
			}
            dt=h/m; 
			t=x;
            for(kk=0; kk<m; kk++)
            {
				FunctionValueGVSI(t,y,d);
                for(ii=0; ii<4; ii++)
                {
					for(j=0; j<n; j++)
                      d[j]=d[j]*hh;
                    for(j=0; j<n; j++)
                    {
						r=(a[ii]*(d[j]-b[ii]*g[j])+y[j])-y[j];
                        y[j]=y[j]+r;
                        s=g[j]+3.0*r;
                        g[j]=s-c[ii]*d[j];
                    }
                    t0=t+e[ii]*hh;
                    FunctionValueGVSI(t0,y,d);
                }
                t += dt;
            }
            p=0.0;
            for(j=0; j<n; j++)
            { 
				qq=Abs(y[j]-v[j]);
                if(qq>p) p=qq;
            }
            hh /= 2.0; 
			m=m+m;
        }
        for(j=0; j<n; j++)
        { 
			q[j]=g[j];
			z(j,i-1)=y[j];
		}
    }
}

//变步长积分默森法
template <class _Ty>
void ODE_MersonVariationalStep(_Ty t, _Ty h, int n, valarray<_Ty>& y,
										_Ty eps, int k, matrix<_Ty>& z)
{
    int j,m,nn;
    _Ty aa(t),bb,x,hh,p,dt,t0,qq;
    valarray<_Ty> a(n), b(n), c(n), u(n), v(n), d(n);
        
    for(int i=0; i<n; i++) z(i,0)=y[i];
    for(i=1; i<k; i++)
    { 
		x=aa+(i-1)*h; 
		nn=1; 
		hh=h;
        for(j=0; j<n; j++) u[j]=y[j];
        p=1.0+eps;
        while(p>eps || FloatEqual(p,eps))
        {
			for(j=0; j<n; j++)
			{
				v[j]=y[j]; 
				y[j]=u[j];
			}
            dt=h/nn; t=x;
            for(m=0; m<nn; m++)
			{ 
				FunctionValueMVS(t,y,d);
                for(j=0; j<n; j++)
                {
					a[j]=d[j];
					y[j]=y[j]+hh*d[j]/3.0;
				}
                t0=t+hh/3.0;
                FunctionValueMVS(t0,y,d);
                for(j=0; j<n; j++)
                {
					b[j]=d[j];
					y[j]=y[j]+hh*(d[j]-a[j])/6.0;
				}
                FunctionValueMVS(t0,y,d);
                for(j=0; j<n; j++)
                {
					b[j]=d[j];
                    bb=(d[j]-4.0*(b[j]+a[j]/4.0)/9.0)/8.0;
                    y[j]=y[j]+3.0*hh*bb;
                }
                t0=t+hh/2.0;
                FunctionValueMVS(t0,y,d);
                for(j=0; j<n; j++)
                { c[j]=d[j];
                    qq=d[j]-15.0*(b[j]-a[j]/5.0)/16.0;
                    y[j]=y[j]+2.0*hh*qq;
                }
                t0=t+hh;
                FunctionValueMVS(t0,y,d);
                for(j=0; j<n; j++)
                { 
					qq=c[j]-9.0*(b[j]-2.0*a[j]/9.0)/8.0;
                    qq=d[j]-8.0*qq;
                    y[j]=y[j]+hh*qq/6.0;
                }
                t += dt;
            }
            p=0.0;
            for(j=0; j<n; j++)
            { qq=Abs(y[j]-v[j]);
                if(qq>p) p=qq;
            }
            hh /= 2.0;
			nn=nn+nn;
        }
        for(j=0; j<n; j++) z(j,i)=y[j];
    }
}

//积分一步连分式法
template <class _Ty>
void ODE_Fraction(_Ty t, _Ty h, valarray<_Ty>& y, _Ty eps)
{
	int i,j,k(1),m,nn(1),it(1);

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
    _Ty x(t),hh(h),dd,q,p,g[10];
	valarray<_Ty> e(n), b(10*n), w(n), u(n), v(n), d(n);

    for(j=0; j<n; j++) v[j]=y[j];
	g[0]=hh;
    Fraction1(x,hh,y,w,d,e);
    for(j=0; j<n; j++)
    { 
		b[j]=y[j];
		u[j]=y[j];
	}

    while(it==1)
    { 
		nn=nn+nn; 
		hh=hh/2.0;
		it=0;
        g[k]=hh;
        for(j=0; j<n; j++) y[j]=v[j];
        t=x;
        for(j=0; j<nn; j++)
        { 
			Fraction1(t,hh,y,w,d,e);
            t=t+hh;
        }
        for(j=0; j<n; j++)
        { 
			dd=y[j]; 
			m=0;
            for(i=0; i<k; i++)
              if(m==0)
              {
				  q=dd-b[i*n+j];
                  if(FloatEqual(q,0.0)) m=1;
                  else dd=(g[k]-g[i])/q;
              }
            b[k*n+j]=dd;
            if(m!=0) b[k*n+j]=1.0e+35;
        }
        for(j=0; j<n; j++)
        {
			dd=0.0;
            for(i=k-1; i>=0; i--)
				dd=-g[i]/(b[(i+1)*n+j]+dd);
            y[j]=dd+b[j];
        }
        p=0.0;
        for(j=0; j<n; j++)
        { 
			q=Abs(y[j]-u[j]);
            if(q>p) p=q;
        }
        if((p>eps || FloatEqual(p,eps))&&(k<7))
        { 
			for(j=0; j<n; j++) u[j]=y[j];
            k=k+1; 
			it=1;
        }
    }
}

//积分一步连分式法ODE_Fraction()辅助函数
template <class _Ty>
void Fraction1(_Ty t, _Ty h, valarray<_Ty>& y, valarray<_Ty>& b, 
									valarray<_Ty>& d, valarray<_Ty>& e)
{
    _Ty a[4],tt;
	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
    
	a[0]=h/2.0;
	a[1]=a[0];
	a[2]=h; 
	a[3]=h;
    FunctionValueF(t,y,d);
    for(int i=0; i<n; i++)
	{
		b[i]=y[i];
		e[i]=y[i];
	}
    for(int k=0; k<=2; k++)
    {
		for(i=0; i<n; i++)
        {
			y[i]=e[i]+a[k]*d[i];
            b[i]=b[i]+a[k+1]*d[i]/3.0;
        }
        tt=t+a[k];
        FunctionValueF(tt,y,d);
    }
    for(i=0; i<n; i++)
      y[i]=b[i]+h*d[i]/6.0;
}


//全区间积分连分式法
template <class _Ty>
void ODE_FractionInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z)
{
	int j,kk,l,m,nn,it;
    _Ty x,hh,dd,tt,q,p,g[10];

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> e(n), b(10*n), w(n), u(n), v(n), d(n);

    for(int i=0; i<n; i++) z(i,0)=y[i];
    for(i=1; i<k; i++)
    { 
		for(j=0; j<n; j++) v[j]=y[j];
        x=t+(i-1)*h; 
		nn=1; 
		hh=h; 
		g[0]=hh;
        Fraction2(x,hh,y,w,d,e);
        for(j=0; j<n; j++)
        {
			b[j]=y[j];
			u[j]=y[j];
		}
        kk=1;
		it=1;
        while(it==1)
        {
			nn=nn+nn; 
			hh=hh/2.0;
			it=0;
            g[kk]=hh;
            for(j=0; j<n; j++) y[j]=v[j];
            tt=x;
            for(j=0; j<nn; j++)
            {
				Fraction2(t,hh,y,w,d,e);
                tt=tt+hh;
            }
            for(j=0; j<n; j++)
            { 
				dd=y[j]; 
				l=0;
                for(m=0; m<kk; m++)
                  if(l==0)
                  {
					  q=dd-b[m*n+j];
                      if(Abs(q)+1.0==1.0) l=1;
                      else dd=(g[kk]-g[m])/q;
                  }
                b[kk*n+j]=dd;
                if(l!=0) b[kk*n+j]=1.0e+35;
            }
            for(j=0; j<n; j++)
            { 
				dd=0.0;
                for(m=kk-1; m>=0; m--)
					dd=-g[m]/(b[(m+1)*n+j]+dd);
                y[j]=dd+b[j];
            }
            p=0.0;
            for(j=0; j<n; j++)
            { 
				q=Abs(y[j]-u[j]);
                if(q>p) p=q;
            }
            if((p>eps || FloatEqual(p,eps))&&(kk<7))
            { for(j=0; j<n; j++) u[j]=y[j];
                kk=kk+1; it=1;
            }
        }
        for(j=0; j<n; j++)
          z(j,i)=y[j];
    }
}

//全区间积分连分式法辅助函数
template <class _Ty>
void Fraction2(_Ty t, _Ty h, valarray<_Ty>& y, valarray<_Ty>& b, 
									valarray<_Ty>& d, valarray<_Ty>& e)
{
    int i,k;
    _Ty a[4],tt;

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
    a[0]=h/2.0; 
	a[1]=a[0];
	a[2]=h;
	a[3]=h;
    FunctionValueFI(t,y,d);
    for(i=0; i<n; i++) 
	{
		b[i]=y[i]; 
		e[i]=y[i];
	}
    for(k=0; k<3; k++)
    { 
		for(i=0; i<n; i++)
        {
			y[i]=e[i]+a[k]*d[i];
            b[i]=b[i]+a[k+1]*d[i]/3.0;
        }
        tt=t+a[k];
        FunctionValueFI(tt,y,d);
    }
    for(i=0; i<n; i++)
      y[i]=b[i]+h*d[i]/6.0;
}

//全区间积分双边法
template <class _Ty>
void ODE_BothSidesInterval(_Ty t, _Ty h, valarray<_Ty>& y,
					   				_Ty eps, int k, matrix<_Ty>& z)
{
    _Ty a(t), qq;

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> d(n), w(n), u(n), v(n), p(n);

    for(int i=0; i<n; i++)
	{
		p[i]=0.0;
		z(i,0)=y[i];
	}

    FunctionValueBSI(t,y,d);
    for(int j=0; j<n; j++) u[j]=d[j];
    BothSidesInterval(t,h,y,eps);		//辅助函数
    t=a+h;
    FunctionValueBSI(t,y,d);
    for(j=0; j<n; j++)
    { 
		z(j,1)=y[j]; v[j]=d[j];}
    for(j=0; j<n; j++)
    { 
		p[j]=-4.0*z(j,1)+5.0*z(j,0)+2.0*h*(2.0*v[j]+u[j]);
        y[j]=p[j];
    }
    t=a+2.0*h;
    FunctionValueBSI(t,y,d);
    for(j=0; j<n; j++)
    { 
		qq=2.0*h*(d[j]-2.0*v[j]-2.0*u[j])/3.0;
        qq=qq+4.0*z(j,1)-3.0*z(j,0);
        z(j,2)=(p[j]+qq)/2.0;
        y[j]=z(j,2);
    }
    for(i=3; i<k; i++)
    { 
		t=a+(i-1)*h;
        FunctionValueBSI(t,y,d);
        for(j=0; j<n; j++)
        {
			u[j]=v[j];
			v[j]=d[j];
		}
        for(j=0; j<n; j++)
        { 
			qq=-4.0*z(j,i-1)+5.0*z(j,i-2);
            p[j]=qq+2.0*h*(2.0*v[j]+u[j]);
            y[j]=p[j];
        }
        t=t+h;
        FunctionValueBSI(t,y,d);
        for(j=0; j<n; j++)
        { 
			qq=2.0*h*(d[j]-2.0*v[j]-2.0*u[j])/3.0;
            qq=qq+4.0*z(j,i-1)-3.0*z(j,i-2);
            y[j]=(p[j]+qq)/2.0;
            z(j,i)=y[j];
        }
    }
}

//全区间积分双边法辅助函数
template <class _Ty>
void BothSidesInterval(_Ty t, _Ty h, valarray<_Ty>& y, _Ty eps)
{ 
    int m(1), j, k;
    _Ty hh(h), dt,x(t), tt, q, a[4];

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> d(n), b(n), c(n), e(n), g(n);
    
    _Ty p=1.0+eps;
    for(int i=0; i<n; i++) c[i]=y[i];
    while(p>eps || FloatEqual(p,eps))
    { 
		a[0]=hh/2.0;
		a[1]=a[0]; 
		a[2]=hh;
		a[3]=hh;
        for(i=0; i<n; i++)
        {
			g[i]=y[i]; 
			y[i]=c[i];
		}
        dt=h/m; 
		t=x;
        for(j=0; j<m; j++)
        {
			FunctionValueBSI(t,y,d);
            for(i=0; i<n; i++) 
            {
				b[i]=y[i];
				e[i]=y[i];
			}
            for(k=0; k<3; k++)
            { 
				for(i=0; i<n; i++)
                {
					y[i]=e[i]+a[k]*d[i];
                    b[i]=b[i]+a[k+1]*d[i]/3.0;
                }
                tt=t+a[k];
                FunctionValueBSI(tt,y,d);
            }
            for(i=0; i<n; i++)
              y[i]=b[i]+hh*d[i]/6.0;
            t=t+dt;
        }
        p=0.0;
        for(i=0; i<n; i++)
        { 
			q=Abs(y[i]-g[i]);
            if(q>p) p=q;
        }
        hh=hh/2.0;
		m=m+m;
    }
}

//全区间积分阿当姆斯预报-校正法
template <class _Ty>
void ODE_AdamsInterval(_Ty t, _Ty h, valarray<_Ty>& y,
									_Ty eps, int k, matrix<_Ty>& z)
{
    int j,m;
    _Ty a(t), q;

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> b(4*n), e(n), s(n), g(n), d(n);
    
    for(int i=0; i<n; i++) z(i,0)=y[i];
    FunctionValueAI(t,y,d);
    for(i=0; i<n; i++) b[i]=d[i];
    for(i=1; i<4; i++)
      if(i<k)
      { 
		  t=a+i*h;
          AdamsInterval(t,h,y,eps);			//辅助函数
          for(j=0; j<n; j++) z(j,i)=y[j];
          FunctionValueAI(t,y,d);
          for(j=0; j<n; j++) b[i*n+j]=d[j];
      }
    for(i=4; i<k; i++)
    { 
		for(j=0; j<n; j++)
        { 
			q=55.0*b[3*n+j]-59.0*b[2*n+j];
            q=q+37.0*b[n+j]-9.0*b[j];
            y[j]=z(j,i-1)+h*q/24.0;
            b[j]=b[n+j];
            b[n+j]=b[n+n+j];
            b[n+n+j]=b[n+n+n+j];
        }
        t=a+i*h;
        FunctionValueAI(t,y,d);
        for(m=0; m<n; m++) b[n+n+n+m]=d[m];
        for(j=0; j<n; j++)
        { 
			q=9.0*b[3*n+j]+19.0*b[n+n+j]-5.0*b[n+j]+b[j];
            y[j]=z(j,i-1)+h*q/24.0;
            z(j,i)=y[j];
        }
        FunctionValueAI(t,y,d);
        for(m=0; m<n; m++) b[3*n+m]=d[m];
    }
}

//全区间积分阿当姆斯预报-校正法辅助函数
template <class _Ty>
void AdamsInterval(_Ty t, _Ty h, valarray<_Ty>& y, _Ty eps)
{
    int m(1), j, k;

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
    _Ty hh(h), dt, x(t), tt, q, a[4];
	valarray<_Ty> g(n), b(n), c(n), d(n), e(n);
    
    _Ty p=1.0+eps;
    for(int i=0; i<n; i++) c[i]=y[i];
    while(p>eps || FloatEqual(p,eps))
    { 
		a[0]=hh/2.0;
		a[1]=a[0]; 
		a[2]=hh;
		a[3]=hh;
        for(i=0; i<n; i++)
        {
			g[i]=y[i]; 
			y[i]=c[i];
		}
        dt=h/m; t=x;
        for(j=0; j<m; j++)
        { 
			FunctionValueAI(t,y,d);
            for(i=0; i<n; i++) 
            {
				b[i]=y[i];
				e[i]=y[i];
			}
            for(k=0; k<3; k++)
            { 
				for(i=0; i<n; i++)
                {
					y[i]=e[i]+a[k]*d[i];
                    b[i]=b[i]+a[k+1]*d[i]/3.0;
                }
                tt=t+a[k];
                FunctionValueAI(tt,y,d);
            }
            for(i=0; i<n; i++)
              y[i]=b[i]+hh*d[i]/6.0;
            t=t+dt;
        }
        p=0.0;
        for(i=0; i<n; i++)
        {
			q=Abs(y[i]-g[i]);
            if(q>p) p=q;
        }
        hh=hh/2.0;
		m=m+m;
    }
}

//全区间积分哈明法
template <class _Ty>
void ODE_HammingInterval(_Ty t, _Ty h, valarray<_Ty>& y,
										_Ty eps, int k, matrix<_Ty>& z)
{
    int j,m;
    _Ty a(t), q;

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> d(n), w(n), u(n), v(n), b(4*n), g(n);
    
    for(int i=0; i<n; i++) z(i,0)=y[i];
    FunctionValueHI(t,y,d);
    for(i=0; i<n; i++) b[i]=d[i];
    for(i=1; i<=3; i++)
      if(i<k)
      { 
		  t=a+i*h;
          HammingInterval(t,h,y,eps);
          for(m=0; m<n; m++) z(m,i)=y[m];
          FunctionValueHI(t,y,d);
          for(m=0; m<n; m++) b[i*n+m]=d[m];
      }
    for(i=0; i<n; i++) u[i]=0.0;
    for(i=4; i<k; i++)
    { 
		for(j=0; j<n; j++)
        {
			q=2.0*b[3*n+j]-b[n+n+j]+2.0*b[n+j];
            y[j]=z(j,i-4)+4.0*h*q/3.0;
        }
        for(j=0; j<n; j++)
			y[j]=y[j]+112.0*u[j]/121.0;
        t=a+i*h;
        FunctionValueHI(t,y,d);
        for(j=0; j<n; j++)
        { 
			q=9.0*z(j,i-1)-z(j,i-3);
            q=(q+3.0*h*(d[j]+2.0*b[3*n+j]-b[n+n+j]))/8.0;
            u[j]=q-y[j];
            z(j,i)=q-9.0*u[j]/121.0;
            y[j]=z(j,i);
            b[n+j]=b[n+n+j];
            b[n+n+j]=b[n+n+n+j];
        }
        FunctionValueHI(t,y,d);
        for(m=0; m<n; m++) b[3*n+m]=d[m];
    }
}

//全区间积分哈明法辅助函数
template <class _Ty>
void HammingInterval(_Ty t, _Ty h, valarray<_Ty>& y, _Ty eps)
{ 
    int m(1), j, k;
    _Ty hh(h), dt, x(t), tt, q, a[4];

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> g(n), b(n), c(n), d(n), e(n);
    
	_Ty p=1.0+eps;
    for(int i=0; i<n; i++) c[i]=y[i];
    while(p>eps || FloatEqual(p,eps))
    { 
		a[0]=hh/2.0; 
		a[1]=a[0]; 
		a[2]=hh;
		a[3]=hh;
        for(i=0; i<n; i++)
        {
			g[i]=y[i]; 
			y[i]=c[i];
		}
        dt=h/m;
		t=x;
        for(j=0; j<m; j++)
        {
			FunctionValueHI(t,y,d);
            for(i=0; i<n; i++) 
            {
				b[i]=y[i];
				e[i]=y[i];
			}
            for(k=0; k<3; k++)
            { 
				for(i=0; i<n; i++)
                {
					y[i]=e[i]+a[k]*d[i];
                    b[i]=b[i]+a[k+1]*d[i]/3.0;
                }
                tt=t+a[k];
                FunctionValueHI(t,y,d);
            }
            for(i=0; i<n; i++)
				y[i]=b[i]+hh*d[i]/6.0;
            t=t+dt;
        }
        p=0.0;
        for(i=0; i<n; i++)
        { 
			q=Abs(y[i]-g[i]);
            if(q>p) p=q;
        }
        hh=hh/2.0; 
		m=m+m;
    }
}

//积分一步特雷纳法
template <class _Ty>
void ODE_Treanor(_Ty t, _Ty h, valarray<_Ty>& y)
{
//extern void gtnr1f();
    _Ty s,aa,bb,dd,g,dy,dy1;

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> w(4*n), q(4*n), r(4*n), d(n), p(n);
    
    for(int j=0; j<n; j++) w[j]=y[j];
    FunctionValueT(t,y,d);
    for(j=0; j<n; j++)
    {
		q[j]=d[j]; 
		y[j]=w[j]+h*d[j]/2.0;
        w[n+j]=y[j];
    }
    s=t+h/2.0;
    FunctionValueT(t,y,d);
    for(j=0; j<n; j++)
    { 
		q[n+j]=d[j];
        y[j]=w[j]+h*d[j]/2.0;
        w[n+n+j]=y[j];
    }
    FunctionValueT(t,y,d);
    for(j=0; j<n; j++) q[n+n+j]=d[j];

    for(j=0; j<n; j++)
    { 
		aa=q[n+n+j]-q[n+j];
        bb=w[n+n+j]-w[n+j];
        if(-aa*bb*h>0.0)
        {
			p[j]=-aa/bb; dd=-p[j]*h;
            r[j]=exp(dd);
            r[n+j]=(r[j]-1.0)/dd;
            r[n+n+j]=(r[n+j]-1.0)/dd;
            r[3*n+j]=(r[n+n+j]-1.0)/dd;
        }
        else p[j]=0.0;
        if(p[j]<0.0||FloatEqual(p[j],0.0)) g=q[n+n+j];
        else
        { 
			g=2.0*(q[n+n+j]-q[j])*r[n+n+j];
            g=g+(q[j]-q[n+j])*r[n+j]+q[n+j];
        }
        w[3*n+j]=w[j]+g*h;
        y[j]=w[3*n+j];
    }
    s=t+h;
    FunctionValueT(t,y,d);
    for(j=0; j<n; j++) q[3*n+j]=d[j];

    for(j=0; j<n; j++)
    { 
		if(p[j]<0.0||FloatEqual(p[j],0.0))
        { 
			dy=q[j]+2.0*(q[n+j]+q[n+n+j]);
            dy=(dy+q[n+n+n+j])*h/6.0;
        }
        else
        {
			dy=-3.0*(q[j]+p[j]*w[j])+2.0*(q[n+j]+p[j]*w[n+j]);
            dy=dy+2.0*(q[n+n+j]+p[j]*w[n+n+j]);
            dy=dy-(q[n+n+n+j]+p[j]*w[n+n+n+j]);
            dy=dy*r[n+n+j]+q[j]*r[n+j];
            dy1=q[j]-q[n+j]-q[n+n+j]+q[n+n+n+j];
            dy1=dy1+(w[j]-w[n+j]-w[n+n+j]+w[n+n+n+j])*p[j];
            dy=(dy+4.0*dy1*r[n+n+n+j])*h;
        }
        y[j]=w[j]+dy;
    }
}

//全区间积分特雷纳法
template <class _Ty>
void ODE_TreanorInterval(_Ty t, _Ty h, valarray<_Ty>& y,
											int k, matrix<_Ty>& z)
{
    _Ty a(t),s,aa,bb,dd,g,dy,dy1;

	int n = y.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> w(4*n), q(4*n), r(4*n), d(n), p(n);
    
    for(int j=0; j<n; j++) z(j,0)=y[j];

    for(int i=1; i<k; i++)
    {
		t=a+(i-1)*h;
        for(j=0; j<n; j++) w[j]=y[j];
        FunctionValueTI(s,y,d);
        for(j=0; j<n; j++)
        { 
			q[j]=d[j]; 
			y[j]=w[j]+h*d[j]/2.0;
            w[n+j]=y[j];
        }
        s=t+h/2.0;
        FunctionValueTI(s,y,d);
        for(j=0; j<n; j++)
        { 
			q[n+j]=d[j];
            y[j]=w[j]+h*d[j]/2.0;
            w[n+n+j]=y[j];
        }
        FunctionValueTI(s,y,d);
        for(j=0; j<n; j++) q[n+n+j]=d[j];

        for(j=0; j<n; j++)
        { 
			aa=q[n+n+j]-q[n+j];
            bb=w[n+n+j]-w[n+j];
            if(-aa*bb*h>0.0)
            {
				p[j]=-aa/bb; dd=-p[j]*h;
                r[j]=exp(dd);
                r[n+j]=(r[j]-1.0)/dd;
                r[n+n+j]=(r[n+j]-1.0)/dd;
                r[3*n+j]=(r[n+n+j]-1.0)/dd;
            }
            else p[j]=0.0;
            if(p[j]<0.0||FloatEqual(p[j],0.0)) g=q[n+n+j];
            else
            { 
				g=2.0*(q[n+n+j]-q[j])*r[n+n+j];
                g=g+(q[j]-q[n+j])*r[n+j]+q[n+j];
            }
            w[3*n+j]=w[j]+g*h;
            y[j]=w[3*n+j];
        }
        s=t+h;
        FunctionValueTI(s,y,d);
        for(j=0; j<n; j++) q[3*n+j]=d[j];
        for(j=0; j<n; j++)
        { 
			if(p[j]<0.0||FloatEqual(p[j],0.0))
            { 
				dy=q[j]+2.0*(q[n+j]+q[n+n+j]);
                dy=(dy+q[n+n+n+j])*h/6.0;
            }
            else
            { 
				dy=-3.0*(q[j]+p[j]*w[j])+2.0*(q[n+j]+p[j]*w[n+j]);
                dy=dy+2.0*(q[n+n+j]+p[j]*w[n+n+j]);
                dy=dy-(q[n+n+n+j]+p[j]*w[n+n+n+j]);
                dy=dy*r[n+n+j]+q[j]*r[n+j];
                dy1=q[j]-q[n+j]-q[n+n+j]+q[n+n+n+j];
                dy1=dy1+(w[j]-w[n+j]-w[n+n+j]+w[n+n+n+j])*p[j];
                dy=(dy+4.0*dy1*r[n+n+n+j])*h;
            }
            y[j]=w[j]+dy;
            z(j,i)=y[j];
        }
    }
}

//积分刚性方程组吉尔法
template <class _Ty>
int ODE_Gear(_Ty a, _Ty b, _Ty hmin, _Ty hmax, _Ty h, _Ty eps, 
		valarray<_Ty>& y0, int k, valarray<_Ty>& t, matrix<_Ty>& z)
{
    int kf(1),jt(0),nn(0),nq(1),m(2),irt1,irt(1),j,nqd,idb;
    int iw,j1,j2,nt,nqw,l;
    _Ty aa[7],hw(h),hd,rm,t0(a),td,r,dd,pr1,pr2,pr3,rr;
    _Ty enq1,enq2,enq3,eup,e,edwn,bnd,r1;

    _Ty pp[7][3] = 
	{
		{2.0,3.0,1.0}, {4.5,6.0,1.0},
        {7.333,9.167,0.5}, {10.42,12.5,0.1667},
        {13.7,15.98,0.04133}, {17.15,1.0,0.008267},
        {1.0,1.0,1.0}
	};
	
	int n = y0.size();	//微分方程组中方程的个数，也是未知函数的个数
	valarray<_Ty> d(n), s(10*n), s02(n), ym(n), er(n);
	valarray<_Ty> yy(n), y(8*n), iis(n), jjs(n);
	matrix<_Ty> p(n,n);

    aa[1]=-1.0;
    for(int i=0; i<8*n; i++) y[i]=0.0;
    for(i=0; i<n; i++) 
    {
		y[i*8]=y0[i];
		yy[i]=y[i*8];
	}
    FunctionValueG(t0,yy,d);
    for(i=0; i<n; i++) y[i*8+1]=h*d[i];
    for(i=0; i<n; i++) ym[i]=1.0;
  l20:
    irt=1;
	kf=1; 
	nn=nn+1;
    t[nn-1]=t0;
    for(i=0; i<n; i++) z(i,nn-1)=y[i*8];
    if((t0>b||FloatEqual(t0,b))||(nn==k))	return(kf);
    
	for(i=0; i<n; i++)
      for(j=0; j<m; j++) s[i*10+j]=y[i*8+j];
    hd=hw;
    if(h!=hd)
    {
		rm=h/hd;
		irt1=0;
        rr=Abs(hmin/hd);
        if(rm<rr) rm=rr;
        rr=Abs(hmax/hd);
        if(rm>rr) rm=rr;
		r=1.0; 
		irt1=irt1+1;
        for(j=1; j<m; j++)
        { 
			r=r*rm;
            for(i=0; i<n; i++)
              y[i*8+j]=s[i*10+j]*r;
        }
        h=hd*rm;
        for(i=0; i<n; i++)
			y[i*8]=s[i*10];
        idb=m;
    }
    nqd=nq;
	td=t0;
	rm=1.0;
    if(jt>0) goto l80;
  l60:
    switch(nq)
    { 
		case 1: aa[0]=-1.0;
				break;
        case 2: aa[0]=-2.0/3.0;
				aa[2]=-1.0/3.0; 
				break;
        case 3: aa[0]=-6.0/11.0; 
				aa[2]=aa[0];
                aa[3]=-1.0/11.0; 
				break;
        case 4: aa[0]=-0.48; 
				aa[2]=-0.7;
				aa[3]=-0.2;
                aa[4]=-0.02; 
				break;
        case 5: aa[0]=-120.0/274.0;
				aa[2]=-225.0/274.0;
                aa[3]=-85.0/274.0; 
				aa[4]=-15.0/274.0;
                aa[5]=-1.0/274.0; 
				break;
        case 6: aa[0]=-720.0/1764.0; 
				aa[2]=-1624.0/1764.0;
                aa[3]=-735.0/1764.0;
				aa[4]=-175.0/1764.0;
                aa[5]=-21.0/1764.0; 
				aa[6]=-1.0/1764.0;
                break;
        default:	return(-2);
    }
    m=nq+1; 
	idb=m;
    enq2=0.5/(nq+1.0);
	enq3=0.5/(nq+2.0);
    enq1=0.5/(nq+0.0);
    eup=pp[nq-1][1]*eps;
	eup=eup*eup;
    e=pp[nq-1][0]*eps;
	e=e*e;
    edwn=pp[nq-1][2]*eps;
	edwn=edwn*edwn;
    if(edwn==0.0)
    {
		for(i=0; i<n; i++)
          for(j=0; j<m; j++)
            y[i*8+j]=s[i*10+j];
        h=hd;
		nq=nqd; 
		jt=nq;
		return(-4);
    }
    bnd=eps*enq3/(n+0.0);
    iw=1;
    if(irt==2)
    { 
		r1=1.0;
        for(j=1; j<m; j++)
        {
			r1=r1*r;
            for(i=0; i<n; i++)
				y[i*8+j]=y[i*8+j]*r1;
        }
        idb=m;
        for(i=0; i<n; i++)
          if(ym[i]<Abs(y[i*8]))
				ym[i]=Abs(y[i*8]);
        jt=nq;
        goto l20;
    }
  l80:
    t0=t0+h;
    for(j=2; j<=m; j++)
      for(j1=j; j1<=m; j1++)
      { 
		  j2=m-j1+j-1;
          for(i=0; i<n; i++)
				y[i*8+j2-1]=y[i*8+j2-1]+y[i*8+j2];
      }
    for(i=0; i<n; i++) er[i]=0.0;
    j1=1; 
	nt=1;
    for(l=0; l<=2; l++)
    {
		if((j1!=0)&&(nt!=0))
        {
			for(i=0; i<n; i++) yy[i]=y[i*8];
            FunctionValueG(t0,yy,d);
            if(iw>=1)
            {
				for(i=0; i<n; i++) yy[i]=y[i*8];
                JacobiMatrix(t0,yy,p);			//计算雅可比矩阵的函数
                r=aa[0]*h;
                for(i=0; i<n; i++)
                  for(j=0; j<n; j++)
                    p(i,j)=p(i,j)*r;
                for(i=0; i<n; i++)
                  p(i,i)=1.0+p(i,i);
                iw=-1;
                jjs[0]=MatrixInversionGS(p);	//全选主元高斯-约当法求矩阵逆
                j1=jjs[0];
            }
            if(jjs[0]!=0)
            {
				for(i=0; i<n; i++)
					s02[i]=y[i*8+1]-d[i]*h;
                for(i=0; i<n; i++)
                { 
					dd=0.0;
                    for(j=0; j<n; j++)
						dd=dd+s02[j]*p(i,j);
                    s[i*10+8]=dd;
                }
                nt=n;
                for(i=0; i<n; i++)
                { 
					y[i*8]=y[i*8]+aa[0]*s[i*10+8];
                    y[i*8+1]=y[i*8+1]-s[i*10+8];
                    er[i]=er[i]+s[i*10+8];
                    if(Abs(s[i*10+8])<=(bnd*ym[i]))
						nt=nt-1;
                }
            }
        }
    }
    if(nt>0)
    {
		t0=td;
        if((h>(hmin*1.00001))||(iw>=0))
        {
			if(iw!=0) rm=0.25*rm;
            iw=1;
			irt1=2;
            rr=Abs(hmin/hd);
            if(rm<rr) rm=rr;
            rr=Abs(hmax/hd);
            if(rm>rr) rm=rr;
            r=1.0;
            for(j=1; j<m; j++)
            {
				r=r*rm;
                for(i=0; i<n; i++)
					y[i*8+j]=s[i*10+j]*r;
            }
            h=hd*rm;
            for(i=0; i<n; i++)
              y[i*8]=s[i*10];
            idb=m;
            goto l80;
        }
        for(i=0; i<n; i++)
          for(j=0; j<m; j++)
            y[i*8+j]=s[i*10+j];
        h=hd; 
		nq=nqd;
		jt=nq;
		return(-3);
    }
    dd=0.0;
    for(i=0; i<n; i++)
		dd=dd+(er[i]/ym[i])*(er[i]/ym[i]);
    iw=0;
    if(dd<=e)
    { 
		if(m>=3)
          for(j=2; j<m; j++)
            for(i=0; i<n; i++)
              y[i*8+j]=y[i*8+j]+aa[j]*er[i];
        kf=1;
		hw=h;
        if(idb>1)
        {
			idb=idb-1;
            if(idb<=1)
              for(i=0; i<n; i++)
                s[i*10+9]=er[i];
            for(i=0; i<n; i++)
	      if(ym[i]<Abs(y[i*8])) ym[i]=Abs(y[i*8]);
            jt=nq;
            goto l20;
        }
    }
    if(dd>e)
    { 
		kf=kf-2;
        if(h<=(hmin*1.00001))
        {
            hw=h;
			jt=nq; 
			return(-1);
        }
        t0=td;
        if(kf<=-5)
        { 
			if(nq==1)
            {
				for(i=0; i<n; i++)
                  for(j=0; j<m; j++)
                    y[i*8+j]=s[i*10+j];
                h=hd;
				nq=nqd;
				jt=nq;
				return(-4);
            }
            for(i=0; i<n; i++) yy[i]=y[i*8];
            FunctionValueG(t0,yy,d);
            r=h/hd;
            for(i=0; i<n; i++)
            { 
				y[i*8]=s[i*10];
                s[i*10+1]=hd*d[i];
                y[i*8+1]=s[i*10+1]*r;
            }
            nq=1;
			kf=1;
			goto l60;
        }
    }
    pr2=log(dd/e);
	pr2=enq2*pr2;
	pr2=exp(pr2);
    pr2=1.2*pr2;
    pr3=1.0e+20;
    if(nq<7)
      if(kf>-1)
      {
		  dd=0.0;
          for(i=0; i<n; i++)
          {
			  pr3=(er[i]-s[i*10+9])/ym[i];
              dd=dd+pr3*pr3;
          }
          pr3=log(dd/eup);
		  pr3=enq3*pr3;
          pr3=exp(pr3);
		  pr3=1.4*pr3;
      }
    pr1=1.0e+20;
    if(nq>1)
    { 
		dd=0.0;
        for(i=0; i<n; i++)
        {
			pr1=y[i*8+m-1]/ym[i];
            dd=dd+pr1*pr1;
        }
        pr1=log(dd/edwn); 
		pr1=enq1*pr1;
        pr1=exp(pr1);
		pr1=1.3*pr1;
    }
    if(pr2<=pr3)
    { 
		if(pr2>pr1)
		{
			r=1.0e+04;
            if(pr1>1.0e-04) r=1.0/pr1;
            nqw=nq-1;
        }
        else
        { 
			nqw=nq;
			r=1.0e+04;
            if(pr2>1.0e-04) r=1.0/pr2;
        }
    }
    else
    { 
		if(pr3<pr1)
        { 
			r=1.0e+04;
            if(pr3>1.0e-04) r=1.0/pr3;
            nqw=nq+1;
        }
        else
        { 
			r=1.0e+04;
            if(pr1>1.0e-04) r=1.0/pr1;
            nqw=nq-1;
        }
    }
    idb=10;
    if(kf==1)
      if(r<1.1)
      { 
		  for(i=0; i<n; i++)
            if(ym[i]<Abs(y[i*8])) ym[i]=Abs(y[i*8]);
          jt=nq; goto l20;
      }
    if(nqw>nq)
      for(i=0; i<n; i++)
        y[i*8+nqw]=er[i]*aa[m-1]/(m+0.0);
    m=nqw+1;
    if(kf==1)
    { 
		irt=2; 
		rr=hmax/Abs(h);
        if(r>rr) r=rr;
        h=h*r; 
		hw=h;
        if(nq==nqw)
        {
			r1=1.0;
            for(j=1; j<m; j++)
            {
				r1=r1*r;
                for(i=0; i<n; i++)
                  y[i*8+j]=y[i*8+j]*r1;
            }
            idb=m;
            for(i=0; i<n; i++)
              if(ym[i]<Abs(y[i*8])) ym[i]=Abs(y[i*8]);
            jt=nq;
			goto l20;
        }
        nq=nqw;
        goto l60;
    }
    rm=rm*r;
	irt1=3;
    rr=Abs(hmin/hd);
    if(rm<rr) rm=rr;
    rr=Abs(hmax/hd);
    if(rm>rr) rm=rr;
    r=1.0;
    for(j=1; j<m; j++)
    {
		r=r*rm;
        for(i=0; i<n; i++)
          y[i*8+j]=s[i*10+j]*r;
    }
    h=hd*rm;
    for(i=0; i<n; i++)
      y[i*8]=s[i*10];
    idb=m;
    if(nqw==nq) goto l80;
    nq=nqw; goto l60;
}

//二阶微分方程边值问题数值解法
template <class _Ty>
void ODE_LinearBoundaryValude(_Ty a, _Ty b, _Ty ya, _Ty yb, int n, 
													valarray<_Ty>& y)
{
    int i,k,nn,m1;
    _Ty h,x;
	valarray<_Ty> g(6*n), d(2*n), z(4);

    h=(b-a)/(n-1.0);
	nn=2*n-1;
    g[0]=1.0;
	g[1]=0.0;
    y[0]=ya; 
	y[n-1]=yb;
    g[3*n-3]=1.0;
	g[3*n-4]=0.0;
    for(i=2; i<n; i++)
    { 
		x=a+(i-1)*h;
        FunctionValueLBV(x,z);
        k=3*(i-1)-1;
        g[k]=z[0]-h*z[1]/2.0;
        g[k+1]=h*h*z[2]-2.0*z[0];
        g[k+2]=z[0]+h*z[1]/2.0;
        y[i-1]=h*h*z[3];
    }
    m1=3*n-2;
	
	valarray<_Ty> gg(m1);
	for(i=0;i<m1;i++)	gg[i] = g[i];
	LE_TridiagonalEquationGauss(gg,y);	//求解三对角方程组

	for(i=0;i<m1;i++)	g[i] = gg[i];

    h=h/2.0;
    g[0]=1.0; 
	g[1]=0.0;
    d[0]=ya; 
	d[nn-1]=yb;
    g[3*nn-3]=1.0;
	g[3*nn-4]=0.0;
    for(i=2; i<nn; i++)
    { 
		x=a+(i-1)*h;
        FunctionValueLBV(x,z);
        k=3*(i-1)-1;
        g[k]=z[0]-h*z[1]/2.0;
        g[k+1]=h*h*z[2]-2.0*z[0];
        g[k+2]=z[0]+h*z[1]/2.0;
        d[i-1]=h*h*z[3];
    }
    m1=3*nn-2;

	valarray<_Ty> ggg(m1), ddd(nn);
	for(i=0;i<m1;i++)	ggg[i] = g[i];
	for(i=0;i<nn;i++)	ddd[i] = d[i];
	LE_TridiagonalEquationGauss(ggg,ddd);	//求解三对角方程组
	
    for(i=0;i<m1;i++)	g[i] = ggg[i];
	for(i=0;i<nn;i++)	d[i] = ddd[i];

	for(i=2; i<n; i++)
    { 
		k=i+i-1;
        y[i-1]=(4.0*d[k-1]-y[i-1])/3.0;
    }
}

#endif		// _ORDINARYDIFFERENTIALEGUATION_INL

