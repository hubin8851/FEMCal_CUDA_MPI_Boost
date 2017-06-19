//FittingApproximation.inl		拟合与逼近头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _FITTINGAPPROXIMATION_INL		//避免多次编译
#define _FITTINGAPPROXIMATION_INL

//最小二乘曲线拟合 
template <class _Ty>
void FitCurveLeastSquares(valarray<_Ty>& x, valarray<_Ty>& y, 
								valarray<_Ty>& a, valarray<_Ty>& dt)
{
	//int i,j,k;
    _Ty g,q,d2,s[20],t[20],b[20];
	int n = x.size();				//给定数据点的个数

	int m = a.size();	//拟合多项式的项数，即拟合多项式的最高次数为m-1
	//要求m<=n且m<=20。若m>n或m>20，本函数自动按m=min{n,20}处理。

    for(int i=0; i<m; i++) a[i]=0.0;
    if(m>n) m = n;
    if(m>20) m=20;
    _Ty z(0), p(0), c(0);;
    for(i=0; i<n; i++) z = z + x[i] / n;
    b[0]=1.0;
	_Ty d1 = n;
    for(i=0; i<n; i++)
    { 
		p=p+(x[i]-z);
		c=c+y[i];}
		c=c/d1; 
		p=p/d1;
		a[0]=c*b[0];
		if(m>1)
		{
			t[1]=1.0;
			t[0]=-p;
			d2=0.0;
			c=0.0; 
			g=0.0;
        for(i=0; i<n; i++)
        { 
			q=x[i]-z-p; d2=d2+q*q;
            c=c+y[i]*q;
            g=g+(x[i]-z)*q*q;
		}
        c=c/d2;
		p=g/d2;
		q=d2/d1;
        d1=d2;
        a[1]=c*t[1]; 
		a[0]=c*t[0]+a[0];
    }
    for(int j=2; j<m; j++)
    { 
		s[j]=t[j-1];
        s[j-1]=-p*t[j-1]+t[j-2];
        if(j>=3)
          for(int k=j-2; k>=1; k--)	s[k]=-p*t[k]+t[k-1]-q*b[k];
        s[0]=-p*t[0]-q*b[0];
        d2=0.0;
		c=0.0;
		g=0.0;
        for(i=0; i<n; i++)
        { 
			q=s[j];
            for(int k=j-1; k>=0; k--)	q=q*(x[i]-z)+s[k];
            d2=d2+q*q; c=c+y[i]*q;
            g=g+(x[i]-z)*q*q;
        }
        c=c/d2; 
		p=g/d2; 
		q=d2/d1;
        d1=d2;
        a[j]=c*s[j]; t[j]=s[j];
        for(int k=j-1; k>=0; k--)
        {
			a[k]=c*s[k]+a[k];
            b[k]=t[k]; t[k]=s[k];
        }
    }
    dt[0]=0.0;
	dt[1]=0.0;
	dt[2]=0.0;
    for(i=0; i<n; i++)
    { 
		q=a[m-1];
        for(int k=m-2; k>=0; k--)	q=a[k]+q*(x[i]-z);
        p=q-y[i];
        if(Abs(p)>dt[2]) dt[2]=fabs(p);
        dt[0]=dt[0]+p*p;
        dt[1]=dt[1]+fabs(p);
    }
}

//切比雪夫曲线拟合
template <class _Ty>
void FitCurveChebyshev(valarray<_Ty>& x, valarray<_Ty>& y, valarray<_Ty>& a)
{
	int ii,k,im,ix[21];
    _Ty h[21],y1,y2,h1,h2,d,hm;
	int n = x.size();	//给定数据点的个数

	int m = a.size()-1;	//拟合多项式的项数，即拟合多项式的最高次数为m-1
	//要求m<n且m<=20。若m>=n或m>20，本函数自动按m=min{n-1,20}处理。

    for(int i=0; i<=m; i++) a[i]=0.0;
    if(m>=n) m=n-1;
    if(m>=20) m=19;
    int m1 = m + 1;
    _Ty ha(0);
    ix[0] = 0; 
	ix[m] = n - 1;
    int l = (n - 1) / m; 
	int j = l;
    for(i=1; i<m; i++)
    {
		ix[i]=j;
		j=j+l;
	}
    while(1)
    {
		_Ty hh=1.0;
        for(i=0; i<=m; i++)
        {
			a[i]=y[ix[i]];
			h[i]=-hh;
			hh=-hh;
		}
        for(int j=1; j<=m; j++)
        {
			ii=m1;
			y2=a[ii-1]; 
			h2=h[ii-1];
            for(i=j; i<=m; i++)
            {
				d=x[ix[ii-1]]-x[ix[m1-i-1]];
                y1=a[m-i+j-1];
                h1=h[m-i+j-1];
                a[ii-1]=(y2-y1)/d;
                h[ii-1]=(h2-h1)/d;
                ii=m-i+j;
				y2=y1;
				h2=h1;
            }
        }
        hh=-a[m]/h[m];
        for(i=0; i<=m; i++)	a[i]=a[i]+h[i]*hh;
        for(j=1; j<m; j++)
        {
			ii=m-j;
			d=x[ix[ii-1]];
            y2=a[ii-1];
            for(k=m1-j; k<=m; k++)
            {
				y1=a[k-1];
				a[ii-1]=y2-d*y1;
                y2=y1; 
				ii=k;
            }
        }
        hm = Abs(hh);
        if(hm < ha || FloatEqual(hm, ha))
		{
			a[m] = -hm;
			return;
		}
        a[m]=hm; 
		ha=hm; 
		im=ix[0]; 
		h1=hh;
        j=0;
        for(i=0; i<n; i++)
        {
			if (i==ix[j])
            { 
				if (j<m) j=j+1;
			}
			else
			{ 
				h2=a[m-1];
				for(k=m-2; k>=0; k--)	h2=h2*x[i]+a[k];
				h2=h2-y[i];
				if(Abs(h2)>hm)
				{
					hm=Abs(h2); 
					h1=h2;
					im=i;
				}
			}
		}
		if(im==ix[0]) return;
		i=0;
		l=1;
		while(l==1)
		{ 
			l=0;
			if(im>=ix[i])
			{
				i=i+1;
				if (i<=m) l=1;
			}
		}
		if(i>m) i=m;
		if(i==(i/2)*2) h2=-hh;
		else h2=hh;
		if(h1*h2>=0.0) ix[i]=im;
		else
		{
			if(im<ix[0])
			{ 
				for(j=m-1; j>=0; j--)	ix[j+1]=ix[j];
				ix[0]=im;
			}
			else
			{
				if(im>ix[m])
				{
					for(j=1; j<=m; j++)	ix[j-1]=ix[j];
					ix[m]=im;
				}
				else 
					ix[i-1]=im;
			}
		}
	}
}

//最佳一致逼近多项式里米兹法
template <class _Ty>
void ApproximationRemez(_Ty a, _Ty b, valarray<_Ty>& p, _Ty eps)
{
    _Ty x[21],g[21],t,s,xx,x0,h,yy;

	int n = p.size()-1;	//n-1次最佳一致逼近多项式的项数
	//要求n<=20; 若n>20，函数自动取n=20

    if(n>20) n=20;
    int m=n+1; 
	_Ty d(1.0e+35);
    for(int k=0; k<=n; k++)
    {
		t = cos((n - k) * 3.1415926 / n);
        x[k] = (b + a + (b-a) * t) / 2.0;
    }
    while(1)
    {
		_Ty u(1);
        for(int i=0; i<m; i++)
        {
			p[i]=FunctionValueAR(x[i]);
            g[i]=-u;
			u=-u;
        }
        for(int j=0; j<n; j++)
        {
			k=m;
			s=p[k-1];
			xx=g[k-1];
            for(i=j; i<n; i++)
            {
				t=p[n-i+j-1];
				x0=g[n-i+j-1];
                p[k-1]=(s-t)/(x[k-1]-x[m-i-2]);
                g[k-1]=(xx-x0)/(x[k-1]-x[m-i-2]);
                k=n-i+j;
				s=t;
				xx=x0;
            }
        }
        u=-p[m-1]/g[m-1];
        for(i=0; i<m; i++)	p[i]=p[i]+g[i]*u;
        for(j=1; j<n; j++)
        {
			k=n-j; 
			h=x[k-1];
			s=p[k-1];
            for(i=m-j; i<=n; i++)
            {
				t=p[i-1];
				p[k-1]=s-h*t;
                s=t; 
				k=i;
            }
        }
        p[m-1]=Abs(u);
		u=p[m-1];
        if(Abs(u-d) <= eps) return;
        d=u; 
		h=0.1*(b-a)/(1.0*n);
        xx=a; 
		x0=a;
        while(x0<=b)
        {
			s=FunctionValueAR(x0);
			t=p[n-1];
            for(i=n-2; i>=0; i--)	t=t*x0+p[i];
            s=Abs(s-t);
            if(s>u) 
			{
				u=s; 
				xx=x0;
			}
            x0=x0+h;
        }
        s=FunctionValueAR(xx); 
		t=p[n-1];
        for(i=n-2; i>=0; i--)	t=t*xx+p[i];
        yy=s-t; 
		i=1;
		j=n+1;
        while((j-i)!=1)
        {
			k=(i+j)/2;
            if(xx<x[k-1]) j=k;
            else i=k;
        }
        if(xx<x[0])
        {
			s=FunctionValueAR(x[0]); 
			t=p[n-1];
            for(k=n-2; k>=0; k--)	t=t*x[0]+p[k];
            s=s-t;
            if(s*yy>0.0) x[0]=xx;
            else
            {
				for(k=n-1; k>=0; k--)	x[k+1]=x[k];
                x[0]=xx;
            }
        }
        else
        {
			if(xx>x[n])
            {
				s=FunctionValueAR(x[n]); 
				t=p[n-1];
                for(k=n-2; k>=0; k--)	t=t*x[n]+p[k];
                s=s-t;
                if(s*yy>0.0) x[n]=xx;
                else
                {
					for(k=0; k<n; k++)	x[k]=x[k+1];
                    x[n]=xx;
                }
            }
            else
            {
				i=i-1;
				j=j-1;
                s=FunctionValueAR(x[i]);
				t=p[n-1];
                for(k=n-2; k>=0; k--)	t=t*x[i]+p[k];
                s=s-t;
                if(s*yy>0.0) x[i]=xx;
                else x[j]=xx;
            }
        }
    }
}

//矩形域的最小二乘曲面拟合
template <class _Ty>
void FitSurfaceLeastSquares(valarray<_Ty>& x, valarray<_Ty>& y, 
				matrix<_Ty>& z, matrix<_Ty>& a,	valarray<_Ty>& dt)
{
	int l,kk;
    _Ty apx[20],apy[20],bx[20],by[20];
    _Ty t[20],t1[20],t2[20],d2,g,g1,g2;
    _Ty x2,dd,y1,x1;

	int n = x.size();			//给定数据点的X坐标个数
	int m = y.size();			//给定数据点的Y坐标个数

	matrix<_Ty> v(20,m), u(20,20);

	int p = a.GetRowNum();		//拟合多项式中变量x的最高次数加1
	//并要求p<=n且p<=20; 否则在本函数中自动取p=min{n,20}处理

	int q = a.GetColNum();		//拟合多项式中变量y的最高次数加1
	//并要求q<=m且q<=20; 否则在本函数中自动取q=min{m,20}处理

    for(int i=0; i<p; i++)
        for(int j=0; j<q; j++) a(i,j)=0.0;

    if(p>n) p=n;
    if(p>20) p=20;
    if(q>m) q=m;
    if(q>20) q=20;
    _Ty xx(0), yy(0), d1(n); 
    for(i=0; i<n; i++)	xx = xx + x[i] / n;
    for(i=0; i<m; i++)  yy = yy + y[i] / m;
    apx[0]=0.0;
    for(i=0; i<n; i++)	apx[0] = apx[0] + x[i] - xx;
    apx[0] = apx[0] / d1;
    for(int j=0; j<m; j++)	//?????
    { 
		v(0,j)=0.0;
        for(i=0; i<n; i++)	v(0,j)=v(0,j)+z(i,j);
        v(0,j)=v(0,j)/d1;
    }
    if(p>1)
    { 
		d2=0.0;
		apx[1]=0.0;
        for(i=0; i<n; i++)
        {
			g=x[i]-xx-apx[0];
            d2=d2+g*g;
            apx[1]=apx[1]+(x[i]-xx)*g*g;
        }
        apx[1]=apx[1]/d2;
        bx[1]=d2/d1;
        for(j=0; j<m; j++)
        { 
			v(1,j)=0.0;
            for(i=0; i<n; i++)
            {
				g=x[i]-xx-apx[0];
                v(1,j)=v(1,j)+z(i,j)*g;
            }
            v(1,j)=v(1,j)/d2;
        }
        d1=d2;
    }
    for(int k=2; k<p; k++)
    { 
		d2=0.0; 
		apx[k]=0.0;
        for(j=0; j<m; j++) v(k,j)=0.0;
        for(i=0; i<n; i++)
        {
			g1=1.0; g2=x[i]-xx-apx[0];
            for(j=2; j<=k; j++)
            {
				g=(x[i]-xx-apx[j-1])*g2-bx[j-1]*g1;
                g1=g2; g2=g;
            }
            d2=d2+g*g;
            apx[k]=apx[k]+(x[i]-xx)*g*g;
            for(j=0; j<m; j++)	v(k,j)=v(k,j)+z(i,j)*g;
        }
        for(j=0; j<m; j++)	v(k,j)=v(k,j)/d2;
        apx[k]=apx[k]/d2;
        bx[k]=d2/d1;
        d1=d2;
    }
    d1=m; 
	apy[0]=0.0;
    for(i=0; i<m; i++)	apy[0]=apy[0]+y[i]-yy;
    apy[0]=apy[0]/d1;
    for(j=0; j<p; j++)
    { 
		u(j,0)=0.0;
        for(i=0; i<m; i++)	u(j,0)=u(j,0)+v(j,i);
		u(j,0)=u(j,0)/d1;
    }
    if(q>1)
    { 
		d2=0.0; 
		apy[1]=0.0;
        for(i=0; i<m; i++)
        {
			g=y[i]-yy-apy[0];
            d2=d2+g*g;
            apy[1]=apy[1]+(y[i]-yy)*g*g;
        }
        apy[1]=apy[1]/d2;
        by[1]=d2/d1;
        for(j=0; j<p; j++)
		{
			u(j,1)=0.0;
            for(i=0; i<m; i++)
            {
				g=y[i]-yy-apy[0];
				u(j,1)=u(j,1)+v(j,i)*g;
            }
			u(j,1)=u(j,1)/d2;
        }
        d1=d2;
    }
    for(k=2; k<q; k++)
    { 
		d2=0.0;
		apy[k]=0.0;
		for(j=0; j<p; j++) u(j,k)=0.0;
        for(i=0; i<m; i++)
        {
			g1=1.0;
            g2=y[i]-yy-apy[0];
            for(j=2; j<=k; j++)
            {
				g=(y[i]-yy-apy[j-1])*g2-by[j-1]*g1;
                g1=g2;
				g2=g;
            }
            d2=d2+g*g;
            apy[k]=apy[k]+(y[i]-yy)*g*g;
            for(j=0; j<p; j++)	u(j,k)=u(j,k)+v(j,i)*g;
        }
        for(j=0; j<p; j++)	u(j,k)=u(j,k)/d2;
        apy[k]=apy[k]/d2;
        by[k]=d2/d1;
        d1=d2;
    }
    v(0,0)=1.0;
	v(1,0)=-apy[0]; 
	v(1,1)=1.0;
    for(i=0; i<p; i++)
      for(j=0; j<q; j++)	a(i,j)=0.0;
    for(i=2; i<q; i++)
    {
		v(i,i)=v(i-1,i-1);
        v(i,i-1)=-apy[i-1]*v(i-1,i-1)+v(i-1,i-2);
        if(i>=3)
          for(k=i-2; k>=1; k--)
            v(i,k)=-apy[i-1]*v(i-1,k)+v(i-1,k-1)-by[i-1]*v(i-2,k);
        v(i,0)=-apy[i-1]*v(i-1,0)-by[i-1]*v(i-2,0);
    }
    for(i=0; i<p; i++)
    { 
		if(i==0)
		{
			t[0]=1.0;
			t1[0]=1.0;
		}
        else
        {
			if(i==1)
            {
				t[0]=-apx[0];
				t[1]=1.0;
                t2[0]=t[0];
				t2[1]=t[1];
            }
            else
            { 
				t[i]=t2[i-1];
                t[i-1]=-apx[i-1]*t2[i-1]+t2[i-2];
                if(i>=3)
                  for(k=i-2; k>=1; k--)
                    t[k]=-apx[i-1]*t2[k]+t2[k-1]
                         -bx[i-1]*t1[k];
                t[0]=-apx[i-1]*t2[0]-bx[i-1]*t1[0];
                t2[i]=t[i];
                for(k=i-1; k>=0; k--)
                { 
					t1[k]=t2[k];
					t2[k]=t[k];
				}
            }
        }
        for(j=0; j<q; j++)
          for(k=i; k>=0; k--)
            for(l=j; l>=0; l--)
				a(k,l)=a(k,l)+u(i,j)*t[k]*v(j,l);
    }
    dt[0]=0.0; 
	dt[1]=0.0;
	dt[2]=0.0;
    for(i=0; i<n; i++)
    {
		x1=x[i]-xx;
        for(j=0; j<m; j++)
        {
			y1=y[j]-yy;
            x2=1.0; 
			dd=0.0;
            for(k=0; k<p; k++)
            {
				g=a(k,q-1);
                for(kk=q-2; kk>=0; kk--)	g=g*y1+a(k,kk);
                g=g*x2;
				dd=dd+g;
				x2=x2*x1;
            }
            dd=dd-z(i,j);
            if(Abs(dd)>dt[2]) dt[2]=Abs(dd);
            dt[0]=dt[0]+dd*dd;
            dt[1]=dt[1]+Abs(dd);
        }
    }
}

#endif		// _FITTINGAPPROXIMATION_INL

