// Extremum.inl		计算极值函数(方法)定义头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _EXTREMUM_INL		//避免多次编译
#define _EXTREMUM_INL

//一维极值连分式法
template <class _Ty>
void ExtremumFraction1D(valarray<_Ty>& x, _Ty eps, 
										int k, valarray<int>& js)
{	
    _Ty xx,h1,dx,y[10],b[10],h2(0);
	valarray<_Ty> z(2);
    int jt = 1;
	js[0] = 0;
    while(jt == 1)
    {
		int j=0;
        while(j<=7)
        {
			if(j <=2) 
				xx=x[0]+0.01*j;
            else xx=h2;
            FunctionValueEF1D(xx,z);
            if(Abs(z[1]) < eps)
            {
				jt=0; 
				j=10;
			}
            else
            {
				h1=z[1]; 
				h2=xx;
                if(j==0)
				{
					y[0]=h1;
					b[0]=h2;
				}
                else
                {
					y[j]=h1;
					int m=0;
					int i=0;
                    while((m==0)&&(i<=j-1))
					{
						if(FloatEqual((h2-b[i]),0)) 
							m=1;
						else
							h2=(h1-y[i])/(h2-b[i]);
                        i=i+1;
					}
                    b[j]=h2;
                    if(m!=0)
						b[j]=1.0e+35;
                    h2=0.0;
                    for(i=j-1; i>=0; i--)
						h2=-y[i]/(b[i+1]+h2);
                    h2=h2+b[0];
                }
                j=j+1;
            }
        }
        x[0]=h2;
        js[0]=js[0]+1;
        if(js[0]==k)  jt=0;
    }
    xx = x[0];
    FunctionValueEF1D(xx,z);
	x[1] = z[0];
    if(Abs(x[0]) < 1.0 || FloatEqual(Abs(x[0]), 1)) dx=1.0e-05;
    else	dx = Abs(x[0] * 1.0e-05);
    xx = x[0] - dx;
    FunctionValueEF1D(xx,z);
	h1 = z[0];
    xx = x[0] + dx;
    FunctionValueEF1D(xx,z);
	h2 = z[0];
    js[1] = -1;
    if((h1 + h2 - 2.0 * x[1]) < 0.0 || FloatEqual((h1 + h2 - 2.0 * x[1]),0)) js[1]=1;
}

//n维极值连分式法
template <class _Ty>
void ExtremumFractionND(valarray<_Ty>& x, _Ty eps, 
										int k, valarray<int>& js)
{
    int i,j,m,l,jt(1),il;
    _Ty y[10],b[10],p,z,t,h1,h2(0),f,dx;

	int n = x.size()-1;			//自变量个数

    js[0] = 0;
	//jt = 1;
	//h2 = 0;
    while(jt == 1)
    {
		t=0.0; 
		js[0]=js[0]+1;
        for(i=1; i<=n; i++)
        {
			f=FunctionValueEFND(x,i);
            t=t+Abs(f);
        }
        if(t<eps) jt=0;
        else
        {
			for(i = 0; i < n; i++)
              { 
				il=5;
                while(il!=0)
                {
					j=0;
					t=x[i];
					il=il-1;
                    while(j<=7)
                    {
						if(j<=2) z=t+j*0.01;
                        else z=h2;
                        x[i]=z;
                        f=FunctionValueEFND(x,i+1);
                        if(FloatEqual(f,0))
                        {
							j=10;
							il=0;
						}
                        else
                        {
							h1=f;
							h2=z;
                            if(j==0)
                            {
								y[0]=h1;
								b[0]=h2;
							}
                            else
                            {
								y[j]=h1;
								m=0;
								l=0;
                                while((m==0)&&(l<=j-1))
                                {
									p=h2-b[l];
                                    if(FloatEqual(p,0)) m=1;
                                    else h2=(h1-y[l])/p;
                                    l=l+1;
                                }
                                b[j]=h2;
                                if(m!=0) b[j]=1.0e+35;
                                h2=0.0;
                                for(l=j-1; l>=0; l--)
                                h2=-y[l]/(b[l+1]+h2);
                                h2=h2+b[0];
                            }
                            j=j+1;
                        }
                    }
                    x[i]=h2; 
                }
                x[i]=z;
            }
            if(js[0]==k) jt=0;
        }
    }
    js[1]=1;
    f=FunctionValueEFND(x,0); 
	x[n]=f;
    dx=0.00001;
	t=x[0];
    x[0]=t+dx;
	h1=FunctionValueEFND(x,0);
    x[0]=t-dx;
	h2=FunctionValueEFND(x,0);
    x[0]=t;
    t=h1+h2-2.0*f;
    if(t>0.0) js[1]=-1;
    j=1; jt=1;
    while(jt==1)
    {
		j=j+1;
		dx=0.00001;
		jt=0;
        t=x[j-1];
        x[j-1]=t+dx;
		h2=FunctionValueEFND(x,0);
        x[j-1]=t-dx;
		h1=FunctionValueEFND(x,0);
        x[j-1]=t;
		t=h1+h2-2.0*f;
        if((t*js[1]<0.0)&&(j<n)) jt=1;
    }
    if(t*js[1]>0.0) js[1]=0;
}

//线性规划
template <class _Ty>
int ExtremumLinePrograming(matrix<_Ty>& a, valarray<_Ty>& b, 
								valarray<_Ty>& c, valarray<_Ty>& x)
{
    _Ty y;
	
	int m = b.size();		//不等式约束条件的个数
	int n = c.size() - m;	//变量个数

	valarray<int> js(m);
	matrix<_Ty> p(m, m);
	matrix<_Ty> d(m, (m+n));
	for(int i = 0; i < m; i++) js[i] = n + i;
    int mn = m + n;
	_Ty s(0);
    while(1==1)
    {
		for(i = 0; i < m; i++)
          for(int j = 0; j < m; j++)
			  p(i,j) = a(i, js[j]);
        int l=MatrixInversionGS(p);
        if(l==0)
        {
			x[n]=s;
            return(l);
        }

		d = p * a;			//矩阵相乘
        
        for(i=0; i<mn; i++) x[i]=0.0;
        for(i=0; i<m; i++)
        {
			s=0.0;
            for(int j=0; j<m; j++)
              s=s + p(i,j) * b[j];
            x[js[i]]=s;
        }
        int k = -1;
		_Ty dd=1.0e-35;
        for(int j = 0; j < mn; j++)
        {
			_Ty z(0);
            for(i=0; i<m; i++)
              z = z + c[js[i]] * d(i, j);
            z = z - c[j];
            if(z > dd)
			{
				dd = z;
				k = j;
			}
        }
        if(k == -1)
          {
			s=0.0;
            for(j=0; j<n; j++)
				s = s + c[j] * x[j];
			x[n] = s;
            return(1);
        }
        j = -1;
        dd = 1.0e+20;
        for(i=0; i<m; i++)
			if(d(i, k) >= 1.0e-20 )
            { 
			  y=x[js[i]] / d(i, k);
              if(y < dd)
			  {
				  dd = y; 
				  j = i;
			  }
			}
        if(j==-1) 
        {
			x[n]=s;
            return(0);
        }
        js[j]=k;
    }
}

//n维极值单形调优法
template <class _Ty>
int ExtremumSimplexND(_Ty d, _Ty u, _Ty v, valarray<_Ty>& x,  
				_Ty eps, int k, matrix<_Ty>& xx, valarray<_Ty>& f)
{
	int r,g,l;
	_Ty fe,ft,ff;
	
	int n = x.size() - 1;	//变量个数

	valarray<_Ty> xt(n);
	valarray<_Ty> xf(n);
	valarray<_Ty> xe(n);

    int kk=0;
	_Ty nn=1.0*n;
    _Ty fr=sqrt(nn+1.0);
    _Ty fl=d*(fr-1.0)/(1.414*nn);
    _Ty fg=d*(fr+nn-1.0)/(1.414*nn);
    for(int i=0; i<n; i++)
    for(int j=0; j<=n; j++)
       xx(i,j) = 0;
    for(i=1; i<=n; i++)
    for(int j=0; j<n; j++)
       xx(j,i) = fl;
    for(i=1; i<=n; i++)
       xx((i-1),i) = fg;
    for(i=0; i<=n; i++)
    { 
		for(int j=0; j<n; j++)	xt[j]=xx(j,i);
        f[i]=FunctionValueESND(xt,n);
    }
    ft=1.0+eps;
    while((kk<k)&&(ft>eps))
    {
		kk=kk+1;
        fr=f[0]; fl=f[0]; r=0; l=0;
        for(i=1; i<=n; i++)
        {
			if(f[i]>fr)
			{
				r=i;
				fr=f[i];
			}
            if(f[i]<fl)
			{
				l=i;
				fl=f[i];
			}
        }
        g=0;
		fg=f[0];
        int j=0;
        if(r==0)
		{
			g=1;
			fg=f[1];
			j=1;
		}
        for(i=j+1; i<=n; i++)
			if((i!=r)&&(f[i]>fg))
            {
				g=i;
				fg=f[i];
			}
        for(j=0; j<n; j++)
        {
			xf[j]=0.0;
            for(i=0; i<=n; i++)
				if(i!=r)	xf[j]=xf[j]+xx(j,i)/nn;
            xt[j]=2.0*xf[j]-xx(j,r);
        }
        ft=FunctionValueESND(xt,n);
        if(ft<f[l])
        {
			for(j=0; j<n; j++)	xf[j]=(1.0+u)*xt[j]-u*xf[j];
            ff=FunctionValueESND(xf,n);
            if(ff<f[l])
            {
				for(j=0; j<n; j++)	xx(j,r)=xf[j];
                f[r]=ff;
            }
            else
			{
				for(j=0; j<n; j++)	xx(j,r)=xt[j];
                f[r]=ft;
            }
        }
        else
        {
			if(ft<=f[g])
            {
				for(j=0; j<=n-1; j++)	xx(j,r)=xt[j];
                f[r]=ft;
            }
            else 
            {
				if(ft<=f[r])
                {
					for(j=0; j<=n-1; j++)	xx(j,r)=xt[j];
                    f[r]=ft;
                }
                for(j=0; j<n; j++)	xf[j]=v*xx(j,r)+(1.0-v)*xf[j];
                ff=FunctionValueESND(xf,n);
                if(ff>f[r])
                  for(i=0; i<=n; i++)
                  {
					  for(j=0; j<=n-1; j++)
					  {
						  xx(j,i)=(xx(j,i)+xx(j,l))/2.0;
                          x[j]=xx(j,i);
						  xe[j]=x[j];
                      }
                      fe=FunctionValueESND(xe,n);
					  f[i]=fe;
                  }
                else
                {
					for(j=0; j<=n-1; j++)	xx(j,r)=xf[j];
                    f[r]=ff;
                }
            }
        }
        ff=0.0;
		ft=0.0;
        for(i=0; i<=n; i++)
        {
			ff=ff+f[i]/(1.0+nn);
            ft=ft+f[i]*f[i];
        }
        ft=(ft-(1.0+n)*ff*ff)/nn;
    }
    for(int j=0; j<n; j++)
    {
		x[j] = 0.0;
        for(i=0; i<=n; i++)	x[j] = x[j] + xx(j,i) / (1.0+nn);
        xe[j] = x[j];
    }
    fe = FunctionValueESND(xe,n);
	x[n] = fe;
    return(kk);
}

//n维极值复形调优法
template <class _Ty>
int ExtremumComplexND(int m, valarray<_Ty>& a, valarray<_Ty>& b,
		_Ty alpha, _Ty eps, valarray<_Ty>& x, matrix<_Ty>& xx, int k)
{
    int r, g, kt, jt, kk;
    _Ty fj, fr, fg, z;
	double rr = 0;				//随机数种子
	int n = a.size();			//变量个数

	valarray<_Ty> c(m);
	valarray<_Ty> d(m);
	valarray<_Ty> w(m);
	valarray<_Ty> xt(m);
	valarray<_Ty> xf(m);

    for(int i=0; i<n; i++)	xx(i,0) = x[i];
    xx(n,0) = FunctionValueTarget(x,n);
    for(int j=1; j<2*n; j++)
    {
		for(i=0; i<n; i++)
		{
			xx(i,j) = a[i] + (b[i] - a[i]) * rand_01_One(rr);
            x[i] = xx(i,j);
        }
        int it = 1;
        while(it==1)
        {
			it = 0;
			r = 0;
			g = 0;
            while((r<n)&&(g==0))
            {
				if((a[r]<=x[r])&&(b[r]>=x[r])) r = r + 1;
                else g = 1;
            }
            if(g==0)
            {
				ConditionValue(n,m,x,c,d,w);
                r = 0;
                while((r<m)&&(g==0))
                {
					if((c[r]<=w[r])&&(d[r]>=w[r])) r=r+1;
                    else g=1;
                }
            }
            if(g!=0)
            {
				for(r=0; r<n; r++)
                {
					z = 0.0;
                    for(g=0; g<j; g++)	z = z + xx(r,g) / (1.0*j);
                    xx(r,j) = (xx(r,j) + z) / 2.0;
                    x[r]=xx(r,j);
                }
                it = 1;
            }
            else xx(n,j) = FunctionValueTarget(x,n);
        }
    }
    kk = 1;
	int it = 1;
    while(it==1)
    {
		it = 0;
        fr = xx(n,0);
		r = 0;
        for(i=1; i<2*n; i++)
          if(xx(n,i)>fr)
          {
			  r = i;
			  fr = xx(n,i);
		  }
        g = 0;
		j = 0;
		fg = xx(n,0);
        if(r==0)
        {
			g = 1;
			j = 1;
			fg = xx(n,1);
		}
        for(i=j+1; i<2*n; i++)
          if(i!=r)
            if(xx(n,i)>fg)
            {
				g = i;
				fg = xx(n,i);
			}
        for(i=0; i<n; i++)
        {
			xf[i] = 0.0;
            for(j=0; j<2*n; j++)
              if(j!=r)	xf[i]=xf[i]+xx(i,j)/(2.0*n-1.0);
            xt[i]=(1.0+alpha)*xf[i]-alpha*xx(i,r);
        }
        jt = 1;
        while(jt==1)
        {
			jt = 0;
            z = FunctionValueTarget(xt,n);
            while(z>fg)
            {
				for(i=0; i<n; i++)	xt[i]=(xt[i]+xf[i])/2.0;
                z = FunctionValueTarget(xt,n);
            }
            j=0;
            for(i=0; i<n; i++)
            { 
				if(a[i]>xt[i])
                {
					xt[i]=xt[i]+0.000001;
					j=1;
				}
				if(b[i]<xt[i])
                {
					xt[i]=xt[i]-0.000001;
					j=1;
				}
            }
            if(j!=0) jt=1;
            else
            {
				ConditionValue(n,m,xt,c,d,w);
                j=0;
				kt=1;
                while((kt==1)&&(j<m))
                {
					if((c[j]<=w[j])&&(d[j]>=w[j])) j=j+1;
                    else kt=0;
                }
                if(j<m)
                { 
					for(i=0; i<n; i++)	xt[i]=(xt[i]+xf[i])/2.0;
                    jt=1;
                }
            }
        }
        for(i=0; i<n; i++)	xx(i,r)=xt[i];
        xx(n,r) = z;
        fr = 0.0; 
		fg = 0.0;
        for(j=0; j<2*n; j++)
        {
			fj = xx(n,j);
            fr = fr + fj / (2.0*n);
            fg = fg + fj * fj;
        }
        fr = (fg - 2.0 * n * fr * fr) / (2.0 * n - 1.0);
        if(fr>=eps)
        {
			kk = kk + 1;
            if(kk<k) it = 1;
        }
    }
    for(i=0; i<n; i++)
    {
		x[i] = 0.0;
        for(j=0; j<2*n; j++)	x[i]=x[i]+xx(i,j)/(2.0*n);
    }
    z = FunctionValueTarget(x,n);
	x[n] = z;
    return(kk);
}

//确定1维函数极小值点所在区间
template <class _Ty>
void MinimizerInterval(_Ty& ax, _Ty& bx, _Ty& cx, _Ty& fa, _Ty& fb, _Ty& fc)
{
    _Ty r,q,dum,gold = 1.0 + GoldNo;
    _Ty u,ulim,fu,tiny = 1e-20;
	int glimit = 100;

    fa = FunctionValue(ax);
    fb = FunctionValue(bx);
    if (fb > fa)
	{
        dum = ax;
        ax = bx;
        bx = dum;
        dum = fb;
        fb = fa;
        fa = dum;
    }
    cx = bx + gold * (bx - ax);
    fc = FunctionValue(cx);
	while (fb >= fc)
	{
        r = (bx - ax) * (fb - fc);
        q = (bx - cx) * (fb - fa);
        dum = q - r;
        if (Abs(dum) < tiny)
		{
			dum = tiny;
		}
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (2 * dum);
        ulim = bx + glimit * (cx - bx);
        if ((bx - u) * (u - cx) > 0)
		{
            fu = FunctionValue(u);
            if (fu < fc)
			{
                ax = bx;
                fa = fb;
                bx = u;
                fb = fu;
                return;
			}
            else
			{
				if (fu > fb)
				{
					cx = u;
					fc = fu;
					return;
				}
            }
            u = cx + gold * (cx - bx);
            fu = FunctionValue(u);
		}
        else
		{
			if ((cx - u) * (u - ulim) > 0)
			{
				fu = FunctionValue(u);
				if (fu < fc)
				{
					bx = cx;
					cx = u;
					u = cx + gold * (cx - bx);
					fb = fc;
					fc = fu;
					fu = FunctionValue(u);
				}
            }
			else
			{
				if ((u - ulim) * (ulim - cx) > 0 || FloatEqual((u - ulim) * (ulim - cx),0))
				{
					u = ulim;
					fu = FunctionValue(u);
				}
				else
				{
					u = cx + gold * (cx - bx);
					fu = FunctionValue(u);
				}
			}
		}
		ax = bx;
		bx = cx;
		cx = u;
		fa = fb;
		fb = fc;
		fc = fu;
	}
}

//确定一维函数极小值点所在区间(重载)
template <class _Ty>
int MinimizerInterval(_Ty& ax, _Ty& bx, _Ty& cx, int& sign, 
									_Ty b, _Ty& fa, _Ty& fb, _Ty& fc)
{
    _Ty r,q,dum,gold = GoldNo + 1.0;
    _Ty u,ulim,fu,tiny = 1e-20;
	int glimit = 100;

    fa = FunctionValue(ax);
    fb = FunctionValue(bx);
    if (fb > fa)
	{
		sign = -1;
		fa *= sign;
		fb *= sign;
    }
    cx = bx + gold * (bx - ax);
	fc = sign * FunctionValue(cx);
	while (fb >= fc)
	{
        r = (bx - ax) * (fb - fc);
        q = (bx - cx) * (fb - fa);
        dum = q - r;
        if (Abs(dum) < tiny)
		{
			dum = tiny;
		}
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (2 * dum);
        ulim = bx + glimit * (cx - bx);
		if(u > b)	//超出f(x)右边界
		{
			cout<<"  The MinimizerInterval exceed right border. Program is over!"<<endl;
			return 0;
		}
        if ((bx - u) * (u - cx) > 0)
		{
            fu = sign * FunctionValue(u);
            if (fu < fc)
			{
                ax = bx;
                fa = fb;
                bx = u;
                fb = fu;
                return 1;
			}
            else
			{
				if (fu > fb)
				{
					cx = u;
					fc = fu;
					return 1;
				}
            }
            u = cx + gold * (cx - bx);
            fu = sign * FunctionValue(u);
		}
        else
		{
			if ((cx - u) * (u - ulim) > 0)
			{
				fu = sign * FunctionValue(u);
				if (fu < fc)
				{
					bx = cx;
					cx = u;
					u = cx + gold * (cx - bx);
					fb = fc;
					fc = fu;
					fu = sign * FunctionValue(u);
				}
            }
			else
			{
				if ((u - ulim) * (ulim - cx) > 0 || FloatEqual((u - ulim) * (ulim - cx),0))
				{
					u = ulim;
					fu = sign * FunctionValue(u);
				}
				else
				{
					u = cx + gold * (cx - bx);
					fu = sign * FunctionValue(u);
				}
			}
		}
		ax = bx;
		bx = cx;
		cx = u;
		fa = fb;
		fb = fc;
		fc = fu;
	}
	return 1;
}

//1维函数极小值的黄金分割法
template <class _Ty>
_Ty ExtremumGold1D(_Ty ax, _Ty bx, _Ty cx, int sign, _Ty tol, _Ty& xmin)
{
    _Ty x0,x1,x2,x3,f0,f1,f2,f3,r = GoldNo;
    _Ty c = 1.0 - GoldNo;
    x0 = ax;
    x3 = cx;
    if (Abs(cx - bx) > Abs(bx - ax))
	{
		x1 = bx;
        x2 = bx + c * (cx - bx);
	}
    else
	{
        x2 = bx;
        x1 = bx - c * (bx - ax);
    }
    f1 = sign * FunctionValue(x1);
    f2 = sign * FunctionValue(x2);
    while (Abs(x3 - x0) > tol * (Abs(x1) + Abs(x2)))
	{
        if (f2 < f1)
		{
            x0 = x1;
            x1 = x2;
            x2 = r * x1 + c * x3;
            f0 = f1;
            f1 = f2;
            f2 = sign * FunctionValue(x2);
		}
        else
		{
            x3 = x2;
            x2 = x1;
            x1 = r * x2 + c * x0;
            f3 = f2;
            f2 = f1;
            f1 = sign * FunctionValue(x1);
        }
    }
    if (f1 < f2)
	{
        xmin = x1;
		if (sign == -1)
			f1 = sign * f1;
        return f1;
	}
    else
	{
        xmin = x2;
		if (sign == -1)
			f2 = sign * f2;
        return f2;
    }
}

//1维函数极小值点的不用导数布伦特(Brent)法
template <class _Ty>
_Ty ExtremumBrentNonDerivative1D(_Ty ax, _Ty bx, _Ty cx, int sign, 
												_Ty tol, _Ty& xmin)
{
    int done,iter,itmax = 100;
    _Ty d,fu,r,q,p,xm,tol1,tol2,a,b,cgold = 1.0 - GoldNo;
    _Ty u,etemp,dum,v,w,x,e,fx,fv1,fw,zeps = DOUBLEERROR;
    a = ax;
    if (cx < ax)
	{
		a = cx;
	}
    b = ax;
    if (cx > ax)
	{
		b = cx;
	}
    v = bx;
    w = v;
    x = v;
    e = 0.0;
    fx = sign * FunctionValue(x);
    fv1 = fx;
    fw = fx;
    for (iter = 1; iter<=itmax; iter++)
	{
        xm = 0.5 * (a + b);
        tol1 = tol * Abs(x) + zeps;
        tol2 = 2.0 * tol1;
        if (Abs(x - xm) <= tol2 - 0.5 * (b - a))
		{
			break;
		}
        done = -1;
        if (Abs(e) > tol1)
		{
            r = (x - w) * (fx - fv1);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0)
			{
				p = -p;
			}
            q = Abs(q);
            etemp = e;
            e = d;
            dum = Abs(0.5 * q * etemp);
            if (Abs(p) < dum && p > q * (a - x) && p < q * (b - x))
			{
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
				{
                    d = Abs(tol1) * Sgn(xm - x);
                }
                done = 0;
            }
        }
        if (done)
		{
            if (x > xm || FloatEqual(x,xm))
			{
                e = a - x;
			}
            else
			{
                e = b - x;
            }
            d = cgold * e;
        }
        if (Abs(d) > tol1 || FloatEqual(d,tol1))
		{
            u = x + d;
		}
        else
		{
            u = x + Abs(tol1) * Sgn(d);
        }
        fu = sign * FunctionValue(u);
        if (fu < fx || FloatEqual(fu,fx))
		{
            if (u > x || FloatEqual(u,x))
			{
                a = x;
			}
            else
			{
                b = x;
            }
            v = w;
            fv1 = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
		}
        else
		{
            if (u < x)
			{
                a = u;
			}
            else
			{
                b = u;
            }
            if (fu < fw || FloatEqual(fu,fw) || FloatEqual(w,x))
			{
                v = w;
                fv1 = fw;
                w = u;
                fw = fu;
			}
            else
			{
				if (fu < fv1 || FloatEqual(fu,fv1) ||  FloatEqual(v,x) ||  FloatEqual(v,w))
				{
					v = u;
					fv1 = fu;
				}
            }
        }
    }
    if (iter > itmax)
	{
		cout<<" ExtremumBrentNonDerivative1D exceed maximum iterations."<<endl;
	}
    xmin = x;
    if(sign == -1) 
	{
		fx *= sign;
	}
	return fx;
}

#endif		// _EXTREMUM_INL

