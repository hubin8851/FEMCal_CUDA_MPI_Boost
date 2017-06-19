//Statistic.inl		数据处理与回归分析头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _STATISTIC_INL		//避免多次编译
#define _STATISTIC_INL

//随机样本分析
template <class _Ty>
void StatRandomSample(valarray<_Ty>& x, _Ty x0, _Ty h, int l, 
			valarray<_Ty>& dt, valarray<int>& g, valarray<int>& q)
{
    char a[50];
	int n = x.size();			//随机样本点数
	int m = g.size();			//直方图中区间总数

    dt[0] = 0;
    for(int i=0; i<n; i++) dt[0]=dt[0]+x[i]/n;
    dt[1] = 0;

    for(i=0; i<n; i++)
      dt[1] = dt[1] + (x[i] - dt[0]) * (x[i] - dt[0]);

    dt[1] = dt[1] / n;
    dt[2] = sqrt(dt[1]);

    for(i=0; i<m; i++)
    {
		q[i]=0;
        _Ty s=x0+(i+0.5)*h-dt[0];
        s=exp(-s*s/(2.0*dt[1]));
        g[i]=n*s*h/(dt[2]*2.5066);
    }
    _Ty s=x0+m*h;
    for(i=0; i<n; i++)
		if((x[i] - x0) > 0 || FloatEqual((x[i]-x0), 0))
        if((s - x[i]) > 0 || FloatEqual((s-x[i]), 0))
        {
			int j = (x[i] - x0) / h;
            q[j] = q[j] + 1;
        }
    if (l==0) return;	//不打印直方图
    cout << endl << "n = " << n << endl;
    cout << endl << "x0 = " << x0 << "\t h = " << h << "\t m = " << m << endl;
    cout << endl << "xa = " << dt[0] << "\t s = " << dt[1];
	cout << "\t t = " << dt[2] << endl << endl;
    int k = 1;
	int z = 0;
    for (i=0; i<m; i++)   if (q[i]>z) z=q[i];

    while(z>50)
    {
		z=z/2;
		k=2*k;
	}
    for(i=0; i<m; i++)
    {
		s=x0+(i+0.5)*h;
        for (int j=0; j<=49; j++) a[j]=' ';
        j=q[i]/k;
        for(z=0; z<j; z++) a[z]='X';
        j = g[i]/k;
        if((j>0) && (j<=50)) a[j] = '*';
		cout << s << "\t" << q[i] << "\t" ;
        for (j=0; j<50; j++)
			cout << a[j];
        cout << endl;
    }
    cout << endl << "k = " << k << endl;
    cout << endl;
	return;		//正常结束程序
}

//一元线性回归分析
template <class _Ty>
void LinearRegression1D(valarray<_Ty>& x, valarray<_Ty>& y, 
								valarray<_Ty>& a, valarray<_Ty>& dt)
{
	_Ty xx(0), yy(0), e(0), f(0), u(0), p(0), umax(0), umin(1.0e+30);

	int n = x.size();			//观测点数

    for(int i=0; i<n; i++)
    {
		xx=xx+x[i]/n;
		yy=yy+y[i]/n;
	}
    
    for(i=0; i<n; i++)
    {
		_Ty q=x[i]-xx;
		e=e+q*q;
        f=f+q*(y[i]-yy);
    }
    a[1]=f/e;
	a[0]=yy-a[1]*xx;
    _Ty q=0.0;
    for(i=0; i<n; i++)
    {
		_Ty s=a[1]*x[i]+a[0];
        q=q+(y[i]-s)*(y[i]-s);
        p=p+(s-yy)*(s-yy);
        e=fabs(y[i]-s);
        if (e>umax) umax=e;
        if (e<umin) umin=e;
        u=u+e/n;
    }
    dt[1]=sqrt(q/n);
    dt[0]=q;
	dt[2]=p;
    dt[3]=umax;
	dt[4]=umin;
	dt[5]=u;
}

//n元线性回归分析
template <class _Ty>
void LinearRegressionND(matrix<_Ty>& x, valarray<_Ty>& y, 
			valarray<_Ty>& a, valarray<_Ty>& dt, valarray<_Ty>& v)
{
	int m = x.GetRowNum();			//自变量个数
	int n = x.GetColNum();			//观测数据的组数
	
	//为调用平方根法求解对称正定方程组函数准备，其两参数须matrix类型
	matrix<_Ty> b((m+1),(m+1));
	matrix<_Ty> aa(m+1,1);

	for(int i=0;i<m+1;i++) aa(i,0) = a[i];

    int mm = m + 1;
	b(m, m) = n;

    for(int j=0; j<m; j++)
    {
		_Ty p(0);
        for(int i=0; i<n; i++)	p = p + x(j,i);
		b(m,j) = p;
		b(j,m) = p;
    }
    for(i=0; i<m; i++)
      for(j=i; j<m; j++)
      {
		  _Ty p(0);
          for(int k=0; k<n; k++)	p=p+x(i,k)*x(j,k);
		  b(j,i) = p;
		  b(i,j) = p;

      }

	aa(m,0) = 0.0;
    for(i=0; i<n; i++) aa(m,0) += y[i];
    for(i=0; i<m; i++)
    {
		aa(i,0) = 0.0;
		for(j=0; j<n; j++)	aa(i,0) = aa(i,0) + x(i,j) * y[j];
    }

	//调用平方根法求解对称正定方程组的函数
    if(LE_SymmetryRegularEuationSquareRoot(b,aa)<1)
	{
		cout << "Matrix is not Symmetry and Regular!" << endl;
		return;
	}
    
	_Ty yy(0);
    for (i=0; i<n; i++)	yy = yy + y[i] / n;
    _Ty q(0), e(0), u(0);
    for(i=0; i<n; i++)
    {
		_Ty p = aa(m,0);
		for(j=0; j<m; j++)	p = p + aa(j,0) * x(j,i);
        q=q+(y[i]-p)*(y[i]-p);
        e=e+(y[i]-yy)*(y[i]-yy);
        u=u+(yy-p)*(yy-p);
    }
    _Ty s = sqrt(q/n);
    _Ty r = sqrt(1.0-q/e);
    for(j=0; j<m; j++)
    {
		_Ty p(0);
        for (i=0; i<n; i++)
        {
			_Ty pp = aa(m,0);
            for(int k=0; k<m; k++)
              if(k!=j) pp = pp + aa(k,0) * x(k,i);
            p = p + (y[i] - pp) * (y[i] - pp);
        }
        v[j] = sqrt(1.0 - q / p);
    }
    dt[0] = q;
	dt[1] = s;
	dt[2] = r;
	dt[3] = u;

	for(i=0; i<m+1; i++) a[i] = aa(i,0);
}

//逐步回归分析
template <class _Ty>
void StepwiseRegression(matrix<_Ty>& x, _Ty f1, _Ty f2, 
		_Ty eps, valarray<_Ty>& xx, valarray<_Ty>& b, valarray<_Ty>& v, 
				valarray<_Ty>& s, valarray<_Ty>& dt, valarray<_Ty>& ye, 
									valarray<_Ty>& yr, matrix<_Ty>& r)
{
	int ii,l;
    _Ty phi,sd,fmi,fmx;
	
	int k = x.GetRowNum();			//观测点数
	int n = x.GetColNum()-1;		//自变量x的个数

    int m = n + 1;
	_Ty q(0);
    for(int j=0; j<=n; j++)
    {
		_Ty z(0);
        for(int i=0; i<k; i++)	z = z + x(i,j) / k;
        xx[j] = z;
    }
    for(int i=0; i<=n; i++)
      for(j=0; j<=i; j++)
      {
		  _Ty z(0);
          for(ii=0; ii<k; ii++)
			  z = z + (x(ii,i) - xx[i]) * (x(ii,j) - xx[j]);
          r(i,j) = z;
      }
    for(i=0; i<=n; i++)	ye[i]=sqrt(r(i,i));
    for (i=0; i<=n; i++)
      for (j=0; j<=i; j++)
      {
		  r(i,j) = r(i,j) / (ye[i] * ye[j]);
          r(j,i) = r(i,j);
      }
    phi = k - 1.0;
    sd = ye[n] / sqrt(k-1.0);
    int it(1);
    while(it==1)
    {
		it = 0;
		_Ty vmi(1.0e+35), vmx(0);
        int imi(-1), imx(-1);
        for(i=0; i<=n; i++)
        {
			v[i]=0.0;
			b[i]=0.0;
			s[i]=0.0;
		}
        for(i=0; i<n; i++)
          if(r(i,i) > eps || FloatEqual(r(i,i),eps))
           {
			  v[i] = r(i,n) * r(n,i) / r(i,i);
              if (v[i] > 0.0 || FloatEqual(v[i],0))
              { 
				  if(v[i]>vmx)
                  {
					vmx=v[i];
					imx=i;
				  }
              }
              else
              {
				  b[i] = r(i,n) * ye[n] / ye[i];
                  s[i] = sqrt(r(i,i)) * sd / ye[i];
                  if(Abs(v[i])<vmi)
                  {
					  vmi=Abs(v[i]);
					  imi=i;
				  }
              }
          }
        if(phi!=n-1.0)
        {
			_Ty z(0);
            for (i=0; i<n; i++)	z = z + b[i] * xx[i];
            b[n] = xx[n] - z;
			s[n] = sd;
			v[n] = q;
        }
        else
        {
			b[n]=xx[n];
			s[n]=sd;
		}
        fmi = vmi * phi / r(n,n);
        fmx = (phi-1.0) * vmx / (r(n,n)-vmx);
        if((fmi<f2)||(fmx>=f1))
        {
			if(fmi<f2)
            {
				phi=phi+1.0;
				l=imi;
			}
            else
            {
				phi=phi-1.0;
				l=imx;
			}
            for(i=0; i<=n; i++)
              if (i!=l)
                for (j=0; j<=n; j++)
                  if (j!=l)
                    r(i,j)=r(i,j)-(r(l,j)/r(l,l))*r(i,l);
            for(j=0; j<=n; j++)
              if(j!=l)
                r(l,j)=r(l,j)/r(l,l);
            for(i=0; i<=n; i++)
              if(i!=l)
                r(i,l)=-r(i,l)/r(l,l);
            r(l,l)=1.0/r(l,l);
            q=r(n,n)*ye[n]*ye[n];
            sd=sqrt(r(n,n)/phi)*ye[n];
            dt[0]=sqrt(1.0-r(n,n));
            dt[1]=(phi*(1.0-r(n,n)))/((k-phi-1.0)*r(n,n));
            it=1;
        }
    }
    for(i=0; i<k; i++)
    {
		_Ty z(0);
        for(j=0; j<n; j++)	z=z+b[j]*x(i,j);
        ye[i]=b[n]+z; 
		yr[i]=x(i,n)-ye[i];
    }
}

//半对数数据相关
template <class _Ty>
void HalfLogarithmCorrelation(valarray<_Ty>& x, valarray<_Ty>& y, 
											_Ty t, valarray<_Ty>& a)
{
	_Ty xx(0), yy(0), dx(0), dxy(0);

	int n = x.size();		//数据点数

    for(int i=0; i<n; i++)
    {
		xx = xx + x[i] / n; 
        yy = yy + log(y[i]) / log(t) / n;
    }
    
    for(i=0; i<n; i++)
    {
		a[2] = x[i] - xx; 
		dx = dx + a[2] * a[2];
        dxy = dxy + a[2] * (log(y[i]) / log(t) - yy);
    }
    a[1] = dxy / dx;
	a[0] = yy - a[1] * xx;
    a[0] = a[0] * log(t); 
	a[0] = exp(a[0]);
    a[2] = a[6] = a[4] = 0.0;
	a[5] = 1.0e+30;
    for(i=0; i<n; i++)
    {
		a[3] = a[1] * x[i] * log(t);
		a[3] = a[0] * exp(a[3]);
        a[2] = a[2] + (y[i] - a[3]) * (y[i] - a[3]);
        dx = Abs(y[i] - a[3]);
        if(dx>a[4]) a[4] = dx;
        if(dx<a[5]) a[5] = dx;
        a[6] = a[6] + dx / n;
    }
    a[3] = sqrt(a[2] / n);
}

//对数数据相关
template <class _Ty>
void LogarithmCorrelation(valarray<_Ty>& x, 
								valarray<_Ty>& y, valarray<_Ty>& a)
{
	_Ty xx(0), yy(0), dx(0), dxy(0);

	int n = x.size();		//数据点数

    for(int i=0; i<n; i++)
    {
		xx = xx + log(x[i]) / n; 
        yy = yy + log(y[i]) / n;
    }
    
    for(i=0; i<n; i++)
    {
		a[2] = log(x[i]) - xx; 
		dx = dx + a[2] * a[2];
        dxy = dxy + a[2] * (log(y[i]) - yy);
    }
	a[1] = dxy / dx;
	a[0] = yy - a[1] * xx;
	a[0] = exp(a[0]);
    a[2] = a[6] = a[4] = 0.0;
	a[5] = 1.0e+30;
    for(i=0; i<n; i++)
    {
		a[3] = a[1] * log(x[i]);
		a[3] = a[0] * exp(a[3]);
        a[2] = a[2] + (y[i] - a[3]) * (y[i] - a[3]);
        dx = Abs(y[i] - a[3]);
        if(dx>a[4]) a[4] = dx;
        if(dx<a[5]) a[5] = dx;
        a[6] = a[6] + dx / n;
    }
    a[3] = sqrt(a[2] / n);
}

#endif		// _STATISTIC_INL
