//Transform.inl		数学变换头文件
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _TRANSFORM_INL		//避免多次编译
#define _TRANSFORM_INL

//傅里叶级数逼近
template <class _Ty>
void FourierSeriesApproach(valarray<_Ty>& f,  
								valarray<_Ty>& a, valarray<_Ty>& b)
{
    _Ty u0, c1(1), s1(0);
	
	int n = a.size() - 1;
    _Ty t = 6.283185306 / (2.0 * n + 1.0);
    _Ty c = cos(t);
	_Ty s = sin(t);
    t = 2.0 / (2.0 * n + 1.0);
		
    for(int i=0; i<=n; i++)
    {
		_Ty u1(0), u2(0);
        for(int j=2*n; j>=1; j--)
        {
			u0 = f[j] + 2.0 * c1 * u1 - u2;
            u2 = u1;
			u1 = u0;
        }
        a[i] = t * (f[0] + u1 * c1 - u2);
        b[i] = t * u1 * s1;
        u0 = c * c1 - s * s1;
		s1 = c * s1 + s * c1;
		c1 = u0;
    }

}

//快速傅里叶变换
template <class _Ty>
void FourierTransform(valarray<complex<_Ty> >& pp, 
						valarray<complex<_Ty> >& ff, int l, int il)
{
	int n = pp.size();		//输入的点数，n>0
	int it, k(0);				
	
	for(it = 1; ; it *= 2)	//求k值, 要满足n等于2的k次方
	{
		if(it == n) break;
		k++;
		if(k > 64) return;
	}

	for(it = 0; it<n; it++)
    {
		int m = it, is(0);
        for(int i=0; i<k; i++)
        {
			int j = m / 2;
			is = 2 * is + (m - 2 * j);
			m = j;
		}
		
		ff[it] = pp[is];		
    }

	pp[0] = complex<_Ty>(1.0, 0.0);
	
	_Ty p = 6.283185306 / (1.0 * n);

    pp[1] = complex<_Ty>(cos(p), -sin(p));
    
	if(l != 0) 
		pp[1] =complex<_Ty>(pp[1].real(), -pp[1].imag());
    
	for(int i = 2; i < n; i ++)
    {
		p = pp[i-1].real() * pp[1].real();
		
		_Ty q = pp[i-1].imag() * pp[1].imag();
        _Ty s = (pp[i-1].real() + pp[i-1].imag()) * (pp[1].real() + pp[1].imag());
        
		pp[i] = complex<_Ty>(p - q, s - p - q);
    }
    for(it = 0; it<n-1; it = it + 2)
    {
		_Ty vr=ff[it].real(); 
		_Ty vi=ff[it].imag();

        ff[it] = complex<_Ty>(vr + ff[it+1].real(), vi + ff[it+1].imag());
		ff[it+1] = complex<_Ty>(vr - ff[it+1].real(), vi - ff[it+1].imag());
		
    }
    
	int m=n/2;
	int nv(2);
    
	for(int l0 = k-2; l0>=0; l0--)
    {
		m = m / 2;
		nv = 2 * nv;
        for(it = 0; it <= (m - 1) * nv; it = it + nv)
          for(int j = 0; j < (nv / 2); j ++)
          {
			  p = pp[m*j].real() * ff[it+j+nv/2].real();
              _Ty q = pp[m*j].imag() * ff[it+j+nv/2].imag();
              _Ty s = pp[m*j].real() + pp[m*j].imag();
              
			  s *= (ff[it+j+nv/2].real() + ff[it+j+nv/2].imag());
              _Ty poddr = p - q;
			  _Ty poddi = s - p - q;
              
			  ff[it+j+nv/2] = complex<_Ty>(ff[it+j].real() - poddr, ff[it+j].imag() - poddi);
              ff[it+j] = complex<_Ty>(ff[it+j].real() + poddr, ff[it+j].imag() + poddi);
          }
    }

    if(l != 0)
      for(i = 0; i < n; i ++)
      {
		  ff[i] = complex<_Ty>(ff[i].real() / (1.0*n),ff[i].imag() / (1.0*n));
      }
    if(il!=0)
      for(i=0; i<n; i++)
      {
		  pp[i] = complex<_Ty>(sqrt(ff[i].real() * ff[i].real() + ff[i].imag() * ff[i].imag()),
			                           pp[i].imag());
          if(Abs(ff[i].real()) < 0.000001 * Abs(ff[i].imag()))
          { 
			  if((ff[i].imag() * ff[i].real()) > 0) pp[i] = complex<_Ty>(pp[i].real(), 90.0);
              else pp[i] = complex<_Ty>(pp[i].real(), -90.0);
          }
          else
			  pp[i] = complex<_Ty>(pp[i].real(), (atan(ff[i].imag() / ff[i].real()) * 360.0 / 6.283185306));
			
      }
}

//快速沃什变换
template <class _Ty>
void WalshTransform(valarray<_Ty>& p, valarray<_Ty>& x)
{
	int i, k(0);
	int n = p.size();			//输入序列的长度，n>0
	for(i=1; ; i*=2)			//n等于2的k次方
	{
		if(i==n) break;
		k++;
		if(k>64) return;
	}

	int l=n, m(1), it=2, ii=n/2;
    x[0]=1;
	x[ii]=2;
    for(i=1; i<k; i++)
    {
		m=m+m;
		l=l/2;
		it=it+it;
        for(int j=0; j<m; j++)	x[j*l+l/2]=it+1-x[j*l];
    }
    for(i=0; i<n; i++)
    {
		ii=x[i]-1;
		x[i]=p[ii];
	}
    l=1;
    for(i=1; i<=k; i++)
     {
		m=n/(2*l)-1;
        for(int j=0; j<=m; j++)
        {
			it=2*l*j;
            for(int is=0; is<l; is++)
            {
				_Ty q=x[it+is]+x[it+is+l];
                x[it+is+l]=x[it+is]-x[it+is+l];
                x[it+is]=q;
            }
        }
        l=2*l;
    }
}

//五点三次平滑(曲线拟合)
template <class _Ty>
void Smooth5_3(valarray<_Ty>& y, valarray<_Ty>& yy)
{
	int n = y.size();	//给定等距观测点的个数。要求n>=5
    if(n<5)
		for(int i=0; i<n; i++) yy[i]=y[i];
    else
    {
		yy[0]=69.0*y[0]+4.0*y[1]-6.0*y[2]+4.0*y[3]-y[4];
        yy[0]=yy[0]/70.0;
        yy[1]=2.0*y[0]+27.0*y[1]+12.0*y[2]-8.0*y[3];
        yy[1]=(yy[1]+2.0*y[4])/35.0;
        for (int i=2; i<n-2; i++)
        {
			yy[i]=-3.0*y[i-2]+12.0*y[i-1]+17.0*y[i];
            yy[i]=(yy[i]+12.0*y[i+1]-3.0*y[i+2])/35.0;
        }
        yy[n-2]=2.0*y[n-5]-8.0*y[n-4]+12.0*y[n-3];
        yy[n-2]=(yy[n-2]+27.0*y[n-2]+2.0*y[n-1])/35.0;
        yy[n-1]=-y[n-5]+4.0*y[n-4]-6.0*y[n-3];
        yy[n-1]=(yy[n-1]+4.0*y[n-2]+69.0*y[n-1])/70.0;
    }
}

//离散随机线性系统的卡尔曼滤波
template <class _Ty>
int SievingKalman(matrix<_Ty>& f, matrix<_Ty>& q, matrix<_Ty>& r, 
					matrix<_Ty>& h, matrix<_Ty>& y, matrix<_Ty>& x, 
									matrix<_Ty>& p, matrix<_Ty>& g)
{
	int i,j,kk,ii,jj,js;

	int n = g.GetRowNum();		//动态系统的维数
	int m = g.GetColNum();		//观测系统的维数
	int k = y.GetRowNum();		//观测序列的长度

	matrix<_Ty> e(m,m);
	int l=m;
    if(l<n) l=n;
	matrix<_Ty> a(l,l);
	matrix<_Ty> b(l,l);
    for(i=0; i<n; i++)
      for(j=0; j<n; j++)
      {
		  a(i,j)=0.0;
          for(kk=0; kk<n; kk++)	a(i,j)=a(i,j)+p(i,kk)*f(j,kk);
      }
    for(i=0; i<n; i++)
    for(j=0; j<n; j++)
    { 
	  p(i,j)=q(i,j);
      for(kk=0; kk<n; kk++) p(i,j)=p(i,j)+f(i,kk)*a(kk,j);
    }
    for(ii=2; ii<=k; ii++)
    { 
		for(i=0; i<n; i++)
	    for(j=0; j<m; j++)
		{
			jj=i*l+j;
			a(i,j)=0.0;
			for(kk=0; kk<n; kk++)
				a(i,j)=a(i,j)+p(i,kk)*h(j,kk);
		}
		for(i=0; i<m; i++)
			for(j=0; j<m; j++)
			{
				jj=i*m+j; 
				e(i,j)=r(i,j);
				for(kk=0; kk<n; kk++)
					e(i,j)=e(i,j)+h(i,kk)*a(kk,j);
			}

		js=MatrixInversionGS(e);	//调用全选主元高斯-约当消去法求逆

		if(js==0)
		{
			cout << "Error! No inversion!" << endl;
			return(js);				//求增益矩阵过程中求逆失败
		}
        
		for(i=0; i<n; i++)
        for(j=0; j<m; j++)
        { 
			jj=i*m+j;
			g(i,j)=0.0;
            for(kk=0; kk<m; kk++)
				g(i,j)=g(i,j)+a(i,kk)*e(j,kk);
        }
        for(i=0; i<n; i++)
        { 
			jj=(ii-1)*n+i; 
			x(ii-1,i)=0.0;
            for(j=0; j<n; j++)
				x(ii-1,i)=x(ii-1,i)+f(i,j)*x(ii-2,j);
        }
        for(i=0; i<m; i++)
        { 
			jj=i*l; 
			b(i,0)=y(ii-1,i);
            for(j=0; j<n; j++)
				b(i,0)=b(i,0)-h(i,j)*x(ii-1,j);
        }
        for(i=0; i<n; i++)
        { 
			jj=(ii-1)*n+i;
            for(j=0; j<m; j++)
				x(ii-1,i)=x(ii-1,i)+g(i,j)*b(j,0);
        }
        if(ii<k)
        { 
			for(i=0; i<n; i++)
            for(j=0; j<n; j++)
            {
				jj=i*l+j;
				a(i,j)=0.0;
                for(kk=0; kk<m; kk++)
					a(i,j)=a(i,j)-g(i,kk)*h(kk,j);
                if(i==j) a(i,j)=1.0+a(i,j);
            }
            for(i=0; i<n; i++)
            for(j=0; j<n; j++)
            { 
				jj=i*l+j; 
				b(i,j)=0.0;
                for(kk=0; kk<n; kk++)
					b(i,j)=b(i,j)+a(i,kk)*p(kk,j);
            }
            for(i=0; i<n; i++)
            for(j=0; j<n; j++)
            { 
				jj=i*l+j; 
				a(i,j)=0.0;
                for(kk=0; kk<n; kk++)
					a(i,j)=a(i,j)+b(i,kk)*f(j,kk);
            }
            for(i=0; i<n; i++)
            for(j=0; j<n; j++)
            { 
				jj=i*n+j;
				p(i,j)=q(i,j);
                for(kk=0; kk<n; kk++)
					p(i,j)=p(i,j)+f(i,kk)*a(j,kk);
          }
      }
  }
  return(js);
}

//a-b-r滤波
template <class _Ty>
void SievingABR(valarray<_Ty>& x, _Ty t, _Ty a, _Ty b,
											_Ty r, valarray<_Ty>& y)
{
    _Ty aa(0), vv(0), ss(0);

	int n = x.size();		//量测数据的点数
    
	for(int i=0; i<n; i++)
    { 
		_Ty s1=ss+t*vv+t*t*aa/2.0;
        _Ty v1=vv+t*aa;
		_Ty a1=aa;
        ss=s1+a*(x[i]-s1);
		y[i]=ss;
        vv=v1+b*(x[i]-s1);
        aa=a1+2.0*r*(x[i]-s1)/(t*t);
    }

}

#endif		// _TRANSFORM_INL
