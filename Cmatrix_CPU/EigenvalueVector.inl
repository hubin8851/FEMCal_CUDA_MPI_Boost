// EigenvalueVector.inl		计算特征值特征向量函数(方法)定义
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002

// 最后修改: 2002.8.10

#ifndef _EIGENVALUEVECTOR_INL
#define _EIGENVALUEVECTOR_INL

//约化对称阵为对称三对角阵的豪斯荷尔德变换法
template <class _Ty>
int HouseholderTransform(matrix<_Ty>& a, matrix<_Ty>& q, 
									valarray<_Ty>& b, valarray<_Ty>& c)
{
	int j, k, stRank;
    _Ty h, f, g, h2;

	if(MatrixSymmetry(a)!=true)
		return int(0);			//矩阵a不是对称阵
	
	stRank = a.GetColNum();		// 矩阵阶数
	
    for(int i=0; i<stRank; i++)
		for(int j=0; j<stRank; j++)
			q(i,j)=a(i,j);

    for(i=stRank-1; i>=1; i--)
    { 
		h=0.0;
        if(i>1)
          for(k=0; k<i; k++)
			  h=h+q(i,k)*q(i,k);
        if(FloatEqual(h,0))
        {
			c[i]=0.0;
            if(i==1) c[i]=q(i,i-1);
            b[i]=0.0;
        }
        else
        { 
			c[i]=sqrt(h);
            if(q(i,i-1)>0.0) c[i]=-c[i];
            h=h-q(i,i-1)*c[i];
            q(i,i-1)=q(i,i-1)-c[i];
            f=0.0;
            for(j=0; j<i; j++)
            { 
				q(j,i) = q(i,j) / h;
                g = 0.0;
                for(k=0; k<=j; k++)
                  g = g + q(j,k) * q(i,k);
                if(j+1<i)
                  for(k=j+1; k<i; k++)
                    g = g + q(k,j) * q(i,k);
                c[j] = g / h;
                f = f + g * q(j,i);
            }
            h2 = f / (h+h);
            for(j=0; j<i; j++)
            { 
				f = q(i,j);
                g = c[j] -h2 * f;
                c[j] = g;
                for(k=0; k<=j; k++)
                    q(j,k)=q(j,k)-f*c[k]-g*q(i,k);
            }
            b[i]=h;
        }
    }
    for(i=0; i<stRank-1; i++) c[i]=c[i+1];
    c[stRank-1]=0.0;
    b[0]=0.0;
    for(i=0; i<stRank; i++)
    { if((b[i]!=0.0)&&(i-1>=0))
          for(j=0; j<i; j++)
          { 
			  g=0.0;
              for(k=0; k<i; k++)
				g=g+q(i,k)*q(k,j);
              for(k=0; k<i; k++)
                q(k,j)=q(k,j)-g*q(k,i);
          }
        b[i]=q(i,i);
		q(i,i)=1.0;
        if(i-1>=0)
          for(j=0; j<i; j++)
          { 
			  q(i,j)=0.0;
			  q(j,i)=0.0;
		  }
    }
	return int(1);	//正常返回
}

//实对称三角阵全部特征值及特征向量QR法
template <class _Ty>
int EigenvalueVectorRealTriangleQR(valarray<_Ty>& b, 
					valarray<_Ty>& c, matrix<_Ty>& q, _Ty eps, int l)
{
	int i, k, m, it, stRank;
    _Ty h, g, p, r, e, s, d(0), f(0);

	stRank = q.GetColNum();		// 矩阵阶数
	
	c[stRank-1]=0.0;

    for(int j=0; j<stRank; j++)
    { 
		it=0;
        h=eps*(Abs(b[j])+Abs(c[j]));
        if(h>d) d=h;
        m=j;
        while((m<stRank)&&(Abs(c[m])>d)) m++;
        if(m!=j)
        {
			do
            {
				if(it==l)
                { 
					cout << "fail" << endl;
                    return(-1);				//出错
                }
                it++;
                g=b[j];
                p=(b[j+1]-g)/(2.0*c[j]);
                r=sqrt(p*p+1.0);
                if(p>0.0 || FloatEqual(p,0)) b[j]=c[j]/(p+r);
                else b[j]=c[j]/(p-r);
                h=g-b[j];
                for(i=j+1; i<stRank; i++)	b[i]=b[i]-h;
                f=f+h; 
				p=b[m]; 
				e=1.0; 
				s=0.0;
                for(i=m-1; i>=j; i--)
                {
					g=e*c[i];
					h=e*p;
                    if(Abs(p)>=Abs(c[i]))
                    {
						e=c[i]/p;
						r=sqrt(e*e+1.0);
                        c[i+1]=s*p*r;
						s=e/r; 
						e=1.0/r;
                    }
                    else
					{
						e=p/c[i]; 
						r=sqrt(e*e+1.0);
                        c[i+1]=s*c[i]*r;
                        s=1.0/r; 
						e=e/r;
                    }
                    p=e*b[i]-s*g;
                    b[i+1]=h+s*(e*g+s*b[i]);
                    for(k=0; k<stRank; k++)
                    { 
                        h=q(k,i+1);
						q(k,i+1)=s*q(k,i)+e*h;
                        q(k,i)=e*q(k,i)-s*h;
                    }
                }
                c[j]=s*p;
				b[j]=e*p;
            }while(Abs(c[j])>d);
        }
        b[j]=b[j]+f;
    }
    for(i=0; i<stRank; i++)
    { 
		k=i;
		p=b[i];
        if(i+1<stRank)
        {
			j=i+1;
            while((j<stRank)&&(b[j]<=p))
            {
				k=j;
				p=b[j];
				j=j+1;
			}
		}
        if(k!=i)
        { 
			b[k]=b[i]; 
			b[i]=p;
            for(j=0; j<stRank; j++)
            {
                p=q(j,i);
				q(j,i)=q(j,k);
				q(j,k)=p;
            }
        }
    }
    return(1);	//正常返回
}

//约化一般实矩阵为赫申伯格阵的初等相似变换法
template <class _Ty>
int HessenbergTransform(matrix<_Ty>& a)
{
	int i, stRank;
    _Ty t;
	
	stRank = a.GetColNum();		// 矩阵阶数
	if(stRank != a.GetRowNum())
		return int(0);			//矩阵a不是方阵

    for(int k=1; k<stRank-1; k++)
    { 
		_Ty d(0);
        for(int j=k; j<stRank; j++)
        { 
			t=a(j,k-1);
            if(Abs(t)>Abs(d))
            {
				d=t; 
				i=j;
			}
        }
        if(FloatNotEqual(d,0))
        { 
			if(i!=k)
            {
				for(j=k-1; j<stRank; j++)
					swap(a(i,j), a(k,j));
                for(j=0; j<stRank; j++)
					swap(a(j,i), a(j,k));
            }
            for(i=k+1; i<stRank; i++)
            { 
				t=a(i,k-1)/d; 
				a(i,k-1)=0.0;
                for(j=k; j<stRank; j++)
                    a(i,j)=a(i,j)-t*a(k,j);
                for(j=0; j<stRank; j++)
                    a(j,k)=a(j,k)+t*a(j,i);
            }
        }
    }
	return int(1);
}

//求赫申伯格阵全部特征值QR法
template <class _Ty>
int EigenvalueVectorHessenbergQR(matrix<_Ty>& a,  
						valarray<complex<_Ty> >& uv, _Ty eps, int jt)
{
	
	int i, j, k, l, it(0), stRank;
    _Ty b,c,w,g,xy,p,q,r,x,s,e,f,z,y;

	stRank = a.GetColNum();		// 矩阵阶数

    while(stRank!=0)
    { 
		l=stRank-1;
        while((l>0)&&(Abs(a(l,l-1))>eps*(Abs(a(l-1,l-1))+Abs(a(l,l))))) l=l-1;
        if(l==stRank-1)
        { 
			uv[stRank-1]=complex<_Ty>(a(stRank-1,stRank-1),0);
            stRank=stRank-1;
			it=0;
        }
        else if(l==stRank-2)
        { 
			b=-(a(stRank-1,stRank-1)+a(stRank-2,stRank-2));
            c=a(stRank-1,stRank-1)*a(stRank-2,stRank-2)-a(stRank-1,stRank-2)*a(stRank-2,stRank-1);
            w=b*b-4.0*c;
            y=sqrt(Abs(w));
            if(w>0.0)
            { 
				xy=1.0;
                if(b<0.0) xy=-1.0;
				uv[stRank-1]=complex<_Ty>((-b-xy*y)/2.0, 0);
				uv[stRank-2]=complex<_Ty>(c/(uv[stRank-1].real()), 0);
            }
            else
            {
				uv[stRank-1]=complex<_Ty>(-b/2.0, y/2.0);
				uv[stRank-2]=complex<_Ty>(uv[stRank-1].real(),-uv[stRank-1].imag());
            }
            stRank=stRank-2; 
			it=0;
        }
        else
        { 
			if(it>=jt)
            { printf("fail\n");
                return(-1);		//出错！
            }
            it=it+1;
            for(j=l+2; j<stRank; j++)	a(j,j-2)=0.0;
            for(j=l+3; j<stRank; j++)    a(j,j-3)=0.0;
            for(k=l; k<stRank-1; k++)
            { 
				if(k!=l)
                {
					p=a(k,k-1);
					q=a(k+1,k-1);
                    r=0.0;
                    if(k!=stRank-2) r=a(k+2,k-1);
                }
                else
                {
					x=a(stRank-1,stRank-1)+a(stRank-2,stRank-2);
                    y=a(stRank-2,stRank-2)*a(stRank-1,stRank-1)-a(stRank-2,stRank-1)*a(stRank-1,stRank-2);
                    p=a(l,l)*(a(l,l)-x)+a(l,l+1)*a(l+1,l)+y;
                    q=a(l+1,l)*(a(l,l)+a(l+1,l+1)-x);
                    r=a(l+1,l)*a(l+2,l+1);
                }
                if(FloatNotEqual((Abs(p)+Abs(q)+Abs(r)),0))
                { 
					xy=1.0;
                    if(p<0.0) xy=-1.0;
                    s=xy*sqrt(p*p+q*q+r*r);
                    if(k!=l) a(k,k-1)=-s;
                    e=-q/s;
					f=-r/s; 
					x=-p/s;
                    y=-x-f*r/(p+s);
                    g=e*r/(p+s);
                    z=-x-e*q/(p+s);
                    for(j=k; j<stRank; j++)
                    {	
                        p=x*a(k,j)+e*a(k+1,j);
                        q=e*a(k,j)+y*a(k+1,j);
                        r=f*a(k,j)+g*a(k+1,j);
                        if(k!=stRank-2)
                        { 
                            p=p+f*a(k+2,j);
                            q=q+g*a(k+2,j);
                            r=r+z*a(k+2,j);
							a(k+2,j)=r;
                        }
                        a(k+1,j)=q; 
						a(k,j)=p;
                    }
                    j=k+3;
                    if(j>=stRank-1) j=stRank-1;
                    for(i=l; i<=j; i++)
                    {	
                        p=x*a(i,k)+e*a(i,k+1);
                        q=e*a(i,k)+y*a(i,k+1);
                        r=f*a(i,k)+g*a(i,k+1);
                        if(k!=stRank-2)
                        { 
                            p=p+f*a(i,k+2);
                            q=q+g*a(i,k+2);
                            r=r+z*a(i,k+2); 
							a(i,k+2)=r;
                        }
                        a(i,k+1)=q; 
						a(i,k)=p;
                    }
                }
            }
        }
    }
    return(1);
}

//实对称阵特征值及特征向量雅可比法
template <class _Ty>
int EigenvalueVectorRealSymmetryJacobi(matrix<_Ty>& a,  
										matrix<_Ty>& v, _Ty eps, int jt)
{
	int j, p, q, l(1), stRank;
    _Ty fm,cn,sn,omega,x,y,d;
	
	if(MatrixSymmetry(a)!=true)	//不是对称阵
		return(0);					

	stRank = a.GetColNum();		// 矩阵阶数

    for(int i=0; i<stRank; i++)
    {
		v(i,i)=1.0;
        for(j=0; j<stRank; j++)
          if(i!=j) v(i,j)=0.0;
    }
    while(1)
    { 
		fm=0.0;
        for(i=1; i<stRank; i++)
        for(j=0; j<i; j++)
        {
			d=Abs(a(i,j));
            if((i!=j)&&(d>fm))
            {
				fm=d; 
				p=i;
				q=j;
			}
        }
        if(fm<eps)  return(1);
        if(l>jt)  return(-1);
        l=l+1;
        x=-a(p,q);
		y=(a(q,q)-a(p,p))/2.0;
        omega=x/sqrt(x*x+y*y);
        if(y<0.0) omega=-omega;
        sn=1.0+sqrt(1.0-omega*omega);
        sn=omega/sqrt(2.0*sn);
        cn=sqrt(1.0-sn*sn);
        fm=a(p,p);
        a(p,p)=fm*cn*cn+a(q,q)*sn*sn+a(p,q)*omega;
        a(q,q)=fm*sn*sn+a(q,q)*cn*cn-a(p,q)*omega;
        a(p,q)=0.0; 
		a(q,p)=0.0;
        for(j=0; j<stRank; j++)
        if((j!=p)&&(j!=q))
        {
            fm=a(p,j);
            a(p,j)=fm*cn+a(q,j)*sn;
            a(q,j)=-fm*sn+a(q,j)*cn;
        }
        for(i=0; i<stRank; i++)
          if((i!=p)&&(i!=q))
          { 
              fm=a(i,p);
              a(i,p)=fm*cn+a(i,q)*sn;
              a(i,q)=-fm*sn+a(i,q)*cn;
          }
        for(i=0; i<stRank; i++)
        {
            fm=v(i,p);
            v(i,p)=fm*cn+v(i,q)*sn;
            v(i,q)=-fm*sn+v(i,q)*cn;
        }
    }
    return(1);
}

//实对称阵特征值及特征向量雅可比过关法
template <class _Ty>
int EigenvalueVectorRealSymmetryJacobiB(matrix<_Ty>& a,  
												matrix<_Ty>& v, _Ty eps)
{
	int j, p, q, stRank;
    _Ty fm, cn, sn, omega, x, y, d;

	if(MatrixSymmetry(a)!=true)
		return(-1);

	stRank = a.GetColNum();		//矩阵阶数

    for(int i=0; i<stRank; i++)
    { v(i,i)=1.0;
        for(j=0; j<stRank; j++)
          if(i!=j) v(i,j)=0.0;
    }
    _Ty ff(0);
    for(i=1; i<stRank; i++)
		for(j=0; j<i; j++)
		{
			d=a(i,j);
			ff=ff+d*d;
		}
    ff=sqrt(2.0*ff);
  loop0:
    ff=ff/(1.0*stRank);
  loop1:
    for(i=1; i<stRank; i++)
        for(j=0; j<i; j++)
        { 
			d=Abs(a(i,j));
            if(d>ff)
            { 
				p=i;
				q=j;
                goto loop;
            }
        }
	if(ff<eps) return(1);	//成功返回
    goto loop0;
  loop:
    x=-a(p,q);
	y=(a(q,q)-a(p,p))/2.0;
    omega=x/sqrt(x*x+y*y);
    if(y<0.0) omega=-omega;
    sn=1.0+sqrt(1.0-omega*omega);
    sn=omega/sqrt(2.0*sn);
    cn=sqrt(1.0-sn*sn);
    fm=a(p,p);
    a(p,p)=fm*cn*cn+a(q,q)*sn*sn+a(p,q)*omega;
    a(q,q)=fm*sn*sn+a(q,q)*cn*cn-a(p,q)*omega;
    a(p,q)=0.0;
	a(q,p)=0.0;
    for(j=0; j<stRank; j++)
        if((j!=p)&&(j!=q))
        {
            fm=a(p,j);
            a(p,j)=fm*cn+a(q,j)*sn;
            a(q,j)=-fm*sn+a(q,j)*cn;
        }
    for(i=0; i<stRank; i++)
        if((i!=p)&&(i!=q))
        {
			fm=a(i,p);
            a(i,p)=fm*cn+a(i,q)*sn;
            a(i,q)=-fm*sn+a(i,q)*cn;
        }
    for(i=0; i<stRank; i++)
    { 
        fm=v(i,p);
        v(i,p)=fm*cn+v(i,q)*sn;
        v(i,q)=-fm*sn+v(i,q)*cn;
    }
    goto loop1;
}

#endif			// _EIGENVALUEVECTOR_INL

