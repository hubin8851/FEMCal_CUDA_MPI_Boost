//polynomials.inl		多项式与连分式头文件
// Ver 1.0.0.0
// 版权所有(C) 2002
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _POLYNOMIALS_INL
#define _POLYNOMIALS_INL

//求一维实(复)多项式值
template<class T, class U>
inline U 	//内联函数
PolyValueOneDim(valarray<T>& dCoff, size_t stNo, U dX)
{		//dCoff为一维多项式系数数组，stNo为多项式次数，dX为自变量x
	U dValue = dCoff[stNo-1];
	for(int st = stNo - 2; st > -1; st--)
		dValue = dValue * dX + dCoff[st];
	return(dValue);			//返回多项式值
}

//求一维多项式组值
template<class T, class V, class U>
inline void 
PolyValueOneDimGroup(valarray<T>& dCoff, valarray<V>& dX, valarray<U>& dValue)
{	//dCoff为一维多项式系数数组，dX存放自变量，dValue存放对应多项式值
  	int i,j,mm,nn,ll,t,s,kk,k;
    V y, z;

    int n = dCoff.size();	//多项式项数，其最高次数为n-1
	int m = dX.size();		//给定自变量个数

	valarray<V> b(2*n);

	//b=malloc(2*n*sizeof(double));
    y=dCoff[n-1];

    for(i=0; i<n; i++) b[i] = dCoff[i]/y;
    k = log(n-0.5) / log(2.0) + 1;
    nn = 1;
    for(i=0; i<k; i++) nn=2*nn;
    for(i=n; i<nn-1; i++) b[i]=0.0;
    b[nn-1] = 1.0;
    t = nn;
	s = 1;
    for(i=1; i<=k-1; i++)
    {
		t = t / 2;
		mm = -t;
        for(j=1; j<=s; j++)
        {
			mm = mm + t + t;
			b[mm-1] -= 1.0; 
            for(kk=2; kk<=t; kk++)
              b[mm-kk] = b[mm-kk] - b[mm-1] * b[mm+t-kk];
        }
        s = s + s;
    }
    for(kk=1; kk<=m; kk++)
      { 
		for(i=0; i<=(nn-2)/2; i++)
           dCoff[i] = dX[kk-1] + b[2*i];
        mm = 1;
		z = dX[kk-1];
        for(i=1; i<=k-1; i++)
        {
			mm = mm + mm; 
			ll = mm + mm;
			z = z * z;
            for(j=0; j<nn; j=j+ll)
              dCoff[j/2]=dCoff[j/2]+dCoff[(j+mm)/2]*(z+b[j+mm-1]);
        }
        z = z * z / dX[kk-1];
        if(nn!=n) dCoff[0] -= z;
        dValue[kk-1] = dCoff[0] * y;
    }
}

//求二维多项式值
template<class T, class U>
inline U	//内联函数
PolyValueTwoDim(matrix<T>& dCoff, U dX, U dY)
{		//dCoff为二维多项式系数数组，dX为自变量x，dY为自变量y
	U dValue = 0, dTemp = 1;

	size_t styNo = dCoff.GetColNum();//二维多项式自变量y的项数，其最高次数为styNo-1
	size_t stxNo = dCoff.GetRowNum();//二维多项式自变量x的项数，其最高次数为stxNo-1

	for(size_t st = 0; st < stxNo; st++)
	{
		U dUS = dCoff(st, styNo-1) * dTemp;
		for(int sr = styNo-2; sr > -1; sr--)
			dUS = dUS * dY + dCoff(st, sr) * dTemp;
		dValue = dValue + dUS;
		dTemp = dTemp * dX;
	}
	return(dValue);		//返回多项式值
}

//两一维多项式相乘
template<class T>
inline void
PolyMultip(valarray<T>& dCoffP, valarray<T>& dCoffQ, valarray<T>& dCoffS)
{	//dCoffP被乘多项式P(x)系数数组，最高次数为stNoP-1
	//dCoffQ乘多项式Q(x)系数数组，最高次数为stNoQ-1
	//dCoffS积多项式S(x)系数数组，最高次数为stNoP+stNoQ-1
	//最后结果得到积多项式S(x)的系数
	size_t stNoP = dCoffP.size();//一维多项式P(x)项数，最高次数为stNoP-1
	size_t stNoQ = dCoffQ.size();//一维多项式Q(x)项数，最高次数为stNoQ-1
	size_t stNoS = dCoffS.size();//一维多项式S(x)项数，最高次数为stNoS-1，stNoS=stNoP+stNoQ-1

	for(size_t st = 0; st < stNoS; st++) dCoffS[st]=0.0;

	for(st = 0; st < stNoP; st++)
		for(size_t sr = 0; sr < stNoQ; sr++)
			dCoffS[st+sr]=dCoffS[st+sr]+dCoffP[st]*dCoffQ[sr];
}

//两一维多项式除法
template<class T>
inline int
PolyDiv(valarray<T>& dCoffP, valarray<T>& dCoffQ, valarray<T>& dCoffS, valarray<T>& dCoffR)
{	//dCoffP被除多项式P(x)系数数组，最高次数为stNoP-1
	//dCoffQ除多项式Q(x)系数数组，最高次数为stNoQ-1
	//dCoffS商多项式S(x)系数数组，最高次数为stNoP-stNoQ
	//dCoffR余多项式R(x)系数数组，最高次数为stNoR-1
	//最后得到的结果是商多项式S(x)和余多项式R(x)的系数
	size_t stNoP = dCoffP.size();//一维多项式P(x)项数，最高次数为stNoP-1
	size_t stNoQ = dCoffQ.size();//一维多项式Q(x)项数，最高次数为stNoQ-1
	size_t stNoS = dCoffS.size();//一维多项式S(x)项数，最高次数为stNoS-stNoQ
	size_t stNoR = dCoffR.size();//一维多项式R(x)项数，最高次数为stNoR-1

	for(size_t st = 0; st < stNoS; st++) dCoffS[st] = 0.0;

	if(FloatEqual(Abs(dCoffQ[stNoQ-1]),0)) return (0);
	
	size_t stk = stNoP - 1;
	
	for(st = stNoS; st > 0; st--)
	{
		dCoffS[st-1] = dCoffP[stk] / dCoffQ[stNoQ-1];
		size_t stm = stk;
		for(size_t sr = 1; sr < stNoQ; sr++)
		{
			dCoffP[stm-1] = dCoffP[stm-1] - dCoffS[st-1] * dCoffQ[stNoQ-sr-1];
			stm = stm-1;
		}
		stk = stk - 1;
	}

	for(st = 0; st < stNoR; st++) dCoffR[st] = dCoffP[st];
	
	return (1);
}

// 计算连分式函数值
template<class T>
inline T
FractionValue(valarray<T>& dXpara, valarray<T>& dCoff, T dX)
{
	size_t stNo = dXpara.size();	//连分式的节数
	T dValue = dCoff[stNo-1];

	for(int st=stNo-2; st>-1; st--)
	{
		if (FloatEqual(dValue,0))
			dValue = 1.0e+35 * (dX-dXpara[st]) / Abs(dX-dXpara[st]);
        else
			dValue = dCoff[st] + (dX-dXpara[st]) / dValue;
    }
    return(dValue);				//返回连分式值
}

#endif //_POLYNOMIALS_INL

