//polynomials.inl		����ʽ������ʽͷ�ļ�
// Ver 1.0.0.0
// ��Ȩ����(C) 2002
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _POLYNOMIALS_INL
#define _POLYNOMIALS_INL

//��һάʵ(��)����ʽֵ
template<class T, class U>
inline U 	//��������
PolyValueOneDim(valarray<T>& dCoff, size_t stNo, U dX)
{		//dCoffΪһά����ʽϵ�����飬stNoΪ����ʽ������dXΪ�Ա���x
	U dValue = dCoff[stNo-1];
	for(int st = stNo - 2; st > -1; st--)
		dValue = dValue * dX + dCoff[st];
	return(dValue);			//���ض���ʽֵ
}

//��һά����ʽ��ֵ
template<class T, class V, class U>
inline void 
PolyValueOneDimGroup(valarray<T>& dCoff, valarray<V>& dX, valarray<U>& dValue)
{	//dCoffΪһά����ʽϵ�����飬dX����Ա�����dValue��Ŷ�Ӧ����ʽֵ
  	int i,j,mm,nn,ll,t,s,kk,k;
    V y, z;

    int n = dCoff.size();	//����ʽ����������ߴ���Ϊn-1
	int m = dX.size();		//�����Ա�������

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

//���ά����ʽֵ
template<class T, class U>
inline U	//��������
PolyValueTwoDim(matrix<T>& dCoff, U dX, U dY)
{		//dCoffΪ��ά����ʽϵ�����飬dXΪ�Ա���x��dYΪ�Ա���y
	U dValue = 0, dTemp = 1;

	size_t styNo = dCoff.GetColNum();//��ά����ʽ�Ա���y������������ߴ���ΪstyNo-1
	size_t stxNo = dCoff.GetRowNum();//��ά����ʽ�Ա���x������������ߴ���ΪstxNo-1

	for(size_t st = 0; st < stxNo; st++)
	{
		U dUS = dCoff(st, styNo-1) * dTemp;
		for(int sr = styNo-2; sr > -1; sr--)
			dUS = dUS * dY + dCoff(st, sr) * dTemp;
		dValue = dValue + dUS;
		dTemp = dTemp * dX;
	}
	return(dValue);		//���ض���ʽֵ
}

//��һά����ʽ���
template<class T>
inline void
PolyMultip(valarray<T>& dCoffP, valarray<T>& dCoffQ, valarray<T>& dCoffS)
{	//dCoffP���˶���ʽP(x)ϵ�����飬��ߴ���ΪstNoP-1
	//dCoffQ�˶���ʽQ(x)ϵ�����飬��ߴ���ΪstNoQ-1
	//dCoffS������ʽS(x)ϵ�����飬��ߴ���ΪstNoP+stNoQ-1
	//������õ�������ʽS(x)��ϵ��
	size_t stNoP = dCoffP.size();//һά����ʽP(x)��������ߴ���ΪstNoP-1
	size_t stNoQ = dCoffQ.size();//һά����ʽQ(x)��������ߴ���ΪstNoQ-1
	size_t stNoS = dCoffS.size();//һά����ʽS(x)��������ߴ���ΪstNoS-1��stNoS=stNoP+stNoQ-1

	for(size_t st = 0; st < stNoS; st++) dCoffS[st]=0.0;

	for(st = 0; st < stNoP; st++)
		for(size_t sr = 0; sr < stNoQ; sr++)
			dCoffS[st+sr]=dCoffS[st+sr]+dCoffP[st]*dCoffQ[sr];
}

//��һά����ʽ����
template<class T>
inline int
PolyDiv(valarray<T>& dCoffP, valarray<T>& dCoffQ, valarray<T>& dCoffS, valarray<T>& dCoffR)
{	//dCoffP��������ʽP(x)ϵ�����飬��ߴ���ΪstNoP-1
	//dCoffQ������ʽQ(x)ϵ�����飬��ߴ���ΪstNoQ-1
	//dCoffS�̶���ʽS(x)ϵ�����飬��ߴ���ΪstNoP-stNoQ
	//dCoffR�����ʽR(x)ϵ�����飬��ߴ���ΪstNoR-1
	//���õ��Ľ�����̶���ʽS(x)�������ʽR(x)��ϵ��
	size_t stNoP = dCoffP.size();//һά����ʽP(x)��������ߴ���ΪstNoP-1
	size_t stNoQ = dCoffQ.size();//һά����ʽQ(x)��������ߴ���ΪstNoQ-1
	size_t stNoS = dCoffS.size();//һά����ʽS(x)��������ߴ���ΪstNoS-stNoQ
	size_t stNoR = dCoffR.size();//һά����ʽR(x)��������ߴ���ΪstNoR-1

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

// ��������ʽ����ֵ
template<class T>
inline T
FractionValue(valarray<T>& dXpara, valarray<T>& dCoff, T dX)
{
	size_t stNo = dXpara.size();	//����ʽ�Ľ���
	T dValue = dCoff[stNo-1];

	for(int st=stNo-2; st>-1; st--)
	{
		if (FloatEqual(dValue,0))
			dValue = 1.0e+35 * (dX-dXpara[st]) / Abs(dX-dXpara[st]);
        else
			dValue = dCoff[st] + (dX-dXpara[st]) / dValue;
    }
    return(dValue);				//��������ʽֵ
}

#endif //_POLYNOMIALS_INL

