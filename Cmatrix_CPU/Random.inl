//random.inl			�������������������
// Ver 1.0.0.0
// ��Ȩ����(C) ����, 2002
// ����޸�: 2002.5.31.

#ifndef _RANDOM_INL
#define _RANDOM_INL

//����һ��[0,1]�����ھ��ȷֲ�α�����
inline double 
rand_01_One(double& seed)
{		//seedΪ���������
	seed = ((unsigned long)seed) % MODUL65536;		//��65536Ϊģ������ȡ��
	seed = RandCoef2053 * (seed) + RandCoef13849;	//�ˡ���ϵ��
	seed = ((unsigned long)seed) % MODUL65536;		//��65536Ϊģ������ȡ��
	double rand = (seed) / (double)MODUL65536;		//ʹ�������[0,1]����
	return( rand );									//�������ֵ
}

//�������[0,1]�����ھ��ȷֲ�α�����
inline void 
rand_01_Series(double& seed, valarray<double>& dp, const size_t stCount)
{	//seedΪ��������ӣ�dp������ɵ�������У�srCountָ��Ҫ���ɵ����������
	for(size_t st=0; st<stCount; st++)
	{
		seed = RandCoef2053 * (seed) + RandCoef13849;
		seed = (unsigned long)seed % MODUL65536;	//��65536Ϊģ������ȡ��
		dp[st] = seed / (double)MODUL65536;
	}
}

//��������[a,b]������һ�����ȷֲ�α�������
inline size_t 
rand_ab_One(size_t a, size_t b, size_t& seed)
{			//a,bΪ�������Ҷ˵㣬seedΪ����
	size_t rand;
	size_t stk = b - a + 1;
	size_t stl = 2;
	while(stl < stk) stl = stl + stl;
	size_t modul = 4 * stl;
	stk = seed;
	size_t sti = 1;
	while(sti <= 1)
	{
		stk = 5 * stk;
		stk = stk % modul;
		stl = stk /4 +a;
		if (stl<=b)
		{
			rand=stl; 
			sti=sti+1;
		}
      }
    seed=stk;
    return(rand);
}

//�����������[a,b]�����ھ��ȷֲ�α�����������
inline void 
rand_ab_Series(size_t a, size_t b, size_t& seed, valarray<size_t>& sp, size_t stCount)
{   //a,bΪ�������Ҷ˵㣬seedΪ���ӣ�sp�����������У�stCountָ�����������
	size_t stk = b - a + 1;
	size_t stl = 2;
	while(stl < stk) stl = stl + stl;
	size_t modul = 4 * stl;
	stk = seed;
	size_t sti = 0;
	while(sti < stCount)
	{
		stk = 5 * stk;
		stk = stk % modul;
		stl = stk / 4 +a;
		if (stl <= b)
		{
			sp[sti]=stl; 
			sti = sti + 1;
		}
    }
    seed=stk;
}

//���������ֵ�뷽���һ����̬�ֲ������
//�����Ĳ�����mu��ro(�̺ͦ�ֵ)Ҫ���ȸ�����Ҳ�ɴ����������л�á�
//����㷨�е�nֵҪȡ�ĸ���(����12)������Ҫ�Լ��޸ġ�
inline double 
rand_NormalDistributing_One(double mu, double ro, double& seed)
{
	double rand = 0.0;
	for (size_t st=0; st<12; st++)
	{
		seed = RandCoef2053 * (seed) + RandCoef13849;	//�ˡ���ϵ��
		seed =(unsigned long)seed % MODUL65536;			//��65536Ϊģ������ȡ��
		rand = rand + (seed) / (double)MODUL65536;
	}
	rand = mu + ro * (rand - 6.0);
	return(rand);
}

//���������ֵ�뷽�����̬�ֲ����������
//�����Ĳ�����mu��ro(�̺ͦ�ֵ)Ҫ���ȸ�����Ҳ�ɴ����������л�á�
//����㷨�е�nֵҪȡ�ĸ���(����12)������Ҫ�Լ��޸ġ�
inline void 
rand_NormalDistributing_Series(double mu, double ro, double seed, valarray<double>& dp, size_t stCount)
{			//seedΪ���ӣ�dp�����������У�stCountָ��Ҫ���ɵ����������
	for (size_t st = 0; st < stCount; st++)
	{
		double rand = 0.0;
		for (size_t sr=1; sr<=12; sr++)
		{
			seed = (seed) * RandCoef2053 + RandCoef13849;
			seed = (unsigned long)seed % MODUL65536;			//��65536Ϊģ������ȡ��
			rand = rand + (seed) / (double)MODUL65536;
		}
		dp[st] = mu + ro * (rand - 6.0);
	}
}

#endif //_RANDOM_INL

