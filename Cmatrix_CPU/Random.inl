//random.inl			随机函数（方法）定义
// Ver 1.0.0.0
// 版权所有(C) 何渝, 2002
// 最后修改: 2002.5.31.

#ifndef _RANDOM_INL
#define _RANDOM_INL

//产生一个[0,1]区间内均匀分布伪随机数
inline double 
rand_01_One(double& seed)
{		//seed为随机数种子
	seed = ((unsigned long)seed) % MODUL65536;		//以65536为模对种子取余
	seed = RandCoef2053 * (seed) + RandCoef13849;	//乘、加系数
	seed = ((unsigned long)seed) % MODUL65536;		//以65536为模对种子取余
	double rand = (seed) / (double)MODUL65536;		//使随机数在[0,1]区间
	return( rand );									//返回随机值
}

//产生多个[0,1]区间内均匀分布伪随机数
inline void 
rand_01_Series(double& seed, valarray<double>& dp, const size_t stCount)
{	//seed为随机数种子，dp存放生成的随机数列，srCount指定要生成的随机数个数
	for(size_t st=0; st<stCount; st++)
	{
		seed = RandCoef2053 * (seed) + RandCoef13849;
		seed = (unsigned long)seed % MODUL65536;	//以65536为模对种子取余
		dp[st] = seed / (double)MODUL65536;
	}
}

//产生任意[a,b]区间内一个均匀分布伪随机整数
inline size_t 
rand_ab_One(size_t a, size_t b, size_t& seed)
{			//a,b为区间左右端点，seed为种子
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

//产生多个任意[a,b]区间内均匀分布伪随机整数序列
inline void 
rand_ab_Series(size_t a, size_t b, size_t& seed, valarray<size_t>& sp, size_t stCount)
{   //a,b为区间左右端点，seed为种子，sp存放随机数序列，stCount指定随机数个数
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

//产生任意均值与方差的一个正态分布随机数
//函数的参数：mu和ro(μ和σ值)要事先给定，也可从其它程序中获得。
//如果算法中的n值要取的更大(超过12)，程序要稍加修改。
inline double 
rand_NormalDistributing_One(double mu, double ro, double& seed)
{
	double rand = 0.0;
	for (size_t st=0; st<12; st++)
	{
		seed = RandCoef2053 * (seed) + RandCoef13849;	//乘、加系数
		seed =(unsigned long)seed % MODUL65536;			//以65536为模对种子取余
		rand = rand + (seed) / (double)MODUL65536;
	}
	rand = mu + ro * (rand - 6.0);
	return(rand);
}

//产生任意均值与方差的正态分布随机数序列
//函数的参数：mu和ro(μ和σ值)要事先给定，也可从其它程序中获得。
//如果算法中的n值要取的更大(超过12)，程序要稍加修改。
inline void 
rand_NormalDistributing_Series(double mu, double ro, double seed, valarray<double>& dp, size_t stCount)
{			//seed为种子，dp存放随机数序列，stCount指定要生成的随机数个数
	for (size_t st = 0; st < stCount; st++)
	{
		double rand = 0.0;
		for (size_t sr=1; sr<=12; sr++)
		{
			seed = (seed) * RandCoef2053 + RandCoef13849;
			seed = (unsigned long)seed % MODUL65536;			//以65536为模对种子取余
			rand = rand + (seed) / (double)MODUL65536;
		}
		dp[st] = mu + ro * (rand - 6.0);
	}
}

#endif //_RANDOM_INL

