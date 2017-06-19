//MatrixToeplitzInversionTrench.cpp  特兰持(Trench)法求托伯利兹(Toeplitz)矩阵逆

#include <iostream>		//输入输出流
#include "Matrix.h"		//矩阵类及相关函数等的定义
using namespace std;	//名字空间


void main()				// 定义控制台应用程序的入口点
{
    const double t[6] = {10.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    const double tuo[6] = { 0.0,-1.0,-2.0,-3.0,-4.0,-5.0};
	valarray<double> vt(t,6), vtuo(tuo,6);

	const double dma[6][6] = 
	{
		{10.0, 5.0, 4.0, 3.0, 2.0, 1.0},
		{-1.0,10.0, 5.0, 4.0, 3.0, 2.0},
		{-2.0,-1.0,10.0, 5.0, 4.0, 3.0},
		{-3.0,-2.0,-1.0,10.0, 5.0, 4.0},
		{-4.0,-3.0,-2.0,-1.0,10.0, 5.0},
		{-5.0,-4.0,-3.0,-2.0,-1.0,10.0}
	};

	matrix<double> matT(&dma[0][0], 6, 6);
	matrix<double> matB(matT);			//用matT初始化matB

	cout << "matT : " << endl;

	MatrixLinePrint(matT);				//按行输出矩阵matT
	
	cout << endl << "matB = Inversion(matT) : " << endl;

	if(MatrixToeplitzInversionTrench(vt, vtuo, matB)>0)	//现在matB是matT的逆
	MatrixLinePrint(matB);				//按行输出矩阵matB

	cout << endl << "matI = matT * matB : " << endl;

	matrix<double> matI(matT * matB);	//matT * matB得到单位阵matI
	MatrixLinePrint(matI);				//按行输出矩阵matI
}

