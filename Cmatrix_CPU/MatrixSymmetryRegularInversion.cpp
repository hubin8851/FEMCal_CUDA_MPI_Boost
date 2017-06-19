//MatrixSymmetryRegularInversion.cpp  用“变量循环重新编号法”法求对称正定矩阵逆

#include <iostream>		//输入输出流
#include "Matrix.h"		//矩阵类及相关函数等的定义
using namespace std;	//名字空间

void main()				// 定义控制台应用程序的入口点
{
	const double dmb[4][4] = 
	{
		{0.2368,0.2471,0.2568,1.2671},
        {1.1161,0.1254,0.1397,0.1490},
        {0.1582,1.1675,0.1768,0.1871},
        {0.1968,0.2071,1.2168,0.2271}
	};

	const double dma[4][4] = 
	{
		{ 5.0,  7.0,  6.0,  5.0},
		{ 7.0, 10.0,  8.0,  7.0},
		{ 6.0,  8.0, 10.0,  9.0},
		{ 5.0,  7.0,  9.0, 10.0} 
	};

	const double dmc[4][4] = 
	{
		{  3.0, -3.0, -2.0,  4.0 },
		{  5.0, -5.0,  1.0,  8.0 },
		{ 11.0,  8.0,  5.0, -7.0 },
		{  5.0, -1.0, -3.0, -1.0 } 
	};

	matrix<double> matA(&dma[0][0], 4, 4);
	matrix<double> matB(matA);
	
	cout << "matA: " << endl;

	MatrixLinePrint(matA);			//按行输出矩阵matA

	int iValue = MatrixSymmetryRegularInversion(matA);
	if(iValue > 1)
	{
		cout << endl << "Inversion(matA) : " << endl;

		MatrixLinePrint(matA);			//按行输出矩阵matA的逆

		cout << endl << "Inversion(matA) * matA : " << endl;

		matA = matA * matB;
		MatrixLinePrint(matA);			//按行输出矩阵matA与逆的乘积
	}
	else
		cout << "matA 不是对称正定阵！" << endl;
}

