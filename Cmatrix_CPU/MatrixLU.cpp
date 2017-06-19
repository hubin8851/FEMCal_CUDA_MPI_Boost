//MatrixLU.cpp			LU分解及求其行列式值

#include <iostream>		//输入输出流
#include "Matrix.h"		//矩阵类及相关函数等的定义
using namespace std;	//名字空间

void main()				// 定义控制台应用程序的入口点
{
	void DislayMatrixInfo(int iRet, matrix<double>& A, 
							matrix<double>& L, matrix<double>& U);
	int iRet;

	const double dma[4][4] = {	{  1.0,  2.0,  3.0,  4.0 },
								{  5.0,  6.0,  7.0,  8.0 },
								{  9.0, 10.0, 11.0, 12.0 },
								{ 13.0, 14.0, 15.0, 16.0 }};

	const double dmb[4][4] = {	{  3.0, -3.0, -2.0, 2.0 },
								{ -3.0,  8.0,  4.0, 1.0 },
								{ -2.0,  4.0,  6.0, 3.0 }, 
								{  2.0,  1.0,  3.0, 9.0 } };

	const double dmc[4][4] = {	{2.0, 4.0,  4.0, 2.0},
								{3.0, 3.0, 12.0, 6.0},
								{2.0, 4.0, -1.0, 2.0},
								{4.0, 2.0,  1.0, 1.0}	};

	matrix<double> matA(&dma[0][0], 4, 4);
	matrix<double> matB(&dmb[0][0], 4, 4);
	matrix<double> matC(&dmc[0][0], 4, 4);

	matrix<double> matA_L(4,4);
	matrix<double> matA_U(4,4);
	
	iRet = MatrixLU(matA, matA_L, matA_U);

	cout << "matA : " << endl;
	DislayMatrixInfo(iRet, matA, matA_L, matA_U); 

	iRet = MatrixLU(matB, matA_L, matA_U);
	cout << endl << "matB : " << endl ;

	DislayMatrixInfo(iRet, matB, matA_L, matA_U); 

	iRet = MatrixLU(matC, matA_L, matA_U);
	cout << endl << "matC : " << endl;

	DislayMatrixInfo(iRet, matC, matA_L, matA_U); 
}

void DislayMatrixInfo(int iRet, matrix<double>& A, 
						matrix<double>& L, matrix<double>& U)
{
	if(iRet==-1)
	{	
		cout << "矩阵不是方阵！" << endl;
	}
	else
		if(iRet==0)
			cout << "矩阵分解失败！" << endl;
		else
		{
			MatrixLinePrint(A);
			cout << endl << "L : " << endl;
			MatrixLinePrint(L);
			cout << endl << "U : " << endl;
			MatrixLinePrint(U);
		}
}

