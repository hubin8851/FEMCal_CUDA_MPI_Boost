//MatrixSymmetryRegular.cpp  判断对称正定阵

//#include <fstream>	//文件流头文件
#include <iostream>		//输入输出流
#include "Matrix.h"		//矩阵类及相关函数等的定义
using namespace std;	//名字空间

void main()				// 定义控制台应用程序的入口点
{
	void DislayMatrixInfo(int sValue);
	int returnValue;

	const double dma[4][4] = 
	{
		{  1.0,  2.0,  3.0,  4.0 },
		{  5.0,  6.0,  7.0,  8.0 },
		{  9.0, 10.0, 11.0, 12.0 },
		{ 13.0, 14.0, 15.0, 16.0 }
	};

	const double dmb[3][3] = 
	{
		{  3.0, -3.0, -2.0 },
		{ -3.0,  8.0,  4.0 },
		{  -2.0, 4.0,  3.0 } 
	};

	const double dmc[4][4] = 
	{
		{ 5.0,  7.0,  6.0,  5.0},
		{ 7.0, 10.0,  8.0,  7.0},
		{ 6.0,  8.0, 10.0, -9.0},
		{ 5.0,  7.0,  9.0, 10.0} 
	};

	const matrix<double> matA(&dma[0][0], 4, 4);
	const matrix<double> matB(&dmb[0][0], 3, 3);
	const matrix<double> matC(&dmc[0][0], 4, 4);

	cout << "matA : " << endl;

	int sym(1);				//判别矩阵是否为对称阵
	returnValue = MatrixSymmetryRegular(matA, sym);
	DislayMatrixInfo(returnValue);

	cout << endl << "matB : " << endl ;

	sym = 1;				//判别矩阵是否为对称阵
	returnValue = MatrixSymmetryRegular(matB, sym);
	DislayMatrixInfo(returnValue);
	
	cout << endl << "matC : "  << endl;

	sym = 0;				//不判别矩阵是否为对称阵
	returnValue = MatrixSymmetryRegular(matC, sym);
	DislayMatrixInfo(returnValue);
}

void DislayMatrixInfo(int sValue)
{
	switch(sValue)
	{
		case 1:		cout << "The matrix is regular." << endl;
					break;

		case 2:		cout << "The matrix is symmetry and regular." << endl;
					break;

		case -1:	cout << "The matrix is not a square." << endl;
					break;

		case -2:	cout << "The matrix is not symmetry." << endl;
					break;

		case -3:	cout << "The matrix is not regular." << endl;
					break;

		case 0:		cout << "The matrix is not regular." << endl;;
					break;
	}
}

