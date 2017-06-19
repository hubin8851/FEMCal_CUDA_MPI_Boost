// MatrixSymmetryRegularCholesky.cpp  
//对称正定阵的乔里斯基(Cholesky)分解及求其行列式值

#include <iostream>		//输入输出流
#include "Matrix.h"		//矩阵类及相关函数等的定义
using namespace std;	//名字空间

void main()				// 定义控制台应用程序的入口点
{
	const double dma[4][4] = {	{  1.0,  2.0,  3.0,  4.0 },
								{  5.0,  6.0,  7.0,  8.0 },
								{  9.0, 10.0, 11.0, 12.0 },
								{ 13.0, 14.0, 15.0, 16.0 }};

	const double dmb[3][3] = {	{  3.0, -3.0, -2.0 },
								{ -3.0,  8.0,  4.0 },
								{  -2.0, 4.0,  3.0 } };

	const double dmc[4][4] = {	{ 5.0,  7.0,  6.0,  5.0},
								{ 7.0, 10.0,  8.0,  7.0},
								{ 6.0,  8.0, 10.0,  9.0},
								{ 5.0,  7.0,  9.0, 10.0} };

	matrix<double> matA(&dma[0][0], 4, 4);
	matrix<double> matB(&dmb[0][0], 3, 3);
	matrix<double> matC(&dmc[0][0], 4, 4);
	
	cout << "matA: " << endl;
	long double mDet = MatrixSymmetryRegularCholesky(matA);
	if(mDet > 0) 
	{
		cout << endl << "Determinant(matA) = " << mDet << endl;
		MatrixLinePrint(matA);
	}
	else cout << "The matA is not symmetry and regular." << endl;

	cout << endl << "matB: " << endl;
	mDet = MatrixSymmetryRegularCholesky(matB);
	if(mDet > 0) 
	{
		cout << "Determinant(matB) = " << mDet << endl;
		MatrixLinePrint(matB);
	}
	else cout << "The matB is not symmetry and regular." << endl;

	cout << endl << "matC: " << endl;
	mDet = MatrixSymmetryRegularCholesky(matC);
	if(mDet > 0) 
	{
		cout << "Determinant(matC) = " << mDet << endl;
		MatrixLinePrint(matC);
	}
	else cout << "The matC is not symmetry and regular." << endl;
}

