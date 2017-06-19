// MatrixRank.cpp		用全选主元高斯消元法求一般矩阵的秩

//#include <fstream>	//文件流头文件
#include <iostream>		//输入输出流
#include "Matrix.h"		//矩阵类及相关函数等的定义
using namespace std;	//名字空间

void main()				// 定义控制台应用程序的入口点
{
	size_t mRank;

	const double dma[4][4] = {	{  1.0,  2.0,  3.0,  4.0 },
								{  5.0,  6.0,  7.0,  8.0 },
								{  9.0, 10.0, 11.0, 12.0 },
								{ 13.0, 14.0, 15.0, 16.0 } };

	const double dmb[3][3] = {	{  3.0, -3.0, -2.0 },
								{ -3.0,  8.0,  4.0 },
								{  -2.0, 4.0, -3.0 } };

	const double dmc[5][4] = {	{ -1.0,  2.0,  3.0,  4.0},
								{  5.0, -6.0,  7.0,  8.0},
								{  9.0, 10.0, 11.0, 12.0},
								{ 13.0, 14.0, 15.0, 16.0},
								{ 17.0, 18.0, 19.0, 20.0} };

	const matrix<double> matA(&dma[0][0], 4, 4);
	const matrix<double> matB(&dmb[0][0], 3, 3);
	const matrix<double> matC(&dmc[0][0], 5, 4);

	mRank = MatrixRank(matA);
	cout << endl << "The rank of matA is : " << mRank << endl;

	mRank = MatrixRank(matB);
	cout << endl << "The rank of matB is : " << mRank << endl;

	mRank = MatrixRank(matC);
	cout << endl << "The rank of matC is : " << mRank << endl;
}

