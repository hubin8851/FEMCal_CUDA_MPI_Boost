//MatrixQR.cpp  
//用豪斯荷尔德(Householder)变换对一般m*n阶的实矩阵进行QR分解

#include <iostream>		//输入输出流
#include "Matrix.h"		//矩阵类及相关函数等的定义
using namespace std;	//名字空间

void main()				// 定义控制台应用程序的入口点
{
	const double dma[4][3] = 
	{
		{ 1.0,  1.0, -1.0},
		{ 2.0,  1.0,  0.0},
		{ 1.0, -1.0,  0.0},
		{-1.0,  2.0,  1.0}
	};

	matrix<double> matA(&dma[0][0], 4, 3);	//生成A阵
	matrix<double> matQ(4, 4);				//定义存放Q阵的矩阵对象

	cout << "matA : " << endl;
	MatrixLinePrint(matA);					//按行输出矩阵matA
	
	cout << endl << "matQ : " << endl;
	
	if(MatrixQR(matA, matQ))				//生成Q阵，且R阵在matA右上三角
	MatrixLinePrint(matQ);					//按行输出矩阵matQ

	cout << endl << "matR : " << endl;
	
	MatrixLinePrint(matA);					//按行输出矩阵matR
}
