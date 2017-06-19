//MatrixTranspose.cpp	矩阵转置与矩阵输出

#include <iostream>		//输入输出流头文件
#include "Comm.h"		//公共部分头文件
#include "Matrix.h"		//矩阵类及相关函数等的定义
using namespace std;	//名字空间

// 定义控制台应用程序的入口点
void main()
{
	double a[5][3] = 
	{
		{1.01, 1.02, 1.03},
		{2.04, 2.05, 2.06},
		{3.07, 3.08, 3.09},
		{4.01, 4.04, 4.08},
		{5.02, 5.05, 5.07}
	};

	matrix<double> ma(&a[0][0],5,3);

	cout << "Matrix ma is : " << endl;
	MatrixLinePrint(ma);			//输出矩阵ma

	matrix<double> mb(3,5);			//生成一新阵，为存放ma转置

	MatrixTranspose(ma, mb);		//生成ma和转置阵mb

	int sR = mb.GetRowNum();		//取mb的行数

	cout << endl << "Matrix mb is : " << endl;
	for(int i=0; i<sR; i++)
		MatrixLinePrint(mb, i);		//一行一行输出矩阵mb
}

