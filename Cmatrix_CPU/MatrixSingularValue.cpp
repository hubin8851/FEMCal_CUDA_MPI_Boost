//MatrixSingularValue.cpp	一般实矩阵的奇异值分解

#include <iostream>			//输入输出流
#include "Matrix.h"			//矩阵类及相关函数等的定义
using namespace std;		//名字空间

void main()					// 定义控制台应用程序的入口点
{
    double a[4][3]={ {1.0,1.0,-1.0},{2.0,1.0,0.0},
                           {1.0,-1.0,0.0},{-1.0,2.0,1.0}};
    double b[3][4]={ {1.0,1.0,-1.0,-1.0},{2.0,1.0,0.0,2.0},
								{1.0,-1.0,0.0,1.0}};

	matrix<double> ma(&a[0][0], 4, 3);
	matrix<double> mb(&b[0][0], 3, 4);	
	matrix<double> mu(4, 4), mv(3, 3), mc(4, 3), md(3, 4);

	double eps(DOUBLEERROR);

	int i = MatrixSingularValue(ma,mu,mv,eps);

	cout << "EXAMPLE(1) " << endl;
	cout << endl << "i = " << i << endl;
	cout << endl << "MAT U IS: " << endl;
	MatrixLinePrint(mu);

	cout << endl << "MAT V IS: " << endl;
	MatrixLinePrint(mv);

	cout << endl << "MAT A IS: " << endl;
	MatrixLinePrint(ma);

	cout << endl << "MAT UAV IS: " << endl;
	mc = mu * ma;
	ma = mc * mv;
	MatrixLinePrint(ma);

	cout << endl << "EXAMPLE(2)" << endl << endl;
	i = MatrixSingularValue(mb,mv,mu,eps);
	cout << "i = " << i << endl << endl;

	cout << "MAT U IS: " << endl;
	MatrixLinePrint(mv);
	
	cout << endl << "MAT V IS: " << endl;
	MatrixLinePrint(mu);

	cout << endl << "MAT B IS: " << endl;
	MatrixLinePrint(mb);

	cout << endl << "MAT UBV IS: " << endl;
	md = mv * mb;
	mb = md * mu;
	MatrixLinePrint(mb);
	cout << endl;
}

