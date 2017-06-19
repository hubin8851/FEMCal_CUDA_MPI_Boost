//GeneralizedInversionSingularValue.cpp  
//广义逆的奇异值分解

#include <iostream>		//输入输出流
#include "Matrix.h"		//矩阵类及相关函数等的定义
using namespace std;	//名字空间

void main()				// 定义控制台应用程序的入口点
{
    double a[5][4]={ {1.0,2.0,3.0,4.0},
                  {6.0,7.0,8.0,9.0},{1.0,2.0,13.0,0.0},
                  {16.0,17.0,8.0,9.0},{2.0,4.0,3.0,4.0}};

	matrix<double> ma(&a[0][0], 5, 4);
	matrix<double> aa(4, 5), mc(5, 4), mu(5, 5), mv(4, 4);

	double eps(DOUBLEERROR);
	int m(5), n(4), ka(6);

	cout << "MAT A IS: " << endl;
	MatrixLinePrint(ma);

	cout << endl << "MAT A+ IS: " << endl;
	int i = GeneralizedInversionSingularValue(ma,aa,eps,mu,mv);
	MatrixLinePrint(aa);

	cout << endl << "MAT A++ IS: " << endl;
	i = GeneralizedInversionSingularValue(aa,mc,eps,mv,mu);
	MatrixLinePrint(mc);
	cout << endl;
}

