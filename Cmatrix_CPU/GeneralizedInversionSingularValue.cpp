//GeneralizedInversionSingularValue.cpp  
//�����������ֵ�ֽ�

#include <iostream>		//���������
#include "Matrix.h"		//�����༰��غ����ȵĶ���
using namespace std;	//���ֿռ�

void main()				// �������̨Ӧ�ó������ڵ�
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

