//MatrixSymmetry.cpp	�жϾ���Գ�

#include <iostream>		//���������
#include "Matrix.h"		//�����༰��غ����ȵĶ���
using namespace std;	//���ֿռ�

void main()				// �������̨Ӧ�ó������ڵ�
{
	const double dma[4][4] = 
	{
		{  1.0,  2.0,  3.0,  4.0 },
		{  5.0,  6.0,  7.0,  8.0 },
		{  9.0, 10.0, 11.0, 12.0 },
		{ 13.0, 14.0, 15.0, 16.0 }
	};

	const int dmb[3][3] = 
	{
		{  3, -3, -2 },
		{ -3,  8,  4 },
		{ -2,  4, -3 } 
	};

	const double dmc[2][3] = 
	{
		{  3.0, -3.0, -2.0 },
		{ -2.0,  4.0, -3.0 } 
	};

	const matrix<double> matA(&dma[0][0], 4, 4);
	const matrix<int> matB(&dmb[0][0], 3, 3);
	const matrix<double> matC(&dmc[0][0], 2, 3);


	if(MatrixSymmetry(matA))
		cout << "The matA is Symmetry." << endl;
	else
		cout << "The matA is not Symmetry." << endl;

	if(MatrixSymmetry(matB))
		cout << "The matB is Symmetry." << endl;
	else
		cout << "The matB is not Symmetry." << endl;

	if(MatrixSymmetry(matC))
		cout << "The matC is Symmetry." << endl;
	else
		cout << "The matC is not Symmetry." << endl;

}

