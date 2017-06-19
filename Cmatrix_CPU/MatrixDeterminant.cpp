//MatrixDeterminant.cpp  ��˹ȫѡ��Ԫ������������ʽ

//#include <fstream>	//�ļ���ͷ�ļ�
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

	const double dmb[4][4] = 
	{
		{  3.8, -3.0, -2.0,  4.0 },
		{  5.0, -5.2,  1.0,  8.0 },
		{ 11.0,  8.0,  5.6, -7.0 },
		{  5.0, -1.0, -3.0, -1.4 } 
	};

	const matrix<double> matA(&dma[0][0], 4, 4);
	const matrix<double> matB(&dmb[0][0], 4, 4);

	cout << "Det(matA) = " << MatrixDeterminant(matA) << endl;

	cout << endl << "Det(matB) = " << MatrixDeterminant(matB) << endl;
	
	const int imc[3][3] = 
	{
		{-1, 2, 3}, 
		{ 9, 0, 2},
		{ 5, 1, 6}
	};

	const matrix<int> matC(&imc[0][0], 3, 3);
	cout << endl << "Det(matC) = " << MatrixDeterminant(matC) << endl;
}

