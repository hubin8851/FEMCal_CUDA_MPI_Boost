// MatrixRank.cpp		��ȫѡ��Ԫ��˹��Ԫ����һ��������

//#include <fstream>	//�ļ���ͷ�ļ�
#include <iostream>		//���������
#include "Matrix.h"		//�����༰��غ����ȵĶ���
using namespace std;	//���ֿռ�

void main()				// �������̨Ӧ�ó������ڵ�
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

