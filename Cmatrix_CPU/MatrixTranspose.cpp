//MatrixTranspose.cpp	����ת����������

#include <iostream>		//���������ͷ�ļ�
#include "Comm.h"		//��������ͷ�ļ�
#include "Matrix.h"		//�����༰��غ����ȵĶ���
using namespace std;	//���ֿռ�

// �������̨Ӧ�ó������ڵ�
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
	MatrixLinePrint(ma);			//�������ma

	matrix<double> mb(3,5);			//����һ����Ϊ���maת��

	MatrixTranspose(ma, mb);		//����ma��ת����mb

	int sR = mb.GetRowNum();		//ȡmb������

	cout << endl << "Matrix mb is : " << endl;
	for(int i=0; i<sR; i++)
		MatrixLinePrint(mb, i);		//һ��һ���������mb
}

