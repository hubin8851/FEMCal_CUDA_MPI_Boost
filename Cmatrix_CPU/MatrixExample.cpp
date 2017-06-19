//MatrixExample.cpp		�������и��ַ�����Ӧ��ʾ��

//#include <fstream>	//�ļ���ͷ�ļ�
#include <iostream>		//���������
#include "Matrix.h"		//������ͷ�ļ�
#include <complex>		//����ģ����ͷ�ļ�
using namespace std;	//���ֿռ�

void main()				// �������̨Ӧ�ó������ڵ�
{
	//ʹ�ù��캯��һ�������������
	matrix<int> matS(2, 1);		//��������2*1����matS
	matrix<int> matR(1, 2);		//��������1*2����matR
	
	//����Ԫ�ظ�ֵ
	matS(0, 0) = 1;
	matS(1, 0) = 2;
	matR(0, 0) = 3;
	matR(0, 1) = 4;

	//���������
	cout << endl << "matS : " << endl;
	
	MatrixLinePrint(matS);			//�����������matS
	cout << endl << "matR : " << endl;
	
	MatrixLinePrint(matR);			//�����������matR

	const double dma[4][4] = {	{  1.2,  2.6,  3.7,  4.8 },
								{  5.3,  6.0,  7.0,  8.0 },
								{  9.4, 10.0, 11.6, 12.0 },
								{  3.5, 14.0, 15.0, 16.8 }};

	const double dmb[4][4] = {	{  3.0, -3.0, -2.0,  4.0 },
								{  5.0, -5.0,  1.0,  8.0 },
								{  1.0,  8.0,  5.0, -7.0 },
								{  5.0, -1.0, -3.0, -1.0 } };

	//ʹ�ù��캯�����������������matA, matB
	const matrix<double> matA(&dma[0][0], 4, 4);
	const matrix<double> matB(&dmb[0][0], 4, 4);
	
	//ʹ�ù��캯�����������������matC
	matrix<double> matC(matA);		
	
	//���������matA, matb
	cout << endl << "matA : " << endl;
	MatrixLinePrint(matA);			//�����������matA

	cout << endl << "matB : " << endl;
	MatrixLinePrint(matB);			//�����������matB

	cout << endl << "matC : " << endl;
	MatrixLinePrint(matC);			//�����������matC

	cout << endl << "Row of matA : " << matA.GetRowNum();	//matA����
	cout << "\t Column of matA : " << matA.GetRowNum() << endl;	//matA����
	
	matC += matA;					//�����ԼӾ���
	cout << endl << "matC += matA : " << endl;
	MatrixLinePrint(matC);			//�����������matC

	matC *= matA;					//�����Գ˾���
	cout << endl << "matC  *= matA : " << endl;
	MatrixLinePrint(matC);			//�����������matC

	matC -= 12.3;					//������-��
	cout << endl << "matC -= 12.3 : " << endl;
	MatrixLinePrint(matC);			//�����������matC

	matC = -matC;					//����-��
	cout << endl << "matC = -matC : " << endl;
	MatrixLinePrint(matC);			//�����������matC

	matC = matA * matB;					//����˾���
	cout << endl << "matC = matA * matB : " << endl;
	MatrixLinePrint(matC);			//�����������matC

	matC = matC / 2.0;					//���������
	cout << endl << "matC = matC / 2.0 : " << endl;
	MatrixLinePrint(matC);			//�����������matC

	matC = 0.5 * matC;					//�����Ծ���
	cout << endl << "matC = 0.5 * matC : " << endl;
	MatrixLinePrint(matC);			//�����������matC

	matC = matA + matC;					//����Ӿ���
	cout << endl << "matC = matA + matC : " << endl;
	MatrixLinePrint(matC);			//�����������matC

	//�Ƚ��������Ƿ���ͬ����ͬ
	cout << endl << "if (matA == matC) : " << (matA == matC) << endl;	//�Ƿ�matA��matC��ͬ
	cout << endl << "if (matA != matB) : " << (matA != matB) << endl;	//�Ƿ�matA��matB����ͬ
	
	matrix< complex<float> > clhs(2, 1);	//���帴��������2*1����clhs
	matrix< complex<float> > crhs(1, 2);	//���帴��������1*2����crhs
	
	//�������;���Ԫ�ظ�
	clhs(0, 0) = complex<float> (1, 2);		//�������;���Ԫ�ظ�ֵ
	clhs(1, 0) = complex<float> (3, 4);
	crhs(0, 0) = complex<float> (5, 6);
	crhs(0, 1) = complex<float> (7, 8);
	
	cout << endl << "ComplexMatrix clhs :" << endl;
	MatrixLinePrint(clhs);		//���������clhs

	cout << endl << "ComplexMatrix crhs :" << endl;
	MatrixLinePrint(crhs);		//���������crhs

	matrix< complex<float> > cmm(clhs * crhs);	//cmm=clhs*crhs
	cout << endl << "ComplexMatrix * ComplexMatrix :" << endl;
	MatrixLinePrint(cmm);	//���������cmm=clhs*crhs

	matrix< complex<float> > cmp(cmm);		//����cmp����cmm��ʼ��
	cmp += cmm;				//cmp += cmm
	cout << endl << "ComplexMatrix * ComplexMatrix :" << endl;
	MatrixLinePrint(cmp);	//���������cmp+=cmm

}

