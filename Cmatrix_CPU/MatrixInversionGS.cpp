// MatrixInversionGS.cpp  
//ȫѡ��Ԫ��˹-Լ��(Gauss-Jordan)���������

#include <iostream>		//���������
#include "Matrix.h"		//�����༰��غ����ȵĶ���
using namespace std;	//���ֿռ�

void main()				// �������̨Ӧ�ó������ڵ�
{
	const double dma[4][4] =
	{
		{0.2368,0.2471,0.2568,1.2671},
		{1.1161,0.1254,0.1397,0.1490},
		{0.1582,1.1675,0.1768,0.1871},
		{0.1968,0.2071,1.2168,0.2271} 
	};

	const double dmc[4][4] = {	{  3.0, -3.0, -2.0,  4.0 },
								{  5.0, -5.0,  1.0,  8.0 },
								{ 11.0,  8.0,  5.0, -7.0 },
								{  5.0, -1.0, -3.0, -1.0 } };

	matrix<double> matA(&dma[0][0], 4, 4);
	matrix<double> matB(matA);
		
	cout << "matA : " << endl;
	MatrixLinePrint(matA);			//�����������matA

	if(MatrixInversionGS(matA)>0)	//�������溯�����ж��Ƿ�>0
	{
		cout << endl << "Inversion(matA) : " << endl;
	
		MatrixLinePrint(matA);		//�����������matA����
	
		cout << endl << "Inversion(matA) * matA : " << endl;

		matA = matA * matB;
		MatrixLinePrint(matA);		//�����������matA������ĳ˻�
	}

	const complex<double> cdma[4][4] = 
	{
		{
			complex<double>(0.2368,0.1345), complex<double>(0.2471,0.1678),
			complex<double>(0.2568,0.1825), complex<double>(1.2671,1.1161)
		},
		{
			complex<double>(1.1161,1.2671), complex<double>(0.1254,0.2617),
			complex<double>(0.1397,0.7024), complex<double>(0.1490,0.2721)
		},
		{
			complex<double>(0.1582,-0.2836), complex<double>(1.1675,-1.1967),
			complex<double>(0.1768,0.3558), complex<double>(0.1871,-0.2078)
		},
		{
			complex<double>(0.1968,0.3567), complex<double>(0.2071,-1.2345),
			complex<double>(1.2168,2.1185), complex<double>(0.2271,0.4773)
		}
	};

	matrix<complex<double> > cmatA(&cdma[0][0], 4, 4);
	matrix<complex<double> > cmatB(cmatA);
	
	cout << endl << "cmatA : " << endl;
	MatrixLinePrint(cmatA);			//�����������cmatA

	if(MatrixInversionGS(cmatA)>0)	//�жϷ���ֵ�Ƿ�>0
	{
		cout << endl <<  "Inversion(cmatA) : " << endl;
		#ifdef FileName
			ofsm << endl << "Inversion(cmatA) : " << endl;
		#endif
	
		MatrixLinePrint(cmatA);		//�����������cmatA����cmaBT
	
		cout << endl << "Inversion(cmatA) * cmatA : " << endl;
		#ifdef FileName
			ofsm << endl << "Inversion(cmatA) * cmatA : " << endl;
		#endif
		cmatA = cmatA * cmatB;
		MatrixLinePrint(cmatA);		//�����������cmatA������ĳ˻�
	}
}

