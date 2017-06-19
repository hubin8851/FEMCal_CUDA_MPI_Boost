// FEMCal_CUDA_MPI_Boost.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "HbxGloFunc.h"

#include "EigenCal.h"
#include "CGSMatCal.h"
#include "CGDMatCal.h"
#include "CGSMatCal_UM.h"
#include "BiCGStabD.h"
#include "PCGSMatCal.h"
#include "PCGDMatCal.h"
#include "SSOR_PCGMatCal.h"
#include "CntToMysql.h"
#include "ImpFrmMatlab.h"
#include "BaseLag.h"
#include "CuInterpolation.h"


#include "DataInput_f.h"
#include "DataInput_d.h"
#include "MatrixAssemble.h"
#include "XMLManage.h"

//���ݶ�ȡ���ȫ�ֱ���
HBXDef::PrecsType_t	g_PrecsType;	//ȫ�ֱ������������ݶ�ȡ�ͼ���ľ���(float or double)
int		g_Dim;	//ȫ��ά�ȣ���������2��3ά��Ͻ�ģ
int		g_NodeNumPerElement;//��Ԫ�ڽڵ���
int		g_SectionNum,g_NodeNum,g_ElementNum,g_MaterialNum;
int		g_SetNum;//�ڵ�͵�Ԫ��������
int		g_LoadingTractionNum,g_LoadingConcentrationNum,g_NodeRestrainNum;	//�غ����ȫ�ֱ���
int		g_RowPerElem;	//��Ԫ������� dim*g_NodeNumPerElement
int		g_nnA,g_nA;
int		g_MatRowNum = 900;//ȫ�־������
int		g_EleNumInMatrix;//�����е�Ԫ��
int		g_WeightMark;	//�Ƿ���������
int		g_PlaneStressOrStrain;	//ƽ��Ӧ��orӦ������

HBXFEMDef::NodeInMap_f		g_vNode_f;		//ȫ�ֽڵ�����,�мǣ���map��Ϊ����������
HBXFEMDef::_vBaseElem_f		g_vBaseElem_f;	//ȫ�ֻ���Ԫ����
HBXFEMDef::_vMaterial_f		g_vMaterial_f;	//ȫ�ֲ�������
HBXFEMDef::_vSection_f		g_vSection_f;	//ȫ�ֽ�������

HBXFEMDef::NodeInMap_d		g_vNode_d;		//ȫ�ֽڵ�����,�мǣ���map��Ϊ����������
HBXFEMDef::_vBaseElem_d		g_vBaseElem_d;	//ȫ�ֻ���Ԫ����
HBXFEMDef::_vMaterial_d		g_vMaterial_d;	//ȫ�ֲ�������
HBXFEMDef::_vSection_d		g_vSection_d;	//ȫ�ֽ�������

HBXFEMDef::_vSet			g_vNSet;//�ڵ㼯������
HBXFEMDef::_vSet			g_vElSet;//��Ԫ��������

//csr��ʽ����������
HBXDef::_CSRInput<double>	g_csrMat;

int		g_Iters = 20000;
float	g_Precision = 1.0e-6f;
const char		*g_ResultFile = (char*)"export_eig.dat";

boost::filesystem::path		g_UFgetPath("..\\data\\UFget");
boost::filesystem::path		g_InPath("..\\data\\source");

CCntToMysql		g_CntToMysql;	//�������ݿ�
CImpFrmMatlab	g_ImpFrmMat;

void	printhelp();	//����̨�ڴ�ӡ�����������������˳�
void	SeachFileTest(int argc, const char **argv);//boost�ļ������������ɹ��������ڷ��ص���������������static����������������ǰ����ջ�Ϸ���ռ䣬�޷�Ԥ���룬�޷����ذ�����
void	ReadElePropTest(int argc, const char **argv);
void	ReadAeroTableTest(int argc, const char **argv);

void	CGSCalTest(int argc, const char **argv);
void	CGDCalTest(int argc, const char **argv);
void	PCGSTest(int argc, const char **argv);
void	PCGDTest(int argc, const char **argv);
void	BiCGStabDTest( int argc, const char **argv );
void	SSOR_PCGDTest(int argc, const char **argv);
void	EigenTest(int argc, const char **argv);//ʹ��Eigen�ĺ���
void	MatReadTest(int argc, const char **argv);

int _tmain(int argc, const char **argv)
{
	if (checkCmdLineFlag(argc, argv, "help"))
	{
		printhelp();
		getchar();
		return 0;
	}

	printf( "argv[0] = %s", argv[0] );

	g_PrecsType = HBXDef::PrecsType_t::_PRECS_DOUBLE_;

	g_CntToMysql.Initial();
	if ( g_CntToMysql.ConnMySQL())//δ���ӳɹ�
	{
		std::cout<<"PCG���Ե�Mysql���ӳɹ�������..."<<std::endl;
	}
	else
	{
		std::cout<<"PCG���Ե�Mysql���Ӳ��ɹ�������..."<<std::endl;
	}

	ReadElePropTest(argc, argv);
	ReadAeroTableTest(argc, argv);

	MatReadTest(argc, argv);
//	SSOR_PCGDTest(argc, argv);
	for( int _iter=0; _iter<1; _iter++)
	{
		Sleep(200);
		EigenTest(argc, argv);
		//CGSCalTest(argc, argv);
		CGDCalTest(argc, argv);
		//PCGSTest(argc, argv);
		//BiCGStabDTest(argc, argv);
		//PCGDTest(argc, argv);
		g_MatRowNum *=4;
 	}
	std::cout<<"�����˳�..."<<std::endl;
	getchar();
	return 0;
}

void ReadElePropTest(int argc, const char **argv)
{
	CXMLManage	g_Xml;//���ڶ�ȡ��Ԫ���Ե�xml
//	g_Xml.SetSourceFilePath( "Trans_Param.xml", "F:\\data from HBX_phd\\VS2010\\boost_sharedMemTestV2\\data");
	g_Xml.SetSourceFilePath();
	g_Xml.ReadElemtProperty();
	HBXFEMDef::_EltPtyInMap g_ElemtPropInMap = g_Xml.GetElemtPtyInMap();

//	˫����double����
	CDataInput_d		g_DataIn_d( g_ElemtPropInMap );
	//g_DataIn_d.SetSourceFilePath("3D8NodeTest.inp", "D:\\AbaqusTemp");	//��ȡinp�ļ�����...
	g_DataIn_d.SetSourceFilePath("C3D8RTest_320E.inp", "F:\\data from HBX_phd\\database\\C3D8RTest");	//��ȡmtx�ļ�����...
	g_DataIn_d.SetInputData();
	MatrixAssemble	g_MatAssemble_d;
	g_MatAssemble_d.AssembleAnalysis(HBXDef::PrecsType_t::_PRECS_DOUBLE_);
	g_MatAssemble_d.ExportResultFile("..\\data\\source");
}

void ReadAeroTableTest(int argc, const char **argv)
{
	using namespace HBXDef;
	CXMLManage	g_Xml;		//���ڶ�ȡ��Ԫ���Ե�xml
	g_Xml.ReadAeroCoef("F:\\data from HBX_phd\\database\\x33��������\\CA_da.xml");

	const unsigned int g_T = 3;
	typedef	HBXDef::CBaseLag<g_T>	_samelag;
	typedef HBXDef::CuInterpolation<g_T, HBXDef::CudaMalloc_t::NORMAL>	_sameCulag;

	_samelag	g_BaseLag( &g_Xml.GetAeroTable("CA_da"), 0 );	//���Բ�ֵ
	HBXDef::CheckUserDefErrors( g_BaseLag.SetBlkId(0) );
	g_BaseLag.Initial();
	HBXDef::UserDefFloat  g_dataIn[] = {-1.5, 0.2, -30.0, -5.0};
	g_BaseLag.ReadTableData(g_dataIn);

	_sameCulag*	g_culag = nullptr;
	checkCudaErrors( cudaMalloc(&g_culag, MAX_GPUITER* sizeof(_sameCulag)) );
	for (int i = 0; i < MAX_GPUITER; i++)
	{
		g_culag[i].Initial(g_BaseLag);
	}
	//_sameCulag	g_CuInterpolation( &g_Xml.GetAeroTable("CA_da"), 0 );
	//g_CuInterpolation.Initial(MAX_GPUITER);
	//HBXDef::_USERDEFT*	_pInput = nullptr;
	//cutestlag1<g_T, HBXDef::CudaMalloc_t::NORMAL><< < GRIDSIZE, BLOCKSIZE, MAX_SHARED, 0 > >>(g_culag);
}

//�ò������������ʵ��ϡ�����
void MatReadTest(int argc, const char **argv)
{
	HBXDef::CheckUserDefErrors ( g_ImpFrmMat.ImpFrmMatlab("UF_Index.mat", "J:\\SuiteSparse-2.1.1\\UFget\\mat") );
	g_csrMat = g_ImpFrmMat.GetStiffMat( true );
	g_csrMat.h_rhs = g_ImpFrmMat.GetRhsVec( true );
	g_csrMat._nnzA = g_ImpFrmMat.GetNoneZeroLgt();
	g_csrMat._nA = g_ImpFrmMat.GetnALgt();
}

void CGSCalTest(int argc, const char **argv)
{
	std::cout<<"��ʼ�����ȹ����ݶȷ��ļ���..."<<std::endl;
	CGSMatCal	g_CGConj(g_MatRowNum);
//	g_CGConj.runtest();
//	g_CGConj.genTridiag(g_MatRowNum, true, "..\\data\\check");
//	g_CGConj.genRhs(g_MatRowNum, true, "..\\data\\check");
	//ͨ��������ȡ
	if ( g_csrMat.h_NoneZeroVal && g_csrMat.h_iColSort && g_csrMat.h_iNonZeroRowSort )
	{
		g_CGConj.SetStiffMat( g_csrMat.h_NoneZeroVal, g_csrMat.h_iColSort, g_csrMat.h_iNonZeroRowSort, g_csrMat._nnzA, g_csrMat._nA, true );
	}
	else//���ļ��ж�ȡ����
	{
		if ( HBXDef::DataAloc_t::INVALID == g_CGConj.SetStiffMat(	"NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath ) )
		{
			std::cout<<"ϵ�������п�ָ��..."<<std::endl;
		}
		//g_EigenCal.genTridiag( g_csrMat._nA-1, true, "\\data\\Eigen" );
	}
	if (g_csrMat.h_rhs)
	{
		g_CGConj.SetLoadVec( g_csrMat.h_rhs, (g_csrMat._nA-1),true );
	}
	else//���ļ��ж�ȡ����
	{
		if ( HBXDef::DataAloc_t::INVALID == g_CGConj.SetLoadVec( "LoadVec.data", g_InPath ) )
		{
			std::cout<<"�Ҷ������������..."<<std::endl;
			g_CGConj.genRhs();
			std::cout<<"�������ɵ�ʽ�Ҷ�..."<<std::endl;
		}
		//g_EigenCal.genRhs( g_csrMat._nA-1, true, "\\data\\Eigen" );
	}
	g_CGConj.DeviceMalloc();
	g_CGConj.GetGPUDeviceQuery();
	g_CGConj.MemCpy(HBXDef::HostToDevice);
	if (true == g_CGConj.InitialDescr())
	{
		std::cout<<"\n��ʼ��CG����������������ɹ���\n"<<std::endl;
	}
	float _UsedTime = g_CGConj.ConjugateWithGPU();
	g_CGConj.ExportEigToFile("cg_f.data");
	//g_CntToMysql.InsertLog("cg_tridiag", g_MatRowNum, g_CGConj.GetIters(), _UsedTime, g_CGConj.Getqaerr()/10 );
}

void CGDCalTest(int argc, const char **argv)
{
	std::cout<<"��ʼ˫���ȹ����ݶȷ��ļ���..."<<std::endl;
	CGDMatCal	g_CGDConj(g_MatRowNum);
	//	g_CGConj.runtest();
	//g_CGDConj.genTridiag(g_MatRowNum, true, "..\\data\\check");
	//g_CGDConj.genRhs(g_MatRowNum, true, "..\\data\\check");
	if ( g_csrMat.h_NoneZeroVal && g_csrMat.h_iColSort && g_csrMat.h_iNonZeroRowSort )
	{
		g_CGDConj.SetStiffMat( g_csrMat.h_NoneZeroVal, g_csrMat.h_iColSort, g_csrMat.h_iNonZeroRowSort, g_csrMat._nnzA, g_csrMat._nA, true );
	}
	else//���ļ��ж�ȡ����
	{
		if ( HBXDef::DataAloc_t::INVALID == g_CGDConj.SetStiffMat(	"NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath ) )
		{
			std::cout<<"ϵ�������п�ָ��..."<<std::endl;
		}
		//g_EigenCal.genTridiag( g_csrMat._nA-1, true, "\\data\\Eigen" );
	}
	g_CGDConj.SetLoadVec("LoadVec.data", g_InPath);
	g_CGDConj.DeviceMalloc();
	g_CGDConj.GetGPUDeviceQuery();
	g_CGDConj.MemCpy(HBXDef::HostToDevice);
	if (true == g_CGDConj.InitialDescr())
	{
		std::cout<<"\n��ʼ��CG����������������ɹ���\n"<<std::endl;
	}
	float _UsedTime = g_CGDConj.ConjugateWithGPU();
	g_CGDConj.ExportEigToFile("cg_d.data", "..\\data\\export");
	//g_CntToMysql.InsertLog("cg_tridiag", g_MatRowNum, g_CGDConj.GetIters(), _UsedTime, g_CGDConj.Getqaerr()/10 );
}

void	BiCGStabDTest( int argc, const char **argv )
{
	std::cout<<"\n��ʼ˫����Ԥ�������ݶȷ��ļ���..."<<std::endl;
	BiCGStabD	g_BiCGStabD( g_MatRowNum );
	HBXDef::DataAloc_t _state1 = g_BiCGStabD.SetStiffMatFromMTX();
	//HBXDef::DataAloc_t _state2 = g_BiCGStabD.SetStiffMat( "NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath );
	g_BiCGStabD.SetLoadVec("x.data",g_InPath); 
	g_BiCGStabD.DeviceMalloc();
	g_BiCGStabD.GetGPUDeviceQuery();
	g_BiCGStabD.MemCpy(HBXDef::HostToDevice);
	if (true == g_BiCGStabD.InitialDescr())
	{
		std::cout<<"\n��ʼ��CG����������������ɹ���\n"<<std::endl;
	}
	float _UsedTime = g_BiCGStabD.ConjugateWithGPU( g_Precision, g_Iters );
}

void	PCGSTest(int argc, const char **argv)
{
	std::cout<<"��ʼ������Ԥ�������ݶȷ��ļ���..."<<std::endl;
	PCGSMatCal	g_PCGConj(g_MatRowNum);
//	g_PCGConj.genLaplace(g_MatRowNum, true, "..\\data\\check");//��������
//	g_PCGConj.genRhs(g_MatRowNum, true, "..\\data\\check");
	HBXDef::DataAloc_t _state = g_PCGConj.SetStiffMat("NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath );
	g_PCGConj.SetLoadVec("LoadVec.data",g_InPath);
	g_PCGConj.DeviceMalloc();
	g_PCGConj.GetGPUDeviceQuery();
	g_PCGConj.MemCpy(HBXDef::HostToDevice);
	if (true == g_PCGConj.InitialDescr())
	{
		std::cout<<"\n��ʼ��CG����������������ɹ���\n"<<std::endl;
	}
	float _UsedTime = g_PCGConj.ConjugateWithGPU();
	g_PCGConj.ExportEigToFile("pcg.data", "..\\data\\export");
	//g_CntToMysql.InsertLog("pcg_tridiag", g_MatRowNum, g_PCGConj.GetIters(), _UsedTime, g_PCGConj.Getqaerr()/10 );
}

void	PCGDTest(int argc, const char **argv)
{
	std::cout<<"��ʼ˫����Ԥ�������ݶȷ��ļ���..."<<std::endl;
	PCGDMatCal	g_PCGConj(g_MatRowNum);
	//g_PCGConj.genLaplace(g_MatRowNum, true, "..\\data\\check");//��������
	//g_PCGConj.genRhs(g_MatRowNum, true, "..\\data\\check");
	HBXDef::DataAloc_t _state = g_PCGConj.SetStiffMat("NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath );
	g_PCGConj.SetLoadVec("LoadVec.data",g_InPath); 
	g_PCGConj.DeviceMalloc();
	g_PCGConj.GetGPUDeviceQuery();
	g_PCGConj.MemCpy(HBXDef::HostToDevice);
	if (true == g_PCGConj.InitialDescr())
	{
		std::cout<<"\n��ʼ��CG����������������ɹ���\n"<<std::endl;
	}
	float _UsedTime = g_PCGConj.ConjugateWithGPU( g_Precision, g_Iters );
	g_PCGConj.ExportEigToFile("pcg.data", "..\\data\\export");
	//g_CntToMysql.InsertLog("pcg_tridiag", g_MatRowNum, g_PCGConj.GetIters(), _UsedTime, g_PCGConj.Getqaerr()/10 );
}

void	SSOR_PCGDTest(int argc, const char **argv)
{
//	SSOR_PCGMatCal	g_SSOR_PCGConj(g_MatRowNum);

//	HBXDef::DataAloc_t _state = g_SSOR_PCGConj.SetStiffMat("NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", "..\\data\\UFget");
//	g_SSOR_PCGConj.SetLoadVec("LoadVec.data","..\\data\\UFget");
//	g_SSOR_PCGConj.HostMalloc(false);
//	g_SSOR_PCGConj.
}

void	EigenTest(int argc, const char **argv)
{
	EigenCal	g_EigenCal(g_MatRowNum);
	
	//ͨ��������ȡ
	if ( g_csrMat.h_NoneZeroVal && g_csrMat.h_iColSort && g_csrMat.h_iNonZeroRowSort )
	{
		g_EigenCal.SetStiffMat( g_csrMat.h_NoneZeroVal, g_csrMat.h_iColSort, g_csrMat.h_iNonZeroRowSort, g_csrMat._nnzA, g_csrMat._nA, true );
	}
	else//���ļ��ж�ȡ����
	{
		if ( HBXDef::DataAloc_t::INVALID == g_EigenCal.SetStiffMat(	"NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath ) )
		{
			std::cout<<"ϵ�������п�ָ��..."<<std::endl;
		}
		//g_EigenCal.genTridiag( g_csrMat._nA-1, true, "\\data\\Eigen" );
	}
	if (g_csrMat.h_rhs)
	{
		g_EigenCal.SetLoadVec( g_csrMat.h_rhs, (g_csrMat._nA-1),true );
	}
	else//���ļ��ж�ȡ����
	{
		if ( HBXDef::DataAloc_t::INVALID == g_EigenCal.SetLoadVec( "LoadVec.data", g_InPath ) )
		{
			std::cout<<"�Ҷ������������..."<<std::endl;
			g_EigenCal.genRhs();
			std::cout<<"�������ɵ�ʽ�Ҷ�..."<<std::endl;
		}
		//g_EigenCal.genRhs( g_csrMat._nA-1, true, "\\data\\Eigen" );
	}
	using namespace Eigen;
	g_EigenCal.EigenCalc( HBXDef::SolveMethod_t::CG_Mtd );
//	g_EigenCal.EigenCalc( HBXDef::SolveMethod_t::BICGSTAB_Mtd );
	g_EigenCal.ExportEigToFile("Eigen_EigVal.data", "..\\data\\export" );
}



void SeachFileTest(int argc, const char **argv)
{
	boost::optional<boost::filesystem::path> r = find_file("F:/","NonZeroVal.data");
	if (r)
	{
		std::cout<<*r<<std::endl;
	}
}
//
void	printhelp()
{
	printf("������CG����PCG������ֵ[OPTION]...\n");
	printf("Test the use of time\n");
	printf("Example:  CG \ PCG����������ά����512-10240֮�䣬��256����������\n");
	printf("./FEMCal_CUDA --matrix-size=1024 --iters=1000 --precision=1e-5 --ResultFile=export_eig.dat increment=256 \n");
	printf("\n");
	printf("Options:\n");
	printf("--help\tDisplay this help menu\n");
	printf("--matrix-size\t��ʼ����ά�ȣ�Ĭ��Ϊ512\n");
	printf("--iters\tÿ�μ���������ֵ����ĵ����������ޣ�Ĭ��Ϊ1000\n");
	printf("--precision\t����Ĭ�����ȣ�Ĭ��Ϊ1e-5\n");
	printf("--increment\t�����ܲ���ʱÿ�ξ���ά�ȵ�����ֵ��Ĭ��Ϊ������\n");
	printf("--ResultFile\t�������ֵ�ļ�����·��Ϊ..//data//export//�£�Ĭ��Ϊexport_eig.dat\n");
}