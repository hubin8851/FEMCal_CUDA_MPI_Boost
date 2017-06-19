// FEMCal_CUDA_MPI_Boost.cpp : 定义控制台应用程序的入口点。
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

//数据读取相关全局变量
HBXDef::PrecsType_t	g_PrecsType;	//全局变量，设置数据读取和计算的精度(float or double)
int		g_Dim;	//全局维度，不适用于2、3维混合建模
int		g_NodeNumPerElement;//单元内节点数
int		g_SectionNum,g_NodeNum,g_ElementNum,g_MaterialNum;
int		g_SetNum;//节点和单元集合总数
int		g_LoadingTractionNum,g_LoadingConcentrationNum,g_NodeRestrainNum;	//载荷相关全局变量
int		g_RowPerElem;	//单元矩阵阶数 dim*g_NodeNumPerElement
int		g_nnA,g_nA;
int		g_MatRowNum = 900;//全局矩阵阶数
int		g_EleNumInMatrix;//矩阵中单元数
int		g_WeightMark;	//是否引入重力
int		g_PlaneStressOrStrain;	//平面应力or应变问题

HBXFEMDef::NodeInMap_f		g_vNode_f;		//全局节点容器,切记！是map！为了搜索更快
HBXFEMDef::_vBaseElem_f		g_vBaseElem_f;	//全局基单元容器
HBXFEMDef::_vMaterial_f		g_vMaterial_f;	//全局材料容器
HBXFEMDef::_vSection_f		g_vSection_f;	//全局截面容器

HBXFEMDef::NodeInMap_d		g_vNode_d;		//全局节点容器,切记！是map！为了搜索更快
HBXFEMDef::_vBaseElem_d		g_vBaseElem_d;	//全局基单元容器
HBXFEMDef::_vMaterial_d		g_vMaterial_d;	//全局材料容器
HBXFEMDef::_vSection_d		g_vSection_d;	//全局截面容器

HBXFEMDef::_vSet			g_vNSet;//节点集合容器
HBXFEMDef::_vSet			g_vElSet;//单元集合容器

//csr格式的三个数组
HBXDef::_CSRInput<double>	g_csrMat;

int		g_Iters = 20000;
float	g_Precision = 1.0e-6f;
const char		*g_ResultFile = (char*)"export_eig.dat";

boost::filesystem::path		g_UFgetPath("..\\data\\UFget");
boost::filesystem::path		g_InPath("..\\data\\source");

CCntToMysql		g_CntToMysql;	//连接数据库
CImpFrmMatlab	g_ImpFrmMat;

void	printhelp();	//控制台内打印帮助，并按键返回退出
void	SeachFileTest(int argc, const char **argv);//boost文件搜索用例，成功。但鉴于返回的是容器，类似于static数据在主函数运行前即在栈上分配空间，无法预编译，无法多重包含。
void	ReadElePropTest(int argc, const char **argv);
void	ReadAeroTableTest(int argc, const char **argv);

void	CGSCalTest(int argc, const char **argv);
void	CGDCalTest(int argc, const char **argv);
void	PCGSTest(int argc, const char **argv);
void	PCGDTest(int argc, const char **argv);
void	BiCGStabDTest( int argc, const char **argv );
void	SSOR_PCGDTest(int argc, const char **argv);
void	EigenTest(int argc, const char **argv);//使用Eigen的函数
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
	if ( g_CntToMysql.ConnMySQL())//未连接成功
	{
		std::cout<<"PCG测试的Mysql连接成功，继续..."<<std::endl;
	}
	else
	{
		std::cout<<"PCG测试的Mysql连接不成功，请检查..."<<std::endl;
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
	std::cout<<"按键退出..."<<std::endl;
	getchar();
	return 0;
}

void ReadElePropTest(int argc, const char **argv)
{
	CXMLManage	g_Xml;//用于读取单元属性的xml
//	g_Xml.SetSourceFilePath( "Trans_Param.xml", "F:\\data from HBX_phd\\VS2010\\boost_sharedMemTestV2\\data");
	g_Xml.SetSourceFilePath();
	g_Xml.ReadElemtProperty();
	HBXFEMDef::_EltPtyInMap g_ElemtPropInMap = g_Xml.GetElemtPtyInMap();

//	双精度double测试
	CDataInput_d		g_DataIn_d( g_ElemtPropInMap );
	//g_DataIn_d.SetSourceFilePath("3D8NodeTest.inp", "D:\\AbaqusTemp");	//读取inp文件测试...
	g_DataIn_d.SetSourceFilePath("C3D8RTest_320E.inp", "F:\\data from HBX_phd\\database\\C3D8RTest");	//读取mtx文件测试...
	g_DataIn_d.SetInputData();
	MatrixAssemble	g_MatAssemble_d;
	g_MatAssemble_d.AssembleAnalysis(HBXDef::PrecsType_t::_PRECS_DOUBLE_);
	g_MatAssemble_d.ExportResultFile("..\\data\\source");
}

void ReadAeroTableTest(int argc, const char **argv)
{
	using namespace HBXDef;
	CXMLManage	g_Xml;		//用于读取单元属性的xml
	g_Xml.ReadAeroCoef("F:\\data from HBX_phd\\database\\x33气动数据\\CA_da.xml");

	const unsigned int g_T = 3;
	typedef	HBXDef::CBaseLag<g_T>	_samelag;
	typedef HBXDef::CuInterpolation<g_T, HBXDef::CudaMalloc_t::NORMAL>	_sameCulag;

	_samelag	g_BaseLag( &g_Xml.GetAeroTable("CA_da"), 0 );	//用以插值
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

//该测试用例仅针对实数稀疏矩阵
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
	std::cout<<"开始单精度共轭梯度法的计算..."<<std::endl;
	CGSMatCal	g_CGConj(g_MatRowNum);
//	g_CGConj.runtest();
//	g_CGConj.genTridiag(g_MatRowNum, true, "..\\data\\check");
//	g_CGConj.genRhs(g_MatRowNum, true, "..\\data\\check");
	//通过向量读取
	if ( g_csrMat.h_NoneZeroVal && g_csrMat.h_iColSort && g_csrMat.h_iNonZeroRowSort )
	{
		g_CGConj.SetStiffMat( g_csrMat.h_NoneZeroVal, g_csrMat.h_iColSort, g_csrMat.h_iNonZeroRowSort, g_csrMat._nnzA, g_csrMat._nA, true );
	}
	else//从文件中读取数据
	{
		if ( HBXDef::DataAloc_t::INVALID == g_CGConj.SetStiffMat(	"NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath ) )
		{
			std::cout<<"系数矩阵有空指针..."<<std::endl;
		}
		//g_EigenCal.genTridiag( g_csrMat._nA-1, true, "\\data\\Eigen" );
	}
	if (g_csrMat.h_rhs)
	{
		g_CGConj.SetLoadVec( g_csrMat.h_rhs, (g_csrMat._nA-1),true );
	}
	else//从文件中读取数据
	{
		if ( HBXDef::DataAloc_t::INVALID == g_CGConj.SetLoadVec( "LoadVec.data", g_InPath ) )
		{
			std::cout<<"右端向量输入错误..."<<std::endl;
			g_CGConj.genRhs();
			std::cout<<"重新生成等式右端..."<<std::endl;
		}
		//g_EigenCal.genRhs( g_csrMat._nA-1, true, "\\data\\Eigen" );
	}
	g_CGConj.DeviceMalloc();
	g_CGConj.GetGPUDeviceQuery();
	g_CGConj.MemCpy(HBXDef::HostToDevice);
	if (true == g_CGConj.InitialDescr())
	{
		std::cout<<"\n初始化CG法的描述器及句柄成功！\n"<<std::endl;
	}
	float _UsedTime = g_CGConj.ConjugateWithGPU();
	g_CGConj.ExportEigToFile("cg_f.data");
	//g_CntToMysql.InsertLog("cg_tridiag", g_MatRowNum, g_CGConj.GetIters(), _UsedTime, g_CGConj.Getqaerr()/10 );
}

void CGDCalTest(int argc, const char **argv)
{
	std::cout<<"开始双精度共轭梯度法的计算..."<<std::endl;
	CGDMatCal	g_CGDConj(g_MatRowNum);
	//	g_CGConj.runtest();
	//g_CGDConj.genTridiag(g_MatRowNum, true, "..\\data\\check");
	//g_CGDConj.genRhs(g_MatRowNum, true, "..\\data\\check");
	if ( g_csrMat.h_NoneZeroVal && g_csrMat.h_iColSort && g_csrMat.h_iNonZeroRowSort )
	{
		g_CGDConj.SetStiffMat( g_csrMat.h_NoneZeroVal, g_csrMat.h_iColSort, g_csrMat.h_iNonZeroRowSort, g_csrMat._nnzA, g_csrMat._nA, true );
	}
	else//从文件中读取数据
	{
		if ( HBXDef::DataAloc_t::INVALID == g_CGDConj.SetStiffMat(	"NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath ) )
		{
			std::cout<<"系数矩阵有空指针..."<<std::endl;
		}
		//g_EigenCal.genTridiag( g_csrMat._nA-1, true, "\\data\\Eigen" );
	}
	g_CGDConj.SetLoadVec("LoadVec.data", g_InPath);
	g_CGDConj.DeviceMalloc();
	g_CGDConj.GetGPUDeviceQuery();
	g_CGDConj.MemCpy(HBXDef::HostToDevice);
	if (true == g_CGDConj.InitialDescr())
	{
		std::cout<<"\n初始化CG法的描述器及句柄成功！\n"<<std::endl;
	}
	float _UsedTime = g_CGDConj.ConjugateWithGPU();
	g_CGDConj.ExportEigToFile("cg_d.data", "..\\data\\export");
	//g_CntToMysql.InsertLog("cg_tridiag", g_MatRowNum, g_CGDConj.GetIters(), _UsedTime, g_CGDConj.Getqaerr()/10 );
}

void	BiCGStabDTest( int argc, const char **argv )
{
	std::cout<<"\n开始双精度预处理共轭梯度法的计算..."<<std::endl;
	BiCGStabD	g_BiCGStabD( g_MatRowNum );
	HBXDef::DataAloc_t _state1 = g_BiCGStabD.SetStiffMatFromMTX();
	//HBXDef::DataAloc_t _state2 = g_BiCGStabD.SetStiffMat( "NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath );
	g_BiCGStabD.SetLoadVec("x.data",g_InPath); 
	g_BiCGStabD.DeviceMalloc();
	g_BiCGStabD.GetGPUDeviceQuery();
	g_BiCGStabD.MemCpy(HBXDef::HostToDevice);
	if (true == g_BiCGStabD.InitialDescr())
	{
		std::cout<<"\n初始化CG法的描述器及句柄成功！\n"<<std::endl;
	}
	float _UsedTime = g_BiCGStabD.ConjugateWithGPU( g_Precision, g_Iters );
}

void	PCGSTest(int argc, const char **argv)
{
	std::cout<<"开始单精度预处理共轭梯度法的计算..."<<std::endl;
	PCGSMatCal	g_PCGConj(g_MatRowNum);
//	g_PCGConj.genLaplace(g_MatRowNum, true, "..\\data\\check");//生成数组
//	g_PCGConj.genRhs(g_MatRowNum, true, "..\\data\\check");
	HBXDef::DataAloc_t _state = g_PCGConj.SetStiffMat("NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath );
	g_PCGConj.SetLoadVec("LoadVec.data",g_InPath);
	g_PCGConj.DeviceMalloc();
	g_PCGConj.GetGPUDeviceQuery();
	g_PCGConj.MemCpy(HBXDef::HostToDevice);
	if (true == g_PCGConj.InitialDescr())
	{
		std::cout<<"\n初始化CG法的描述器及句柄成功！\n"<<std::endl;
	}
	float _UsedTime = g_PCGConj.ConjugateWithGPU();
	g_PCGConj.ExportEigToFile("pcg.data", "..\\data\\export");
	//g_CntToMysql.InsertLog("pcg_tridiag", g_MatRowNum, g_PCGConj.GetIters(), _UsedTime, g_PCGConj.Getqaerr()/10 );
}

void	PCGDTest(int argc, const char **argv)
{
	std::cout<<"开始双精度预处理共轭梯度法的计算..."<<std::endl;
	PCGDMatCal	g_PCGConj(g_MatRowNum);
	//g_PCGConj.genLaplace(g_MatRowNum, true, "..\\data\\check");//生成数组
	//g_PCGConj.genRhs(g_MatRowNum, true, "..\\data\\check");
	HBXDef::DataAloc_t _state = g_PCGConj.SetStiffMat("NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath );
	g_PCGConj.SetLoadVec("LoadVec.data",g_InPath); 
	g_PCGConj.DeviceMalloc();
	g_PCGConj.GetGPUDeviceQuery();
	g_PCGConj.MemCpy(HBXDef::HostToDevice);
	if (true == g_PCGConj.InitialDescr())
	{
		std::cout<<"\n初始化CG法的描述器及句柄成功！\n"<<std::endl;
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
	
	//通过向量读取
	if ( g_csrMat.h_NoneZeroVal && g_csrMat.h_iColSort && g_csrMat.h_iNonZeroRowSort )
	{
		g_EigenCal.SetStiffMat( g_csrMat.h_NoneZeroVal, g_csrMat.h_iColSort, g_csrMat.h_iNonZeroRowSort, g_csrMat._nnzA, g_csrMat._nA, true );
	}
	else//从文件中读取数据
	{
		if ( HBXDef::DataAloc_t::INVALID == g_EigenCal.SetStiffMat(	"NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data", g_InPath ) )
		{
			std::cout<<"系数矩阵有空指针..."<<std::endl;
		}
		//g_EigenCal.genTridiag( g_csrMat._nA-1, true, "\\data\\Eigen" );
	}
	if (g_csrMat.h_rhs)
	{
		g_EigenCal.SetLoadVec( g_csrMat.h_rhs, (g_csrMat._nA-1),true );
	}
	else//从文件中读取数据
	{
		if ( HBXDef::DataAloc_t::INVALID == g_EigenCal.SetLoadVec( "LoadVec.data", g_InPath ) )
		{
			std::cout<<"右端向量输入错误..."<<std::endl;
			g_EigenCal.genRhs();
			std::cout<<"重新生成等式右端..."<<std::endl;
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
	printf("用例：CG法、PCG求特征值[OPTION]...\n");
	printf("Test the use of time\n");
	printf("Example:  CG \ PCG测试用例，维度在512-10240之间，以256递增做测试\n");
	printf("./FEMCal_CUDA --matrix-size=1024 --iters=1000 --precision=1e-5 --ResultFile=export_eig.dat increment=256 \n");
	printf("\n");
	printf("Options:\n");
	printf("--help\tDisplay this help menu\n");
	printf("--matrix-size\t初始矩阵维度，默认为512\n");
	printf("--iters\t每次计算求特征值矩阵的迭代次数上限，默认为1000\n");
	printf("--precision\t迭代默认误差精度，默认为1e-5\n");
	printf("--increment\t作性能测试时每次矩阵维度递增数值，默认为不递增\n");
	printf("--ResultFile\t输出特征值文件名，路径为..//data//export//下，默认为export_eig.dat\n");
}