#pragma once

#include "stdafx.h"

extern	int		g_Dim;	//全局维度，不适用于2、3维混合建模
extern	int		g_NodeNumPerElement;//单元内节点数
extern	int		g_SectionNum,g_NodeNum,g_ElementNum,g_MaterialNum, g_SetNum;
extern	int		g_LoadingTractionNum,g_LoadingConcentrationNum,g_NodeRestrainNum;	//载荷相关全局变量
extern	int		g_RowPerElem;	//单元矩阵阶数 dim*g_NodeNumPerElement
extern	int		g_MatRowNum;//全局矩阵阶数
extern	int		g_EleNumInMatrix;//矩阵中单元数
extern	int		g_WeightMark;	//是否引入重力
extern	int		g_PlaneStressOrStrain;	//平面应力or应变问题

extern HBXFEMDef::_vMaterial_f		g_vMaterial_f;	//全局材料容器
extern HBXFEMDef::_vSection_f		g_vSection_f;	//全局截面容器
extern HBXFEMDef::NodeInMap_f		g_vNode_f;		//全局节点容器,切记！是map！为了搜索更快
extern HBXFEMDef::_vBaseElem_f		g_vBaseElem_f;	//全局基单元容器
extern HBXFEMDef::_vSet				g_vNSet;//节点集合容器
extern HBXFEMDef::_vSet				g_vElSet;//单元集合容器

extern HBXFEMDef::NodeInMap_d		g_vNode_d;		//全局节点容器,切记！是map！为了搜索更快
extern HBXFEMDef::_vBaseElem_d		g_vBaseElem_d;	//全局基单元容器
extern HBXFEMDef::_vMaterial_d		g_vMaterial_d;	//全局材料容器
extern HBXFEMDef::_vSection_d		g_vSection_d;	//全局截面容器


class MatrixAssemble
{
public:
	MatrixAssemble(void);
	~MatrixAssemble(void);

	void	AssembleAnalysis( HBXDef::PrecsType_t _precs = HBXDef::PrecsType_t::_PRECS_FLOAT_ );		//总体分析程序，对外接口
	void	PreCondit( HBXDef::PrecsType_t _precs );//预处理相关数据，比如赋单元的材料属性
	void	StiffMatrixCal( HBXDef::PrecsType_t _precs );	//总刚矩阵的组装
	void	MassMatrixCal( HBXDef::PrecsType_t _precs );		//总质量矩阵的组装
	void	LoadingVecCal();	//总体荷载矢量组装
	void	EssentialConditionEnforce();	//边界条件加入到总刚矩阵和总体荷载矢量中
	//计算每个单元的单元刚度矩阵
	//int		GetElemStiff(bool _gpuflag = false);

	//结果输出到文件
	void	ExportResultFile(  boost::filesystem::path _path = "..\\data\\source", HBXDef::SaveType_t _save_t = HBXDef::SaveType_t::CSR, HBXDef::PrecsType_t _precs= HBXDef::PrecsType_t::_PRECS_DOUBLE_ );		

private:
	HBXDef::PrecsType_t m_precs;

	boost::filesystem::path	m_SavePath;
	//CSR存储格式
	//CSR第一个数组表示非零元的值
	std::vector<double>  _vNoneZeroVal;
	//m_iColSort表示的是CSR压缩格式的第二个数组，记录原矩阵非零对应的列号
	std::vector<int> _viColSort;
	//m_iNonZeroRowSort表示的是CSR压缩格式的第三个数组，表示原矩阵第一个非零元在第一个数组中的索引值
	std::vector<int> _viNonZeroRowSort;


	void ElementStiffCal();		//单元刚度矩阵计算函数
	void ElementMassCal();	//单元质量矩阵计算函数
	void AssembleStiffMat( HBXDef::PrecsType_t _precs= HBXDef::PrecsType_t::_PRECS_FLOAT_ );	//组装单元刚度矩阵成为总装矩阵

	HBXDef::CMatrix<float>	m_SumStiffMat_f;
	HBXDef::CMatrix<double>	m_SumStiffMat_d;

	HBXDef::CMatrix<float>	m_SumMassMat_f;
	HBXDef::CMatrix<double>	m_SumMassMat_d;
};

