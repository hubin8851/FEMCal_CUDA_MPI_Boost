#pragma once

#include "stdafx.h"

extern	int		g_Dim;	//ȫ��ά�ȣ���������2��3ά��Ͻ�ģ
extern	int		g_NodeNumPerElement;//��Ԫ�ڽڵ���
extern	int		g_SectionNum,g_NodeNum,g_ElementNum,g_MaterialNum, g_SetNum;
extern	int		g_LoadingTractionNum,g_LoadingConcentrationNum,g_NodeRestrainNum;	//�غ����ȫ�ֱ���
extern	int		g_RowPerElem;	//��Ԫ������� dim*g_NodeNumPerElement
extern	int		g_MatRowNum;//ȫ�־������
extern	int		g_EleNumInMatrix;//�����е�Ԫ��
extern	int		g_WeightMark;	//�Ƿ���������
extern	int		g_PlaneStressOrStrain;	//ƽ��Ӧ��orӦ������

extern HBXFEMDef::_vMaterial_f		g_vMaterial_f;	//ȫ�ֲ�������
extern HBXFEMDef::_vSection_f		g_vSection_f;	//ȫ�ֽ�������
extern HBXFEMDef::NodeInMap_f		g_vNode_f;		//ȫ�ֽڵ�����,�мǣ���map��Ϊ����������
extern HBXFEMDef::_vBaseElem_f		g_vBaseElem_f;	//ȫ�ֻ���Ԫ����
extern HBXFEMDef::_vSet				g_vNSet;//�ڵ㼯������
extern HBXFEMDef::_vSet				g_vElSet;//��Ԫ��������

extern HBXFEMDef::NodeInMap_d		g_vNode_d;		//ȫ�ֽڵ�����,�мǣ���map��Ϊ����������
extern HBXFEMDef::_vBaseElem_d		g_vBaseElem_d;	//ȫ�ֻ���Ԫ����
extern HBXFEMDef::_vMaterial_d		g_vMaterial_d;	//ȫ�ֲ�������
extern HBXFEMDef::_vSection_d		g_vSection_d;	//ȫ�ֽ�������


class MatrixAssemble
{
public:
	MatrixAssemble(void);
	~MatrixAssemble(void);

	void	AssembleAnalysis( HBXDef::PrecsType_t _precs = HBXDef::PrecsType_t::_PRECS_FLOAT_ );		//����������򣬶���ӿ�
	void	PreCondit( HBXDef::PrecsType_t _precs );//Ԥ����������ݣ����縳��Ԫ�Ĳ�������
	void	StiffMatrixCal( HBXDef::PrecsType_t _precs );	//�ܸվ������װ
	void	MassMatrixCal( HBXDef::PrecsType_t _precs );		//�������������װ
	void	LoadingVecCal();	//�������ʸ����װ
	void	EssentialConditionEnforce();	//�߽��������뵽�ܸվ�����������ʸ����
	//����ÿ����Ԫ�ĵ�Ԫ�նȾ���
	//int		GetElemStiff(bool _gpuflag = false);

	//���������ļ�
	void	ExportResultFile(  boost::filesystem::path _path = "..\\data\\source", HBXDef::SaveType_t _save_t = HBXDef::SaveType_t::CSR, HBXDef::PrecsType_t _precs= HBXDef::PrecsType_t::_PRECS_DOUBLE_ );		

private:
	HBXDef::PrecsType_t m_precs;

	boost::filesystem::path	m_SavePath;
	//CSR�洢��ʽ
	//CSR��һ�������ʾ����Ԫ��ֵ
	std::vector<double>  _vNoneZeroVal;
	//m_iColSort��ʾ����CSRѹ����ʽ�ĵڶ������飬��¼ԭ��������Ӧ���к�
	std::vector<int> _viColSort;
	//m_iNonZeroRowSort��ʾ����CSRѹ����ʽ�ĵ��������飬��ʾԭ�����һ������Ԫ�ڵ�һ�������е�����ֵ
	std::vector<int> _viNonZeroRowSort;


	void ElementStiffCal();		//��Ԫ�նȾ�����㺯��
	void ElementMassCal();	//��Ԫ����������㺯��
	void AssembleStiffMat( HBXDef::PrecsType_t _precs= HBXDef::PrecsType_t::_PRECS_FLOAT_ );	//��װ��Ԫ�նȾ����Ϊ��װ����

	HBXDef::CMatrix<float>	m_SumStiffMat_f;
	HBXDef::CMatrix<double>	m_SumStiffMat_d;

	HBXDef::CMatrix<float>	m_SumMassMat_f;
	HBXDef::CMatrix<double>	m_SumMassMat_d;
};

