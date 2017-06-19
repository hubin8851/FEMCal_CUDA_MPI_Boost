#pragma once
#include "stdafx.h"

extern int g_MatRowNum;

class CImpFrmMatlab
{
public:
	CImpFrmMatlab(void);
	~CImpFrmMatlab(void);

	ImpMatError_t ImpFrmMatlab( const char* _matname,  boost::filesystem::path _searchpath = "C:\\Users\\hbx\\Desktop\\SuiteSparse-2.1.1\\UFget");
//	ImpMatError_t ImpFrmMatlab( boost::filesystem::path _searchpath, const char* _matname = "UF_Index.mat" );

	//���CSR��ʽ��ϵ������ģ�庯����Ĭ�ϲ���double
	HBXDef::_CSRInput<double>& GetStiffMat(bool _bsave = false, boost::filesystem::path _SavePath="..\\data\\UFget\\" );
	//���������ʽ���Ҷ�����
	double*	GetRhsVec(bool _bsave = false);
	//���nnA�ĳ���
	int& GetNoneZeroLgt();
	//���rowsort�ĳ��ȣ���rownum+1
	int& GetnALgt();
	//��ȡ��ǰ�����ID�ţ���ΪĳЩ�����CSR��ʽ��rowsort������
	int& GetID();	
private:
	//���ڶ�ȡ��������������Ϣ
	void ReadObject( const mxArray *const _Array );
	//���ڶ�ȡUF_Index�ڵ���Ϣ
	void ReadMsg( const mxArray* const _Array );
	void ReadSparseMat( const mxArray *const _Array );
	void ReadRhs( const mxArray *const _Array );
	void ReadCell( const mxArray* const _Array );

	//���ⲿ��ȡ��������ָ��õ��նȾ����CSR��ʽ
	void	SetStiffMat(const double* _srcVal, const size_t* _srcCol, size_t* _srcRow, size_t _nnA, size_t _nA);
	//�����ⲿ����غ�����ָ��
	void	SetLoadVec(const double* _LoadVec, int _RowNum);

	//��������ĵ�ʽ�Ҷ�����
	void	genRhs( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "..\\data\\UFget" );

	boost::filesystem::path	m_SavePath;

	MATFile*	m_pmatFile;    
	mxArray*	m_pMxArray;

	std::map<std::string, mxArray*>	_mxArray_Map;
	typedef std::map<std::string, mxArray*>::iterator	mx_Map_itr;

	mwSize m_RowNum, m_ColNum;//�����������
	size_t m_id;	//�����ڱ�׼���е�ID��

	HBXDef::_CSRInput<double>  m_MatMsg;
	//CSR��ʽ�ĸ�����
//	double *_pNoneZeroVal;
//	size_t *_piColSort;
//	size_t *_piNonZeroRowSort;
	//���x
	double	*_px;
	//�Ҷ�����
//	double	*_prhs;
	//CSR��ʽ�ĸ����鳤��
	int		m_nnzA;	//CSR��ʽ��һ���͵ڶ�������ĳ���;
	int		m_nA;	//CSR��ʽ����������ĳ��ȣ����ھ���ĳһά�ȳ���+1;

	std::vector<double>	_vNoneZeroVal;
	std::vector<int>	_viColSort;
	std::vector<int>	_viNonZeroRowSort;
	std::vector<double>	_vlhsX;
	std::vector<double>	_vRhsY;
	
};

