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

	//输出CSR形式的系数矩阵模板函数，默认参数double
	HBXDef::_CSRInput<double>& GetStiffMat(bool _bsave = false, boost::filesystem::path _SavePath="..\\data\\UFget\\" );
	//输出向量形式的右端向量
	double*	GetRhsVec(bool _bsave = false);
	//输出nnA的长度
	int& GetNoneZeroLgt();
	//输出rowsort的长度，即rownum+1
	int& GetnALgt();
	//获取当前矩阵的ID号，因为某些矩阵的CSR格式其rowsort有问题
	int& GetID();	
private:
	//用于读取单个矩阵的相关信息
	void ReadObject( const mxArray *const _Array );
	//用于读取UF_Index内的信息
	void ReadMsg( const mxArray* const _Array );
	void ReadSparseMat( const mxArray *const _Array );
	void ReadRhs( const mxArray *const _Array );
	void ReadCell( const mxArray* const _Array );

	//从外部获取三个数组指针得到刚度矩阵的CSR格式
	void	SetStiffMat(const double* _srcVal, const size_t* _srcCol, size_t* _srcRow, size_t _nnA, size_t _nA);
	//从类外部获得载荷数组指针
	void	SetLoadVec(const double* _LoadVec, int _RowNum);

	//生成随机的等式右端向量
	void	genRhs( size_t _genRowNum = g_MatRowNum, bool _bsave = false, boost::filesystem::path _savepath = "..\\data\\UFget" );

	boost::filesystem::path	m_SavePath;

	MATFile*	m_pmatFile;    
	mxArray*	m_pMxArray;

	std::map<std::string, mxArray*>	_mxArray_Map;
	typedef std::map<std::string, mxArray*>::iterator	mx_Map_itr;

	mwSize m_RowNum, m_ColNum;//矩阵的行列数
	size_t m_id;	//矩阵在标准库中的ID号

	HBXDef::_CSRInput<double>  m_MatMsg;
	//CSR格式的各数组
//	double *_pNoneZeroVal;
//	size_t *_piColSort;
//	size_t *_piNonZeroRowSort;
	//左端x
	double	*_px;
	//右端数组
//	double	*_prhs;
	//CSR格式的各数组长度
	int		m_nnzA;	//CSR格式第一个和第二个数组的长度;
	int		m_nA;	//CSR格式第三个数组的长度，等于矩阵某一维度长度+1;

	std::vector<double>	_vNoneZeroVal;
	std::vector<int>	_viColSort;
	std::vector<int>	_viNonZeroRowSort;
	std::vector<double>	_vlhsX;
	std::vector<double>	_vRhsY;
	
};

