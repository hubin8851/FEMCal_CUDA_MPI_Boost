#ifndef _DATAINPUT_D_H
#define _DATAINPUT_D_H


class CDataInput_d
{
public:
	CDataInput_d(HBXFEMDef::_EltPtyInMap & _EltPty);
	~CDataInput_d(void);

	bool	SetSourceFilePath( const std::string &_SourceFile, boost::filesystem::path _savepath = "F:\\data from HBX_phd\\VS2010\\FEMCal_CUDA_MPI_Boost\\data\\DataIn" );
	bool	SetInputData();//完成inp文件的导入，填充各容器并计算总刚矩阵的维度等参数
private:
	bool	ReadInpFile();	//读取数据文件
	bool	ReadMtxFile();	//读取Abaqus导出的mtx文件
	void	SetElementSectionAndMaterial();//对所有的单元的截面和材料进行填充

private:

	std::string	m_SrcFileName;
	boost::filesystem::path	m_path;

	HBXFEMDef::_EltPtyInMap m_ElemtPropInMap;
};

#endif



