#ifndef _DATAINPUT_D_H
#define _DATAINPUT_D_H


class CDataInput_d
{
public:
	CDataInput_d(HBXFEMDef::_EltPtyInMap & _EltPty);
	~CDataInput_d(void);

	bool	SetSourceFilePath( const std::string &_SourceFile, boost::filesystem::path _savepath = "F:\\data from HBX_phd\\VS2010\\FEMCal_CUDA_MPI_Boost\\data\\DataIn" );
	bool	SetInputData();//���inp�ļ��ĵ��룬���������������ܸվ����ά�ȵȲ���
private:
	bool	ReadInpFile();	//��ȡ�����ļ�
	bool	ReadMtxFile();	//��ȡAbaqus������mtx�ļ�
	void	SetElementSectionAndMaterial();//�����еĵ�Ԫ�Ľ���Ͳ��Ͻ������

private:

	std::string	m_SrcFileName;
	boost::filesystem::path	m_path;

	HBXFEMDef::_EltPtyInMap m_ElemtPropInMap;
};

#endif



