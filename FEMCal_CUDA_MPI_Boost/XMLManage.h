#pragma once
#include <tinyxml/tinyxml.h>
#include <tinyxml/tinystr.h>

using namespace std;


class CXMLManage
{
public:
	CXMLManage( const std::string& _SourceFile = "EltProperty.xml", boost::filesystem::path _savepath = "F:\\data from HBX_phd\\VS2010\\FEMCal_CUDA_MPI_Boost\\data\\DataIn" );	//构造函数
	int	TestXML(const char* _filename);
	bool	SetSourceFilePath( const std::string& _SourceFile = "EltProperty.xml", boost::filesystem::path _savepath = "F:\\data from HBX_phd\\VS2010\\FEMCal_CUDA_MPI_Boost\\data\\DataIn" );
	void	ReadElemtProperty();		//获取所有的表号编码
	HBXFEMDef::_EltPtyInMap&	GetElemtPtyInMap();		//获取[单元名称-(属性)]-对象名称的映射

	//读取气动力参数表相关函数
	void	ReadAeroCoef( const char* _filename );				
	bool	ReadTableAttri( TiXmlElement* _inputT, boost::shared_ptr<HBXDef::_AeroTable> _outTable );
	bool	ReadBlock( TiXmlElement* _inputT, HBXDef::_AeroBlock* _outBlock  );
	bool	ReadBlockAttri( TiXmlElement* _inputT, HBXDef::_AeroBlock* _outBlock );
	//@_inputT：输入的xml树形结构的某个节点
	//@_outBlock：输出的_AeroBlock指针
	//@_idx：_AeroBlock下的维度表的索引号
	//返回值：当前维度下的插值点个数
	size_t	ReadDemtion( TiXmlElement* _inputT, HBXDef::_AeroBlock* _outBlock, size_t _idx );
	bool	ReadInterpolationData( TiXmlElement* _inputT, HBXDef::_AeroBlock* _outBlock, size_t _lgt );

	HBXDef::_AEROTABLE&	GetAeroTable( std::string _str );	//获取气动数据表


	~CXMLManage(void);

private:
	boost::filesystem::path	m_path;
	std::string		m_SrcFileName;

	HBXFEMDef::_EltPtyInMap	m_ElemtPtyInMap;

	std::map < std::string, boost::shared_ptr<HBXDef::_AeroTable> >	m_AeroTableInMap;
};

