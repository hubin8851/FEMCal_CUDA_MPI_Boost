#include "stdafx.h"
#include "XMLManage.h"
#include "AeroTable.h"


CXMLManage::CXMLManage( const std::string& _SourceFile, boost::filesystem::path _savepath )
{
	m_SrcFileName = _SourceFile;
	m_path = _savepath;
}

CXMLManage::~CXMLManage(void)
{
}

bool CXMLManage::SetSourceFilePath( const std::string& _SourceFile, boost::filesystem::path _savepath )
{
	m_SrcFileName = _SourceFile;
	m_path = _savepath;
	return true;
}

// 获取所有的表号编码
// void	CXMLManage::ReadElemtProperty()
// {
// 	using namespace std;
// 	using namespace boost;
// 	using boost::property_tree::ptree;
// 
// 	std::string	_tmppath;
// 	_tmppath = m_path.string()+"\\";
// 	_tmppath.append(m_SrcFileName);
// 
// 	size_t	_iCount = 0;	//单元总和计数器
// 	ptree	pt;
// 	ptree	root;
// 	HBXFEMDef::ElemProperty	_tmpElPro;	//单元属性结构体临时变量
// 	string	_ElemtName;//单元名称，临时变量
// 	try
// 	{
// 		read_xml(_tmppath, pt);
// 		std::cout<<"打开节点属性配置文件成功"<<pt.data()<<std::endl;
// 		root = pt.get_child("root");
// 	}
// 	catch (std::exception& ec)
// 	{
// 		std::cout << "Error: " << ec.what() << std::endl;  
// 		return;  
// 	}
// 	std::cout << "root的长度即所有表的个数：" << root.size() << std::endl;
// 	for (BOOST_AUTO(pos, root.begin()); pos!=root.end(); ++pos)	//使用了boost中的auto
// 	{
// 		if (pos->first == "ElementObj")
// 		{
// 			ptree	_table = pos->second;	
// 			_ElemtName = _table.get<string>("<xmlattr>.name");
// 			_tmpElPro._type = _table.get<unsigned int>("<xmlattr>.Type");
// 			for (BOOST_AUTO(_itr, _table.begin()); _itr!=_table.end(); ++_itr)
// 			{
// 				if (_itr->first == "Property")
// 				{
// 					ptree	_data = _itr->second;
// 					if ( iequals("Stress", _data.get<string>("<xmlattr>.StressOrStrain")) )
// 					{
// 						_tmpElPro.StressOrStrain = 1;
// 					}
// 					else if ( iequals("strain", _data.get<string>("<xmlattr>.StressOrStrain")) )
// 					{
// 						_tmpElPro.StressOrStrain = 2;
// 					}
// 					else _tmpElPro.StressOrStrain = 0;
// 				}
// 			}
// 			m_ElemtPtyInMap.insert( make_pair(_ElemtName, _tmpElPro) );
// 			_iCount++;	//计数器只在《数据》子树下，无论配对成功与否均+1
// 		}
// 		
// 	}
// }

// HBXFEMDef::_EltPtyInMap&	CXMLManage::GetElemtPtyInMap()
// {
// 	return m_ElemtPtyInMap;
// }

//读取气动力参数表
void	CXMLManage::ReadAeroCoef( const char* _filename )
{
	using namespace std;
	using namespace boost;
	using boost::property_tree::ptree;

	TiXmlDocument _tmpdoc( _filename );
	TiXmlElement*	_tmpElmt = nullptr;
	HBXDef::_AEROTABLE*	_tmpTable = new HBXDef::_AEROTABLE();

	bool	_loadOkey = _tmpdoc.LoadFile();
	if (!_loadOkey)
	{
		std::cerr<<"载入XML文件错误，该文件路径为"<<_filename<<std::endl;
		std::cerr<<"当前错误码为"<<_tmpdoc.ErrorDesc()<<std::endl;
		exit(1);
	}
	_tmpElmt = _tmpdoc.RootElement();
	std::string _tmpstr(_tmpElmt->FirstAttribute()->Value());
	boost::to_lower(_tmpstr);
	m_AeroTableInMap.insert( make_pair(_tmpstr, _tmpTable) );

	if (!ReadTableAttri(_tmpElmt, _tmpTable))
	{
		std::cerr<<"载入表头属性出错..."<<std::endl;
		exit(1);
	}

	TiXmlElement* _blocks = _tmpElmt->FirstChildElement();
	for (size_t i=0; i<_tmpTable->_blocknum; i++)
	{
		if ( iequals( _blocks->Value(), "block" ) )
		{
			for (size_t j=0; j<_tmpTable->_blocknum; j++)
			{
				if (!ReadBlock(_blocks, _tmpTable->_blocks+j ))
				{
					std::cerr<<"读入BLOCK属性错误"<<std::endl;
					exit(1);
				}
			}
		}
	}
	
}

bool	CXMLManage::ReadTableAttri( TiXmlElement* _inputT, HBXDef::_AeroTable* _outTable )
{
	using namespace boost;
	if (nullptr == _inputT)
	{
		std::cerr<<"读入TinyXml的表指针错误..."<<std::endl;
		return false;
	}
	TiXmlAttribute* _TmpAttri = _inputT->FirstAttribute();
	while (_TmpAttri)//对block做循环
	{
		std::string _name = _TmpAttri->Name();
		if ( iequals(_name, "name") )
		{
			strcmp( _outTable->_name, _TmpAttri->Value() );
		}
		if ( iequals(_name, "block") )
		{
			_outTable->_blocknum = lexical_cast<int>( _TmpAttri->Value() );
			_outTable->_blocks = new HBXDef::_AeroBlock[_outTable->_blocknum];
		}

		_TmpAttri = _TmpAttri->Next();
	}
	delete _TmpAttri;
	return true;
}

bool	CXMLManage::ReadBlock( TiXmlElement* _inputT, HBXDef::_AeroBlock* _outBlock  )
{
	using namespace boost;
	if (nullptr == _inputT)
	{
		std::cerr<<"读入TinyXml的表指针错误..."<<std::endl;
		return false;
	}
	if ( !ReadBlockAttri( _inputT, _outBlock ) )
	{
		std::cerr<<"读入TinyXml的块属性错误..."<<std::endl;
		return false;
	}
	TiXmlElement* _tmpNode = _inputT->FirstChildElement();

	size_t _tmpdataNum = 1;
	if ( iequals( _tmpNode->Value(), "demtion") )
	{
		for (size_t i=0; i<_outBlock->_dim; i++)
		{
			_tmpdataNum *= ReadDemtion( _tmpNode, _outBlock, i );
			if ( _tmpdataNum < 0 )
			{
				std::cerr<<"读入TinyXml的维度属性错误..."<<std::endl;
				return false;
			}
			_tmpNode = _tmpNode->NextSiblingElement();
		}
	}
	//TiXmlNode*  _tmpNode = _tmpNode;
	if ( equals( _tmpNode->Value(), "Data") )
	{
		_outBlock->_data = new HBXDef::UserDefFloat[_tmpdataNum];	//在此时所有参数已知
		ReadInterpolationData( _tmpNode, _outBlock, _tmpdataNum );
	}
	return true;
}

bool	CXMLManage::ReadBlockAttri( TiXmlElement* _inputT, HBXDef::_AeroBlock* _outBlock )
{
	using namespace boost;
	if (nullptr == _inputT)
	{
		std::cerr<<"读入TinyXml的块属性指针错误..."<<std::endl;
		return false;
	}
	TiXmlAttribute* _TmpAttri = _inputT->FirstAttribute();
	while (_TmpAttri)//对block内的属性做循环
	{
		std::string _name = _TmpAttri->Name();
		if ( iequals(_name, "name") )
		{
			strcmp( _outBlock->_name, _TmpAttri->Value() );
		}
		if ( iequals(_name, "demention") )
		{
			_outBlock->_dim = lexical_cast<int>( _TmpAttri->Value() );
			if (_outBlock->_dim>0 && _outBlock->_dim<100)	//这个参数范围合适么？...
			{
				_outBlock->_numperdim = new unsigned int[_outBlock->_dim];
				_outBlock->_corddata = new HBXDef::UserDefFloat*[_outBlock->_dim];//对二维数组的外层维度分配空间
				//_outBlock->_beg = new unsigned int[_outBlock->_dim];
			}
			else
			{
				std::cout<<"当前块内维度不满足要求..."<<std::endl;
			}
		}

		_TmpAttri = _TmpAttri->Next();
	}
	delete _TmpAttri;
	return true;
}

//@_inputT：输入的xml树形结构的某个节点
//@_outBlock：输出的_AeroBlock指针
//@_idx：_AeroBlock下的维度表的索引号
//返回值：当前维度下的插值点个数
size_t	CXMLManage::ReadDemtion( TiXmlElement* _inputT, HBXDef::_AeroBlock* _outBlock, size_t _idx )
{
	using namespace boost;
	HBXDef::IFNULL(_inputT, "读入TinyXml的表指针错误...");

	TiXmlAttribute* _TmpAttri = _inputT->FirstAttribute();
	std::string _strText;	//最长4G，足够
	vector<string> vstrText;

	while ( _TmpAttri )
	{
		if ( iequals( _TmpAttri->Name(), "name" ) )
		{
			//此处留下
		}
		if ( iequals( _TmpAttri->Name(), "demention" ) )
		{
			_outBlock->_numperdim[_idx] = _TmpAttri->IntValue();
			//为当前维度插值点数组分配内存
#ifdef DEBUG
			_outBlock->_corddata[_idx] = new HBXDef::UserDefFloat[_TmpAttri->IntValue()];
#endif // DEBUG
			_outBlock->_corddata[_idx] = new HBXDef::UserDefFloat[_TmpAttri->IntValue()];
		}
		_TmpAttri = _TmpAttri->Next();
	}
	_strText = _inputT->GetText();

	// 定义分割方式为英文逗号，中文逗号和空格，构造一个分词器，  
	boost::char_separator<char> sep(" ,，*");  
	typedef boost::tokenizer<boost::char_separator<char> >  CustonTokenizer;  
	CustonTokenizer tok(_strText, sep);
	// 输出分割结果
	vstrText.clear();
	size_t _idy = 0; 
	for(CustonTokenizer::iterator beg=tok.begin(); beg!=tok.end();++beg, _idy++)  
	{  
		vstrText.push_back(*beg);	//验证用
		_outBlock->_corddata[_idx][_idy] = lexical_cast<double>(*beg);
	}
	if (vstrText.size() != _outBlock->_numperdim[_idx])//判断属性内的参数与实际数据个数是否相符
	{
		std::cerr<<"维度内属性与插值点个数不相符..."<<std::endl;
	}

	return _outBlock->_numperdim[_idx];
}

bool	CXMLManage::ReadInterpolationData( TiXmlElement* _inputT, HBXDef::_AeroBlock* _outBlock, size_t _lgt )
{
	using namespace boost;
	HBXDef::IFNULL( _inputT, "插值数据表导入错误..." )

	std::string _strText;
	vector<string> vstrText;
//	std::strstream	thestream;
// 	thestream << _inputT->GetText();
// 	thestream << '\t';	
// 	for(size_t i=0; i<_lgt; i++)  
// 	{  
// 		if (thestream.eof()) return 0;
// 		thestream >> _outBlock->_data[i];
// 	}
	_strText = _inputT->GetText();
	HBXDef::IFNULL( _strText.c_str(), "当前块内数据为空" )
	// 定义分割方式为英文逗号，中文逗号和空格，构造一个分词器，  
	boost::char_separator<char> sep(" ,，*");  
	typedef boost::tokenizer<boost::char_separator<char> >  CustonTokenizer;  
	CustonTokenizer tok(_strText, sep);
	// 输出分割结果
	vstrText.clear();
	size_t _idx = 0; 
	for(CustonTokenizer::iterator beg=tok.begin(); beg!=tok.end();++beg, _idx++)  
	{  
		vstrText.push_back(*beg);	//验证用
		_outBlock->_data[_idx] = lexical_cast<double>(*beg);
	}
	if ( vstrText.size() != _lgt )//判断属性内的参数与实际数据个数是否相符
	{
		std::cerr<<"维度内属性与插值点个数不相符..."<<std::endl;
	}

	return true;
}


HBXDef::_AEROTABLE&	CXMLManage::GetAeroTable( std::string _str )
{
	_AeroTableMapIter _iter;
	boost::to_lower(_str);
	_iter = m_AeroTableInMap.find( _str );
	if (_iter == m_AeroTableInMap.end())
	{
		cout<<_str<<"关联的气动数据表未找到"<<endl;
	}
	return *_iter->second;
}

