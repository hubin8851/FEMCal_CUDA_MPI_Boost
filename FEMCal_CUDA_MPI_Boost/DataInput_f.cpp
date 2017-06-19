#include "stdafx.h"
#include "HbxFEMFunc.h"
#include "HbxGloFunc.h"
#include "XMLManage.h"
#include "DataInput_f.h"

extern	int		g_Dim;	//全局维度，不适用于2、3维混合建模
extern	int		g_SectionNum,g_NodeNum,g_ElementNum,g_MaterialNum;
extern	int		g_SetNum;//节点和单元集合总数
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


CDataInput_f::CDataInput_f(HBXFEMDef::_EltPtyInMap & _EltPty):m_ElemtPropInMap(_EltPty)
{
}


CDataInput_f::~CDataInput_f(void)
{
}

bool	CDataInput_f::SetSourceFilePath( const std::string &_SourceFile, boost::filesystem::path _savepath )
{
	m_SrcFileName = _SourceFile;
	m_path = _savepath;
	return true;
}


bool	CDataInput_f::SetInputData()
{
	ReadInpFile();
	g_NodeNum	= GetTotalNodeNum_f();
	g_ElementNum	= GetTotalElementNum_f();
	g_MaterialNum	= GetTotalMaterialNum_f();
	g_SectionNum	= GetTotalSectionNum_f();
	g_SetNum = GetTotalSetNum(false) + GetTotalSetNum(true);

	g_MatRowNum = g_NodeNum * g_Dim;		//节点数乘以维度。目前来说没问题，但不适合二、三维混合网格
	std::cout<<"总结点数为："<<g_NodeNum<<" 总单元数为："<<g_ElementNum<<std::endl;
	std::cout<<"总材料数为："<<g_MaterialNum<<"总截面数为："<<g_SectionNum<<std::endl;
	std::cout<<"包括单元和节点的集合数为："<<g_SetNum<<std::endl<<std::endl;
	SetElementSectionAndMaterial();
	return true;
}

bool	CDataInput_f::ReadInpFile()
{
	using namespace std;
	using namespace boost;
	using namespace HBXDef;
	using namespace HBXFEMDef;
	int		_Dim;
	unsigned short _stressFlag;
	unsigned short _nodenum ;

	float	rdens = 0.0;//重度
	std::string stringLine, _MatName, _SectionName;	//读取的当前行内容，材料名称，截面名称 
	std::string _Elset, _Nset;	//单元集合，节点集合
	vector<string> vstrLine;
	unsigned int _tmptype;//单元类型号字符串
	int _markloop;
	int	 marknumber;//数据类型的判定
	ControlMark_t ipControlMark;//控制流标记
	std::string _SetName;
	std::string	_tmppath;
	_tmppath = m_path.string()+"\\";
	_tmppath.append(m_SrcFileName);
	std::ifstream inFile;
	inFile.open( _tmppath );
	if(inFile.bad())
	{
		std::cerr<<"Error1： Can not open input data file"<<std::endl;
		exit(1);
	}
	ipControlMark = RESET;
	boost::shared_ptr<_Node<float>>	_tmpNode;
	boost::shared_ptr<_BaseElem<float>> _tmpElem;
	boost::shared_ptr< _Material<float> > _tmpMaterial;
	boost::shared_ptr< _BaseSection<float> > _tmpSection;
	boost::shared_ptr< _Cset > _tmpSet;
	while(!inFile.eof())
	{
		getline(inFile, stringLine);
		// 定义分割方式为英文逗号，中文逗号和空格，构造一个分词器，  
		boost::char_separator<char> sep(" ,，*");  
		typedef boost::tokenizer<boost::char_separator<char> >  CustonTokenizer;  
		CustonTokenizer tok(stringLine,sep);
		// 输出分割结果
		vstrLine.clear();
		for(CustonTokenizer::iterator beg=tok.begin(); beg!=tok.end();++beg)  
		{  
			vstrLine.push_back(*beg);  
		}  
		marknumber = 0;
		if ( "" == vstrLine[0]   )
		{
			_markloop = 0;
			continue;
		}
		if ( iequals("Node", vstrLine[0]) )//判断是否进入节点段
		{
			_markloop = 10;
			ipControlMark = NODE;
		}
		else if ( iequals("Element", vstrLine[0]) && icontains(stringLine, "type=") )//判断是否进入单元段
		{
			ierase_first( vstrLine[1], "type=" );
			g_PlaneStressOrStrain = m_ElemtPropInMap[vstrLine[1]].StressOrStrain;
			_tmptype = lexical_cast<unsigned int>( to_string( m_ElemtPropInMap[vstrLine[1]]._type) );//获得单元编号
			_Dim = lexical_cast<int>( to_string( m_ElemtPropInMap[vstrLine[1]]._type ).substr(0,1) );
			_stressFlag = lexical_cast<unsigned short>( to_string( m_ElemtPropInMap[vstrLine[1]]._type ).substr(1,2) );
			_nodenum = lexical_cast<unsigned short>( to_string( m_ElemtPropInMap[vstrLine[1]]._type ).substr(3,3) );
			//待续
			_markloop = 20;
			ipControlMark = ELEMENT;
		}
		//读取集合相关名称
		else if ( iequals("Elset", vstrLine[0]) )
		{
			ierase_first( vstrLine[1], "elset=" );
			_SetName = vstrLine[1];
			_markloop = 51;
			ipControlMark = ELSET;
		}
		else if ( iequals("Nset", vstrLine[0]) )
		{
			ierase_first( vstrLine[1], "nset=" );
			_SetName = vstrLine[1];
			_markloop = 52;
			ipControlMark = NSET;
		}
		//读取截面相关名称、所用集合
		else if ( iequals("section:", vstrLine[0]) )
		{
			_SectionName = vstrLine[1];
			_markloop = 0;
			ipControlMark = SECTION;
		}
		else if ( ipControlMark ==SECTION && iequals("Section", vstrLine[1]) )
		{
			ierase_first( vstrLine[2], "elset=" );
			_Elset = vstrLine[2];
			ierase_first( vstrLine[3], "material=" );
			_MatName = vstrLine[3];
			_markloop = 42;
		}
		//读取材料相关名称、所用集合
		else if ( iequals("Material", vstrLine[0]) )
		{
			ierase_first( vstrLine[1], "name=" );
			_MatName = vstrLine[1];
			_markloop = 0;
			ipControlMark = MATERIAL;
		}
		else if( ipControlMark ==MATERIAL && iequals("Density", vstrLine[0]) ) _markloop=ControlMark_t::DENSITY;
		else if( ipControlMark ==MATERIAL && iequals("Elastic", vstrLine[0]) ) _markloop=ControlMark_t::ELASTIC;
		else if ( all( vstrLine[0].substr(0,1) , is_digit()) )//all：（检测一个字符串中的所有元素是否满足给定的判断规则）boost提供的一些判定规则可以到boost官网进行
		{
			marknumber = _markloop;
		}
		switch ( marknumber )
		{
		case ControlMark_t::NODE:
			_tmpNode = boost::make_shared<_Node<float>>( lexical_cast<float>(vstrLine[1]), lexical_cast<float>(vstrLine[2]), lexical_cast<float>(vstrLine[3]) );//三维最多，如果两维，最后的自动置0
			g_vNode_f.insert( make_pair( lexical_cast<int>(vstrLine[0])-1, _tmpNode) );
			break;
		case ControlMark_t::ELEMENT:
			switch(_tmptype)
			{
			case 200003:
				//_tmpElem = boost::make_shared<_Tri2D3NElemt<float>>( (unsigned char)_Dim, _stressFlag, vstrLine.begin()+1, vstrLine.size()-1 );//因为不知道一个单元有几个节点，只能传入迭代器
				g_vBaseElem_f.push_back(_tmpElem);
				break;
			case 300008:
				//_tmpElem = boost::make_shared<_C3D8RElemt<float>>( (unsigned char)_Dim, _stressFlag, vstrLine.begin()+1, vstrLine.size()-1 );//因为不知道一个单元有几个节点，只能传入迭代器
				g_vBaseElem_f.push_back(_tmpElem);
				break;
			default:
				break;
			}
			break;
		case ControlMark_t::DENSITY:
			rdens = lexical_cast<float>(vstrLine[0]) * 9.8f;	//重度 = 密度*9.8
			g_WeightMark = 1;
			break;
		case ControlMark_t::ELASTIC:
			_tmpMaterial = boost::make_shared< _Material<float>>( lexical_cast<float>(vstrLine[0]), lexical_cast<float>(vstrLine[1]), rdens, _MatName.c_str() );
			g_vMaterial_f.push_back( _tmpMaterial );
			break;
		case 42:
			_tmpSection = boost::make_shared< _BaseSection<float>>( lexical_cast<float>(vstrLine[0]), _SectionName.c_str(), _Elset.c_str(), _MatName.c_str() );
			g_vSection_f.push_back( _tmpSection );
			break;
		case ControlMark_t::ELSET:
			_tmpSet = boost::make_shared< _Cset>( vstrLine.begin(), vstrLine.size(), true );
			g_vElSet.insert( make_pair( _SetName, _tmpSet ) );
			break;
		case ControlMark_t::NSET:
			_tmpSet = boost::make_shared< _Cset>( vstrLine.begin(), vstrLine.size(), false );
			g_vNSet.insert( make_pair( _SetName, _tmpSet ) );
			break;
		default:
			break;
		}
	}

	return true;
}

//对所有的单元的截面和材料进行填充
void	CDataInput_f::SetElementSectionAndMaterial()
{
	if ( 0 == g_SetNum && 0 == g_SectionNum)
	{
		std::cout<<"无集合无截面，选择材料0作为所有单元材料"<<std::endl;
		HBXFEMDef::_vBaseElemIter_f _ElemIter_f;
		for (auto _iter = g_vBaseElem_f.begin(); _iter!=g_vBaseElem_f.end(); ++_iter)
		{
			(*_iter)->SetMaterial( make_pair(g_vMaterial_f[0], nullptr) );
		}
		return;
	}
	for (auto _iter = g_vSection_f.begin(); _iter!=g_vSection_f.end(); ++_iter)
	{
		unsigned int _bg = g_vElSet[ (*_iter)->elsetName ]->_begin;
		unsigned int _ed = g_vElSet[ (*_iter)->elsetName ]->_end;
		for (unsigned int i=_bg; i<=_ed; i++)
		{
			g_vBaseElem_f[i]->SetMaterial( make_pair(g_vMaterial_f[0], *_iter) );
		}
	}

	return;
}
