#include "stdafx.h"
#include "HbxFEMFunc.h"
#include "HbxGloFunc.h"
#include "XMLManage.h"
#include "DataInput_f.h"

extern	int		g_Dim;	//ȫ��ά�ȣ���������2��3ά��Ͻ�ģ
extern	int		g_SectionNum,g_NodeNum,g_ElementNum,g_MaterialNum;
extern	int		g_SetNum;//�ڵ�͵�Ԫ��������
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

	g_MatRowNum = g_NodeNum * g_Dim;		//�ڵ�������ά�ȡ�Ŀǰ��˵û���⣬�����ʺ϶�����ά�������
	std::cout<<"�ܽ����Ϊ��"<<g_NodeNum<<" �ܵ�Ԫ��Ϊ��"<<g_ElementNum<<std::endl;
	std::cout<<"�ܲ�����Ϊ��"<<g_MaterialNum<<"�ܽ�����Ϊ��"<<g_SectionNum<<std::endl;
	std::cout<<"������Ԫ�ͽڵ�ļ�����Ϊ��"<<g_SetNum<<std::endl<<std::endl;
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

	float	rdens = 0.0;//�ض�
	std::string stringLine, _MatName, _SectionName;	//��ȡ�ĵ�ǰ�����ݣ��������ƣ��������� 
	std::string _Elset, _Nset;	//��Ԫ���ϣ��ڵ㼯��
	vector<string> vstrLine;
	unsigned int _tmptype;//��Ԫ���ͺ��ַ���
	int _markloop;
	int	 marknumber;//�������͵��ж�
	ControlMark_t ipControlMark;//���������
	std::string _SetName;
	std::string	_tmppath;
	_tmppath = m_path.string()+"\\";
	_tmppath.append(m_SrcFileName);
	std::ifstream inFile;
	inFile.open( _tmppath );
	if(inFile.bad())
	{
		std::cerr<<"Error1�� Can not open input data file"<<std::endl;
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
		// ����ָʽΪӢ�Ķ��ţ����Ķ��źͿո񣬹���һ���ִ�����  
		boost::char_separator<char> sep(" ,��*");  
		typedef boost::tokenizer<boost::char_separator<char> >  CustonTokenizer;  
		CustonTokenizer tok(stringLine,sep);
		// ����ָ���
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
		if ( iequals("Node", vstrLine[0]) )//�ж��Ƿ����ڵ��
		{
			_markloop = 10;
			ipControlMark = NODE;
		}
		else if ( iequals("Element", vstrLine[0]) && icontains(stringLine, "type=") )//�ж��Ƿ���뵥Ԫ��
		{
			ierase_first( vstrLine[1], "type=" );
			g_PlaneStressOrStrain = m_ElemtPropInMap[vstrLine[1]].StressOrStrain;
			_tmptype = lexical_cast<unsigned int>( to_string( m_ElemtPropInMap[vstrLine[1]]._type) );//��õ�Ԫ���
			_Dim = lexical_cast<int>( to_string( m_ElemtPropInMap[vstrLine[1]]._type ).substr(0,1) );
			_stressFlag = lexical_cast<unsigned short>( to_string( m_ElemtPropInMap[vstrLine[1]]._type ).substr(1,2) );
			_nodenum = lexical_cast<unsigned short>( to_string( m_ElemtPropInMap[vstrLine[1]]._type ).substr(3,3) );
			//����
			_markloop = 20;
			ipControlMark = ELEMENT;
		}
		//��ȡ�����������
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
		//��ȡ����������ơ����ü���
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
		//��ȡ����������ơ����ü���
		else if ( iequals("Material", vstrLine[0]) )
		{
			ierase_first( vstrLine[1], "name=" );
			_MatName = vstrLine[1];
			_markloop = 0;
			ipControlMark = MATERIAL;
		}
		else if( ipControlMark ==MATERIAL && iequals("Density", vstrLine[0]) ) _markloop=ControlMark_t::DENSITY;
		else if( ipControlMark ==MATERIAL && iequals("Elastic", vstrLine[0]) ) _markloop=ControlMark_t::ELASTIC;
		else if ( all( vstrLine[0].substr(0,1) , is_digit()) )//all�������һ���ַ����е�����Ԫ���Ƿ�����������жϹ���boost�ṩ��һЩ�ж�������Ե�boost��������
		{
			marknumber = _markloop;
		}
		switch ( marknumber )
		{
		case ControlMark_t::NODE:
			_tmpNode = boost::make_shared<_Node<float>>( lexical_cast<float>(vstrLine[1]), lexical_cast<float>(vstrLine[2]), lexical_cast<float>(vstrLine[3]) );//��ά��࣬�����ά�������Զ���0
			g_vNode_f.insert( make_pair( lexical_cast<int>(vstrLine[0])-1, _tmpNode) );
			break;
		case ControlMark_t::ELEMENT:
			switch(_tmptype)
			{
			case 200003:
				//_tmpElem = boost::make_shared<_Tri2D3NElemt<float>>( (unsigned char)_Dim, _stressFlag, vstrLine.begin()+1, vstrLine.size()-1 );//��Ϊ��֪��һ����Ԫ�м����ڵ㣬ֻ�ܴ��������
				g_vBaseElem_f.push_back(_tmpElem);
				break;
			case 300008:
				//_tmpElem = boost::make_shared<_C3D8RElemt<float>>( (unsigned char)_Dim, _stressFlag, vstrLine.begin()+1, vstrLine.size()-1 );//��Ϊ��֪��һ����Ԫ�м����ڵ㣬ֻ�ܴ��������
				g_vBaseElem_f.push_back(_tmpElem);
				break;
			default:
				break;
			}
			break;
		case ControlMark_t::DENSITY:
			rdens = lexical_cast<float>(vstrLine[0]) * 9.8f;	//�ض� = �ܶ�*9.8
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

//�����еĵ�Ԫ�Ľ���Ͳ��Ͻ������
void	CDataInput_f::SetElementSectionAndMaterial()
{
	if ( 0 == g_SetNum && 0 == g_SectionNum)
	{
		std::cout<<"�޼����޽��棬ѡ�����0��Ϊ���е�Ԫ����"<<std::endl;
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
