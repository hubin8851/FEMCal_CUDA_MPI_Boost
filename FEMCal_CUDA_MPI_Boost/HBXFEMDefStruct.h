#ifndef _HBXFEMDEFSTRUCT_
#define _HBXFEMDEFSTRUCT_

#include "CBaseSection.hpp"
#include "CBaseMaterial.hpp"
#include "CBaseElemt.hpp"

#define TRUSS 1000
#define FRAME 1001

extern	HBXDef::PrecsType_t	g_PrecsType;	//全局变量，设置数据读取和计算的精度(float or double)

//2015-1-28, HBX, 定义单元类型宏//
//2015-9-30,HBX,吐血修改...囧...

namespace HBXFEMDef
{
	//单元属性结构体，包括单元编号，平面（三维）应力、应变问题，便于程序修改
	struct ElemProperty
	{
		unsigned int		_type;			//单元类型的编号
		char	StressOrStrain;	//应力、应变问题
	};


	//节点定义结构体
	template <class T=float>
	class _Node
	{
	public:
		_Node(T _x=0, T _y=0, T _z=0):_X(_x),_Y(_y),_Z(_z)
		{
		}
		short iType;	//节点类型
		T	_X, _Y, _Z;	//节点坐标
		int iaDOFIndex[3];	//节点自由度编号
	};


	//节点或单元集合结构体
	class _Cset
	{
	public:
		_Cset( std::vector<std::string>::iterator _lhs, size_t _num, bool _whichset=1 ):_ElOrN(_whichset)
		{
			_vSort.clear();
			std::vector<std::string>::iterator _iter=_lhs;
			for ( unsigned int i=0; i<(unsigned int)_num;i++, _iter++)
			{
				_vSort.push_back(boost::lexical_cast<unsigned int>(*_iter)-1);
			}
			_begin = _vSort[0];
			if (1 == _vSort.size())
			{
				_end = _begin;
			}
			else
			{
				_end = _vSort[1];
			}
		}

	public:
		std::vector<unsigned int>	_vSort;
		unsigned int _begin;	//起始单元号
		unsigned int _end;		//结束单元号
		//目前inp文件内对应的generate尚未赋值
		bool _ElOrN;	//节点还是单元集合标志位，false表示节点，true表示单元
	};

	template< class T=float >
	class _BaseLoad
	{
	public:
		char	iType;			//载荷类型
		unsigned int	iLoadedElem;	//载荷作用的单元号
	};
	//2015-1-30，HBX
	//2016-1-11,hbx,修改.取消模板类;2016-11-10,恢复模板类
	//边界条件，单元面载荷
	template<typename T>
	class _TractionLoad:
		public _BaseLoad<T>
	{
		int		iNodeNum[3]; // 节点编号
		T		m_qx[3];
		T		m_qy[3];
		T		m_qz[3];
	};
	//2015-1-30，HBX
	//边界条件，集中节点力
	template<typename T=float>
	class _ConcentrationLoad:
		public _BaseLoad<T>
	{
		int iLoadNode;		//载荷作用的节点编号
		T	m_q[3];			//载荷三个方向的大小
	};

	//2015-1-30，HBX
	//边界条件，受约束节点
	template <class T>
	class _NodeRestrain
	{
		int iRNode;	//受约束节点号
		int kRDerecIp[3];
		// kRDerecIp=1：x derection restrain; kRDerecIp=2：y derection restrain;
		// kRDerecIp=3：x and y derection restrain;
		T uxyRest[3]; // x,y方向的位移
	};

	/************************************************************************/
	/*                     有限元相关容器								                                                */
	/************************************************************************/
//为了程序的拓展性，可从XML里读取每种单元对应的属性
	typedef std::map < std::string, ElemProperty >	_EltPtyInMap;
	typedef _EltPtyInMap::iterator					_EltPtyMapIter;

//	typedef std::vector < boost::shared_ptr<Node<float>> > _vNode;
	typedef	std::map< int, boost::shared_ptr< HBXFEMDef::_Node<float>> >	NodeInMap_f;//映射关系下的Node，索引更快
	typedef std::map< int, boost::shared_ptr< HBXFEMDef::_Node<float>> >::iterator	 _NodeInMapIter_f;	//迭代器
	typedef	std::map< int, boost::shared_ptr< HBXFEMDef::_Node<double>> >	NodeInMap_d;//映射关系下的Node，索引更快
	typedef std::map< int, boost::shared_ptr< HBXFEMDef::_Node<double>> >::iterator	 _NodeInMapIter_d;	//迭代器

	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseElem<float>> >	_vBaseElem_f;	
	typedef	_vBaseElem_f::iterator											_vBaseElemIter_f;//单元容器的迭代器(float)
	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseElem<double>> >	_vBaseElem_d;	
	typedef	_vBaseElem_d::iterator											_vBaseElemIter_d;//单元容器的迭代器(double)

	typedef std::vector < boost::shared_ptr<HBXFEMDef::_Material<float>> >	_vMaterial_f;//材料容器(float)
	typedef std::vector < boost::shared_ptr<HBXFEMDef::_Material<double>> >	_vMaterial_d;//材料容器(double)

	typedef std::map < std::string, boost::shared_ptr<HBXFEMDef::_Cset> > _vSet;//节点集合容器
	typedef std::pair < std::string, boost::shared_ptr<HBXFEMDef::_Cset> > _Set_pair_;
	typedef std::map < std::string, boost::shared_ptr<HBXFEMDef::_Cset> >::iterator		_vSetIter;

	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseSection<float>> >	_vSection_f;//截面容器(float)
	typedef	_vSection_f::iterator												_vSectionIter_f;	//截面容器迭代器(float)
	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseSection<double>> >	_vSection_d;//截面容器(double)
	typedef	_vSection_d::iterator												_vSectionIter_d;	//截面容器迭代器(double)

	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseLoad<float>> >	_vBaseLoad_f;//边界条件容器
	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseLoad<double>> > _vBaseLoad_d;//边界条件容器

//	typedef	std::map< int, boost::shared_ptr< HBXFEMDef::BaseElem<float>> >	ElemInMap;//映射关系下的Element，索引更快
}



//2016-1-8,HBX,已删除...





#endif