#ifndef _HBXFEMDEFSTRUCT_
#define _HBXFEMDEFSTRUCT_

#include "CBaseSection.hpp"
#include "CBaseMaterial.hpp"
#include "CBaseElemt.hpp"

#define TRUSS 1000
#define FRAME 1001

extern	HBXDef::PrecsType_t	g_PrecsType;	//ȫ�ֱ������������ݶ�ȡ�ͼ���ľ���(float or double)

//2015-1-28, HBX, ���嵥Ԫ���ͺ�//
//2015-9-30,HBX,��Ѫ�޸�...��...

namespace HBXFEMDef
{
	//��Ԫ���Խṹ�壬������Ԫ��ţ�ƽ�棨��ά��Ӧ����Ӧ�����⣬���ڳ����޸�
	struct ElemProperty
	{
		unsigned int		_type;			//��Ԫ���͵ı��
		char	StressOrStrain;	//Ӧ����Ӧ������
	};


	//�ڵ㶨��ṹ��
	template <class T=float>
	class _Node
	{
	public:
		_Node(T _x=0, T _y=0, T _z=0):_X(_x),_Y(_y),_Z(_z)
		{
		}
		short iType;	//�ڵ�����
		T	_X, _Y, _Z;	//�ڵ�����
		int iaDOFIndex[3];	//�ڵ����ɶȱ��
	};


	//�ڵ��Ԫ���Ͻṹ��
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
		unsigned int _begin;	//��ʼ��Ԫ��
		unsigned int _end;		//������Ԫ��
		//Ŀǰinp�ļ��ڶ�Ӧ��generate��δ��ֵ
		bool _ElOrN;	//�ڵ㻹�ǵ�Ԫ���ϱ�־λ��false��ʾ�ڵ㣬true��ʾ��Ԫ
	};

	template< class T=float >
	class _BaseLoad
	{
	public:
		char	iType;			//�غ�����
		unsigned int	iLoadedElem;	//�غ����õĵ�Ԫ��
	};
	//2015-1-30��HBX
	//2016-1-11,hbx,�޸�.ȡ��ģ����;2016-11-10,�ָ�ģ����
	//�߽���������Ԫ���غ�
	template<typename T>
	class _TractionLoad:
		public _BaseLoad<T>
	{
		int		iNodeNum[3]; // �ڵ���
		T		m_qx[3];
		T		m_qy[3];
		T		m_qz[3];
	};
	//2015-1-30��HBX
	//�߽����������нڵ���
	template<typename T=float>
	class _ConcentrationLoad:
		public _BaseLoad<T>
	{
		int iLoadNode;		//�غ����õĽڵ���
		T	m_q[3];			//�غ���������Ĵ�С
	};

	//2015-1-30��HBX
	//�߽���������Լ���ڵ�
	template <class T>
	class _NodeRestrain
	{
		int iRNode;	//��Լ���ڵ��
		int kRDerecIp[3];
		// kRDerecIp=1��x derection restrain; kRDerecIp=2��y derection restrain;
		// kRDerecIp=3��x and y derection restrain;
		T uxyRest[3]; // x,y�����λ��
	};

	/************************************************************************/
	/*                     ����Ԫ�������								                                                */
	/************************************************************************/
//Ϊ�˳������չ�ԣ��ɴ�XML���ȡÿ�ֵ�Ԫ��Ӧ������
	typedef std::map < std::string, ElemProperty >	_EltPtyInMap;
	typedef _EltPtyInMap::iterator					_EltPtyMapIter;

//	typedef std::vector < boost::shared_ptr<Node<float>> > _vNode;
	typedef	std::map< int, boost::shared_ptr< HBXFEMDef::_Node<float>> >	NodeInMap_f;//ӳ���ϵ�µ�Node����������
	typedef std::map< int, boost::shared_ptr< HBXFEMDef::_Node<float>> >::iterator	 _NodeInMapIter_f;	//������
	typedef	std::map< int, boost::shared_ptr< HBXFEMDef::_Node<double>> >	NodeInMap_d;//ӳ���ϵ�µ�Node����������
	typedef std::map< int, boost::shared_ptr< HBXFEMDef::_Node<double>> >::iterator	 _NodeInMapIter_d;	//������

	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseElem<float>> >	_vBaseElem_f;	
	typedef	_vBaseElem_f::iterator											_vBaseElemIter_f;//��Ԫ�����ĵ�����(float)
	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseElem<double>> >	_vBaseElem_d;	
	typedef	_vBaseElem_d::iterator											_vBaseElemIter_d;//��Ԫ�����ĵ�����(double)

	typedef std::vector < boost::shared_ptr<HBXFEMDef::_Material<float>> >	_vMaterial_f;//��������(float)
	typedef std::vector < boost::shared_ptr<HBXFEMDef::_Material<double>> >	_vMaterial_d;//��������(double)

	typedef std::map < std::string, boost::shared_ptr<HBXFEMDef::_Cset> > _vSet;//�ڵ㼯������
	typedef std::pair < std::string, boost::shared_ptr<HBXFEMDef::_Cset> > _Set_pair_;
	typedef std::map < std::string, boost::shared_ptr<HBXFEMDef::_Cset> >::iterator		_vSetIter;

	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseSection<float>> >	_vSection_f;//��������(float)
	typedef	_vSection_f::iterator												_vSectionIter_f;	//��������������(float)
	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseSection<double>> >	_vSection_d;//��������(double)
	typedef	_vSection_d::iterator												_vSectionIter_d;	//��������������(double)

	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseLoad<float>> >	_vBaseLoad_f;//�߽���������
	typedef std::vector < boost::shared_ptr<HBXFEMDef::_BaseLoad<double>> > _vBaseLoad_d;//�߽���������

//	typedef	std::map< int, boost::shared_ptr< HBXFEMDef::BaseElem<float>> >	ElemInMap;//ӳ���ϵ�µ�Element����������
}



//2016-1-8,HBX,��ɾ��...





#endif