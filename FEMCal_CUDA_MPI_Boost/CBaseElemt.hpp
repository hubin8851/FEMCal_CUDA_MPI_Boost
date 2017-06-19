#pragma once
#include "stdafx.h"
#include "HBXFEMDefStruct.h"


extern	HBXDef::PrecsType_t	g_PrecsType;	//全局变量，设置数据读取和计算的精度(float or double)
extern	int		g_Dim;	//全局维度，不适用于2、3维混合建模

namespace HBXFEMDef
{
#pragma region 单元代码段

	template<class T=double>	//单元基类
	class  _BaseElem
	{
	protected:
		//short _iMaterial;	//单元材料索引号
	public:
		//@_Num_:当前单元编号
		//@_type_:单元类型
		//@_lhs：节点编号向量迭代器
		//@_nodenum：节点总数
		_BaseElem( unsigned int _Num_, unsigned short _type_, std::vector<std::string>::iterator _lhs, size_t _nodenum, unsigned char _dim_ );

		virtual void SetMaterial( std::pair< boost::shared_ptr<_Material<T>>, boost::shared_ptr<_BaseSection<T>> > _MatAndSec ) = 0;
		virtual void SetStiffMat( HBXDef::CMatrix<T> &_rhsmat )
		{
			_Stiffmat = _rhsmat;
		}
		virtual void StiffMatCal() = 0;	//纯虚函数，完成单元刚度矩阵的计算
		virtual T VolumeCalc() = 0;	//纯虚函数，完成单元体积的计算
		virtual void MassMatCal( MassType_t _MassType = MassType_t::LUMPED ) = 0;	//纯虚函数，完成单元质量矩阵的计算

	public:
		unsigned int	_iNum;	//当前单元编号
		unsigned char	_idim;	//每个节点自由度
		unsigned short	_iType;		//应力亦或应变问题
		std::vector<unsigned int>	_vNode;//单元内各节点号
		HBXDef::CMatrix<T>	_Stiffmat;	//将质量矩阵和刚度矩阵分别放在派生类和基类中，是因为若没有计算重度的要求，质量矩阵可不用，但刚度矩阵必须有，故为节省空间，分开
		HBXDef::CMatrix<T>	_Massmat;	//质量矩阵,还是置于基类中吧...
		std::vector<int>	_ipmark;	//单元刚度矩阵与总刚矩阵之间的映射表,即单元刚度阵中的i行(列)在总刚阵中的位置
		boost::shared_ptr<_Material<T>>	m_Matreial_ptr;
	};

#pragma  endregion 结束单元代码段

}

#include "CBaseElemt.inl"