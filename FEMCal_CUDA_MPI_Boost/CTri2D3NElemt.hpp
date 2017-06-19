#pragma once
#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//全局变量，设置数据读取和计算的精度(float or double)

namespace HBXFEMDef
{
	//平面3节点三角单元结构体
	template <class T>
	class _Tri2D3NElemt:
		public _BaseElem<T>
	{
	protected:
		HBXDef::CMatrix<T>	_B;	//单元刚度矩阵计算过程中的临时变量
		HBXDef::CMatrix<T>	_D;	//单元刚度矩阵计算过程中的临时变量
		boost::shared_ptr<_BaseSection<T>>	m_Section_ptr;
		short iSection;	//单元截面索引号
		int iaNode[3];	
		T Area;//面积

	public:
		_Tri2D3NElemt( unsigned int _Num_, short _type_, std::vector<std::string>::iterator _lhs, size_t _nodenum, unsigned char _dim_):
			_BaseElem( _Num_, _type_, _lhs, _nodenum, _dim_ )
		{
			_Stiffmat.ResetSize( _vNode.size()*_dim_, _vNode.size()*_dim_);
		};
		void StiffMatCal()
		{
			T xi, yi, xj, yj, xm, ym;
			if ( HBXDef::PrecsType_t::_PRECS_FLOAT_ == g_PrecsType )
			{
				xi = g_vNode_f[_vNode[0]]->_X;
				yi = g_vNode_f[_vNode[0]]->_Y;
				xj = g_vNode_f[_vNode[1]]->_X;
				yj = g_vNode_f[_vNode[1]]->_Y;
				xm = g_vNode_f[_vNode[2]]->_X;
				ym = g_vNode_f[_vNode[2]]->_Y;
			}
			else if ( HBXDef::PrecsType_t::_PRECS_DOUBLE_ == g_PrecsType )
			{
				xi = g_vNode_d[_vNode[0]]->_X;
				yi = g_vNode_d[_vNode[0]]->_Y;
				xj = g_vNode_d[_vNode[1]]->_X;
				yj = g_vNode_d[_vNode[1]]->_Y;
				xm = g_vNode_d[_vNode[2]]->_X;
				ym = g_vNode_d[_vNode[2]]->_Y;
			}
			else std::cout<<"单元刚度矩阵计算节点位置读取出错"<<std::endl;

			T thick = m_Section_ptr->sThick;
			T E = m_Matreial_ptr->dE;
			T Mu = m_Matreial_ptr->dMu;
			_D.ResetSize(3, 3);
			_D(0, 1) = _D(1, 0) = m_Matreial_ptr->dMu;
			if ( HBXDef::SimpleType_t::PLANESTRESS == _iType )
			{
				_D(0, 0) = _D(1, 1) = 1;
				_D(2, 2) = (1-Mu)/2;
				_D(0, 1) = _D(1, 0) = Mu;
				_D *= (E/(1-Mu*Mu));
			}
			else if( HBXDef::SimpleType_t::PLANESTRAIN == _iType )
			{
				_D(0, 0) = _D(1, 1) = 1-Mu;
				_D(2, 2) = (1-2*Mu)/2;
				_D*=E/(1+Mu)/(1-2*Mu);
			}
			//正式运算
			Area = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2;
			_B.ResetSize(3,6);
			_B(0, 0) = _B(2, 1) = yj-ym;
			_B(0, 2) = _B(2, 3) = ym-yi;
			_B(0, 4) = _B(2, 5) = yi-yj;
			_B(1, 1) = _B(2, 0) = xm-xj;
			_B(1, 3) = _B(2, 2) = xi-xm;
			_B(1, 5) = _B(2, 4) = xj-xi;
			_B/=(2*Area);
			_Stiffmat = thick*Area*(!_B)*_D*_B;
		};

		T VolumeCalc()
		{
			return (Area*m_Section_ptr->sThick);
		}

		void MassMatCal( MassType_t _MassType = MassType_t::LUMPED )
		{
			if ( m_Matreial_ptr->dMaGama < 1e-2)
			{
				std::cerr<<"无重度标志无法进行质量矩阵的计算..."<<std::endl;
				return;
			}
			_Massmat.ResetSize( _vNode.size()*_idim, _vNode.size()*_idim);//质量矩阵分配内存
			T __mass = VolumeCalc()*m_Matreial_ptr->dMaGama/9.8;
			if ( MassType_t::LUMPED == _MassType )
			{
				_Massmat(0,0) = _Massmat(1,1) = _Massmat(2,2) = _Massmat(3,3) = _Massmat(4,4) = _Massmat(5,5) = __mass/3;
				return;
			}	
			else if ( MassType_t::CONSISTENT == _MassType )
			{
				_Massmat(0,0) = _Massmat(1,1) = _Massmat(2,2) = _Massmat(3,3) = _Massmat(4,4) = _Massmat(5,5) = __mass/6;
				_Massmat(0,2) = _Massmat(0,4) = _Massmat(1,3) = _Massmat(1,5) = _Massmat(2,0) = _Massmat(2,4) = __mass/12;
				_Massmat(3,1) = _Massmat(3,5) = _Massmat(4,0) = _Massmat(4,2) = _Massmat(5,1) = _Massmat(5,3) = __mass/12;
				return;
			}
			else
			{
				std::cerr<<"输入质量矩阵算法不正确..."<<std::endl;
				return;
			}
		}

		void SetMaterial( std::pair< boost::shared_ptr<_Material<T>>, boost::shared_ptr<_BaseSection<T>> > _MatAndSec )
		{
			m_Matreial_ptr = _MatAndSec.first;
			m_Section_ptr = _MatAndSec.second;
		};

	};

	typedef boost::shared_ptr<_Tri2D3NElemt<float>> _Tri2D3NElemt_f_ptr;
	typedef boost::shared_ptr<_Tri2D3NElemt<double>> _Tri2D3NElemt_d_ptr;


}

