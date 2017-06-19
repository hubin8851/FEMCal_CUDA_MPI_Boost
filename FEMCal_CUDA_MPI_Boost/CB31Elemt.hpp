#ifndef _CB31ELEMT_HPP_
#define _CB31ELEMT_HPP_


#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//全局变量，设置数据读取和计算的精度(float or double)

namespace HBXFEMDef
{

	//三维2节点梁单元
	template <class T>
	class _B31Elemt:
		public _BaseElem<T>
	{
	protected:
		double	m_Lgth;	//梁单元长度
		double  m_kappay, m_kappaz;	//考虑剪切变形的影响时梁单元刚度矩阵y、z方向的修正系数
		double	m_Iy, m_Iz, m_Ix;	//惯性矩和极惯性矩
		boost::shared_ptr<_BaseSection<T>>	m_Sec_ptr;	//截面智能指针
		HBXDef::CMatrix<T>	_B;	//单元刚度矩阵计算过程中的临时变量
		HBXDef::CMatrix<T>	_D;	//单元刚度矩阵计算过程中的临时变量,弹性矩阵	
		HBXDef::CMatrix<T>  _Loc;//内部用于高斯积分点的小矩阵
	private:
		//计算当前单元的惯性矩等参数
		UserStatusError_t InertialCal()
		{
			m_Iy = (double)( m_Sec_ptr->dWide * pow(m_Sec_ptr->dH, 3) / 12. );
			m_Iz = (double)( m_Sec_ptr->dH * pow(m_Sec_ptr->dWide, 3) / 12. );
			m_Ix = (double)( 0.229 * m_Sec_ptr->dWide * pow(m_Sec_ptr->dH, 3) );
		}
		//计算考虑剪切变形的影响时梁单元刚度矩阵修正系数
		UserStatusError_t KappayCoefCal()
		{
			// m_kappay = (12*E*Iy)/(k*G*A*l^2),参见《结构分析中的有限元法与程序设计》P90
			if ( m_kappay < 0.0 ) 
			{
				m_kappay = 12. * m_Matreial_ptr->dE * / (m_Lgth*m_Lgth);
			}

			return UserStatusError_t::USER_STATUS_SUCCESS;
		}
		UserStatusError_t KappazCoefCal()
		{
			return UserStatusError_t::USER_STATUS_SUCCESS;
		}
		UserStatusError_t LengthCal()
		{
			double _tmpx,_tmpy, _tmpz;
			if (m_Lgth<=0)
			{
				HBXFEMDef::_Node _nA,_nB;
				_nA = g_vNode_d[_vNode[0]];
				_nB = g_vNode_d[_vNode[1]];
				_tmpx = _nA._X - _nB._X;
				_tmpy = _nA._Y - _nB._Y;
				_tmpz = _nA._Z - _nB._Z;
				m_Lgth = sqrt( _tmpx*_tmpx + _tmpy*_tmpy + _tmpz*_tmpz );
			}
			return UserStatusError_t::USER_STATUS_SUCCESS;
		}
	public:
		//构造函数No.1
		//@_iNum_:当前单元的编号
		//@_type_:单元类型
		//@_lhs：节点编号向量迭代器
		//@_nodenum：节点总数
		_B31Elemt( unsigned int _Num_, short _type_, std::vector<std::string>::iterator _lhs, size_t _nodenum, unsigned char _idim_):
			_BaseElem( _Num_, _type_, _lhs, _nodenum, _idim_ )
		{
			m_Lgth = -1.0f;
			m_kappay = -1.0f;
			m_kappaz = -1.0f;
		}

		
		void SetMaterial( std::pair< boost::shared_ptr<_Material<T>>, boost::shared_ptr<_BaseSection<T>> > _MatAndSec )
		{
			m_Matreial_ptr = _MatAndSec.first;
			m_Section_ptr = _MatAndSec.second;
		};

		//求解局部坐标系下的单元刚度矩阵,参照《结构分析中的...》P90，2017-4-28，hbx
		void LocalStiffMatCal();
		//求解全局坐标系下的Ke单元刚度矩阵
		void StiffMatCal();
		//求解局部坐标系下的单元质量矩阵
		void LocalMassMatCal();
		//求解质量矩阵，默认集中质量法
		void MassMatCal( MassType_t _MassType = MassType_t::LUMPED );

	}
	typedef boost::shared_ptr<_B31Elemt<float>> __B31Elemt_f_ptr;
	typedef boost::shared_ptr<_B31Elemt<double>> __B31Elemt_d_ptr;
}

#include "CB31Elemt.inl"

#endif // !_CB31ELEMT_HPP_