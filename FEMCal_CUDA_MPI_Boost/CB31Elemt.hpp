#ifndef _CB31ELEMT_HPP_
#define _CB31ELEMT_HPP_


#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//ȫ�ֱ������������ݶ�ȡ�ͼ���ľ���(float or double)

namespace HBXFEMDef
{

	//��ά2�ڵ�����Ԫ
	template <class T>
	class _B31Elemt:
		public _BaseElem<T>
	{
	protected:
		double	m_Lgth;	//����Ԫ����
		double  m_kappay, m_kappaz;	//���Ǽ��б��ε�Ӱ��ʱ����Ԫ�նȾ���y��z���������ϵ��
		double	m_Iy, m_Iz, m_Ix;	//���Ծغͼ����Ծ�
		boost::shared_ptr<_BaseSection<T>>	m_Sec_ptr;	//��������ָ��
		HBXDef::CMatrix<T>	_B;	//��Ԫ�նȾ����������е���ʱ����
		HBXDef::CMatrix<T>	_D;	//��Ԫ�նȾ����������е���ʱ����,���Ծ���	
		HBXDef::CMatrix<T>  _Loc;//�ڲ����ڸ�˹���ֵ��С����
	private:
		//���㵱ǰ��Ԫ�Ĺ��ԾصȲ���
		UserStatusError_t InertialCal()
		{
			m_Iy = (double)( m_Sec_ptr->dWide * pow(m_Sec_ptr->dH, 3) / 12. );
			m_Iz = (double)( m_Sec_ptr->dH * pow(m_Sec_ptr->dWide, 3) / 12. );
			m_Ix = (double)( 0.229 * m_Sec_ptr->dWide * pow(m_Sec_ptr->dH, 3) );
		}
		//���㿼�Ǽ��б��ε�Ӱ��ʱ����Ԫ�նȾ�������ϵ��
		UserStatusError_t KappayCoefCal()
		{
			// m_kappay = (12*E*Iy)/(k*G*A*l^2),�μ����ṹ�����е�����Ԫ���������ơ�P90
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
		//���캯��No.1
		//@_iNum_:��ǰ��Ԫ�ı��
		//@_type_:��Ԫ����
		//@_lhs���ڵ�������������
		//@_nodenum���ڵ�����
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

		//���ֲ�����ϵ�µĵ�Ԫ�նȾ���,���ա��ṹ�����е�...��P90��2017-4-28��hbx
		void LocalStiffMatCal();
		//���ȫ������ϵ�µ�Ke��Ԫ�նȾ���
		void StiffMatCal();
		//���ֲ�����ϵ�µĵ�Ԫ��������
		void LocalMassMatCal();
		//�����������Ĭ�ϼ���������
		void MassMatCal( MassType_t _MassType = MassType_t::LUMPED );

	}
	typedef boost::shared_ptr<_B31Elemt<float>> __B31Elemt_f_ptr;
	typedef boost::shared_ptr<_B31Elemt<double>> __B31Elemt_d_ptr;
}

#include "CB31Elemt.inl"

#endif // !_CB31ELEMT_HPP_