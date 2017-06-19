#pragma once
#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//ȫ�ֱ������������ݶ�ȡ�ͼ���ľ���(float or double)

namespace HBXFEMDef
{
	//��ά8�ڵ������嵥Ԫ
	template <class T>
	class _C3D8RElemt:
		public _BaseElem<T>
	{
	protected:
		HBXDef::CMatrix<T>	_B;	//��Ԫ�նȾ����������е���ʱ����
		HBXDef::CMatrix<T>	_D;	//��Ԫ�նȾ����������е���ʱ����,���Ծ���	
		HBXDef::CMatrix<T>  _Loc;//�ڲ����ڸ�˹���ֵ��С����

	private:
		//���ݲ��ɱȺ�����ģ����_D��ֵ
		void SetD()
		{
			if ( nullptr == m_Matreial_ptr )
			{
				std::cerr<<"��Ԫ����δ��ֵ"<<std::endl;
			}
			T _Mu = m_Matreial_ptr->dMu;
			T _E  = m_Matreial_ptr->dE;
			_D( 0, 0 ) = 1-_Mu;
			_D( 0, 1 ) = _Mu;
			_D( 0, 2 ) = _Mu;
			_D( 1, 0 ) = _Mu;
			_D( 1, 1 ) = 1-_Mu;
			_D( 1, 2 ) = _Mu;
			_D( 2, 0 ) = _Mu;
			_D( 2, 1 ) = _Mu;
			_D( 2, 2 ) = 1-_Mu;
			_D( 3, 3 ) = 0.5-_Mu;
			_D( 4, 4 ) = 0.5-_Mu;
			_D( 5, 5 ) = 0.5-_Mu;
			_D *= ( _E/((1+_Mu)*(1-2*_Mu)) ); 
		}

		//����˸��ڵ��xyz������
		void SetLoc()
		{
			if ( HBXDef::PrecsType_t::_PRECS_FLOAT_ == g_PrecsType )
			{
				_Loc( 0, 0 ) = g_vNode_f[_vNode[0]]->_X;
				_Loc( 0, 1 ) = g_vNode_f[_vNode[0]]->_Y;
				_Loc( 0, 2 ) = g_vNode_f[_vNode[0]]->_Z;
				_Loc( 1, 0 ) = g_vNode_f[_vNode[1]]->_X;
				_Loc( 1, 1 ) = g_vNode_f[_vNode[1]]->_Y;
				_Loc( 1, 2 ) = g_vNode_f[_vNode[1]]->_Z;
				_Loc( 2, 0 ) = g_vNode_f[_vNode[2]]->_X;
				_Loc( 2, 1 ) = g_vNode_f[_vNode[2]]->_Y;
				_Loc( 2, 2 ) = g_vNode_f[_vNode[2]]->_Z;
				_Loc( 3, 0 ) = g_vNode_f[_vNode[3]]->_X;
				_Loc( 3, 1 ) = g_vNode_f[_vNode[3]]->_Y;
				_Loc( 3, 2 ) = g_vNode_f[_vNode[3]]->_Z;
				_Loc( 4, 0 ) = g_vNode_f[_vNode[4]]->_X;
				_Loc( 4, 1 ) = g_vNode_f[_vNode[4]]->_Y;
				_Loc( 4, 2 ) = g_vNode_f[_vNode[4]]->_Z;
				_Loc( 5, 0 ) = g_vNode_f[_vNode[5]]->_X;
				_Loc( 5, 1 ) = g_vNode_f[_vNode[5]]->_Y;
				_Loc( 5, 2 ) = g_vNode_f[_vNode[5]]->_Z;
				_Loc( 6, 0 ) = g_vNode_f[_vNode[6]]->_X;
				_Loc( 6, 1 ) = g_vNode_f[_vNode[6]]->_Y;
				_Loc( 6, 2 ) = g_vNode_f[_vNode[6]]->_Z;
				_Loc( 7, 0 ) = g_vNode_f[_vNode[7]]->_X;
				_Loc( 7, 1 ) = g_vNode_f[_vNode[7]]->_Y;
				_Loc( 7, 2 ) = g_vNode_f[_vNode[7]]->_Z;
			}
			else if ( HBXDef::PrecsType_t::_PRECS_DOUBLE_ == g_PrecsType )
			{
				_Loc( 0, 0 ) = g_vNode_d[_vNode[0]]->_X;
				_Loc( 0, 1 ) = g_vNode_d[_vNode[0]]->_Y;
				_Loc( 0, 2 ) = g_vNode_d[_vNode[0]]->_Z;
				_Loc( 1, 0 ) = g_vNode_d[_vNode[1]]->_X;
				_Loc( 1, 1 ) = g_vNode_d[_vNode[1]]->_Y;
				_Loc( 1, 2 ) = g_vNode_d[_vNode[1]]->_Z;
				_Loc( 2, 0 ) = g_vNode_d[_vNode[2]]->_X;
				_Loc( 2, 1 ) = g_vNode_d[_vNode[2]]->_Y;
				_Loc( 2, 2 ) = g_vNode_d[_vNode[2]]->_Z;
				_Loc( 3, 0 ) = g_vNode_d[_vNode[3]]->_X;
				_Loc( 3, 1 ) = g_vNode_d[_vNode[3]]->_Y;
				_Loc( 3, 2 ) = g_vNode_d[_vNode[3]]->_Z;
				_Loc( 4, 0 ) = g_vNode_d[_vNode[4]]->_X;
				_Loc( 4, 1 ) = g_vNode_d[_vNode[4]]->_Y;
				_Loc( 4, 2 ) = g_vNode_d[_vNode[4]]->_Z;
				_Loc( 5, 0 ) = g_vNode_d[_vNode[5]]->_X;
				_Loc( 5, 1 ) = g_vNode_d[_vNode[5]]->_Y;
				_Loc( 5, 2 ) = g_vNode_d[_vNode[5]]->_Z;
				_Loc( 6, 0 ) = g_vNode_d[_vNode[6]]->_X;
				_Loc( 6, 1 ) = g_vNode_d[_vNode[6]]->_Y;
				_Loc( 6, 2 ) = g_vNode_d[_vNode[6]]->_Z;
				_Loc( 7, 0 ) = g_vNode_d[_vNode[7]]->_X;
				_Loc( 7, 1 ) = g_vNode_d[_vNode[7]]->_Y;
				_Loc( 7, 2 ) = g_vNode_d[_vNode[7]]->_Z;
			}
		}

	public:
		_C3D8RElemt( unsigned int _Num_, short _type_, std::vector<std::string>::iterator _lhs, size_t _nodenum, unsigned char _dim_ ):
			_BaseElem( _Num_, _type_, _lhs, _nodenum, _dim_ )
		{
			//���캯��
			_Stiffmat.ResetSize( _nodenum*_dim_, _nodenum*_dim_);
			_Massmat.ResetSize( _nodenum*_dim_, _nodenum*_dim_);
			_Loc.ResetSize( 8, 3 );
			_B.ResetSize( 6, 24 );
			_D.ResetSize( 6, 6 );
		}
		
		//����ÿ�����ֵ��ϵľ���˻�
		HBXDef::CMatrix<T> BDcalc( T _gridx, T _blkx, T _thrdx )
		{
			HBXDef::CMatrix<T>  dNsnt( 8, 3 );
			
			dNsnt( 0, 0 ) = - (1-_blkx) * (1+_thrdx)/8;
			dNsnt( 0, 1 ) = -(1-_gridx) * (1+_thrdx)/8;
			dNsnt( 0, 2 ) = (1-_gridx) * (1-_blkx)/8;

			dNsnt( 1, 0 ) = (1-_blkx) * (1+_thrdx)/8;
			dNsnt( 1, 1 ) = -(1+_gridx) * (1+_thrdx)/8;
			dNsnt( 1, 2 ) = (1+_gridx) * (1-_blkx)/8;
			
			dNsnt( 2, 0 ) = (1-_blkx) * (1-_thrdx)/8;
			dNsnt( 2, 1 ) = -(1+_gridx) * (1-_thrdx)/8;
			dNsnt( 2, 2 ) = -(1+_gridx) * (1-_blkx)/8;
			
			dNsnt( 3, 0 ) = -(1-_blkx) * (1-_thrdx)/8;
			dNsnt( 3, 1 ) = -(1-_gridx) * (1-_thrdx)/8;
			dNsnt( 3, 2 ) = -(1-_gridx) * (1-_blkx)/8;

			dNsnt( 4, 0 ) = -(1+_blkx) * (1+_thrdx)/8;
			dNsnt( 4, 1 ) = (1-_gridx) * (1+_thrdx)/8;
			dNsnt( 4, 2 ) = (1-_gridx) * (1+_blkx)/8;

			dNsnt( 5, 0 ) = (1+_blkx) * (1+_thrdx)/8;
			dNsnt( 5, 1 ) = (1+_gridx) * (1+_thrdx)/8;
			dNsnt( 5, 2 ) = (1+_gridx) * (1+_blkx)/8;

			dNsnt( 6, 0 ) = (1+_blkx) * (1-_thrdx)/8;
			dNsnt( 6, 1 ) = (1+_gridx) * (1-_thrdx)/8;
			dNsnt( 6, 2 ) = -(1+_gridx) * (1+_blkx)/8;

			dNsnt( 7, 0 ) = -(1+_blkx) * (1-_thrdx)/8;
			dNsnt( 7, 1 ) = (1-_gridx) * (1-_thrdx)/8;
			dNsnt( 7, 2 ) = -(1-_gridx) * (1+_blkx)/8;

			HBXDef::CMatrix<T> _TransdNsnt(!dNsnt);
			HBXDef::CMatrix<T> _TmpJ = _TransdNsnt * _Loc;	//3*3
			long double	_detJ = MatrixDetrminant( _TmpJ );

			HBXDef::MatrixInversionGS(_TmpJ);
			HBXDef::CMatrix<T> _dNxyz = _TmpJ * _TransdNsnt;//3*8

			for (int i = 0; i < 8; i++)
			{
				_B( 0, 3*i )   = _dNxyz( 0, i );
				_B( 1, 3*i+1 ) = _dNxyz( 1, i );
				_B( 2, 3*i+2 ) = _dNxyz( 2, i );
				_B( 3, 3*i )   = _dNxyz( 1, i );
				_B( 3, 3*i+1 ) = _dNxyz( 0, i );
				_B( 4, 3*i+1 ) = _dNxyz( 2, i );
				_B( 4, 3*i+2 ) = _dNxyz( 1, i );
				_B( 5, 3*i )   = _dNxyz( 2, i );
				_B( 5, 3*i+2 ) = _dNxyz( 0, i );
			}
			SetD();
			HBXDef::CMatrix<T> _TransB(!_B);
			HBXDef::CMatrix<T> _MatOut = _detJ * _TransB * _D * _B;
			return _MatOut;
		}

		//��˹�������Ke��Ԫ�նȾ���
		void StiffMatCal()
		{
			double gridx, gridw, blkx, blkw, thrdx, thrdw;
			SetLoc();
			double gsx[3] = { -0.7745966692, 0, 0.7745966692 };
			double gsw[3] = { 0.55555555556, 0.888888888889, 0.55555555556 };
			for ( int _indexi=0; _indexi<3; _indexi++)
			{
				gridx = gsx[_indexi];
				gridw = gsw[_indexi];
				for ( int _indexj=0; _indexj<3; _indexj++)
				{
					blkx = gsx[_indexj];
					blkw = gsw[_indexj];
					for ( int _indexk=0; _indexk<3; _indexk++)
					{
						thrdx = gsx[_indexk];
						thrdw = gsw[_indexk];
						_Stiffmat += BDcalc( gridx, blkx, thrdx ) * gridw * blkw * thrdw; 
					}
				}
			}
		}
		T VolumeCalc(){return 0;};
		void MassMatCal( MassType_t _MassType ){};
		void SetMaterial( std::pair< boost::shared_ptr<_Material<T>>, boost::shared_ptr<_BaseSection<T>> > _MatAndSec )
		{
			m_Matreial_ptr = _MatAndSec.first;
		};

	};
	typedef boost::shared_ptr<_C3D8RElemt<float>> _C3D8RElemt_f_ptr;
	typedef boost::shared_ptr<_C3D8RElemt<double>> _C3D8RElemt_d_ptr;
}