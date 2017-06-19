namespace HBXFEMDef
{
	template<class T>
	void _B31Elemt<T>::LocalStiffMatCal()
	{
		double l2, l3, EA, eiy, eiz, c1y, c1z;
		CheckUserDefErrors( InertialCal() );	//计算转动惯量
		CheckUserDefErrors( KappayCoefCal() );	//计算y方向修正系数
		CheckUserDefErrors( KappazCoefCal() );	//计算z方向修正系数
		l2 = m_Lgth * m_Lgth;	//l^2
		l3 = l2 * m_Lgth;		//l^3
		EA = m_Matreial_ptr->dE * m_Sec_ptr->dA;
		c1y = 1 + m_kappay;
		c1z = 1 + m_kappaz;
		//局部坐标系下矩阵赋值
		_Stiffmat( 0, 0 ) = _Stiffmat( 6, 6 ) = EA / m_Lgth;
		_Stiffmat( 0, 6 ) = -EA / m_Lgth;

		_Stiffmat( 1, 1 ) = 12 * eiz / (l3 * c1y);
		_Stiffmat( 1, 5 ) = 6 * eiz / (l2 * c1y);
		_Stiffmat( 1, 7 ) = -12 * eiz / (l3 * c1y);
		_Stiffmat( 1, 11 ) = 6 * eiz / (l2 * c1y);

		_Stiffmat( 2, 2 ) = 12 * eiy / (l3 * c1z);
		_Stiffmat( 2, 4 ) = -6 * eiy / (l2 * c1z);
		_Stiffmat( 2, 8 ) = -12 * eiy / (l3 * c1z);
		_Stiffmat( 2, 10 ) = -6 * eiy / (l2 * c1z);

		_Stiffmat( 3, 3 ) = _Stiffmat( 9, 9 ) = m_Matreial_ptr->dG * m_Ix / m_Lgth;
		_Stiffmat( 3, 9 ) = -_Stiffmat( 3, 3 );

		_Stiffmat( 4, 4 ) = (4 + m_kappaz) * eiy / ( m_Lgth * c1z );
		_Stiffmat( 4, 8 ) = 6 * eiy / (l2 * c1z);
		_Stiffmat( 4, 10) = ( 2 - m_kappaz) * eiy /(m_Lgth * c1z);

		_Stiffmat( 5, 5 ) = (4 + m_kappay) * eiz / (m_Lgth * c1y);
		_Stiffmat( 5, 7 ) = -6 * eiz / (l2 * c1y);
		_Stiffmat( 5, 11) = ( 2 - m_kappay) * eiz /(m_Lgth * c1y);

		_Stiffmat( 7, 7 ) = 12 * eiz / (l3 * c1y);
		_Stiffmat( 7, 11 ) = -6 * eiz / (l2 * c1y);

		_Stiffmat( 8, 8 ) = 12 * eiy / (l3 * c1z);
		_Stiffmat( 8, 10 ) = 6 * eiy / (l2 * c1z);

		_Stiffmat( 10, 10 ) = (4 + m_kappaz) * eiy / ( m_Lgth * c1z );
		_Stiffmat( 11, 11 ) = (4 + m_kappay) * eiz / ( m_Lgth * c1y );
		_Stiffmat.Symmetrized();

		return;
	}

	template<class T>
	void _B31Elemt<T>::StiffMatCal()
	{
		this->LocalStiffMatCal();

	};

	template<class T>
	void _B31Elemt<T>::MassMatCal( MassType_t _MassType )
	{

	}


}