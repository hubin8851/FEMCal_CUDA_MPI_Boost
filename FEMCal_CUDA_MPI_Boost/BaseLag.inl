
namespace HBXDef
{
	template< unsigned int T >
	__host__ CBaseLag<T>::CBaseLag( HBXDef::_AEROTABLE* _IptTable, size_t _blkId )
	{
		Assert( T>0 )
		IFNULL( _IptTable, "�����������ݱ����...")
		m_isInit = false;
		m_table_p = _IptTable;
		this->SetBlkId( _blkId );

		memset( m_lag_cordinate, 0, T * sizeof(unsigned int) );
	}

	template< unsigned int T >
	__host__ CBaseLag<T>::CBaseLag( const _SameBaseLag& _rhs )
	{
		Assert( T>0 )
		memcpy( m_numperdim, _rhs.m_numperdim, T * sizeof(unsigned int) );

		memset( m_lag_cordinate, 0, T * sizeof(unsigned int) );

		m_cordials = _rhs.m_cordials;
		m_data = _rhs.m_data;

		m_isInit = true;
	}

	template< unsigned int T >
	CBaseLag<T>::~CBaseLag(void)
	{
		if (m_table_p)
		{
			delete	m_table_p;
			m_table_p = nullptr;
		}

		if (m_cordials)
		{
			delete	[]m_cordials;
			m_cordials = nullptr;
		}
	}

	template< unsigned int T >
	UserStatusError_t CBaseLag<T>::SetBlkId( size_t _blkId )
	{
		if (_blkId >= m_table_p->_blocknum)
		{
			return UserStatusError_t::USER_STATUS_INVALID_VALUE;
		}
		m_slct_block_num = _blkId;
		m_dim = m_table_p->_blocks[_blkId]._dim;
		if ( m_dim <= 0 || T != m_dim )
		{
			return UserStatusError_t::USER_STATUS_INVALID_VALUE;
		}
		
		return UserStatusError_t::USER_STATUS_SUCCESS;
	}

	template< unsigned int T >
	bool CBaseLag<T>::Initial()
	{
		memcpy( m_numperdim, m_table_p->_blocks[m_slct_block_num]._numperdim, sizeof(unsigned int)*T );

		unsigned int _cordlgth = 0;	//��ֵ�ڵ�����ĳ���
		for (int i = 0; i < m_dim; i++)
		{
			_cordlgth += m_numperdim[i];
		}
		m_cordials = new HBXDef::UserDefFloat[_cordlgth];

		_cordlgth = 0;	//�ڴ�����άתһά�������ʱ����
		for (int i = 0; i < m_dim; i++)
		{
			for (int j = 0; j < m_numperdim[i]; j++)
			{
				m_cordials[_cordlgth] = m_table_p->_blocks[m_slct_block_num]._corddata[i][j];
				_cordlgth++;
			}
		}
		
		m_data = m_table_p->_blocks[m_slct_block_num]._data;

		m_isInit = true;
		return UserStatusError_t::USER_STATUS_SUCCESS;
	}


	template< unsigned int T >
	HBXDef::UserDefFloat	CBaseLag<T>::ReadTableData( HBXDef::UserDefFloat* _ArrayIn )
	{
		unsigned int _tmp_position[ TPOW2(T) ];	//�ݹ���ʱ����
		HBXDef::UserDefFloat	_tmp_reducdata[ TPOW2(T) ];	//�ݹ�������ʱ����
		HBXDef::UserDefFloat	_tmp_proper[ T ];	//�ݹ����ϵ��

		memset((void*)&_tmp_position, 0, TPOW2(T) * sizeof(unsigned int) );

		findcordi( _tmp_position, _tmp_proper, _ArrayIn );
		//�����ݹ��㷨
		unsigned int _tmpidx = TPOW2(T)/2;
		for (unsigned int idx = 0; idx < TPOW2(T); idx++)//����ʱ������
		{
			_tmp_reducdata[idx] = m_data[ _tmp_position[idx] ];
		}
		for (unsigned int idx = T-1; idx > 0 ; idx--)
		{
			for (int idx2 = 0; idx2 < _tmpidx; idx2++)
			{//ѭ���ڱ��ʽ�����֮��Ϊb��i�� = a��i+offset��+coef*(a(i+offset)-a(i)),�����ֵ
				_tmp_reducdata[idx2] = _tmp_reducdata[idx2 + _tmpidx] - _tmp_proper[idx] * (_tmp_reducdata[idx2 + _tmpidx] - _tmp_reducdata[idx2]);
			}
			_tmpidx >> 1;//ͨ�����Ʋ�������һά�ȵ����ݽ��в�ֵ
		}//�ڴ�û�н����һ������ѭ���У�����δ�Ż����������ٶ�
		return _tmp_reducdata[1] - _tmp_proper[0] * (_tmp_reducdata[1] - _tmp_reducdata[0]);
	}

	template< unsigned int T >
	void	CBaseLag<T>::findcordi( unsigned int _tmp_position[TPOW2(T)], HBXDef::UserDefFloat _tmp_proper[T], HBXDef::UserDefFloat _vecIn[T] )
	{
		unsigned int tempmul = 1;		//ÿ��ά�Ȳ�ֵ�����ϵ�����ĳ˻�
		unsigned int temppow = 1;		//���ڿ��Ƶ���������vecÿ��ά�ȵ�ƫ����
		unsigned int _tmpoffset = 0;	//��ά����ʼ������һά�����е�����
		unsigned int		_length = 0;
		HBXDef::UserDefFloat*	_cordials = nullptr;
		unsigned int		_cordial = 0;

		for (unsigned int idx1 = 0; idx1 < T; idx1++)
		{
			_length = m_numperdim[idx1];	//Ѱ�ҵ�ǰά����ÿ��ά�ȵĽڵ���Ŀ

			_cordials = &m_cordials[_tmpoffset];	
			_cordial = m_lag_cordinate[idx1];

			//�жϵ�ǰ��ֵ������һ����������ҷ�
			if (_vecIn[idx1] > _cordials[_cordial])
			{
				for (int idx2 = _cordial; idx2 < _length; idx2++)
				{
					if ((_length - 1) == idx2)	//�ж��Ƿ����ұ߽�
					{
						_tmp_proper[idx1] = 0;	//�ݹ�������0
						m_lag_cordinate[idx1] = _length - 1; 
					}
					else
					{
						if (_vecIn[idx1] < _cordials[idx2+1])
						{
							m_lag_cordinate[idx1] = idx2 + 1;
							_tmp_proper[idx1] = (_cordials[idx2 + 1] - _vecIn[idx1]) / (_cordials[idx2 + 1] - _cordials[idx2]);
							break;
						}
					}
				}
			}
			else
			{//��
				for (int idx2 = _cordial; idx2 >= 0; idx2--)
				{
					if (idx2 == 0)//���ò�ֵ�㿿����˱߽�
					{
						_tmp_proper[idx1] = 1;				//�ݹ�����
						m_lag_cordinate[idx1] = 1;
						break;
					}
					else
					{
						if ( _vecIn[idx1] > _cordials[idx2-1] )	//����_cordials[idx2-1]��_cordials[idx2]֮��
						{
							m_lag_cordinate[idx1] = idx2;
							_tmp_proper[idx1] = (_cordials[idx2] - _vecIn[idx1]) / (_cordials[idx2] - _cordials[idx2-1]);
							break;
						}
					}
				}
			}
			for (int idx2 = 0; idx2 < TPOW2(T); idx2++)
			{
				if (0 != (idx2 / temppow) % 2)
					_tmp_position[idx2] += m_lag_cordinate[idx1] * tempmul;
				else
					_tmp_position[idx2] += (m_lag_cordinate[idx1]-1) * tempmul;
			}
			tempmul *= _length;
			_tmpoffset += _length;
			temppow = temppow << 1;
		}

	}

}


