#ifndef _MATRIXCAL_INL
#define _MATRIXCAL_INL

namespace HBXDef
{
	//����˷�����, �������mOut��
	template <class _TOut, class _Tlhs, class _Trhs>
	CMatrix<_TOut>& MatrixMultiply(CMatrix<_TOut>& matOut, const CMatrix<_Tlhs>& matlhs, const CMatrix<_Trhs>& matrhs)
	{
		//�ж���߾�����������ұ߾�����������
		Assert(matlhs.GetColNum() == matrhs.GetRowNum());
		//���ɾ����¶�����lhs������Ϊ�������������rhs��������Ϊ���������
		CMatrix<_TOut> _matTemp(matlhs.GetRowNum(), matrhs.GetColNum());

		for(size_t i = 0; i < _matTemp.GetRowNum(); i ++)
		{
			for(size_t j = 0; j < _matTemp.GetColNum(); j ++)
			{
				_matTemp(i, j) = _TOut(0);		//����ֵ0
				for(size_t k = 0; k < matlhs.GetColNum(); k ++)
				{
					_matTemp(i, j) += matlhs(i, k) * matrhs(k, j);
				}
			}
		}
		matOut = _matTemp;		//�������ת����mOut������
		return matOut;		//���ؽ������mOut
	}


	//ȫѡ��Ԫ�����������ʽ
	template<class _Ty>
	long double MatrixDetrminant( const CMatrix<_Ty>& _rhs )
	{
		long double MaxValue, tmp;
		size_t stRank = _rhs.GetColNum();
		//����Ϊ�����򷵻�0
		if ( stRank != _rhs.GetRowNum() )
		return long double(0);
		
		CMatrix< long double > m( stRank, stRank );

		for(size_t i=0; i<stRank; i++)
			for(size_t j=0; j<stRank; j++)
				m(i, j) = (long double)_rhs(i, j);	//��ʼ��

		size_t iSign, jSign;				// ��Ԫ��λ�ñ�־

		long double Det(1);					// ����ʽ��ֵ
		int	nSgn = 1;						// ����

		for ( size_t k=0; k<stRank-1; k++ )
		{
			MaxValue = 0.0;
			for(size_t i = k; i < stRank; i ++)
			{
				for(size_t j = k; j < stRank; j ++)//�Խ�������
				{
					tmp = fabsl(m(i, j));		//��m(i,j)����ֵ
					if(tmp > MaxValue)
					{
						MaxValue = tmp;
						iSign = i;			//������Ԫλ��
						jSign = j;
					}
				}
			}
			//�ж���Ԫ�Ƿ�Ϊ0
			if ( LongDoubleEqual(MaxValue, 0) )
				return long double(0);
			
			if ( iSign != k )//����ֵ���Ԫ�ز��ڵ�ǰ��
			{
				nSgn = -nSgn;
				for (size_t j = k; j < stRank; j++)
					swap( m(k,j), m(iSign,j) );
			}

			if(jSign != k)			//����ֵ���Ԫ�ز��ڵ�ǰ��
			{
				nSgn = -nSgn;		//�任����ʽ����
				for(size_t i = k; i < stRank; i ++)		//������
					swap(m(i, jSign), m(i, k));
			}

			Det *= m(k, k);					//�Խ�Ԫ���
			for(size_t i = k + 1; i < stRank; i ++)
			{
				long double d(m(i, k) / m(k, k));		//��Ԫ����
				for(size_t j = k + 1; j < stRank; j ++)	//����Ԫ�·�Ԫ����Ϊ0
					m(i, j) -= d * m(k, j);	//��ǰ��Ԫ����������Ԫ�����任
			}
		}

		return Det * nSgn * m(stRank - 1, stRank - 1);	//��������ʽ��ֵ
	}


	//ȫѡ��Ԫ��˹-Լ��(Gauss-Jordan)���������
	template <class _Ty>
	int MatrixInversionGS( CMatrix<_Ty>& _rhs )
	{
		size_t stRank = _rhs.GetColNum();	// �������
		if(stRank != _rhs.GetRowNum())
			return int(-1);					//rhs���Ƿ���

		std::valarray<size_t> is(stRank);		//�н�����Ϣ
		std::valarray<size_t> js(stRank);		//�н�����Ϣ

		CMatrix<_Ty> m(_rhs);					//����һmatrix����

		for(size_t k = 0; k < stRank; k++)
		{	// ȫѡ��Ԫ
			long double MaxValue(0);
			for(size_t i = k; i < stRank; i ++)
			{
				for(size_t j = k; j < stRank; j ++)
				{
					long double tmp(fabsl(m(i, j)));	//��m(i,j)����ֵ
					if(tmp > MaxValue)				//��Ԫ���ڶԽ�����
					{
						MaxValue = tmp;
						is[k] = i;					//��¼��Ԫ����
						js[k] = j;					//��¼��Ԫ����
					}
				}
			}

			if( LongDoubleEqual(MaxValue, 0) ) 
				return int(0);						//��ԪΪ0����������

			if(is[k] != k)							//��Ԫ���ڵ�ǰ��
			{
				for(size_t j = 0; j < stRank; j ++)	//������Ԫ��
					swap(m(k, j), m(is[k], j));
			}

			if(js[k] != k)							//��Ԫ���ڵ�ǰ��
			{
				for(size_t i = 0; i < stRank; i ++)	//������Ԫ��
					swap(m(i, k), m(i, js[k]));
			}

			m(k, k) = 1.0 / m(k, k);				//��Ԫ����
			for(size_t j = 0; j < stRank; j ++)
				if(j != k)
					m(k, j) *= m(k, k);

			for(size_t i = 0; i < stRank; i ++)
				if(i != k)
					for(size_t j = 0; j < stRank; j ++)	
						if(j != k)
							m(i, j) = m(i, j) - m(i, k) * m(k, j);

			for(size_t i = 0; i < stRank; i ++)
				if(i != k)
					m(i, k) = -m(i, k) * m(k, k);
		}

		for(int r = stRank - 1; r >= 0; r--)
		{
			if(js[r] != r)
				for(size_t j = 0; j < stRank; j ++)
					swap(m(r, j), m(js[r], j));
			if(is[r] != r)
				for(size_t i = 0; i < stRank; i ++)
					swap(m(i, r), m(i, is[r]));
		}

		_rhs = m;

		return int(1);
	}

}//namespace hbxdef

#endif //����_MATRIXCAL_INL