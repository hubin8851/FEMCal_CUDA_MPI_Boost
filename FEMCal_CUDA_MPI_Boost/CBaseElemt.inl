namespace HBXFEMDef
{
	//@_Num_:��ǰ��Ԫ���
	//@_type_:��Ԫ����
	//@_lhs���ڵ�������������
	//@_nodenum���ڵ�����
	template<class T>
	_BaseElem<T>::_BaseElem( unsigned int _Num_, unsigned short _type_, std::vector<std::string>::iterator _lhs, size_t _nodenum, unsigned char _dim_ ):
		_iNum(_Num_),_iType(_type_),_idim(_dim_)
	{
		_vNode.clear();
		std::vector<std::string>::iterator _iter=_lhs;
		for ( size_t i=0; i<_nodenum; i++, _iter++)
		{
			_vNode.push_back(boost::lexical_cast<unsigned int>(*_iter)-1);
		}
		//��Ԫ�նȾ������ܸվ���ӳ���ϵ����
		int i,j,i1;
		_ipmark.resize( _nodenum*_dim_ );
		for(i=0,i1=0; i<_nodenum; i++)
		{
			for(j=0; j<_idim; j++,i1++) _ipmark[i1] = _dim_*(this->_vNode[i]) + j;
		}
 	}

}
