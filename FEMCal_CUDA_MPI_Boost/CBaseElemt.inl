namespace HBXFEMDef
{
	//@_Num_:当前单元编号
	//@_type_:单元类型
	//@_lhs：节点编号向量迭代器
	//@_nodenum：节点总数
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
		//单元刚度矩阵与总刚矩阵映射关系计算
		int i,j,i1;
		_ipmark.resize( _nodenum*_dim_ );
		for(i=0,i1=0; i<_nodenum; i++)
		{
			for(j=0; j<_idim; j++,i1++) _ipmark[i1] = _dim_*(this->_vNode[i]) + j;
		}
 	}

}
