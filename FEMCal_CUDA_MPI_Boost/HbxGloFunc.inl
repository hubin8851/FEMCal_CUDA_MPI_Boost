

//�����ռ��ڲ��ֺ���
namespace HBXDef
{

	//���ļ��а��ж�ȡ���ݲ���������ָ��
	template<class T>
	size_t ReadVecByRowFromFile(std::vector<T> &_OutputPtr, const char *_FileName ,const boost::filesystem::path& _dir )
	{
		_OutputPtr.clear();
		char str[255];
		boost::filesystem::path _tmppath(_dir);
		_tmppath.append(_FileName);
		std::ifstream inFile(_tmppath.string());
		if (!inFile)
		{
			std::cout<<"��ȡ�ļ�����������!\n";
			return false;
		}
		int _offset = 0;
		while(!inFile.eof())//��������ļ�ĩβ����һ������ֵ����û����ʱ����0����0�ͼ���~
		{
			inFile.getline(str,255);
			if (0 != strlen(str))
			{
				_OutputPtr.push_back((T)atof(str));
			}		
			_offset++;
		}
		return _OutputPtr.size();
	}

	//���ļ��а��ж�ȡ���ݲ���������ָ��
	template<class T>
	size_t	ReadVecByRowFromFile(thrust::host_vector<T> &_OutputPtr, const char *_FileName ,const boost::filesystem::path& _dir )
	{
		_OutputPtr.clear();
		char str[255];
		boost::filesystem::path _tmppath(_dir);
		_tmppath.append(_FileName);
		std::ifstream inFile(_tmppath.string());
		if (!inFile)
		{
			std::cout<<"File open error!\n";
			return false;
		}
		int _offset = 0;
		while(!inFile.eof())//��������ļ�ĩβ����һ������ֵ����û����ʱ����0����0�ͼ���~
		{
			inFile.getline(str,255);
			//strcpy_s(bstr,str);
			if (0 != strlen(str))
			{
				_OutputPtr.push_back((T)atof(str));
			}		
			_offset++;
		}
		return _OutputPtr.size();
	}

	template<class T>
	HBXDef::CMatrix<T>& GetMatFrmMtx( std::ifstream& _file, std::vector<std::string> &_vstrLine, int &_subMatDim )
	{
		using namespace std;
		using namespace boost;
		using namespace HBXDef;
		using namespace HBXFEMDef;

		std::string _stringLine;
		size_t _sum = 0;	//�Խ���ʽ��Ԫ�ؼ���
		std::vector<double> _vElement;
		HBXDef::CMatrix<T>	_tmpMat( _subMatDim, _subMatDim );

		do
		{
			_sum += _vstrLine.size();
			for ( int i=0; i<_vstrLine.size();i++ )
			{
				_vElement.push_back( lexical_cast<double>( _vstrLine[i] ) );
			}
			getline(_file, _stringLine);
			boost::char_separator<char> sep(" ,��*=");  
			typedef boost::tokenizer<boost::char_separator<char> >  CustonTokenizer;  
			CustonTokenizer tok(_stringLine, sep);
			// ����ָ���
			_vstrLine.clear();
			for(CustonTokenizer::iterator beg=tok.begin(); beg!=tok.end();++beg)  
			{  
				_vstrLine.push_back(*beg);  
			}  
		}while ( _sum<300 );

		_sum = 0;
		for (int i = 0; i<_subMatDim; i++)
		{
			for (int j = 0; j<i; j++)
			{
				_tmpMat(i,j) = _vElement[_sum];
				_sum++;
			}

		}
		return _tmpMat;
	}

	//Vector���������ڴ濽��
	template<class T>
	void		MemcpyArrayTransToVec( T* _SourceArray, size_t _length, std::vector<T>& _DestVec)
	{
		_DestVec.resize(_length);
		for (size_t i=0; i<_length; i++)
		{
			_DestVec[0] = ( *(_SourceArray+i) );
		}
	}

	//������Vector���ڴ濽��
	template<class T>
	int MemcpyVecTransToArray(std::vector<T> &_vecIn,  T  *_arrayOut)
	{
		if (nullptr == _arrayOut)
		{
			std::cout<<"\n ����δ�����ڴ�! \n"<<endl;
			return false;
		}
		size_t iMax = _vecIn.size();
		memcpy(_arrayOut, _vecIn.data(), sizeof(T)*iMax);
		return true;
	}


	//�������������
	template<class T>
	UserStatusError_t	OutputVectorData(std::vector<T> _OutputVec, const char *_FileName)
	{
		std::ofstream	fout(_FileName);
		if (!_FileName) return USER_EMPTY_POINT;
		if (!exists(_dir))//�ж�·���Ƿ����
		{
			std::cout<<"�洢·����Ч..."<<std::endl;
			return UserStatusError_t::USER_STATUS_INVALID_VALUE;
		}
		for (int i=0; i<_OutputVec.size(); i++)
		{
			fout<<_OutputVec[i]<<std::endl;
		}
		fout.close();
		return USER_STATUS_SUCCESS;
	}

	//������������
	template<class T>
	UserStatusError_t	OutputPointData(T*	_OutputPoint, size_t _length, const char *_ExportFileName, boost::filesystem::path _dir )
	{
		if (_OutputPoint == nullptr)
		{
			std::cout<<"����ָ��Ϊ��..."<<std::endl;
			return UserStatusError_t::USER_EMPTY_POINT;
		}
		if (!exists(_dir))//�ж�·���Ƿ����
		{
			std::cout<<"�洢·����Ч..."<<std::endl;
			return UserStatusError_t::USER_STATUS_INVALID_VALUE;
		}

		_dir.append(_ExportFileName);
		std::ofstream	fout(_dir.string());
		if (!_ExportFileName) return USER_EMPTY_POINT;
		for ( size_t i=0; i<_length; i++ )
		{
			fout<<*(_OutputPoint+i)<<std::endl;
		}
		fout.close();
		return UserStatusError_t::USER_STATUS_SUCCESS;
	}

	//��CSR��ʽ����������ֱ��¼�������ļ�
	//@_Val����������
	//@_ColSort����������������
	//@_ValRowSort����������������
	//@_nnA��_Val����
	//@_N�������size
	template<class T1, class T2>
	void		_ExportCSRdataToFile(T1 *_Val,			
									 T2 *_ColSort,		
									 T2 *_ValRowSort,
									 int _nnA, int _N,
									 const char*	_NoneZeroValFileName,
									 const char*	_ColSortFileName,
									 const char*	_NoneZeroRowSortFileName,
									 const boost::filesystem::path& _dir )
	{
		HBXDef::CheckUserDefErrors( HBXDef::OutputPointData( _Val, _nnA, _NoneZeroValFileName, _dir ) );
		HBXDef::CheckUserDefErrors( HBXDef::OutputPointData( _ColSort, _nnA, _ColSortFileName, _dir ) );
		HBXDef::CheckUserDefErrors( HBXDef::OutputPointData( _ValRowSort, (_N+1), _NoneZeroRowSortFileName, _dir ) );
	}

	//��CSR��ʽ����������ֱ��¼�������ļ�
	//@_vecVal��������������
	//@_vecColSort����������������
	//@_vecValRowSort����������������
	template<class T1, class T2>
	void _ExportCSRdataToFile(std::vector<T1>	&_vecVal,			
							  std::vector<T2> &_vecColSort,		
							  std::vector<T2> &_vecValRowSort,
							  const char*		_NoneZeroValFileName,
							  const char*		_ColSortFileName,
							  const char*		_NoneZeroRowSortFileName,
							  const boost::filesystem::path& _dir)
	{
		HBXDef::CheckUserDefErrors( HBXDef::OutputVectorData(_vecVal,		_NoneZeroValFileName) );
		HBXDef::CheckUserDefErrors( HBXDef::OutputVectorData(_vecColSort,		_ColSortFileName) );
		HBXDef::CheckUserDefErrors( HBXDef::OutputVectorData(_vecValRowSort,	_NoneZeroRowSortFileName) );
	}

	//����CSR��ʽ�ı������ݲ�������������
	//@_vecVal��������������
	//@_vecColSort����������������
	//@_vecValRowSort����������������
	template<class T1,class T2>
	size_t	_ImportCSRdataFromFile(std::vector<T1>		&_vecVal,			
								   std::vector<T2>	&_vecColSort,		
								   std::vector<T2>	&_vecValRowSort,
								   const char*		_NoneZeroValFileName,
								   const char*		_ColSortFileName,
								   const char*		_NoneZeroRowSortFileName,
								   const boost::filesystem::path& _dir)
	{
		_vecVal.clear();
		_vecColSort.clear();
		_vecValRowSort.clear();
		size_t _nz1 = HBXDef::ReadVecByRowFromFile(_vecVal,		_NoneZeroValFileName, _dir);
		size_t _nz2 = HBXDef::ReadVecByRowFromFile( _vecColSort,	_ColSortFileName, _dir);
		size_t _N = HBXDef::ReadVecByRowFromFile(_vecValRowSort,	_NoneZeroRowSortFileName, _dir);
		if (_nz1 == _nz2) return _nz1;
		else return 0;
	}

	//����CSR��ʽ�ı������ݲ�������������,thrust��
	//@_vecVal��������������
	//@_vecColSort����������������
	//@_vecValRowSort����������������
	template<class T1,class T2>
	size_t	_ImportCSRdataFromFile(	thrust::host_vector<T1> &_vecVal,			
									thrust::host_vector<T2> &_vecColSort,		
									thrust::host_vector<T2> &_vecValRowSort,
									const char*	_NoneZeroValFileName,
									const char*	_ColSortFileName,
									const char*	_NoneZeroRowSortFileName,
									const boost::filesystem::path& _dir)
	{
		_vecVal.clear();
		_vecColSort.clear();
		_vecValRowSort.clear();
		size_t _nz1 = HBXDef::ReadVecByRowFromFile(_vecVal,		_NoneZeroValFileName, _dir);
		size_t _nz2 = HBXDef::ReadVecByRowFromFile( _vecColSort,	_ColSortFileName, _dir);
		size_t _N = HBXDef::ReadVecByRowFromFile(_vecValRowSort,	_NoneZeroRowSortFileName, _dir);
		if (_nz1 == _nz2) return _nz1;
		else return 0;
	}


	//�������ľ���ת����CSR��ʽ
	//@_InCMatric	���������
	//@_dNonZeroVal CSRģʽ��һ��<T>����
	//@_iColSort	CSRģʽ�ڶ���int����
	//@_iNonZeroRowSort	CSRģʽ������int����
	//�μ�����CUDA�Ĵ��ģϡ������PCG�㷨�Ż�-֣��γ-�廪��ѧѧ��
	template<class _iT, class _oT>
	bool	_FullMatricToCSR( HBXDef::CMatrix<_iT>& _InCMatric, 
		std::vector<_oT>& _dNonZeroVal, 
		std::vector<int>& _iColSort, 
		std::vector<int>& _iNonZeroRowSort)
	{
		//У��CSR�����������Ƿ�Ϊ�գ���Ϊ�����������Ϊ�������
		if (false == _dNonZeroVal.empty())
		{
			_dNonZeroVal.clear();
		}
		if (false == _iColSort.empty())
		{
			_iColSort.clear();
		}
		if (false == _iNonZeroRowSort.empty())
		{
			_iNonZeroRowSort.clear();
		}
		//��ʱ������������
		int _tempRow,_tempCol;
		//��i��j��ֵ����ʱ����
		float _tempVal;
		//Ϊ��һ�еĵ�һ������Ԫ�ı�־λ
		bool _FlagNonZero;
		_InCMatric.GetSize(_tempRow, _tempCol);
		for (int i=0; i<_tempRow; i++)
		{
			_FlagNonZero = true;
			for (int j=0; j<_tempCol; j++)
			{
				_tempVal = _InCMatric(i,j);
				if (_tempVal !=0 )
				{
					_dNonZeroVal.push_back(_tempVal);
					_iColSort.push_back(j);
					if (true == _FlagNonZero)
					{
						_iNonZeroRowSort.push_back((int)_iColSort.size()-1);
						_FlagNonZero = false;
					}
				}
			}
		}
		if (_iNonZeroRowSort.size() == _tempRow)
		{
			//A_R���һ��������ʾ
			_iNonZeroRowSort.push_back((int)_iColSort.size());
			return true;
		}
		else return false;
	}
	
}

