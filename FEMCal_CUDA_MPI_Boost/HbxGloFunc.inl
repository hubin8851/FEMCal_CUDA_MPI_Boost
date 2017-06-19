

//命名空间内部分函数
namespace HBXDef
{

	//从文件中按行读取数据并放入数组指针
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
			std::cout<<"读取文件列向量出错!\n";
			return false;
		}
		int _offset = 0;
		while(!inFile.eof())//如果到达文件末尾返回一个非零值，即没到达时返回0，非0就继续~
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

	//从文件中按行读取数据并放入数组指针
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
		while(!inFile.eof())//如果到达文件末尾返回一个非零值，即没到达时返回0，非0就继续~
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
		size_t _sum = 0;	//对角形式的元素计数
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
			boost::char_separator<char> sep(" ,，*=");  
			typedef boost::tokenizer<boost::char_separator<char> >  CustonTokenizer;  
			CustonTokenizer tok(_stringLine, sep);
			// 输出分割结果
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

	//Vector至向量的内存拷贝
	template<class T>
	void		MemcpyArrayTransToVec( T* _SourceArray, size_t _length, std::vector<T>& _DestVec)
	{
		_DestVec.resize(_length);
		for (size_t i=0; i<_length; i++)
		{
			_DestVec[0] = ( *(_SourceArray+i) );
		}
	}

	//向量至Vector的内存拷贝
	template<class T>
	int MemcpyVecTransToArray(std::vector<T> &_vecIn,  T  *_arrayOut)
	{
		if (nullptr == _arrayOut)
		{
			std::cout<<"\n 数组未分配内存! \n"<<endl;
			return false;
		}
		size_t iMax = _vecIn.size();
		memcpy(_arrayOut, _vecIn.data(), sizeof(T)*iMax);
		return true;
	}


	//输出容器的数据
	template<class T>
	UserStatusError_t	OutputVectorData(std::vector<T> _OutputVec, const char *_FileName)
	{
		std::ofstream	fout(_FileName);
		if (!_FileName) return USER_EMPTY_POINT;
		if (!exists(_dir))//判断路径是否存在
		{
			std::cout<<"存储路径无效..."<<std::endl;
			return UserStatusError_t::USER_STATUS_INVALID_VALUE;
		}
		for (int i=0; i<_OutputVec.size(); i++)
		{
			fout<<_OutputVec[i]<<std::endl;
		}
		fout.close();
		return USER_STATUS_SUCCESS;
	}

	//输出数组的数据
	template<class T>
	UserStatusError_t	OutputPointData(T*	_OutputPoint, size_t _length, const char *_ExportFileName, boost::filesystem::path _dir )
	{
		if (_OutputPoint == nullptr)
		{
			std::cout<<"输入指针为空..."<<std::endl;
			return UserStatusError_t::USER_EMPTY_POINT;
		}
		if (!exists(_dir))//判断路径是否存在
		{
			std::cout<<"存储路径无效..."<<std::endl;
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

	//将CSR格式的三个数组分别记录进三个文件
	//@_Val：非零向量
	//@_ColSort：非零向量列索引
	//@_ValRowSort：非零向量行索引
	//@_nnA：_Val长度
	//@_N：矩阵的size
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

	//将CSR格式的三个数组分别记录进三个文件
	//@_vecVal：非零向量数组
	//@_vecColSort：非零向量列索引
	//@_vecValRowSort：非零向量行索引
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

	//导入CSR格式文本的数据并存入三个数组
	//@_vecVal：非零向量数组
	//@_vecColSort：非零向量列索引
	//@_vecValRowSort：非零向量行索引
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

	//导入CSR格式文本的数据并存入三个数组,thrust版
	//@_vecVal：非零向量数组
	//@_vecColSort：非零向量列索引
	//@_vecValRowSort：非零向量行索引
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


	//将完整的矩阵转换成CSR形式
	//@_InCMatric	输入矩阵类
	//@_dNonZeroVal CSR模式第一个<T>容器
	//@_iColSort	CSR模式第二个int容器
	//@_iNonZeroRowSort	CSR模式第三个int容器
	//参见基于CUDA的大规模稀疏矩阵的PCG算法优化-郑经纬-清华大学学报
	template<class _iT, class _oT>
	bool	_FullMatricToCSR( HBXDef::CMatrix<_iT>& _InCMatric, 
		std::vector<_oT>& _dNonZeroVal, 
		std::vector<int>& _iColSort, 
		std::vector<int>& _iNonZeroRowSort)
	{
		//校验CSR的三个容器是否为空，若为空则继续，不为空则清空
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
		//临时变量，行列数
		int _tempRow,_tempCol;
		//第i行j列值的临时变量
		float _tempVal;
		//为第一行的第一个非零元的标志位
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
			//A_R会多一个数，表示
			_iNonZeroRowSort.push_back((int)_iColSort.size());
			return true;
		}
		else return false;
	}
	
}

