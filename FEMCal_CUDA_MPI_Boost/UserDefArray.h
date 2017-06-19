#include "stdafx.h"

namespace HBXDef
{

template< class T= double>
class UserDefArray
{
	typedef UserDefArray<T>	_SameArray;
protected:
	int size;
	T*	values;
public:
	//默认构造函数
	UserDefArray(int = 0);	
	//拷贝构造
	UserDefArray( const UserDefArray& );
	virtual ~UserDefArray(){ if (values) { FreeVal(values) }; }
	//清除并拷贝
	UserDefArray& operator=(const UserDefArray& _rhs){};
	UserDefArray& operator=(const UserDefArray&& _rhs ){} ;

	//查找某项数据
	const T& operator()(int i) const
	{
#ifdef DEBUG
		if (i>=size)
		{
			CheckUserDefErrors(UserStatusError::USER_STATUS_INVALID_VALUE);
		}
#endif
		return values[i];
	}

	//重新再堆上分配空间。在此处非常重要，因为不能再在栈上分配，因为总装矩阵可能很大，导致栈溢出
};

}