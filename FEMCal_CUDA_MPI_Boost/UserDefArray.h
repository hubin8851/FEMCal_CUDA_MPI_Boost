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
	//Ĭ�Ϲ��캯��
	UserDefArray(int = 0);	
	//��������
	UserDefArray( const UserDefArray& );
	virtual ~UserDefArray(){ if (values) { FreeVal(values) }; }
	//���������
	UserDefArray& operator=(const UserDefArray& _rhs){};
	UserDefArray& operator=(const UserDefArray&& _rhs ){} ;

	//����ĳ������
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

	//�����ٶ��Ϸ���ռ䡣�ڴ˴��ǳ���Ҫ����Ϊ��������ջ�Ϸ��䣬��Ϊ��װ������ܴܺ󣬵���ջ���
};

}