#pragma once
#include <mysql/mysql.h>

class CCntToMysql
{
public:
	CCntToMysql(void);
	~CCntToMysql(void);

#ifdef _DEBUG
	friend bool runTest(int argc, char **argv); 
#endif

	//初始化主机，端口，数据库，所要填的表等信息
	//@host:	MYSQL服务器IP
	//@port：	数据库端口
	//@DataBase:数据库名称
	//@user：	数据库用户
	//@passwd：	数据库用户的密码
	//@charset：希望使用的字符集
	//@Msg：		返回的消息，包括错误消息
	void Initial(const std::string _host="127.0.0.1", std::string _port="3306", std::string _Database="femconjugate", 
				std::string _user="root",std::string _passwd="admin", std::string _charset="GBK");
	//主要功能：连接数据库
	//连接是否成功看返回的bool值
	bool	ConnMySQL();
	//插入计算的时间，矩阵维度，迭代次数，收敛精度，计算耗时,返回值是mysql语句的错误码
	//@_matsz:	矩阵维度
	//@_iters:		迭代次数
	//@_usedT:	计算耗时
	int		InsertLog( const std::string& _Table, size_t _matsz, size_t _iters, float _usedT, float _prec = 1e-6 );
	//遍历表内所有数据并放入容器
	//@SQL:		查询的SQL语句
	//@Cnum:	查询的列数
	//@Msg:		返回的消息
	//@int		为返回的数据
	int		SelectData(std::string _Table="cgconjugate");
	//获取最大Id
	bool	GetMaxId(std::string _strDatabase, int& _MaxID);
	//根据当前时间获取该表起始索引ID
	//@_DatabaseName:	数据库名称
	//@_SortNum:		所得到的索引号
	//@_SortTime:		所需索引的“当前时间”，不写则默认当前local时间
	bool GetStartID(const char*	const _DatabaseName, int&	_SortNum, char* _SortTime = nullptr);
	//根据KEY值轮询数据
	//@_Database:	数据库名称
	//@_MatSize:		主键号
	//@_OutFaultDig:输出的FAULTDIG结构体
	int SelectDataById(std::string			_Database,
		int&			_Id,
		HBXDef::CalRecord_t&	_CalRecord);
	//修改数据，暂时不用
	//@_Table：查询的数据库表名称
	//@_fault_key：测试用例ID
	//@_CalRecord 测试用例记录
	int	ModifyData(std::string _Table, const int _fault_key,const HBXDef::CalRecord_t&	_CalRecord);
	//删除数据
	//@_Table：查询的数据库中的表
	//@_Id：故障序列ID
	bool	DeleteData(std::string _Table, const int _Id);
	//关闭数据库连接
	void	CloseMySQL();
	//返回检索状态，多线程使用
	bool	IsSorting();
	//返回连接成功标志位
	bool	IsConnected();
	//重置标志位为false
	void	ResetConnect();
private:
	bool	m_isConnect;
	bool	m_isSorting;//是否在检索的标志位,多线程备用
	std::string		m_OutStrTime;//输出到数据库的时间

	MYSQL	mysql;
	MYSQL_FIELD	*fd;

	std::string m_host;	//IP地址
	std::string m_port;		//端口
	std::string m_Database;	//数据库名称
	std::string m_user;			//用户，root
	std::string m_passwd;	//密码，admin
	std::string m_charset;	//“GBK”
	std::string m_Msg;		//错误消息记录

	std::vector<HBXDef::CalRecord_t>	m_vCalRecord;
};

