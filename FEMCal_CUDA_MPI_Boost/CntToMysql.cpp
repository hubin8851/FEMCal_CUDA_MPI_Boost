#include "StdAfx.h"
#include "CntToMysql.h"


CCntToMysql::CCntToMysql(void)
{
}


CCntToMysql::~CCntToMysql(void)
{
}

//��ʼ���������˿ڣ����ݿ⣬��Ҫ��ı����Ϣ
void CCntToMysql::Initial(const std::string _host, std::string _port, 
										std::string _Database, std::string _user, std::string _passwd, std::string _charset)
{
	m_host = _host;
	m_port = _port;
	m_Database = _Database;
	m_user = _user;
	m_passwd = _passwd;
	m_charset = _charset;

	m_isConnect = false;
	m_isSorting = false;
}

//��ʼ������
bool	CCntToMysql::ConnMySQL()
{
	//����������ֱ�ӷ���
	if (IsConnected())
	{
		std::cout<<"�Ѿ�����"<<std::endl;
		return true;
	}
	//string��char*ת��
	char *host = (char*)m_host.c_str();
	char *user = (char*)m_user.c_str();
	char *port = (char*)m_port.c_str();
	char *passwd = (char*)m_passwd.c_str();
	char *dbname = (char*)m_Database.c_str();
	char *charset = (char*)m_charset.c_str();
	if( NULL == mysql_init(&mysql) )
	{
		m_Msg ="Error: Initial mysql handle error";
		return false;
	}
	if( NULL == mysql_real_connect(&mysql, host, user, passwd, dbname, 0, NULL, 0) )
	{
		m_Msg ="Error: Failed to connect to database";
		return false;
	}

	if( 0 != mysql_set_character_set(&mysql,"GBK") )
	{
		m_Msg = "Error: mysql_set_character_set";
		return false;
	}
	m_isConnect = true;
	return true;
}

//��������ʱ�䣬����ά�ȣ������������������ȣ������ʱ
//@_matsz:	����ά��
//@_iters:		��������
//@_usedT:	�����ʱ
int		CCntToMysql::InsertLog( const std::string&	_Table, size_t _matsz, size_t _iters, float _usedT, float _prec )
{
	using namespace boost;
	using namespace std;
	int _err_c;
	boost::format	_fmt( "insert into %s (mat_size, iters, precs,usedT)values(%s, %s, %s,%s)");
	_fmt%_Table % _matsz % _iters % _prec % _usedT;
	_err_c = mysql_query(&mysql,  (char*)_fmt.str().c_str());
	if(0 != _err_c)
	{
		std::cout<<"Error: ���붯������ʧ��... "<<std::endl;	//��ʼ�����ﲻ�����жϵģ���������жϡ�ֻ��Ϊ�˷����Ժ�boost::log��
		return _err_c;
	}	
	return _err_c;
}

//��ѯ����
int CCntToMysql::SelectData(std::string _Table)
{
	MYSQL_ROW	m_row;
	MYSQL_RES*	m_res;
	char column[32][32];
	int i,j;
	std::string _tempstringAll = "SELECT * from ";
	_tempstringAll.append(_Table);
	if(mysql_query(&mysql, (char*)_tempstringAll.c_str()) != 0)
	{
		m_Msg = "Error: SelectData��������...";
		return false;
	}
	//�����ѯ�����������ӵ�m_res
	m_res = mysql_store_result(&mysql);

	if(m_res != NULL)
	{
		std::cout<<"number of result: "<<(unsigned long)mysql_num_rows(m_res)<<std::endl;
		//��ȡ����
		for (i=0; fd=mysql_fetch_field(m_res);i++)
		{
			strcpy(column[i],fd->name);
		}
		j = mysql_num_fields(m_res);//��ȡ��¼�����ֶε�����
		for (i=0; i<j; i++)
		{
			std::cout<<column[i]<<" ";
		}
		std::cout<<std::endl;

		//��ȡ��������
		HBXDef::CalRecord_t _pCalRecord;
		while(m_row = mysql_fetch_row(m_res))
		{
			_pCalRecord._matsize = boost::lexical_cast<int>(m_row[2]);
			_pCalRecord._iters = boost::lexical_cast<int>(m_row[3]);
			_pCalRecord._tol = boost::lexical_cast<float>(m_row[4]);
			_pCalRecord._eigfilename = m_row[4];
			m_vCalRecord.push_back(_pCalRecord);
		}
	}
	else
	{
		m_Msg = "Error: store to m_res error! ";
		ResetConnect();
		return false;
	}
	mysql_free_result(m_res);
	return true;
}

//��ȡ���Id
bool	CCntToMysql::GetMaxId( std::string _strDatabase, int& _MaxID)
{
	MYSQL_ROW	m_row;
	MYSQL_RES*	m_res;
	char SqlCommand[128];
	std::sprintf(SqlCommand, "select max(Id) from %s", (char*)_strDatabase.c_str());
	if(mysql_query(&mysql, SqlCommand) != 0)
	{
		std::cout<<"Error: ��ȡ���Id��ʧ��,������ID-(Id)��Сд��������... "<<std::endl;
		return false;
	}
	m_res = mysql_store_result(&mysql);
	if (0 == m_res->row_count)
	{
		std::cout<<"Error: ��ȡ��ʼ����IDʧ�ܣ�û���ҵ������¼... "<<std::endl;
		_MaxID = 0;
		return true;
	}
	m_row = mysql_fetch_row(m_res);
	_MaxID = boost::lexical_cast<int>(m_row[0]);
	return true;
}

//���ݵ�ǰʱ���ȡ�ñ���ʼ����ID��֮ǰ�ļ�¼���ӡ�
bool CCntToMysql::GetStartID(const char*	const	_DatabaseName,
								int&				_SortNum,
								char*				_SortTime)
{
	using namespace boost::gregorian;
	using namespace boost::posix_time;
	m_isSorting = true;	//��ʼ������־λ�����߳�ʹ�ã��ݲ���mutex
	char SqlCommand[200];
	MYSQL_ROW	m_row;
	MYSQL_RES*	m_res;
	ptime _NowTime = second_clock::local_time();//��ȡ��ǰʱ��
	int _tempi;
	if (nullptr == _SortTime)
	{
		_SortTime = (char*)to_simple_string(_NowTime).c_str();
	}
	std::sprintf(SqlCommand, "select ID from %s where CalTime = '%s' ", 
		_DatabaseName,			_SortTime);
	if(mysql_query(&mysql, SqlCommand) != 0)	//���ɹ�
	{
		std::cout<<"Error: ��ȡ����ʼ����IDʧ��... "<<std::endl;
		m_isSorting = false;
		return false;
	};
	m_res = mysql_store_result(&mysql);
	//���ü�¼��Ϊ��ʱѡ������һ��KEYֵ��Ϊ�������ΪKEY��������
	if (0 == m_res->row_count)
	{
		std::cout<<"Error: ��ȡ��ʼ����IDʧ�ܣ�û���ҵ������¼... "<<std::endl;
		std::cout<<"Error: Ĭ��ѡ�����һ����¼��ID��Ϊ������ʼֵ... "<<std::endl;
		std::sprintf(SqlCommand, "select ID from %s", _DatabaseName);
		if(mysql_query(&mysql, SqlCommand) != 0)
		{
			std::cout<<"Error: ��ȡ����ʼ����IDʧ��... "<<std::endl;
			m_isSorting = false;
			return false;
		};
		m_res = mysql_store_result(&mysql);
		while (m_row = mysql_fetch_row(m_res))
		{
			_tempi = boost::lexical_cast<int>(m_row[0]);
		}
		if(0 >= _tempi)
		{
			_SortNum = 0;
		}
		else _SortNum = _tempi;
		std::cout<<"�õ���ʼ����ֵΪ:"<<_SortNum<<std::endl;
		m_isSorting = false;
		return true;
	}
	else
	{
		m_row = mysql_fetch_row(m_res);
		_SortNum = boost::lexical_cast<int>(m_row[0]);
		std::cout<<"�õ���ʼ����ֵΪ:"<<_SortNum<<std::endl;
		m_isSorting = false;
		return true;
	}
}

//����KEYֵ��ѯ����
int CCntToMysql::SelectDataById(std::string			_Database,
									int&			_Id,
									HBXDef::CalRecord_t&	_CalRecord)
{
	m_isSorting = true;
	char SQL[200];
	MYSQL_ROW	m_row;
	MYSQL_RES*	m_res;
	std::sprintf(SQL, "select * from %s where ID = %d ",
		(char*)_Database.c_str(),			_Id);
	if(mysql_query(&mysql, SQL) != 0)
	{
		std::cout<<"Error: ��Idѡ������ʧ��... "<<std::endl;
		m_isSorting = false;
		return _Id;
	};
	//�����ѯ�����������ӵ�m_res
	m_res = mysql_store_result(&mysql);
	if (0 == m_res->row_count)
	{
		m_isSorting = false;
		return _Id;
	}
	m_row = mysql_fetch_row(m_res);
	_CalRecord._matsize = boost::lexical_cast<int>(m_row[2]);
	_CalRecord._iters = boost::lexical_cast<int>(m_row[3]);
	_CalRecord._tol = boost::lexical_cast<float>(m_row[4]);
	m_isSorting = false;
	return 0;
}

//ɾ������
//@_Table����ѯ�����ݿ��еı�
//@_Id����������ID
bool		CCntToMysql::DeleteData(std::string _Table, const int _Id)
{
	char SqlCommand[200];
	std::sprintf(SqlCommand, "delete from %s where FAULT_KEY = %d",
		(char*)_Table.c_str(), _Id );
	if(mysql_query(&mysql, SqlCommand) != 0)
	{
		m_Msg = "Error: delete Data... ";
		return false;
	}
	return true;
}

void	CCntToMysql::CloseMySQL()
{
	mysql_close(&mysql);
}

bool	CCntToMysql::IsSorting()
{
	return m_isSorting;
};

bool	CCntToMysql::IsConnected()
{
	return m_isConnect;
}

void	CCntToMysql::ResetConnect()
{
	if (!m_isConnect)
	{
		std::cout<<"�ڴ�֮ǰ�Ѿ�������..."<<std::endl;
		return;
	}
	m_isConnect = false;
}