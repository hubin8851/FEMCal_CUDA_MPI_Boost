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

	//��ʼ���������˿ڣ����ݿ⣬��Ҫ��ı����Ϣ
	//@host:	MYSQL������IP
	//@port��	���ݿ�˿�
	//@DataBase:���ݿ�����
	//@user��	���ݿ��û�
	//@passwd��	���ݿ��û�������
	//@charset��ϣ��ʹ�õ��ַ���
	//@Msg��		���ص���Ϣ������������Ϣ
	void Initial(const std::string _host="127.0.0.1", std::string _port="3306", std::string _Database="femconjugate", 
				std::string _user="root",std::string _passwd="admin", std::string _charset="GBK");
	//��Ҫ���ܣ��������ݿ�
	//�����Ƿ�ɹ������ص�boolֵ
	bool	ConnMySQL();
	//��������ʱ�䣬����ά�ȣ������������������ȣ������ʱ,����ֵ��mysql���Ĵ�����
	//@_matsz:	����ά��
	//@_iters:		��������
	//@_usedT:	�����ʱ
	int		InsertLog( const std::string& _Table, size_t _matsz, size_t _iters, float _usedT, float _prec = 1e-6 );
	//���������������ݲ���������
	//@SQL:		��ѯ��SQL���
	//@Cnum:	��ѯ������
	//@Msg:		���ص���Ϣ
	//@int		Ϊ���ص�����
	int		SelectData(std::string _Table="cgconjugate");
	//��ȡ���Id
	bool	GetMaxId(std::string _strDatabase, int& _MaxID);
	//���ݵ�ǰʱ���ȡ�ñ���ʼ����ID
	//@_DatabaseName:	���ݿ�����
	//@_SortNum:		���õ���������
	//@_SortTime:		���������ġ���ǰʱ�䡱����д��Ĭ�ϵ�ǰlocalʱ��
	bool GetStartID(const char*	const _DatabaseName, int&	_SortNum, char* _SortTime = nullptr);
	//����KEYֵ��ѯ����
	//@_Database:	���ݿ�����
	//@_MatSize:		������
	//@_OutFaultDig:�����FAULTDIG�ṹ��
	int SelectDataById(std::string			_Database,
		int&			_Id,
		HBXDef::CalRecord_t&	_CalRecord);
	//�޸����ݣ���ʱ����
	//@_Table����ѯ�����ݿ������
	//@_fault_key����������ID
	//@_CalRecord ����������¼
	int	ModifyData(std::string _Table, const int _fault_key,const HBXDef::CalRecord_t&	_CalRecord);
	//ɾ������
	//@_Table����ѯ�����ݿ��еı�
	//@_Id����������ID
	bool	DeleteData(std::string _Table, const int _Id);
	//�ر����ݿ�����
	void	CloseMySQL();
	//���ؼ���״̬�����߳�ʹ��
	bool	IsSorting();
	//�������ӳɹ���־λ
	bool	IsConnected();
	//���ñ�־λΪfalse
	void	ResetConnect();
private:
	bool	m_isConnect;
	bool	m_isSorting;//�Ƿ��ڼ����ı�־λ,���̱߳���
	std::string		m_OutStrTime;//��������ݿ��ʱ��

	MYSQL	mysql;
	MYSQL_FIELD	*fd;

	std::string m_host;	//IP��ַ
	std::string m_port;		//�˿�
	std::string m_Database;	//���ݿ�����
	std::string m_user;			//�û���root
	std::string m_passwd;	//���룬admin
	std::string m_charset;	//��GBK��
	std::string m_Msg;		//������Ϣ��¼

	std::vector<HBXDef::CalRecord_t>	m_vCalRecord;
};

