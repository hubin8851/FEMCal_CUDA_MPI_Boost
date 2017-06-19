#ifndef _HBXFEMFUNC_H_
#define _HBXFEMFUNC_H_

extern	HBXDef::PrecsType_t	g_PrecsType;	//ȫ�ֱ������������ݶ�ȡ�ͼ���ľ���(float or double)

extern	HBXFEMDef::NodeInMap_f		g_vNode_f;		//ȫ�ֽڵ�����,�мǣ���map��Ϊ����������
extern	HBXFEMDef::_vBaseElem_f		g_vBaseElem_f;	//ȫ�ֻ���Ԫ����
extern	HBXFEMDef::_vMaterial_f		g_vMaterial_f;	//ȫ�ֲ�������
extern	HBXFEMDef::_vSection_f		g_vSection_f;	//ȫ�ֽ�������
extern	HBXFEMDef::_vSet			g_vNSet;//�ڵ㼯������
extern	HBXFEMDef::_vSet			g_vElSet;//��Ԫ��������

extern	HBXFEMDef::NodeInMap_d		g_vNode_d;		//ȫ�ֽڵ�����,�мǣ���map��Ϊ����������
extern	HBXFEMDef::_vBaseElem_d		g_vBaseElem_d;	//ȫ�ֻ���Ԫ����
extern	HBXFEMDef::_vMaterial_d		g_vMaterial_d;	//ȫ�ֲ�������
extern	HBXFEMDef::_vSection_d		g_vSection_d;	//ȫ�ֽ�������


#pragma region ����Ԫ���ݶ�ȡ���

/************************************************************************/
/*           �ڵ���صĺ���                                             */
/************************************************************************/
//������нڵ���Ŀ
int GetTotalNodeNum_f();
int GetTotalNodeNum_d();
//���ýڵ����������ֵ��Ϊһ���ڵ������������
//@_datainput:�������(����)��ָ��(�׵�ַ)
//	int gAdd2NodeVec(float *_datainput);

/************************************************************************/
/*           ��Ԫ��صĺ���                                             */
/************************************************************************/
//������е�Ԫ��Ŀ
int GetTotalElementNum_f();
int GetTotalElementNum_d();
//���õ�Ԫ����
//@kNum:�������е�����
//@cset:���õ�ǰ��Ԫ�Ľ���
//@cm:���õ�Ԫ����
int gSet2D3NodeElemParam(int kNum,int cset,int cm);
//��ȡ��Ԫ����
//@kNum:�������е�����
//@kEle:��ȡ��ǰ��Ԫ�������ڵ��
//@ksection:��ȡ��ǰ��Ԫ����
//@kme:��ȡ��ǰ��Ԫ����
int gGet2D3NodeElemParam(int kNum,int *kEle,int &ksection,int &kme);
//��ȡ��Ԫ�ڽڵ�����
//@_outXY������ڵ������һάָ��
template <class T>
int gGetNodeCoorXInElem(int kNodeNum,T *_outXY);
//���õ�Ԫ�Ľڵ�źͽ���ŷ�������
//@iaNode:�������(�ڵ�)��ָ��(�׵�ַ)
//@_Section:�����
int gAdd2ElementVec(int *iaNode, int _Section);
//���õ�Ԫ�Ľڵ�źͽ���ŷ���ȫ������
//@iaNode:�������(�ڵ�)��ָ��(�׵�ַ)
//@_Section:�����
//��ʼ�ڵ��Ϊ1...
int gAdd2ElementVecStartFrom1(int *iaNode, int _Section);
/************************************************************************/
/*           ������صĺ���                                             */
/************************************************************************/

//��ȡ������Ŀ
int GetTotalMaterialNum_f();
int GetTotalMaterialNum_d();
//���ò������Բ���ӽ�ȫ������
//@sE����ǰ���ϵ�����ģ��
//@sMu:��ǰ���ϵĲ��ɱ�
//@sGama:�ض�
//@ch����������
//@_dAlpha��������ϵ����Ĭ��Ϊ0�ɲ���
int gAdd2MaterialVec(float &sE, float &sMu, float &sGama, char *ch, float _dAlpha=0);
//��ȡ��ǰ��������������
//@kNum:�������е�����
//@sE����ǰ���ϵ�����ģ��
//@sMu:��ǰ���ϵĲ��ɱ�
template <class T>
int		gGetMaterialInfo(int kNum,T &sE,T &sMu);
//��ȡ��ǰ��������������
//@kNum:�������е�����
//@ch:��ǰ�����Ĳ�������
char*	gGetMaterialName(int kNum,char *ch);

int GetTotalSetNum(bool _ElOrN = true);
/************************************************************************/
/*           ������صĺ���                                             */
/************************************************************************/
//��ȡ���н�����Ŀ
int	GetTotalSectionNum_f();
int	GetTotalSectionNum_d();
//��ȡ������
//@kNum:�������е�����
//@st:��ǰ�����Ľ�����
template <class T>
T gGetSectionThick(int kNum,T &st);
//��ȡ������
//@kNum:�������е�����
//@scname����������
//@maname:�����������
int gGetSectionData(int kNum, char *scname, char *maname);
//�Խ��������ʼ������ֵ
//@scname:��������
//@maname:��������
int gAdd2SectionVec(char *SectionName, char *MaterialName);
//���ý����ȣ���������������һ��Ԫ����
//@st:������
int gSetSectionThick(float st);

/************************************************************************/
/* �նȾ��������                                                       */
/************************************************************************/
//��Ԫ�նȾ���洢
//@m_line
//@m_row
//@kij
//@mstiffness
template <class T>
int gElementStiffMatrixStorage(int m_line,int m_row,T *kij,T *mstiffness);
//HBX,2015-2-4
//����ĳ����Ԫ�ĸնȾ���
//@kNum:����ĵ�Ԫ�ڵĵ�kNum���ڵ�
//@mstiffness������նȾ���ָ��
//@ipmark�������ܸվ����е�λ��
int gElementStiffMatrixNew(int kNum, float *mstiffness,int *ipmark);



#pragma endregion ����Ԫ���ݶ�ȡ���


#endif