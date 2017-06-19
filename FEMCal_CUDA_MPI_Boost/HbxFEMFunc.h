#ifndef _HBXFEMFUNC_H_
#define _HBXFEMFUNC_H_

extern	HBXDef::PrecsType_t	g_PrecsType;	//全局变量，设置数据读取和计算的精度(float or double)

extern	HBXFEMDef::NodeInMap_f		g_vNode_f;		//全局节点容器,切记！是map！为了搜索更快
extern	HBXFEMDef::_vBaseElem_f		g_vBaseElem_f;	//全局基单元容器
extern	HBXFEMDef::_vMaterial_f		g_vMaterial_f;	//全局材料容器
extern	HBXFEMDef::_vSection_f		g_vSection_f;	//全局截面容器
extern	HBXFEMDef::_vSet			g_vNSet;//节点集合容器
extern	HBXFEMDef::_vSet			g_vElSet;//单元集合容器

extern	HBXFEMDef::NodeInMap_d		g_vNode_d;		//全局节点容器,切记！是map！为了搜索更快
extern	HBXFEMDef::_vBaseElem_d		g_vBaseElem_d;	//全局基单元容器
extern	HBXFEMDef::_vMaterial_d		g_vMaterial_d;	//全局材料容器
extern	HBXFEMDef::_vSection_d		g_vSection_d;	//全局截面容器


#pragma region 有限元数据读取相关

/************************************************************************/
/*           节点相关的函数                                             */
/************************************************************************/
//获得所有节点数目
int GetTotalNodeNum_f();
int GetTotalNodeNum_d();
//将该节点的三个坐标值作为一个节点参数放入容器
//@_datainput:输入参数(坐标)的指针(首地址)
//	int gAdd2NodeVec(float *_datainput);

/************************************************************************/
/*           单元相关的函数                                             */
/************************************************************************/
//获得所有单元数目
int GetTotalElementNum_f();
int GetTotalElementNum_d();
//设置单元属性
//@kNum:在容器中的索引
//@cset:设置当前单元的截面
//@cm:设置单元材料
int gSet2D3NodeElemParam(int kNum,int cset,int cm);
//获取单元属性
//@kNum:在容器中的索引
//@kEle:获取当前单元的三个节点号
//@ksection:获取当前单元截面
//@kme:获取当前单元材料
int gGet2D3NodeElemParam(int kNum,int *kEle,int &ksection,int &kme);
//获取单元内节点坐标
//@_outXY：输出节点坐标的一维指针
template <class T>
int gGetNodeCoorXInElem(int kNodeNum,T *_outXY);
//将该单元的节点号和截面号放入容器
//@iaNode:输入参数(节点)的指针(首地址)
//@_Section:截面号
int gAdd2ElementVec(int *iaNode, int _Section);
//将该单元的节点号和截面号放入全局容器
//@iaNode:输入参数(节点)的指针(首地址)
//@_Section:截面号
//起始节点号为1...
int gAdd2ElementVecStartFrom1(int *iaNode, int _Section);
/************************************************************************/
/*           材料相关的函数                                             */
/************************************************************************/

//获取材料数目
int GetTotalMaterialNum_f();
int GetTotalMaterialNum_d();
//设置材料属性并添加进全局容器
//@sE：当前材料的杨氏模量
//@sMu:当前材料的泊松比
//@sGama:重度
//@ch：材料名称
//@_dAlpha：线膨胀系数，默认为0可不填
int gAdd2MaterialVec(float &sE, float &sMu, float &sGama, char *ch, float _dAlpha=0);
//获取当前材料索引的名称
//@kNum:在容器中的索引
//@sE：当前材料的杨氏模量
//@sMu:当前材料的泊松比
template <class T>
int		gGetMaterialInfo(int kNum,T &sE,T &sMu);
//获取当前材料索引的名称
//@kNum:在容器中的索引
//@ch:当前索引的材料名称
char*	gGetMaterialName(int kNum,char *ch);

int GetTotalSetNum(bool _ElOrN = true);
/************************************************************************/
/*           截面相关的函数                                             */
/************************************************************************/
//获取所有截面数目
int	GetTotalSectionNum_f();
int	GetTotalSectionNum_d();
//获取截面厚度
//@kNum:在容器中的索引
//@st:当前索引的截面厚度
template <class T>
T gGetSectionThick(int kNum,T &st);
//获取截面厚度
//@kNum:在容器中的索引
//@scname：截面名称
//@maname:截面材料名称
int gGetSectionData(int kNum, char *scname, char *maname);
//对截面参数初始化并赋值
//@scname:截面名称
//@maname:材料名称
int gAdd2SectionVec(char *SectionName, char *MaterialName);
//设置截面厚度，添加至容器的最后一个元素内
//@st:截面厚度
int gSetSectionThick(float st);

/************************************************************************/
/* 刚度矩阵计算类                                                       */
/************************************************************************/
//单元刚度矩阵存储
//@m_line
//@m_row
//@kij
//@mstiffness
template <class T>
int gElementStiffMatrixStorage(int m_line,int m_row,T *kij,T *mstiffness);
//HBX,2015-2-4
//创建某个单元的刚度矩阵
//@kNum:输入的单元内的第kNum个节点
//@mstiffness：输出刚度矩阵指针
//@ipmark：所在总刚矩阵中的位置
int gElementStiffMatrixNew(int kNum, float *mstiffness,int *ipmark);



#pragma endregion 有限元数据读取相关


#endif