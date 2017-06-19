/************************************************************************/
/*								 HBX,2016-11-10															           */
/************************************************************************/
#include "stdafx.h"
#include "HbxFEMFunc.h"

//获得所有节点数目
int GetTotalNodeNum_f()
{ 
	return((int)g_vNode_f.size()); 
}
//获得所有单元数目
int GetTotalElementNum_f()
{ 
	return((int)g_vBaseElem_f.size()); 
}
//获取材料数目
int GetTotalMaterialNum_f()
{ 
	return((int)g_vMaterial_f.size()); 
}
//获取所有截面数目
int	GetTotalSectionNum_f()
{ 
	return((int)g_vSection_f.size()); 
}

//获得所有节点数目
int GetTotalNodeNum_d()
{ 
	return((int)g_vNode_d.size()); 
}
//获得所有单元数目
int GetTotalElementNum_d()
{ 
	return((int)g_vBaseElem_d.size()); 
}
//获取材料数目
int GetTotalMaterialNum_d()
{ 
	return((int)g_vMaterial_d.size()); 
}
//获取所有截面数目
int	GetTotalSectionNum_d()
{ 
	return((int)g_vSection_d.size()); 
}

//true为单元，false为节点
 int GetTotalSetNum(bool _ElOrN)
 {
	 if (true ==_ElOrN)
	 {
		 return((int)g_vElSet.size()); 
	 }
	 else return((int)g_vNSet.size()); 
 }
