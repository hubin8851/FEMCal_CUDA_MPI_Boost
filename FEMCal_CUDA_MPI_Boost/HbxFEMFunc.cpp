/************************************************************************/
/*								 HBX,2016-11-10															           */
/************************************************************************/
#include "stdafx.h"
#include "HbxFEMFunc.h"

//������нڵ���Ŀ
int GetTotalNodeNum_f()
{ 
	return((int)g_vNode_f.size()); 
}
//������е�Ԫ��Ŀ
int GetTotalElementNum_f()
{ 
	return((int)g_vBaseElem_f.size()); 
}
//��ȡ������Ŀ
int GetTotalMaterialNum_f()
{ 
	return((int)g_vMaterial_f.size()); 
}
//��ȡ���н�����Ŀ
int	GetTotalSectionNum_f()
{ 
	return((int)g_vSection_f.size()); 
}

//������нڵ���Ŀ
int GetTotalNodeNum_d()
{ 
	return((int)g_vNode_d.size()); 
}
//������е�Ԫ��Ŀ
int GetTotalElementNum_d()
{ 
	return((int)g_vBaseElem_d.size()); 
}
//��ȡ������Ŀ
int GetTotalMaterialNum_d()
{ 
	return((int)g_vMaterial_d.size()); 
}
//��ȡ���н�����Ŀ
int	GetTotalSectionNum_d()
{ 
	return((int)g_vSection_d.size()); 
}

//trueΪ��Ԫ��falseΪ�ڵ�
 int GetTotalSetNum(bool _ElOrN)
 {
	 if (true ==_ElOrN)
	 {
		 return((int)g_vElSet.size()); 
	 }
	 else return((int)g_vNSet.size()); 
 }
