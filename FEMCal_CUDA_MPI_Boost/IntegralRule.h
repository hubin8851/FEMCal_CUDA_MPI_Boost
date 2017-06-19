#pragma once

#include "HBXFEMDefStruct.h"

namespace HBXFEMDef
{
	/* ��Բ�ͬ��Ԫ���ж��ֻ�����ʽ����ȫ/��������˹/ţ��-����/��������ɭ �ȣ�������������
	�㷨�ڣ���������ȷ�Ļ��ֵ��Ȩֵ�����ڶ���ѧ�ڲ��ϵ����Իᷢ���仯�������ò��ϵ�shared_ptr������
	ÿ����Ԫ���Ӧ
	*/
	class CIntegralRule
	{
	private:
	protected:
		boost::shared_ptr<HBXFEMDef::_BaseElem<double>> _elem;
		bool	isDynamic;//����ѧ��־λ������Ƕ���ѧ����¼֮ǰ�Ļ��ֵ㣻����ѧ�����¼��˹���ֵ�����������еı仯���

		virtual int SetUpPointOnLine( int iPoint, MaterialType_t _material_t )=0;
	public:
		CIntegralRule(void);
		virtual ~CIntegralRule(void);

		//
		int SetUpIntegrationPoint( Integral_t _Int_type  );
	};
}


