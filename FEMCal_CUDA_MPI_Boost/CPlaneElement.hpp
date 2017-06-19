#pragma once
#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//ȫ�ֱ������������ݶ�ȡ�ͼ���ľ���(float or double)

namespace HBXFEMDef
{
	//�˵�Ԫ
	template <class T>
	class _PlaneElement:
		public _BaseElem<T>
	{
	protected:
		int iSection;	//��Ԫ����������
		T dLength;	//��Ԫ����
		T dSin, dCos;	//��Ԫ�ֲ�����x������������x��ļнǵ�������
		T daEndInterForce[6];	//��Ԫ�˶�������

	public:
		void StiffMatCal(){};
		void MassMatCal( MassType_t _MassType ){};
		void SetMaterial( std::pair< boost::shared_ptr<_Material<T>>, boost::shared_ptr<_BaseSection<T>> > _MatAndSec ){ };
		T VolumeCalc(){return 0;};

	};

	typedef boost::shared_ptr<_PlaneElement<float>> _PlaneElement_f_ptr;
	typedef boost::shared_ptr<_PlaneElement<double>> _PlaneElement_d_ptr;
}
