#pragma once
#include "stdafx.h"
#include "HBXFEMDefStruct.h"

extern	HBXDef::PrecsType_t	g_PrecsType;	//ȫ�ֱ������������ݶ�ȡ�ͼ���ľ���(float or double)

namespace HBXFEMDef
{
	//��ܵ�Ԫ�ṹ��
	template <class T>
	class _FrameElement:
		public _BaseElem<T>
	{
	protected:
		short iSection;	//��Ԫ����������
		T dLength;	//��Ԫ����
		T dSin, dCos;	//��Ԫ�ֲ�����x������������x��ļнǵ�������
		T daEndInterForce[6];	//��Ԫ�˶�������

	public:
		void StiffMatCal(){};
		T VolumeCalc(){return 0;};
		void MassMatCal( MassType_t _MassType = MassType_t::LUMPED ){};
		void SetMaterial( std::pair< boost::shared_ptr<_Material<T>>, boost::shared_ptr<_BaseSection<T>> > _MatAndSec ){ };

	};

	typedef boost::shared_ptr<_FrameElement<float>> _FrameElement_f_ptr;
	typedef boost::shared_ptr<_FrameElement<double>> _FrameElement_d_ptr;
}