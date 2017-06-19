#pragma once
#include "stdafx.h"
#include "HBXFEMDefStruct.h"


extern	HBXDef::PrecsType_t	g_PrecsType;	//ȫ�ֱ������������ݶ�ȡ�ͼ���ľ���(float or double)
extern	int		g_Dim;	//ȫ��ά�ȣ���������2��3ά��Ͻ�ģ

namespace HBXFEMDef
{
#pragma region ��Ԫ�����

	template<class T=double>	//��Ԫ����
	class  _BaseElem
	{
	protected:
		//short _iMaterial;	//��Ԫ����������
	public:
		//@_Num_:��ǰ��Ԫ���
		//@_type_:��Ԫ����
		//@_lhs���ڵ�������������
		//@_nodenum���ڵ�����
		_BaseElem( unsigned int _Num_, unsigned short _type_, std::vector<std::string>::iterator _lhs, size_t _nodenum, unsigned char _dim_ );

		virtual void SetMaterial( std::pair< boost::shared_ptr<_Material<T>>, boost::shared_ptr<_BaseSection<T>> > _MatAndSec ) = 0;
		virtual void SetStiffMat( HBXDef::CMatrix<T> &_rhsmat )
		{
			_Stiffmat = _rhsmat;
		}
		virtual void StiffMatCal() = 0;	//���麯������ɵ�Ԫ�նȾ���ļ���
		virtual T VolumeCalc() = 0;	//���麯������ɵ�Ԫ����ļ���
		virtual void MassMatCal( MassType_t _MassType = MassType_t::LUMPED ) = 0;	//���麯������ɵ�Ԫ��������ļ���

	public:
		unsigned int	_iNum;	//��ǰ��Ԫ���
		unsigned char	_idim;	//ÿ���ڵ����ɶ�
		unsigned short	_iType;		//Ӧ�����Ӧ������
		std::vector<unsigned int>	_vNode;//��Ԫ�ڸ��ڵ��
		HBXDef::CMatrix<T>	_Stiffmat;	//����������͸նȾ���ֱ����������ͻ����У�����Ϊ��û�м����ضȵ�Ҫ����������ɲ��ã����նȾ�������У���Ϊ��ʡ�ռ䣬�ֿ�
		HBXDef::CMatrix<T>	_Massmat;	//��������,�������ڻ����а�...
		std::vector<int>	_ipmark;	//��Ԫ�նȾ������ܸվ���֮���ӳ���,����Ԫ�ն����е�i��(��)���ܸ����е�λ��
		boost::shared_ptr<_Material<T>>	m_Matreial_ptr;
	};

#pragma  endregion ������Ԫ�����

}

#include "CBaseElemt.inl"