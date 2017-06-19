#pragma once

//ģ�����������޷�ʵ�֣���Ϊ�麯����ģ�������ڱ���֮ǰ����֪���麯����Ĵ�С����ģ�����ڱ���ǰ�ǲ�֪���麯����Ĵ�С�ġ����޳���
#include "StdAfx.h"
#include "HbxGloFunc.h"

extern int g_MatRowNum;

template < typename _PrescType >
class CBaseXConjugate
{
	typedef CBaseXConjugate<T>	_SameClass;
//ģ�����Ȱѱ���д��ǰ�棬tips...
private:

protected:
	typedef boost::shared_ptr< HBXDef::_CSRInput<_PrescType> > CSRInput_ptr;

	bool		m_bCudaFree;	//�ж��Ƿ����ͷ�CUDA��Դ
	bool		m_bGenTridiag;	//�ж������Ƿ��Ѿ�malloc�������malloc����hostmalloc�����ﵱ���ָ�벻Ϊ0ʱ����delete[]�����û����delete
	bool		m_bGPUFlag;		//�Ƿ�GPU����
	int			m_DevID;		//��ѡ�����õ�GPU��ID��
	bool		m_bSave;		//�Ƿ����
	boost::filesystem::path		m_SavePath;
	size_t		m_iters;	//��������

	cudaError_t			_cuError_id;		//cuda�ڴ濽���ȴ��󷵻�ֵ
	cublasStatus_t		_cublasError_id;		//cublas����ֵ
	cusparseStatus_t	_cusparseError_id;		//cusparse�÷���ֵ
	HBXDef::DataAloc_t		m_DataAloc;

	int		m_nnzA;	//CSR��ʽ��һ���͵ڶ�������ĳ���;
	int		m_nA;	//CSR��ʽ����������ĳ��ȣ����ھ���ĳһά�ȳ���+1;
	int		m_RowNum;	//�����������������Ϊ�����ĳһά�ĳ���

	//�����ݶȷ�����ر���
	cublasHandle_t		cublasHandle;	//CUBLAS����
	cusparseHandle_t	cusparseHandle;	//CUSPARSE����
	cusparseMatDescr_t	Matdescr;		//A�����������
	cusparseIndexBase_t	m_CSRIndexBase;	//�����������ʼλ��Ϊ0����1

	_PrescType	msecUsed;//ʱ��ͳ��double�ͱ���
	_PrescType	m_qaerr1;		//�в�

	//������buffer
	int			*h_iColSort, *h_iNonZeroRowSort;
	_PrescType	*h_NoneZeroVal;
	_PrescType	*h_x;
	_PrescType	*h_rhs;
	//�Դ�����ָ��
	_PrescType*	d_NonZeroVal;
	int*		d_iColSort;
	int*		d_iNonZeroRowSort;
	_PrescType*	d_x;	//��������ֵ
	_PrescType*	d_r;		//��ʽ�ұ�����

	std::vector<_PrescType>	h_vNoneZeroVal;
	std::vector<int>		h_viColSort;
	std::vector<int>		h_viNonZeroRowSort;
	std::vector<_PrescType>	h_vRhs;


public:
	CBaseXConjugate( size_t _RowNum = g_MatRowNum )
	{
		m_nA = (int)m_RowNum+1;
		m_bGenTridiag = false;
		m_bCudaFree = true;
		m_bGPUFlag = false;
		m_DevID = -1;
		m_DataAloc = HBXDef::INVALID;
		m_bSave = false;				//�Ƿ����

		msecUsed = (_PrescType)0.0f;
		m_iters = (_PrescType)0.0f;

		h_NoneZeroVal = nullptr;
		h_iColSort = nullptr;
		h_iNonZeroRowSort = nullptr;
		h_x = nullptr;
		h_rhs = nullptr;
	}
	virtual ~CBaseXConjugate()
	{
		if (nullptr == cusparseHandle)
		{
			return;
		}
		else FreeGPUResource();
		//�ͷ�CPU��Դ��todo something
	}

	virtual void	HostMalloc()
	{

	}
	virtual void	DeviceMalloc(){};	




};

