#pragma once

//模板类做基类无法实现，因为虚函数在模板类内在编译之前就需知道虚函数表的大小，而模板类在编译前是不知道虚函数表的大小的。故剔除。
#include "StdAfx.h"
#include "HbxGloFunc.h"

extern int g_MatRowNum;

template < typename _PrescType >
class CBaseXConjugate
{
	typedef CBaseXConjugate<T>	_SameClass;
//模板类先把变量写在前面，tips...
private:

protected:
	typedef boost::shared_ptr< HBXDef::_CSRInput<_PrescType> > CSRInput_ptr;

	bool		m_bCudaFree;	//判断是否在释放CUDA资源
	bool		m_bGenTridiag;	//判断数组是否已经malloc过，如果malloc过，hostmalloc函数里当检测指针不为0时，用delete[]；如果没有用delete
	bool		m_bGPUFlag;		//是否GPU可用
	int			m_DevID;		//所选择的最好的GPU的ID号
	bool		m_bSave;		//是否存盘
	boost::filesystem::path		m_SavePath;
	size_t		m_iters;	//迭代次数

	cudaError_t			_cuError_id;		//cuda内存拷贝等错误返回值
	cublasStatus_t		_cublasError_id;		//cublas返回值
	cusparseStatus_t	_cusparseError_id;		//cusparse得返回值
	HBXDef::DataAloc_t		m_DataAloc;

	int		m_nnzA;	//CSR格式第一个和第二个数组的长度;
	int		m_nA;	//CSR格式第三个数组的长度，等于矩阵某一维度长度+1;
	int		m_RowNum;	//矩阵的行数，方阵则为矩阵的某一维的长度

	//共轭梯度法的相关变量
	cublasHandle_t		cublasHandle;	//CUBLAS对象
	cusparseHandle_t	cusparseHandle;	//CUSPARSE对象
	cusparseMatDescr_t	Matdescr;		//A矩阵的描述器
	cusparseIndexBase_t	m_CSRIndexBase;	//矩阵的索引起始位，为0或者1

	_PrescType	msecUsed;//时间统计double型变量
	_PrescType	m_qaerr1;		//残差

	//主机端buffer
	int			*h_iColSort, *h_iNonZeroRowSort;
	_PrescType	*h_NoneZeroVal;
	_PrescType	*h_x;
	_PrescType	*h_rhs;
	//显存内裸指针
	_PrescType*	d_NonZeroVal;
	int*		d_iColSort;
	int*		d_iNonZeroRowSort;
	_PrescType*	d_x;	//待求特征值
	_PrescType*	d_r;		//等式右边向量

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
		m_bSave = false;				//是否存盘

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
		//释放CPU资源，todo something
	}

	virtual void	HostMalloc()
	{

	}
	virtual void	DeviceMalloc(){};	




};

