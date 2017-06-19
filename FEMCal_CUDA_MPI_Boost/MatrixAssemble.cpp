#include "StdAfx.h"
#include "HbxGloFunc.h"
#include "MatrixAssemble.h"


MatrixAssemble::MatrixAssemble(void)
{
}


MatrixAssemble::~MatrixAssemble(void)
{
}

// CMatrixAssemble类函数
void MatrixAssemble::AssembleAnalysis( HBXDef::PrecsType_t _precs )
{
	if ( HBXDef::PrecsType_t::_PRECS_FLOAT_ == _precs )
	{
		m_SumStiffMat_f.ResetSize( g_Dim * g_vNode_f.size(), g_Dim * g_vNode_f.size() );
		m_SumMassMat_f.ResetSize( g_Dim * g_vNode_f.size(), g_Dim * g_vNode_f.size() );
	}
	else if ( HBXDef::PrecsType_t::_PRECS_DOUBLE_ == _precs )
	{
		m_SumStiffMat_d.ResetSize( g_Dim * g_vNode_d.size(), g_Dim * g_vNode_d.size() );
		m_SumMassMat_d.ResetSize( g_Dim * g_vNode_d.size(), g_Dim * g_vNode_d.size() );
	}
//	PreCondit( _precs );
	boost::timer::auto_cpu_timer	StiffCalTimer;
	StiffMatrixCal( _precs );
	std::cout<<StiffCalTimer.format(6)<<std::endl;
	MassMatrixCal( _precs );

}

//总体分析程序，单元刚度和质量矩阵计算的预处理函数
void	MatrixAssemble::PreCondit( HBXDef::PrecsType_t _precs )
{
	return;
}

void MatrixAssemble::StiffMatrixCal( HBXDef::PrecsType_t _precs )
{
	if ( HBXDef::PrecsType_t::_PRECS_FLOAT_ == _precs )
	{
		for (auto _iter=g_vBaseElem_f.begin(); _iter!=g_vBaseElem_f.end(); _iter++)
		{
			(*_iter)->StiffMatCal();
			//直接存储方法
			int i,j,k,ic,jc;
			int EleMatRowNum = ((*_iter)->_idim)*((*_iter)->_vNode.size());
			for (i=0; i<EleMatRowNum; i++)
			{
				ic = (*_iter)->_ipmark[i];
				for(j=0; j<EleMatRowNum; j++)
				{
					jc = (*_iter)->_ipmark[j];
					m_SumStiffMat_d(ic,jc) += (*_iter)->_Stiffmat(i,j);
				}
			}
		}
		return;
	}
	else if ( HBXDef::PrecsType_t::_PRECS_DOUBLE_ == _precs )
	{
		for (auto _iter=g_vBaseElem_d.begin(); _iter!=g_vBaseElem_d.end(); _iter++)
		{
			(*_iter)->StiffMatCal();
			//直接存储方法
			int i,j,k,ic,jc;
			int EleMatRowNum = ((*_iter)->_idim)*((*_iter)->_vNode.size());
			for (i=0; i<EleMatRowNum; i++)
			{
				ic = (*_iter)->_ipmark[i];
				for(j=0; j<EleMatRowNum; j++)
				{
					jc = (*_iter)->_ipmark[j];
					m_SumStiffMat_d(ic,jc) += (*_iter)->_Stiffmat(i,j);
				}
			}
		}
		return;
	}
	else
	{
		std::cerr<<"刚度矩阵计算函数传入无效精度标识"<<std::endl;
		return;
	}
}

void MatrixAssemble::MassMatrixCal( HBXDef::PrecsType_t _precs )
{
	using namespace HBXFEMDef;
	if ( HBXDef::PrecsType_t::_PRECS_FLOAT_ == _precs )
	{
		for (auto _iter=g_vBaseElem_f.begin(); _iter!=g_vBaseElem_f.end(); _iter++)
		{
			(*_iter)->MassMatCal( MassType_t::LUMPED );
			//直接存储方法
			int i,j,k,ic,jc;
			int EleMatRowNum = ((*_iter)->_idim)*((*_iter)->_vNode.size());
			for (i=0; i<EleMatRowNum; i++)
			{
				ic = (*_iter)->_ipmark[i];
				for(j=0; j<EleMatRowNum; j++)
				{
					jc = (*_iter)->_ipmark[j];
					m_SumMassMat_f(ic,jc) += (*_iter)->_Massmat(i,j);
				}
			}
		}
		return;
	}
	else if ( HBXDef::PrecsType_t::_PRECS_DOUBLE_ == _precs )
	{
		for (auto _iter=g_vBaseElem_d.begin(); _iter!=g_vBaseElem_d.end(); _iter++)
		{
			(*_iter)->MassMatCal( MassType_t::LUMPED );
			//直接存储方法
			int i,j,k,ic,jc;
			int EleMatRowNum = ((*_iter)->_idim)*((*_iter)->_vNode.size());
			for (i=0; i<EleMatRowNum; i++)
			{
				ic = (*_iter)->_ipmark[i];
				for(j=0; j<EleMatRowNum; j++)
				{
					jc = (*_iter)->_ipmark[j];
					m_SumMassMat_d(ic,jc) += (*_iter)->_Massmat(i,j);
				}
			}
		}
		return;
	}
	else
	{
		std::cerr<<"质量矩阵计算函数传入无效精度标识"<<std::endl;
		return;
	}
}

//组装单元刚度矩阵成为总装矩阵
void MatrixAssemble::AssembleStiffMat( HBXDef::PrecsType_t _precs )
{

}

//结果输出到文件
void	MatrixAssemble::ExportResultFile( boost::filesystem::path _path,  HBXDef::SaveType_t _save_t , HBXDef::PrecsType_t _precs )
{
	using namespace boost::gregorian;
	using namespace boost::posix_time;

	m_SavePath = _path;
	//		m_SavePath.append( to_iso_string(_pathtime) );
	if (exists(m_SavePath))
	{
		if ( boost::filesystem::is_empty(m_SavePath) )
		{
			remove(m_SavePath);
		}
	}
	else
	{
		remove_all(m_SavePath);
	}
	assert(!exists(m_SavePath));
	create_directories(m_SavePath);
	if ( HBXDef::PrecsType_t::_PRECS_DOUBLE_ == _precs )
	{
		HBXDef::_FullMatricToCSR( this->m_SumStiffMat_d, _vNoneZeroVal, _viColSort,_viNonZeroRowSort );
	}
	else if ( HBXDef::PrecsType_t::_PRECS_FLOAT_ == _precs )
	{
		HBXDef::_FullMatricToCSR( this->m_SumStiffMat_f, _vNoneZeroVal, _viColSort,_viNonZeroRowSort );
	}

	if ( HBXDef::SaveType_t::CSR == _save_t )
	{
		HBXDef::_ExportCSRdataToFile( _vNoneZeroVal.data(), _viColSort.data(), _viNonZeroRowSort.data(),
			_vNoneZeroVal.size(), _viNonZeroRowSort.size()-1, "NoneZeroVal.data", "ColSort.data", "NoneZeroRowSort.data",
			m_SavePath );
	}
	
}