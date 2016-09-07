#ifndef COMMON_H
#define COMMON_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>

#define USE_PETSC_LA

namespace LA
{
#ifdef USE_PETSC_LA
  using namespace dealii::LinearAlgebraPETSc;
#else
  using namespace dealii::LinearAlgebraTrilinos;
#endif
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_dgq.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <dirent.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <unistd.h>
/*
  Variables here are uniform within all problems.
*/

#define   Ndim 3

//e.g., Utilities::int_to_string(step=999,MAX_STEPS)-> 00999
#define MAX_STEPS 5
//e.g., Utilities::int_to_string(nprc=999,MAX_STEPS)-> 00999
#define MAX_PROCS 5

//#define HYPER_ELASTICITY
#define IN_ELASTICITY
//#define LINEAR_ELASTICITY
#define DPSI1

#define NMNDRCSPLIT 100

//define v solver

//IC
#define zeroU
#define FisX

//#define PETSC_SPARSE_DIRECT

//#define VSVR_PETSC_CG 
//#define ASVR_PETSC_CG 
//#define XSVR_PETSC_CG
#define FSVR_PETSC_CG

#define VSVR_PETSC_BICG_STAB
#define ASVR_PETSC_BICG_STAB
#define XSVR_PETSC_BICG_STAB
//#define FSVR_PETSC_BICG_STAB

#define ASSERT(condition, message) \
  do { \
  if (! (condition)) { \
  std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
  << " line " << __LINE__ << ": " << message << std::endl; \
  std::exit(EXIT_FAILURE); \
  } \
  } while (false)
				  
enum Commands {
  INITL,
  SOLVE,
  DSOLVE,
  CHECK,
  PRSRV,
  BCKWD,
  DPRSRV,
  DBCKWD,
  UPDATE,
  CLEAR,
  GETDT,
  UHIST,
  SETSF,
  SSCRV
} ;

// For a cube:
// reaction11(xmax)/A(xmax) vs strain_11 
// reaction22(ymax)/A(ymax) vs strain_22
// reaction33(zmax)/A(zmax) vs strain_33
// reaction11(ymax)/A(ymax) vs strain_12 
// reaction22(zmax)/A(zmax) vs strain_23
// reaction33(xmax)/A(xmax) vs strain_23
enum SScurveType {
  XX, 
  YY,
  ZZ,
  XY,
  YZ,
  XZ
};

#endif

