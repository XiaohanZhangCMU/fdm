#ifndef DISPL_EQUATION_H
#define DISPL_EQUATION_H

#include "fdm.h"

namespace FDM
{
  using namespace dealii;

  template<int dim> class HeatEquation;
  template<int dim> class PlasticEquation;

  template<int dim>
    class DisplEquation //: public Equation<dim> 
    {
    public:
      DisplEquation(Fdm<dim>* fdm)
	:
      fdm(fdm),
	//mapping(2),
	fe(FE_Q<dim>(1),dim),
	dof_handler(fdm->triangulation)
	  {
	    fe_values = NULL;
	    fe_face_values = NULL;
	  }
  
      ~DisplEquation()
	{	  	 	  
	  dof_handler.clear();	  	  
	}

      void set_coupled(PlasticEquation<dim-1>* peqn)
	{
	  this->peqn = peqn;
	}

      
    public:
      
      Fdm<dim>* fdm;
      PlasticEquation<dim-1>* peqn;
      HeatEquation<dim>* heqn;
      
      double clock;
      double dclock;
      int tstep;
      int static_step;
      std::vector<unsigned int> fix_rigid_dofs;
      std::vector<unsigned int> fix_translation_dofs;

      //const MappingQ<dim> mapping;    
      FESystem<dim> fe;      
      DoFHandler<dim> dof_handler;
      /* QGauss<dim> quadrature_formula; */
      /* QGauss<dim-1> face_quadrature_formula; */

      FEValues<dim>* fe_values;
      FEFaceValues<dim>* fe_face_values;

      IndexSet locally_owned_dofs;
      IndexSet locally_relevant_dofs;
      ConstraintMatrix constraints;
      LA::MPI::SparseMatrix system_matrix;
      LA::MPI::Vector system_rhs;
      LA::MPI::Vector locally_owned_u;
      LA::MPI::Vector locally_owned_du;
      LA::MPI::Vector static_residual;
      LA::MPI::Vector local_residual;
                  
    public:
      
      /* 
	 distribute dofs to all sovlers once for all!
	 set up solution, rhs and system matrix for all solvers.
      */

      int drive(Commands cmd);     
      int setup_system();
      int static_solve();
      int solve();
      void reinit_fe_values();
      int newton_solve();
      int static_assemble_system();
      void assemble_preconditioner(){}                      
      void output_results();
      bool dof_belongs_to_fix_rigid(const Point<dim>& p,
				    unsigned int d);
      bool dof_belongs_to_fix_translation(const Point<dim>&, unsigned int);

      void clear();
      /* helper functions */
      int print_fe_mesh();
      void initialize();
      void constraintMatrix();
      void set_rigid_dof();

      void
	static_element(
		       typename DoFHandler<dim>::active_cell_iterator&,  
		       unsigned int,
		       unsigned int,
		       unsigned int,
		       std::vector<types::global_dof_index>,
		       FullMatrix<double>&, 
		       Vector<double>&);

      Tensor<1,dim> 
	neumann(const Point<dim>& P, const Point<dim>& normal);

    };
}
#endif
