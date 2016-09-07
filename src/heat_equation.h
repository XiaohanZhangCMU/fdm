#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H

#include "fdm.h"
//#include "equation.h"

#include "displ_equation.h"

namespace FDM
{
  using namespace dealii;

  template<int dim> class DisplEquation;

  template<int dim>
    class HeatEquation //: public Equation<dim> 
    {
    public:
      HeatEquation(Fdm<dim>* fdm)
	:
      fdm(fdm),
	//mapping(2),
	fe(FE_Q<dim>(1),1),
	dof_handler(fdm->triangulation)
	/* quadrature_formula(3), */
	/* face_quadrature_formula(3), */
	  {
	    fe_values = NULL;
	    fe_face_values = NULL;
	  }
  
      ~HeatEquation()
	{	  	 	  
	  dof_handler.clear();	  	  
	}

      void set_coupled(DisplEquation<dim>* heqn,
		       PlasticEquation<dim>* peqn)
	{
	  this->heqn = heqn;
	  this->peqn = peqn;
	}

      
    public:
      
      Fdm<dim>* fdm;
      PlasticEquation<dim>* peqn;
      DisplEquation<dim>* deqn;

      double clock;
      double dclock;

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

      //current solve
      LA::MPI::Vector locally_owned_theta;

      LA::MPI::Vector locally_owned_dtheta;
      
      //for restore
      LA::MPI::Vector locally_owned_p_theta;
      LA::MPI::Vector locally_owned_pp_theta;
                  
    public:
      
      /* 
	 distribute dofs to all sovlers once for all!
	 set up solution, rhs and system matrix for all solvers.
      */

      int drive(Commands cmd);      
      int setup_system();
      void preserve();
      void restore();
      void clear();
      int solve();
      void reinit_fe_values();
      int assemble_system();
      void assemble_preconditioner(){}                      
      void output_results();
      void set_rigid_dof();
      
      bool dof_belongs_to_fix_translation(const Point<dim>& p,
					  unsigned int d);

      /* helper functions */
      void initialize();
      void constraintMatrix(int);

      void 
	element(typename DoFHandler<dim>::active_cell_iterator&,  
		unsigned int,
		unsigned int,
		unsigned int,
		FullMatrix<double>&, 
		Vector<double>&);

    };
}
#endif
