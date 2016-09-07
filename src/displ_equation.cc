#include "displ_equation.h"

//#define plastic_Debug

namespace FDM
{
using namespace dealii;

template<int dim>
int DisplEquation<dim>::drive(Commands cmd) {

  this->dclock = fdm->dclock;
  reinit_fe_values();

  if (cmd == INITL) {
    initialize();
    print_fe_mesh();
  }
  else if (cmd == SOLVE) {
    newton_solve();
  }
  else {
    std::cerr<<"fsvr::not implemented command"<<std::endl;
    exit(1);
  }
  return 0;
}

template<int dim>
void DisplEquation<dim>::clear()
{
  if (fe_values != NULL)
    delete fe_values;
  if (fe_face_values != NULL)
    delete fe_face_values;

}

template<int dim>
void DisplEquation<dim>::reinit_fe_values()
{
  clear();

  fe_values = new
  FEValues<dim> (fe, fdm->quadrature_formula,
                 update_values    |  update_gradients |
                 update_quadrature_points |
                 update_JxW_values);

  fe_face_values = new
  FEFaceValues<dim> (fe, fdm->face_quadrature_formula,
                     update_values   |
                     update_quadrature_points |
                     update_normal_vectors |
                     update_JxW_values);
}

template<int dim>
int DisplEquation<dim>::newton_solve() {
  double norm = 1e50;
  double Eps = 1e-7;
  int static_step = 0;

  fdm->pcout << "f static solve starts.....";
  while (norm > Eps) {
    static_assemble_system();
    norm = static_residual.linfty_norm();
#if defined plastic_Debug
    fdm->pcout << "static_step = " <<
                static_step << ", norm = " << norm << std::endl;
#endif
    static_solve();
    if (static_step == 150) {
      fdm->pcout << "Newton f static does not converge." << std::endl;
      return -1;
    }
    static_step++;
  }

  fdm->pcout << "f static solve ends with " <<
               static_step << " iters"<< std::endl;
  return 0;
}

template<int dim>
int DisplEquation<dim>::setup_system() {

  TimerOutput::Scope t(fdm->computing_timer, "plastic setup_system");
  fdm->pcout << "f setup starts......";
  dof_handler.distribute_dofs(fe);

  locally_owned_dofs = dof_handler.locally_owned_dofs();

  DoFTools::extract_locally_relevant_dofs(dof_handler,
                                          locally_relevant_dofs);
  fdm->U.reinit(locally_owned_dofs,
                locally_relevant_dofs,
                fdm->mpi_communicator);

  static_residual.reinit(locally_owned_dofs,
                         locally_relevant_dofs,
                         fdm->mpi_communicator);

  system_rhs.reinit(locally_owned_dofs,
                    fdm->mpi_communicator);

  system_rhs = 0;

  constraintMatrix();

  fdm->pcout << "ends" << std::endl;

  return 0;
}

template<int dim>
void DisplEquation<dim>::constraintMatrix()
{
  constraints.clear();
  constraints.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler,
      constraints);

  for (unsigned int dof_index = 0;
       dof_index < fix_rigid_dofs.size(); dof_index++)
    constraints.add_line(fix_rigid_dofs[dof_index]);

  constraints.close();

  DynamicSparsityPattern csp (locally_relevant_dofs);

  DoFTools::make_sparsity_pattern (dof_handler, csp,
                                   constraints, false);
  SparsityTools::distribute_sparsity_pattern (csp,
      dof_handler.n_locally_owned_dofs_per_processor(),
      fdm->mpi_communicator,
      locally_relevant_dofs);
  system_matrix.reinit (locally_owned_dofs,
                        locally_owned_dofs,
                        csp,
                        fdm->mpi_communicator);
}

template<int dim>
int DisplEquation<dim>::solve()
{

  TimerOutput::Scope t(fdm->computing_timer, "plastic solve");
  fdm->pcout<<"f solve starts......";

  SolverControl solver_control (dof_handler.n_dofs(), 1e-10);

#ifdef FSVR_PETSC_CG
  LA::SolverCG solver(solver_control, fdm->mpi_communicator);
  LA::MPI::PreconditionAMG preconditioner;

  LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
  data.symmetric_operator = true;
#else
  //trilinos defaults are good
#endif
  preconditioner.initialize(system_matrix, data);

  solver.solve (system_matrix,
                locally_owned_du,
                system_rhs,
                preconditioner);
  //if(tstep %nprint == 0)
  // fdm->pcout << "   step: " << fdm->tstep <<
  //            "\tplastic Solved in " << solver_control.last_step()
  //            << " iterations." << std::endl;

#elif defined FSVR_PETSC_SPARSE_DIRECT

  dealii::PETScWrappers::SparseDirectMUMPS
  solver(solver_control, fdm->mpi_communicator);

  solver.solve (system_matrix,
                locally_owned_du,
                system_rhs);
  // if(tstep % nprint == 0)
  //   fdm->pcout<<"   Step: "<<tstep<<
  //  "\tfsvr Solved with SparseDirectMumps "<<std::endl;
#elif defined FSVR_PETSC_BICG_STAB
  dealii::PETScWrappers::PreconditionNone preconditioner;
  preconditioner.initialize(system_matrix);
  dealii::PETScWrappers::SolverBicgstab
  solver(solver_control, fdm->mpi_communicator);

  solver.solve (system_matrix,
                locally_owned_du,
                system_rhs,
                preconditioner);

  //if(tstep % nprint == 0)
  //fdm->pcout << "   Step: " << fdm->tstep <<
  //           "\tplastic Solved with PETSc_Bicgstab" << std::endl;
#else

  throw ("No Solver Specified Error");

#endif

  constraints.distribute (locally_owned_du);

  //locally_owned_udot =locally_owned_du;
  // locally_owned_udot *= (1.0/dclock);
  /*
     follow two lines are only for debug
     turn them off to do real simulation
  */
  locally_owned_du *= dclock;
  locally_owned_u += locally_owned_du;

  fdm->U = locally_owned_u;
  fdm->pcout<<"ends"<<std::endl;
  return 0;
}

template<int dim>
int DisplEquation<dim>::static_assemble_system()
{
  TimerOutput::Scope t(fdm->computing_timer, "plastic static assembly");
  //fdm->pcout << "f static assemble starts......";

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = fdm->quadrature_formula.size();
  const unsigned int   n_face_q_points = fdm->face_quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  constraintMatrix();
  system_rhs = 0;
  system_matrix = 0;

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for (; cell != endc; ++cell)
    if (cell->is_locally_owned())
    {
      cell_matrix = 0;
      cell_rhs = 0;

      //std::cerr<<"start to assemble f static element..."<<fe++<<std::endl;

      try {
        cell->get_dof_indices (local_dof_indices);
        static_element(cell,
                       dofs_per_cell,
                       n_q_points, n_face_q_points,
		       local_dof_indices,
                       cell_matrix, cell_rhs);

        // for(int i = 0;i<dofs_per_cell;i++)
        //   std::cerr<<"cell_rhs("<<i<<") = "<<cell_rhs(i)
        //         <<std::endl;
      }
      catch (...)
      {
        std::cerr << "Exception on static 1" << std::endl
                  << "Aborting!" << std::endl
                  << std::endl;
      }
      /*
        assemble routine, see step-40.cc
      */

      try {

        //cell->distribute_local_to_global(cell_rhs,
        //local_residual);
        constraints.
        distribute_local_to_global (cell_matrix,
                                    cell_rhs,
                                    local_dof_indices,
                                    system_matrix,
                                    system_rhs);
      }
      catch (...)
      {
        std::cerr << "Exception on static 2" << std::endl

                  << "Aborting!" << std::endl
                  << std::endl;
      }
    }
  system_matrix.compress (VectorOperation::add);
  system_rhs.compress (VectorOperation::add);
  //local_residual.compress(VectorOperation::add);

#ifdef Debug
  const std::string filename = ("f" + Utilities::int_to_string(fdm->tstep, MAX_STEPS));
  std::ofstream out(filename.c_str());
  system_rhs.print(out, true, false);
#endif

  //static_residual = local_residual;
  static_residual = system_rhs;

  //fdm->pcout << "ends" << std::endl;

  return 0;
}


template<int dim>
int DisplEquation<dim>::static_solve()
{
  //fdm->pcout<<"plastic static solve starts......";
  TimerOutput::Scope t(fdm->computing_timer, "plastic static_solve");

  SolverControl solver_control (dof_handler.n_dofs(), 1e-10);

#ifdef FSVR_PETSC_CG
  LA::SolverCG solver(solver_control, fdm->mpi_communicator);
  LA::MPI::PreconditionAMG preconditioner;

  LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
  data.symmetric_operator = true;
#else
  //trilinos defaults are good
#endif
  preconditioner.initialize(system_matrix, data);

  solver.solve (system_matrix,
                locally_owned_du,
                system_rhs,
                preconditioner);
  // if(tstep %nprint == 0)
  // fdm->pcout << "   step: " << tstep <<
  //            "\tplastic static Solved in " << solver_control.last_step()
  //            << " iterations." << std::endl;

#elif defined FSVR_PETSC_SPARSE_DIRECT

  dealii::PETScWrappers::SparseDirectMUMPS
  solver(solver_control, fdm->mpi_communicator);

  solver.solve (system_matrix,
                locally_owned_du,
                system_rhs);
  // if(tstep % nprint == 0)
  //   fdm->pcout<<"   Step: "<<tstep<<
  //  "\tplastic static Solved with SparseDirectMumps "<<std::endl;
#elif defined FSVR_PETSC_BICG_STAB
  dealii::PETScWrappers::PreconditionNone preconditioner;
  preconditioner.initialize(system_matrix);
  dealii::PETScWrappers::SolverBicgstab
  solver(solver_control, fdm->mpi_communicator);

  solver.solve (system_matrix,
                locally_owned_du,
                system_rhs,
                preconditioner);

  // if(tstep % nprint == 0)
  // fdm->pcout<<"   Step: "<<tstep<<
  //   "\tplastic static Solved with PETSc_Bicgstab"<<std::endl;
#else

  throw ("No Solver Specified Error");

#endif

  constraints.distribute (locally_owned_du);

  locally_owned_u += locally_owned_du;

  //constraints.distribute (locally_owned_u);

  fdm->U = locally_owned_u;

  // locally_owned_udot = locally_owned_du;
  // locally_owned_udot *= (1.0/dclock);
  //fdm->Udot = locally_owned_udot;

  //fdm->pcout<<"plastic static solve ends"<<std::endl;
  //    static_residual = system_rhs;
  return 0;
}


//return the max plastic strain rate of all Gaussian points
template <int dim>
void DisplEquation<dim>::static_element(typename DoFHandler<dim>::active_cell_iterator& cell,
					  unsigned int dofs_per_cell,  
					  unsigned int n_q_points,
					  unsigned int n_face_q_points,
					  std::vector<types::global_dof_index> dof_indices,
					  FullMatrix<double>& cell_matrix,
					  Vector<double>& cell_rhs)
{
  assert(dof_indices == dof_indices);
  unsigned int i, j, q_point, face, 
    comp_i, comp_j;
  unsigned int mid = cell->material_id();
  //shape_grad = dphi^0/dx_i, i = 1,2,3 ;
  //shape = phi^0, phi^1, phi^2. 
  Tensor<1,dim> traction;
    
  Tensor<2,dim> eX, temp, dest;

  // std::cout<<"lambda = "<<this->lambda<<std::endl;
  // std::cout<<"mu = "<<this->mu<<std::endl;
  this->fe_values->reinit (cell);

  for(q_point = 0;q_point<n_q_points;q_point++)
    {    	
      Point<dim> xx = fe_values->quadrature_point(q_point);
	       
      //fdm->sLeX(xx,mid,eX,dest);

      for(i = 0;i<dofs_per_cell;i++)
	{
	  comp_i = this->fe.system_to_component_index(i).first;
	  cell_rhs(i)+= dest[comp_i]
	    *this->fe_values->shape_grad(i,q_point)
	    *this->fe_values->JxW(q_point);
	
	  for(j = 0;j<dofs_per_cell;j++)
	    {
	      comp_j = this->fe.system_to_component_index(j).first;
	      cell_matrix(i,j)
		+=
		(
		 (this->fe_values->shape_grad(i,q_point)[comp_i] *
		  this->fe_values->shape_grad(j,q_point)[comp_j] *
		  this->fdm->get_lambda(xx,mid))
		 +
		 (this->fe_values->shape_grad(i,q_point)[comp_j] *
		  this->fe_values->shape_grad(j,q_point)[comp_i] *
		  this->fdm->get_mu(xx,mid))
		 +
		 ((comp_i == comp_j) ?
		  (this->fe_values->shape_grad(i,q_point) *
		   this->fe_values->shape_grad(j,q_point) *
		   this->fdm->get_mu(xx,mid))  :
		  0)
		 )
		*
		this->fe_values->JxW(q_point);
	    }
	}
    }

  for( face = 0;face< GeometryInfo<dim>::faces_per_cell;++face) {
    //unsigned int boundary_indicator = cell->face(face)->boundary_indicator();
    unsigned int boundary_indicator = cell->face(face)->boundary_id();
    if((cell->face(face)->at_boundary())&&
       boundary_indicator > NMNDRCSPLIT )
      {
	this->fe_face_values->reinit(cell,face);
	   
	//std::cout<<"traction = "<<traction<<std::endl;
	for( i = 0;i<dofs_per_cell;i++)		 
	  {
	    comp_i = this->fe.system_to_component_index(i).first;
	    for( q_point = 0;q_point<n_face_q_points;++q_point)
	      {

		traction = 
		  (fdm->neumannFunction).vector_value(this->fe_face_values->quadrature_point(q_point), this->fe_face_values->normal_vector(q_point), boundary_indicator);
	   
		cell_rhs(i) += 
		  traction[comp_i] * 
		  this->fe_face_values->shape_value(i,q_point) * 
		  this->fe_face_values->JxW(q_point); 	       
	      }
	  }
      }     
  }
  
}//element


template <int dim>
void DisplEquation<dim>::output_results ()
{
  //fdm->pcout << "plastic output starts......";
  //MPI_Barrier(this->mpi_communicator);
  //assert(dim == 3);
  reinit_fe_values();
  DataOut<dim> data_out;
  data_out.attach_dof_handler (this->dof_handler);
  data_out.add_data_vector (fdm->U, "f");
  
  double local_sum_stress11,  local_sum_stress12,  local_sum_stress13, 
    local_sum_stress22,  local_sum_stress23,  local_sum_stress33;
  local_sum_stress11 = local_sum_stress12 = local_sum_stress13
    = local_sum_stress22 = local_sum_stress23 = local_sum_stress33 = 0;
  
  // data_out.add_data_vector (this->vsvr_total_displacement, "u");

  /* cauchy  stress -> vtk. the averaged stress of each cell */

  Vector<double> stress11(fdm->triangulation.n_active_cells());
  Vector<double> stress12(fdm->triangulation.n_active_cells());
  Vector<double> stress13(fdm->triangulation.n_active_cells());

  Vector<double> stress22(fdm->triangulation.n_active_cells());
  Vector<double> stress23(fdm->triangulation.n_active_cells());
  Vector<double> stress33(fdm->triangulation.n_active_cells());


  // Vector<double> temperature(fdm->triangulation.n_active_cells());

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points  = fdm->quadrature_formula.size();

  std::vector<std::vector< Tensor< 1, dim > > >
  dfdx(n_q_points, std::vector<Tensor<1, dim> >(dim));

  std::vector<Vector<double> >
  eX_values(n_q_points, Vector<double>(dim * dim));

  Vector<double>   f_values(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = this->dof_handler.begin_active(),
    endc = this->dof_handler.end();

  for (unsigned int index = 0; cell != endc; ++cell, ++index)
    if (cell->is_locally_owned())
    {

      this->fe_values->reinit (cell);     
      this->fe_values->get_function_gradients(fdm->U, dfdx);
      cell->get_dof_values(fdm->U, f_values);

      Tensor<2, dim> accumulated_stress, stress_of_this_qpt, strain;
      // double  temperature_of_this_qpt, accumulated_temperature= 0;
      Tensor<2, dim> G, F, eX, temp;
      Tensor<4, dim> d2psi;
      //unsigned int mid = cell->material_id();
      PointHistory<dim>* ph
        = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

      for (unsigned int q = 0;
           q < fdm->quadrature_formula.size(); ++q)
      {
        stress_of_this_qpt = ph[q].dpsi;

        accumulated_stress += stress_of_this_qpt;
        // accumulated_temperature += temperature_of_this_qpt;
      }

      stress11(index) = accumulated_stress[0][0] /
                        (fdm->quadrature_formula.size());
      stress12(index) = accumulated_stress[0][1] /
                        (fdm->quadrature_formula.size());
      stress13(index) = accumulated_stress[0][2] /
                        (fdm->quadrature_formula.size());
      stress22(index) = accumulated_stress[1][1] /
                        (fdm->quadrature_formula.size());
      stress23(index) = accumulated_stress[1][2] /
                        (fdm->quadrature_formula.size());
      stress33(index) = accumulated_stress[2][2] /
                        (fdm->quadrature_formula.size());


      local_sum_stress11 += stress11(index);
      local_sum_stress12 += stress12(index);
      local_sum_stress13 += stress13(index);
      local_sum_stress22 += stress22(index);
      local_sum_stress23 += stress23(index);
      local_sum_stress33 += stress33(index);
      
    }

  /*
    say the mesh is 3x3x5 = 45, then
    stressII(index) is a vector of length = fdm->n_active_cells = 45
    local_sum_stressII sums over all cells locally owned by this process
    add all local_sum_stressII across all processors then divide by 45
   */

  local_sum_stress11 =  Utilities::MPI::sum(local_sum_stress11,fdm->mpi_communicator);
  local_sum_stress12 =  Utilities::MPI::sum(local_sum_stress12,fdm->mpi_communicator);
  local_sum_stress13 =  Utilities::MPI::sum(local_sum_stress13,fdm->mpi_communicator);
  local_sum_stress22 =  Utilities::MPI::sum(local_sum_stress22,fdm->mpi_communicator);
  local_sum_stress23 =  Utilities::MPI::sum(local_sum_stress23,fdm->mpi_communicator);
  local_sum_stress33 =  Utilities::MPI::sum(local_sum_stress33,fdm->mpi_communicator);

  local_sum_stress11 /= fdm->triangulation.n_active_cells(); local_sum_stress12 /= fdm->triangulation.n_active_cells();
  local_sum_stress13 /= fdm->triangulation.n_active_cells(); local_sum_stress22 /= fdm->triangulation.n_active_cells();
  local_sum_stress23 /= fdm->triangulation.n_active_cells(); local_sum_stress33 /= fdm->triangulation.n_active_cells();
  
  data_out.add_data_vector(stress11, "stress11");
  data_out.add_data_vector(stress12, "stress12");
  data_out.add_data_vector(stress13, "stress13");

  data_out.add_data_vector(stress22, "stress22");
  data_out.add_data_vector(stress23, "stress23");
  data_out.add_data_vector(stress33, "stress33");

  // data_out.add_data_vector(temperature,"Theta");

  Vector<float> subdomain (fdm->triangulation.n_active_cells());
  //std::cout<<"subdomain.size = "<<subdomain.size()<<std::endl;
  for (unsigned int i = 0; i < subdomain.size(); ++i)
  {
    subdomain(i) = fdm->triangulation.locally_owned_subdomain();
    //std::cout<<"subdomain("<<i<<") ="<< subdomain(i)<<std::endl;
  }

  data_out.add_data_vector (subdomain, "subdomain");

  data_out.build_patches ();

  const std::string filename = ("fsvr_solution-" +
                                Utilities::int_to_string (fdm->tstep, MAX_STEPS) +
                                "." +
                                Utilities::int_to_string
                                (fdm->triangulation.locally_owned_subdomain(), MAX_PROCS));
  std::ofstream output ((filename + ".vtu").c_str());
  data_out.write_vtu (output);


  if (Utilities::MPI::this_mpi_process(fdm->mpi_communicator) == 0)
  {
    std::vector<std::string> filenames;
    for (unsigned int i = 0;
         i < Utilities::MPI::n_mpi_processes(fdm->mpi_communicator);
         ++i)
      filenames.push_back ("fsvr_solution-" +
                           Utilities::int_to_string (fdm->tstep, MAX_STEPS) +
                           "." +
                           Utilities::int_to_string (i, MAX_PROCS) +
                           ".vtu");

    std::ofstream master_output ((filename + ".pvtu").c_str());
    data_out.write_pvtu_record (master_output, filenames);
  }

  //fdm->pcout << "ends" << std::endl;
}

template<int dim>
void DisplEquation<dim>::set_rigid_dof()
{
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::map<types::global_dof_index, bool> dof_touched;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  for (; cell != endc; ++cell)
    if (cell->is_locally_owned())
    {
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; i++)
        dof_touched[local_dof_indices[i]] = false;
    }

  for (cell = dof_handler.begin_active(); cell != endc; ++cell)
    if (cell->is_locally_owned())
    {
      cell_rhs = 0;

      cell->get_dof_indices(local_dof_indices);

      for (unsigned int i = 0; i < dofs_per_cell; i++) {

        if (dof_touched[local_dof_indices[i]] == false)
        {
          // std::cerr<<"local_dof_indices["<<i<<"] = "<<local_dof_indices[i] <<"is false"<<std::endl;
          dof_touched[local_dof_indices[i]] = true;

          for (unsigned int v = 0;
               v < GeometryInfo<dim>::vertices_per_cell; ++v)
          {
            for (unsigned int d = 0; d < dim; d++)
            {
              if (local_dof_indices[i] ==
                  cell->vertex_dof_index(v, d))
              {
                Point<dim> the_vertex = cell->vertex(v);

                if (dof_belongs_to_fix_rigid(the_vertex, d))
                  fix_rigid_dofs.push_back(local_dof_indices[i]);

                if (dof_belongs_to_fix_translation(the_vertex, d))
                  fix_translation_dofs.push_back(local_dof_indices[i]);
              }
            }
          }
        }
      }
    }

  fdm->pcout << "Domain 0: f fix_rigid_dofs = ";
  for (unsigned int i = 0; i < fix_rigid_dofs.size(); i++)
    fdm->pcout << fix_rigid_dofs[i] << " ";
  fdm->pcout << "......"<<std::endl;
}


template<int dim>
bool DisplEquation<dim>:: dof_belongs_to_fix_rigid(const Point<dim>& p,
    unsigned int d)
{
  /*
     set fix_free_traction_dofs
     3d: Fix x,y,z degree for [xmax,ymin,zmin]
     Fix x,z degree for [xmax,ymax,zmin]
     Fix z degree for [xmin,ymax,zmin]
     2d: Fix x,y degree for [xmin,ymin]
     Fix y degree for [xmax,ymin]
  */

  if (dim == 2)
  {
    double xmin, xmax, ymin; //ymax;
    xmin = fdm->limits[0][0];
    ymin = fdm->limits[1][0];
    xmax = fdm->limits[0][1];

    if (fdm->mesher == 1 || fdm->mesher == 2 || fdm->mesher == 3) {
      // std::cerr<<"xmin = "<<xmin<<std::endl;
      // std::cerr<<"ymin = "<<ymin<<std::endl;
      // std::cerr<<"xmax = "<<xmax<<std::endl;

      //ymax = limits[1][1];

      if ((fabs(p[0] - xmin) < 1e-12 &&
           fabs(p[1] - ymin) < 1e-12 )

          || (fabs(p[0] - xmax) < 1e-12 &&
              fabs(p[1] - ymin) < 1e-12 &&
              d == 1))
      {
        // std::cerr<<"xmin = "<<xmin<<std::endl;
        // std::cerr<<"ymin = "<<ymin<<std::endl;
        // std::cerr<<"xmax = "<<xmax<<std::endl;
        // std::cerr<<"p[0] = "<<p[0]<<std::endl;
        // std::cerr<<"p[1] = "<<p[1]<<std::endl;
        // std::cerr<<"(fabs(p[0]-xmin) = "<<(fabs(p[0]-xmin))<<std::endl;
        // std::cerr<<"(fabs(p[1]-ymin) = "<<(fabs(p[1]-ymin))<<std::endl;
        // std::cerr<<"(fabs(p[0]-xmax) = "<<(fabs(p[0]-xmax))<<std::endl;
        // std::cerr<<"p = "<<p<<std::endl;
        return true;
      }
    }
  }
  else if (dim == 3)
  {
    double xmin, xmax, ymin, zmin, zmax, ymax;
    xmin = fdm->limits[0][0];
    ymin = fdm->limits[1][0];
    zmin = fdm->limits[2][0];
    xmax = fdm->limits[0][1];
    ymax = fdm->limits[1][1];
    zmax = fdm->limits[2][1];

    if (fdm->mesher == 1 || fdm->mesher == 2 || fdm->mesher == 3 || fdm->mesher == 4) {
      if ((fabs(p[0] - xmax) < 1e-12 &&
           fabs(p[1] - ymin) < 1e-12 &&
           fabs(p[2] - zmin) < 1e-12)
          // ||(fabs(p[0]-xmax)<1e-12 &&
          //   fabs(p[1]-ymax)<1e-12 &&
          //   fabs(p[2]-zmin)<1e-12 &&
          //   d != 1)
          || (fabs(p[0] - xmax) < 1e-12 &&
              fabs(p[1] - ymin) < 1e-12 &&
              fabs(p[2] - zmax) < 1e-12 &&
              d != 2)

          // ||(fabs(p[0]-xmin)<1e-12 &&
          //   fabs(p[1]-ymax)<1e-12 &&
          //   fabs(p[2]-zmin)<1e-12 &&
          //   d == 2))
          || (fabs(p[0] - xmin) < 1e-12 &&
              fabs(p[1] - ymin) < 1e-12 &&
              fabs(p[2] - zmax) < 1e-12 &&
              d == 1))

      {
        return true;
      }
    }
    else if (fdm->mesher == 4)
    {
      Point<dim> p0(xmin, 0, zmin);
      Point<dim> p1(xmin, ymax, 0);
      Point<dim> p2(xmin, ymin, 0);

      if (p.distance(p0) < 1e-12 ||
          (p.distance(p1) < 1e-12 && d != 1) ||
          (p.distance(p2) < 1e-12 && d == 2))
        return true;
      // if((fabs(p[0]-xmin)<1e-12 &&
      //  fabs(p[1]-ymin)<1e-12 &&
      //  (p[2])<0)
      //    ||(fabs(p[0]-xmin)<1e-12 &&
      //    fabs(p[1]-ymin)<1e-12 &&
      //    (p[2])>0 &&
      //    d != 1)
      //    ||(fabs(p[0]-xmin)<1e-12 &&
      //    fabs(p[1]-ymax)<1e-12 &&
      //    (p[2])<0 &&
      //    d == 2))
      //   {
      //  return true;
      //   }

    }
    //zmax = limits[2][1];

    // std::cerr<<"xmin = "<<xmin<<std::endl;
    // std::cerr<<"ymin = "<<ymin<<std::endl;
    // std::cerr<<"zmin = "<<zmin<<std::endl;
    // std::cerr<<"xmax = "<<xmax<<std::endl;
    // std::cerr<<"ymax = "<<ymax<<std::endl;

  }

  return false;
}


template<int dim>
bool DisplEquation<dim>:: dof_belongs_to_fix_translation(const Point<dim>& p,
    unsigned int d)
{
  assert(d<Ndim);
  /*
     set fix_translation_dofs
     3d: Fix x,y,z degree for [xmax,ymin,zmin]
     2d: Fix x,y degree for [xmin,ymin]
  */

  if (dim == 2)
  {
    double xmin, ymin; //ymax;
    xmin = fdm->limits[0][0];
    ymin = fdm->limits[1][0];

    // std::cerr<<"xmin = "<<xmin<<std::endl;
    // std::cerr<<"ymin = "<<ymin<<std::endl;
    // std::cerr<<"xmax = "<<xmax<<std::endl;

    //ymax = limits[1][1];

    if ((fabs(p[0] - xmin) < 1e-12 &&
         fabs(p[1] - ymin) < 1e-12 ))
      return true;
  }
  else if (dim == 3)
  {
    double xmax, ymin, zmin; //,zmax;

    ymin = fdm->limits[1][0];
    zmin = fdm->limits[2][0];
    xmax = fdm->limits[0][1];

    //zmax = limits[2][1];

    // std::cerr<<"xmin = "<<xmin<<std::endl;
    // std::cerr<<"ymin = "<<ymin<<std::endl;
    // std::cerr<<"zmin = "<<zmin<<std::endl;
    // std::cerr<<"xmax = "<<xmax<<std::endl;
    // std::cerr<<"ymax = "<<ymax<<std::endl;

    if ((fabs(p[0] - xmax) < 1e-12 &&
         fabs(p[1] - ymin) < 1e-12 &&
         fabs(p[2] - zmin) < 1e-12))
    {
      return true;
    }
  }

  return false;
}



template<int dim>
void DisplEquation<dim>::initialize()
{
  //fdm->pcout << "f initialize f......" << std::endl;
  VectorTools::interpolate(//this->mapping,
    dof_handler,
    ZeroFunction<dim>(dim),
    this->locally_owned_u);

  fdm->U = this->locally_owned_u;
  set_rigid_dof();
  //fdm->pcout << "ends" << std::endl;
}




template<int dim>
int DisplEquation<dim>::print_fe_mesh()
{
#define PrintFeMesh

#ifdef PrintFeMesh
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  //   std::vector<bool> vertex_touched (triangulation.n_vertices(),false);

  const std::string filename = "./Output/femesh." + Utilities::int_to_string(fdm->triangulation.locally_owned_subdomain(), 4);
  std::ofstream out(filename.c_str());

  unsigned int cell_id = 0;

  out << "Localy owned finite elements:" << std::endl;

  for (; cell != endc; ++cell)
  {
    /* vsvr */
    if (cell -> is_locally_owned())
    {
      out << "elem\t" << cell_id++ << std::endl;

      for (unsigned int v = 0;
           v < GeometryInfo<dim>::vertices_per_cell; ++v)
      {
        out << "Vertex\t" << cell->vertex_index(v) <<
            ":\t" << cell->vertex(v) << std::endl;
        out << "\tv_dof:\t";
        for (unsigned int d = 0; d < dim; ++d)
        {
          out << cell->vertex_dof_index(v, d) << "\t";
        }
        out << "\n\tx_dof:\t";
      }
    }
  }
  
  out << "Ghost elements on this process" << std::endl;

#endif
  //fdm->pcout << "fe mesh printed." << std::endl;
  return 0;
}


}
