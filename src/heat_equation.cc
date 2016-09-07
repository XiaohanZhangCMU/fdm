#include "heat_equation.h"

#define Debug

/* err: 
   -1: newton solve does not converge
*/

namespace FDM
{
  using namespace dealii;

  template<int dim>
  int HeatEquation<dim>::drive(Commands cmd) {

    this->dclock = fdm->dclock;
    reinit_fe_values();

    if (cmd == INITL) {
      initialize();
    }
    else if (cmd == SOLVE){
      assemble_system();
      solve();
    }
    else if (cmd == BCKWD) {
      restore();
    } else if (cmd == CHECK) {

    } else if (cmd == CLEAR) {
      clear();          
    } else if (cmd == PRSRV) {
      preserve();
    } else 
      return -1;
    return 0;
  }

  template<int dim>
  void HeatEquation<dim>::clear()
  { 
    if (fe_values != NULL)
      delete fe_values;
    if (fe_face_values != NULL)
      delete fe_face_values;

  }

  template<int dim>
  void HeatEquation<dim>::reinit_fe_values()
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
			 update_normal_vectors|
			 update_JxW_values);
  }


  template<int dim>
  int HeatEquation<dim>::setup_system() {

    TimerOutput::Scope t(fdm->computing_timer, "theta setup_system");
    fdm->pcout<<"t setup starts......";
    dof_handler.distribute_dofs(fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    
    DoFTools::extract_locally_relevant_dofs(dof_handler,
					    locally_relevant_dofs);
    fdm->theta.reinit(locally_owned_dofs,
		      locally_relevant_dofs,
		      fdm->mpi_communicator);
    system_rhs.reinit(locally_owned_dofs,
		      fdm->mpi_communicator);

    system_rhs = 0;

    locally_owned_theta.reinit(locally_owned_dofs,
				fdm->mpi_communicator);
    locally_owned_dtheta.reinit(locally_owned_dofs,
				fdm->mpi_communicator);
    locally_owned_p_theta.reinit(locally_owned_dofs,
			     fdm->mpi_communicator);    
    locally_owned_pp_theta.reinit(locally_owned_dofs,
			      fdm->mpi_communicator);    

    constraintMatrix(0);

    fdm->pcout<<"ends"<<std::endl;

    return 0;
  }

  template<int dim>
  void HeatEquation<dim>::constraintMatrix(int caseid)
  {    
    assert(caseid>=0);
    constraints.clear();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, 
					     constraints);

    // for(unsigned int dof_index = 0;
    // 	dof_index<fix_translation_dofs.size();dof_index++)	
    //   constraints.add_line(fix_translation_dofs[dof_index]);
    
    constraints.close();

    //CompressedSimpleSparsityPattern csp (locally_relevant_dofs);
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
  void HeatEquation<dim>::preserve()
  {
    if (fdm->tstep == 0)
      {
	locally_owned_p_theta = locally_owned_theta;
	locally_owned_pp_theta = locally_owned_p_theta;
      }
    else 
      {
	locally_owned_pp_theta = locally_owned_p_theta;
	locally_owned_p_theta = locally_owned_theta;
      }
  }

  template<int dim>
  void HeatEquation<dim>::restore()
  {
    locally_owned_theta = locally_owned_pp_theta;
    locally_owned_p_theta = locally_owned_pp_theta;
    fdm->theta = locally_owned_theta;   
  }

  template<int dim>
  int HeatEquation<dim>::solve()
  {

    TimerOutput::Scope t(fdm->computing_timer, "plastic solve");
    fdm->pcout<<"t solve starts......";
    
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
		  locally_owned_dtheta, 
		  system_rhs,
		  preconditioner);
    //if(tstep %nprint == 0)
    // fdm->pcout << "   step: "<<fdm->tstep<<
    //   "\tplastic Solved in " << solver_control.last_step()
    // 	       << " iterations." << std::endl;

#elif defined FSVR_PETSC_SPARSE_DIRECT

    dealii::PETScWrappers::SparseDirectMUMPS 
      solver(solver_control, fdm->mpi_communicator);
    
    solver.solve (system_matrix, 
		  locally_owned_dtheta, 
		  system_rhs);
    // if(tstep % nprint == 0)
    //   fdm->pcout<<"   Step: "<<tstep<<
    // 	"\tfsvr Solved with SparseDirectMumps "<<std::endl;
#elif defined FSVR_PETSC_BICG_STAB
    dealii::PETScWrappers::PreconditionNone preconditioner;
    preconditioner.initialize(system_matrix); 
    dealii::PETScWrappers::SolverBicgstab 
      solver(solver_control, fdm->mpi_communicator);
   
    solver.solve (system_matrix, 
		  locally_owned_dtheta, 
		  system_rhs,
		  preconditioner);

    //if(tstep % nprint == 0)
    // fdm->pcout<<"   Step: "<<fdm->tstep<<
    //   "\ttheta Solved with PETSc_Bicgstab"<<std::endl;
#else

    throw("No Solver Specified Error");

#endif

    constraints.distribute (locally_owned_dtheta);    

    //locally_owned_fdot =locally_owned_dtheta;
    // locally_owned_fdot *= (1.0/dclock);
    /* 
       follow two lines are only for debug
       turn them off to do real simulation 
    */
    locally_owned_dtheta *= this-> dclock;
    locally_owned_theta += locally_owned_dtheta;

    fdm->theta = locally_owned_theta;    
    fdm->pcout<<"ends"<<std::endl;
    return 0;
  }


  template<int dim>
  int HeatEquation<dim>::assemble_system()
  {
    TimerOutput::Scope t(fdm->computing_timer, "plastic assembly");
    fdm->pcout<<"t assemble starts......";

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = fdm->quadrature_formula.size();
    const unsigned int   n_face_q_points = fdm->face_quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    constraintMatrix(1);

    system_rhs = 0;
    system_matrix = 0;

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
	{
	  cell_matrix = 0;
	  cell_rhs = 0;

	  //std::cerr<<"start to assemble velement"<<std::endl;
	  try{
	    element(cell,
		    dofs_per_cell, n_q_points, n_face_q_points,
		    cell_matrix, cell_rhs);	    

	  }
	  catch(...)
	    {	      
	      std::cerr<<"Exception on plastic 1"<<std::endl     
		       <<"Aborting!"<<std::endl
		       <<std::endl;
	      
	    }

	  try{
	    cell->get_dof_indices (local_dof_indices);
	  
	    constraints.
	      distribute_local_to_global (cell_matrix,
					  cell_rhs,
					  local_dof_indices,
					  system_matrix,
					  system_rhs);
	  }
	  catch(...)
	    {
	      std::cerr<<"Exception on plastic 2"<<std::endl
		       <<"Aborting!"<<std::endl
		       <<std::endl;
	    }
	}

    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);       
    fdm->pcout<<"ends......";
    return 0;
  }


  template <int dim>
  void HeatEquation<dim>::element(typename DoFHandler<dim>::active_cell_iterator& cell,  
				  unsigned int dofs_per_cell,
				  unsigned int n_q_points,
				  unsigned int n_face_q_points,
				  FullMatrix<double>& cell_matrix, 
				  Vector<double>& cell_rhs)
  {
    assert(n_face_q_points == n_face_q_points);
    /* declare the following for reusing */
    unsigned int i, j, d, m, q_point, 
      comp_i, comp_j;
    double  G;    
    Tensor<1,dim> dtheta, V;
    Tensor<2,dim> alpha, temp, Lp, dpsi;    
    unsigned int mid = cell->material_id();
    std::vector<Vector<double> > 
      theta_values(n_q_points, Vector<double>(1));
    std::vector<std::vector<Tensor<1,dim> > > 
      dtheta_values(n_q_points, std::vector<Tensor<1,dim> >(1));

    std::vector<Vector<double> > 
      alpha_values(n_q_points, Vector<double>(dim*dim));

    this->fe_values->reinit (cell);

    this->fe_values->get_function_values(fdm->theta,
					 theta_values);
    this->fe_values->get_function_gradients(fdm->theta,
					    dtheta_values);
    PointHistory<dim>* ph 
      =reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    /* quadrature loop */

    for(q_point = 0;q_point < n_q_points;q_point++)
      {
	Point<dim> xx = fe_values->quadrature_point(q_point);
	double k = fdm->get_ThetaDiffusion(xx,mid);
	double kai = fdm->get_ThetaWorkConvs(xx,mid);

	for(d = 0;d<dim;d++) {
	  dtheta[d] = dtheta_values[q_point][0][d];
	  for(m = 0;m<dim;m++)
	    alpha[d][m] = alpha_values[q_point][d*dim+m];	     
	}	

	for(i = 0;i<dofs_per_cell;i++)
	  {
	    comp_i = this->fe.
	      system_to_component_index(i).first;
	    assert(comp_i < dim*dim);
	    
	    dpsi = ph[q_point].dpsi;
	    Lp = ph[q_point].Lp;
	    V  = ph[q_point].V;
	    temp = fdm->AxV(alpha,V) + Lp;
	    G = double_contract(dpsi,temp);

	    cell_rhs(i)+=
	      (
	       -k*dtheta* (fe_values->shape_grad(i,q_point))
	       + kai* G *this->fe_values->shape_value(i,q_point))
	      *this->fe_values->JxW(q_point);
	  
	    for(j = 0;j<dofs_per_cell;j++)
	      {
		comp_j = this->fe.
		  system_to_component_index(j).first;

		assert(comp_j < dim*dim);

		if(comp_j == comp_i)
		  {		   
		    cell_matrix(i,j) += 
		      this->fe_values->shape_value(i,q_point)*
		      this->fe_values->shape_value(j,q_point)*
		      this->fe_values->JxW(q_point);	      
		  }
	      }	      
	  }
      }    

  }

  template <int dim>
  void HeatEquation<dim>::output_results ()
  {
    //fdm->pcout<<"heat output starts......";
    //MPI_Barrier(this->mpi_communicator);
    //assert(dim == 3);
    reinit_fe_values();
    DataOut<dim> data_out;
    data_out.attach_dof_handler (this->dof_handler);
    data_out.add_data_vector (fdm->theta, "theta");

    // data_out.add_data_vector (this->vsvr_total_displacement, "u");

    /* cauchy  stress -> vtk. the averaged stress of each cell */


    Vector<double> temperature(fdm->triangulation.n_active_cells());
   
    typename DoFHandler<dim>::active_cell_iterator
      cell = this->dof_handler.begin_active(),
      endc = this->dof_handler.end();
    
    for (unsigned int index = 0; cell!=endc; ++cell, ++index)
      if (cell->is_locally_owned())
	{	  
	  this->fe_values->reinit (cell);

	  double  temperature_of_this_qpt, accumulated_temperature= 0;

	  PointHistory<dim>* ph 
	    =reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

	  for(unsigned int q = 0;
	      q<fdm->quadrature_formula.size();++q)
	    {
	      temperature_of_this_qpt = ph[q].theta;
	      accumulated_temperature += temperature_of_this_qpt;
	    }

	  temperature(index) = accumulated_temperature/
	    (fdm->quadrature_formula.size());	  
	}

    data_out.add_data_vector(temperature,"thetaGauss");

    Vector<float> subdomain (fdm->triangulation.n_active_cells());
    //std::cout<<"subdomain.size = "<<subdomain.size()<<std::endl;
    for (unsigned int i=0; i<subdomain.size(); ++i)
      {
	subdomain(i) = fdm->triangulation.locally_owned_subdomain();
	//std::cout<<"subdomain("<<i<<") ="<< subdomain(i)<<std::endl;
      }

    data_out.add_data_vector (subdomain, "subdomain");

    data_out.build_patches ();
   
    const std::string filename = ("hsvr_solution-" +
                                  Utilities::int_to_string (fdm->tstep, MAX_STEPS) +
                                  "." +
                                  Utilities::int_to_string
                                  (fdm->triangulation.locally_owned_subdomain(), MAX_PROCS));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);


    if (Utilities::MPI::this_mpi_process(fdm->mpi_communicator) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i=0;
             i<Utilities::MPI::n_mpi_processes(fdm->mpi_communicator);
             ++i)
          filenames.push_back ("hsvr_solution-" +
                               Utilities::int_to_string (fdm->tstep, MAX_STEPS) +
                               "." +
                               Utilities::int_to_string (i, MAX_PROCS) +
                               ".vtu");

        std::ofstream master_output ((filename + ".pvtu").c_str());
        data_out.write_pvtu_record (master_output, filenames);
      }

    //fdm->pcout<<"ends"<<std::endl;
  }
  


  template<int dim>
  void HeatEquation<dim>::set_rigid_dof()
  {
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::map<types::global_dof_index, bool> dof_touched;
    
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    for(; cell!=endc;++cell)
      if(cell->is_locally_owned())	
	{
	  cell->get_dof_indices(local_dof_indices);
	  for(unsigned int i= 0;i<dofs_per_cell;i++)	 
	    dof_touched[local_dof_indices[i]] = false;	     
	}

    for(cell = dof_handler.begin_active(); cell!=endc; ++cell)
      if(cell->is_locally_owned())
	{
	  cell->get_dof_indices(local_dof_indices);

	  for(unsigned int i = 0;i<dofs_per_cell;i++)	{       
	  	   
	    if(dof_touched[local_dof_indices[i]] == false)
	      {
		// std::cerr<<"local_dof_indices["<<i<<"] = "<<local_dof_indices[i] <<"is false"<<std::endl;
		dof_touched[local_dof_indices[i]] = true;	  

		for (unsigned int v = 0;
		     v<GeometryInfo<dim>::vertices_per_cell;++v)
		  {
		    for(unsigned int d = 0;d<1;d++)
		      {		
			if(local_dof_indices[i] == 
			   cell->vertex_dof_index(v,d))
			  {			    		
			    Point<dim> the_vertex = cell->vertex(v);

			    if(dof_belongs_to_fix_translation(the_vertex,d))
			       fix_translation_dofs.push_back(local_dof_indices[i]);
			  }
		      }
		  }	
	      }	   
	  }	 
	}

    fdm->pcout<<"Domain 0: h fix_rigid_dofs = ";
    for(unsigned int i = 0;i<fix_rigid_dofs.size();i++)
      fdm->pcout<<fix_rigid_dofs[i]<<" ";
    fdm->pcout<<"......"<<std::endl;
  }


  template<int dim>
  bool HeatEquation<dim>:: dof_belongs_to_fix_translation(const Point<dim>& p,
						 unsigned int d)
  {
    assert(d==d);
    /* 
       set fix_translation_dofs 
       3d: Fix x,y,z degree for [xmax,ymin,zmin]
       2d: Fix x,y degree for [xmin,ymin]    
    */                  

    if(dim == 2)   
      {
	double xmin,ymin;//ymax;
	xmin = fdm->limits[0][0]; 
	ymin = fdm->limits[1][0]; 	

	if((fabs(p[0]-xmin)<1e-12 && 
	    fabs(p[1]-ymin)<1e-12 ))	  
	  return true;
      }
    else if(dim == 3)
      {			     
	double xmax,ymin,zmin;//,zmax;

	ymin = fdm->limits[1][0]; 
	zmin = fdm->limits[2][0];
	xmax = fdm->limits[0][1]; 

	if((fabs(p[0]-xmax)<1e-12 && 
	    fabs(p[1]-ymin)<1e-12 &&
	    fabs(p[2]-zmin)<1e-12))
	  {
	    return true;	   
	  }
      }

    return false;
  }


  template<int dim>
  void HeatEquation<dim>::initialize()
  {            
    //fdm->pcout<<"f initialize f......"<<std::endl;
    VectorTools::interpolate(//this->mapping,
			     dof_handler,
			     fdm->initHeatFunction,
			     this->locally_owned_theta);
    set_rigid_dof();    
    fdm->theta = this->locally_owned_theta;      
    //fdm->pcout<<"ends"<<std::endl;
  }
}
