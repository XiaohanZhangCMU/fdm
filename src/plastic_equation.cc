#include "plastic_equation.h"

namespace FDM
{
  using namespace dealii;

  template<int dim>
  int PlasticEquation<dim>::drive(Commands cmd) {    
    switch(cmd)
      {
      case INITL:
	initialize();   
	break;      
      case SOLVE:
	solve();
	break;
      default:
	std::cerr<<"fsvr::not implemented command"<<std::endl;
	exit(1);
      }
    return 0;
  }

  //Setup layer grids and dofs. Initialize phi and omega values.  
  template<int dim>
  int PlasticEquation<dim>::initialize()
  {
    setup_system();

    distribute_memories();

    initialize_alpha();

    return 0;
  }

  template<int dim>
  int PlasticEquation<dim>::setup_system()
  {
    //set up a slip plane on y = 0 
    double Lx = fdm->limits[0][1]-fdm->limits[0][0];
    double Lz = fdm->limits[2][1]-fdm->limits[2][0];
    int nNx = fdm->repetitions[0];
    int nNz = fdm->repetitions[2];
    
    for (int j = 0; j< nNz;j++)
      for (int i = 0; i< nNx; i++) {
	imagNodes[i+j*nNx][0] = fdm->limits[1][0]+(i+0.5) * Lx/nNx;
	imagNodes[i+j*nNx][1] = fdm->limits[2][0]+(j+0.5) * Lz/nNz;
      }
    return 0;
  }

  template<int dim>
  int PlasticEquation<dim>::solve()
  {

    //interpolate layer stress
    if (tau_g_t()!=0) 
      { 
	std::cout<<"upwind->tau_g_t() fails at step " << step << std::endl;
	return -1;
      }

    double* extra = new double[n];
    double* alpha = new double[n];
    double* vr = new double[n];

    step = 0;

    std::vector<double> vec_phi_dot;
    std::vector<double> phi_t_last;
  
    for(int node = 0;node<n;node++)  phi_t_last.push_back(0.0);       

    evaluateLoad();
    
    //get velocity term vl[]
    for(int node =0;node<n;node++)
      {	  	  
	phi_differentiate(node); 
	
	F_t[node] = fabs(phi_x[node]);
	  
	double s = SHIFT;
	double p = phi_t[node];
	double e = ela_s[node];
	
	if(abs(s)>=Eps)            
	  extra[node] = eta_over_phi(s,p,e);

	vr[node] = 1.0/B * (tau_t[node]+extra[node]);
	vl[node] = 1.0/B * (B*vr[node]+epsilon*phi_xx[node]); 
      }
  
    semi_implicit_check_dtime();

    for(int node = 0;node<n;node++)
      {
	F_phi_x_c (node);  //upwinding, get F_t = |phi_x|	  
	double Fv = pow(F_t[node],1+v);	
	phi_dot[node] = Fv*vr[node]*dtime + phi_t[node]; //rhs
	alpha[node] = Fv*epsilon*dtime/B/dx/dx;
      }
	  
    solve_phi_t(alpha,phi_dot);

    for(int node = 0;node<n;node++)
      {
	if(fdm_node_at_boundary(node)==0){ 
	  phi_t[node] = P[node];	      

	  //APPLY Neumann BC 
	  if(fdm_node_at_boundary(node)==-1){
	    phi_t[node] = P[node+1];
	  }
	  if(fdm_node_at_boundary(node)==1){
	    phi_t[node] = P[node-1];
	  }
	}
      }

    vec_phi_dot.clear(); 
      
    for(int node = 0;node<n;node++)
      {	  
	vec_phi_dot.push_back( (phi_t[node] - phi_t_last[node])/dtime );
      }      
     
    phi_t_last.clear();

    for(int node = 0;node<n;node++)
      phi_t_last.push_back(phi_t[node]);
  
    //double nrm = norm(vec_phi_dot);
    nrm = norm(vec_phi_dot);     

    //pump.close();vump.close();zump.close();
    if(extra!=NULL)   delete[] extra;
    if(alpha!=NULL)   delete[] alpha;
    if(vr!=NULL)      delete[] vr;
    return 0;
  }


  template<int dim>
  int PlasticEquation<dim>::tau_g_t()
  {
    double x, y;
    for (int i= 0;i<n;i++)
      {
	x = imagNodes[i][0]; y = imagNodes[i][1];
	tau_t[i] = interpolate_layer_stress(x,y);
      }
    return 0;
  } 

  //Find all elements along the thickness. 
  //Get stress[0][1], stress[3][2] of each element
 
  template<int dim>
  int PlasticEquation<dim>::interpolate_layer_stress(double x, double y)
  {
    assert(x==y);
    return 0;
  }
 
#if 0
  template<int dim>
  void PlasticEquation<dim>::
  solve_phi_t(double* alpha, double* phi_dot)
  {
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> PM;

    Vector<double>       PV;

    DynamicSparsityPattern dsp(n);
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);
    PM.reinit (sparsity_pattern);
    P.reinit (n);
    PV.reinit (n);

    for(int i = 0;i<n;i++)
      {
	if(i == 0)
	  {
	    double pm[3][3] = { { 1,0,0 },{0,0,0},{0,0,0}};
	    double pv[] = {phi_dot[i],0,0};
	    int pm_dof[] = {i,i+1,i+2};
	    //PM.assemble(3,pm_dof,(double*)pm,PV,pv);	  
	  }
	else if(i == n-1)
	  {
	    double pm[3][3] = {{0,0,0},{0,0,0},{0,0, 1} };
	    double pv[] = {0,0,phi_dot[i]};
	    int pm_dof[] = {i-2,i-1,i};
	    //PM.assemble(3,pm_dof,(double*)pm,PV,pv);	  
	  }      
	else
	  {
	    if(i==n-2)
	      {
		double pm[2][2] = {{1+2*alpha[i], -alpha[i]}, {0, 0}};
		double pv[] = {phi_dot[i],0};
		int pm_dof[] = {i,i+1};
		//PM.assemble(2, pm_dof,(double*)pm, PV, pv);		  
	      }
	    else
	      {
		double pm[2][2] = {{1+2*alpha[i], -alpha[i]}, {-alpha[i+1], 0} };
		double pv[] = {phi_dot[i],0};
		int pm_dof[] = {i,i+1};
		//PM.assemble(2, pm_dof,(double*)pm, PV, pv);		  
	      }	  
	  }
      } 

    SolverControl           solver_control (1000, 1e-12);
    SolverCG<>              solver (solver_control);
    solver.solve (PM, P, PV,
		  PreconditionIdentity());
  }
#endif


  template<int dim>
  void PlasticEquation<dim>::
  solve_phi_t(double* alpha, double* phi_dot)
  {
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> PM;

    Vector<double>       PV;

    SparsityPattern dsp(n,n,3);
    PM.reinit (dsp);
    P.reinit (n);
    PV.reinit (n);

    for(int i = 0;i<n;i++)
      {
	if(i == 0)
	  {
	    // double pm[3][3] = { { 1,0,0 },{0,0,0},{0,0,0}};
	    // double pv[] = {phi_dot[i],0,0};
	    // int pm_dof[] = {i,i+1,i+2};
	    //PM.assemble(3,pm_dof,(double*)pm,PV,pv);	  
	    PM.add(i,i,1);
	    PV(i) = PV(i)+phi_dot[i];
	  }
	else if(i == n-1)
	  {
	    // double pm[3][3] = {{0,0,0},{0,0,0},{0,0, 1} };
	    // double pv[] = {0,0,phi_dot[i]};
	    // int pm_dof[] = {i-2,i-1,i};
	    //PM.assemble(3,pm_dof,(double*)pm,PV,pv);	  
	    PM.add(i,i,1);
	    PV(i) = PV(i)+phi_dot[i];
	  }      
	else
	  {
	    if(i==n-2)
	      {
		// double pm[2][2] = {{1+2*alpha[i], -alpha[i]}, {0, 0}};
		// double pv[] = {phi_dot[i],0};
		// int pm_dof[] = {i,i+1};
		//PM.assemble(2, pm_dof,(double*)pm, PV, pv);		  
		PM.add(i,i, 1+2*alpha[i]);
		PM.add(i,i+1, -alpha[i]);
		PV(i) = PV(i)+ phi_dot[i];
	      }
	    else
	      {
		// double pm[2][2] = {{1+2*alpha[i], -alpha[i]}, {-alpha[i+1], 0} };
		// double pv[] = {phi_dot[i],0};
		// int pm_dof[] = {i,i+1};
		PM.add(i,i,1+2*alpha[i]);
		PM.add(i,i+1,-alpha[i]);
		PM.add(i+1,i,-alpha[i+1]);
		PV(i) = PV(i)+phi_dot[i];
		//PM.assemble(2, pm_dof,(double*)pm, PV, pv);		  
	      }	  
	  }
      } 

    SolverControl           solver_control (1000, 1e-12);
    SolverCG<>              solver (solver_control);
    solver.solve (PM, P, PV,
		  PreconditionIdentity());
  }

  template<int dim>
  double PlasticEquation<dim>::evaluateLoad()
  {
    return 0;
  }

  //GL energy: -(d/d_imag * mu * s^2) / (4 * pi^2) * cos(2pi p/(s/num_pile_up))
  template<int dim>
  double PlasticEquation<dim>::eta_over_phi(double s, double p, double e)
  {
    assert(e==e);
    double val;
    double MU = 1;
    double d = 2;
    double d_imag = 2;
    double PI = 3.1415926;
    //double d_imag = d; //imaginary layer thickness
    //epsilon = 0.1;// MU/2.0/PI *s/2.0;
    val = (-1.0)*MU/2.0/PI * s*d/d_imag * (sin(2*PI*p/(s)));
    //std::cout<<"epsilon = "<<epsilon<<std::endl;
    return val;
  }


  template<int dim>
  double PlasticEquation<dim>::eta_over_phi_prime(double s, double p, double e)
  {
    assert(e==e);
    double val;
    double MU = 1;
    double d = 2;
    double d_imag = 2;
    double PI = 3.1415926;
    
    //double d_imag = d; //imaginary layer thickness
    //epsilon = 0.1;// MU/2.0/PI *s/2.0;
    val = (-1.0)*MU *d/d_imag * (cos(2*PI*p/(s)));
    //std::cout<<"epsilon = "<<epsilon<<std::endl;
    return val;
  }

  template<int dim>
  int PlasticEquation<dim>::fdm_node_at_boundary(int node)
  {
    assert(node<n && node>=0);
    // std::cerr<<"node=  "<<node<<std::endl;
    // std::cerr<<"dist left=  "<<fabs(imagNodes[node]-dx-xmin)<<std::endl;
    // std::cerr<<"dist right=  "<<fabs(imagNodes[node]+dx-xmax)<<std::endl;
    if(fabs(imagNodes[node][0]-dx/2.0 -xmin)<Eps)
      return -1;
    if(fabs(imagNodes[node][0]+dx/2.0 -xmax)<Eps)    
      return 1;    
    else
      return 0;
  }

  template<int dim>
  int PlasticEquation<dim>::phi_differentiate(int node)
  { 
    double h1, h2, alpha, gamma, beta, r, s, t;

    if(fdm_node_at_boundary(node) == 0)
      {
	h1 = imagNodes[node][0] - imagNodes[node-1][0];
	h2 = imagNodes[node+1][0]-imagNodes[node][0];

	if(h1<=0 || h2<=0) std::cout<<"node = "<<node<<"h1 = "<<h1 <<"h2 = "<<h2<<std::endl; 
	//assert(h1<=0);assert(h2<=0);
	//throw("dx is negative.");
      }
    else if(fdm_node_at_boundary(node)==-1)
      {
	h1 = imagNodes[node+1][0] - imagNodes[node][0];
	h2 = imagNodes[node+2][0]-imagNodes[node+1][0];   
      }
    else //(node == n-1)
      {
	h1 = imagNodes[node-1][0] - imagNodes[node-2][0];
	h2 = imagNodes[node][0]-imagNodes[node-1][0];   
      }
 
    alpha = -h2/(h1*(h1+h2));
    gamma = h1/(h2*(h1+h2));
    beta = -alpha-gamma;
      
    r = 2.0/(h1*(h1+h2));
    t = 2.0/(h2*(h1+h2));
    s = -r-t;
      
    if(fdm_node_at_boundary(node)==0)
      {
	phi_x[node] = alpha*phi_t[node-1]+beta*phi_t[node]+gamma*phi_t[node+1];
	phi_xx[node] = r*phi_t[node-1]+s*phi_t[node] + t*phi_t[node+1];
      }
    else if(fdm_node_at_boundary(node)==-1)
      {
	phi_x[node] = alpha*phi_t[node]+beta*phi_t[node+1]+gamma*phi_t[node+2];
	phi_xx[node] = r*phi_t[node]+s*phi_t[node+1] + t*phi_t[node+2];
      }
    else//(node==n-1)
      {
	phi_x[node] = alpha*phi_t[node-2]+beta*phi_t[node-1]+gamma*phi_t[node];
	phi_xx[node] = r*phi_t[node-2]+s*phi_t[node-1] + t*phi_t[node];
      }
    return 0;
  } 

  template<int dim>
  int PlasticEquation<dim>::semi_implicit_check_dtime()
  {
    double st = 0.5;
    double pr_val;
   
    double c_max = 0.0;

    double dtime_min = st;
    double sign_F_dash = 0;

    double dtime_node[3];
  
    for(int node = 0;node<n;node++)
      {
	if(phi_x[node]>0)
	  sign_F_dash = 1.0;
	if(phi_x[node]<0)
	  sign_F_dash = -1.0;
	if(phi_x[node]==0)
	  sign_F_dash = 0.0;
      
	double F_v = pow(fabs(phi_x[node]),v);
     
	c[node] = -((sign_F_dash))*F_v*(1+v)*(vl[node]);

	//std::cout<<"c["<<node<<"]="<<c[node]<<"; dx = "<<dx <<"; F_v = "<<F_v<<"; phi_x["<<node<<"]="<<phi_x[node]<<std::endl;
      
	if(fabs(c[node])>fabs(c_max))
	  c_max = fabs(c[node]);
	
	dtime_node[0] = dx/c[node];
      
	//dtime_node[1] = (dx*dx)*B/(2.0*epsilon*F_t[node]);
      
	pr_val = fabs(eta_over_phi_prime(SHIFT,phi_t[node],0));
	//pr_val+=norm(tau_t);
	dtime_node[1] = B/(F_t[node] * pr_val );

	dtime_node[2] = dx/C1;
	// for(int i = 0;i<2;i++){
	//  	std::cout<<"dx= "<<dx<<" c_max = "<<c_max<<std::endl;
	//  	std::cout<<"dtime_node["<<i<<"] = "<<dtime_node[i]<<std::endl;
	// }

	double dtime_node_min = st;
	int sz = 2; 

	for(int i =0;i<sz;i++)
	  {
	    if(fabs(dtime_node[i])<dtime_node_min)
	      dtime_node_min = fabs(dtime_node[i]);
	  }
      
	if(fabs(dtime_node_min)<dtime_min)
	  dtime_min = fabs(dtime_node_min);
      }

    dtime = 0.5* dtime_min;

    for(int node = 0;node<n;node++)
      {
	if(fabs(c[node])/fabs(c_max) <= 1e-6)
	  c[node] = 0;
      }
    return 0;
  }


  template<int dim>
  int PlasticEquation<dim>::F_phi_x_c( int node)
  {    
    if(fdm_node_at_boundary(node)==0)
      {
	if(c[node]==0)
	  {
	    F_t[node] = fabs((phi_t[node+1]-phi_t[node-1])/(2.0*dx));
	    //F_t[node] = 0;
	    //	  std::cout<<"c[node] == 0"<<std::endl;
	  }

	if(c[node]<0)
	  {
	    F_t[node] = fabs((phi_t[node+1]-phi_t[node])/dx);
	  }
	if(c[node]>0)
	  {
	    F_t[node] = fabs((phi_t[node]-phi_t[node-1])/dx);
	  }	
      }
    else
      F_t[node] = 0.0;

    return 0;
  }

  template<int dim>
  int PlasticEquation<dim>::distribute_memories()
  {
    F_t = new double[n];  tau_t =  new double[n];  tau_dash_t =  new double[n];  phi =  new double[n]; prev_phi = new double[n];

    phi_t =  new double[n];  phi_dot =  new double[n];  phi_x =  new double[n];  phi_xx = new double[n]; 

    c = new double[n]; ela_s = new double[n]; vl = new double[n]; slip = new double[n];

    imagNodes = new double*[n];

    for(int i = 0;i<n;i++)
      {
	F_t[i] = 0; tau_t[i] = 0; tau_dash_t[i] = 0; phi[i] = 0;
      
	phi_t[i] = 0; phi_dot[i] = 0; phi_x[i] = 0;phi_xx[i] = 0;

	c[i] = 0; ela_s[i] = 0; vl[i] = 0; slip[i] = 0; prev_phi[i] = 0;

	imagNodes[i] = new double[2]; imagNodes[i][0] = imagNodes[i][1]=0.0;
      }

    //check the pointers
    if(F_t !=NULL && tau_t!=NULL && tau_dash_t!=NULL && phi!=NULL
       && phi_t !=NULL && phi_dot!=NULL && phi_x!=NULL&& phi_xx!=NULL &&c!=NULL && ela_s!=NULL && slip!=NULL && vl!=NULL && prev_phi!=NULL)  return 0;

    return -1;
  }

  template<int dim>
  int PlasticEquation<dim>::initialize_alpha()
  {
    double s = SHIFT;
    double MU = 1;
    double a = sqrt(MU/0.01/epsilon);
    double xcc = 0;
    int nNx = fdm->repetitions[0];
    int nNz = fdm->repetitions[2];
    for (int j = 0;j<nNz;j++)
      for (int i = 0;i<nNx;i++)
	phi_t[i+j*nNx] = tanh(a*(imagNodes[i][0]-xcc))*s/2.0+s/2.0;
    return 0;
  }

  template<int dim>
  double PlasticEquation<dim>::norm(std::vector<double> vec_phi_dot)
  {
    return vec_phi_dot.size();
  }
}
