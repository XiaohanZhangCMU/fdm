#ifndef PLASTIC_EQUATION_H
#define PLASTIC_EQUATION_H

#include "fdm.h"
#include "displ_equation.h"
//#include "heat_equation.h"

namespace FDM
{
  using namespace dealii;

  template<int dim> class DisplEquation;
  template<int dim> class HeatEquation;

  template<int dim>
    class PlasticEquation //: public Equation<dim> 
    {
    public:
      PlasticEquation(Fdm<dim+1>* fdm)	
	:
      fdm(fdm)
      {

      }
  
      ~PlasticEquation()
	{	  	  
	  if(F_t!=NULL) delete[] F_t;
	  if(tau_t!=NULL) delete[] tau_t;
	  if(tau_dash_t!=NULL) delete[] tau_dash_t;
	  if(phi!=NULL) delete[] phi;
	  if(phi_t!=NULL) delete[] phi_t;
	  if(phi_dot!=NULL) delete[] phi_dot;
	  if(phi_x!=NULL) delete[] phi_x;
	  if(phi_xx!=NULL) delete[] phi_xx;
	  if(c!=NULL) delete[] c;
	  if(slip!=NULL) delete[] slip;
	  if(ela_s!=NULL) delete[] ela_s;
	  if(vl!=NULL) delete[] vl;
	  if(prev_phi!=NULL) delete[] prev_phi;	 	  
	}

      void set_coupled(DisplEquation<dim+1>* ueqn)
	{
	  this->ueqn = ueqn;
	}

      
    public:
      
      Fdm<dim+1>* fdm;
      DisplEquation<dim+1>* ueqn;
      HeatEquation<dim+1>* heqn;

      /* Triangulation<2>     triangulation; */
      /* FESystem<2>           fe; */
      /* DoFHandler<2>        dof_handler; */

      Vector<double>       P; //solution
      
      double* F_t;  double* tau_t;   double* tau_dash_t;  double* phi; double* ela_s;
      double* phi_t;  double* phi_dot;  double* phi_x;   double* phi_xx; 
      double* c; double* vl; double* slip; double* prev_phi;
      double** imagNodes;

      double B;
      double v;
      double nrm;
      double dtime;
      double Eps;
      double xmin, xmax, zmin, zmax, dx;
      double clock;
      double dclock;
      int step;
      int n;
      
      const double SHIFT = 0.5;
      const double epsilon = 0.025;
      const double C1 = 0;
      
    public:
      
      int drive(Commands cmd);      
      double evaluateLoad();
      double eta_over_phi(double,double,double);
      double eta_over_phi_prime(double,double,double);
      int initialize();
      double norm(std::vector<double>);
      int initialize_alpha();
      int setup_system();
      int solve();
      void solve_phi_t(double* alpha, double* phi_dot);
      double eta_over_phi6(double s, double p, double e);
      double eta_over_phi_prime6(double s, double p, double e);
      int fdm_node_at_boundary(int node);
      int phi_differentiate(int node);
      int semi_implicit_check_dtime();
      int F_phi_x_c(int node);
      int distribute_memories();
      int tau_g_t();
      int interpolate_layer_stress(double, double);
    };
}
#endif
