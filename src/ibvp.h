#ifndef IBVP_H
#define IBVP_H

/*
 * Dislocation free.
 * put simple shear traction.
 * traction history is load and hold
 * total time 10
 */

namespace FDM
{
  using namespace dealii;

  template<int dim> class Fdm;

  template<int dim> 
    class NeumannLoadProfileFunc : public Function<dim> {
  public:  
  NeumannLoadProfileFunc() : Function<dim>(dim) { }
    inline double evaluate_load(double clock, Fdm<dim>* fdm) {
      assert(clock >= 0);
      double targetLoad = fdm->neumann_target_load;
      double dLT = fdm->neumann_load_period;    
      if (clock < dLT)  return clock * (targetLoad / dLT);
      else return 0;
    }
    inline double evaluate_load_rate(double clock, Fdm<dim>* fdm) {
      assert(clock >= 0);
      double targetLoad = fdm->neumann_target_load;
      double dLT = fdm->neumann_load_period;
      if (clock < dLT) return targetLoad / dLT;
      else return 0;
    }

  private:
    double targetLoad;
    double dLT;
  };


  template<int dim>
    class DirichletLoadProfileFunc : public Function<dim> {
  public:
  DirichletLoadProfileFunc() : Function<dim>(dim) { }
    inline double evaluate_load(double clock, Fdm<dim>* fdm) {
      assert(clock >= 0);
      double targetLoad = fdm->dirichlet_target_load;
      double dLT = fdm->dirichlet_load_period;
      if (clock < dLT) return clock * (targetLoad / dLT);
      else return targetLoad;
    }
    inline double evaluate_load_rate(double clock, Fdm<dim>* fdm) {
      assert(clock >= 0);
      double targetLoad = fdm->dirichlet_target_load;
      double dLT = fdm->dirichlet_load_period;
      if (clock < dLT) return targetLoad / dLT ;
      else return 0;
    }
  private:
    double targetLoad;
    double dLT;
  };


  /******************************************************/

  template<int dim>
    class DirichletDisplFunction : public Function<dim> {
  public:
  DirichletDisplFunction(Fdm<dim>* fdm, unsigned int bd) : Function<dim>(dim) {
      this->fdm = fdm;
      this->load_profile = new DirichletLoadProfileFunc<dim>();
      this->bd = bd;
    }

    inline std::vector<bool> set_component_mask(unsigned int bd) {    
      std::vector<bool> cm;
      for (int i = 0; i < dim; i++)
	cm.push_back(true); //fix dirichlet on all components

#if 0  //simple tensile
      if (bd == 5 || bd == 6) {
	cm[0] = false;
	cm[1] = false;
      }
#endif

#if 1  //simple shear in x-y
      if (bd == 5 || bd == 6) {
	cm[0] = false;
	cm[1] = false;
      }
#endif
      
      return cm;
    }

    inline void vector_value(const Point<dim> &p, Vector<double> &value) const {
      assert(p==p);
      std::fill(value.begin(), value.end(), 0.0);

#if 0  // simple tensile
      if ( bd == 6 )
	value(2) = load_profile->evaluate_load_rate(fdm->clock,fdm);	    
#endif

#if 1  // simple shear
      if (bd == 1 || bd == 2 || bd == 3 || bd == 4) {
	double v_max = load_profile->evaluate_load_rate(fdm->clock, fdm);
	value[0]  = v_max * (p[1]  - fdm->Point2[1]) / (fdm->Point1[1] - fdm->Point2[1]);

	/* std::cerr<<"p = "<<p<<std::endl; */
	/* std::cerr<<"value = "<<value<<std::endl; */
      }

#endif
    }
    
  private:
    Fdm<dim>* fdm;
    DirichletLoadProfileFunc<dim>* load_profile;
    unsigned int bd;
  };
  

  template<int dim>
    class NeumannFunction : public Function<dim> {
  public:
  NeumannFunction(Fdm<dim>* fdm) : Function<dim>(dim) {
      this->fdm = fdm;
      this->load_profile = new NeumannLoadProfileFunc<dim>();
    }

    inline Tensor<1, dim> vector_value(const Point<dim> &p, const Point<dim>& normal) const {
      throw ("this should not be called!");
    }

    inline Tensor<1, dim> vector_value(const Point<dim> &p, const Point<dim>& normal, unsigned int bd) const {
      assert(p == p && normal == normal && bd == bd);
      Tensor<2,dim> du; 
      /* du[dim-2][dim-1] = gamma; */
      /* du[dim-1][dim-2] = gamma; */
      return du*normal; 
      // Tensor<1, dim> t;
      // return t;
      //return t;
    }

    inline Tensor<1, dim> vector_value_rate(const Point<dim> &p,const Point<dim>& normal, unsigned int bd) const {
      assert(p == p && bd == bd && normal == normal);
      Tensor<1, dim> t;
      return t;
    }

  private:
    Fdm<dim>* fdm;
    NeumannLoadProfileFunc<dim>* load_profile;
  };



  /*************************************************

  Set initial conditions of alpha, f, v


  **************************************************/


  template<int dim>
    class InitDisplFunction : public Function<dim> {
  public:

  InitDisplFunction(Fdm<dim>* fdm) : Function<dim>(dim) {
      this->fdm = fdm;
    }

    inline void vector_value(const Point<dim> &p, Vector<double> &value) const {
      assert(p == p);
      for (unsigned int i = 0; i < value.size(); i++)
	value(i) = 0.0;
    }
  private:
    Fdm<dim>* fdm;
  };


  template<int dim>
    class InitHeatFunction : public Function<dim>
    {
    public:

    InitHeatFunction(Fdm<dim>* fdm)
      : Function<dim>(1)
	{
	  this->fdm = fdm;
	}

      inline double value(const Point<dim> &p, const unsigned int i ) const
      {
	assert(p == p &&
	       i==i);

	return fdm->T0[0];
      }
    private:
      Fdm<dim>* fdm;
    };

  template<int dim>
    class InitDislocationFunction : public Function<dim> {
  public:
  InitDislocationFunction(Fdm<dim>* fdm) : Function<dim>(dim * dim) {
      this->fdm = fdm;
    }

    inline void vector_value(const Point<dim> &p, Vector<double> &value) const {
      assert(p == p);
      value(2) = fdm->NeyNorm;

      double r = sqrt(p[0]*p[0]+p[1]*p[1]);
      std::cerr<<"r = "<<r<<std::endl;
      for (unsigned int i = 0; i < value.size(); i++) {
        if (i == 2) //alpha[1,3]
          if (r < 1) {
            value(i) = (fdm->NeyNorm) * (1 - tanh(r) * tanh(r));
	    std::cerr<<"value("<<i<<")="<<value(i)<<std::endl;
	  }
          else
            value(i) = 0.0;
      }
    }
  private:
    Fdm<dim>* fdm;
  };
 
}
#endif



