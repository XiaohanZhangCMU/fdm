/*
 * FDM contains FDM problem constitutive equaitons and standard tensors
 */

#ifndef FDM_H
#define FDM_H

#include "common.h"
#include "ibvp.h"

namespace FDM
{
  using namespace dealii;

  template<int dim>
    struct PointHistory
    {
      Tensor<2, dim> dpsi, p_dpsi, pp_dpsi;
    };


  template<int dim> 
    class Fdm {
   
  public:
  Fdm(int refine_cycle, std::ostream& loger)
    :
    mpi_communicator (MPI_COMM_WORLD),
      pcout (loger, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
      computing_timer (pcout,
		       TimerOutput::summary,
		       TimerOutput::wall_times),
      triangulation (mpi_communicator,
		     typename Triangulation<dim>::MeshSmoothing
		     (Triangulation<dim>::smoothing_on_refinement
		      | Triangulation<dim>::smoothing_on_coarsening)),

      quadrature_formula(3),
      face_quadrature_formula(3),
      neumannFunction(NeumannFunction<dim>(this))
	{
	  read_in_sim_consts();

	  limits = new double*[dim];
	  sum_rf = new double[dim];
	  for (int d = 0; d < dim; d++)
	    limits[d] = new double[2];

	  grid(refine_cycle);
	  set_b_c();
	  set_material_id();
	  setup_quadrature_point_history ();
	  set_standard_tensors();

	  for (unsigned int i = 0; i < b_c_codes.size(); i++) {
	    if (b_c_codes[i] < NMNDRCSPLIT) {

	      //V
	      DirichletDisplFunction<dim>* dfv = new DirichletDisplFunction<dim>(this, b_c_codes[i]);
	      velocity_dirichlet_boundary_functions[b_c_codes[i]] = dfv;

	      velocity_component_masks[b_c_codes[i]] =
		dfv->set_component_mask(b_c_codes[i]);
	    }

	  }
	}


    ~Fdm()
      {
	for (unsigned int bidcnt = 0; bidcnt < b_c_codes.size(); bidcnt++) {
	  if (b_c_codes[bidcnt] < NMNDRCSPLIT) {
	    if (velocity_dirichlet_boundary_functions[b_c_codes[bidcnt]] != NULL)
	      delete velocity_dirichlet_boundary_functions[b_c_codes[bidcnt]];
	  }
	}

	for (int d = 0; d < dim; d++)
	  if (limits[d] != NULL)
	    delete[] limits[d];
	if (limits != NULL)
	  delete[] limits;
	if (sum_rf != NULL)
	  delete sum_rf;
      }

    /* Grid for the problem, every equation uses the same grid. */
    int mesher;
    int mesoscale;
    int grid(int);

    /* containers for boundary limits */
    double** limits;
    double* sum_rf;
    double element_size;

    /* shared the same mpi world & I/O/Timer */

    MPI_Comm mpi_communicator;
    ConditionalOStream    pcout;
    TimerOutput computing_timer;

    /* share the same triangulation */
    parallel::distributed::Triangulation<dim> triangulation;

    /* Solution data */
    LA::MPI::Vector U;
    LA::MPI::Vector RF; //ghosted reaction force sumover all time steps
    LA::MPI::Vector theta;

    QGauss<dim> quadrature_formula;
    QGauss < dim - 1 > face_quadrature_formula;

    //std::vector<double> step_records;

  public:
    SymmetricTensor<2, dim> I;
    SymmetricTensor<4, dim> IxI;
    SymmetricTensor<4, dim> II;
    Tensor<3, dim> eE;

    std::vector<Tensor<4, dim> > sL;

    //DirichletDisplFunction<dim> dirichletDisplFunction;
    typename FunctionMap<dim>::type velocity_dirichlet_boundary_functions;

    std::map<unsigned int, std::vector<bool> > velocity_component_masks;

    /* std::vector<DirichletFunction<dim> > dirichletFunctions; */
    NeumannFunction<dim> neumannFunction;

    /* read following params from const.in */
    std::vector<int> b_c_codes;

    std::vector<double> Point1;
    std::vector<double> Point2;
    std::vector<unsigned int> repts;
    std::vector<PointHistory<dim> > quadrature_point_history;
    std::vector<std::vector<Tensor<1, dim> > >  FCC;
    std::vector<std::vector<Tensor<1, dim> > >  BCC;

    inline double get_lambda(const Point<dim>& P,
			     const unsigned int mid)
    {
      assert(P == P&&mid==mid);

      return lambda[mid];
    }

    inline double get_mu(const Point<dim>& P,
			 const unsigned int mid)
    {
      assert(P == P&&mid==mid);

      return mu[mid];
    }


    inline double get_rho(const Point<dim>& P,
			  const unsigned int mid)
    {
      assert(P == P&&mid==mid);

      return Rho[mid];
    }

    inline double get_max_vp()
    {
      assert(Vp.size() > 0);
      double vp = Vp[0];
      for (int i = 1; i < Vp.size(); i++)
	if (Vp[i] > vp)
	  vp = Vp[i];
      return vp;
    }

    inline double get_ThetaDiffusion(const Point<dim>& P,
				     const unsigned int mid)
    {
      assert(P == P&&mid==mid);

      return ThetaDiffusion[mid];
    }

    inline double get_ThetaWorkConvs(const Point<dim>& P,
				     const unsigned int mid)
    {
      assert(P == P&&mid==mid);

      return ThetaWorkConvs[mid];
    }

    inline double getMaxDiffusion()
    {
      assert(ThetaDiffusion.size() > 0);
      double thmx = ThetaDiffusion[0];
      for (unsigned int i = 0; i < ThetaDiffusion.size(); i++)
	if (thmx < ThetaDiffusion[i])
	  thmx = ThetaDiffusion[i];
      return thmx;
    }

    /* constitutive equations */
    inline void dPsi(const Point<dim>& P,
		     const unsigned int mid,
		     const Tensor<2, dim>& F,
		     const bool d2psiFlag,
		     Tensor<2, dim>& dpsi,
		     Tensor<4, dim>& d2psi)
    {
      assert(P == P&& mid==mid);

      unsigned int i, k, l, r, s, j, m, n;

      Tensor<2, dim> E, dest;
      dpsi = 0;
      d2psi = 0;

      E = 0.5 * (transpose(F) * F - I);

      double_contract(dest, sL[mid], E);
      dpsi = F * dest * transpose(F);

      if (!d2psiFlag)
	return;

      for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	  for (m = 0; m < dim; m++)
	    for (n = 0; n < dim; n++)

	      for (l = 0; l < dim; l++)
		for (r = 0; r < dim; r++)
		  for (s = 0; s < dim; s++)
		    for (k = 0; k < dim; k++)
		      {
			d2psi[i][j][m][n] += I[i][m] *
			  sL[mid][n][l][r][s] * E[r][s] * F[j][l] +
			  F[i][k] * sL[mid][k][n][r][s] * E[r][s] * I[j][m] +
			  F[i][k] * 0.5 * sL[mid][k][l][r][s] *
			  (F[m][s] * I[r][n] + F[m][r] * I[n][s]) * F[j][l];
		      }
    }
  
    inline
      Tensor<2, dim>
      AxV(Tensor<2, dim> A, Tensor<1, dim> V)
      {
	Tensor<2, dim> temp;
	for (int i = 0; i < dim; i++)
	  for (int m = 0; m < dim; m++)

	    for (int j = 0; j < dim; j++)
	      for (int k = 0; k < dim; k++)
		temp[i][m] += eE[m][j][k] * A[i][j] * V[k];
	return temp;
      }

    inline
      double dot_Up(Tensor<2, dim> A, Tensor<1, dim> V)
    {
      Tensor<2, dim> temp;
      for (int i = 0; i < dim; i++)
	for (int m = 0; m < dim; m++)

	  for (int j = 0; j < dim; j++)
	    for (int k = 0; k < dim; k++)
	      temp[i][m] += eE[m][j][k] * A[i][j] * V[k];
      return temp.norm();
    }


    /* helper for constitutive equations */
    void set_b_c();
    void set_material_id();
    void setup_quadrature_point_history ();
    void set_standard_tensors();
    void set_predefined_slip_system();

    /* Tensor<2,dim> tensor_cross_vector(Tensor<2,dim> A,  */
    /* 				      Tensor<1,dim> V); */

    unsigned int ShearLockSoluble;
    unsigned int b_c_split_id;
    double dclock;
    double dt; //sub evolve mesh
   
    double p_dclock;
    double pp_dclock;
    double max_dclock;
    double beta;
    double duration;

    unsigned int matid;
    std::vector<unsigned int > repetitions;
    std::vector<double> lambda;
    std::vector<double> mu;
    std::vector<double> Vp; //primary wave velocity
    std::vector<double> Vs;

    std::vector<double> youngs;
    std::vector<double> poisson;
    std::vector<double> drag;

    std::vector<double> LL ;
    std::vector<double> cc ;
    std::vector<double> mm ;
    std::vector<double> Gm0 ;
    std::vector<double> Gs ;
    std::vector<double> G0 ;
    std::vector<double> Theta0 ;
    std::vector<double> eta ;
    std::vector<double> K0 ;
    std::vector<double> Burger ;
    std::vector<double> Rho ;
    std::vector<double> T0;
    std::vector<double> TK0;
    std::vector<double> ThetaDiffusion;
    std::vector<double> ThetaWorkConvs;
    std::vector<double> DELTA;


    double msfe;
    double NeyNorm;
    double max_d_psr;
    double neumann_target_load;
    double dirichlet_target_load;
    double neumann_load_period;
    double dirichlet_load_period;
    SScurveType sstype;

    int tstep;
    int cutback_tstep;
    double clock;
    double meshsize;

    double inner_radius;
    double outer_radius;
    double extrusion;
    double layer_elm_sz;

    double left_lim;
    double right_lim;

    void read_in_sim_consts();
    void display();

    bool AllTraction;


  };



  template<int dim> void Fdm<dim>::read_in_sim_consts()
    {
      std::ifstream file( "../const_params.in" );
      if (!file.is_open())
	{
	  std::cerr << "cannot find const_params.in" << std::endl;
	  exit(1);
	}

      std::string line;
      while ( std::getline( file, line ) ) {
	std::istringstream iss( line );

	std::string result;

	if (line.find_first_of("#") != std::string::npos)
	  continue;

	if ( std::getline( iss, result , '=')) {
	  result.erase(std::remove(result.begin(), result.end(), ' '), result.end());

	  //algorithm control
	  if ( result == "ShearLockSoluble" ) {
	    std::string token;
	    std::getline( iss, token );
	    ShearLockSoluble = atoi(token.c_str());
	  }

	  //TEMPORAL SPECIFICS
	  if ( result == "dclock" ) {
	    std::string token;
	    std::getline( iss, token );
	    max_dclock = atof(token.c_str());
	    dclock = max_dclock;
	  }

	  if ( result == "duration" ) {
	    std::string token;
	    std::getline( iss, token );
	    duration = atof(token.c_str());
	  }

	  if ( result == "tstep" ) {
	    std::string token;
	    std::getline( iss, token );
	    tstep = atoi(token.c_str());
	    cutback_tstep = tstep;
	  }

	  if ( result == "mesoscale" ) {
	    std::string token;
	    std::getline( iss, token );
	    mesoscale = atoi(token.c_str());
	  }

	  if ( result == "clock" ) {
	    std::string token;
	    std::getline( iss, token );
	    clock = atof(token.c_str());
	  }

	  if ( result == "beta" ) {
	    std::string token;
	    std::getline( iss, token );
	    beta = atof(token.c_str());
	  }
	  if ( result == "NeyNorm" ) {
	    std::string token;
	    std::getline( iss, token );
	    NeyNorm = atof(token.c_str());
	  }

	  if ( result == "max_d_psr" ) {
	    std::string token;
	    std::getline( iss, token );
	    max_d_psr = atof(token.c_str());
	  }

	  //MATERIAL PARAMS
	  if ( result == "matid" ) {
	    std::string token;
	    std::getline( iss, token );
	    matid = atoi(token.c_str());
	  }

	  if ( result == "YoungsM" ) {
	    std::string token;
	    std::getline( iss, token );
	    youngs.push_back(atof(token.c_str()));
	    assert(youngs.size() == matid);
	  }
	  if ( result == "Poisson" ) {
	    std::string token;
	    std::getline( iss, token );
	    poisson.push_back(atof(token.c_str()));
	    assert(poisson.size() == matid);
	  }

	  if ( result == "drag" ) {
	    std::string token;
	    std::getline( iss, token );
	    drag.push_back(atof(token.c_str()));
	    assert(drag.size() == matid);
	  }

	  if ( result == "LL" ) {
	    std::string token;
	    std::getline( iss, token );
	    LL.push_back(atof(token.c_str()));
	    assert(LL.size() == matid);
	  }
	  if ( result == "cc" ) {
	    std::string token;
	    std::getline( iss, token );
	    cc.push_back(atof(token.c_str()));
	    assert(cc.size() == matid);
	  }
	  if ( result == "mm" ) {
	    std::string token;
	    std::getline( iss, token );
	    mm.push_back(atof(token.c_str()));
	    assert(mm.size() == matid);
	  }
	  if ( result == "Gm0" ) {
	    std::string token;
	    std::getline( iss, token );
	    Gm0.push_back(atof(token.c_str()));
	    assert(Gm0.size() == matid);
	  }
	  if ( result == "G0" ) {
	    std::string token;
	    std::getline( iss, token );
	    G0.push_back(atof(token.c_str()));
	    assert(G0.size() == matid);
	  }

	  if ( result == "Gs" ) {
	    std::string token;
	    std::getline( iss, token );
	    Gs.push_back(atof(token.c_str()));
	    assert(Gs.size() == matid);
	  }
	  if ( result == "Theta0" ) {
	    std::string token;
	    std::getline( iss, token );
	    Theta0.push_back(atof(token.c_str()));
	    assert(Theta0.size() == matid);
	  }
	  if ( result == "eta" ) {
	    std::string token;
	    std::getline( iss, token );
	    eta.push_back(atof(token.c_str()));
	    assert(eta.size() == matid);
	  }
	  if ( result == "K0" ) {
	    std::string token;
	    std::getline( iss, token );
	    K0.push_back(atof(token.c_str()));
	    assert(K0.size() == matid);
	  }
	  if ( result == "Burger" ) {
	    std::string token;
	    std::getline( iss, token );
	    Burger.push_back(atof(token.c_str()));
	    assert(Burger.size() == matid);
	  }
	  if ( result == "Rho" ) {
	    std::string token;
	    std::getline( iss, token );
	    Rho.push_back(atof(token.c_str()));
	    assert(Rho.size() == matid);
	  }
	  if ( result == "DELTA" ) {
	    std::string token;
	    std::getline( iss, token );
	    DELTA.push_back(atof(token.c_str()));
	    assert(DELTA.size() == matid);
	  }
	  if ( result == "T0" ) {
	    std::string token;
	    std::getline( iss, token );
	    T0.push_back(atof(token.c_str()));
	    assert(T0.size() == matid);
	  }
	  if ( result == "TK0" ) {
	    std::string token;
	    std::getline( iss, token );
	    TK0.push_back(atof(token.c_str()));
	    assert(TK0.size() == matid);
	  }
	  if ( result == "ThetaDiffusion" ) {
	    std::string token;
	    std::getline( iss, token );
	    ThetaDiffusion.push_back(atof(token.c_str()));
	    assert(ThetaDiffusion.size() == matid);
	  }
	  if ( result == "ThetaWorkConvs" ) {
	    std::string token;
	    std::getline( iss, token );
	    ThetaWorkConvs.push_back(atof(token.c_str()));
	    assert(ThetaWorkConvs.size() == matid);
	  }


	  //LOAD PARAMS
	  if ( result == "NeumannTargetLoad" ) {
	    std::string token;
	    std::getline( iss, token );
	    //neumann_target_load = stod(token);
	    neumann_target_load = atof(token.c_str());
	  }
	  if ( result == "DirichletTargetLoad" ) {
	    std::string token;
	    std::getline( iss, token );
	    //dirichlet_target_load = stod(token);
	    dirichlet_target_load = atof(token.c_str());
	  }

	  if ( result == "NeumannLoadPeriod" ) {
	    std::string token;
	    std::getline( iss, token );
	    // neumann_load_period = stod(token);
	    neumann_load_period = atof(token.c_str());
	  }
	  if ( result == "DirichletLoadPeriod" ) {
	    std::string token;
	    std::getline( iss, token );
	    //dirichlet_load_period = stod(token);
	    dirichlet_load_period = atof(token.c_str());
	  }
	  if ( result == "sstype" ) {
	    std::string token;
	    std::getline( iss, token );
	    sstype = (SScurveType)atoi(token.c_str());
	  }

	  //MESH SPECIFICS
	  if ( result == "mesher" ) {
	    std::string token;
	    std::getline( iss, token );
	    //mesher = stoi(token);
	    mesher = atoi(token.c_str());
	  }
	  if ( result == "meshsize" ) {
	    std::string token;
	    std::getline( iss, token );
	    //meshsize = stod(token);
	    meshsize = atof(token.c_str());
	  }
	  if ( result == "inner_radius" ) {
	    std::string token;
	    std::getline( iss, token );
	    //meshsize = stod(token);
	    inner_radius = atof(token.c_str());
	  }

	  if ( result == "outer_radius" ) {
	    std::string token;
	    std::getline( iss, token );
	    //meshsize = stod(token);
	    outer_radius = atof(token.c_str());
	  }

	  if ( result == "extrusion" ) {
	    std::string token;
	    std::getline( iss, token );
	    //meshsize = stod(token);
	    extrusion = atof(token.c_str());
	  }
	  if ( result == "layer_elm_sz" ) {
	    std::string token;
	    std::getline( iss, token );
	    //meshsize = stod(token);
	    layer_elm_sz = atof(token.c_str());
	  }

	  if ( result == "left_lim" ) {
	    std::string token;
	    std::getline( iss, token );
	    //meshsize = stod(token);
	    left_lim = atof(token.c_str());
	  }

	  if ( result == "right_lim" ) {
	    std::string token;
	    std::getline( iss, token );
	    //meshsize = stod(token);
	    right_lim = atof(token.c_str());
	  }


	  if (result == "b_c_codes") {
	    std::string token;
	    std::getline( iss, token );
	    std::stringstream ss(token);
	    int number;
	    while (ss >> number)
	      b_c_codes.push_back(number);
	  }
	  if (result == "Point1") {
	    std::string token;
	    std::getline( iss, token );
	    std::stringstream ss(token);
	    double number;
	    while (ss >> number)
	      Point1.push_back(number);
	  }
	  if (result == "Point2") {
	    std::string token;
	    std::getline( iss, token );
	    std::stringstream ss(token);
	    double number;
	    while (ss >> number)
	      Point2.push_back(number);
	  }
	  if (result == "repts") {
	    std::string token;
	    std::getline( iss, token );
	    std::stringstream ss(token);
	    unsigned int number;
	    while (ss >> number)
	      repts.push_back(number);
	  }
	}
      }

      assert(matid > 0);
      assert(youngs.size() == matid);
      assert(poisson.size() == matid);
      assert(drag.size() == matid);

      assert(LL.size() == matid);
      assert(cc.size() == matid);
      assert(mm.size() == matid);
      assert(Gm0.size() == matid);
      assert(Gs.size() == matid);
      assert(G0.size() == matid);
      assert(Theta0.size() == matid);
      assert(eta.size() == matid);
      assert(K0.size() == matid);
      assert(Burger.size() == matid);
      assert(Rho.size() == matid);

    }

  template<int dim> void Fdm<dim>::display() {
    pcout << "clock = " << clock << "; ";
    pcout << "tstep = " << tstep << "; ";
    pcout << "dclock = " << dclock << "; ";
    pcout << "duration = " << duration << "; ";

    pcout << std::endl;
    for (unsigned int i = 0; i < lambda.size(); i++)
      pcout << "lambda = " << lambda[i] << "; ";

    for (unsigned int i = 0; i < mu.size(); i++)
      pcout << "mu = " << mu[i] << "; ";
    pcout << std::endl;

    pcout << "mesher = " << mesher << "; ";
    pcout << "b_c_codes = ";
    for (std::vector<int>::iterator it = b_c_codes.begin();
	 it != b_c_codes.end(); it++)
      pcout << *it << " ";
    pcout << std::endl;
    pcout << "point 1 = ";
    for (unsigned int i = 0; i < Point1.size(); i++)
      pcout << Point1[i] << " ";

    pcout << "point 2 = ";
    for (unsigned int i = 0; i < Point2.size(); i++)
      pcout << Point2[i] << " ";
    pcout << std::endl;

    pcout << "rept.size = " << repts.size() << "; ";
    pcout << "element_size = " << element_size << "; ";
    pcout << "b_c_split = " << b_c_split_id << "; ";
    pcout << "wave speed Vp = " << Vp[0] << "; ";
    /*     pcout<<"LL = "<<LL<<"; cc = "<<"; mm"<<mm<<"; lambda = "<< */
    /* mm = 0.1 */
    /* Gm0 = 1 */
    /* Gs = 210 */
    /* G0 = 50 */
    /* Theta0 = 205 */
    /* eta = 1/3 */
    /* K0 = 20.0 */
    /* Burger = 0.00025 */

  }



  template<int dim> int Fdm<dim>::grid(int refine_cycle) {
    element_size = 1e50;

    switch (mesher) {
    case 1: {
      /* Point<dim> p1(20,20,20); */
      /* Point<dim> p2(-20,-20,-20); */

      Point<dim> p1(Point1[0], Point1[1], Point1[2]);
      Point<dim> p2(Point2[0], Point2[1], Point2[2]);

      for (int d = 0; d < dim; d++) {
	limits[d][0] = p1[d] > p2[d] ? p2[d] : p1[d];
	limits[d][1] = p1[d] > p2[d] ? p1[d] : p2[d];
	double h = (limits[d][1] - limits[d][0] ) / pow(2, refine_cycle) ;
	if (h < element_size)
	  element_size = h;
      }
      GridGenerator::hyper_rectangle(triangulation, p1, p2);
      triangulation.refine_global(refine_cycle);
      element_size = -1.0 * Utilities::MPI::max(-element_size,
						mpi_communicator);
      break;
    }
    case 2: {
      GridGenerator::hyper_cube(triangulation, -meshsize, meshsize);
      triangulation.refine_global(refine_cycle);
      element_size = 1e50;
      for (int d = 0; d < dim; d++) {
	limits[d][0] = -1.0 * meshsize;
	limits[d][1] = meshsize;
	double h = (limits[d][1] - limits[d][0] ) / pow(2, refine_cycle) ;
	if (h < element_size)
	  element_size = h;
      }
      element_size = -1.0 * Utilities::MPI::max(-element_size,
						mpi_communicator);
      break;
    }
    case 3: {
      Point<dim> p1(Point1[0], Point1[1], Point1[2]);
      Point<dim> p2(Point2[0], Point2[1], Point2[2]);

      for (int d = 0; d < dim; d++)
	repetitions.push_back(repts[d]);

      for (int d = 0; d < dim; d++) {
	limits[d][0] = p1[d] > p2[d] ? p2[d] : p1[d];
	limits[d][1] = p1[d] > p2[d] ? p1[d] : p2[d];
	double h = (limits[d][1] - limits[d][0] ) /
	  pow(repetitions[d], refine_cycle + 1);

	if (h < element_size)
	  element_size = h;
      }
      GridGenerator::subdivided_hyper_rectangle(triangulation,
						repetitions,
						p1, p2);
      triangulation.refine_global(refine_cycle);
      element_size = -1.0 * Utilities::MPI::max(-element_size,
						mpi_communicator);

      break;
    }
    case 4: {
      GridGenerator::hyper_cube_with_cylindrical_hole(triangulation,
						      inner_radius,
						      outer_radius,
						      extrusion,
						      repts[2],
						      true);
      triangulation.refine_global(refine_cycle);
      element_size = 1e50;
      for (int d = 0; d < 2; d++) {
	limits[d][0] = -outer_radius;
	limits[d][1] =  outer_radius;
	double h = (limits[d][1] - limits[d][0] ) / pow(2, refine_cycle) ;
	if (h < element_size)
	  element_size = h;
      }
      limits[2][0] = 0;
      limits[2][1] = extrusion;
      element_size = -1.0 * Utilities::MPI::max(-element_size,
						mpi_communicator);
      break;
    }
    case 5: {
      GridGenerator::hyper_cube_slit(triangulation,
				     left_lim,
				     right_lim,
				     false);
      triangulation.refine_global(refine_cycle);
      element_size = 1e50;
      for (int d = 0; d < 2; d++) {
	limits[d][0] = left_lim;
	limits[d][1] =  right_lim;
	double h = (limits[d][1] - limits[d][0] ) / pow(2, refine_cycle) ;
	if (h < element_size)
	  element_size = h;
      }
      limits[2][0] = left_lim / 2.0;
      limits[2][1] = right_lim / 2.0;
      element_size = -1.0 * Utilities::MPI::max(-element_size,
						mpi_communicator);
      break;
    }

    case 6: {

      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream file("../asb_plate.ucd");
      if (!file.is_open())
	{
	  std::cerr << "cannot find asb_plate.ucd" << std::endl;
	  exit(1);
	}

      grid_in.read_ucd(file);

      Point<dim> p1(Point1[0], Point1[1], Point1[2]);
      Point<dim> p2(Point2[0], Point2[1], Point2[2]);

      for (int d = 0; d < dim; d++) {
	limits[d][0] = p1[d] > p2[d] ? p2[d] : p1[d];
	limits[d][1] = p1[d] > p2[d] ? p1[d] : p2[d];
      }

      element_size = layer_elm_sz;

      break;
    }

    default:
      return -1;
    }

    //Point1 and Point2 are required to track the global deformation approximately, so they need to be specified here.
    assert(Point1 != Point2);
    return 0;
  }



  template<int dim>
    void Fdm<dim>::set_material_id()
    {
      typename Triangulation<dim>::cell_iterator
	cell = this->triangulation.begin(),
	endc = this->triangulation.end();
      for (; cell != endc; ++cell) {
	//if set extra material properties

	/* const Point<dim> cell_center = cell->center(); */
	/* if (fabs(cell_center[1])<=0.5 && cell_center[0] >=9) */
	/*   cell->set_material_id(1); */
	/* else */

	//end
	cell->set_material_id(0);
      }
    }

  template<int dim>
    void Fdm<dim>::set_b_c() {

    this->AllTraction = true;
    //b_c_split_id = b_c_codes.size()+1;

    for (unsigned int i = 0; i < b_c_codes.size(); i++)
      if (b_c_codes[i] < NMNDRCSPLIT) {
	this->AllTraction = false;
	break;
      }
    pcout << "fdm->Alltraction = " << AllTraction << "; ";

    //set b.c for vsvr
    typename Triangulation<dim>::cell_iterator
      cell = this->triangulation.begin(),
      endc = this->triangulation.end();
    for (; cell != endc; ++cell)
      for (unsigned int face  = 0;
	   face < GeometryInfo<dim>::faces_per_cell; ++face) {
	for (unsigned int d = 0; d < dim; d++) {
	  double negative_limit_this_dimension = limits[d][0];
	  double positive_limit_this_dimension = limits[d][1];

	  if (std::fabs(cell->face(face)->center()(d)
			- negative_limit_this_dimension) < 1e-12)
	    cell->face(face)->set_boundary_id(b_c_codes[d * 2]);
	  if (std::fabs(cell->face(face)->center()(d)
			- positive_limit_this_dimension) < 1e-12)
	    cell->face(face)->set_boundary_id(b_c_codes[d * 2 + 1]);
	}
      }

    //if with slit/inerhole, mark slit/hole with last 2nd id in bccodes
    if (mesher == 5 || mesher == 4 || mesher == 6) {
      assert( b_c_codes.size() > 6 );
      for (cell = this->triangulation.begin(); cell != endc; ++cell)
	for (unsigned int face  = 0;
	     face < GeometryInfo<dim>::faces_per_cell; ++face)
	  if (cell->face(face)->boundary_id() == 0) {
	    cell->face(face)->set_boundary_id(b_c_codes[6]);
	  }
    }

    //if set extra boundary indicator
    if (b_c_codes.size() > 7)
      for (cell = this->triangulation.begin(); cell != endc; ++cell)
	for (unsigned int face  = 0;
	     face < GeometryInfo<dim>::faces_per_cell; ++face) {
	  if (cell->face(face)->boundary_id() == 102)
	    {
	      Point<dim> P = cell->face(face)->center();
	      assert(fabs(P[0] - limits[0][1]) < 1e-12);
	      if (P[1] < 0)
		cell->face(face)->set_boundary_id(b_c_codes[7]);

	      /* if (fabs(P[1]-9)<1.5) */
	      /*   cell->face(face)->set_boundary_id(b_c_codes[8]); */
	    }

	  /* if (cell->face(face)->boundary_id() == 101) */
	  /*   { */
	  /*     Point<dim> P = cell->face(face)->center(); */
	  /*     assert(fabs(P[0]-limits[0][0]) < 1e-12); */
	  /*     if (fabs(P[1]-9)<1.5) */
	  /* 	cell->face(face)->set_boundary_id(b_c_codes[9]); */

	  /*     if (fabs(P[1]-5)<1.5) */
	  /* 	cell->face(face)->set_boundary_id(b_c_codes[10]); */
	  /*   } */
	}
  }

  template <int dim>
    void Fdm<dim>::setup_quadrature_point_history ()
    {
      unsigned int our_cells = 0;
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = triangulation.begin_active();
	   cell != triangulation.end(); ++cell)
	//if (cell->subdomain_id() == this_mpi_process)
	if (cell->is_locally_owned())
	  ++our_cells;

      triangulation.clear_user_data();
      {
	std::vector<PointHistory<dim> > tmp;
	tmp.swap (quadrature_point_history);
      }
      quadrature_point_history.resize (our_cells *
				       quadrature_formula.size());
      unsigned int history_index = 0;
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = triangulation.begin_active();
	   cell != triangulation.end(); ++cell)
	//if (cell->subdomain_id() == this_mpi_process)
	if (cell->is_locally_owned())
	  {
	    cell->set_user_pointer (&quadrature_point_history[history_index]);
	    history_index += quadrature_formula.size();
	  }
      Assert (history_index == quadrature_point_history.size(),
	      ExcInternalError());
    }

  template<int dim>
    void Fdm<dim>::set_standard_tensors() {
    int i, j, k, l;

    I = unit_symmetric_tensor<dim>();
    IxI = outer_product(I, I);
    II = identity_tensor<dim>();

    //    if (lambda <1e-15 && mu <1e-15)
    assert(mu.size() == 0);
    assert(lambda.size() == 0);
    assert(sL.size() == 0);

    for (unsigned int mid = 0; mid < matid; mid++) {
      Tensor<4, dim> elast_tensor;

      lambda.push_back(youngs[mid]*poisson[mid] /
		       (1 + poisson[mid]) / (1 - 2 * poisson[mid]));
      mu.push_back(youngs[mid] / 2.0 / (1 + poisson[mid]));

      Vp.push_back(sqrt((lambda[mid] + 2 * mu[mid]) / Rho[mid]));
      Vs.push_back(sqrt(mu[mid] / Rho[mid]));

      for ( i = 0; i < Ndim; i++)
	for ( j = 0; j < Ndim; j++)
	  for ( k = 0; k < Ndim; k++)
	    for ( l = 0; l < Ndim; l++)
	      elast_tensor[i][j][k][l] = lambda[mid] * I[i][j] * I[k][l]
		+ mu[mid] * (I[i][k] * I[j][l] + I[i][l] * I[j][k]);
      sL.push_back(elast_tensor);
    }

    assert(mu.size() == matid);
    assert(lambda.size() == matid);
    assert(sL.size() == matid);

    eE = 0;

    eE[0][1][2] = 1.0;
    eE[0][2][1] = -1.0;
    eE[1][2][0] = 1.0;
    eE[1][0][2] = -1.0;
    eE[2][0][1] = 1.0;
    eE[2][1][0] = -1.0;
  }
}
#endif
