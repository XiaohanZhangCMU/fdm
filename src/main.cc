#include "common.h"
#include "fdm.h"
#include "displ_equation.cc"
#include "plastic_equation.cc"
#include "heat_equation.cc"


int main(int argc, char * argv[]) {
  //default refinement grade
  unsigned int n_cycles = 2;
  std::string folder = "outputs";

  /* prepare output folder */
  try
    {
      if (argc >= 2)
	{
	  std::istringstream ist(argv[1]);
	  ist >> n_cycles;
	}
      if (argc >= 3)
	{
	  std::istringstream ist(argv[2]);
	  ist >> folder;
	}

      DIR* dir = opendir(folder.c_str());

      if (!dir) {
	std::string cmd = "mkdir -p " + folder;
	if (system(cmd.c_str()) != 0)
	  return -1;
      }
      if (chdir(folder.c_str()) != 0) {
	std::cerr <<
	  "cannot cd dir " << folder << std::endl;
	return -1;
      }
    }
  catch (std::exception &excpt)
    {
      std::cerr << "Exception on parsing arguments: " << std::endl
		<< excpt.what() << std::endl
		<< "Aborting!" << std::endl
		<< std::endl;
      return 1;
    }

  try
    {
      using namespace dealii;
      using namespace FDM;

      Utilities::MPI::
	MPI_InitFinalize mpi_initialization(argc, argv, 1);

      deallog.depth_console (0);
      {
	if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
	  // Make sure you change ncpyfiles
	  unsigned int ncpyfiles = 3;
	  std::string str_arr[] = { "src", "const_params.in", "inputs"};
	  std::vector<std::string> cpfiles(str_arr,str_arr+ncpyfiles);
	  assert(cpfiles.size() == ncpyfiles);
	  for(unsigned int i = 0;i<ncpyfiles;i++) 
	    if (system(("cp -r ../" + cpfiles[i] + " .").c_str()) != 0) {
	      if (system(("rm -rf ./"+cpfiles[i]).c_str()) != 0) {
		std::cerr << "cannot copy " << cpfiles[i]<<std::endl;
		return -1;
	      }
	      else
		system(("cp -r ../" + cpfiles[i] + " .").c_str()); 
	    }	
	}
	  
#if 1
	Fdm<Ndim> *fdm = new Fdm<Ndim>(n_cycles, std::cout); 
#else
	std::ofstream out("log.out");
	Fdm<Ndim> *fdm = new Fdm<Ndim>(n_cycles, out);		       
#endif

	DisplEquation<Ndim>* ueqn    = new DisplEquation<Ndim>(fdm);
	PlasticEquation<Ndim-1>* peqn  = new PlasticEquation<Ndim-1>(fdm);
	//HeatEquation<Ndim>* heqn        = new HeatEquation<Ndim>(fdm);

	//All equations have access to the other eqns local information
	ueqn->set_coupled(peqn);
	peqn->set_coupled(ueqn);
	//heqn->set_coupled(peqn, veqn);
			
	ueqn->setup_system();
	peqn->setup_system();
	//heqn->setup_system();

	fdm->tstep = -1;

	while (fdm->clock <= fdm->duration) {

	  fdm->pcout << "\n\n..........SOLUTION OF TSTEP.........."
		     << fdm->tstep
		     << "\nCLOCK = " << fdm->clock << "; DCLOCK = "
		     << fdm->dclock
		     << "; Cutback since " << fdm->cutback_tstep <<"\n"
	             <<"; normRF = "<<fdm->RF.linfty_norm()
		     <<"; normDisplacement= "<<"??"
		     << "\n\n" << std::endl;
	    
	  if (fdm->tstep == -1) {

	    peqn->drive ( INITL ); 
	    //heqn->drive ( INITL ); 
	    ueqn->drive ( INITL ); 
	    fdm->tstep ++;
	    continue;

	  } else if (fdm->tstep == 0) { //ECDD solve on x(t=0) configuration.

	    ueqn->drive ( SOLVE );
	    ueqn->drive ( SSCRV );

	  } else {

	    peqn->drive ( SOLVE );
	    ueqn->drive ( SOLVE );
	    //heqn->drive ( SOLVE );

	  }	  
	  
	  if (fdm->cutback_tstep>0 && fdm->tstep % 1 == 0) {
	    ueqn->output_results();
	    //heqn->output_results();
	  }
	  
	  fdm->clock = fdm->clock + fdm->dclock;
	  fdm->tstep ++;
	  fdm->cutback_tstep ++;
	}

	peqn->drive ( CLEAR );	ueqn->drive ( CLEAR ); //heqn->drive ( CLEAR );
			
	if (ueqn != NULL) delete ueqn;
	fdm->pcout << "ueqn removed" << std::endl;
	if (peqn != NULL) delete peqn;
	fdm->pcout << "peqn removed" << std::endl;
	//if (heqn != NULL) delete heqn;	
	//fdm->pcout<<"heqn removed"<<std::endl;
	fdm->pcout << "About to exit fdm.run() normally...\nBut I am not sure if fdm is freed :)" << std::endl;
	delete fdm;
      }
    }

  catch (std::exception &excpt)
    {
      std::cerr << "Exception on processing: " << std::endl
		<< excpt.what() << std::endl
		<< "Aborting!" << std::endl
		<< std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << "Exception on processing: " << std::endl
		<< "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< std::endl;
      return 1;
    }


  return 0;
} 
