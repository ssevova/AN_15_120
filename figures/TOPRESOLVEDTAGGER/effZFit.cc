#include "CEffZFitter.hh"
//#include "CEffPlotter.hh"
#include "KStyle.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings and constants
  //==============================================================================================================

  // handle input arguments
  //MT

  const std::string infname  = argv[1];         // input ROOT file of probes
  const std::string outdir   = argv[2];         // output directory
  const std::string temfname1 = argv[3];        // ROOT file for generating MC-based signal templates
  const std::string temfname2 = argv[4];        // ROOT file for generating MC-based bg1 (Wrong comb tt1l) templates
  const std::string temfname3 = argv[5];        // ROOT file for generating MC-based (nontt1l bg) templates
  const float         topmvacut  = atof(argv[6]); //topmvacut to split inclusive sample in pass and fail  
  const bool         maketemplates  = atoi(argv[7]); //decide whether to make templates or not (they may already exist)  
  // other settings
  const double       massLo    = 0;//60;
  const double       massHi    = 800;//120;
  const double       fitMassLo = massLo;
  const double       fitMassHi = massHi;
  const unsigned int runNumLo  = 0;
  const unsigned int runNumHi  = 999999;
  
  std::cout << std::endl;
  std::cout << " <> Processing probes file: " << infname << std::endl;
  std::cout << " outdir: " << outdir << std::endl;
  //std::cout << " doPU: " << doPU << std::endl;
  std::cout << " outdir: " << outdir << std::endl;
  std::cout << " temfname1, temfname2, temfname3: " << temfname1 <<" "<<temfname2<<" "<<temfname3<< std::endl;
  std::cout << "topmvacut "<<topmvacut<<" argv[6] "<<argv[6]<<std::endl;
  std::cout << "maketemplates "<<maketemplates<<" argv[7] "<<argv[7]<<std::endl;
  std::cout << std::endl;

  KStyle();  
  
  CEffZFitter fitter;

 fitter.initialize(infname, outdir, temfname1,temfname2,temfname3,
                    massLo, massHi, fitMassLo, fitMassHi, 
		   runNumLo, runNumHi,
		   topmvacut, maketemplates);
  fitter.computeEff();
  
  std::cout << " <> Output saved in " << outdir << std::endl;






  return 0;
}
