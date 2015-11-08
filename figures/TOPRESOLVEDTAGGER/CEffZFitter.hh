#ifndef CEFFZFITTER_HH
#define CEFFZFITTER_HH

//================================================================================================
//
// Signal Extraction
//-------------------
//  0: probe counting
//  1: Breit-Wigner convolved with Crystal Ball function
//  2: MC template convolved with Gaussian
//  3: Phil's Crystal Ball based "Voigtian" shape
//  4: Unbinned MC data convolved with Gaussian
//
// Background Model
//------------------
//  0: no background
//  1: exponential model
//  2: erfc*exp model
//  3: double exponential model
//  4: linear*exp model
//  5: quadratic*exp model
//
//________________________________________________________________________________________________

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <TLorentzVector.h>

class TTree;
class TCanvas;
class TGraphAsymmErrors;
class TH1D;
class TH2D;
class RooFitResult;

class CEffZFitter
{
public:
  CEffZFitter();
  ~CEffZFitter();
  
  
 
void initialize(const std::string infname, const std::string outdir, const std::string temfname1, const std::string temfname2, const std::string temfname3,
		const double massLo, const double massHi, const double fitMassLo, const double fitMassHi, 
		const unsigned int runNumLo, const unsigned int runNumHi,
		const double TOPMVACUT, const bool MAKETEMPLATES);
  void computeEff();
  
    
protected:

  void makeBinnedTemplates(const std::string temfname, std::string name); //MT
  
  
  
  
  void performFit(double &resEff, double &resErrl, double &resErrh,
		  TTree *passTree, TTree *failTree);


  
  float calcmass(TLorentzVector *obj1, TLorentzVector *obj2, TLorentzVector *obj3);  
  ///// data members /////
  
  bool fIsInitialized;
  bool maketemplates;  
  
  int fSigPass, fBkgPass, fSigFail, fBkgFail;
  
  double fMassLo, fMassHi;        // signal extraction mass window  
  double fFitMassLo, fFitMassHi;  // fit mass window
  
  double topmvacut;


  
  // output directory for results
  std::string fOutputDir;
      
  TTree *fPassTree, *fFailTree;
  
  
};

#endif
