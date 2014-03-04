#ifndef SPECTRA_H_
#define SPECTRA_H_

#include "TObject.h"
#include "TSpline.h"
#include <string>
#include "TRandom3.h"


using namespace std;

class Spectra {

private:
  TRandom3 Rand3;
  double maxflux;   // max flux value
//  static const int NSPECTRA_MAX=300;  // why need this??
  static const int E_bin_max = 50;
  int E_bin;   // initialize # of energy bins (max = 50)

  void GetFlux(string filename);    // read neutrino flux EdNdEdAdt (in GeV) from filename file

//--------------------------------------------------
//   TGraph *gEdNdEdAdt;   //graph for EdNdEdAdt flux
//   TGraph *gE2dNdEdAdt;  //graph for E2dNdEdAdt flux
// 
//   TSpline3 *sEdNdEdAdt; //spline of EdNdEdAdt
//   TSpline3 *sE2dNdEdAdt;    //spline of E2dNdEdAdt
//-------------------------------------------------- 

  //int EXPONENT_model; // set flux model
  double EXPONENT_model; // set flux model

  double pnu_EXPONENT;  // if mono energy from EXPONENT, pnu_EXPONENT = log10(pnu), constant pnu for all events.

public:  


  double energy[E_bin_max]; // energies that correspond to the fluxes in the previous array  
  double EdNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E*dN/dE/dA/dt
  double E2dNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E^2*dN/dE/dA/dt
  
  Spectra();    // default constructor
  //Spectra(int EXPONENT); // constructor  
  Spectra(double EXPONENT); // constructor  
  ~Spectra();   // destructor
  
  double GetNuEnergy_bin(); // get the neutrino energy which follows neutrino flux. (bin step)

  double GetNuEnergy(); // get the neutrino energy which follows neutrino flux. (interpolated flux)


  double SimpleLinearInterpolation_value(int n1, double *x1, double *y1, double x2 );    // reads n1 array values x1, y1 and do simple linear interpolated value at x2 and return it


//--------------------------------------------------
//   TGraph *GetGEdNdEdAdt();
//   TGraph *GetGE2dNdEdAdt();
// 
//   TSpline3 *GetSEdNdEdAdt();
//   TSpline3 *GetSE2dNdEdAdt();
//-------------------------------------------------- 

  double *Getenergy();
  double *GetEdNdEdAdt();
  double *GetE2dNdEdAdt();
  double GetEdNdEdAdt(double E_val);    // return flux value from TSpline
  double GetE2dNdEdAdt(double E_val);   // return flux value from TSpline

  double Getmaxflux();
  
  int GetE_bin();   // return energy bin number

  int IsSpectrum(); // return 1 or 0 depend on EXPONENT value
  int IsMonoenergetic();    // return 1 or 0 depend of EXPONENT value

  // destructor


  ClassDef(Spectra,1);

}; //class Spectra

#endif
