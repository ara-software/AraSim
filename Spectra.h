#ifndef SPECTRA_H_
#define SPECTRA_H_

#include "TObject.h"
#include "TSpline.h"
#include <string>
#include "TRandom3.h"

#include<vector>

using namespace std;

class Settings;


class Spectra 
{

  private:
    TRandom3 Rand3;
    double maxflux;   // max flux value
    int E_bin;   // initialize # of energy bins (max = 50)

    void GetFlux(string filename);    // read neutrino flux EdNdEdAdt (in GeV) from filename file

    double EXPONENT_model; // set flux model
    double EXPONENT_min;
    double EXPONENT_max;

    double pnu_EXPONENT;  // if mono energy from EXPONENT, pnu_EXPONENT = log10(pnu), constant pnu for all events.

  public:  


    vector<double> energy; // energies that correspond to the fluxes in the previous array  
    vector<double> EdNdEdAdt; //flux of incident neutrinos vs. energy E*dN/dE/dA/dt
    vector<double> E2dNdEdAdt; //flux of incident neutrinos vs. energy E^2*dN/dE/dA/dt
    vector<double> lgEdNdEdAdt; 
    vector<double> lgE2dNdEdAdt; 
   
    Spectra();    // default constructor
    Spectra(Settings *settings1); // constructor  
    ~Spectra();   // destructor
    
    double GetNuEnergy(); // get the neutrino energy which follows neutrino flux. (interpolated flux)

    double SimpleLinearInterpolation_value(int n1, double *x1, double *y1, double x2 );    // reads n1 array values x1, y1 and do simple linear interpolated value at x2 and return it

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
