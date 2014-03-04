#ifndef SIGNAL_H_
#define SIGNAL_H_
////////////////////////////////////////////////////////////////////////////////////////////////
//class Signal:
////////////////////////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>
#include <vector>

class Settings;
class Event;

using std::cout;

class Signal {

protected:
  // properties of ice
  double x0ice; 
  //double X0ICE; // radiation length of ice (meters)
  double ecice;                     // critical energy in ice (MeV)
  //const static double ECICE; // critical energy in ice (MeV)
  double nice;                      // index of refraction of ice
  double nfirn;                   // index of refraction at the very surface - Peter
  double invnfirn; 
  double invnice;
  double rhoice;                     // density of ice (kg/m**3)
  double kelvins_ice;          // temperature in Kelvin (ice+system)
  double changle_ice;
  double aex_ice;  //efficiency for producing charge asymmetry relative to ice.  1 by definition
  //double n_depth;  // index of refraction at the interaction depth
  double alphaice; // exponent that goes into cutting off the spectrum at high frequencies
  double rm_ice; // moliere radius, in g/cm^2
  double ke_ice; // const staticant in jaime's parameterization, in V/cm/MHz
  double kl_ice; //const staticant in jaime's parameterization
  double kdelta_ice; // const staticant in jaime's parameterization
  double kr_ice; // const staticant in jaime's parameterization
  double betaice; // exponent, in jaime's parameterization
  double nu0_modified; // nu_0 modified for a specific medium
  double nu_r;// parameter for signal parameterization
  int WHICHPARAMETERIZATION;
  double vmmhz1m_reference; // reference value for V/m/MHz at f=1 MHz and pnu=10^18 eV
  double freq_reference; // reference frequency in MHz
  double pnu_reference; // reference energy in eV

    double CalpulserVmMHz1m[60];



  double KR_MEDIUM; // constant in jaime's parameterization
  double RM_MEDIUM; // moliere radius, in g/cm^2
double KL_MEDIUM; //constant in jaime's parameterization
  
  double KE_MEDIUM; // constant in jaime's parameterization, in V/cm/MHz
  double ECMEDIUM; // radiation length of medium
  
  
  double ALPHAMEDIUM;// exponent that goes into cutting off the spectrum at high frequencies
  double AEXMEDIUM; // efficiency for making charge asymmetry
  double KDELTA_MEDIUM; // constant in jaime's parameterization
double BETAMEDIUM; // exponent, in jaime's parameterization


  double JAIME_FACTOR; // factor to multiply Jaime's parameterization for error analysis
  int MEDIUM;
  int LPM;

  static const double RHOSALT;                 // density of salt (kg/m**3)
  
  

  static const double RM_ICE; // moliere radius, in g/cm^2
  static const double RM_SALT; // moliere radius, in g/cm^2
  static const double KR_SALT; // constant in jaime's parameterization
  static const double KR_ICE; // constant in jaime's parameterization
  //const double X0SALT=0.1024;                // radiation length of salt (meters)
  static const double X0SALT;                // radiation length of salt (meters)
  //const double ECSALT=40.;                   // critical energy in salt (MeV)
  static const double ECSALT;                   // critical energy in salt (MeV)
  static const double X0ICE; 
  // //const double X0ICE=0.392; // radiation length of ice (meters)
  static const double ECICE;                     // critical energy in ice (MeV)
  // //const double ECICE=73.0; // critical energy in ice (MeV)
static const double AEX_ICE;  //efficiency for producing charge asymmetry relative to ice.  1 by definition
  


  static const double ALPHAICE; // exponent that goes into cutting off the spectrum at high frequencies
  static const double AEX_SALT;  // efficiency for producing charge asymmetry relative to ice
  static const double ALPHASALT; // exponent that goes into cutting off the spectrum at high frequencies
  static const double KE_SALT; // constant in jaime's parameterization, in V/cm/MHz
  static const double KL_SALT; //constant in jaime's parameterization
  static const double KDELTA_SALT; // constant in jaime's parameterization
  static const double KE_ICE; // constant in jaime's parameterization, in V/cm/MHz
  static const double KL_ICE; //constant in jaime's parameterization
  static const double KDELTA_ICE; // constant in jaime's parameterization
  
  static const double KELVINS_ICE;          // temperature in Kelvin (ice+system)
  static const double KELVINS_SALT;            // temperature in salt (350) + receiver temp (150)
  static const double BETAICE; // exponent, in jaime's parameterization
  // double NU0_MODIFIED=0.; // nu_0 modified for a specific medium
  // double NU_R;// parameter for signal parameterization
  static const double BETASALT; // exponent, in jaime's parameterization 
  



public:
  Signal();
  Signal(Settings *settings1);
  ~Signal();

  void TaperVmMHz(double viewangle,double deltheta_em,double deltheta_had,double emfrac,double hadfrac,
		double& vmmhz1m,
	       double& vmmhz_em); // returns 1 if viewangle-changle<20*width for both em and had showers
  double GetVmMHz1m(double pnu,double freq); // constructor
    double GetVmMHz1mCalPulser(int bin);
    void ReadCalPulserSpectrum();
  void GetVmMHz(double vmmhz_max,double vmmhz1m_max,double pnu,double *freq,double notch_min,double notch_max,double *vmmhz,int nfreq);
  void Initialize();
  void Initialize(Settings *settings1);
  
  void SetParameterization(int whichparameterization);
  double vmmhz1m_max; // V/m/MHz at 1m
  int GetLPM();  // lpm
double GetELPM();  // elpm
  void GetSpread(double pnu,
		 double emfrac,
		 double hadfrac,
		 double freq,
		 //		 double n_depth,
		 // double X0DEPTH,
		 
		 double& deltheta_em_max,
		 double& deltheta_had_max);
 

  //
  // functions for parameterized t-domain signal calculation
  //
  void GetShowerProfile(double E_shower, // energy of shower
                            int shower_mode, // 0 : EM shower, 1 : HAD shower 
                            double shower_step_m, // shower step in meters
                            int param_model, // 0 : Jaime's fit, 1 : Carl's fit
                            //int q_excess_model, // 0 : const 25% (only option now)
                            std::vector <double> &depth, // shower depth array, output
                            std::vector <double> &Q_shower, // shower charge excess profile, output
                            double &LQ // integrated Q_shower
        );

  double GaisserHillas(double x_in, double *par);
  double Greisen(double x_in, double *par);


  //void GetVm_FarField_Tarray( Event *event, Settings *settings1, double viewangle, double atten_factor, int outbin, double *Tarray, double *Earray );
  void GetVm_FarField_Tarray( Event *event, Settings *settings1, double viewangle, double atten_factor, int outbin, double *Tarray, double *Earray, int &skip_bins );

  double Param_RE_Tterm(double Tterm, double *par);
  double Param_RE_Tterm_approx(double Tterm, double *par); // use approx expansion if possible




double X0MEDIUM; // radiation length of medium
double KELVINS;  // temperature of medium + system
static const double RHOICE;                     // density of ice (kg/m**3)
  static const double RHOAIR;          // density of air (kg/m**3)
  static const double RHOH20;          // density of water (kg/m**3) 
  double N_DEPTH;  // index of refraction at the interaction depth
  double RHO_DEPTH;  // density at the interaction depth
  double X0_DEPTH;  // density at the interaction depth
  double NMEDIUM_RECEIVER; // index of refraction at receiver
double changle; // cherenkov angle  
  double RHOMEDIUM; // density of medium
  double logscalefactor_taper;
  static const double N_AIR;             // index of refr for air
  static const double NICE;                 // index of refraction of ice
  static const double NSALT;                 // index of refraction of salt
  static const double CHANGLE_ICE; // cherenkov angle in ice
  void SetMedium(int medium) {
    MEDIUM=medium;
    if (MEDIUM!=0) {
      std::cout << "Medium is " << MEDIUM << "\n";
      std::cout << "Non-default setting:  Not ice!\n";
    }
    InitializeMedium();
  }
  static const double VIEWANGLE_CUT;
  void InitializeMedium();

  void SetNMediumReceiver(double nmedium_receiver) {
    NMEDIUM_RECEIVER=nmedium_receiver;
  }
void SetLPM(double lpm) {
  LPM=(int)lpm;
  }
  void SetKelvins(double kelvins) {
    KELVINS=kelvins;
  }
  void SetBetaMedium(double betamedium) {
    BETAMEDIUM=betamedium;
  }


  void SetRhoMedium(double rhomedium) {
    RHOMEDIUM=rhomedium;
  }
  void SetKrMedium(double kr_medium) {
    KR_MEDIUM=kr_medium;
  }
  void SetKlMedium(double kl_medium) {
    KL_MEDIUM=kl_medium;
  }

  void SetRmMedium(double rm_medium) {
    RM_MEDIUM=rm_medium;
  }
  void SetNDepth(double n_depth) {
    N_DEPTH=n_depth;
    SetChangle(acos(1/N_DEPTH));

    SetrhoDepth((N_DEPTH-1.)/0.86*1000.);
    //SetrhoDepth(RHOICE);

    SetX0Depth(X0MEDIUM); // no dependence on rho


  }
  void SetX0Depth(double x0_depth) {
    X0_DEPTH=x0_depth;
  }
  void SetrhoDepth(double rho_depth) {
    RHO_DEPTH=rho_depth;
  }
  void SetKeMedium(double ke_medium) {
    KE_MEDIUM=ke_medium;
  }
  void SetEcMedium(double ecmedium) {
    ECMEDIUM=ecmedium;
  }
  void SetX0Medium(double x0medium) {
    X0MEDIUM=x0medium;
  }
  void SetChangle(double thischangle) {
    changle=thischangle;
  }
  void SetAlphaMedium(double alphamedium) {
    ALPHAMEDIUM=alphamedium;
  }
  void SetAexMedium(double aexmedium) {
    AEXMEDIUM=aexmedium;
  }
  void SetKdelta_Medium(double kdelta_medium) {
    KDELTA_MEDIUM=kdelta_medium;
  }
  void SetJaime_Factor(double jaime_factor) {
    JAIME_FACTOR=jaime_factor;
  if (JAIME_FACTOR!=1)
    std::cout << "Non-default setting:  JAIME_FACTOR= " << JAIME_FACTOR << "\n";
  }


}; //class Position

#endif
