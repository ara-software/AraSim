////////////////////////////////////////////////////////////////////////////////////////////////
//class Counting:
////////////////////////////////////////////////////////////////////////////////////////////////

class Vector;
class Event;

class Counting {
public:

  Counting();
  ~Counting();

  void initializeEachRun();
  
  int npass[2]; // count events that pass
  int npassestrigger[2]; // incremented if passes trigger
  int nchanceinhell2[2]; // based on direction of ray and thickness in Cerenkov cone,
  // signal has a chance to pass after accounting for 
  // angle, ice attenuation and 1/r 
  int nviewanglecut[2];
  int nchanceinhell[2]; // based on depth, 
  // signal has a chance to pass after accounting for ice attenuation and 1/r
  int nchanceinhell_1overr[2]; // after 1/r
// signal has chance of passing
  int nchanceinhell_fresnel[2]; // after fresnel coefficients
  // signal has chance of passing

  int nconverges[2]; // ray tracing converges to within 10m at ice surface

  int nacceptablerf[2]; // ray leaves where there is ice, and not where there is water
  int nraywithincontinent1[2]; // reality check:  exiting ray is within 30 degrees of south pole
  int nraywithincontinent2[2]; // same, after next iteration.
  int nraypointsup1[2]; // ray from exit point to balloon does not intersect earth
  int nnottoosmall[2]; // based on neutrino position, 
  int nraypointsup2[2]; // same, after next iteration to get refracted ray.
  int nviewangle_lt_90[2]; // viewing angle lt 90
  int ngoodfracs[2]; // for debugging
  int nbadfracs[2]; // for debugging
  int nnottir[2]; // not totally internally reflected
  int nentersice[2]; // Reality check:  place where neutrino enters ice (from below) within 30 deg of south pole
  int nabsorbed[2]; //  Event has more than 0.001 probability of surviving trip through earth 
  int noway[2]; // no way the event will see any ice give its earth entrance point and its direction
  int wheredoesitleave_err[2]; // wheredoesitleave gives error
  int neverseesice[2];  // determined that the neutrino path never sees ice
  int iceinteraction[2];  // there is an interaction in the ice
  int inhorizon[2];  // there is an interaction in the ice
  int wheredoesitenterice_err[2];  // there is an interaction in the ice
  int toohigh[2];  // there is an interaction in the ice
  int toolow[2];  // there is an interaction in the ice


// variables for counting neutrinos and reporting results.
  int nnu_e;  //counting the number of e,mu,tau neutrinos
  int nnu_mu;  
  int nnu_tau;

  static const int NCOSTHETA=180;
  static const int NPHI=180;
  static constexpr double COSTHETAMAX=1.0;
  static constexpr double COSTHETAMIN=0.0;
  static constexpr double PHIMAX=2*3.14159;
  static constexpr double PHIMIN=0.;
 double weights_rin[NCOSTHETA][NPHI];

  void findCosthetaPhiBins(Position r_in,int &icostheta,int &iphi);
  void IncrementWeights_r_in(Position r_in,double weight) ;


  // added for increment weight from icemc
  static const int NBINS=10; // keep track of the number of events found, binned
  // by weights
  static constexpr double MIN_LOGWEIGHT=-3;
  //static const double MAX_LOGWEIGHT=-1;
  static constexpr double MAX_LOGWEIGHT=0;

  double eventsfound_binned[NBINS]; 
  double eventsfound_binned_e[NBINS];
  double eventsfound_binned_mu[NBINS];
  double eventsfound_binned_tau[NBINS];
  double sum[3]; // sum of weight for events found for 3 flavors  

  static void findErrorOnSumWeights(double *eventsfound_binned,double &error_plus,double &error_minus);
  //void incrementEventsFound(double weight,Interaction *interaction1);
  void incrementEventsFound(double weight,Event *event);
  static int findWeightBin(double logweight);
  //


protected:

};

