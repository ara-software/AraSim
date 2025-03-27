////////////////////////////////////////////////////////////////////////////////////////////////
//class Primaries:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PRIMARIES_H_
#define PRIMARIES_H_

#include "TRandom3.h" 
#include <iostream>
#include "TObject.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TF3.h"
#include "TH2D.h"
#include "Vector.h"
#include <vector>
#include "Position.h"

class IceModel;
class Counting;
class Settings;
class Detector;
class Ray;
class Primaries;
class Spectra;
class Secondaries;
class Signal;
class RaySolver;
class Report;

using namespace std;

class Y {

	//Code from Connolly Calc 2011.
  private:
    TF1* ffrac; // fraction of events in the low y region.
    TF1* fC1_high[2][2]; // parameterization of C1 in high y region

    TF1 *fC1_low; // parameterization of C1 in the low y region.

    TF1 *fC2; // parameterization of C2 in the low y region.

    TF3* fy0_low; // for picking y0. 
    TF2* fy0_high;  // for picking y0.
    TRandom3 Rand3;

  public:
    Y();
    ~Y();

    double pickY(int NU,int CURRENT,double e); // pick an inelasticity
    // NU=0: nubar, NU=1: nu
    // CURRENT=0: CC, CURRENT-1: NC

    static const double miny_low;//2.E-5
    static const double maxy_low;
    static const double miny_high;
    static const double maxy_high;


  ClassDef(Y,1);

};//Y

class Primaries {

  private:
    TRandom3 Rand3;

    TH2D *m_hsigma;
    TCanvas *m_csigma;
    Y *m_myY;
    int run_old_code;

  public:

    Primaries();//constructor //default
    ~Primaries();//destructor //default

    double ymin_low, ymax_low, ymin_high, ymax_high;
    double A_low[4];//same for any nu_nubar and current type.
    double A0_high[2][2];
    double A1_high[2][2];
    double A2_high[2][2];
    double A3_high[2][2];
    double b0, b1; 

    TF1* m_fy[2][2];
    TF1* m_fsigma[2][2];

    TF1* m_fsigma_upper[2][2];
    TF1* m_fsigma_lower[2][2];

    double c0[2][2];
    double c1[2][2];
    double c2[2][2];
    double c3[2][2];
    double c4[2][2];

    double c0_upper[2][2];
    double c1_upper[2][2];
    double c2_upper[2][2];
    double c3_upper[2][2];
    double c4_upper[2][2];

    double c0_lower[2][2];
    double c1_lower[2][2];
    double c2_lower[2][2];
    double c3_lower[2][2];
    double c4_lower[2][2];

    static const int NSIGMAS=2;// number of possible cross section models
    // 0=Gandhi et al.
    // 1=Connolly et al. 2011
    double mine[NSIGMAS];// minimum energy for cross section parametrizations, in eV
    double maxe[NSIGMAS]; //max

    double col1[900]; // array for air density
    void GetAir(double *col1);
    double GetThisAirColumn(Settings* settings1, Position r_in,Vector nnu,Position posnu, double& cosalpha,double& mytheta, double& cosbeta0,double& mybeta);

    // GetAnyDirection for selecting nnu
    Vector GetAnyDirection(double phi, double d_phi, int nnu_this_phi);   // get random direction for nnu. Added this function instead of removing PickAnyDirection in Interaction to prevent confliction
    Vector GetAnyDirection();
    Vector GetThatDirection( double theta, double d_theta);   // get direction for nnu near theta angle with d_theta variation. Added this function instead of removing PickAnyDirection in Interaction to prevent confliction
    Vector GetThatDirection( double theta, double d_theta, double phi, double d_phi, int nnu_this_phi);
    double costheta_nutraject; //theta of nnu with earth center to balloon as z axis 
    double phi_nutraject; //phi of nnu with earth center to balloon as z axis

    int GetSigma(double pnu,double& sigma,double &len_int_kgm2,Settings *settings1,int nu_nubar,int currentint);//not static
    int GetSigma(double pnu,double& sigma,double &len_int_kgm2,Settings *settings1,int nu_nubar,int currentint, double &len_int_kgm2_total);// len_int_kgm2_total is interaction length with (CC + NC cross section)
    double Gety(Settings *settings1,double pnu,int nu_nubar,int currentint);
    double Getyweight(double pnu,double y,int nu_nubar,int currentint);
    string GetCurrent();
    string GetCurrent(Settings *settings1);
    string GetNuFlavor();
    string GetNuFlavor(Settings *settings1);

    int GetNuNuBar( string nuflavor );
    int GetNuNuBar( string nuflavor, Settings *settings1);


    int IsCalpulser;

  ClassDef(Primaries,1);

};//Primaries

class Interaction  {

  private:


  public:

    // default constructor
    Interaction();
    // constructor for setting posnu, y, emfrac, hadfrac, vmmhz1m at cherenkov angle, etc
    Interaction (double pnu, string nuflavor, int nu_nubar, int &n_interactions, IceModel *antarctica, Detector *detector, Settings *settings1, Primaries *primary1, Signal *signal, Secondaries *sec1 );
    // constructor for setting posnu, y, emfrac, hadfrac, vmmhz1m at cherenkov angle, etc
    Interaction (IceModel *antarctica, Detector *detector, Settings *settings1, Primaries *primary1, Signal *signal, Secondaries *sec1 );
    // constructor for setting posnu, y, emfrac, hadfrac, vmmhz1m at cherenkov angle, etc 
    Interaction (Settings *settings1, Detector *detector, IceModel *antarctica, Primaries *primary1, Signal *signal); 

    ~Interaction();

    void Initialize (); // initialize values for Interaction class.

    // variables for GetSignal
    //
    string taudecay;
    double emfrac, hadfrac;
    double elast_y;
    int ray_sol_cnt;   // counting number of solutions from Solve_Ray
    vector <double> vmmhz1m;
    vector <double> vmmhz1m_em;
    vector <double> d_theta_em;
    vector <double> d_theta_had;
    double vmmhz1m_tmp;
    // end variables for GetSignal

    // values related to index of refraction
    double indexN;     // index of refraction at posnu depth
    double changle;    // cherenkov angle at posnu depth

    void PickAnyDirection();

    int noway;
    int wheredoesitleave_err;
    int neverseesice;
    int wheredoesitenterice_err;
    int toohigh;
    int toolow;

    int pickposnu;    // : 0 for fail picking posnu, nnu  : 1 sucess picking posnu, nnu
      
    int ray_solver_toggle;    // : 0 no solution, : 1 solution exists

    int PickUnbiased (IceModel *antarctica);
    int WhereDoesItLeave( const Position &posnu, const Vector &ntemp, IceModel *antarctica, Position &r_out);
    int WhereDoesItEnterIce ( const Position &posnu, const Vector &nnu, double stepsize, Position &r_enterice_output, IceModel *antarctica);
    int WhereDoesItExitIce ( const Position &posnu, const Vector &nnu, double stepsize, Position &r_enterice_output, IceModel *antarctica);
    int WhereDoesItExitIceForward ( const Position &posnu, const Vector &nnu, double stepsize, Position &r_enterice_output, IceModel *antarctica);
    void FlattoEarth ( IceModel *antarctica, double X, double Y, double D);
    void FlattoEarth_AboveIce ( IceModel *antarctica, double X, double Y, double D, double height);
    void FlattoEarth_Near_Surface ( IceModel *antarctica, double X, double Y, double D, double max_depth);
    void FlattoEarth_Spherical ( IceModel *antarctica, double X, double Y, double Z);

    void PosNuFromAntennaCenter (Detector *detector); ///< re-calculate Neutrino position (x, y, z, r, theta, phi) from antenna center point of view.

    void PickNear_Cylinder (IceModel *antarctica, Detector *detector, Settings *settings1, double energy = 0);
    double PickNear_Sphere (IceModel *antarctica, Detector *detector, Settings *settings1);

    void PickExact(IceModel *antarctica, Detector *detector, Settings *settings1, double R, double Theta, double Phi);
    void PickNear_Cylinder_AboveIce (IceModel *antarctica, Detector *detector, Settings *settings1);

    bool Does_Interact(double x, double y, double z,
    double theta, double phi, double r,
    double &newx, double &newy, double &newz, double& l);

    void PickExactGlobal(IceModel *antarctica, Detector *detector, Settings *settings1, double thisLat, double thisLong, double thisDepth);

    double pathlength_inice;

    Vector nnu;  // direction of neutrino (+z in south pole direction)
    Vector cone_axis;  // cone axis which can be same with neutrino direction or not

    double costheta_nutraject; //theta of nnu with earth center to balloon as z axis 
    double phi_nutraject; //phi of nnu with earth center to balloon as z axis

    Position r_in; // position where neutrino enters the earth
    Position r_enterice; // position where neutrino enters the ice
    Position nuexit; // place where neutrino would have left the earth
    Position nuexitice; // place where neutrino would have left the ice
    double chord;  // chord in m from earth entrance to rock-ice boundary
    double logchord; // log_10 of chord length earth entrance to where it enters ice
    double weight_bestcase; // what weight1 would be if whole earth had density of crust - for quick and dirty calculation of best case scenario
    int sigma_err;    // 0 if GetSigma is un-successful (!), otherwise 1
    double sigma;

    // input information for Getchord
    double len_int_kgm2;
    double len_int_kgm2_total;
    double weight;
    double probability; // weight * probability to interact inside the ice
    double nearthlayers;
    double myair;
    double total_kgm2;
    int crust_entered;
    int mantle_entered;
    int core_entered;
    // input information for Getchord

    void setCurrent(Primaries *primary1, Settings *settings1);
    Position posnu;

    Position posnu_from_antcen; ///< Nu position (x,y,z,r,theta,phi) from antenna center.

    Position posnu_down;
    string  current;                    //  CC or NC?
    int currentint;                 // Ditto - Stephen

    // values required for t-domain signal mode
    vector <double> shower_depth_m; // shower depth array in meters
    vector <double> shower_Q_profile; // shower charge excess array
    double LQ; // integrated charge excess
    vector <double> EM_shower_depth_m; // EM shower depth array in meters
    vector <double> EM_shower_Q_profile; // EM shower charge excess array
    double EM_LQ; // integrated charge excess

    vector <double> HAD_shower_depth_m;
    vector <double> HAD_shower_Q_profile;
    double HAD_LQ;
    double pnuenergy;

    int primary_shower; // 0 or EM, 1 for HAD. set by emfrac, hadfrac

    void clear_useless(Settings *settings1);
  
  ClassDef(Interaction,3);


};//Interaction

#endif
