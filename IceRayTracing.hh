/*
This is the IceRayTracing namespace. Author: Uzair Latif 
released under GPL3.
*/
#ifndef IRT_HEAD
#define IRT_HEAD

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <chrono>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <sys/time.h>
#include <gsl/gsl_integration.h>

using namespace std;

namespace IceRayTracing{

  /********Stuff for Interpolation**********/
  static vector <double> GridPositionX;
  static vector <double> GridPositionZ;
  static vector <double> GridZValue[4];

  static double GridStepSizeX_O=0.2;
  static double GridStepSizeZ_O=0.2;
  static double GridWidthX=20;
  static double GridWidthZ=20;

  static int GridPoints=100;////just set a non-zeronumber for now
  static int TotalStepsX_O=100;////just set a non-zeronumber for now
  static int TotalStepsZ_O=100;////just set a non-zeronumber for now
  static double GridStartX=40;////just set a non-zeronumber for now
  static double GridStopX=60;////just set a non-zeronumber for now
  static double GridStartZ=-20;////just set a non-zeronumber for now
  static double GridStopZ=0;
  
  /* Set the value of pi */
  static constexpr double pi=3.14159265359;
  /* Set the value of the speed of light in m/s */ 
  static constexpr double c_light_ms=299792458;
  /* Set the value of the asymptotic parameter of the refractive index model */
  static constexpr double A_ice=1.78;
  static constexpr double TransitionBoundary=0;
  // const double A_ice=1.775;
  // const double TransitionBoundary=14.9;
 
  /* Get the value of the B parameter for the refractive index model */
  double GetB(double z);
  
  /* Get the value of the C parameter for the refractive index model */
  double GetC(double z);

  /* Get the value of refractive index model for a given depth  */
  double Getnz(double z);

  /* E-feild Power Fresnel coefficient for S-polarised wave which is perpendicular to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R */
  double Refl_S(double thetai);
  
/* E-feild Power Fresnel coefficient for P-polarised wave which is parallel to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R */
  double Refl_P(double thetai);
    
  /* The temperature and attenuation model has been taken from AraSim which also took it from here http://icecube.wisc.edu/~araproject/radio/ . This is basically Matt Newcomb's icecube directory which has alot of information, plots and codes about South Pole Ice activities. Please read it if you find it interesting. */

  /* Temperature model:The model takes in value of depth z in m and returns the value of temperature in Celsius.*/
  double GetIceTemperature(double z);

  /* Ice Attenuation Length model: Takes in value of frequency in Ghz and depth z and returns you the value of attenuation length in m */
  double GetIceAttenuationLength(double z, double frequency);

  /* Setup the integrand to calculate the attenuation */
  double AttenuationIntegrand (double x, void * params);
  
/* Integrate over the integrand to calculate the attenuation */
  double IntegrateOverLAttn (double A0, double Frequency, double z0, double z1, double Lvalue);
  
/* Calculate the total attenuation for each type of ray */
  double GetTotalAttenuationDirect (double A0, double frequency, double z0, double z1, double Lvalue);
  
  double GetTotalAttenuationReflected (double A0, double frequency, double z0, double z1, double Lvalue);
  
  double GetTotalAttenuationRefracted (double A0, double frequency, double z0, double z1, double zmax, double Lvalue);
  
  /* Use GSL minimiser which relies on calculating function deriavtives. This function uses GSL's Newton's algorithm to find root for a given function. */
  double FindFunctionRootFDF(gsl_function_fdf FDF,double x_lo, double x_hi);

  /* Use GSL minimiser which uses GSL's false position algorithm to find root for a given function. */
  double FindFunctionRoot(gsl_function F,double x_lo, double x_hi);

  /* Use GSL minimiser which uses GSL's false position algorithm to find root for a given function. */
  double FindFunctionRootZmax(gsl_function F,double x_lo, double x_hi);

  /* Define the function that will be minimized to get the value of the depth of the turning point for a given refracted ray. This function basically requires the value of the L parameter to find the depth of the turning point.  This comes from the part of the fDnfR function where sqrt( n(z) - L ). This imposes the constraint then n(z)=L at the turning point of the ray from which we can find zmax. */
  struct Minnz_params { double a,l; };
  double GetMinnz(double x,void *params);

  /* Get the value of the depth of the turning point for the refracted ray */
  double GetZmax(double A, double L);

  /* Analytical solution describing ray paths in ice as function of depth */
  struct fDnfR_params { double a, b, c, l; };
  double fDnfR(double x,void *params);

  /* Analytical solution describing the ray path in ice as a function of the L parameter */
  struct fDnfR_L_params { double a, b, c, z; };
  double fDnfR_L(double x,void *params);

  /* The function used to calculate ray propogation time in ice */
  struct ftimeD_params { double a, b, c, speedc,l; };
  double ftimeD(double x,void *params);

  /* The function used to calculate ray geometric path in ice */
  double fpathD(double x,void *params);
  
  /* The set of functions starting with the name "fDa" are used in the minimisation procedure to find the launch angle (or the L parameter) for the direct ray */
  struct fDanfRa_params { double a, z0, x1, z1; };
  double fDa(double x,void *params);

  double fDa_df(double x,void *params);

  void fDa_fdf (double x, void *params,double *y, double *dy);

  /* The set of functions starting with the name "fRa" are used in the minimisation procedure to find the launch angle (or the L parameter) for the reflected ray */
  double fRa(double x,void *params);

  double fRa_df(double x,void *params);

  void fRa_fdf (double x, void *params,double *y, double *dy);

  /* This function is minimised to find the launch angle (or the L parameter) for the refracted ray */
  double fRaa(double x,void *params);

  double fRaa_df(double x,void *params);

  void fRaa_fdf (double x, void *params,double *y, double *dy);

  /* This functions works for the Direct ray and gives you back the launch angle, receive angle and propagation time of the ray together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
  double* GetDirectRayPar(double z0, double x1, double z1);
  
  /* This functions works for the Reflected ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
  double *GetReflectedRayPar(double z0, double x1 ,double z1);

  /* This functions works for the Refracted ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. It requires the launch angle of the reflected ray as an input. */
  double *GetRefractedRayPar(double z0, double x1 ,double z1, double LangR, double RangR, double checkzeroD, double checkzeroR);

  
  /* This function returns the x and z values for the full Direct ray path and prints out the ray path in a text file */
  void GetFullDirectRayPath(double z0, double x1, double z1,double lvalueD, vector <double> &x, vector <double> &z);

  /* This function returns the x and z values for the full Reflected ray path and prints out the ray path in a text file */
  void GetFullReflectedRayPath(double z0, double x1, double z1,double lvalueR, vector <double> &x, vector <double> &z);
  
  /* This function returns the x and z values for the full Refracted ray path and prints out the ray path in a text file */
  void GetFullRefractedRayPath(double z0, double x1, double z1, double zmax, double lvalueRa, vector <double> &x, vector <double> &z,int raynumber);

  /* function for plotting and storing all the rays */
  void PlotAndStoreRays(double x0,double z0, double z1, double x1, double zmax[2], double lvalues[4], double checkzeroes[4]);

  double *IceRayTracing(double x0, double z0, double x1, double z1);

  /* Analytical solution describing ray paths in ice as function of depth for constant refractive index*/
  double fDnfR_Cnz(double x,void *params);

  /* Analytical solution describing the ray path in ice as a function of the L parameter for constant refractive index*/
  double fDnfR_L_Cnz(double x,void *params);

  /* This function is minimised to find the launch angle (or the L parameter) for the reflected ray for constant refractive index*/
  double fRa_Cnz(double x,void *params);

  double fRa_Cnz_df(double x,void *params);

  void fRa_Cnz_fdf (double x, void *params,double *y, double *dy);

  /* This functions works for the Direct ray and gives you back the launch angle, receive angle and propagation time of the ray together with values of the L parameter. This for constant refractive index*/
  double* GetDirectRayPar_Cnz(double z0, double x1, double z1, double A_ice_Cnz);

  /* This functions works for the Reflected ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter. This is for constant refractive index*/
  double *GetReflectedRayPar_Cnz(double z0, double x1 , double z1, double A_ice_Cnz);

  /* This function returns the x and z values for the full Direct ray path and prints out the ray path in a text file. This is for a constant refractive index. */
  void GetFullDirectRayPath_Cnz(double z0, double x1, double z1, double lvalueD, double A_ice_Cnz,vector <double> &x, vector <double> &z);

  /* This function returns the x and z values for the full Reflected ray path and prints out the ray path in a text file. This is for a constant refractive index. */
  void GetFullReflectedRayPath_Cnz(double z0, double x1, double z1, double lvalueR, double A_ice_Cnz,vector <double> &x, vector <double> &z);

  /* function for plotting and storing all the rays. This is for constant refractive index. */
  void PlotAndStoreRays_Cnz(double x0,double z0, double z1, double x1, double lvalues[2], double A_ice_Cnz);

  /* This is the main raytracing function. x0 always has to be zero. z0 is the Tx depth in m and z1 is the depth of the Rx in m. Both depths are negative. x1 is the distance between them. This functions works for a constant refractive index */
  double *IceRayTracing_Cnz(double x0, double z0, double x1, double z1, double A_ice_Cnz); 

 /* The set of functions starting with the name "fDa" are used in the minimisation procedure to find the launch angle (or the L parameter) for the direct ray */
  double fDa_Air(double x,void *params);

  /* This functions works for the Direct ray and gives you back the launch angle, receive angle and propagation time of the ray together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
  double* GetDirectRayPar_Air(double z0, double x1, double z1);

  double *GeantRayTracer(double xT, double yT, double zT, double xR, double yR, double zR);
  
  /* Function that makes interpolation tables for raytracing */
  void MakeTable(double ShowerHitDistance,double zT);

  /* Function that calculates the interpolated value for raytracing. The rt parameter: 0 is for launch angle, 1 is for recieve angle, 2 is for propagation time, 3 is for distance */
  
  double GetInterpolatedValue(double xR, double zR, int rtParameter);
			      
  void GetRayTracingSolutions(double RxDepth, double Distance, double TxDepth, double TimeRay[2], double PathRay[2], double LaunchAngle[2], double RecieveAngle[2], int IgnoreCh[2], double IncidenceAngleInIce[2],vector <double> xRay[2], vector <double> zRay[2]);
  
}
#endif
