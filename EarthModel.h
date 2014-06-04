#ifndef EARTHMODEL_H
#define EARTHMODEL_H

#include <string>
#include "TRandom3.h"
#include "TObject.h"
class Position;
class Vector;

class Primaries;
class Settings;
class IceModel;
class Secondaries;
//class TRandom3;
using std::string;




////////////////////////////////////////////////////////////////////////////////////////////////
//class EarthModel
//This class contains a model of the physical structure of the Earth, with information on shape,
//density at all points inside, and depth of water and ice on the surface.
//
//   methods:
//
//IceThickness : Returns the thickness of ice in meters at a given Position or lon/lat.
//WaterDepth  : Returns the depth of water in meters at a given Position or lon/lat.
//Surface   : Returns the distance in meters from the center of the Earth to the surface,
//            i.e. to the start of air.
//SurfaceAboveGeoid : Returns the distance in meters from the local geoid to the start of air.
//Geoid  : Returns the height in meters of the geoid at a given Position or lon/lat.
//RockSurface  : Returns the distance in meters from the center of the Earth to the end of rock
//               (the begninning of the ice/water/air)
//GetSurfaceNormal  : Returns a unit vector pointing in the direction of the surface normal at 
//                    a given Position.
//Getchord
//
////////////////////////////////////////////////////////////////////////////////////////////////

class EarthModel {

private:
  TRandom3 Rand3;


protected:
  int EARTH_MODEL;
  int CONSTANTICETHICKNESS;
  int CONSTANTCRUST;
  int FIXEDELEVATION;
  int FLATSURFACE;
  int weightabsorption;
  

  // pick method to step neutrino through the earth and get attenuation length
  static const int getchord_method=2; // 1=core,mantle,crust; 2=crust 2.0 model  

  static const double GEOID_MAX; // parameters of geoid model
  static const double GEOID_MIN; 
  static const double COASTLINE; // if the rf leaves from beyond this "coastline" (in degrees of latitude relative to south pole) it's not in antarctica.  Real coastline is always closer than this. 

  // parameters of the Crust 2.0 earth model
  static const int NLON=180; // number of bins in longitude for crust 2.0
  static const int NLAT=90;  // number of bins in latitude
  static const int NPHI=180; // bins in longitude for visible ice in horizon
  static const double MAXTHETA; // maximum value of theta in degrees in polar coordinates
  double thetastep; // how big do you step in theta-> always 2deg with Crust 2.0
  double phistep; // how big do you step in phi->always 2deg in Crust 2.0
  static const int ILAT_MAX;
  double surfacer[NLON][NLAT]; // elevation at the surface (top of ice) compared to geoid (in meters)
  double icer[NLON][NLAT];  // elevation at the *bottom* of ice layer (in meters)
  double waterr[NLON][NLAT]; // elevation at the bottom of water layer (in meters)
  double softsedr[NLON][NLAT]; // elevation at the bottom of soft set layer (in meters)
  double hardsedr[NLON][NLAT]; // elev at bottom of hard sed layer (in meters)
  double uppercrustr[NLON][NLAT]; // elev at bottom of upper crust layer (in meters)
  double middlecrustr[NLON][NLAT]; // elev at bottom of middle crust layer (in meters)
  double lowercrustr[NLON][NLAT]; // elev at bottom of lower crust layer (in meters)
  double geoid[NLAT];     // realistic shape of earth-radius at each latitude (in meters)
  double MIN_ALTITUDE_CRUST; // maximum depth of crust- determines when to start stepping
  //double MAX_VOL; // maximum volume of ice in a bin in Crust 2.0 - not used
  double elevationarray[NLON][NLAT]; // If no water, measures the elevation (relative to geoid, in meters) of the top of the ice or rock (i.e., air interface).  If there is water, measures elevation to bottom of water. (There may or may not be ice on top of the water.) 
  double waterthkarray[NLON][NLAT];  // thickness of water layer (in km)
  double icethkarray[NLON][NLAT]; // thickness of ice layer (in km)
  double softsedthkarray[NLON][NLAT]; // thickness of soft sed layer (in km)
  double hardsedthkarray[NLON][NLAT]; // thickness of hard sed layer (in km)
  double uppercrustthkarray[NLON][NLAT]; // thickness of upper crust layer (in km)
  double middlecrustthkarray[NLON][NLAT]; // thickness of middle crust layer (in km)
  double lowercrustthkarray[NLON][NLAT]; // thickness of lower crust layer (in km)
  double crustthkarray[NLON][NLAT]; // total thickness of crust (in km)
  double waterdensityarray[NLON][NLAT]; // density of water layer bin by bin
  double icedensityarray[NLON][NLAT]; // density of ice layer bin by bin
  double softseddensityarray[NLON][NLAT]; // density of soft sed layer
  double hardseddensityarray[NLON][NLAT]; // density of hard sed layer
  double uppercrustdensityarray[NLON][NLAT]; // density of upper crust layer
  double middlecrustdensityarray[NLON][NLAT]; // density of middle crust layer
  double lowercrustdensityarray[NLON][NLAT]; // density of lower crust layer
  double area[NLAT]; // area of a bin at a given latitude- calculated once
  double average_iceth; // average ice thickness over the continent-calculated once

  /////////////////////////////////////
  //methods
  void ReadCrust(string);
  double SmearPhi(int ilon);
  double SmearTheta(int ilat) ;
  double dGetTheta(int itheta) const;
  double dGetPhi(int ilon) const;
  void GetILonILat(const Position&,int& ilon,int& ilat) const;
  double GetLat(double theta) const;
  double GetLon(double phi) const;



public:
  EarthModel(int model = 0,int WEIGHTABSORPTION_SETTING=1);
  virtual ~EarthModel();
// properties of the Earth
//  static const double R_EARTH=6.378140E6;        // radius of Earth in m at bulge
  static const double R_EARTH;        // radius of Earth in m at bulge

  double radii[3];
  // = {1.2e13,(EarthModel::R_EARTH-4.0E4)*(EarthModel::R_EARTH-4.0E4),EarthModel::R_EARTH*EarthModel::R_EARTH}; // average radii of boundaries between earth layers

    double volume; // sums the volume of medium (ice or salt)
  virtual double Geoid(double latitude) const;
  virtual double Geoid(const Position &pos) const;
  virtual double IceThickness(double lon,double lat) const;
  virtual double IceThickness(const Position& pos) const;
  virtual double Surface(double lon,double lat) const;
  virtual double Surface(const Position& pos) const;
  virtual int InFirn(const Position& pos) const;
  virtual double SurfaceDeepIce(const Position& pos) const;
  virtual double SurfaceAboveGeoid(double lon,double lat) const;
  virtual double SurfaceAboveGeoid(const Position& pos) const;
  virtual double WaterDepth(double lon,double lat) const;
  virtual double WaterDepth(const Position& pos) const;
  virtual double RockSurface(double lon,double lat) const;
  virtual double RockSurface(const Position& pos) const;
  int Getchord(double len_int_kgm2,
		       const Position &r_in,
		       const Position &posnu,
		       int inu,

		       double& chord, 
		       double& weight1,
		       double& nearthlayers,
		       double myair,
		       double& total_kgm2,
		       int& crust_entered,
		       int& mantle_entered,
		       int& core_entered) ; 


  int Getchord(Primaries *primary1, Settings *settings1,IceModel *antarctica1, Secondaries *sec1,
			 double len_int_kgm2,
			 const Position &earth_in, // place where neutrino entered the earth
			 const Position &r_enterice,
			 const Position &nuexitice,
			 
			 const Position &posnu, // position of the interaction
			 int inu,
			 double& chord, // chord length
			 double& probability_tmp, // weight
			 double& weight1_tmp,
			 double& nearthlayers, // core, mantle, crust
			 double myair,
			 double& total_kgm2, // length in kg m^2
			 int& crust_entered, // 1 or 0
			 int& mantle_entered, // 1 or 0
			 int& core_entered, 
                         string thisnuflavor,double pnu, double Etau_final,
			 int nu_nubar, int currentint,int taumodes1, double *myxarray, double *myEarray, double *myyweightarray,
			 double *mytausurvarray, double& tauweight, double& tauchord, double *avgdensityarray, double *densityarray);


  int Getchord(Primaries *primary1, Settings *settings1,IceModel *antarctica1, Secondaries *sec1,
			 double len_int_kgm2,
			 const Position &earth_in, // place where neutrino entered the earth
			 const Position &r_enterice,
			 const Position &nuexitice,
			 
			 const Position &posnu, // position of the interaction
			 int inu,
			 double& chord, // chord length
			 double& probability_tmp, // weight
			 double& weight1_tmp,
			 double& nearthlayers, // core, mantle, crust
			 double myair,
			 double& total_kgm2, // length in kg m^2
			 int& crust_entered, // 1 or 0
			 int& mantle_entered, // 1 or 0
			 int& core_entered, 
                         string thisnuflavor,double pnu, double Etau_final,
			 int nu_nubar, int currentint,int taumodes1
			 );


  int Getchord(Primaries *primary1, Settings *settings1,IceModel *antarctica1, Secondaries *sec1,
			 double len_int_kgm2,
			 const Position &earth_in, // place where neutrino entered the earth
			 const Position &r_enterice,
			 const Position &nuexitice,
			 
			 const Position &posnu, // position of the interaction
			 int inu,
			 double& chord, // chord length
			 double& probability_tmp, // weight
			 double& weight1_tmp,
			 double& nearthlayers, // core, mantle, crust
			 double myair,
			 double& total_kgm2, // length in kg m^2
			 int& crust_entered, // 1 or 0
			 int& mantle_entered, // 1 or 0
			 int& core_entered, 
                         string thisnuflavor,double pnu, double Etau_final,
                         int nu_nubar, int currentint,int taumodes1, double lpath
			 );




  double GetDensity(double altitude, const Position earth_in, const Position posnu,
			      int& crust_entered, // 1 or 0
			      int& mantle_entered, // 1 or 0
			      int& core_entered);


  double GetDensity1(double altitude, const Position earth_in, const Position posnu,
			      int& crust_entered, // 1 or 0
			      int& mantle_entered, // 1 or 0
			       int& core_entered, double& abovesurface);



  Vector GetSurfaceNormal(const Position &r_out) const;
  static double LongtoPhi_0isPrimeMeridian(double longitude); // convert longitude to phiwith 0 longitude being the prime meridian
  static double LongtoPhi_0is180thMeridian(double longitude); // convert longitude to phi with 0 longitude being at the 180th meridian
  void EarthCurvature(double *array,double depth_temp); // adjusts coordinates within the mine to account for the curvature of the earth.

  static double GetGeoid(double latitude);
  static double GetCOASTLINE();

  ClassDef(EarthModel,1);


}; //class EarthModel


const double densities[3]={14000.,3400.,2900.}; // average density of each earth layer

#endif
