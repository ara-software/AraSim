#ifndef ICEMODEL_H
#define ICEMODEL_H

#include <cmath>
#include "EarthModel.h"
#include "Constants.h"
#include "Vector.h"
#include "Position.h"
//#include "Primaries.h"

class Interaction;
class Settings;

//Constants relating to all ice models
const double FIRNDEPTH=-150.;                // depth of the firn, in meters: currently a constant over all ice

class IceModel : public EarthModel {

protected:
  int ice_model;
  int DEPTH_DEPENDENT_N;
  int mooreBayFlag;

  //Information on horizons - what ice the balloon can see at each position along its path.
 
 
//   double volume_inhorizon_average; // average volume of ice seen by balloon
//   vector<int> ilon_inhorizon[NBNPOSITIONS_MAX]; // indices in lon and lat for bins in horizon for NPHI balloon positions along 80 deg latitude line.
//   vector<int> ilat_inhorizon[NBNPOSITIONS_MAX];
//   vector<int> easting_inhorizon[NBNPOSITIONS_MAX]; //indicies in easting and northing for bins in horizon for NPHI balloon positions along 80 deg latitude line.
//   vector<int> northing_inhorizon[NBNPOSITIONS_MAX];
//   double maxvol_inhorizon[NBNPOSITIONS_MAX]; // maximum volume of ice for a bin 

  //BEDMAP utility methods
  double Area(double latitude) const;

void ENtoLonLat(int e_coord, 
		  int n_coord,
		  double xLowerLeft,
		  double yLowerLeft,
		  
		  double& lon, 
		  double& lat) const;



  void WaterENtoLonLat(int e,
		       int n,
		       
		       double& lon,
		       double& lat) const;
  void LonLattoEN(double lon, 
		  double lat,
		  double xLowerLeft,
		  double yLowerLeft,
		
		  int& e_coord, 
		  int& n_coord) const;
 
  void GroundLonLattoEN(double lon, 
			double lat,
			
			int& e_coord, 
			int& n_coord) const;
  void WaterLonLattoEN(double lon,
		       double lat,
		       
		       int& e_coord,
		       int& n_coord) const;

  //BEDMAP data input methods
  void ReadIceThickness();
  void ReadGroundBed();
  void ReadWaterDepth();

private:

  const static int N_sheetup=2810;
  double d_sheetup[N_sheetup], l_sheetup[N_sheetup];
  const static int N_shelfup=420;
  double d_shelfup[N_shelfup], l_shelfup[N_shelfup];
  const static int N_westlandup=420;
  double d_westlandup[N_westlandup],l_westlandup[N_westlandup];

  const static int N_sheetdown=2810;
  double d_sheetdown[N_sheetup], l_sheetdown[N_sheetdown];
  const static int N_shelfdown=420;
  double d_shelfdown[N_shelfdown], l_shelfdown[N_shelfdown];
  const static int N_westlanddown=420;
  double d_westlanddown[N_westlanddown],l_westlanddown[N_westlanddown];



public:


  //BEDMAP data
  double ice_thickness_array[1200][1000];  //thickness of the ice
  double ground_elevation[1068][869]; //elevation above geoid at which ice starts
  double water_depth[1200][1000]; //depth of water under ice

  void IceENtoLonLat(int e,
		     int n,
		     
		     double& lon,
		     double& lat) const;  
  void GroundENtoLonLat(int e,
			int n,
			
			double& lon,
			double& lat) const;

  //  const static int NBNPOSITIONS_MAX=26000;
  //double volume_inhorizon[NBNPOSITIONS_MAX]; // volume of ice within horizon for each balloon phi position 
//  IceModel();   //default constructor
  IceModel(int model=0,int earth_model=0,int mooreBay=0);
  ~IceModel();

  void setUpIceModel(int model=0);
  double IceThickness(double lon,double lat) const;
  double IceThickness(const Position& pos) const;
  double Surface(double lon,double lat) const;
  double Surface(const Position& pos) const;
  double SurfaceAboveGeoid(double lon,double lat) const;
  double SurfaceAboveGeoid(const Position& pos) const;
  double WaterDepth(double lon,double lat) const;
  double WaterDepth(const Position& pos) const;

  Position PickBalloonPosition() const;
  //void GetMAXHORIZON(double bn_altitude); // get upper limit on the horizon wrt the balloon.
  int RossIceShelf(const Position &position) const; 
  int IceOnWater(const Position &postition) const;
  int RossExcept(const Position &position) const;
  int RonneIceShelf(const Position &position) const;
  int WestLand(const Position &pos) const; 
  int AcceptableRfexit(const Vector &nsurf_rfexit,const Position &rfexit,const Vector &n_exit2rx) const; 
  //double GetBalloonPositionWeight(int ibnpos) const;
  int OutsideAntarctica(const Position &pos) const;
  int OutsideAntarctica(double lat) const;
  Position WhereDoesItEnterIce(const Position &posnu,
			       const Vector &nnu,
			       double stepsize) const;


  int WhereDoesItEnter_sphere(const Position &sphere_in, const Vector &nnu, Position &r_in ) const;


  Position WhereDoesItEnter(const Position &posnu,const Vector &nnu) const;
  Position WhereDoesItLeave(const Position &posnu,const Vector &nnu) const;

  //void CreateHorizons(int whichpath,int reduceballoonpositions,Balloon *bn1,double theta_bn,double phi_bn,double altitude_bn,ofstream &foutput);
  Vector GetSurfaceNormal(const Position &r_out) const; //overloaded from EarthModel to include procedures for new ice models.
  double GetN(double depth) const;
  double GetN(const Position &pos) const;
  double EffectiveAttenuationLength(const Position &pos, const int &whichray) const;
  double EffectiveAttenuationLength(Settings *settings1, const Position &pos, const int &whichray) const;
  
  void IceLonLattoEN(double lon, 
		     double lat,
		     
		     int& e_coord, 
		     int& n_coord) const;
  int Getice_model();

  void GetFresnel_slappy (double i_ang, double n1, double n2, double &r, double &t);
  void GetFresnel_pokey (double i_ang, double n1, double n2, double &r, double &t);
void GetFresnel (
        double launch_angle, double rec_angle,
        double refl_angle, Position &posnu, Vector &launch_vector, Vector &rec_vector, Settings *settings1, double &fresnel, double &mag,
        Vector &Pol // will read the polarization at the source and return polarization at the target antenna
        );

  // new ARA ice attenuation measurement values (at 300 MHz)
  //
  double ARA_IceAtten_Depth[100]; // maximum 100 bins
  double ARA_IceAtten_Length[100];
  int ARA_IceAtten_bin;

  double GetARAIceAttenuLength(double depth);

  //  void FillArraysforTree(double lon_ground[1068][869],double lat_ground[1068][869],double lon_ice[1200][1000],double lat_ice[1200][1000],double lon_water[1200][1000],double lat_water[1200][1000]);

  // below three members are copied from icemc icemodel.
//--------------------------------------------------
//   int PickUnbiased(int inu, Interaction *interaction1, IceModel *antarctica);
//   int PickNear(int inu, Interaction *interaction1, IceModel *antarctica);
//   int WhereDoesItEnterIce(const Position &posnu,
// 			       const Vector &nnu,
// 			       double stepsize,
// 			       Position &r_enterice);
// 
//   int WhereDoesItExitIce(int inu,const Position &posnu,
// 			 const Vector &nnu,
// 			 double stepsize,
// 			 Position &r_enterice);
//-------------------------------------------------- 
  // end three copied members from icemc icemodel.


  ClassDef(IceModel,1);



}; //class IceModel
// input files for Crust 2.0
const string crust20_in="data/outcr"; // Crust 2.0 data
const string crust20_out="altitudes.txt"; // output file for plotting


#endif //ICEMODEL_H
