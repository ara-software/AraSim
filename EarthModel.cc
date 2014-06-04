#include "Constants.h"
//#include "Signal.h"
#include "signal.hh"
//#include "earthmodel.hh"
#include "EarthModel.h"
//#include "icemodel.hh"
#include "IceModel.h"
#include <cmath>
#include "Tools.h"
//#include "vector.hh"
#include "Vector.h"
//#include "position.hh"
#include "Position.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "signal.hh"

#include "Primaries.h"
#include "Settings.h"
#include "secondaries.hh"

ClassImp(EarthModel);


using std::cout;
using std::endl;
using std::ios;
using std::fstream;


const double EarthModel::COASTLINE(30.);
const double EarthModel::MAXTHETA(180.);
const int EarthModel::ILAT_MAX((int)((COASTLINE/MAXTHETA)*(double)NLAT+0.00001)); // corresponding latitude bin to "coastline"
const double EarthModel::GEOID_MAX(6.378137E6); // parameters of geoid model
const double EarthModel::GEOID_MIN(6.356752E6); // from Geodetic Reference System 1980, Bulletin Geodesique, Vol 54:395,1980. // The previous reference gave issue number 3 instead of page number 395

// test ClassDef works with static const double
const double EarthModel::R_EARTH(6.378140E6);

double EarthModel::GetCOASTLINE() {
    return COASTLINE;
}

EarthModel::EarthModel(int model,int WEIGHTABSORPTION_SETTING) {

  radii[0]=1.2e13;
  radii[1]=(EarthModel::R_EARTH-4.0E4)*(EarthModel::R_EARTH-4.0E4);
  radii[2]=(EarthModel::R_EARTH*EarthModel::R_EARTH); // average radii of boundaries between earth layers


  //  cout << "In EarthModel, model is " << model << "\n";
  weightabsorption= WEIGHTABSORPTION_SETTING;

  CONSTANTICETHICKNESS = (int) (model / 1000);
  model -= CONSTANTICETHICKNESS * 1000;

  CONSTANTCRUST = (int) (model / 100);
  model -= CONSTANTCRUST * 100;

  FIXEDELEVATION = (int) (model / 10);
  model -= FIXEDELEVATION * 10;

  EARTH_MODEL = model;
  //cout<<"CONSTICETHK = "<<CONSTANTICETHICKNESS<<", CNSTCRST = "<<CONSTANTCRUST<<", FIXDELV = "<<FIXEDELEVATION<<", EARTH_MODEL = "<<EARTH_MODEL<<endl;
  for (int i=0;i<NLON;i++) {
    
    Tools::Zero(elevationarray[i],NLAT);
    Tools::Zero(waterthkarray[i],NLAT);
    Tools::Zero(icethkarray[i],NLAT);
    Tools::Zero(softsedthkarray[i],NLAT);
    Tools::Zero(hardsedthkarray[i],NLAT);
    Tools::Zero(uppercrustthkarray[i],NLAT);
    Tools::Zero(middlecrustthkarray[i],NLAT);
    Tools::Zero(lowercrustthkarray[i],NLAT);
    Tools::Zero(crustthkarray[i],NLAT);
    
    
    Tools::Zero(surfacer[i],NLAT);
    Tools::Zero(icer[i],NLAT);
    Tools::Zero(waterr[i],NLAT);
    Tools::Zero(softsedr[i],NLAT);
    Tools::Zero(hardsedr[i],NLAT);
    Tools::Zero(uppercrustr[i],NLAT);
    Tools::Zero(middlecrustr[i],NLAT);
    Tools::Zero(lowercrustr[i],NLAT);
    
    Tools::Zero(waterdensityarray[i],NLAT);
    Tools::Zero(icedensityarray[i],NLAT);
    Tools::Zero(softseddensityarray[i],NLAT);
    Tools::Zero(hardseddensityarray[i],NLAT);
    Tools::Zero(uppercrustdensityarray[i],NLAT);
    Tools::Zero(middlecrustdensityarray[i],NLAT);
    Tools::Zero(lowercrustdensityarray[i],NLAT);
        
  } //Zero Earth model arrays

  // see monte carlo note #17
  for (int i=0;i<NLAT;i++) {
    geoid[i]=GEOID_MIN*GEOID_MAX/sqrt(pow(GEOID_MIN,2.)-(pow(GEOID_MIN,2.)-pow(GEOID_MAX,2.))*pow(cos(dGetTheta(i)),2.));
  } //for


  // Crust 2.0 is binned in 2deg x 2deg bins, area of bin depends on latitude.
  // calculating surface area of bins
  phistep=2*PI/(double)NLON;
  thetastep=(MAXTHETA*RADDEG)/NLAT;
  for (int i=0;i<NLAT;i++) {
    area[i]=phistep*(cos(dGetTheta(i))-cos(dGetTheta(i+1)))*pow(geoid[i],2.);
  } //for


  if (EARTH_MODEL == 0)
    ReadCrust(crust20_in);
  else {
    cout<<"Error!  Unknown Earth model requested!  Defaulting to Crust 2.0 model.\n";
    ReadCrust(crust20_in);
  } //else

} //EarthModel constructor (int mode)

EarthModel::~EarthModel() {} //EarthModel destructor - no dynamic variables, nothing to delete


 double EarthModel::LongtoPhi_0isPrimeMeridian(double longitude) {

  double phi;
  // convert longitude (-180 to 180) to phi (0 to 2pi) wrt +x
  // in radians
  phi=(90-longitude); 
  if (phi<0.)
    phi+=360.;

  phi=phi*RADDEG;

  return phi;
}
 double EarthModel::LongtoPhi_0is180thMeridian(double longitude) {

  double phi;
  // convert longitude (0 to 360) to phi (0 to 2pi) wrt +x
  
  phi=(270.-longitude); 
  if (phi<0.)
    phi+=360.;

  phi=phi*RADDEG;
  // in radians

  return phi;
}

 double EarthModel::GetGeoid(double latitude) {

  return (GEOID_MIN*GEOID_MAX/
	  sqrt(GEOID_MIN*GEOID_MIN-(GEOID_MIN*GEOID_MIN-GEOID_MAX*GEOID_MAX)*
	       cos(latitude*RADDEG)*cos(latitude*RADDEG)));
 }


 double EarthModel::Geoid(double latitude) const {
  // latitude here is 0 at the south pole and 180 at the north pole
 
  return (GEOID_MIN*GEOID_MAX/
	  sqrt(GEOID_MIN*GEOID_MIN-(GEOID_MIN*GEOID_MIN-GEOID_MAX*GEOID_MAX)*
	       cos(latitude*RADDEG)*cos(latitude*RADDEG)));
} //Geoid(lat)

 double EarthModel::Geoid(const Position &pos) const {
  return Geoid(pos.Lat());
} //Geoid(Position)

 double EarthModel::IceThickness(double lon,double lat) const {
  return icethkarray[(int)(lon/2)][(int)(lat/2)]*1000.;
} //IceThickness(lon,lat)

 double EarthModel::IceThickness(const Position& pos) const {
  return IceThickness(pos.Lon(),pos.Lat());
} //IceThickness(Position)
 int EarthModel::InFirn(const Position& pos) const {
  if (pos.Mag()-Surface(pos)<FIRNDEPTH)
    return 0;
  return 1;
} //InFirn(Position)
 double EarthModel::SurfaceDeepIce(const Position& pos) const { // surface of the deep ice (where you reach the firn)
  return  surfacer[(int)(pos.Lon()/2)][(int)(pos.Lat()/2)] + geoid[(int)(pos.Lat()/2)] + FIRNDEPTH;
} //Surface(lon,lat)

 double EarthModel::Surface(double lon,double lat) const {
  return surfacer[(int)(lon/2)][(int)(lat/2)] + geoid[(int)(lat/2)];
} //Surface(lon,lat)

 double EarthModel::Surface(const Position& pos) const {
  return surfacer[(int)(pos.Lon()/2)][(int)(pos.Lat()/2)] + geoid[(int)(pos.Lat()/2)];
} //Surface(Position)

 double EarthModel::RockSurface(double lon,double lat) const {
  return (Surface(lon,lat) - IceThickness(lon,lat) - WaterDepth(lon,lat));
} //RockSurface(lon,lat)

 double EarthModel::RockSurface(const Position& pos) const {
  return RockSurface(pos.Lon(),pos.Lat());
} //RockSurface(lon,lat)

 double EarthModel::SurfaceAboveGeoid(double lon,double lat) const {
  return surfacer[(int)(lon/2)][(int)(lat/2)];
} //SurfaceAboveGeoid(lon,lat)

 double EarthModel::SurfaceAboveGeoid(const Position& pos) const {
  return surfacer[(int)(pos.Lon()/2)][(int)(pos.Lat()/2)];
} //SurfaceAboveGeoid(Position)

 double EarthModel::WaterDepth(double lon,double lat) const {
  return waterthkarray[(int)(lon/2)][(int)(lat/2)]*1000;
} //WaterDepth(lon,lat)

 double EarthModel::WaterDepth(const Position& pos) const {
  return WaterDepth(pos.Lon(),pos.Lat());
} //WaterDepth(Position)

 double EarthModel::GetLat(double theta) const {
  return theta*DEGRAD;
} //GetLat

 double EarthModel::GetLon(double phi) const {
  // input is phi in radians wrt +x
  double phi_deg = phi*DEGRAD; 
  if (phi_deg > 270)   
    phi_deg = phi_deg - 360.; // this puts it from -90 to 270

  return (360.*3./4. - phi_deg); // returns 0 to 360 degrees (going from -180 to 180 deg longitude like Crust 2.0 does)
} //GetLon



double EarthModel::GetDensity(double altitude, const Position earth_in, const Position posnu,
			      int& crust_entered, // 1 or 0
			      int& mantle_entered, // 1 or 0
			      int& core_entered){
  
        Position where = earth_in;
	//cout<<"where is "<<where<<"\n";
  	double x = 0; //where.Mag();
	double lon = where.Lon();
	double lat = where.Lat();
	//cout <<"Lon and Lat are "<<lon<<","<<lat<<"\n";

	int ilon = (int)(lon/2);
	int ilat = (int)(lat/2);

	double ddensity =0; //initilize ddensity

	double surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point

	double local_icethickness = this->IceThickness(lon,lat);
	double local_waterdepth = WaterDepth(lon,lat);

	//altitude=altitude-Geoid(lat); // what is the altitude of the entrance point

	if(altitude>surface_elevation+0.1){ // if it is above the surface, it's messed up
	  //cout << "neutrino entrance point is above the surface. Density0 \n";
	  
	  // cout <<"altitude is "<<altitude<<"\n";
	    }
	if(altitude>surface_elevation+0.1){
	  ddensity=1.25;
	  //cout <<"density is air! \n";
	}
	if (altitude<=surface_elevation+0.1 && altitude>(surface_elevation-local_icethickness)) // the 0.1 is just to take care of precision issues.   It could have been 0.01 or 0.001.
		ddensity=icedensityarray[ilon][ilat]*1000;
	  
	else if (altitude<=(surface_elevation-local_icethickness) && altitude>(surface_elevation-local_icethickness-local_waterdepth))
		ddensity=waterdensityarray[ilon][ilat]*1000;
	else if (altitude<=(surface_elevation-local_icethickness-local_waterdepth) && altitude>softsedr[ilon][ilat]) {
		ddensity=softseddensityarray[ilon][ilat]*1000;
	       	crust_entered=1; //Switch that lets us know we've penetrated into the crust
		 } //end if
	else if (altitude<=softsedr[ilon][ilat] && altitude>hardsedr[ilon][ilat])
		ddensity=hardseddensityarray[ilon][ilat]*1000;
	else if (altitude<=hardsedr[ilon][ilat] && altitude>uppercrustr[ilon][ilat])
		ddensity=uppercrustdensityarray[ilon][ilat]*1000;
	else if (altitude<=uppercrustr[ilon][ilat] && altitude>middlecrustr[ilon][ilat])
		ddensity=middlecrustdensityarray[ilon][ilat]*1000;
	else if (altitude<=middlecrustr[ilon][ilat] && altitude>lowercrustr[ilon][ilat])
		ddensity=lowercrustdensityarray[ilon][ilat]*1000;
	else if (altitude<=lowercrustr[ilon][ilat])
		ddensity=densities[1];
	
	    return ddensity;

}//Get Density



double EarthModel::GetDensity1(double altitude, const Position earth_in, const Position posnu,
			      int& crust_entered, // 1 or 0
			      int& mantle_entered, // 1 or 0
			       int& core_entered, double& abovesurface){
  
        Position where = earth_in;
	//cout<<"where is "<<where<<"\n";
  	double x = 0; //where.Mag();
	double lon = where.Lon();
	double lat = where.Lat();
	//cout <<"Lon and Lat are "<<lon<<","<<lat<<"\n";

	int ilon = (int)(lon/2);
	int ilat = (int)(lat/2);

	double ddensity =0; //initilize ddensity

	double surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point

	double local_icethickness = this->IceThickness(lon,lat);
	double local_waterdepth = WaterDepth(lon,lat);

	//altitude=altitude-Geoid(lat); // what is the altitude of the entrance point

	if(altitude>surface_elevation+0.1){ // if it is above the surface, it's messed up
	  //cout << "neutrino entrance point is above the surface density1.\n";
	 	  //cout <<"altitude is "<<altitude<<"\n";
	  abovesurface=1;
	    }
	if(altitude>surface_elevation+0.1){
	  ddensity=1.25;
	  // cout<<"density is air. \n";
	}
	if (altitude<=surface_elevation+0.1 && altitude>(surface_elevation-local_icethickness)) // the 0.1 is just to take care of precision issues.   It could have been 0.01 or 0.001.
		ddensity=icedensityarray[ilon][ilat]*1000;
	  
	else if (altitude<=(surface_elevation-local_icethickness) && altitude>(surface_elevation-local_icethickness-local_waterdepth))
		ddensity=waterdensityarray[ilon][ilat]*1000;
	else if (altitude<=(surface_elevation-local_icethickness-local_waterdepth) && altitude>softsedr[ilon][ilat]) {
		ddensity=softseddensityarray[ilon][ilat]*1000;
	       	crust_entered=1; //Switch that lets us know we've penetrated into the crust
		 } //end if
	else if (altitude<=softsedr[ilon][ilat] && altitude>hardsedr[ilon][ilat])
		ddensity=hardseddensityarray[ilon][ilat]*1000;
	else if (altitude<=hardsedr[ilon][ilat] && altitude>uppercrustr[ilon][ilat])
		ddensity=uppercrustdensityarray[ilon][ilat]*1000;
	else if (altitude<=uppercrustr[ilon][ilat] && altitude>middlecrustr[ilon][ilat])
		ddensity=middlecrustdensityarray[ilon][ilat]*1000;
	else if (altitude<=middlecrustr[ilon][ilat] && altitude>lowercrustr[ilon][ilat])
		ddensity=lowercrustdensityarray[ilon][ilat]*1000;
	else if (altitude<=lowercrustr[ilon][ilat])
		ddensity=densities[1];
	
	    return ddensity;

}//Get Density




 int EarthModel::Getchord(double len_int_kgm2,
				const Position &earth_in, // place where neutrino entered the earth
				const Position &posnu, // position of the interaction
				int inu,

				double& chord, // chord length
				double& weight1, // weight
				double& nearthlayers, // core, mantle, crust
				double myair,
				double& total_kgm2, // length in kg m^2
				int& crust_entered, // 1 or 0
				int& mantle_entered, // 1 or 0
				int& core_entered)  {

  Vector chord3;
  Vector nchord;
  double x=0;
  double lat,lon;
  int ilon,ilat;


  total_kgm2 = 0; //Initialize column density
  nearthlayers=0; // this counts the number of earth layers the neutrino traverses.
  crust_entered=0;
  mantle_entered=0;
  core_entered=0;
  // Want to find probability that the neutrino survives its trip
  // through the earth.
  
  //Find the chord, its length and its unit vector.
  chord3 = posnu - earth_in;
  chord=chord3.Mag();
  nchord = chord3 / chord;
  
  if (chord<=1) {
    cout << "short chord " << chord << "\n";
    return 0;
  }
  if (chord>2.*R_EARTH+1000) {
    cout << "bad chord" << " " << chord << ".  Event is " << inu << "\n";
  }

  Position where=earth_in;

  // the sin of the angle between the neutrino path and the 
  // radial vector to its earth entrance point determines
  // if it will get to the next layer down.
  double costh=(where*nchord)/where.Mag();
  double sinth=sqrt(1-costh*costh);
  double distance=0;
  double halfchord=0;
  
  if (getchord_method<1 || getchord_method>3)
    cout << "Bogus method!\n";
  
  
  // we are really focusing on method 2 - method 1 has not been maintenanced in a long time!!
  // use at your own risk.
  if (getchord_method==1) {
    double L=0;
    weight1=0;
        
    if (sinth>sqrt(radii[1]/radii[2])) {
      nearthlayers++;
      
      // these only skim the first layer.
      L=len_int_kgm2/densities[2];
      
      weight1=exp(-posnu.Distance(where)/L);
    }
    else {
      nearthlayers++;
      
      // these get to the second layer down.
      L=len_int_kgm2/densities[2];
      // compute distance before it gets to the next layer.
      halfchord=sqrt(radii[1]-radii[2]*sinth*sinth);
      distance=sqrt(radii[2])*costh-halfchord;
      
      weight1=exp(-distance/L);

      // get position where it enters the second layer.
      where = earth_in + distance*nchord;

      // determine if it enters the core or not.
      costh=(where*nchord)/where.Mag();
      sinth=sqrt(1-costh*costh);
      
      if (sinth>sqrt(radii[0]/radii[1])) {
	

	halfchord=sqrt(radii[1])*costh;
	nearthlayers++;
	
	// these do not enter the core.
	L=len_int_kgm2/densities[1];

	
	weight1 *= exp(-2*halfchord/L); 

	L=len_int_kgm2/densities[2];
	// this is where it exits the second layer and enters the crust again.
	where = where + 2*halfchord*nchord;
	weight1*=exp(-where.Distance(posnu)/L);
      }
      else {
	nearthlayers++;
	// these enter the core.
	L=len_int_kgm2/densities[1];
	
	// compute distance before entering the core.
	halfchord=sqrt(radii[0]-radii[1]*sinth*sinth);
	distance=sqrt(radii[1])*costh-halfchord;
	weight1*=exp(-distance/L);

	// go through the core.
	L=len_int_kgm2/densities[0];
	weight1*=exp(-2*halfchord/L);

	// go through 2nd layer again.
	L=len_int_kgm2/densities[1];
	weight1*=exp(-distance/L);
	
	// through the crust and end up at posnu.
	L=len_int_kgm2/densities[2];
	
	where = where + (2*distance+2*halfchord)*nchord;
	weight1*=exp(-where.Distance(posnu)/L);
      } //else
    } //else
  } //if getchord_method==1
  if (getchord_method==2) {

    x=0; // x is the distance you move as you step through the earth.
    
    lon = where.Lon();
    lat = where.Lat();
    ilon = (int)(lon/2);
    ilat = (int)(lat/2);
    
    double surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
    double local_icethickness = this->IceThickness(lon,lat);
    double local_waterdepth = WaterDepth(lon,lat);
    double altitude=0;
    weight1=1;
    double step=Tools::dMin(len_int_kgm2/densities[1]/10,500.); //how big is the step size
    //double step=Tools::dMin(len_int_kgm2/densities[1]/10,5.); //how big is the step size
//--------------------------------------------------
//     cout<<"len_int_kgm2 = "<<len_int_kgm2<<endl;
//     cout<<"densities[1] = "<<densities[1]<<endl;
//     cout<<"Getchord.step = "<<step<<endl;
//-------------------------------------------------- 
    // either 1/10 of an interaction length in the mantle or 500 m, whichever is smaller.
    // 500 m is approximately the feature size in Crust 2.0.
    //------------------added on Dec 8------------------------
    weight1*=exp(-myair/len_int_kgm2);//add atmosphere attenuation // fenfang's atten. due to atmosphere
    //------------------added on Dec 8------------------------
    total_kgm2+=myair;

    double L=0;

    double ddensity=Signal::RHOAIR;
    nearthlayers=1;
    
    if (where*nchord>0.)  { // look at direction of neutrino where it enters the earth.
      cout << "This one's trouble.  Neutrino exit point looks more like an entrance point.  Event is " << inu << "\n";
      cout << "where is " << where[0] << " " << where[1] << " " << where[2] << "\n";
      cout << "nchord is " << nchord[0] << " " << nchord[1] << " " << nchord[2] << "\n";
      cout << "dot product is " << where*nchord/sqrt(where*where) << "\n";
      cout << "posnu is " << posnu[0] << " " << posnu[1] << " " << posnu[2] << "\n";
      cout << "Length of chord is : "<<chord<<endl;
    } //end if

    altitude=where.Mag()-Geoid(lat); // what is the altitude of the entrance point
      
    if(altitude>surface_elevation+0.1) // if it is above the surface, it's messed up
      cout << "neutrino entrance point is above the surface.  Event is " << inu << "\n";

    while(altitude>MIN_ALTITUDE_CRUST && x<posnu.Distance(earth_in)) { // starting at earth entrance point, step toward interaction position until you've reached the interaction or you are below the crust.
      //    while(altitude>MIN_ALTITUDE_CRUST && x<dDistance(enterice,earth_in)) {
      
      // find the density of the earth at this altitude
      ddensity=Signal::RHOAIR;
      if (altitude<=surface_elevation+0.1 && altitude>(surface_elevation-local_icethickness)) // the 0.1 is just to take care of precision issues.   It could have been 0.01 or 0.001.
	ddensity=icedensityarray[ilon][ilat]*1000;
      else if (altitude<=(surface_elevation-local_icethickness) && altitude>(surface_elevation-local_icethickness-local_waterdepth))
	ddensity=waterdensityarray[ilon][ilat]*1000;
      else if (altitude<=(surface_elevation-local_icethickness-local_waterdepth) && altitude>softsedr[ilon][ilat]) {
	ddensity=softseddensityarray[ilon][ilat]*1000;
	crust_entered=1; //Switch that lets us know we've penetrated into the crust
      } //end if
      else if (altitude<=softsedr[ilon][ilat] && altitude>hardsedr[ilon][ilat]) {
	ddensity=hardseddensityarray[ilon][ilat]*1000;
	crust_entered=1; //Switch that lets us know we've penetrated into the crust
      } //end if
      else if (altitude<=hardsedr[ilon][ilat] && altitude>uppercrustr[ilon][ilat]) {
	ddensity=uppercrustdensityarray[ilon][ilat]*1000;
	crust_entered=1; //Switch that lets us know we've penetrated into the crust
      } //end if
      else if (altitude<=uppercrustr[ilon][ilat] && altitude>middlecrustr[ilon][ilat]) {
	ddensity=middlecrustdensityarray[ilon][ilat]*1000;
	crust_entered=1; //Switch that lets us know we've penetrated into the crust
      } //end if
      else if (altitude<=middlecrustr[ilon][ilat] && altitude>lowercrustr[ilon][ilat]) {
	ddensity=lowercrustdensityarray[ilon][ilat]*1000;
	crust_entered=1; //Switch that lets us know we've penetrated into the crust
      } //end if
      else if (altitude<=lowercrustr[ilon][ilat]) {
	ddensity=densities[1];
	crust_entered=1; //Switch that lets us know we've penetrated into the crust
      } //end if

      // sometimes altitude will not satisfy any of these because it is above the surface.
      // the neutrino is skimming the surface and will fly through the air for a while.

      L=len_int_kgm2/ddensity; // get the interaction length for that density
      weight1*=exp(-step/L);  // adjust the weight accordingly
      total_kgm2+=ddensity*step; //increase column density accordingly
      if (exp(-step/L) > 1)
	cout<<"Oops! len_int_kgm2, ddensity, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<exp(-step/L)<<endl;
      x+=step; // distance you have stepped through the earth so far.

      where += step*nchord;// find where you are now along the neutrino's path 

      lon = where.Lon();
      lat = where.Lat();
      ilon = (int)(lon/2);
      ilat = (int)(lat/2);
      altitude=where.Mag()-Geoid(lat); //what is the altitude
      surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
      local_icethickness = this->IceThickness(lon,lat);
      local_waterdepth = WaterDepth(lon,lat);
      
    } //end while

    if (x>posnu.Distance(earth_in) && weightabsorption) // if you left the loop because you have already stepped the whole distance from the entrance point to the neutrino interaction position
      return 1;

    // if you left the loop because you're not in the crust anymore
    if (altitude<=MIN_ALTITUDE_CRUST) {
      
      mantle_entered=1; //Switch that lets us know we're into the mantle

      // determine if it enters the core or not.
      sinth=sin(where.Angle(nchord));
      costh=sqrt(1-sinth*sinth);

      if (sinth>sqrt(radii[0]/radii[1])) { // it does not enter the core, just the mantle
	
	nearthlayers++;  // count the mantle as a layer traversed.
       
	L=len_int_kgm2/densities[1]; // interaction length in the mantle
	halfchord=sqrt(radii[1])*costh; // 1/2 chord the neutrino goes through in the mantle
	
	weight1 *= exp(-2*halfchord/L); // adjust the weight for path through mantle
	total_kgm2+= 2*halfchord*densities[1];  //add column density for path through mantle 
	where += (2*halfchord)*nchord; // neutrino's new position once it reaches the crust again
	
      } //end if (not entering core)
      // these enter the core
      else {
	core_entered=1; //Switch that lets us know we've entered the core

	nearthlayers+=2; // count the mantle and core as a layer traversed.
	       
	L=len_int_kgm2/densities[1]; // interaction length in mantle
	
	// compute distance before entering the core.
	halfchord=sqrt(radii[0]-radii[1]*sinth*sinth); // find distance it travels in the mantle
	distance=sqrt(radii[1])*costh-halfchord;
	weight1*=exp(-distance/L); // adjust the weight
	total_kgm2 += 2*distance*densities[1]; //Add column density for trip through mantle
	// go through the core.
	L=len_int_kgm2/densities[0]; // interaction length in core
	weight1*=exp(-2*halfchord/L); // adjust the weight
	total_kgm2 += 2*halfchord*densities[0]; //Add column density for trip through core
	// go through 2nd layer again.
	L=len_int_kgm2/densities[1];
	weight1*=exp(-distance/L);
 	if (exp(-distance/L) > 1)
	  cout<<"Oops2! len_int_kgm2, ddensity, distance, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<distance<<" , "<<exp(-distance/L)<<endl;
	where += (2*distance+2*halfchord)*nchord;  // neutrino's new position once it reaches the crust again
	
      } //end else(enter core)
    } //end if(left crust)

    lon = where.Lon();
    lat = where.Lat();
    ilon = (int)(lon/2);
    ilat = (int)(lat/2);
    altitude=where.Mag()-Geoid(lat); //what is the altitude
    surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
    local_icethickness = this->IceThickness(lon,lat);
    local_waterdepth = WaterDepth(lon,lat);

    double distance_remaining=where.Distance(posnu); // how much farther you need to travel before you reach the neutrino interaction point
    
    x=0; // this keeps track of how far you've stepped along the neutrino path, starting at the crust entrance.
    while(x<=distance_remaining) { // keep going until you have reached the interaction position

      double ddensity=Signal::RHOAIR;

      // which layer does it go through
      if (altitude<=surface_elevation && altitude>(surface_elevation-local_icethickness))
	ddensity=icedensityarray[ilon][ilat]*1000;
      if (altitude<=(surface_elevation-local_icethickness) && altitude>(surface_elevation-local_icethickness-local_waterdepth))
	ddensity=waterdensityarray[ilon][ilat]*1000;
      if (altitude<=(surface_elevation-local_icethickness-local_waterdepth) && altitude>softsedr[ilon][ilat])
	ddensity=softseddensityarray[ilon][ilat]*1000;
      if (altitude<=softsedr[ilon][ilat] && altitude>hardsedr[ilon][ilat])
	ddensity=hardseddensityarray[ilon][ilat]*1000;
      if (altitude<=hardsedr[ilon][ilat] && altitude>uppercrustr[ilon][ilat])
	ddensity=uppercrustdensityarray[ilon][ilat]*1000;
      if (altitude<=uppercrustr[ilon][ilat] && altitude>middlecrustr[ilon][ilat])
	ddensity=middlecrustdensityarray[ilon][ilat]*1000;
      if (altitude<=middlecrustr[ilon][ilat] && altitude>lowercrustr[ilon][ilat])
	ddensity=lowercrustdensityarray[ilon][ilat]*1000;
      if (altitude<=lowercrustr[ilon][ilat])
	ddensity=3.4*1000;

      L=len_int_kgm2/ddensity; // get the interaction length for that density
      weight1*=exp(-step/L);  // adjust the weight accordingly
      total_kgm2 += step*ddensity;

      if (exp(-step/L) > 1)
	  cout<<"Oops3! len_int_kgm2, ddensity, step, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<step<<" , "<<exp(-step/L)<<endl;
      x+=step; // increment how far you've stepped through crust


      // possible for a neutrino to go through the air but not likely because they aren't the most extreme skimmers (they went through the mantle)
      where += step*nchord; // where you are now along neutrino's path
          
      lon = where.Lon();
      lat = where.Lat();
      ilon = (int)(lon/2);
      ilat = (int)(lat/2);
      altitude=where.Mag()-Geoid(lat); //what is the altitude
      surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
      local_icethickness = this->IceThickness(lon,lat);
      local_waterdepth = WaterDepth(lon,lat);
    } //while
   
  } //if (getchord_method == 2)


  if (weightabsorption==0) {
    if (Rand3.Rndm()>weight1) { 
   
      weight1=0.;
      return 0;
    }
    else {
     
      weight1=1.;
      return 1;
    }
  }
  else 
    return 1;
      
  cout << "made it this far.\n";

  return 1;
} //end Getchord





// new Getchord from icemc (org version)
// new tau weight calculation, probability to interact inside the ice
int EarthModel::Getchord(Primaries *primary1, Settings *settings1,IceModel *antarctica1, Secondaries *sec1,
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
			 double *mytausurvarray, double& tauweight, double& tauchord, double *avgdensityarray, double *densityarray)  {
    
    Vector chord3;
    Vector nchord;
    double x=0;
    double lat,lon;
    int ilon,ilat;
    
    
    total_kgm2 = 0; //Initialize column density
    nearthlayers=0; // this counts the number of earth layers the neutrino traverses.
    // Want to find probability that the neutrino survives its trip
    // through the earth.
    
    //Find the chord, its length and its unit vector.
    chord3 = posnu - earth_in;
    chord=chord3.Mag();
    nchord = chord3 / chord;
    
    if (chord<=1) {
	cout << "short chord " << chord << "\n";
	return 0;
    }
    if (chord>2.*R_EARTH+1000) {
	cout << "bad chord" << " " << chord << ".  Event is " << inu << "\n";
    }
    
    Position where=earth_in;
    //cout <<"where(1) is "<<where;
    // the sin of the angle between the neutrino path and the 
    // radial vector to its earth entrance point determines
    // if it will get to the next layer down.
    double costh=(where*nchord)/where.Mag();
    double sinth=sqrt(1-costh*costh);
    double distance=0;
    double halfchord=0;
    double taumodes = settings1->taumodes;


    if (getchord_method<1 || getchord_method>3)
	cout << "Bogus method!\n";
    if (thisnuflavor =="nutau" && taumodes1==1 ){
      // cout <<"nuflavor is nutau, \n";
      sec1->GetTauWeight(primary1, settings1,antarctica1,pnu, nu_nubar,currentint,Etau_final,posnu, earth_in,
				       crust_entered, mantle_entered, core_entered, myxarray, myEarray, myyweightarray, 
				       mytausurvarray,tauchord,avgdensityarray,densityarray,inu,weight1_tmp,probability_tmp);
      weight1_tmp *=2; //only going to get half the normal number, so multiply by 2 to compensate.
      tauweight = weight1_tmp;
      //cout <<"weight1_tmp(tau) is "<<weight1_tmp<<".\n";
      return 1;
    }//thisnuflavor
    //cout <<"nuflavor is not nutau. \n";
    // we are really focusing on method 2 - method 1 has not been maintenanced in a long time!!
    // use at your own risk.
    if (getchord_method==1) {
	double L=0;
	weight1_tmp=0;
	
	if (sinth>sqrt(radii[1]/radii[2])) {
	    nearthlayers++;
	    
	    // these only skim the first layer.
	    L=len_int_kgm2/densities[2];
	    
	    weight1_tmp=exp(-posnu.Distance(where)/L);
	}
	else {
	    nearthlayers++;
	    
	    // these get to the second layer down.
	    L=len_int_kgm2/densities[2];
	    // compute distance before it gets to the next layer.
	    halfchord=sqrt(radii[1]-radii[2]*sinth*sinth);
	    distance=sqrt(radii[2])*costh-halfchord;
	    
	    weight1_tmp=exp(-distance/L);
	    
	    // get position where it enters the second layer.
	    where = earth_in + distance*nchord;
	    
	    // determine if it enters the core or not.
	    costh=(where*nchord)/where.Mag();
	    sinth=sqrt(1-costh*costh);
	    
	    if (sinth>sqrt(radii[0]/radii[1])) {
		
		
		halfchord=sqrt(radii[1])*costh;
		nearthlayers++;
		
		// these do not enter the core.
		L=len_int_kgm2/densities[1];
		
		
		weight1_tmp *= exp(-2*halfchord/L); 
		
		L=len_int_kgm2/densities[2];
		// this is where it exits the second layer and enters the crust again.
		where = where + 2*halfchord*nchord;
		weight1_tmp*=exp(-where.Distance(posnu)/L);
	    }
	    else {
		nearthlayers++;
		// these enter the core.
		L=len_int_kgm2/densities[1];
		
		// compute distance before entering the core.
		halfchord=sqrt(radii[0]-radii[1]*sinth*sinth);
		distance=sqrt(radii[1])*costh-halfchord;
		weight1_tmp*=exp(-distance/L);
		
		// go through the core.
		L=len_int_kgm2/densities[0];
		weight1_tmp*=exp(-2*halfchord/L);
		
		// go through 2nd layer again.
		L=len_int_kgm2/densities[1];
		weight1_tmp*=exp(-distance/L);
		
		// through the crust and end up at posnu.
		L=len_int_kgm2/densities[2];
		
		where = where + (2*distance+2*halfchord)*nchord;
		weight1_tmp*=exp(-where.Distance(posnu)/L);
	    } //else
	} //else
    } //if getchord_method==1
    if (getchord_method==2) {
	
	x=0; // x is the distance you move as you step through the earth.
	
	lon = where.Lon();
	lat = where.Lat();
	ilon = (int)(lon/2);
	ilat = (int)(lat/2);
	
	double surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	
	double local_icethickness = this->IceThickness(lon,lat);
	double local_waterdepth = WaterDepth(lon,lat);
	double altitude=0;
	weight1_tmp=1;
	probability_tmp=1;
	double step=Tools::dMin(len_int_kgm2/densities[1]/10,500.); //how big is the step size
	// either 1/10 of an interaction length in the mantle or 500 m, whichever is smaller.
	// 500 m is approximately the feature size in Crust 2.0.
	//------------------added on Dec 8------------------------
	weight1_tmp*=exp(-myair/len_int_kgm2);//add atmosphere attenuation // fenfang's atten. due to atmosphere
	//------------------added on Dec 8------------------------
	total_kgm2+=myair;
	
	
	
	double L_ice=len_int_kgm2/Signal::RHOICE;
	
	if (settings1->UNBIASED_SELECTION)
	    probability_tmp*=1.-exp(-1.*(r_enterice.Distance(nuexitice)/L_ice)); // probability it interacts in ice along its path
	
	double L=0;
	
	double ddensity=Signal::RHOAIR;
	nearthlayers=1;
	
	if (where*nchord>0.)  { // look at direction of neutrino where it enters the earth.
	    cout << "This one's trouble.  Neutrino exit point looks more like an entrance point.  Event is " << inu << "\n";
	    cout << "where is " << where[0] << " " << where[1] << " " << where[2] << "\n";
	    cout << "nchord is " << nchord[0] << " " << nchord[1] << " " << nchord[2] << "\n";
	    cout << "dot product is " << where*nchord/sqrt(where*where) << "\n";
	    cout << "posnu is " << posnu[0] << " " << posnu[1] << " " << posnu[2] << "\n";
	    cout << "Length of chord is : "<<chord<<endl;
	} //end if
	
	altitude=where.Mag()-Geoid(lat); // what is the altitude of the entrance point
	
	if(altitude>surface_elevation+0.1) // if it is above the surface, it's messed up
	    cout << "neutrino entrance point is above the surface.  Event is " << inu << "\n";
	//cout <<"altitude is "<<altitude<<"\n";
	
	while(altitude>MIN_ALTITUDE_CRUST && x<posnu.Distance(earth_in)) { // starting at earth entrance point, step toward interaction position until you've reached the interaction or you are below the crust.
	    //    while(altitude>MIN_ALTITUDE_CRUST && x<dDistance(enterice,earth_in)) {
	  double abovesurface =0;
	  ddensity = this->GetDensity(altitude,where,posnu,crust_entered,mantle_entered,core_entered);
	  //  if (abovesurface==1)
	  //  cout <<"altitude is "<<altitude<<"\n";
	    // sometimes altitude will not satisfy any of these because it is above the surface.
	    // the neutrino is skimming the surface and will fly through the air for a while.
	    
	    L=len_int_kgm2/ddensity; // get the interaction length for that density
	    weight1_tmp*=exp(-step/L);  // adjust the weight accordingly
	    
	    
	    total_kgm2+=ddensity*step; //increase column density accordingly
	    if (exp(-step/L) > 1)
		cout<<"Oops! len_int_kgm2, ddensity, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<exp(-step/L)<<endl;
	    x+=step; // distance you have stepped through the earth so far.
	    
	    where += step*nchord;// find where you are now along the neutrino's path 
	    
	    lon = where.Lon();
	    lat = where.Lat();
	    ilon = (int)(lon/2);
	    ilat = (int)(lat/2);
	    altitude=where.Mag()-Geoid(lat); //what is the altitude
	    surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	    local_icethickness = this->IceThickness(lon,lat);
	    local_waterdepth = WaterDepth(lon,lat);
	    
	} //end while
	
	if (x>posnu.Distance(earth_in) && weightabsorption) { // if you left the loop because you have already stepped the whole distance from the entrance point to the neutrino interaction position
	  if (taumodes1==1 && thisnuflavor =="nutau")
	    weight1_tmp *=2; //factor of two for only having half the regular neutrinos.
	    probability_tmp*=weight1_tmp;
	    return 1;
	}
	// if you left the loop because you're not in the crust anymore
	if (altitude<=MIN_ALTITUDE_CRUST) {
	    
	    mantle_entered=1; //Switch that lets us know we're into the mantle
	    
	    // determine if it enters the core or not.
	    sinth=sin(where.Angle(nchord));
	    costh=sqrt(1-sinth*sinth);
	    
	    if (sinth>sqrt(radii[0]/radii[1])) { // it does not enter the core, just the mantle
		
		nearthlayers++;  // count the mantle as a layer traversed.
		
		L=len_int_kgm2/densities[1]; // interaction length in the mantle
		halfchord=sqrt(radii[1])*costh; // 1/2 chord the neutrino goes through in the mantle
		
		weight1_tmp *= exp(-2*halfchord/L); // adjust the weight for path through mantle
		total_kgm2+= 2*halfchord*densities[1];  //add column density for path through mantle 
		where += (2*halfchord)*nchord; // neutrino's new position once it reaches the crust again
		
	    } //end if (not entering core)
	    // these enter the core
	    else {
		core_entered=1; //Switch that lets us know we've entered the core
		
		nearthlayers+=2; // count the mantle and core as a layer traversed.
		
		L=len_int_kgm2/densities[1]; // interaction length in mantle
		
		// compute distance before entering the core.
		halfchord=sqrt(radii[0]-radii[1]*sinth*sinth); // find distance it travels in the mantle
		distance=sqrt(radii[1])*costh-halfchord;
		weight1_tmp*=exp(-distance/L); // adjust the weight
		total_kgm2 += 2*distance*densities[1]; //Add column density for trip through mantle
		// go through the core.
		L=len_int_kgm2/densities[0]; // interaction length in core
		weight1_tmp*=exp(-2*halfchord/L); // adjust the weight
		total_kgm2 += 2*halfchord*densities[0]; //Add column density for trip through core
		// go through 2nd layer again.
		L=len_int_kgm2/densities[1];
		weight1_tmp*=exp(-distance/L);
		if (exp(-distance/L) > 1)
		    cout<<"Oops2! len_int_kgm2, ddensity, distance, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<distance<<" , "<<exp(-distance/L)<<endl;
		where += (2*distance+2*halfchord)*nchord;  // neutrino's new position once it reaches the crust again
		
	    } //end else(enter core)
	} //end if(left crust)
	
	lon = where.Lon();
	lat = where.Lat();
	ilon = (int)(lon/2);
	ilat = (int)(lat/2);
	altitude=where.Mag()-Geoid(lat); //what is the altitude
	surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	local_icethickness = this->IceThickness(lon,lat);
	local_waterdepth = WaterDepth(lon,lat);
	
	double distance_remaining=where.Distance(posnu); // how much farther you need to travel before you reach the neutrino interaction point
	
	x=0; // this keeps track of how far you've stepped along the neutrino path, starting at the crust entrance.
	while(x<=distance_remaining) { // keep going until you have reached the interaction position
	    
	  ddensity=this->GetDensity(altitude,where,posnu,crust_entered,mantle_entered,core_entered);
	  
	    L=len_int_kgm2/ddensity; // get the interaction length for that density
	    weight1_tmp*=exp(-step/L);  // adjust the weight accordingly
	    total_kgm2 += step*ddensity;
	    
	    if (exp(-step/L) > 1)
		cout<<"Oops3! len_int_kgm2, ddensity, step, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<step<<" , "<<exp(-step/L)<<endl;
	    x+=step; // increment how far you've stepped through crust
	    
	    
	    // possible for a neutrino to go through the air but not likely because they aren't the most extreme skimmers (they went through the mantle)
	    where += step*nchord; // where you are now along neutrino's path
	    
	    lon = where.Lon();
	    lat = where.Lat();
	    ilon = (int)(lon/2);
	    ilat = (int)(lat/2);
	    altitude=where.Mag()-Geoid(lat); //what is the altitude
	    surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	    local_icethickness = this->IceThickness(lon,lat);
	    local_waterdepth = WaterDepth(lon,lat);
	} //while
	
    } //if (getchord_method == 2)
    if (taumodes==1&& thisnuflavor =="nutau")
      weight1_tmp *=2; //compensation factor
    probability_tmp*=weight1_tmp;
    //cout <<"probability_tmp(non-tau) is "<<probability_tmp<<".\n";
    
    if (weightabsorption==0) {
	if (gRandom->Rndm()>weight1_tmp) { 
	    
	    weight1_tmp=0.;
	    return 0;
	}
	else {
	    
	    weight1_tmp=1.;
	    return 1;
	}
    }
    else 
	return 1;
    
    cout << "made it this far.\n";
    
    return 1;
} //end Getchord




// new Getchord from icemc (modified version)
// new tau weight calculation, probability to interact inside the ice
int EarthModel::Getchord(Primaries *primary1, Settings *settings1,IceModel *antarctica1, Secondaries *sec1,
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
			 )  {
    
    Vector chord3;
    Vector nchord;
    double x=0;
    double lat,lon;
    int ilon,ilat;
    
    
    total_kgm2 = 0; //Initialize column density
    nearthlayers=0; // this counts the number of earth layers the neutrino traverses.
    // Want to find probability that the neutrino survives its trip
    // through the earth.
    
    //Find the chord, its length and its unit vector.
    chord3 = posnu - earth_in;
    chord=chord3.Mag();
    nchord = chord3 / chord;
    
    if (chord<=1) {
	cout << "short chord " << chord << "\n";
	return 0;
    }
    if (chord>2.*R_EARTH+1000) {
	cout << "bad chord" << " " << chord << ".  Event is " << inu << "\n";
    }
    
    Position where=earth_in;
    //cout <<"where(1) is "<<where;
    // the sin of the angle between the neutrino path and the 
    // radial vector to its earth entrance point determines
    // if it will get to the next layer down.
    double costh=(where*nchord)/where.Mag();
    double sinth=sqrt(1-costh*costh);
    double distance=0;
    double halfchord=0;
    double taumodes = settings1->taumodes;


    if (getchord_method<1 || getchord_method>3)
	cout << "Bogus method!\n";
    if (thisnuflavor =="nutau" && taumodes1==1 ){
      // cout <<"nuflavor is nutau, \n";
      /*
      sec1->GetTauWeight(primary1, settings1,antarctica1,pnu, nu_nubar,currentint,Etau_final,posnu, earth_in,
				       crust_entered, mantle_entered, core_entered, myxarray, myEarray, myyweightarray, 
				       mytausurvarray,tauchord,avgdensityarray,densityarray,inu,weight1_tmp,probability_tmp);
      */
      weight1_tmp=0.; // initialization necessary (added by K.M. 2013.10.31)

      sec1->GetTauWeight(primary1, settings1,antarctica1,pnu, nu_nubar,currentint,Etau_final,posnu, earth_in,
				       crust_entered, mantle_entered, core_entered, 
				       weight1_tmp,probability_tmp);

      weight1_tmp *=2; //only going to get half the normal number, so multiply by 2 to compensate.
      //tauweight = weight1_tmp;
      //cout <<"weight1_tmp(tau) is "<<weight1_tmp<<".\n";
      return 1;
    }//thisnuflavor
    //cout <<"nuflavor is not nutau. \n";
    // we are really focusing on method 2 - method 1 has not been maintenanced in a long time!!
    // use at your own risk.
    if (getchord_method==1) {
	double L=0;
	weight1_tmp=0;
	
	if (sinth>sqrt(radii[1]/radii[2])) {
	    nearthlayers++;
	    
	    // these only skim the first layer.
	    L=len_int_kgm2/densities[2];
	    
	    weight1_tmp=exp(-posnu.Distance(where)/L);
	}
	else {
	    nearthlayers++;
	    
	    // these get to the second layer down.
	    L=len_int_kgm2/densities[2];
	    // compute distance before it gets to the next layer.
	    halfchord=sqrt(radii[1]-radii[2]*sinth*sinth);
	    distance=sqrt(radii[2])*costh-halfchord;
	    
	    weight1_tmp=exp(-distance/L);
	    
	    // get position where it enters the second layer.
	    where = earth_in + distance*nchord;
	    
	    // determine if it enters the core or not.
	    costh=(where*nchord)/where.Mag();
	    sinth=sqrt(1-costh*costh);
	    
	    if (sinth>sqrt(radii[0]/radii[1])) {
		
		
		halfchord=sqrt(radii[1])*costh;
		nearthlayers++;
		
		// these do not enter the core.
		L=len_int_kgm2/densities[1];
		
		
		weight1_tmp *= exp(-2*halfchord/L); 
		
		L=len_int_kgm2/densities[2];
		// this is where it exits the second layer and enters the crust again.
		where = where + 2*halfchord*nchord;
		weight1_tmp*=exp(-where.Distance(posnu)/L);
	    }
	    else {
		nearthlayers++;
		// these enter the core.
		L=len_int_kgm2/densities[1];
		
		// compute distance before entering the core.
		halfchord=sqrt(radii[0]-radii[1]*sinth*sinth);
		distance=sqrt(radii[1])*costh-halfchord;
		weight1_tmp*=exp(-distance/L);
		
		// go through the core.
		L=len_int_kgm2/densities[0];
		weight1_tmp*=exp(-2*halfchord/L);
		
		// go through 2nd layer again.
		L=len_int_kgm2/densities[1];
		weight1_tmp*=exp(-distance/L);
		
		// through the crust and end up at posnu.
		L=len_int_kgm2/densities[2];
		
		where = where + (2*distance+2*halfchord)*nchord;
		weight1_tmp*=exp(-where.Distance(posnu)/L);
	    } //else
	} //else
    } //if getchord_method==1
    if (getchord_method==2) {
	
	x=0; // x is the distance you move as you step through the earth.
	
	lon = where.Lon();
	lat = where.Lat();
	ilon = (int)(lon/2);
	ilat = (int)(lat/2);
	
	double surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	
	double local_icethickness = this->IceThickness(lon,lat);
	double local_waterdepth = WaterDepth(lon,lat);
	double altitude=0;
	weight1_tmp=1;
	probability_tmp=1;
	double step=Tools::dMin(len_int_kgm2/densities[1]/10,500.); //how big is the step size
	// either 1/10 of an interaction length in the mantle or 500 m, whichever is smaller.
	// 500 m is approximately the feature size in Crust 2.0.
	//------------------added on Dec 8------------------------
	weight1_tmp*=exp(-myair/len_int_kgm2);//add atmosphere attenuation // fenfang's atten. due to atmosphere
	//------------------added on Dec 8------------------------
	total_kgm2+=myair;
	
	
	
	double L_ice=len_int_kgm2/Signal::RHOICE;
	
	if (settings1->UNBIASED_SELECTION)
	    probability_tmp*=1.-exp(-1.*(r_enterice.Distance(nuexitice)/L_ice)); // probability it interacts in ice along its path
	
	double L=0;
	
	double ddensity=Signal::RHOAIR;
	nearthlayers=1;
	
	if (where*nchord>0.)  { // look at direction of neutrino where it enters the earth.
	    cout << "This one's trouble.  Neutrino exit point looks more like an entrance point.  Event is " << inu << "\n";
	    cout << "where is " << where[0] << " " << where[1] << " " << where[2] << "\n";
	    cout << "nchord is " << nchord[0] << " " << nchord[1] << " " << nchord[2] << "\n";
	    cout << "dot product is " << where*nchord/sqrt(where*where) << "\n";
	    cout << "posnu is " << posnu[0] << " " << posnu[1] << " " << posnu[2] << "\n";
	    cout << "Length of chord is : "<<chord<<endl;
	} //end if
	
	altitude=where.Mag()-Geoid(lat); // what is the altitude of the entrance point
	
	if(altitude>surface_elevation+0.1) // if it is above the surface, it's messed up
	    cout << "neutrino entrance point is above the surface.  Event is " << inu << "\n";
	//cout <<"altitude is "<<altitude<<"\n";
	
	while(altitude>MIN_ALTITUDE_CRUST && x<posnu.Distance(earth_in)) { // starting at earth entrance point, step toward interaction position until you've reached the interaction or you are below the crust.
	    //    while(altitude>MIN_ALTITUDE_CRUST && x<dDistance(enterice,earth_in)) {
	  double abovesurface =0;
	  ddensity = this->GetDensity(altitude,where,posnu,crust_entered,mantle_entered,core_entered);
	  //  if (abovesurface==1)
	  //  cout <<"altitude is "<<altitude<<"\n";
	    // sometimes altitude will not satisfy any of these because it is above the surface.
	    // the neutrino is skimming the surface and will fly through the air for a while.
	    
	    L=len_int_kgm2/ddensity; // get the interaction length for that density
	    weight1_tmp*=exp(-step/L);  // adjust the weight accordingly
	    
	    
	    total_kgm2+=ddensity*step; //increase column density accordingly
	    if (exp(-step/L) > 1)
		cout<<"Oops! len_int_kgm2, ddensity, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<exp(-step/L)<<endl;
	    x+=step; // distance you have stepped through the earth so far.
	    
	    where += step*nchord;// find where you are now along the neutrino's path 
	    
	    lon = where.Lon();
	    lat = where.Lat();
	    ilon = (int)(lon/2);
	    ilat = (int)(lat/2);
	    altitude=where.Mag()-Geoid(lat); //what is the altitude
	    surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	    local_icethickness = this->IceThickness(lon,lat);
	    local_waterdepth = WaterDepth(lon,lat);
	    
	} //end while
	
	if (x>posnu.Distance(earth_in) && weightabsorption) { // if you left the loop because you have already stepped the whole distance from the entrance point to the neutrino interaction position
	  if (taumodes1==1 && thisnuflavor =="nutau")
	    weight1_tmp *=2; //factor of two for only having half the regular neutrinos.
	    probability_tmp*=weight1_tmp;
	    return 1;
	}
	// if you left the loop because you're not in the crust anymore
	if (altitude<=MIN_ALTITUDE_CRUST) {
	    
	    mantle_entered=1; //Switch that lets us know we're into the mantle
	    
	    // determine if it enters the core or not.
	    sinth=sin(where.Angle(nchord));
	    costh=sqrt(1-sinth*sinth);
	    
	    if (sinth>sqrt(radii[0]/radii[1])) { // it does not enter the core, just the mantle
		
		nearthlayers++;  // count the mantle as a layer traversed.
		
		L=len_int_kgm2/densities[1]; // interaction length in the mantle
		halfchord=sqrt(radii[1])*costh; // 1/2 chord the neutrino goes through in the mantle
		
		weight1_tmp *= exp(-2*halfchord/L); // adjust the weight for path through mantle
		total_kgm2+= 2*halfchord*densities[1];  //add column density for path through mantle 
		where += (2*halfchord)*nchord; // neutrino's new position once it reaches the crust again
		
	    } //end if (not entering core)
	    // these enter the core
	    else {
		core_entered=1; //Switch that lets us know we've entered the core
		
		nearthlayers+=2; // count the mantle and core as a layer traversed.
		
		L=len_int_kgm2/densities[1]; // interaction length in mantle
		
		// compute distance before entering the core.
		halfchord=sqrt(radii[0]-radii[1]*sinth*sinth); // find distance it travels in the mantle
		distance=sqrt(radii[1])*costh-halfchord;
		weight1_tmp*=exp(-distance/L); // adjust the weight
		total_kgm2 += 2*distance*densities[1]; //Add column density for trip through mantle
		// go through the core.
		L=len_int_kgm2/densities[0]; // interaction length in core
		weight1_tmp*=exp(-2*halfchord/L); // adjust the weight
		total_kgm2 += 2*halfchord*densities[0]; //Add column density for trip through core
		// go through 2nd layer again.
		L=len_int_kgm2/densities[1];
		weight1_tmp*=exp(-distance/L);
		if (exp(-distance/L) > 1)
		    cout<<"Oops2! len_int_kgm2, ddensity, distance, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<distance<<" , "<<exp(-distance/L)<<endl;
		where += (2*distance+2*halfchord)*nchord;  // neutrino's new position once it reaches the crust again
		
	    } //end else(enter core)
	} //end if(left crust)
	
	lon = where.Lon();
	lat = where.Lat();
	ilon = (int)(lon/2);
	ilat = (int)(lat/2);
	altitude=where.Mag()-Geoid(lat); //what is the altitude
	surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	local_icethickness = this->IceThickness(lon,lat);
	local_waterdepth = WaterDepth(lon,lat);
	
	double distance_remaining=where.Distance(posnu); // how much farther you need to travel before you reach the neutrino interaction point
	
	x=0; // this keeps track of how far you've stepped along the neutrino path, starting at the crust entrance.
	while(x<=distance_remaining) { // keep going until you have reached the interaction position
	    
	  ddensity=this->GetDensity(altitude,where,posnu,crust_entered,mantle_entered,core_entered);
	  
	    L=len_int_kgm2/ddensity; // get the interaction length for that density
	    weight1_tmp*=exp(-step/L);  // adjust the weight accordingly
	    total_kgm2 += step*ddensity;
	    
	    if (exp(-step/L) > 1)
		cout<<"Oops3! len_int_kgm2, ddensity, step, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<step<<" , "<<exp(-step/L)<<endl;
	    x+=step; // increment how far you've stepped through crust
	    
	    
	    // possible for a neutrino to go through the air but not likely because they aren't the most extreme skimmers (they went through the mantle)
	    where += step*nchord; // where you are now along neutrino's path
	    
	    lon = where.Lon();
	    lat = where.Lat();
	    ilon = (int)(lon/2);
	    ilat = (int)(lat/2);
	    altitude=where.Mag()-Geoid(lat); //what is the altitude
	    surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	    local_icethickness = this->IceThickness(lon,lat);
	    local_waterdepth = WaterDepth(lon,lat);
	} //while
	
    } //if (getchord_method == 2)
    if (taumodes==1&& thisnuflavor =="nutau")
      weight1_tmp *=2; //compensation factor
    probability_tmp*=weight1_tmp;
    //cout <<"probability_tmp(non-tau) is "<<probability_tmp<<".\n";
    
    if (weightabsorption==0) {
	if (gRandom->Rndm()>weight1_tmp) { 
	    
	    weight1_tmp=0.;
	    return 0;
	}
	else {
	    
	    weight1_tmp=1.;
	    return 1;
	}
    }
    else 
	return 1;
    
    cout << "made it this far.\n";
    
    return 1;
} //end Getchord








// new Getchord from icemc (modified version)
// new tau weight calculation, probability to interact inside the ice
// additional L0 input to account the probability to interact inside the sphere (INTERACTION_MODE=0)
int EarthModel::Getchord(Primaries *primary1, Settings *settings1,IceModel *antarctica1, Secondaries *sec1,
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
			 int nu_nubar, int currentint,int taumodes1, double L0
			 )  {
    
    Vector chord3;
    Vector nchord;
    double x=0;
    double lat,lon;
    int ilon,ilat;
    
    
    total_kgm2 = 0; //Initialize column density
    nearthlayers=0; // this counts the number of earth layers the neutrino traverses.
    // Want to find probability that the neutrino survives its trip
    // through the earth.
    
    //Find the chord, its length and its unit vector.
    chord3 = posnu - earth_in;
    chord=chord3.Mag();
    nchord = chord3 / chord;
    
    if (chord<=1) {
	cout << "short chord " << chord << "\n";
	return 0;
    }
    if (chord>2.*R_EARTH+1000) {
	cout << "bad chord" << " " << chord << ".  Event is " << inu << "\n";
    }
    
    Position where=earth_in;
    //cout <<"where(1) is "<<where;
    // the sin of the angle between the neutrino path and the 
    // radial vector to its earth entrance point determines
    // if it will get to the next layer down.
    double costh=(where*nchord)/where.Mag();
    double sinth=sqrt(1-costh*costh);
    double distance=0;
    double halfchord=0;
    double taumodes = settings1->taumodes;


    if (getchord_method<1 || getchord_method>3)
	cout << "Bogus method!\n";
    if (thisnuflavor =="nutau" && taumodes1==1 ){
      // cout <<"nuflavor is nutau, \n";
      /*
      sec1->GetTauWeight(primary1, settings1,antarctica1,pnu, nu_nubar,currentint,Etau_final,posnu, earth_in,
				       crust_entered, mantle_entered, core_entered, myxarray, myEarray, myyweightarray, 
				       mytausurvarray,tauchord,avgdensityarray,densityarray,inu,weight1_tmp,probability_tmp);
      */
      weight1_tmp=0.; // initialization necessary (added by K.M. 2013.10.31)

      //cout <<"weight1_tmp(tau before) is "<<weight1_tmp<<".\n";
      sec1->GetTauWeight(primary1, settings1,antarctica1,pnu, nu_nubar,currentint,Etau_final,posnu, earth_in,
				       crust_entered, mantle_entered, core_entered, 
				       weight1_tmp,probability_tmp);

      weight1_tmp *=2; //only going to get half the normal number, so multiply by 2 to compensate.
      //tauweight = weight1_tmp;
      //cout <<"weight1_tmp(tau after) is "<<weight1_tmp<<".\n";
      return 1;
    }//thisnuflavor
    //cout <<"nuflavor is not nutau. \n";
    // we are really focusing on method 2 - method 1 has not been maintenanced in a long time!!
    // use at your own risk.
    if (getchord_method==1) {
	double L=0;
	weight1_tmp=0;
	
	if (sinth>sqrt(radii[1]/radii[2])) {
	    nearthlayers++;
	    
	    // these only skim the first layer.
	    L=len_int_kgm2/densities[2];
	    
	    weight1_tmp=exp(-posnu.Distance(where)/L);
	}
	else {
	    nearthlayers++;
	    
	    // these get to the second layer down.
	    L=len_int_kgm2/densities[2];
	    // compute distance before it gets to the next layer.
	    halfchord=sqrt(radii[1]-radii[2]*sinth*sinth);
	    distance=sqrt(radii[2])*costh-halfchord;
	    
	    weight1_tmp=exp(-distance/L);
	    
	    // get position where it enters the second layer.
	    where = earth_in + distance*nchord;
	    
	    // determine if it enters the core or not.
	    costh=(where*nchord)/where.Mag();
	    sinth=sqrt(1-costh*costh);
	    
	    if (sinth>sqrt(radii[0]/radii[1])) {
		
		
		halfchord=sqrt(radii[1])*costh;
		nearthlayers++;
		
		// these do not enter the core.
		L=len_int_kgm2/densities[1];
		
		
		weight1_tmp *= exp(-2*halfchord/L); 
		
		L=len_int_kgm2/densities[2];
		// this is where it exits the second layer and enters the crust again.
		where = where + 2*halfchord*nchord;
		weight1_tmp*=exp(-where.Distance(posnu)/L);
	    }
	    else {
		nearthlayers++;
		// these enter the core.
		L=len_int_kgm2/densities[1];
		
		// compute distance before entering the core.
		halfchord=sqrt(radii[0]-radii[1]*sinth*sinth);
		distance=sqrt(radii[1])*costh-halfchord;
		weight1_tmp*=exp(-distance/L);
		
		// go through the core.
		L=len_int_kgm2/densities[0];
		weight1_tmp*=exp(-2*halfchord/L);
		
		// go through 2nd layer again.
		L=len_int_kgm2/densities[1];
		weight1_tmp*=exp(-distance/L);
		
		// through the crust and end up at posnu.
		L=len_int_kgm2/densities[2];
		
		where = where + (2*distance+2*halfchord)*nchord;
		weight1_tmp*=exp(-where.Distance(posnu)/L);
	    } //else
	} //else
    } //if getchord_method==1
    if (getchord_method==2) {
	
	x=0; // x is the distance you move as you step through the earth.
	
	lon = where.Lon();
	lat = where.Lat();
	ilon = (int)(lon/2);
	ilat = (int)(lat/2);
	
	double surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	
	double local_icethickness = this->IceThickness(lon,lat);
	double local_waterdepth = WaterDepth(lon,lat);
	double altitude=0;
	weight1_tmp=1;
	probability_tmp=1;
	double step=Tools::dMin(len_int_kgm2/densities[1]/10,500.); //how big is the step size
	// either 1/10 of an interaction length in the mantle or 500 m, whichever is smaller.
	// 500 m is approximately the feature size in Crust 2.0.
	//------------------added on Dec 8------------------------
	weight1_tmp*=exp(-myair/len_int_kgm2);//add atmosphere attenuation // fenfang's atten. due to atmosphere
	//------------------added on Dec 8------------------------
	total_kgm2+=myair;
	
	
	
	double L_ice=len_int_kgm2/Signal::RHOICE;
	
	if (settings1->UNBIASED_SELECTION){

	  if(settings1->INTERACTION_MODE!=0){
	    probability_tmp*=1.-exp(-1.*(r_enterice.Distance(nuexitice)/L_ice)); // probability it interacts in ice along its path  
	  }
          else{  // if INTERACTION_MODE==0 (using sphere area)
	    probability_tmp*=1.-exp(-1.*(L0/L_ice)); // probability it interacts in ice along its path

	    //cout << "Yeah!!!" << endl;
	    //cout << L0 << " " << r_enterice.Distance(nuexitice) 
		 //<< " " << L_ice << endl;
	  }
	}
	
	double L=0;
	
	double ddensity=Signal::RHOAIR;
	nearthlayers=1;
	
	if (where*nchord>0.)  { // look at direction of neutrino where it enters the earth.
	    cout << "This one's trouble.  Neutrino exit point looks more like an entrance point.  Event is " << inu << "\n";
	    cout << "where is " << where[0] << " " << where[1] << " " << where[2] << "\n";
	    cout << "nchord is " << nchord[0] << " " << nchord[1] << " " << nchord[2] << "\n";
	    cout << "dot product is " << where*nchord/sqrt(where*where) << "\n";
	    cout << "posnu is " << posnu[0] << " " << posnu[1] << " " << posnu[2] << "\n";
	    cout << "Length of chord is : "<<chord<<endl;
	} //end if
	
	altitude=where.Mag()-Geoid(lat); // what is the altitude of the entrance point
	
	if(altitude>surface_elevation+0.1) // if it is above the surface, it's messed up
	    cout << "neutrino entrance point is above the surface.  Event is " << inu << "\n";
	//cout <<"altitude is "<<altitude<<"\n";
	
	while(altitude>MIN_ALTITUDE_CRUST && x<posnu.Distance(earth_in)) { // starting at earth entrance point, step toward interaction position until you've reached the interaction or you are below the crust.
	    //    while(altitude>MIN_ALTITUDE_CRUST && x<dDistance(enterice,earth_in)) {
	  double abovesurface =0;
	  ddensity = this->GetDensity(altitude,where,posnu,crust_entered,mantle_entered,core_entered);
	  //  if (abovesurface==1)
	  //  cout <<"altitude is "<<altitude<<"\n";
	    // sometimes altitude will not satisfy any of these because it is above the surface.
	    // the neutrino is skimming the surface and will fly through the air for a while.
	    
	    L=len_int_kgm2/ddensity; // get the interaction length for that density
	    weight1_tmp*=exp(-step/L);  // adjust the weight accordingly
	    
	    
	    total_kgm2+=ddensity*step; //increase column density accordingly
	    if (exp(-step/L) > 1)
		cout<<"Oops! len_int_kgm2, ddensity, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<exp(-step/L)<<endl;
	    x+=step; // distance you have stepped through the earth so far.
	    
	    where += step*nchord;// find where you are now along the neutrino's path 
	    
	    lon = where.Lon();
	    lat = where.Lat();
	    ilon = (int)(lon/2);
	    ilat = (int)(lat/2);
	    altitude=where.Mag()-Geoid(lat); //what is the altitude
	    surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	    local_icethickness = this->IceThickness(lon,lat);
	    local_waterdepth = WaterDepth(lon,lat);
	    
	} //end while
	
	if (x>posnu.Distance(earth_in) && weightabsorption) { // if you left the loop because you have already stepped the whole distance from the entrance point to the neutrino interaction position
	  if (taumodes1==1 && thisnuflavor =="nutau")
	    weight1_tmp *=2; //factor of two for only having half the regular neutrinos.
	    probability_tmp*=weight1_tmp;
	    return 1;
	}
	// if you left the loop because you're not in the crust anymore
	if (altitude<=MIN_ALTITUDE_CRUST) {
	    
	    mantle_entered=1; //Switch that lets us know we're into the mantle
	    
	    // determine if it enters the core or not.
	    sinth=sin(where.Angle(nchord));
	    costh=sqrt(1-sinth*sinth);
	    
	    if (sinth>sqrt(radii[0]/radii[1])) { // it does not enter the core, just the mantle
		
		nearthlayers++;  // count the mantle as a layer traversed.
		
		L=len_int_kgm2/densities[1]; // interaction length in the mantle
		halfchord=sqrt(radii[1])*costh; // 1/2 chord the neutrino goes through in the mantle
		
		weight1_tmp *= exp(-2*halfchord/L); // adjust the weight for path through mantle
		total_kgm2+= 2*halfchord*densities[1];  //add column density for path through mantle 
		where += (2*halfchord)*nchord; // neutrino's new position once it reaches the crust again
		
	    } //end if (not entering core)
	    // these enter the core
	    else {
		core_entered=1; //Switch that lets us know we've entered the core
		
		nearthlayers+=2; // count the mantle and core as a layer traversed.
		
		L=len_int_kgm2/densities[1]; // interaction length in mantle
		
		// compute distance before entering the core.
		halfchord=sqrt(radii[0]-radii[1]*sinth*sinth); // find distance it travels in the mantle
		distance=sqrt(radii[1])*costh-halfchord;
		weight1_tmp*=exp(-distance/L); // adjust the weight
		total_kgm2 += 2*distance*densities[1]; //Add column density for trip through mantle
		// go through the core.
		L=len_int_kgm2/densities[0]; // interaction length in core
		weight1_tmp*=exp(-2*halfchord/L); // adjust the weight
		total_kgm2 += 2*halfchord*densities[0]; //Add column density for trip through core
		// go through 2nd layer again.
		L=len_int_kgm2/densities[1];
		weight1_tmp*=exp(-distance/L);
		if (exp(-distance/L) > 1)
		    cout<<"Oops2! len_int_kgm2, ddensity, distance, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<distance<<" , "<<exp(-distance/L)<<endl;
		where += (2*distance+2*halfchord)*nchord;  // neutrino's new position once it reaches the crust again
		
	    } //end else(enter core)
	} //end if(left crust)
	
	lon = where.Lon();
	lat = where.Lat();
	ilon = (int)(lon/2);
	ilat = (int)(lat/2);
	altitude=where.Mag()-Geoid(lat); //what is the altitude
	surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	local_icethickness = this->IceThickness(lon,lat);
	local_waterdepth = WaterDepth(lon,lat);
	
	double distance_remaining=where.Distance(posnu); // how much farther you need to travel before you reach the neutrino interaction point
	
	x=0; // this keeps track of how far you've stepped along the neutrino path, starting at the crust entrance.
	while(x<=distance_remaining) { // keep going until you have reached the interaction position
	    
	  ddensity=this->GetDensity(altitude,where,posnu,crust_entered,mantle_entered,core_entered);
	  
	    L=len_int_kgm2/ddensity; // get the interaction length for that density
	    weight1_tmp*=exp(-step/L);  // adjust the weight accordingly
	    total_kgm2 += step*ddensity;
	    
	    if (exp(-step/L) > 1)
		cout<<"Oops3! len_int_kgm2, ddensity, step, factor : "<<len_int_kgm2<<" , "<<ddensity<<" , "<<step<<" , "<<exp(-step/L)<<endl;
	    x+=step; // increment how far you've stepped through crust
	    
	    
	    // possible for a neutrino to go through the air but not likely because they aren't the most extreme skimmers (they went through the mantle)
	    where += step*nchord; // where you are now along neutrino's path
	    
	    lon = where.Lon();
	    lat = where.Lat();
	    ilon = (int)(lon/2);
	    ilat = (int)(lat/2);
	    altitude=where.Mag()-Geoid(lat); //what is the altitude
	    surface_elevation = this->SurfaceAboveGeoid(lon,lat); // altitude of surface relative to geoid at earth entrance point
	    local_icethickness = this->IceThickness(lon,lat);
	    local_waterdepth = WaterDepth(lon,lat);
	} //while
	
    } //if (getchord_method == 2)
    if (taumodes==1&& thisnuflavor =="nutau")
      weight1_tmp *=2; //compensation factor
    probability_tmp*=weight1_tmp;
    //cout <<"probability_tmp(non-tau) is "<<probability_tmp<<".\n";
    
    if (weightabsorption==0) {
	if (gRandom->Rndm()>weight1_tmp) { 
	    
	    weight1_tmp=0.;
	    return 0;
	}
	else {
	    
	    weight1_tmp=1.;
	    return 1;
	}
    }
    else 
	return 1;
    
    cout << "made it this far.\n";
    
    return 1;
} //end Getchord










 Vector EarthModel::GetSurfaceNormal(const Position &r_out) const
{
  Vector n_surf = r_out.Unit();
  if (FLATSURFACE)
    return n_surf;

  double theta=r_out.Theta();
  
  int ilon,ilat;
  GetILonILat(r_out,ilon,ilat);
  
  int ilon_previous=ilon-1;
  if (ilon_previous<0)
    ilon_previous=NLON-1;
  
  int ilon_next=ilon+1;
  if (ilon_next==NLON)
    ilon_next=0;
  
  double r=(geoid[ilat]+surfacer[ilon][ilat])*sin(theta);
  
  double slope_phi=(surfacer[ilon_next][ilat]-surfacer[ilon_previous][ilat])/(r*2*phistep);
  
  int ilat_previous=ilat-1;
  if (ilat_previous<0)
    ilat_previous=0;
  
  int ilat_next=ilat+1;
  if (ilat_next==NLAT)
    ilat_next=NLAT-1;
     
  double slope_costheta=(surfacer[ilon][ilat_next]-surfacer[ilon][ilat_previous])/((geoid[ilat]+surfacer[ilon][ilat])*2*thetastep);
  
  // first rotate n_surf according to tilt in costheta and position on continent - rotate around the y axis.
  double angle=atan(slope_costheta);
  
  n_surf = n_surf.RotateY(angle);
  
  // now rotate n_surf according to tilt in phi - rotate around the z axis.
  angle=atan(slope_phi);
  
  n_surf = n_surf.RotateZ(angle);
  
  return n_surf;
    
} //method GetSurfaceNormal

 double EarthModel::SmearPhi(int ilon)  {


  double phi=((double)(360.*3./4.-((double)ilon+Rand3.Rndm())*360/180))*RADDEG;
  if (phi<0 && phi>-1*PI/2)
    phi+=2*PI;

  
  return phi;
} //SmearPhi

 double EarthModel::SmearTheta(int ilat)  {

  // remember that we should smear it evenly in cos(theta).
  // first get the cos(theta)'s at the boundaries.
 
  double theta1=dGetTheta(ilat)-PI/(double)NLAT/2.;
  double theta2=dGetTheta(ilat+1)-PI/(double)NLAT/2.;
 
  double costheta1=cos(theta1);
  double costheta2=cos(theta2);



  double costheta=Rand3.Rndm()*(costheta2-costheta1)+costheta1;

  double theta=acos(costheta);

  return theta;
} //SmearTheta

void EarthModel::ReadCrust(string test) {

  // reads in altitudes of 7 layers of crust, ice and water
  // puts data in arrays
  
  fstream infile(test.c_str(),ios::in);

  string thisline; // for reading in file
  string slon; //longitude as a string
  string slat; // latitude as a string
  string selev; // elevation (km relative to geoid)
  string sdepth; // depth (km)
  string sdensity; // density (g/cm^3)
  double dlon,dlat; // longitude, latitude as double
  int endindex; // index along thisline for parsing
  int beginindex; // same

  int indexlon=0; // 180 bins in longitude
  int indexlat=0; // 90 bins in latitude
   
  string layertype; // water, ice, etc.

  while(!infile.eof()) {
    getline(infile,thisline,'\n'); 
    
    int loc=thisline.find("type, latitude, longitude,"); 
    
    if (loc!=(int)(string::npos)) {      
      
      beginindex=thisline.find_first_not_of(" ",57);
      
      endindex=thisline.find_first_of(" ",61);
      
      slat=thisline.substr(beginindex,endindex-beginindex);
      dlat=(double)atof(slat.c_str());

      beginindex=thisline.find_first_not_of(" ",68);
      endindex=thisline.find_first_of(" ",72);

      slon=thisline.substr(beginindex,endindex-beginindex);
      dlon=(double)atof(slon.c_str());
      

      indexlon=(int)((dlon+180)/2);
      indexlat=(int)((90+dlat)/2);

      beginindex=thisline.find_first_not_of(" ",78);
      endindex=thisline.find_first_of(" ",83);

      selev=thisline.substr(beginindex,endindex-beginindex);
      elevationarray[indexlon][indexlat]=(double)atof(selev.c_str());

    } //if

    for (int i=0;i<4;i++) {
      getline(infile,thisline,'\n');
    } //for
       
    for (int i=0;i<7;i++) {
      getline(infile,thisline,'\n');
      
      endindex=thisline.length()-1;
      beginindex=thisline.find_last_of("0123456789",1000);
      layertype=thisline.substr(beginindex+3,endindex-beginindex);

     
      beginindex=thisline.find_first_not_of(" ",0);
      endindex=thisline.find_first_of(" ",beginindex);
     
      sdepth=thisline.substr(beginindex,endindex-beginindex-1);
      

      // fills arrays of thicknesses of each layer
      if (layertype.substr(0,5)=="water") 
	waterthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str()); 
      if (layertype.substr(0,3)=="ice") 
	icethkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
      if (layertype.substr(0,8)=="soft sed") 
	softsedthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
      if (layertype.substr(0,8)=="hard sed") 
	hardsedthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
      if (layertype.substr(0,11)=="upper crust") 
	uppercrustthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
      if (layertype.substr(0,12)=="middle crust") 
	middlecrustthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());
      if (layertype.substr(0,11)=="lower crust") 
	lowercrustthkarray[indexlon][indexlat]=(double)atof(sdepth.c_str());

      //      cout << "indexlon, indexlat, icethkarray " << indexlon << " " << indexlat << " " << icethkarray[indexlon][indexlat] << "\n";

      // region where Ross Ice Shelf was not accounted for in Crust 2.0
      // add it in by hand
      if (indexlat==5 && (indexlon<=5 || indexlon>=176)) // Ross Ice Shelf
	icethkarray[indexlon][indexlat]=0.5;

      beginindex=thisline.find_first_not_of(" ",endindex);
      endindex=thisline.find_first_of(" ",beginindex);
      

      beginindex=thisline.find_first_not_of(" ",endindex);
      endindex=thisline.find_first_of(" ",beginindex);
      
      beginindex=thisline.find_first_not_of(" ",endindex);
      endindex=thisline.find_first_of(" ",beginindex);

     
      sdensity=thisline.substr(beginindex,endindex-beginindex);

      double ddensity=(double)atof(sdensity.c_str());
      

      // fills arrays of densities of each layer
      if (layertype.substr(0,5)=="water") 
	waterdensityarray[indexlon][indexlat]=ddensity; 
      if (layertype.substr(0,3)=="ice") 
	icedensityarray[indexlon][indexlat]=ddensity;
      if (layertype.substr(0,8)=="soft sed") 
	softseddensityarray[indexlon][indexlat]=ddensity;
      if (layertype.substr(0,8)=="hard sed") 
	hardseddensityarray[indexlon][indexlat]=ddensity;
      if (layertype.substr(0,11)=="upper crust") 
	uppercrustdensityarray[indexlon][indexlat]=ddensity;
      if (layertype.substr(0,12)=="middle crust")
	middlecrustdensityarray[indexlon][indexlat]=ddensity;
      if (layertype.substr(0,11)=="lower crust") 
	lowercrustdensityarray[indexlon][indexlat]=ddensity;
    } //for (reading all lines for one location given in Crust 2.0 input file)

    if (CONSTANTCRUST) {
      softsedthkarray[indexlon][indexlat]=40.;
      hardsedthkarray[indexlon][indexlat]=0;
      uppercrustthkarray[indexlon][indexlat]=0;
      middlecrustthkarray[indexlon][indexlat]=0;
      lowercrustthkarray[indexlon][indexlat]=0;
      crustthkarray[indexlon][indexlat]=0;
      softseddensityarray[indexlon][indexlat]=2.9;
    } //if (set crust thickness to constant everywhere)
    if (CONSTANTICETHICKNESS) {
      icethkarray[indexlon][indexlat]=3.;
      waterthkarray[indexlon][indexlat]=0.;
    } //if (set ice thickness to constant everywhere)

    // adds up total thickness of crust
    crustthkarray[indexlon][indexlat]=softsedthkarray[indexlon][indexlat]+
      hardsedthkarray[indexlon][indexlat]+
      uppercrustthkarray[indexlon][indexlat]+
      middlecrustthkarray[indexlon][indexlat]+
      lowercrustthkarray[indexlon][indexlat];
    
    if (indexlon==179 && indexlat==0)
      break;
  }  // done reading file
  
  for (int i=0;i<NLON;i++) {
    for (int j=0;j<NLAT;j++) {

      if (FIXEDELEVATION) 	
	elevationarray[i][j]=icethkarray[i][j]*1000;

      if (waterthkarray[i][j] != 0) 
	surfacer[i][j]=elevationarray[i][j]+waterthkarray[i][j]*1000+icethkarray[i][j]*1000;	  
      else
	surfacer[i][j]=elevationarray[i][j];

      if (fabs(surfacer[i][j])<1.E-10)
	surfacer[i][j] = 0;	

      // reminder- waterr is elevation at *bottom* of water layer, etc. 
      // in units of m
      waterr[i][j]=surfacer[i][j]-(icethkarray[i][j]+waterthkarray[i][j])*1000;
      if ((double)fabs(waterr[i][j])<1.E-10)
	waterr[i][j]=0;
      icer[i][j]=waterr[i][j]+
	waterthkarray[i][j]*1000;
      softsedr[i][j]=waterr[i][j]-
	softsedthkarray[i][j]*1000;
      hardsedr[i][j]=waterr[i][j]-
	(softsedthkarray[i][j]+
	 hardsedthkarray[i][j])*1000;
      uppercrustr[i][j]=waterr[i][j]-
	(softsedthkarray[i][j]+
	 hardsedthkarray[i][j]+
	 uppercrustthkarray[i][j])*1000;
      middlecrustr[i][j]=waterr[i][j]-
	(softsedthkarray[i][j]+
	 hardsedthkarray[i][j]+
	 uppercrustthkarray[i][j]+
	 middlecrustthkarray[i][j])*1000;
      lowercrustr[i][j]=waterr[i][j]-
	(softsedthkarray[i][j]+
	 hardsedthkarray[i][j]+
	 uppercrustthkarray[i][j]+
	 middlecrustthkarray[i][j]+
	 lowercrustthkarray[i][j])*1000;
    } //for (latitude bins)
  } //for (longitude bins)

  // calculate ice volume for comparison with expectation
  volume=0;
  double sumarea=0; // sum of surface area of ice
  average_iceth = 0;

  for (int j=0;j<ILAT_MAX;j++) {
   for (int i=0;i<NLON;i++) {
     volume+=icethkarray[i][j]*1000.*area[j];

     /*
      // j=6 corresponds to 80deg S
      if (j==6) {
	// fill output file, just for plotting
	outfile << surfacer[i][j] << "\t" << waterr[i][j] << "\t" << icer[i][j] << "\t" << icethkarray[i][j] << " " << waterr[i][j]-icer[i][j] << " " << (waterr[i][j]-icer[i][j])*area[j] << " " << area[j] << " " << volume << "\n";
      }//if
     */

      // find average ice thickness
      average_iceth+=(surfacer[i][j]-icer[i][j])*area[j];
      sumarea+=area[j];
   } //for
  } //for
  average_iceth=average_iceth/sumarea; 

  // find the place where the crust is the deepest.
  // for finding where to start stepping in Getchord
  MIN_ALTITUDE_CRUST=1.E6;
  //MAX_VOL=-1.E6;
  for (int i=0;i<NLON;i++) {
    for (int j=0;j<NLAT;j++) {
      if (elevationarray[i][j]-(crustthkarray[i][j])*1000<MIN_ALTITUDE_CRUST) {
	if (waterthkarray[i][j]==0)
	  MIN_ALTITUDE_CRUST=elevationarray[i][j]-(crustthkarray[i][j]+icethkarray[i][j])*1000;
	else
	  MIN_ALTITUDE_CRUST=elevationarray[i][j]-crustthkarray[i][j]*1000;
      }//if
      //if (icethkarray[i][j]*1000.*area[j]>MAX_VOL) 
      //MAX_VOL=icethkarray[i][j]*1000.*area[j];      
    }//for
  }//for
  
  //record depth of crust-mantle interface
  radii[1]=(Geoid(0.)+MIN_ALTITUDE_CRUST)*(Geoid(0.)+MIN_ALTITUDE_CRUST);

}//ReadCrust

 double EarthModel::dGetTheta(int ilat) const {
  return (((double)ilat+0.5)/(double)NLAT*MAXTHETA)*RADDEG;
} //dGetTheta(int)

 double EarthModel::dGetPhi(int ilon) const {
  // this takes as an input the crust 2.0 index 0=-180 deg longitude to 179=+180 deg longitude
  // its output is phi in radians
  // from ~ -pi/2 to 3*pi/2 
  return (double)(-1*((double)ilon+0.5)+(double)NLON)*2*PI/(double)NLON-PI/2;
} //dGetPhi(int)

 void EarthModel::GetILonILat(const Position &p,int& ilon,int& ilat) const {
  // Phi function outputs from 0 to 2*pi wrt +x
  double phi_deg=p.Phi()*DEGRAD;

  if (phi_deg>270)
    phi_deg=phi_deg-360;
  // now it's from -90 to 270
  
  ilon=(int)((360.*3./4.-phi_deg)*180./360.); // ilon is from 0 (at -180 longitude) to 180 (at 180 longitude)
  
  ilat=(int)((p.Theta()*DEGRAD)/2.);
  
} //method GetILonILat
void EarthModel::EarthCurvature(double *array,double depth_temp) {

  Position parray;
  parray.SetXYZ(array[0],array[1],array[2]);

  // adjust array coordinates so that it fits to a curved earth surface at a specific depth 
  double length=Surface(parray)-depth_temp; // length=distance from center of earth

  double rxposx=array[0]; // x coordinate of antenna position
  double rxposy=array[1]; // y coordinate of antenna position
  double rxdr=sqrt(rxposx*rxposx+rxposy*rxposy); // distance in horizontal plane from the center of the detector to the antenna
  if (Tools::dSquare(array)==0) cout << "Attempt of divide by zero in Earth curvature!!\n";
  double rxdtheta=asin(rxdr/sqrt(Tools::dSquare(array)));
  double rxdphi=atan2(rxposy,rxposx);

  array[0]=length*sin(rxdtheta)*cos(rxdphi);// have the array sit on a sphere of radius "length"
  array[1]=length*sin(rxdtheta)*sin(rxdphi);
  array[2]=length*cos(rxdtheta);

}
