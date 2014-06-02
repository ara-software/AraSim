#include "TRandom3.h"
#include "Constants.h"
#include "Primaries.h"

#include "IceModel.h"
#include "EarthModel.h"
#include "Vector.h"
#include "Ray.h"
#include "Settings.h"
//#include "icemodel.hh"
//#include "earthmodel.hh"
//#include "vector.hh"
//#include "ray.hh"

#include "Tools.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

ClassImp(IceModel);

using namespace std;

//class Interaction;



//Parameters of the BEDMAP ice model. (See http://www.antarctica.ac.uk/aedc/bedmap/download/)
int nCols_ice=1200; //number of columns in data, set by header file (should be 1200)
int nRows_ice=1000; //number of rows in data, set by header file (should be 1000)
int cellSize=5000; //in meters, set by header file (should be 5000) - same for both files
int xLowerLeft_ice=-3000000; 
int yLowerLeft_ice=-2500000;
int nCols_ground=1068;
int nRows_ground=869;
int xLowerLeft_ground=-2661600;
int yLowerLeft_ground=-2149967;
int nCols_water=1200;
int nRows_water=1000;
int xLowerLeft_water=-3000000;
int yLowerLeft_water=-2500000;
int NODATA=-9999;

//Variables for conversion between BEDMAP polar stereographic coordinates and lat/lon.  Conversion equations from ftp://164.214.2.65/pub/gig/tm8358.2/TM8358_2.pdf
const double scale_factor=0.97276901289;  //scale factor at pole corresponding to 71 deg S latitude of true scale (used in BEDMAP)
const double ellipsoid_inv_f = 298.257223563; //of Earth
const double ellipsoid_b = EarthModel::R_EARTH*(1-(1/ellipsoid_inv_f));
const double eccentricity = sqrt((1/ellipsoid_inv_f)*(2-(1/ellipsoid_inv_f)));
const double bedmap_a_bar = pow(eccentricity,2)/2 + 5*pow(eccentricity,4)/24 + pow(eccentricity,6)/12 + 13*pow(eccentricity,8)/360;
const double bedmap_b_bar = 7*pow(eccentricity,4)/48 + 29*pow(eccentricity,6)/240 + 811*pow(eccentricity,8)/11520;
const double bedmap_c_bar = 7*pow(eccentricity,6)/120 + 81*pow(eccentricity,8)/1120;
const double bedmap_d_bar = 4279*pow(eccentricity,8)/161280;
const double bedmap_c_0 = (2*EarthModel::R_EARTH / sqrt(1-pow(eccentricity,2))) * pow(( (1-eccentricity) / (1+eccentricity) ),eccentricity/2);
double bedmap_R = scale_factor*bedmap_c_0 * pow(( (1 + eccentricity*sin(71*RADDEG)) / (1 - eccentricity*sin(71*RADDEG)) ),eccentricity/2) * tan((PI/4) - (71*RADDEG)/2); //varies with latitude, defined here for 71 deg S latitude
const double bedmap_nu = bedmap_R / cos(71*RADDEG);

/*
IceModel::IceModel() {
    //default constructor
}
*/

IceModel::IceModel(int model,int earth_model,int moorebay) : EarthModel(earth_model),mooreBayFlag(moorebay) {

  
  setUpIceModel(model);

 


 }


IceModel::~IceModel () {

}


void IceModel::setUpIceModel(int model) {
  
  DEPTH_DEPENDENT_N = (int) (model / 10);
  model -= DEPTH_DEPENDENT_N * 10;
  ice_model=model;

  if (ice_model != 0 && ice_model != 1) {
    cout<<"Error!  Unknown ice model requested!  Defaulting to Crust 2.0.\n";
    ice_model = 0;
  } //if
  else if (ice_model==1) {
    ReadIceThickness();
    ReadWaterDepth();
    ReadGroundBed();
  } //else if (BEDMAP)
 //read in attenuation length data for direct signals
  int i=0;
  ifstream sheetup("data/icesheet_attenlength_up.txt");
  if(sheetup.fail())
    {
      cerr << "Failed to open icesheet_attenlength_up.txt" << endl;
      exit(1);
    }
  
  i=0;
  while(sheetup>>d_sheetup[i]>>l_sheetup[i])
    {
      i++;
    }
  sheetup.close();
   
  ifstream shelfup("data/iceshelf_attenlength_up.txt");
  if(shelfup.fail())
    {
      cerr << "Failed to open iceshelf_attenlength_up.txt" << endl;
      exit(1);
    }
  
  i=0;
  while(shelfup>>d_shelfup[i]>>l_shelfup[i])
    {
      i++;
    }
  shelfup.close();
   
  ifstream westlandup("data/westland_attenlength_up.txt");
  if(westlandup.fail())
    {cerr << "Failed to open westland_attenlength_up.txt";
      exit(1);
    }
  i=0;
  while(westlandup>>d_westlandup[i]>>l_westlandup[i])
    {
      i++;
    }
  westlandup.close();
   
  //read in attenuation length for downgoing signals
  ifstream sheetdown("data/icesheet_attenlength_down.txt");
  if(sheetdown.fail())
    {
      cerr << "Failed to open icesheet_attenlength_down.txt" << endl;
      exit(1);
    }
  
  i=0;
  while(sheetdown>>d_sheetdown[i]>>l_sheetdown[i])
    {
      i++;
    }
  sheetdown.close();

  
  ifstream shelfdown("data/iceshelf_attenlength_down.txt");
  if(shelfdown.fail())
    {
      cerr << "Failed to open iceshelf_attenlength_down.txt" << endl;
      exit(1);
    }
  
  i=0;
  while(shelfdown>>d_shelfdown[i]>>l_shelfdown[i])
    {
      i++;
    }
  shelfdown.close();
   
  ifstream westlanddown("data/westland_attenlength_down.txt");
  if(westlanddown.fail())
    {cerr << "Failed to open westland_attenlength_down.txt";
      exit(1);
    }
  i=0;
  while(westlanddown>>d_westlanddown[i]>>l_westlanddown[i])
    {
      i++;
    }
  westlanddown.close();


  // new ARA ice attenuation measurement values (at 300 MHz)
  //
  /*
  ifstream file( "data/ARA_IceAttenL.txt" );

  string line;
  string line2;

  int N=-1;

  int skipline = 1;
  int first_time = 1;


if ( file.is_open() ) {
    while (file.good() ) {
        
        if ( first_time == 1 ) {
            for (int sl=0; sl<skipline; sl++) {
                getline (file, line);
            }
            first_time = 0;
        }
                

        getline (file, line);

        N++;


        //ARA_IceAtten_Depth[N] = atof( line.substr(0, line.find_first_of(",")).c_str() );
        cout<< line.substr(0, line.find_first_of(",")).c_str()<<endl;

        cout<<"ARA_IceAtten Depth"<<N<<" : "<<ARA_IceAtten_Depth[N]<<"\t";

        line2 = line.substr( line.find_first_of(",")+1);

        ARA_IceAtten_Length[N] = atof( line2.substr(0).c_str() );

        cout<<"ARA_IceAtten L"<<N<<" : "<<ARA_IceAtten_Length[N]<<"\n";

    }
    file.close();
}

    ARA_IceAtten_bin = N;
    cout<<"ARA_IceAtten total N : "<<ARA_IceAtten_bin<<"\n";

    // done reading ARA ice attenuation info
    */

  //
  double ARA_IceAtten_Depth_tmp[53] = { 72.7412,   76.5697,    80.3982,    91.8836,    95.7121,    107.198,    118.683,    133.997,    153.139,    179.939,    206.738,    245.023,    298.622,    356.049,    405.819,    470.904,    516.845,    566.616,    616.386,    669.985,    727.412,    784.839,    838.438,    899.694,    949.464,    1003.06,    1060.49,    1121.75,    1179.17,    1236.6,    1297.86,    1347.63,    1405.05,    1466.31,    1516.08,    1565.85,    1611.79,    1657.73,    1699.85,    1745.79,    1791.73,    1833.84,    1883.61,    1929.56,    1990.81,    2052.07,    2109.49,    2170.75,    2232.01,    2304.75,    2362.17,    2431.09,    2496.17 };

  double ARA_IceAtten_Length_tmp[53] = { 1994.67,   1952,    1896,    1842.67,    1797.33,    1733.33,    1680,    1632,    1586.67,    1552,    1522.67,    1501.33,    1474.67,    1458.67,    1437.33,    1416,    1392,    1365.33,    1344,    1312,    1274.67,    1242.67,    1205.33,    1168,    1128,    1090.67,    1048,    1008,    965.333,    920,    874.667,    834.667,    797.333,    752,    714.667,    677.333,    648,    616,    589.333,    557.333,    530.667,    506.667,    477.333,    453.333,    418.667,    389.333,    362.667,    333.333,    309.333,    285.333,    264,    242.667,    221.333 };

  ARA_IceAtten_bin = 53;
  for (int bin=0; bin<ARA_IceAtten_bin; bin++) {

      ARA_IceAtten_Depth[bin] = ARA_IceAtten_Depth_tmp[bin];
      ARA_IceAtten_Length[bin] = ARA_IceAtten_Length_tmp[bin];
  }


 
}



// read depth in positive value and return attenuation length (m) at the depth
double IceModel::GetARAIceAttenuLength(double depth) {

    double AttenL;

    // check if depth is positive value
    if ( depth < 0. ) {// whether above the ice or wrong value!

        cerr<<"depth negative! "<<depth<<endl;
    }
    else {

        AttenL = Tools::SimpleLinearInterpolation_extend_Single(ARA_IceAtten_bin, ARA_IceAtten_Depth, ARA_IceAtten_Length, depth );
    }

    return AttenL;

}




int IceModel::Getice_model() {
    return ice_model;
}

 //constructor IceModel(int model)
Position IceModel::PickBalloonPosition() const {
  Vector temp;
  return temp;

}


//--------------------------------------------------
// int IceModel::PickUnbiased(int inu, Interaction *interaction1, IceModel *antarctica) {
//     
// 
//   interaction1->PickAnyDirection(); // first pick the neutrino direction
// 
//   double mincos=cos(COASTLINE*RADDEG);
//   double maxcos=cos(0.);
//   double minphi=0.;
//   double maxphi=2.*PI;
//   double thisphi,thiscos,thissin;
//   double theta=0.;
//   double phi=0.;
// 
//   int ilon,ilat;    
//   int e_coord,n_coord;
//   double vol_thisbin=0.;
//   double lon=0.;
//   double lat=0.;
//   
//  
//     thisphi=gRandom->Rndm()*(maxphi-minphi)+minphi;
//     thiscos=gRandom->Rndm()*(maxcos-mincos)+mincos;
//     thissin=sqrt(1.-thiscos*thiscos);
//     Position thisr_in;// entrance point
//     Position thisr_enterice;
// Position thisr_enterice_tmp;
//     Position thisnuexitearth;
//     Position thisnuexitice;
//     Position thisr_exitice;
//     interaction1->noway=0;
//     interaction1->wheredoesitleave_err=0;
//     interaction1->neverseesice=0;
//     interaction1->wheredoesitenterice_err=0;
//     interaction1->toohigh=0;
//     interaction1->toolow=0;
// 
//     thisr_in.SetXYZ(R_EARTH*thissin*cos(thisphi),R_EARTH*thissin*sin(thisphi),R_EARTH*thiscos);
//     if (thisr_in.Dot(interaction1->nnu)>0)
//       interaction1->nnu=-1.*interaction1->nnu;
//     // does this intersect any ice
//     //cout << "lat, coastline, cos are " << thisr_in.Lat() << " " << COASTLINE << " " << cos(interaction1->nnu.Theta()) << "\n";
//     if (thisr_in.Lat()>COASTLINE && cos(interaction1->nnu.Theta())<0) {
//       interaction1->noway=1;
// 
//       interaction1->pickunbiased=0;
//       return 0; // there is no way it's going through the ice
//     }
// 
//     int count1=0;
//     int count2=0;
// 
//    
//     if (Ray::WhereDoesItLeave(0,thisr_in,interaction1->nnu,antarctica,thisnuexitearth)) { // where does it leave Earth
//       // really want to find where it leaves ice
//       int err;
//       // Does it leave in an ice bin
//       if (IceThickness(thisnuexitearth) && thisnuexitearth.Lat()<COASTLINE) { // if this is an ice bin in the Antarctic
// 	//cout << "inu is " << inu << " it's in ice.\n";
// 	//cout << "this is an ice bin.\n";
// 	thisnuexitice=thisnuexitearth;
// 	thisr_exitice=thisnuexitearth;
// 	if (thisnuexitice.Mag()>Surface(thisnuexitice)) { // if the exit point is above the surface
// 	  if ((thisnuexitice.Mag()-Surface(thisnuexitice))/cos(interaction1->nnu.Theta())>5.E3) { 
// 	    WhereDoesItExitIce(inu,thisnuexitearth,interaction1->nnu,5.E3, // then back up and find it more precisely
// 			       thisr_exitice);
// 	    thisnuexitice=(5000.)*interaction1->nnu;
// 	    thisnuexitice+=thisr_exitice;
// 	    count1++;
// 	  }
// 	  if ((thisnuexitice.Mag()-Surface(thisnuexitice))/cos(interaction1->nnu.Theta())>5.E2) {
// 	    
// 	    WhereDoesItExitIce(inu,thisnuexitice,interaction1->nnu,5.E2, // then back up and find it more precisely
// 			       thisr_exitice);
// 	    thisnuexitice=5.E2*interaction1->nnu;
// 	    thisnuexitice+=thisr_exitice;
// 	    count1++;
// 	  }
// 	  if ((thisnuexitice.Mag()-Surface(thisnuexitice))/cos(interaction1->nnu.Theta())>50.) {
// 
// 	    WhereDoesItExitIce(inu,thisnuexitice,interaction1->nnu,50., // then back up and find it more precisely
// 			     thisr_exitice);
// 	    count1++;
// 	  } // end third wheredoesitexit
// 	  thisnuexitice=thisr_exitice;
// 	} // if the exit point overshoots
// 	else
// 	  thisnuexitice=thisnuexitearth;
// 
// 	// should also correct for undershooting
//     if (count1>10)
//       cout << "count1 is " << count1 << "\n";	  
//       } // if it's an Antarctic ice bin
//       else { // it leaves a rock bin so back up and find where it leaves ice
// 	//cout << "inu is " << inu << " it's in rock.\n";
// 	if (thisr_in.Distance(thisnuexitearth)>5.E4) {
// 	  count2++;
// 	  if (WhereDoesItExitIce(inu,thisnuexitearth,interaction1->nnu,5.E4, // then back up and find it more precisely
// 				 thisr_exitice)) {
// 	    
// 	    thisnuexitice=(5.E4)*interaction1->nnu;
// 	    thisnuexitice+=thisr_exitice;
// 	    //cout << "inu is " << inu << " I'm here 1.\n";
// 
// 	  }
// 	  else {
// 	    interaction1->neverseesice=1;
//             interaction1->pickunbiased = 0;
// 	    return 0;
// 	  }
// 	}
// 	else
// 	  thisnuexitice=thisnuexitearth;
// 	//   WhereDoesItExitIce(inu,thisnuexit,interaction1->nnu,5.E4, // then back up and find it more precisely
// // 			     thisr_exitice);
// // 	  thisnuexit=5.E4*interaction1->nnu;
// // 	  thisnuexit+=thisr_exitice;
// 	if (thisr_in.Distance(thisnuexitice)>5.E3) {
// 
// 	  
// 	  if (WhereDoesItExitIce(inu,thisnuexitice,interaction1->nnu,5.E3, // then back up and find it more precisely
// 				  thisr_exitice)) {
// 	    count2++;
// 	    //interaction1->neverseesice=1;
// 	    thisnuexitice=5.E3*interaction1->nnu;
// 	    thisnuexitice+=thisr_exitice;
// 	    //cout << "inu is " << inu << " I'm here 2\n";
// 	    //return 0;
// 	    
// 	  }
// 	}
// 	if (thisr_in.Distance(thisnuexitice)>5.E2) {
// 
// 
// 	  if (WhereDoesItExitIce(inu,thisnuexitice,interaction1->nnu,5.E2, // then back up and find it more precisely
// 				  thisr_exitice)) {
// 	    count2++;
// 	    //interaction1->neverseesice=1;
// 
// 	    thisnuexitice=5.E2*interaction1->nnu;
// 	    thisnuexitice+=thisr_exitice;
// 	    //cout << "inu is " << inu << " I'm here 3\n";
// 	    //return 0;
// 	  }
// 	
// 	}
// 	if (thisr_in.Distance(thisnuexitice)>50.) {
// 
// 
// 	  if (WhereDoesItExitIce(inu,thisnuexitice,interaction1->nnu,50., // then back up and find it more precisely
// 				  thisr_exitice)) {
// 	    //interaction1->neverseesice=1;
// 	    count2++;
// 	    //cout << "inu is " << inu << " I'm here 4\n";
// 	    //return 0;
// 	  }
// 	}
// 	thisnuexitice=thisr_exitice;
// 	if (count2>10)
// 	  cout << "count1 is " << count2 << "\n";
// 	//	else return 0;  // never reaches any ice or is it because our step is too big
//       } // if the nu leaves a rock bin
//     } // end wheredoesitleave
//     else {
//       interaction1->wheredoesitleave_err=1;
//       interaction1->pickunbiased = 0;
//       return 0;
//     }
//     // end finding where it leaves ice
// 
// // 	if (thisnuexit.Mag()<Surface(thisnuexit)) { // if the exit point is below the surface
// // 	  WhereDoesItExitIceForward(thisnuexit,interaction1->nnu,20., // then find it more finely
// // 			     thisr_exitice);
// // 	  thisnuexit=thisr_enterice;
// // 	  // then back up and find it more precisely
// // 	}
// 
//     if (WhereDoesItEnterIce(thisnuexitearth,interaction1->nnu,5.E3, // first pass with sort of course binning
// 			    thisr_enterice)) {
//       thisr_enterice_tmp=thisr_enterice+5.E3*interaction1->nnu;
//       //cout << "inu is " << inu << " thisr_enterice is ";thisr_enterice.Print();
//       if (WhereDoesItEnterIce(thisr_enterice_tmp,interaction1->nnu,20., // second pass with finer binning
// 			      thisr_enterice)) {
// 	//cout << "inu is " << inu << " thisr_enterice is ";thisr_enterice.Print();
// 	//cout << "entersice is ";thisr_enterice.Print();
// 	//cout << "thisnuexitice is ";thisnuexitice.Print();
// 	interaction1->pathlength_inice=thisr_enterice.Distance(thisnuexitice);
// 	//cout << "distance is " << distance << "\n";
// 	//cout << "inu " << inu << " thisr_enterice, thisnuexitice are ";thisr_enterice.Print();thisnuexitice.Print();
// 	interaction1->posnu=interaction1->pathlength_inice*gRandom->Rndm()*interaction1->nnu;
// 	interaction1->posnu=interaction1->posnu+thisr_enterice;
// 	//cout << "inu" << inu << " thisr_enterice, thisnuexitice are ";thisr_enterice.Print();thisnuexitice.Print();
// 	//cout << "inu " << inu << " distance is " << distance << "\n";
//       }
//     }
//     else {
//       thisr_enterice=thisr_in;
//       interaction1->wheredoesitenterice_err=1;
//       interaction1->pickunbiased = 0;
//       return 0;
//     }
//     interaction1->nuexitice=thisnuexitice;
//     interaction1->r_enterice=thisr_enterice;
//     
//     if (interaction1->posnu.Mag()-Surface(interaction1->posnu)>0) {
//       interaction1->toohigh=1;
//       //cout << "inu, toohigh is " << inu << " " << interaction1->toohigh << "\n";
//       interaction1->pickunbiased = 0;
//       return 0;
//     }
//     if (interaction1->posnu.Mag()-Surface(interaction1->posnu)+IceThickness(interaction1->posnu)<0) {
//       interaction1->toolow=1;
//       //cout << "inu, toolow is " << inu << " " << interaction1->toolow << "\n";
//       interaction1->pickunbiased = 0;
//       return 0;
//     }    
//     interaction1->pickunbiased = 1;
//     return 1;
// 
// }
//-------------------------------------------------- 


//--------------------------------------------------
// int IceModel::PickNear() {
// }
//-------------------------------------------------- 




Vector IceModel::GetSurfaceNormal(const Position &r_out) const {
  Vector n_surf = r_out.Unit();
  if (FLATSURFACE) 
    return n_surf;

  if (ice_model==0) {
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
  } //end if(Crust 2.0)
  else if (ice_model==1) {
    double dist_to_check = 7500; //check elevations at this distance north, south, east and west of event
    double lon,lat;
    double lon_prev,lon_next;
    double lat_prev,lat_next;
    lon = r_out.Lon();
    lat = r_out.Lat(); //longitude and latitude of interaction
    double local_surface_elevation = Surface(lon,lat);

    lat_next = lat + dist_to_check * (180 / (local_surface_elevation * PI)); //the latitude 7.5 km south of the interaction
    lat_prev = lat - dist_to_check * (180 / (local_surface_elevation * PI)); //the latitude 7.5 km north of the interaction

    lon_next = lon + dist_to_check * (180 / (sin(lat*RADDEG) * local_surface_elevation * PI)); 
    lon_prev = lon - dist_to_check * (180 / (sin(lat*RADDEG) * local_surface_elevation * PI)); 

    if (lat_next > 90) {
      //cout<<"lat_next is > 90"<<endl;
      lat_next = 90 - (lat_next - 90);  //if we went past the pole, set coordinates for the other side
      lon_next += 180;
      lon_prev += 180;
    } //end if
    //cout<<"lon, lat: "<<lon<<" , "<<lat<<endl;
    //correct any out of range longitudes
    if (lon_next > 360) {
      //cout<<"lon_next > 360\n";
      lon_next -= 360;
    }
    else if (lon_next < 0) {
      //cout<<"lon_next < 0\n";
      lon_next += 360;
    }
    if (lon_prev > 360) {
      //cout<<"lon_prev > 360\n";
      lon_prev -= 360;
    }
    else if (lon_prev < 0) {
      //cout << "lon_prev < 0";
      lon_prev += 360;
    }
   
    double slope_phi=(SurfaceAboveGeoid(lon_next,lat)-SurfaceAboveGeoid(lon_prev,lat))/(2*dist_to_check);

    double slope_costheta=(SurfaceAboveGeoid(lon,lat_next)-SurfaceAboveGeoid(lon,lat_prev))/(2*dist_to_check);
    
    // first rotate n_surf according to tilt in costheta - rotate around the y axis.
    double angle=atan(slope_costheta);

    n_surf = n_surf.RotateY(angle);
   
    // now rotate n_surf according to tilt in phi - rotate around the z axis.
    angle=atan(slope_phi);

    n_surf = n_surf.RotateZ(angle);
  } //end if(BEDMAP)
    
  return n_surf;
    
} //method GetSurfaceNormal

Position IceModel::WhereDoesItEnterIce(const Position &posnu,
				       const Vector &nnu,
				       double stepsize) const {
  // now get exit point...
  //   see my geometry notes.
  // parameterize the neutrino trajectory and just see where it
  // crosses the earth radius.

  Position r_enterice;
  double distance=0;
  int left_edge=0;
  Position x = posnu;
  double x2;
  
  Position x_previous = posnu;

  double x_previous2= x_previous * x_previous;
  x2=x_previous2;
  
  double lon = x.Lon(),lat = x.Lat();
  double lon_old = lon,lat_old = lat;
  double local_surface = Surface(lon,lat);
  double rock_previous2= pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
  double surface_previous2=pow(local_surface,2);

  double rock2=rock_previous2;
  double surface2=surface_previous2;


  while (distance<2*local_surface+1000) {

    distance+=stepsize;

    x -= stepsize*nnu;
    x2=x*x;

    lon = x.Lon();
    lat = x.Lat();

    if (lon!=lon_old || lat!=lat_old) {
      local_surface = Surface(lon,lat);

      if (lat>COASTLINE) 
	left_edge=1;
      
      rock2=pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
      surface2=pow(local_surface,2);    

      if (ice_model==0) {
	if ((int)(lat)==COASTLINE && rock_previous2 < x2 && surface2 > x2)
	  left_edge=1;
      } //if (Crust 2.0)
    } //if (neutrino has stepped into new lon/lat bin)

    if ((x_previous2>rock_previous2 && x2<rock2)
	|| (x_previous2<surface_previous2 && x2>surface2)
	|| left_edge) {
       
      r_enterice = x;
      // this gets you out of the loop.
      distance=3*Geoid(lat);
    } //if

    x_previous = x;
    x_previous2 = x2;

    if (lon!=lon_old || lat!=lat_old) {
      rock_previous2 = rock2;
      surface_previous2 = surface2;
      lat_old = lat;
      lon_old = lon;
    } //if

  } //while

  return r_enterice;
}//WhereDoesItEnterIce




Position IceModel::WhereDoesItEnter(const Position &posnu,const Vector &nnu) const {
    // now get neutrino entry point...
//    cout << posnu.GetX() << " : " << posnu.GetY() << " : " << posnu.GetZ() << endl;
//    cout << posnu.R() << " : " << posnu.Theta() << " : " << posnu.Phi() << endl;
//    cout << posnu.Lon() << " : " << posnu.Lat() << endl;

    double p = posnu.Mag(); // radius of interaction
    double costheta = (nnu*posnu) / p; // theta of neutrino at interaction position
    double sintheta = sqrt(1-costheta*costheta);
    
    double lon = posnu.Lon();
    double lat = posnu.Lat();
    
    double a=0; // length of chord
    
    double R = Surface(lon,lat);
    double delta = R - p; // depth of the interaction
    // if interaction occurs below surface, as it should
    
    if (delta>-0.001) {
	a=p*costheta+sqrt(R*R*costheta*costheta+2*delta*R*sintheta*sintheta); // chord length
	if (a<0) {
	    cout << "Negative chord length: " << a << "\n";
	} //end if
    } //end if (interaction below surface)  
    else if (delta<=-0.001) {
	
	//cout << "Error in interaction position.  whichray is " << whichray << "\n";
	cout << "lon, lat from WhereDoesItEnter is " << " " << lon << " " << lat << "\n";
	cout << "geoid, surface, p, surface-p are " << Geoid(lat) << " " << Surface(lon,lat) << ", " << p << " , "<<(Surface(lon,lat)-p)<<"\n";
	
    } //else if: error: interaction takes place above the surface
    
    // first approx
    Position r_in = posnu - a*nnu;

    int iter = 0;
    // now do correction 3 times
    //for (iter=0; iter<3; iter++) {
    //    delta = r_in.Mag() - Surface( r_in );
    //    r_in = r_in + (delta * nnu);
    //}
    
    delta = r_in.Mag() - Surface( r_in );
    while ( fabs(delta) >= 0.1 ) {
        r_in = r_in + (delta * nnu);
        delta = r_in.Mag() - Surface( r_in );
        iter++;
        if ( iter > 10 ) {
            //cout<<"\n r_in iteration more than 10 times!!! delta : "<<delta<<". now set r_in as before."<<endl;
            r_in = Surface( r_in ) * r_in.Unit();   // the way before
            delta = r_in.Mag() - Surface( r_in );
        }
    }

    
    //lon = r_in.Lon();
    //lat = r_in.Lat();
    
    //r_in = Surface(lon,lat) * r_in.Unit();
    
    return r_in;
} //method WhereDoesItEnter






int IceModel::WhereDoesItEnter_sphere(const Position &sphere_in, const Vector &nnu, Position &r_in ) const {

    double p = sphere_in.Mag(); // radius of interaction
    double costheta = (nnu*sphere_in) / p; // theta of neutrino at interaction position
    double sintheta = sqrt(1-costheta*costheta);
    
    double lon = sphere_in.Lon();
    double lat = sphere_in.Lat();
    
    double a=0; // length of chord
    
    double R = Surface(lon,lat);
    double delta = R - p; // depth of the interaction
    // if interaction occurs below surface, as it should
    


    // if sphere_in is inside the earth
    if (delta>-0.001) {
	//a=p*costheta+sqrt(R*R*costheta*costheta+2*delta*R*sintheta*sintheta); // chord length
	a=p*costheta + sqrt(R*R-p*p*sintheta*sintheta); // chord length
	if (a<0) {
	    cout << "Negative chord length: " << a << "\n";
	} //end if

        // first approx
        r_in = sphere_in - a*nnu;

    } //end if (sphere_in below surface)  

    // if sphere_in is outside the earth
    else if (delta<=-0.001) {
	
        // D : shortest distance between earth center and neutrino trajectory
        Position D = sphere_in + (p*costheta)*nnu;

        // neutrino pass through the earth
        if ( D.Mag() < Surface(D) ) {
        
            // first approx
            r_in = sphere_in - (sqrt(R*R-D.Mag()*D.Mag())+costheta*p)*nnu;
        }
        // neutrino don't pass through the earth
        else {
            return 0;
        }
	
    } //else if (sphere_in above surface)
    

    int iter = 0;
    // now do correction 3 times
    //for (iter=0; iter<3; iter++) {
    //    delta = r_in.Mag() - Surface( r_in );
    //    r_in = r_in + (delta * nnu);
    //}
    
    delta = r_in.Mag() - Surface( r_in );
    while ( fabs(delta) >= 0.1 ) {
        r_in = r_in + (delta * nnu);
        delta = r_in.Mag() - Surface( r_in );
        iter++;
        if ( iter > 10 ) {
            //cout<<"\n r_in iteration more than 10 times!!! delta : "<<delta<<". now set r_in as before."<<endl;
            r_in = Surface( r_in ) * r_in.Unit();   // the way before
            delta = r_in.Mag() - Surface( r_in );
        }
    }

    // we found r_in properly
    return 1;

} //method WhereDoesItEnter_new








Position IceModel::WhereDoesItLeave(const Position &posnu,const Vector &nnu) const {
    // now get neutrino entry point...
    double p = posnu.Mag(); // radius of interaction
    double costheta = (nnu*posnu) / p; // theta of neutrino at interaction position
    double sintheta = sqrt(1-costheta*costheta);
    
    double lon = posnu.Lon();
    double lat = posnu.Lat();
    
    double a=0; // length of chord
    
    double R = Surface(lon,lat);
    double delta = R - p; // depth of the interaction
    // if interaction occurs below surface, as it should
    
    if (delta>-0.001) {
	a=sqrt(R*R*costheta*costheta+2*delta*R*sintheta*sintheta) - p*costheta; // chord length
	if (a<0) {
	    cout << "Negative chord length: " << a << "\n";
	} //end if
    } //end if (interaction below surface)  
    else if (delta<=-0.001) {
	
	//cout << "Error in interaction position.  whichray is " << whichray << "\n";
	cout << "lon, lat from WhereDoesItLeave is " << " " << lon << " " << lat << "\n";
	cout << "geoid, surface, p, surface-p are " << Geoid(lat) << " " << Surface(lon,lat) << " " << p << " , "<<(Surface(lon,lat)-p)<<"\n";
	
    } //else if: error: interaction takes place above the surface
    
    Position r_in = posnu + a*nnu;
    
    lon = r_in.Lon();
    lat = r_in.Lat();
    
    r_in = Surface(lon,lat) * r_in.Unit();
    
    return r_in;
} //method WhereDoesItLeave





// Below WhereDoesItEnterIce is from icemodel in icemc.
//--------------------------------------------------
// int IceModel::WhereDoesItEnterIce(const Position &posnu,
// 				       const Vector &nnu,
// 				       double stepsize,
// 				       Position &r_enterice) {
//   // now get exit point...
//   //   see my geometry notes.
//   // parameterize the neutrino trajectory and just see where it
//   // crosses the earth radius.
// 
//   //  Position r_enterice;
//   double distance=0;
//   int left_edge=0;
//   Position x = posnu;
//   double x2;
//   
//   Position x_previous = posnu;
// 
//   double x_previous2= x_previous * x_previous;
//   x2=x_previous2;
//   
//   double lon = x.Lon(),lat = x.Lat();
//   double lon_old = lon,lat_old = lat;
//   double local_surface = Surface(lon,lat);
//   double rock_previous2= pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
//   double surface_previous2=pow(local_surface,2);
// 
//   double rock2=rock_previous2;
//   double surface2=surface_previous2;
//   int foundit=0;  // keeps track of whether you found an ice entrance point
// 
//   //  cout << "lon, lat are " << posnu.Lon() << " " << posnu.Lat() << "\n";
//   //cout << "x2 at start is " << x2 << "\n";
//   while (distance<2*local_surface+1000) {
// 
//     distance+=stepsize;
// 
//     x -= stepsize*nnu;
//     x2=x*x;
//     //cout << "x2 is " << x2 << "\n";
//     lon = x.Lon();
//     lat = x.Lat();
// 
//       double ice_thickness=IceThickness(lon,lat);
//     if (lon!=lon_old || lat!=lat_old) {
//       local_surface = Surface(lon,lat);
// 
//       //if (lat>COASTLINE) 
//       //left_edge=1;
// 
//       rock2=pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
//       surface2=pow(local_surface,2);    
// 
//       if (ice_model==0) {
// 	if ((int)(lat)==COASTLINE && rock_previous2 < x2 && surface2 > x2)
// 	  left_edge=1;
//       } //if (Crust 2.0)
//     } //if (neutrino has stepped into new lon/lat bin)
// 
//     if ((((x_previous2>rock_previous2 && x2<rock2) // crosses rock boundary from above
// 	 || (x_previous2<surface_previous2 && x2>surface2)) && ice_thickness>0 && lat<COASTLINE) // crosses surface boundary from below
// 	|| left_edge) {
//       //  cout << "lat, COASTLINE, left_edge is " << lat << " " << COASTLINE<< " " << left_edge << "\n";
//       //cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
//       r_enterice = x;
//       // this gets you out of the loop.
//       //continue;
//       distance=3*Geoid(lat);
//       foundit=1;
//       //cout << "foundit is " << foundit << "\n";
//       //cout << "r_enterice is ";r_enterice.Print();
//       //continue;
//     } //if
// 
//     x_previous = x;
//     x_previous2 = x2;
//     //cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";
// 
//     if (lon!=lon_old || lat!=lat_old) {
//       rock_previous2 = rock2;
//       surface_previous2 = surface2;
//       lat_old = lat;
//       lon_old = lon;
//     } //if
// 
//   } //while
// 
//   return foundit;
// }//WhereDoesItEnterIce
// 
// 
// 
// 
// int IceModel::WhereDoesItExitIce(int inu,const Position &posnu,
// 				       const Vector &nnu,
// 				       double stepsize,
// 				       Position &r_enterice) {
//   // now get exit point...
//   //   see my geometry notes.
//   // parameterize the neutrino trajectory and just see where it
//   // crosses the earth radius.
// 
//   //  Position r_enterice;
//   double distance=0;
//   int left_edge=0;
//   Position x = posnu;
//   double x2;
//   
//   if (inu==1491) {
//    
//     cout << "posnu is";posnu.Print();
//     cout << "nnu is ";nnu.Print();
//   }
//    
// 
//   Position x_previous = posnu;
// 
//   double x_previous2= x_previous * x_previous;
//   x2=x_previous2;
//   
//   double lon = x.Lon(),lat = x.Lat();
//   double lon_old = lon,lat_old = lat;
//   double local_surface = Surface(lon,lat);
//   double rock_previous2= pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
//   double surface_previous2=pow(local_surface,2);
// 
//   double rock2=rock_previous2;
//   double surface2=surface_previous2;
//   int foundit=0;  // keeps track of whether you found an ice entrance point
// 
//  
// 
//   //  cout << "lon, lat are " << posnu.Lon() << " " << posnu.Lat() << "\n";
//   //cout << "x2 at start is " << x2 << "\n";
//   int nsteps=0;
//   while (distance<2*local_surface+1000) {
//     //cout << "another step.\n";
//     distance+=stepsize;
//     nsteps++;
//     //    cout << "inu, nsteps is " << inu << " " << nsteps << "\n";
//     x -= stepsize*nnu;
//     x2=x*x;
//     //cout << "x2 is " << x2 << "\n";
//     lon = x.Lon();
//     lat = x.Lat();
// 
//       double ice_thickness=IceThickness(lon,lat);
//     if (lon!=lon_old || lat!=lat_old) {
//       local_surface = Surface(lon,lat);
// 
//       //if (lat>COASTLINE) 
//       //left_edge=1;
// 
//       rock2=pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
//       surface2=pow(local_surface,2);    
// 
//       if (ice_model==0) {
// 	if ((int)(lat)==COASTLINE && rock_previous2 < x2 && surface2 > x2)
// 	  left_edge=1;
//       } //if (Crust 2.0)
//     } //if (neutrino has stepped into new lon/lat bin)
// 
//     if (inu==1491 && nsteps<10)
//       cout << "inu, x_previous2, rock_previous2, x2, rock2 are " << inu << " " << x_previous2 << " " << rock_previous2 << " " << x2 << " " << rock2 << "\n";
// 
//     if ((((x_previous2<rock_previous2 && x2>rock2) // crosses rock boundary from above
// 	 || (x_previous2>surface_previous2 && x2<surface2)) && ice_thickness>0 && lat<COASTLINE) // crosses surface boundary from above
// 	|| left_edge) {
//       //  cout << "lat, COASTLINE, left_edge is " << lat << " " << COASTLINE<< " " << left_edge << "\n";
//       //cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
//       r_enterice = x;
//       // this gets you out of the loop.
//       //continue;
//       distance=3*Geoid(lat);
//       foundit=1;
//       //cout << "foundit is " << foundit << "\n";
//       //continue;
//     } //if
// 
//     x_previous = x;
//     x_previous2 = x2;
//     //cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";
// 
//     if (lon!=lon_old || lat!=lat_old) {
//       rock_previous2 = rock2;
//       surface_previous2 = surface2;
//       lat_old = lat;
//       lon_old = lon;
//     } //if
// 
//   } //while
//   if (inu==0) {
//     cout << "r_enterice is ";r_enterice.Print();}
//   return foundit;
// }//WhereDoesItExitIce
//-------------------------------------------------- 





double IceModel::IceThickness(double lon, double lat) const {
  //This method returns the thickness of the ice in meters at a location specified by a latitude and longitude (in degrees).  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  double ice_thickness=0;
  //cout << "ice_model is " << ice_model << "\n";
  //cout << "icethkarray is " << icethkarray[(int)(lon/2)][(int)(lat/2)]*1000. << "\n";
  if (ice_model==1) {
    int e_coord=0;
    int n_coord=0;
    IceLonLattoEN(lon,lat,e_coord,n_coord);
    if (e_coord <= 1200 && e_coord >= 0 && n_coord <= 1000 && n_coord > 0)
      ice_thickness = ice_thickness_array[e_coord][n_coord]; //if this region has BEDMAP data, use it.
    else
      ice_thickness = icethkarray[(int)(lon/2)][(int)(lat/2)]*1000.; //if the location given is not covered by BEDMAP, use Crust 2.0 data
  } //BEDMAP ice thickness
  else if (ice_model==0) {
    ice_thickness = icethkarray[(int)(lon/2)][(int)(lat/2)]*1000.;
    //cout << "ilon, ilat are " << (int)(lon/2) << " " << (int)(lat/2) << "\n";
  } //Crust 2.0 ice thickness

  return ice_thickness;
} //method IceThickness
double IceModel::IceThickness(const Position &pos) const {
  //This method returns the thickness of the ice in meters at a location under a given position vector.  Code by Stephen Hoover.

  return IceThickness(pos.Lon(),pos.Lat());
} //method IceThickness(position)

double IceModel::SurfaceAboveGeoid(double lon, double lat) const {
  //This method returns the elevation above the geoid of the surface of the ice (or bare ground, if no ice is present) in meters, at a location specified by a latitude and longitude (in degrees).  In areas covered by water where no ice present, the method returns 0.  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  // lon must be 0 to 360
  double surface=0;

  if (ice_model==1) {
    int e_coord_ice=0;
    int n_coord_ice=0;
    int e_coord_ground=0;
    int n_coord_ground=0;
    IceLonLattoEN(lon,lat,e_coord_ice,n_coord_ice);
    GroundLonLattoEN(lon,lat,e_coord_ground,n_coord_ground);
    if (e_coord_ground <= 1068 && e_coord_ground >= 0 && n_coord_ground <= 869 && n_coord_ground >= 0 && e_coord_ice <= 1200 && e_coord_ice >= 0 && n_coord_ice <= 1000 && n_coord_ice >= 0)
      surface = ground_elevation[e_coord_ground][n_coord_ground] + ice_thickness_array[e_coord_ice][n_coord_ice] + water_depth[e_coord_ice][n_coord_ice];
    else
      surface = surfacer[(int)(lon/2)][(int)(lat/2)]; //If the position requested is outside the bounds of the BEDMAP data, use the Crust 2.0 data, regardless of the ice_model flag.
  } //Elevation of surface above geoid according to BEDMAP
  else if (ice_model==0) {
    surface = surfacer[(int)(lon/2)][(int)(lat/2)];
  } //Elevation of surface above geoid according to Crust 2.0

  return surface;
} //method SurfaceAboveGeoid

double IceModel::SurfaceAboveGeoid(const Position &pos) const {
  //This method returns the elevation above the geoid of the surface of the ice (or bare ground, if no ice is present) in meters, at a location specified by a position vector.  Code by Stephen Hoover.

  return SurfaceAboveGeoid(pos.Lon(),pos.Lat());
} //method SurfaceAboveGeoid(position)

double IceModel::Surface(double lon,double lat) const {
  return (SurfaceAboveGeoid(lon,lat) + Geoid(lat)); // distance from center of the earth to surface
} //Surface

double IceModel::Surface(const Position& pos) const {
  return Surface(pos.Lon(),pos.Lat());
} //Surface

double IceModel::WaterDepth(double lon, double lat) const {
  //This method returns the depth of water beneath ice shelves in meters, at a location specified by a latitude and longitude (in degrees).  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
  double water_depth_value=0;

  if (ice_model==0) {
      water_depth_value = waterthkarray[(int)(lon/2)][(int)(lat/2)]*1000;
  } //if(Crust 2.0)
  else if (ice_model==1) {
    int e_coord=0;
    int n_coord=0;
    WaterLonLattoEN(lon,lat,e_coord,n_coord);
    if (e_coord <= 1200 && e_coord >= 0 && n_coord <= 1000 && n_coord >= 0)
      water_depth_value = water_depth[e_coord][n_coord];
    else
      water_depth_value = waterthkarray[(int)(lon/2)][(int)(lat/2)]*1000;
  } //else if(BEDMAP)

  return water_depth_value;
} //method WaterDepth(longitude, latitude)
double IceModel::WaterDepth(const Position &pos) const {
  //This method returns the depth of water beneath ice shelves in meters, at a location specified by a position vector.  Code by Stephen Hoover.

  return WaterDepth(pos.Lon(),pos.Lat());
} //method WaterDepth(position)

int IceModel::IceOnWater(const Position &pos) const
{
  if(IceThickness(pos)>0.&&WaterDepth(pos)>0.)
    return 1;
  else return 0;

}
int IceModel::RossIceShelf(const Position &pos) const {
  int ilon,ilat;

  GetILonILat(pos,ilon,ilat);

  if ((ilat==2 && ilon>=5 && ilon<=14) ||
      (ilat==3 && (ilon>=168 || ilon<=14)) ||
      (ilat==4 && (ilon>=168 || ilon<=13)) ||
      (ilat==5 && (ilon>=168 || ilon<=14)))
    return 1;
  else
    return 0;
}//RossIceShelf

int IceModel::RossExcept(const Position &pos) const{
  int ilon,ilat;
  GetILonILat(pos,ilon,ilat);
if(ilon<=178&&ilon>=174&&ilat>=4&&ilat<=5)
    return 1;
  else 
    return 0;
}


int IceModel::RonneIceShelf(const Position &pos) const {
  int ilon,ilat;

  GetILonILat(pos,ilon,ilat);

  if ((ilat==4 && ilon>=52 && ilon<=74) ||
      (ilat==5 && ilon>=50 && ilon<=71) ||
      (ilat==6 && ilon>=55 && ilon<=64))
    return 1;
  else
    return 0;

}//RonneIceShelf

int IceModel::WestLand(const Position &pos) const {
  double lon = pos.Lon() , lat = pos.Lat();

  if((lat>=4&&lat<=26)&&((lon>=0&&lon<=180)||lon>=336))
    return 1;
  else return 0;

}//WestLand


int IceModel::OutsideAntarctica(const Position &pos) const {
  return (pos.Lat() >= COASTLINE);
} //OutsideAntarctica(Position)

int IceModel::OutsideAntarctica(double lat) const {
  return (lat >= COASTLINE);
} //OutsideAntarctica(double lat)

int IceModel::AcceptableRfexit(const Vector &nsurf_rfexit,const Position &rfexit,const Vector &n_exit2rx) const {

  //Make sure there's actually ice where the ray leaves
  if (rfexit.Lat()>COASTLINE || IceThickness(rfexit)<0.0001) {
    cout << "latitude is " << rfexit.Lat() << " compared to COASTLINE at " << COASTLINE << "\n";
    cout << "ice thickness is " << IceThickness(rfexit) << "\n";
    return 0;

  } //if

  if (nsurf_rfexit*n_exit2rx<0) {
    cout << "dot product is " << nsurf_rfexit*n_exit2rx << "\n";
    return 0;
  } //if

  return 1;
} //AcceptableRfexit

double IceModel::GetN(double altitude) const {
  // these are Peter's fit parameters
  double a1=0.463251;
  double b1=0.0140157;
  double n=0;

  if (altitude < FIRNDEPTH) 
    n=NICE;
  else if (altitude >= FIRNDEPTH && altitude <=0 && DEPTH_DEPENDENT_N) 
    //    N_DEPTH=NFIRN-(4.6198+13.62*(altitude_int/1000.))*
    //(altitude_int/1000.);   // Besson's equation for n(z)
    n=NFIRN+a1*(1.0-exp(b1*altitude));   // Peter's equation for n(z)
  else if (altitude > 0)
    cout<<"Error!  N requested for position in air!\n";
  else if (!DEPTH_DEPENDENT_N)
    n = NFIRN;

  return n;
} //GetN(altitude)

double IceModel::GetN(const Position &pos) const{
  return GetN(pos.Mag() - Surface(pos.Lon(),pos.Lat()));
} //GetN(Position)

double IceModel::EffectiveAttenuationLength(const Position &pos,const int &whichray) const {
  double localmaxdepth = IceThickness(pos);
  double depth = Surface(pos) - pos.Mag();
  
  int depth_index=0;
  double attenuation_length=0.0;
//   if (inu<10) {
//     cout << "pos is ";pos.Print();
//     cout << "surface is " << Surface(pos) << "\n";
//   }
  if(WestLand(pos) && !CONSTANTICETHICKNESS) 
    {
      depth_index=int(depth*419.9/localmaxdepth);//use 420 m ice shelf attenuation length data as the standard, squeeze or stretch if localmaxdepth is longer or shorter than 420m.
      if(RossIceShelf(pos) || RonneIceShelf(pos)) 
	{	  
	  if(whichray==0)
	    attenuation_length=l_shelfup[depth_index];
	  else if(whichray==1)
	    attenuation_length=l_shelfdown[depth_index];
	  else
	    cerr << " wrong attenuation length " <<endl;
	  
	  //for sanity check
	  if((depth_index+0.5)!=d_shelfup[depth_index])
	    {
	      cerr << "the index of the array l_iceshelfup is wrong!" << endl;
	      exit(1);
	    }
	}
      else //in ice sheet of westland
	{
	  if(whichray==0)
	    attenuation_length=l_westlandup[depth_index]; 
	  else if(whichray==1)
	    attenuation_length=l_westlanddown[depth_index];
	  else
	    cerr << " wrong attenuation length " <<endl;
      	}
       
      if(mooreBayFlag)//if use Moore's Bay measured data for the west land
	attenuation_length*=1.717557; //about 450 m (field attenuation length) for one whole way when assuming -3dB for the power loss at the bottom
    }
  else //in east antarctica or constant ice thickness
     { 
//        if (inu<10) {
//        cout << "localmaxdepth is " << localmaxdepth << "\n";
//        cout << "depth is " << depth << "\n";
//        }
       depth_index =int(depth*(2809.9/localmaxdepth));
       //if (inu<10)
	 //       cout << "depth_index is " << depth_index << "\n";

       if(whichray==0)
	 attenuation_length =l_sheetup[depth_index];
       else if(whichray==1)
	 attenuation_length =l_sheetdown[depth_index];
       else
	 cerr << " wrong attenuation length " <<endl;
     } //else

  return attenuation_length;
} //EffectiveAttenuationLengthUp


double IceModel::EffectiveAttenuationLength(Settings *settings1, const Position &pos,const int &whichray) const {
  double localmaxdepth = IceThickness(pos);
  double depth = Surface(pos) - pos.Mag();
  
  int depth_index=0;
  double attenuation_length=0.0;
//   if (inu<10) {
//     cout << "pos is ";pos.Print();
//     cout << "surface is " << Surface(pos) << "\n";
//   }
  if(WestLand(pos) && !CONSTANTICETHICKNESS) 
    {
      depth_index=int(depth*419.9/localmaxdepth);//use 420 m ice shelf attenuation length data as the standard, squeeze or stretch if localmaxdepth is longer or shorter than 420m.
      if(RossIceShelf(pos) || RonneIceShelf(pos)) 
	{	  
	  if(whichray==0)
	    attenuation_length=l_shelfup[depth_index];
	  else if(whichray==1)
	    attenuation_length=l_shelfdown[depth_index];
	  else
	    cerr << " wrong attenuation length " <<endl;
	  
	  //for sanity check
	  if((depth_index+0.5)!=d_shelfup[depth_index])
	    {
	      cerr << "the index of the array l_iceshelfup is wrong!" << endl;
	      exit(1);
	    }
	}
      else //in ice sheet of westland
	{
	  if(whichray==0)
	    attenuation_length=l_westlandup[depth_index]; 
	  else if(whichray==1)
	    attenuation_length=l_westlanddown[depth_index];
	  else
	    cerr << " wrong attenuation length " <<endl;
      	}
       
      //if(mooreBayFlag)//if use Moore's Bay measured data for the west land
      if(settings1->MOOREBAY)//if use Moore's Bay measured data for the west land
	attenuation_length*=1.717557; //about 450 m (field attenuation length) for one whole way when assuming -3dB for the power loss at the bottom
    }
  else //in east antarctica or constant ice thickness
     { 
//        if (inu<10) {
//        cout << "localmaxdepth is " << localmaxdepth << "\n";
//        cout << "depth is " << depth << "\n";
//        }
       depth_index =int(depth*(2809.9/localmaxdepth));
       //if (inu<10)
	 //       cout << "depth_index is " << depth_index << "\n";

       if(whichray==0)
	 attenuation_length =l_sheetup[depth_index];
       else if(whichray==1)
	 attenuation_length =l_sheetdown[depth_index];
       else
	 cerr << " wrong attenuation length " <<endl;
     } //else

  return attenuation_length;
} //EffectiveAttenuationLengthUp





double IceModel::Area(double latitude) const {
  //Returns the area of one square of the BEDMAP data at a given latitude. 
  double lat_rad = (90 - latitude) * RADDEG;

  return (pow(cellSize* ((1 + sin(71*RADDEG)) / (1 + sin(lat_rad))),2));
} //method Area

void IceModel::LonLattoEN(double lon, double lat, double xLowerLeft, double yLowerLeft, int& e_coord, int& n_coord) const {
  //takes as input a latitude and longitude (in degrees) and converts to indicies for BEDMAP matricies. Needs a location for the corner of the matrix, as not all the BEDMAP files cover the same area.  Code by Stephen Hoover.

  double easting=0;
  double northing=0;

  double lon_rad = (lon - 180) * RADDEG; //convert to radians, and shift origin to conventional spot
  double lat_rad = (90 - lat) * RADDEG;

  bedmap_R = scale_factor*bedmap_c_0 * pow(( (1 + eccentricity*sin(lat_rad)) / (1 - eccentricity*sin(lat_rad)) ),eccentricity/2) * tan((PI/4) - lat_rad/2);

  easting = bedmap_R * sin(lon_rad);
  northing = bedmap_R * cos(lon_rad);

  //  cout << "bedmap_R is " << bedmap_R << "\n";
  //cout << "easting, northing are " << easting << " " << northing << "\n";

  e_coord = (int)((easting - xLowerLeft) / cellSize);
  n_coord = (int)((-1*northing - yLowerLeft) / cellSize);

  return;
} //method LonLattoEN

void IceModel::IceLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
  //Converts a latitude and longitude (in degrees) to indicies for BEDMAP ice thickness data.  Code by Stephen Hoover.
  LonLattoEN(lon, lat, xLowerLeft_ice, yLowerLeft_ice, e_coord, n_coord);
}//IceLonLattoEN
void IceModel::GroundLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
  //Converts a latitude and longitude (in degrees) to indicies for BEDMAP ground elevation data.  Code by Stephen Hoover.
  LonLattoEN(lon, lat, xLowerLeft_ground, yLowerLeft_ground, e_coord, n_coord);
}//GroundLonLattoEN
void IceModel::WaterLonLattoEN(double lon, double lat, int& e_coord, int& n_coord) const {
  //Converts a latitude and longitude (in degrees) to indicies for BEDMAP water depth data.  Code by Stephen Hoover.
  LonLattoEN(lon, lat, xLowerLeft_water, yLowerLeft_water, e_coord, n_coord);
}//WaterLonLattoEN

void IceModel::ENtoLonLat(int e_coord, int n_coord, double xLowerLeft, double yLowerLeft, double& lon, double& lat) const {
  //Takes as input the indicies from a BEDMAP data set, and turns them into latitude and longitude coordinates.  Information on which data set (surface data, ice depth, water depth) is necessary, in the form of coordinates of a corner of the map.  Code by Stephen Hoover.

  double isometric_lat=0;
  double easting = xLowerLeft+(cellSize*(e_coord+0.5)); //Add offset of 0.5 to get coordinates of middle of cell instead of edges.
  double northing = -1*(yLowerLeft+(cellSize*(n_coord+0.5)));

  //  cout << "easting, northing are " << easting << " " << northing << "\n";

  //first set longitude

  if (northing!=0)
    lon = atan(easting/northing);
  else
    lon = 90*RADDEG;

  // this puts lon between -pi and pi
  if (easting > 0 && lon < 0) //adjust sign of longitude
    lon += PI;
  else if (easting < 0 && lon > 0)
    lon -= PI;
  else if (easting == 0 && northing < 0)
    lon += PI;

  //  now find latitude

  if (easting != 0)
    bedmap_R = fabs(easting/sin(lon));
  else if (easting == 0 && northing != 0)
    bedmap_R = fabs(northing);
  else {
    lat = 0; //at the pole, set lat=0 degrees
    lon = lon*DEGRAD; // now put lon between 180 and 180 (only at pol)
    return;
  } //else

  isometric_lat = (PI/2) - 2*atan(bedmap_R/(scale_factor*bedmap_c_0));

  lat = isometric_lat + bedmap_a_bar*sin(2*isometric_lat) + bedmap_b_bar*sin(4*isometric_lat) + bedmap_c_bar*sin(6*isometric_lat) + bedmap_d_bar*sin(8*isometric_lat);

  lon = lon * DEGRAD + 180;  //convert to degrees, shift 0 to line up with bin 0 of Crust 2.0
  lat = 90 - lat*DEGRAD; //convert to degrees, with 0 degrees at the south pole

  //  if (lon>160 && lon<165)
  //cout << "e_coord, n_coord, easting, northing, lon are " << e_coord << " " << n_coord << " " << easting << " " << northing << " " << lon << "\n";
  return;
  
} //method ENtoLonLat

void IceModel::IceENtoLonLat(int e, int n, double& lon, double& lat) const {
  //Converts indicies of the BEDMAP ice thickness matrix into longitude and latitude.  Code by Stephen Hoover.
  // cout << "I'm inside IceENtoLonLat.\n";
  ENtoLonLat(e,n,xLowerLeft_ice,yLowerLeft_ice,lon,lat);
}//IceENtoLonLat
void IceModel::GroundENtoLonLat(int e, int n, double& lon, double& lat) const {
  //Converts indicies of the BEDMAP ground elevation matrix into longitude and latitude.  Code by Stephen Hoover.
  ENtoLonLat(e,n,xLowerLeft_ground,yLowerLeft_ground,lon,lat);
}//GroundENtoLonLat
void IceModel::WaterENtoLonLat(int e, int n, double& lon, double& lat) const {
  //Converts indicies of the BEDMAP water depth matrix into longitude and latitude.  Code by Stephen Hoover.
  ENtoLonLat(e,n,xLowerLeft_water,yLowerLeft_water,lon,lat);
}//WaterENtoLonLat



void IceModel::ReadIceThickness() {
  //Reads the BEDMAP ice thickness data.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover

  ifstream IceThicknessFile("data/icethic.asc");
  if(!IceThicknessFile) {
    cerr << "Couldn't open: data/icethic.asc" << endl;
    exit(1);
  }

  cout<<"Reading in BEDMAP data on ice thickness.\n";

  string tempBuf1;
  string tempBuf2;
  string tempBuf3;
  string tempBuf4;
  string tempBuf5;
  string tempBuf6;
  int temp1,temp2,temp3,temp4,temp5,temp6;
  
  IceThicknessFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
		   >> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
		   >> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
  
  if(tempBuf1 == string("ncols")) {
    nCols_ice=temp1;
  }
  if(tempBuf2 == string("nrows")) {
    nRows_ice=temp2;
  }
  if(tempBuf3 == string("xllcorner")) {
    xLowerLeft_ice=temp3;
  }
  if(tempBuf4 == string("yllcorner")) {
    yLowerLeft_ice=temp4;
  }
  if(tempBuf5 == string("cellsize")) {
    cellSize=temp5;
  }
  if(tempBuf6 == string("NODATA_value")) {
    NODATA=temp6;
  }
  //cout<<"nCols_ice, nRows_ice "<<nCols_ice<<" , "<<nRows_ice<<endl;
  //cout<<"xLL_ice, yLL_ice, cellsize "<<xLowerLeft_ice<<" , "<<yLowerLeft_ice<<" , "<<cellSize<<endl<<endl;
  
  double theValue;
  for(int rowNum=0;rowNum<nRows_ice;rowNum++) {
    for(int colNum=0;colNum<nCols_ice;colNum++) {
      IceThicknessFile >> theValue;
      if(theValue==NODATA)
	theValue=0; //Set ice depth to 0 where we have no data.
      ice_thickness_array[colNum][rowNum] = double(theValue); //This stores data as ice_thickness_array[easting][northing]
    }//for
  }//for
  
  IceThicknessFile.close();
  return;
} //method ReadIceThickness

void IceModel::ReadGroundBed() {
  //Reads the BEDMAP data on the elevation of the ground beneath the ice.  If there is water beneath the ice, the ground elevation is given the value 0.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
  ifstream GroundBedFile("data/groundbed.asc");
  if(!GroundBedFile) {
    cerr << "Couldn't open: data/groundbed.asc" << endl;
    exit(1);
  }

  cout<<"Reading in BEDMAP data on elevation of ground.\n";

  string tempBuf1;
  string tempBuf2;
  string tempBuf3;
  string tempBuf4;
  string tempBuf5;
  string tempBuf6;
  int temp1,temp2,temp3,temp4,temp5,temp6;
  
  GroundBedFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
		>> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
		>> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
  
  if(tempBuf1 == string("ncols")) {
    nCols_ground=temp1;
  }
  if(tempBuf2 == string("nrows")) {
    nRows_ground=temp2;
  }
  if(tempBuf3 == string("xllcorner")) {
    xLowerLeft_ground=temp3;
  }
  if(tempBuf4 == string("yllcorner")) {
    yLowerLeft_ground=temp4;
  }
  if(tempBuf5 == string("cellsize")) {
    cellSize=temp5;
  }
  if(tempBuf6 == string("NODATA_value")) {
    NODATA=temp6;
  }

  //cout<<"nCols_ground, nRows_ground "<<nCols_ground<<" , "<<nRows_ground<<endl;
  //cout<<"xLL_ground, yLL_ground, cellsize "<<xLowerLeft_ground<<" , "<<yLowerLeft_ground<<" , "<<cellSize<<endl<<endl;
  
  double theValue;
  for(int rowNum=0;rowNum<nRows_ground;rowNum++) {
    for(int colNum=0;colNum<nCols_ground;colNum++) {
      GroundBedFile >> theValue;
      
      if(theValue==NODATA)
	theValue=0; //Set elevation to 0 where we have no data.
      ground_elevation[colNum][rowNum] = double(theValue);
      //if (theValue != -96 && theValue != 0)
      //cout<<"ground_elevation: "<<theValue<<endl;
    }//for
  }//for
  
  GroundBedFile.close();
  return;
} //method ReadGroundBed

void IceModel::ReadWaterDepth() {
  //Reads BEDMAP data on the depth of water beneath the ice.  Where no water is present, the value 0 is entered.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
  ifstream WaterDepthFile("data/water.asc");
  if(!WaterDepthFile) {
    cerr << "Couldn't open: data/water.asc" << endl;
    exit(1);
  }

  cout<<"Reading in BEDMAP data on water depth.\n";

  string tempBuf1;
  string tempBuf2;
  string tempBuf3;
  string tempBuf4;
  string tempBuf5;
  string tempBuf6;
  int temp1,temp2,temp3,temp4,temp5,temp6;
  
  WaterDepthFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
		 >> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
		 >> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
  
  if(tempBuf1 == string("ncols")) {
    nCols_water=temp1;
  }
  if(tempBuf2 == string("nrows")) {
    nRows_water=temp2;
  }
  if(tempBuf3 == string("xllcorner")) {
    xLowerLeft_water=temp3;
  }
  if(tempBuf4 == string("yllcorner")) {
    yLowerLeft_water=temp4;
  }
  if(tempBuf5 == string("cellsize")) {
    cellSize=temp5;
  }
  if(tempBuf6 == string("NODATA_value")) {
    NODATA=temp6;
  }

  //cout<<"nCols_water, nRows_water "<<nCols_water<<" , "<<nRows_water<<endl;
  //cout<<"xLL_water, yLL_water, cellsize "<<xLowerLeft_water<<" , "<<yLowerLeft_water<<" , "<<cellSize<<endl<<endl;
  
  double theValue;
  for(int rowNum=0;rowNum<nRows_water;rowNum++) {
    for(int colNum=0;colNum<nCols_water;colNum++) {
      WaterDepthFile >> theValue;
      
      if(theValue==NODATA)
	theValue=0; //Set depth to 0 where we have no data.
      water_depth[colNum][rowNum] = double(theValue);
    }//for
  }//for
  
  WaterDepthFile.close();
  return;
} //method ReadWaterDepth



void IceModel::GetFresnel_slappy (double i_ang, double n1, double n2, double &r, double &t ) {

    //calculate the fresnel coefficient at i_ang incident angle, n1 (incident medium index of refraction), n2.
    //r for reflection coefficient,
    //t for transmitted coefficient.
    //

    if (n1/n2 * sin(i_ang) >= 1.) { // total internal reflection
        r = 1.;
        t = 0.;
    }
    else {
        double t_ang = asin( n1/n2 * sin(i_ang) );

        r = sin (i_ang - t_ang) / sin (i_ang + t_ang);

        t = 2.*sin(t_ang) * cos(i_ang) / (sin(i_ang + t_ang));
    }

    return;
}


void IceModel::GetFresnel_pokey (double i_ang, double n1, double n2, double &r, double &t ) {

    //calculate the fresnel coefficient at i_ang incident angle, n1 (incident medium index of refraction), n2.
    //r for reflection coefficient,
    //t for transmitted coefficient.
    //

    if (n1/n2 * sin(i_ang) >= 1.) { // total internal reflection
        r = 1.;
        t = 0.;
    }
    else {
        double t_ang = asin( n1/n2 * sin(i_ang) );

        r = tan (i_ang - t_ang) / tan (i_ang + t_ang);

        t = 2.*sin(t_ang) * cos(i_ang) / (sin(i_ang + t_ang) * cos(i_ang - t_ang) );
    }

    return;
}


void IceModel::GetFresnel (
        double launch_angle, double rec_angle,
        double refl_angle, Position &posnu, Vector &launch_vector, Vector &rec_vector, Settings *settings1, double &fresnel, double &mag,
        Vector &Pol // will read the polarization at the source and return polarization at the target antenna
        ) {

    double n1 = 1.35;   // index of refraction at the firn
    double n2 = 1.;     // index of refraction at the air

    // calculate the magnification factor for plane / spherical wave case
    if (settings1->WAVE_TYPE == 0) { // plane wave
        mag = sqrt( abs(tan(rec_angle)) / abs(tan (launch_angle)) );
    }
    else if (settings1->WAVE_TYPE == 1) { // spherical wave
        mag = sqrt( abs(tan (launch_angle)) / abs(tan (rec_angle)) );
    }

    Vector perp = launch_vector.Cross( posnu ).Unit();    // perp unit vector it should be same in both src and trg
    Vector src_parallel = perp.Cross( launch_vector ).Unit();
    Vector trg_parallel = perp.Cross( rec_vector ).Unit();

    double pol_perp_src = Pol * perp;
    double pol_parallel_src = Pol * src_parallel;
    double pol_perp_trg=0, pol_parallel_trg=0;

    // check if ray is reflected or not
    if (refl_angle < PI/2.) {   // the ray is reflected at the surface

        double r_coeff_pokey, r_coeff_slappy;

        if (n1/n2 * sin(launch_angle) >= 1.) {  // total internal reflection case
            r_coeff_pokey = 1.;
            r_coeff_slappy = 1.;
        }

        else {  // there is refracted ray to air

            double t_angle = asin( n1/n2 * sin(launch_angle) ); // transmitted ray (which we don't care) angle

            r_coeff_pokey = tan(launch_angle - t_angle) / tan(launch_angle + t_angle);  // only reflected ray can be a signal
            r_coeff_slappy = sin(launch_angle - t_angle) / sin(launch_angle + t_angle);

        }

        pol_parallel_trg = r_coeff_pokey * pol_parallel_src;
        pol_perp_trg = r_coeff_slappy * pol_perp_src;

    }
    else {      // ray didn't relfected at the surface; no fresnel coeff need to be applied
        pol_parallel_trg = pol_parallel_src;
        pol_perp_trg = pol_perp_src;
    }

    fresnel = sqrt( pow(pol_perp_trg,2) + pow(pol_parallel_trg,2) );

    Pol = (pol_perp_trg * perp + pol_parallel_trg * trg_parallel).Unit();


}


//void IceModel::FillArraysforTree(double icethck[1200][1000],double elev[1068][869],double lon_ground[1068][869],double lat_ground[1068][869],double lon_ice[1200][1000],double lat_ice[1200][1000],double h20_depth[1200][1000],double lon_water[1200][1000],double lat_water[1200][1000]) {
// void IceModel::FillArraysforTree(double lon_ground[1068][869],double lat_ground[1068][869],double lon_ice[1200][1000],double lat_ice[1200][1000],double lon_water[1200][1000],double lat_water[1200][1000]) {
 
  //  for (int rowNum=0;rowNum<nRows_ice;rowNum++) {
  //for (int colNum=0;colNum<nCols_ice;colNum++) {
//   for (int i=0;i<nRows_ice;i++) {
//     for (int j=0;j<nCols_ice;j++) {
//       //     this->IceENtoLonLat(colNum,rowNum,lon_ice[colNum][rowNum],lat_ice[colNum][rowNum]); //Recall that the e / n coordinates in horizon were picked from the ground bed array.
//       double test1,test2;
//       cout << "rowNum, colNum are " << i << " " << j << "\n";
//       //      cout << "lon_ice is " << rowNum << " " << colNum << " " << lon_ice[rowNum][colNum] << "\n";
//       //this->IceENtoLonLat(i,j,test1,test2); //Recall that the e / n coordinates in horizon were picked from the ground bed array.
//       //icethck[colNum][rowNum]=IceThickness(lon_ice[colNum][rowNum],lat_ice[colNum][rowNum]);
     
//       //      WaterENtoLonLat(colNum,rowNum,lon_water[colNum][rowNum],lat_water[colNum][rowNum]);
//       //h20_depth[colNum][rowNum]=WaterDepth(lon_water[colNum][rowNum],lat_water[colNum][rowNum]);

//     }
//   }h
//   for (int n_coord=0;n_coord<nRows_ground;n_coord++) {
//     for (int e_coord=0;e_coord<nCols_ground;e_coord++) {
//       GroundENtoLonLat(e_coord,n_coord,lon_ground[e_coord][n_coord],lat_ground[e_coord][n_coord]);
//       //     elev[e_coord][n_coord] = SurfaceAboveGeoid(lon_ice[e_coord][n_coord],lat_ice[e_coord][n_coord]);

//     }
//   }





//}
