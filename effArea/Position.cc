#include "TRandom3.h"
#include "Position.h"
#include "EarthModel.h"
#include "Constants.h"
#include <cmath>

ClassImp(Position);

// Vector Constants
static Vector x_axis = Vector(1,0,0);
static Vector y_axis = Vector(0,1,0);
static Vector z_axis = Vector(0,0,1);

Position::Position() : Vector() 
{
  //This method intentionally left blank.
} //Position default constructor

Position::Position(Vector vec) : Vector(vec[0],vec[1],vec[2]) 
{
  //This method intentionally left blank.
} //Position constructor from Vector

Position::Position(double longitude, double latitude, double altitude) {
  Vector location = z_axis;
  theta = latitude * RADDEG;

  phi=EarthModel::LongtoPhi_0isPrimeMeridian(longitude); // convert longitude (-180 to 180) to phi (0 to 2pi wrt 90E, counter-clockwise)

  location = location.RotateY(theta);
  location = location.RotateZ(phi);
  location = altitude * location;

  x = location[0];
  y = location[1];
  z = location[2];
} //Position constructor from longitude and latitude

Position::Position(double theta_inp, double phi_inp) : Vector(theta_inp,phi_inp) 
{
  //This method intentionally left blank.
} //Constructor Position(theta,phi)

double Position::Distance(const Position &second) const {
  return sqrt((x - second.x)*(x-second.x) 
	      + (y - second.y)*(y-second.y) 
	      + (z - second.z)*(z-second.z)); // it saves time to multiply them like this rather than use pow
} //Position::Distance

double Position::SurfaceDistance(const Position &second, double local_surface) const {
  return  this->Angle(second) * local_surface;
} //Position::SurfaceDistance

double Position::Lat() const {
  return theta*DEGRAD;
} //Position::Lat

double Position::Lon() const {
  double phi_deg = phi*DEGRAD;

  if (phi_deg > 270.)
    phi_deg = phi_deg-360.;
  
  return (270. - phi_deg);
} //Position::Lon
