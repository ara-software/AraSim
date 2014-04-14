#ifndef POSITION_H
#define POSITION_H

#include "Vector.h"
//#include "TObject.h"

////////////////////////////////////////////////////////////////////////////////////////////////
//class Position:
//This class is a 3-vector that represents a position on the Earth's surface.
//
//   Methods:
//
// Lat()  : Returns latitude of this position.
//
// Lon()  : Returns longitude of this position.
//
// Distance(second position)  : Returns distance between this position and a second.  Takes
//                              a Position as input.
//
// SurfaceDistance(second position, surface elevation)  : Returns distance over the surface of
//              the Earth between the spot on the surface under/above this position and the spot
//              under/above a second position.
//              Input: a Position, and the distance from the center of the Earth to the surface at
//                     this Position.
//
////////////////////////////////////////////////////////////////////////////////////////////////

//class Position : public Vector, public TObject {
class Position : public Vector {
public:
  Position();
  //Default constructor: calls default constructor of Vector

  Position(Vector vec);

  Position(double theta_inp, double phi_inp);
  //Identical to the Vector constructor with the same inputs.

  Position(double longitude, double latitude, double altitude);
  //Constructs a position vector given a longitude, latitude,
  //and distance from the center of the Earth.

  double Lat() const;
  //Returns latitude, where the +z direction is at 0 latitude.

  double Lon() const;
  //Returns longitude, where 0 degrees longitude corresponds to phi = 270 deg (-y axis),
  //and increases clockwise.
  //Hence, lon=0 -> phi=270,
  //lon=90 -> phi=180,
  //lon=180 -> phi=90,
  //lon=270 -> phi=0,

  double Distance(const Position &second) const;
  //Returns chord distance (direct distance between two vectors)

  double SurfaceDistance(const Position &second, double local_surface) const;
  //Returns "surface distance" between two positions.  The surface distance
  //is the the length of arc between two positions.
  //Altitude (i.e. length of the position vector) is irrelevant; only angle
  //between the two position vectors is considered.

  ClassDef(Position,1);

}; //class Position


#endif //POSITION_H
