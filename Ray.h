////////////////////////////////////////////////////////////////////////////////////////////////
//class Ray:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef RAY_H
#define RAY_H

#include <cmath>
#include <iostream>

class Event;
class Position;
class Vector;
class IceModel;
class Detector;
class Interaction;

using std::vector;
using std::cout;


class Ray {

 private:
//--------------------------------------------------
//   int WhereDoesItEnterIce(const Position &posnu,
// 			       const Vector &nnu,
// 			       double stepsize,
// 			       Position &r_enterice,
//                                IceModel *antarctica);
// 
//   int WhereDoesItExitIce(const Position &posnu,
// 			 const Vector &nnu,
// 			 double stepsize,
// 			 Position &r_enterice,
//                          IceModel *antarctica);
//   void FlattoEarth (Interaction *interaction1, IceModel *antarctica, double X, double Y, double D);
// 
//-------------------------------------------------- 

 public:
  Ray();
  ~Ray();

//--------------------------------------------------
//   vector<double> getEField(Event *event,Position pos); // Get E field from the shower at a particular location pos
//-------------------------------------------------- 


  //
  // below WhereDoesItLeave member is copied from ray.hh in icemc.
//--------------------------------------------------
//   static int WhereDoesItLeave(const Position &posnu,
// 			      const Vector &ntemp,IceModel *antarctica,
// 			      Position &r_out); 
//   int PickUnbiased(Interaction *interaction1, IceModel *antarctica);
//-------------------------------------------------- 
//--------------------------------------------------
//   void PickNear(Interaction *interaction1, IceModel *antarctica, Detector *detector);
//-------------------------------------------------- 
    

};

#endif //RAY_H

