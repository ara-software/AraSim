#include "Vector.h"
#include "Constants.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRotation.h"
#include <cmath>
#include <iostream>

ClassImp(Vector);


using namespace std;

 Vector::Vector(double x_inp,double y_inp,double z_inp) {
 x = x_inp; 
  y = y_inp; 
  z = z_inp;
  UpdateThetaPhi();
} //Constructor Vector(double,double,double)
 Vector::Vector(double* xarray) {
  x=xarray[0];
  y=xarray[1];
  z=xarray[2];
  UpdateThetaPhi();
}
 Vector::Vector(double theta_inp, double phi_inp) {

  if (theta_inp < 0.0 || theta_inp > PI) {

    cout<<"Error!  Attempt to construct Vector from invalid theta!\n";

    x = 0.; 
    y = 0.; 
    z = 1.; 
    UpdateThetaPhi();
  } //check to see if theta is valid

  x = sin(theta_inp) * cos(phi_inp);
  y = sin(theta_inp) * sin(phi_inp);
  z = cos(theta_inp);

  UpdateThetaPhi();
} //Constructor Vector(theta,phi)


 Vector::Vector() {
  x = 0.; 
  y = 0.; 
  z = 1.; 
  UpdateThetaPhi();
} //Vector default constructor

 Vector Vector::Cross(const Vector &vec) const {
  return Vector(y * vec.z - z * vec.y,
		-x * vec.z + z * vec.x,
		x * vec.y - y * vec.x);
} //Vector::Cross

double Vector::Dot(const Vector &vec) const
{return x * vec.x +y * vec.y + z * vec.z ;}
//Takes the dot product  this x vec.

 double Vector::Mag() const {
  return sqrt(x*x + y*y + z*z);
} //Vector::Mag

 double Vector::Angle(const Vector &vec) const {
  return acos(((*this)*vec) / (this->Mag() * vec.Mag()));
} //Vector::Angle

 Vector Vector::ChangeCoord(const Vector &new_x_axis,const Vector &new_y_axis) const {
  

  TVector3 temp;
  temp.SetX(this->GetX());
  temp.SetY(this->GetY());
  temp.SetZ(this->GetZ());
  
  TVector3 tnew_x_axis;
  tnew_x_axis.SetX(new_x_axis.GetX());
  tnew_x_axis.SetY(new_x_axis.GetY());
  tnew_x_axis.SetZ(new_x_axis.GetZ());
  
  TVector3 tnew_y_axis;
  tnew_y_axis.SetX(new_y_axis.GetX());
  tnew_y_axis.SetY(new_y_axis.GetY());
  tnew_y_axis.SetZ(new_y_axis.GetZ());
  
  TRotation r;
  r.SetXAxis(tnew_x_axis);
  r.SetYAxis(tnew_y_axis);
  r.SetZAxis(tnew_x_axis.Cross(tnew_y_axis));
  
  
  //r.SetToIdentity();
  //r.RotateAxes(newX,newY,newZ);
  //r.Invert();
  temp.Transform(r);
  
  Vector new_vector;
  new_vector.SetX(temp.X());
  new_vector.SetY(temp.Y());
  new_vector.SetZ(temp.Z());
  
  
  
  //Vector new_vector = this->RotateY(new_z_axis.theta);
  //new_vector = new_vector.RotateZ(new_z_axis.phi);

  return new_vector;
} //Vector::ChangeCoord

 Vector Vector::ChangeCoord(const Vector &new_z_axis) const {

  Vector new_vector = this->RotateY(new_z_axis.theta);
  new_vector = new_vector.RotateZ(new_z_axis.phi);

  return new_vector;

}

 Vector Vector::Unit() const {
  return (*this) / this->Mag();
} //Vector::Unit

void Vector::Print() const {
  cout << x << " " << y << " " << z << "\n";
}


void Vector::PrintSp() const {
  cout << r << " " << theta << " " << phi << "\n";
}

 double Vector::GetX() const {
  return x;
} //Vector::GetX

 double Vector::GetY() const {
  return y;
} //Vector::GetY

 double Vector::GetZ() const {
  return z;
} //Vector::GetZ


 double Vector::R() const {
  return r;
} //Vector::R

 double Vector::Theta() const {
  return theta;
} //Vector::Theta

 double Vector::Phi() const {
  return phi;
} //Vector::Phi

 void Vector::SetX(double inp) {
  x = inp;
  UpdateThetaPhi();
} //Vector::SetX

 void Vector::SetY(double inp) {
  y = inp;
  UpdateThetaPhi();
} //Vector::SetY

 void Vector::SetZ(double inp) {
  z = inp;
  UpdateThetaPhi();
} //Vector::SetZ

 void Vector::SetXYZ(double inpx,double inpy,double inpz) {
  x = inpx;
  y = inpy;
  z = inpz;
  UpdateThetaPhi();
} //Vector::SetXYZ


 void Vector::SetR(double inpr) {
  r = inpr;
  UpdateXYZ();
} //Vector::SetRThetaPhi


 void Vector::SetThetaPhi(double inptheta,double inpphi) {
  theta = inptheta;
  phi = inpphi;
  UpdateXYZ();
} //Vector::SetRThetaPhi


 void Vector::SetRThetaPhi(double inpr,double inptheta,double inpphi) {
  r = inpr;
  theta = inptheta;
  phi = inpphi;
  UpdateXYZ();
} //Vector::SetRThetaPhi


 void Vector::Reset(double x_inp, double y_inp, double z_inp) {
  x = x_inp;
  y = y_inp;
  z = z_inp;
  UpdateThetaPhi();
} //Vector::Reset

 Vector Vector::RotateX(double angle) const {
  double new_x = x;
  double new_y = cos(angle)*y - sin(angle)*z;
  double new_z = sin(angle)*y + cos(angle)*z;
  Vector rotated_vector(new_x,new_y,new_z);
  return rotated_vector;
} //RotateX

 Vector Vector::RotateY(double angle) const {
  double new_x = cos(angle)*x + sin(angle)*z;
  double new_y = y;
  double new_z = -sin(angle)*x + cos(angle)*z;
  Vector rotated_vector(new_x,new_y,new_z);
  return rotated_vector;
} //RotateY

 Vector Vector::RotateZ(double angle) const {
  double new_x = cos(angle)*x - sin(angle)*y;
  double new_y = sin(angle)*x + cos(angle)*y;
  double new_z = z;
  Vector rotated_vector(new_x,new_y,new_z);
  return rotated_vector;
} //RotateZ

 Vector Vector::Rotate(const double angle,const Vector& axis) const {
  //Code blatently stolen from Root's TRotation::Rotate method
  //Example: If you rotate the vector (0,0,1) by 90 degrees around the vector (0,1,0), the result is (1,0,0).
  double length = axis.Mag();

  double s = sin(angle);
  double c = cos(angle);
  double dx = axis.x / length;
  double dy = axis.y / length;
  double dz = axis.z / length;

  double newx = (c+(1-c)*dx*dx) * x + ((1-c)*dx*dy-s*dz) * y + ((1-c)*dx*dz+s*dy) * z;
  double newy = ((1-c)*dy*dx+s*dz) * x + (c+(1-c)*dy*dy) * y + ((1-c)*dy*dz-s*dx) * z;
  double newz = ((1-c)*dz*dx-s*dy) * x + ((1-c)*dz*dy+s*dx) * y + (c+(1-c)*dz*dz) * z;

  return Vector(newx,newy,newz);
} //Vector::Rotate

 void Vector::UpdateThetaPhi() {
  //This is a private method that will calculate values of theta and phi from x,y,z coordinates,
  //and store the results in the class variables.
  //double transverse = hypot(x,y);
  double transverse = sqrt(x*x+y*y);

  // atan2 outputs in the range -pi to pi
  theta = atan2(transverse,z);
 
  phi=atan2(y,x);

  r = sqrt(x*x + y*y + z*z);

  if (phi<0)
    phi+=2*PI;
  // phi is now from 0 to 2*pi wrt +x
  
} //UpdateThetaPhi

void Vector::UpdateXYZ() {
    
  x = r * sin(theta) * cos(phi);
  y = r * sin(theta) * sin(phi);
  z = r * cos(theta);


}


 Vector Vector::Zero() {
  //Zero the vector
   
  x=0;
  y=0;
  z=0;
  return Vector(x,y,z);
} // Zero


 void Vector::RotateUz(const Vector &NewUzVector) {
   // Change original z-axis to match the new vector
   // Note the new vector must be normalized
   // stolen from TVector3

   double u1 = NewUzVector.GetX();
   double u2 = NewUzVector.GetY();
   double u3 = NewUzVector.GetZ();
   double up = u1*u1 + u2*u2;
   
   if(up){
     up = sqrt(up);
     double px = x;
     double py = y;
     double pz = z;
     x = (u1*u3*px - u2*py + u1*up*pz)/up;
     y = (u2*u3*px + u1*py + u2*up*pz)/up;
     z = (u3*u3*px -    px + u3*up*pz)/up;
   }else if(u3 < 0.){ // phi = 0 theta = pi
     x = -x;
     z = -z;
   }
   else{};
} // RotateUz



 double Vector::operator [](int i) const {
  //code taken from ROOT's TVector3 class
  switch(i) {
  case 0:
    return x;
  case 1:
    return y;
  case 2:
    return z;
  default:
    printf("Vector::operator[](i) has been given a bad index: %d.  Returning zero.",i);
  } //end switch

  return 0.;
} //operator[]

//Overloaded operators - friends of Vector, not member functions.

 Vector operator +(const Vector& vector1, const Vector& vector2) {
  return Vector(vector1.x + vector2.x , vector1.y + vector2.y , vector1.z + vector2.z);
} //Vector overloaded operator +

 void operator +=(Vector& vector1, const Vector& vector2) {
  vector1 = vector1 + vector2;
} //Vector overloaded operator +=

 Vector operator -(const Vector& vector1, const Vector& vector2) {
  return Vector(vector1.x - vector2.x , vector1.y - vector2.y , vector1.z - vector2.z);
} //Vector overloaded operator -

 void operator -=(Vector& vector1, const Vector& vector2) {
  vector1 = vector1 - vector2;
} //Vector overloaded operator -=

 Vector operator -(const Vector& vec) {
  return Vector(-vec.x,-vec.y,-vec.z);
} //Vector overloaded operator - (negative)

 double operator *(const Vector& vector1, const Vector& vector2) {
  return ((vector1.x * vector2.x) + (vector1.y * vector2.y) + (vector1.z * vector2.z));
} //Vector overloaded operator * (dot product)

 Vector operator /(const Vector &v, const double& a) {
  return Vector(v.x / a, v.y / a, v.z / a);
} //Vector overloaded operator /

 Vector operator *(const double& a, const Vector& v) {
  return Vector(a*v.x , a*v.y , a*v.z);
} //Vector overloaded operator * (scalar multiplication)

 ostream& operator <<(ostream& outs, const Vector& vec) {
  cout<<vec.x<<","<<vec.y<<","<<vec.z; 
  return outs;
} //Vector overloaded operator <<
