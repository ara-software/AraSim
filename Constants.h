#ifndef CONSTANTS_H
#define CONSTANTS_H


// constants in math
const double TWOPI=6.2831852;
const double PI=3.141592654;
const double ALOG2=0.693147;        // natural log of 2
const double INV_E=0.36787944;      // 1/e
const double sr=4*PI;


// conversion constants
const double CMINCH=2.54;          // inches to cm
const double RADDEG=0.017453292;   // radians/degree  
const double DEGRAD=57.2957795;    // degree/rad

// physical constants
const double CLIGHT=3.0E8;            // speed of light m/s
const double KBOLTZ=1.38E-23;         // Boltzmann constant J/K
const double Z0=377.;                 // resistivity of free space
const double Zr=50.;  // radiation resistance (50 Ohms?)
const double M_NUCL=1.66E-27;         // amu mass in kg
const double MTAU=1.777E9;            // mass of the tau
const double TAUDECAY_TIME=290.6E-15; // lifetime of tau

// properties of water
const double X0H20=0.361;          // radiation length of water (meters)


// properties of air


const double Z_AIR=377;            // resistance of air = sqrt(epsilon/mu)
  const double RHOAIR=1.25;          // density of air (kg/m**3)
// // properties of ice


const double NFIRN=1.3250;                   // index of refraction at the very surface - Peter
const double NICE=1.79;                      // index of refraction of ice

enum {kNC, kCC};    // neutrino interaction constant nc : 0, cc :1


const double poissonerror_minus[21] = {0.-0.00,1.-0.37,2.-0.74,3.-1.10,4.-2.34,5.-2.75,6.-3.82,7.-4.25,8.-5.30,9.-6.33,10.-6.78,11.-7.81,12.-8.83,13.-9.28,14.-10.30,15.-11.32,16.-12.33,17.-12.79,18.-13.81,19.-14.82,20.-15.83};
const double poissonerror_plus[21] = {1.29-0.,2.75-1.,4.25-2.,5.30-3.,6.78-4.,7.81-5.,9.28-6.,10.30-7.,11.32-8.,12.79-9.,13.81-10.,14.82-11.,16.29-12.,17.30-13.,18.32-14.,19.32-15.,20.80-16.,21.81-17.,22.82-18.,23.82-19.,25.30-20.};

#endif //CONSTANTS_H
