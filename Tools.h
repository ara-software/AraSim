////////////////////////////////////////////////////////////////////////////////////////////////
//class Tools:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TOOLS_H_
#define TOOLS_H_


//#include "TSpline.h"
//#include "TH1.h"
/* #include "TGraph.h" */
/* #include "TH2F.h" */
#include "TRandom3.h" 
#include <iostream>
#include <fstream>

using namespace std;
using std::string;
using std::vector;

class TSpline5;
class TH1;
class TGraph;
class TH2F;
class TRandom3;


class Tools {


public:

  static double dMax(double,double);
  static double dMax(const double*,int);
  static double dvMax(const vector<double>);
  static double dsMax(TSpline5 *sp);
  static double dMin(const double*,int);
  static double dMinNotZero(const double*,int);
  static double dMin(double,double);
  static double getMaxMagnitude(vector<double> v);
  static int Getifreq(double freq,double freq_low,double freq_high,int n);
  static void InterpolateComplex(double *array, const int n);

  static void four1(double *data, const int isign,int nsize);
  static void realft(double *data, const int isign, int nsize);

  static void SWAP(double &a, double &b) // swaps two numbers
  {double dum=a; a=b; b=dum;}
  static void NormalTimeOrdering(const int n,double *volts);
  static void NormalTimeOrdering_InvT(const int n,double *volts);
  static void ShiftLeft(double *x,const int n,int ishift);  
  static void ShiftRight(double *x,const int n,int ishift);  
  static double GetFWHM(TH1 *h1);
  static void  MakeGraph(int n,double *time,double *volts,TGraph *&mygraph,TH2F *&h2, double scalex,double scaley,string xaxistitle,string yaxistitle);
  // Function declarations
  static double dDot(double*,double*,int);
  static void dCross(double*,double*,double*);
  static double dSquare(double*);
  static double Step(double x);
  static double dGetTheta(double*);
  static double dGetPhi(double*);
  static int WhichIsMax(double *x,int n);
  static int WhichIsMin(double *x,int n);

  static double dSum(double*,int);
  static int iSum(int*,int);
  static void Print(double*,int);
  static void Print(int*,int);
  static void Zero(double *anarray,int n);
  static void Zero(int *anarray,int n);
  static void GetNumbersAsStringArray(ifstream& fin, ofstream& fout,vector<string>& vnumbers, int nelements);
  static void GetNext2NumbersAsString(ifstream& fin,ofstream& fout,string& number1,string& number2, string& stherest);
  static void GetNextNumberAsString(ifstream& fin,ofstream& fout,string& number);

  static void SimpleLinearInterpolation(int n1, double *x1, double *y1, int n2, double *x2, double *y2 );

  static void SimpleLinearInterpolation_OutZero(int n1, double *x1, double *y1, int n2, double *x2, double *y2 );


  static void SimpleLinearInterpolation_extend_180cut(int n1, double *x1, double *y1, int n2, double *x2, double *y2 ); 


  static void SimpleLinearInterpolation_extend_PIcut(int n1, double *x1, double *y1, int n2, double *x2, double *y2 ); 

  static void SimpleLinearInterpolation_OutZero_180cut(int n1, double *x1, double *y1, int n2, double *x2, double *y2 ); 

  static void SimpleLinearInterpolation_OutZero_PIcut(int n1, double *x1, double *y1, int n2, double *x2, double *y2 ); 

  static double SimpleLinearInterpolation_extend_Single(int n1, double *x1, double *y1, double x2 );


  static  void get_random_rician(double signal, double signal_phase, double sigma, double& amplitude, double &phase);
  static  void get_circular_bivariate_normal_random_variable(double& rand_gauss_a, double& rand_gauss_b);

  static void Exchange(double &a, double &b);

protected:

};
#endif
