#include "Vector.h"
#include "TRandom3.h"
#include "Settings.h"
#include "Position.h"
#include "Primaries.h"
#include "secondaries.hh"
#include "Tools.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

#include "IceModel.h"

#include "TH1F.h"
#include "Constants.h"
#include "Settings.h"

using std::cout;
using std::stringstream;
using std::setprecision;
using std::accumulate;
using std::max_element;
using std::partial_sum;
using std::max;


// Vector Constants
static Vector x_axis = Vector(1,0,0);
static Vector y_axis = Vector(0,1,0);
static Vector z_axis = Vector(0,0,1);

Secondaries::~Secondaries() {

}



 Secondaries::Secondaries(Settings *settings1) {
	//For Total Tau Survival probability equation
	//n.b. not in SI units.
	////from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
	//Equation 16  &  used in Equation 30. 
	/////////////////////|Units////|Description////////////////////////
	B0=1.2*pow(10.,-6.); //| cm^2/g  |\,
	B1=0.16*pow(10.,-6.);//| cm^2/g  | }parameterization using a logarithmic dependence on energy for B,
	E0=pow(10.,10.);     //| GeV     |/   the tau elecromagnetic energy loss parameter.
	p=2.65;            //| g/cm^3  |Density of Standard Rock
	mT=1.777;	   //| GeV     |Mass of Tau
	cT=0.008693; 	   //| cm      |Tau Decay length (86.93 microMeters)
	                   //|         |
	Mn=1.672622E-24;   //| g       |nucleon/ proton mass in grams,also equal to 0.938 GeV. 
	A=1.;              //| none    |constant that sets the total probability to unity
	//these last two constanst from Connolly Calc 2011, used in d_dzPsurvNu().
  
  flavors[0]="nue";
  flavors[1]="numu";
  flavors[2]="nutau"; // the gps path of the anita-lite flight

  //SECONDARIES=1; // include secondary interactions
  //TAUDECAY=1; // include secondary interactions
  SECONDARIES = settings1->SECONDARIES;
  TAUDECAY = settings1->TAUDECAY;


    // reading in tauola data file for tau decays
  InitTauola();
  
  TAUFRAC=.5; //fraction of tau neutrino-cc current events where the primare interaction point is the first bang   

  tauolainfile.open("data/tau_decay_tauola.dat",ifstream::in);


  count_nfb=0;
  secondary_e_noncons=0;

  for (int i=0;i<7;i++) {
    Tools::Zero(dsdy_muon_brems[i],NPROB_MAX);
    Tools::Zero(dsdy_muon_epair[i],NPROB_MAX);
    Tools::Zero(dsdy_muon_pn[i],NPROB_MAX);
    
    Tools::Zero(y_muon_brems[i],NPROB_MAX);
    Tools::Zero(y_muon_epair[i],NPROB_MAX);
    Tools::Zero(y_muon_pn[i],NPROB_MAX);
    
    Tools::Zero(dsdy_tauon_brems[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_epair[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_pn[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_hadrdecay[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_edecay[i],NPROB_MAX);
    Tools::Zero(dsdy_tauon_mudecay[i],NPROB_MAX);
  


    Tools::Zero(y_tauon_brems[i],NPROB_MAX);
    Tools::Zero(y_tauon_epair[i],NPROB_MAX);
    Tools::Zero(y_tauon_pn[i],NPROB_MAX);
    Tools::Zero(y_tauon_hadrdecay[i],NPROB_MAX);
    Tools::Zero(y_tauon_edecay[i],NPROB_MAX);
    Tools::Zero(y_tauon_mudecay[i],NPROB_MAX);
  } //for (Tools::Zeroing)

  Tools::Zero(int_muon_brems,7);
  Tools::Zero(int_muon_epair,7);
  Tools::Zero(int_muon_pn,7);
  
  Tools::Zero(int_tauon_brems,7);
  Tools::Zero(int_tauon_epair,7);
  Tools::Zero(int_tauon_pn,7);
  Tools::Zero(int_tauon_hadrdecay,7);
  Tools::Zero(int_tauon_edecay,7);
  Tools::Zero(int_tauon_mudecay,7);

  // Read probability distributions for secondary interactions

  ReadSecondaries();

	
	
}//Secondaries Constructor

 void Secondaries::readData(string nuflavor,string secndryType, double (*y)[NPROB_MAX], double (*dsdy)[NPROB_MAX])
{
  
  stringstream senergy;
  
  ifstream ifile;
  string suffix=".vec";
  if(nuflavor=="tauon")
    suffix="_tau.vec";
  
  for(int index=0;index<7;index++)
    {senergy.str("");
      double energy=18+0.5*index;
      int precision=(index%2==0)?2:3;
      senergy << setprecision(precision) << energy;
      string path="secondary/"+nuflavor+"/dsdy_"+secndryType+"_1e"+senergy.str()+suffix;
      //cout << "openning file " << path.c_str() << endl;
      ifile.open(path.c_str());
      NPROB=0;
      while(!ifile.eof())
	{
	  ifile >> y[index][NPROB] >> dsdy[index][NPROB];
	  NPROB++;
	  if(NPROB>=NPROB_MAX)
	    {
	      // cerr << " ERROR in reading in y_muon_brem. \n";
	      break;
	    }
	 
	}
      ifile.close();
    }
  
}

 void Secondaries::ReadSecondaries() {
  // reading in data for secondary interactions
  
  cout<<"Reading in data on secondary interactions.\n";

  readData("muons","brems",y_muon_brems,dsdy_muon_brems);
  readData("muons","epair",y_muon_epair,dsdy_muon_epair);
  readData("muons","pn",y_muon_pn,dsdy_muon_pn);
  readData("tauon","brems",y_tauon_brems,dsdy_tauon_brems);
  readData("tauon","epair",y_tauon_epair,dsdy_tauon_epair);
  readData("tauon","pn",y_tauon_pn,dsdy_tauon_pn);
  readData("tauon","hadrdecay",y_tauon_hadrdecay,dsdy_tauon_hadrdecay);
  readData("tauon","edecay",y_tauon_edecay,dsdy_tauon_edecay);
  readData("tauon","mudecay",y_tauon_mudecay,dsdy_tauon_mudecay);
  //cout << "NPROB=" << NPROB << ",  NPROB_MAX=" << NPROB_MAX << endl;
 for(int j=0;j<7;j++) {
    // integrating prob. distributions.
    int_muon_brems[j]=accumulate(dsdy_muon_brems[j],dsdy_muon_brems[j]+NPROB_MAX,0.);//very important to keep the initial value the same type as the elements type
    int_muon_epair[j]=accumulate(dsdy_muon_epair[j],dsdy_muon_epair[j]+NPROB_MAX,0.);
    int_muon_pn[j]=accumulate(dsdy_muon_pn[j],dsdy_muon_pn[j]+NPROB_MAX,0.);
    int_tauon_brems[j]=accumulate(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+NPROB_MAX,0.);
    int_tauon_epair[j]=accumulate(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+NPROB_MAX,0.);
    int_tauon_pn[j]=accumulate(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+NPROB_MAX,0.);
    int_tauon_hadrdecay[j]=accumulate(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+NPROB_MAX,0.);
    int_tauon_edecay[j]=accumulate(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+NPROB_MAX,0.);
    int_tauon_mudecay[j]=accumulate(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+NPROB_MAX,0.);
    
    // maximum value of prob. dist.
    max_muon_brems=*max_element(dsdy_muon_brems[j],dsdy_muon_brems[j]+NPROB_MAX);
    //cout << "max_muon_brems=" << max_muon_brems << endl;//fenfang
    max_muon_epair=*max_element(dsdy_muon_epair[j],dsdy_muon_epair[j]+NPROB_MAX);
    max_muon_pn=*max_element(dsdy_muon_pn[j],dsdy_muon_pn[j]+NPROB_MAX);   
    max_tauon_brems=*max_element(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+NPROB_MAX);
    max_tauon_epair=*max_element(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+NPROB_MAX);
    max_tauon_pn=*max_element(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+NPROB_MAX);
    max_tauon_hadrdecay=*max_element(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+NPROB_MAX);
    max_tauon_edecay=*max_element(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+NPROB_MAX);
    max_tauon_mudecay=*max_element(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+NPROB_MAX);
     
    // minimum value of prob. dist.
    min_muon_brems=Tools::dMinNotZero(dsdy_muon_brems[j],NPROB_MAX);
    min_muon_epair=Tools::dMinNotZero(dsdy_muon_epair[j],NPROB_MAX);
    min_muon_pn=Tools::dMinNotZero(dsdy_muon_pn[j],NPROB_MAX);   
    min_tauon_brems=Tools::dMinNotZero(dsdy_tauon_brems[j],NPROB_MAX);
    min_tauon_epair=Tools::dMinNotZero(dsdy_tauon_epair[j],NPROB_MAX);
    min_tauon_pn=Tools::dMinNotZero(dsdy_tauon_pn[j],NPROB_MAX);
    min_tauon_hadrdecay=Tools::dMinNotZero(dsdy_tauon_hadrdecay[j],NPROB_MAX);
    min_tauon_edecay=Tools::dMinNotZero(dsdy_tauon_edecay[j],NPROB_MAX);
    min_tauon_mudecay=Tools::dMinNotZero(dsdy_tauon_mudecay[j],NPROB_MAX);
     
    if (min_muon_brems<=0)
      cout << "Minimum probability is <=0!\n";
    
    partial_sum(dsdy_muon_brems[j],dsdy_muon_brems[j]+NPROB_MAX,y_cumulative_muon_brems[j]);
    partial_sum(dsdy_muon_epair[j],dsdy_muon_epair[j]+NPROB_MAX,y_cumulative_muon_epair[j]);
    partial_sum(dsdy_muon_pn[j],dsdy_muon_pn[j]+NPROB_MAX,y_cumulative_muon_pn[j]);
    partial_sum(dsdy_tauon_brems[j],dsdy_tauon_brems[j]+NPROB_MAX,y_cumulative_tauon_brems[j]);
    partial_sum(dsdy_tauon_epair[j],dsdy_tauon_epair[j]+NPROB_MAX,y_cumulative_tauon_epair[j]);
    partial_sum(dsdy_tauon_pn[j],dsdy_tauon_pn[j]+NPROB_MAX,y_cumulative_tauon_pn[j]);
    partial_sum(dsdy_tauon_hadrdecay[j],dsdy_tauon_hadrdecay[j]+NPROB_MAX,y_cumulative_tauon_hadrdecay[j]);
    partial_sum(dsdy_tauon_mudecay[j],dsdy_tauon_mudecay[j]+NPROB_MAX,y_cumulative_tauon_mudecay[j]);
    partial_sum(dsdy_tauon_edecay[j],dsdy_tauon_edecay[j]+NPROB_MAX,y_cumulative_tauon_edecay[j]);
     
    for (int i=0;i<NPROB_MAX;i++) {
       y_cumulative_muon_brems[j][i]      /= y_cumulative_muon_brems[j][NPROB_MAX-1];
       y_cumulative_muon_epair[j][i]      /= y_cumulative_muon_epair[j][NPROB_MAX-1];
       y_cumulative_muon_pn[j][i]         /= y_cumulative_muon_pn[j][NPROB_MAX-1];
       y_cumulative_tauon_brems[j][i]     /= y_cumulative_tauon_brems[j][NPROB_MAX-1];
       y_cumulative_tauon_epair[j][i]     /= y_cumulative_tauon_epair[j][NPROB_MAX-1];
       y_cumulative_tauon_pn[j][i]        /= y_cumulative_tauon_pn[j][NPROB_MAX-1];
       y_cumulative_tauon_hadrdecay[j][i] /= y_cumulative_tauon_hadrdecay[j][NPROB_MAX-1];
       y_cumulative_tauon_mudecay[j][i]   /= y_cumulative_tauon_mudecay[j][NPROB_MAX-1];
       y_cumulative_tauon_edecay[j][i]    /= y_cumulative_tauon_edecay[j][NPROB_MAX-1];
    } //for

  }
  cout<<"Finished reading secondary interaction data.\n"; 
  
  /*string istring;
  char buffer[50];
  int n;  // counter
  

  ifstream ifile;
  int index;
  cout<<"Reading in data on secondary interactions.\n";
  for (int j=0;j<2;j++) {
    // for each energy
    for (int i=18;i<=21;i++) {
      if (!(i==21 && j==1)) {
      if (j==0)
	n=sprintf (buffer, "secondary/muons/dsdy_brems_1e%d.vec", i);
      if (j==1)
	n=sprintf (buffer, "secondary/muons/dsdy_brems_1e%d.5.vec", i);

      istring=buffer;

      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_muon_brems[index][NPROB] >> dsdy_muon_brems[index][NPROB];
	NPROB++;
      }
      ifile.close();   
      }
    }
  }

  for (int j=0;j<2;j++) {
    for (int i=18;i<=21;i++) {
      if (!(i==21 && j==1)) {
      if (j==0)
	n=sprintf (buffer, "secondary/muons/dsdy_epair_1e%d.vec", i);
      if (j==1)
	n=sprintf (buffer, "secondary/muons/dsdy_epair_1e%d.5.vec", i);

      istring=buffer;
      
      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_muon_epair[index][NPROB] >> dsdy_muon_epair[index][NPROB];
	NPROB++;
      }
      ifile.close();   
      }
    }
  }

  for (int j=0;j<2;j++) {
    for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {      
      if (j==0)
	n=sprintf (buffer, "secondary/muons/dsdy_pn_1e%d.vec", i);
      if (j==1)
	n=sprintf (buffer, "secondary/muons/dsdy_pn_1e%d.5.vec", i);

      istring=buffer;
      
      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_muon_pn[index][NPROB] >> dsdy_muon_pn[index][NPROB];
	NPROB++;
      }
      ifile.close();   
       }
    }
  }

   for (int j=0;j<2;j++) {
    for (int i=18;i<=21;i++) {
      if (!(i==21 && j==1)) {
      if (j==0)
	n=sprintf (buffer, "secondary/tauon/dsdy_brems_1e%d_tau.vec", i);
      if (j==1)
	n=sprintf (buffer, "secondary/tauon/dsdy_brems_1e%d.5_tau.vec", i);

      istring=buffer;
      
      ifstream ifile;
      ifile.open(istring.c_str());
      NPROB=0;
      index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
      while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	ifile >> y_tauon_brems[index][NPROB] >> dsdy_tauon_brems[index][NPROB];
	NPROB++;
      }
      ifile.close();   
      }
    }
   }

   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {       

       if (j==0)
	 n=sprintf (buffer, "secondary/tauon/dsdy_epair_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "secondary/tauon/dsdy_epair_1e%d.5_tau.vec", i);
       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_epair[index][NPROB] >> dsdy_tauon_epair[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {  
       if (j==0)
	 n=sprintf (buffer, "secondary/tauon/dsdy_pn_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "secondary/tauon/dsdy_pn_1e%d.5_tau.vec", i);
       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_pn[index][NPROB] >> dsdy_tauon_pn[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {

       if (j==0)
	 n=sprintf (buffer, "secondary/tauon/dsdy_hadrdecay_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "secondary/tauon/dsdy_hadrdecay_1e%d.5_tau.vec", i);       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_hadrdecay[index][NPROB] >> dsdy_tauon_hadrdecay[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       
       if (!(i==21 && j==1)) {
       if (j==0)
	 n=sprintf (buffer, "secondary/tauon/dsdy_edecay_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "secondary/tauon/dsdy_edecay_1e%d.5_tau.vec", i);       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_edecay[index][NPROB] >> dsdy_tauon_edecay[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
   for (int j=0;j<2;j++) {
     for (int i=18;i<=21;i++) {
       if (!(i==21 && j==1)) {
       if (j==0)
	 n=sprintf (buffer, "secondary/tauon/dsdy_mudecay_1e%d_tau.vec", i);
       if (j==1)
	 n=sprintf (buffer, "secondary/tauon/dsdy_mudecay_1e%d.5_tau.vec", i);       
       istring=buffer;
       
       ifstream ifile;
       ifile.open(istring.c_str());
       NPROB=0;
       index=2*(i-18)+j;
      if (index>=7)
	cout << "ERROR in reading in y_muon_brems.\n";
       while (!ifile.eof()) {
	if (NPROB>=NPROB_MAX)
	  cout << "ERROR in reading in y_muon_brems.\n";
	 ifile >> y_tauon_mudecay[index][NPROB] >> dsdy_tauon_mudecay[index][NPROB];
	 NPROB++;
       }
       ifile.close();   
       }
     }
   }
 
   // for filling vectors with y values distributed so that they follow
   // dsdy distributions.

   for (int j=0;j<7;j++) {

     int_muon_brems[j]=0;
     int_muon_epair[j]=0;
     int_muon_pn[j]=0;   
     int_tauon_brems[j]=0;
     int_tauon_epair[j]=0;
     int_tauon_pn[j]=0;
     int_tauon_hadrdecay[j]=0;
     int_tauon_edecay[j]=0;
     int_tauon_mudecay[j]=0;
     
     // integrating prob. distributions.
     for (int i=0;i<NPROB_MAX;i++) {
       int_muon_brems[j]+=dsdy_muon_brems[j][i];
       //cout << "int_muon_brems is " << int_muon_brems[j] << "\n";
       int_muon_epair[j]+=dsdy_muon_epair[j][i];
       int_muon_pn[j]+=dsdy_muon_pn[j][i];
       int_tauon_brems[j]+=dsdy_tauon_brems[j][i];
       int_tauon_epair[j]+=dsdy_tauon_epair[j][i];
       int_tauon_pn[j]+=dsdy_tauon_pn[j][i];
       int_tauon_hadrdecay[j]+=dsdy_tauon_hadrdecay[j][i];
       int_tauon_edecay[j]+=dsdy_tauon_edecay[j][i];
       int_tauon_mudecay[j]+=dsdy_tauon_mudecay[j][i];
     }

     // maximum value of prob. dist. 
     max_muon_brems=dMax(dsdy_muon_brems[j],NPROB_MAX);
     max_muon_epair=dMax(dsdy_muon_epair[j],NPROB_MAX);
     max_muon_pn=dMax(dsdy_muon_pn[j],NPROB_MAX);   
     max_tauon_brems=dMax(dsdy_tauon_brems[j],NPROB_MAX);
     max_tauon_epair=dMax(dsdy_tauon_epair[j],NPROB_MAX);
     max_tauon_pn=dMax(dsdy_tauon_pn[j],NPROB_MAX);
     max_tauon_hadrdecay=dMax(dsdy_tauon_hadrdecay[j],NPROB_MAX);
     max_tauon_edecay=dMax(dsdy_tauon_edecay[j],NPROB_MAX);
     max_tauon_mudecay=dMax(dsdy_tauon_mudecay[j],NPROB_MAX);
     
     // minimum value of prob. dist.
     min_muon_brems=Tools::dMinNotZero(dsdy_muon_brems[j],NPROB_MAX);
     min_muon_epair=Tools::dMinNotZero(dsdy_muon_epair[j],NPROB_MAX);
     min_muon_pn=Tools::dMinNotZero(dsdy_muon_pn[j],NPROB_MAX);   
     min_tauon_brems=Tools::dMinNotZero(dsdy_tauon_brems[j],NPROB_MAX);
     min_tauon_epair=Tools::dMinNotZero(dsdy_tauon_epair[j],NPROB_MAX);
     min_tauon_pn=Tools::dMinNotZero(dsdy_tauon_pn[j],NPROB_MAX);
     min_tauon_hadrdecay=Tools::dMinNotZero(dsdy_tauon_hadrdecay[j],NPROB_MAX);
     min_tauon_edecay=Tools::dMinNotZero(dsdy_tauon_edecay[j],NPROB_MAX);
     min_tauon_mudecay=Tools::dMinNotZero(dsdy_tauon_mudecay[j],NPROB_MAX);
     
     if (min_muon_brems<=0)
       cout << "Minimum probability is <=0!\n";

     
     // for each y bin in dsdy curve, fill vector y_muon_brem with as
     // many of y's as you need to get the right distribution.
     for (int i=0;i<NPROB_MAX;i++) {
       y_cumulative_muon_brems[j][i]      = dSum(dsdy_muon_brems[j],i+1);
       y_cumulative_muon_epair[j][i]      = dSum(dsdy_muon_epair[j],i+1);
       y_cumulative_muon_pn[j][i]         = dSum(dsdy_muon_pn[j],i+1);
       y_cumulative_tauon_brems[j][i]     = dSum(dsdy_tauon_brems[j],i+1);
       y_cumulative_tauon_epair[j][i]     = dSum(dsdy_tauon_epair[j],i+1);
       y_cumulative_tauon_pn[j][i]        = dSum(dsdy_tauon_pn[j],i+1);
       y_cumulative_tauon_hadrdecay[j][i] = dSum(dsdy_tauon_hadrdecay[j],i+1);
       y_cumulative_tauon_mudecay[j][i]   = dSum(dsdy_tauon_mudecay[j],i+1);
       y_cumulative_tauon_edecay[j][i]    = dSum(dsdy_tauon_edecay[j],i+1);
     } //for

     // normalize the distributions
     for (int i=0;i<NPROB_MAX;i++) {
       y_cumulative_muon_brems[j][i]      /= y_cumulative_muon_brems[j][NPROB_MAX-1];
       y_cumulative_muon_epair[j][i]      /= y_cumulative_muon_epair[j][NPROB_MAX-1];
       y_cumulative_muon_pn[j][i]         /= y_cumulative_muon_pn[j][NPROB_MAX-1];
       y_cumulative_tauon_brems[j][i]     /= y_cumulative_tauon_brems[j][NPROB_MAX-1];
       y_cumulative_tauon_epair[j][i]     /= y_cumulative_tauon_epair[j][NPROB_MAX-1];
       y_cumulative_tauon_pn[j][i]        /= y_cumulative_tauon_pn[j][NPROB_MAX-1];
       y_cumulative_tauon_hadrdecay[j][i] /= y_cumulative_tauon_hadrdecay[j][NPROB_MAX-1];
       y_cumulative_tauon_mudecay[j][i]   /= y_cumulative_tauon_mudecay[j][NPROB_MAX-1];
       y_cumulative_tauon_edecay[j][i]    /= y_cumulative_tauon_edecay[j][NPROB_MAX-1];
     } //for
     
   }

   cout<<"Finished reading secondary interaction data.\n";
  */
} //end method ReadSecondaries


 void Secondaries::GetSecondaries(Settings *settings1,string nuflavor,double plepton,double &em_secondaries_max,double &had_secondaries_max,int &n_interactions,TH1F *hy) {


  em_secondaries_max=0.;
  had_secondaries_max=0.;

  int i=(int)((log10(plepton)-18.)*2.);
  if (i>6)
    i=6;
  if (i<0)
    i=0;

  int n_brems,n_epair,n_pn; // number of interactions of each type.
  int index_y; // index along the horizontal axis of ped's plots
  double rnd1=1000.;
  double rnd2=1000.;  // random numbers for throwing at dart board
  double y = 0; // inelasticity
 
  string whichtype; // which type of interaction corresponds to that index
  


  if (nuflavor=="numu") {   
    n_brems=gRandom->Poisson(int_muon_brems[i]); // pick number of brem interactions
    n_epair=gRandom->Poisson(int_muon_epair[i]); // # of pair production
    n_pn=gRandom->Poisson(int_muon_pn[i]); // # photonuclear interactions   
    
    n_interactions+=(n_brems+n_epair+n_pn);	


    for (int j=0;j<n_brems+n_epair+n_pn;j++) {
      rnd1=gRandom->Rndm();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn))
	whichtype="brems";
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn))
	whichtype="epair";
      else
	whichtype="pn";



      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;

      if (whichtype=="brems") {	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_brems[i],NPROB,rnd1,y);
      }
      else if (whichtype=="epair") {	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_epair[i],NPROB,rnd1,y);	
      }
      else if (whichtype=="pn") {
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_pn[i],NPROB,rnd1,y);
      }
     
      if (y*plepton>max(em_secondaries_max,had_secondaries_max)) {  // if this is the largest interaction for this event so far
	if (whichtype=="brems" || whichtype=="epair") {  // save it
	  em_secondaries_max=y*plepton;

	}
	if (whichtype=="pn") 
	  had_secondaries_max=y*plepton;
	 
		
      }
    } // loop over secondary interactions
  } // end if it was a muon neutrino
  if (nuflavor=="nutau") {
    n_brems=gRandom->Poisson(int_tauon_brems[i]);
    n_epair=gRandom->Poisson(int_tauon_epair[i]);
    n_pn=gRandom->Poisson(int_tauon_pn[i]);

    n_interactions+=(n_brems+n_epair+n_pn); // increment number of secondary interactions.

    for (int j=0;j<n_brems+n_epair+n_pn;j++) { // loop over secondary interactions. 
      
      rnd1=gRandom->Rndm();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn))
	whichtype="brems";
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn))
	whichtype="epair";
      else
	whichtype="pn";
  
      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;

      if (whichtype=="brems") {  // bremstrahlung interaction
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_brems[i],NPROB,rnd1,y);
      }
      if (whichtype=="epair") { // pair production
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_epair[i],NPROB,rnd1,y);
      }
      if (whichtype=="pn") {
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_pn[i],NPROB,rnd1,y);
      }

      if (settings1->HIST==1 && !settings1->ONLYFINAL && hy->GetEntries()<settings1->HIST_MAX_ENTRIES)
	hy->Fill(y);
      if (y*plepton>max(em_secondaries_max,had_secondaries_max)) { // if this is the biggest secondary signal yet,
	if (whichtype=="brems" || whichtype=="epair") // save it.
	  em_secondaries_max=y*plepton;
	if (whichtype=="pn")
	  had_secondaries_max=y*plepton;
      }
    }
   

    if (TAUDECAY) {
      n_interactions++; // increment number of interactions, for plotting.

      rnd1=gRandom->Rndm();
      if (rnd1<0.65011)  // pick which type of decay it is.
	whichtype="hadrdecay";
      if (rnd1>=0.65011 && rnd1<0.8219)
	whichtype="mudecay";
      if (rnd1>=0.8219)
	whichtype="edecay";
           
      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;     
      
      if (whichtype=="hadrdecay") { // hadronic decay
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_hadrdecay[i],NPROB,rnd1,y);	
      }
      else if (whichtype=="edecay") { // e decay	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_edecay[i],NPROB,rnd1,y);
      }
      else if (whichtype=="mudecay") { // mu decay
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_mudecay[i],NPROB,rnd1,y);
      }
      
     
      if (y*plepton>max(em_secondaries_max, had_secondaries_max)) {  // if this is the biggest interaction yet,    
	if (whichtype=="edecay") // save it.
	  em_secondaries_max=y*plepton;
	if (whichtype=="hadrdecay")
	  had_secondaries_max=y*plepton;
      } //if     
    } //if (TAUDECAY)
  } //if (nutau)

} //GetSecondaries


// Eugene version. (without histogram hy)
 void Secondaries::GetSecondaries(Settings *settings1,string nuflavor,double plepton,double &em_secondaries_max,double &had_secondaries_max,int &n_interactions ) {


  em_secondaries_max=0.;
  had_secondaries_max=0.;

  int i=(int)((log10(plepton)-18.)*2.);
  if (i>6)
    i=6;
  if (i<0)
    i=0;

  int n_brems,n_epair,n_pn; // number of interactions of each type.
  int index_y; // index along the horizontal axis of ped's plots
  double rnd1=1000.;
  double rnd2=1000.;  // random numbers for throwing at dart board
  double y = 0; // inelasticity
 
  string whichtype; // which type of interaction corresponds to that index
  


  if (nuflavor=="numu") {   
    n_brems=gRandom->Poisson(int_muon_brems[i]); // pick number of brem interactions
    n_epair=gRandom->Poisson(int_muon_epair[i]); // # of pair production
    n_pn=gRandom->Poisson(int_muon_pn[i]); // # photonuclear interactions   
    
    n_interactions+=(n_brems+n_epair+n_pn);	


    for (int j=0;j<n_brems+n_epair+n_pn;j++) {
      rnd1=gRandom->Rndm();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn))
	whichtype="brems";
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn))
	whichtype="epair";
      else
	whichtype="pn";



      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;

      if (whichtype=="brems") {	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_brems[i],NPROB,rnd1,y);
      }
      else if (whichtype=="epair") {	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_epair[i],NPROB,rnd1,y);	
      }
      else if (whichtype=="pn") {
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_muon_pn[i],NPROB,rnd1,y);
      }
     
      if (y*plepton>max(em_secondaries_max,had_secondaries_max)) {  // if this is the largest interaction for this event so far
	if (whichtype=="brems" || whichtype=="epair") {  // save it
	  em_secondaries_max=y*plepton;

	}
	if (whichtype=="pn") 
	  had_secondaries_max=y*plepton;
	 
		
      }
    } // loop over secondary interactions
  } // end if it was a muon neutrino
  if (nuflavor=="nutau") {
    n_brems=gRandom->Poisson(int_tauon_brems[i]);
    n_epair=gRandom->Poisson(int_tauon_epair[i]);
    n_pn=gRandom->Poisson(int_tauon_pn[i]);

    n_interactions+=(n_brems+n_epair+n_pn); // increment number of secondary interactions.

    for (int j=0;j<n_brems+n_epair+n_pn;j++) { // loop over secondary interactions. 
      
      rnd1=gRandom->Rndm();
      if (rnd1<=(double)n_brems/(double)(n_brems+n_epair+n_pn))
	whichtype="brems";
      else if (rnd1<=(double)(n_brems+n_epair)/(double)(n_brems+n_epair+n_pn))
	whichtype="epair";
      else
	whichtype="pn";
  
      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;

      if (whichtype=="brems") {  // bremstrahlung interaction
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_brems[i],NPROB,rnd1,y);
      }
      if (whichtype=="epair") { // pair production
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_epair[i],NPROB,rnd1,y);
      }
      if (whichtype=="pn") {
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_pn[i],NPROB,rnd1,y);
      }

//--------------------------------------------------
//       if (settings1->HIST==1 && !settings1->ONLYFINAL && hy->GetEntries()<settings1->HIST_MAX_ENTRIES)
// 	hy->Fill(y);
//-------------------------------------------------- 
      if (y*plepton>max(em_secondaries_max,had_secondaries_max)) { // if this is the biggest secondary signal yet,
	if (whichtype=="brems" || whichtype=="epair") // save it.
	  em_secondaries_max=y*plepton;
	if (whichtype=="pn")
	  had_secondaries_max=y*plepton;
      }
    }
   

    if (TAUDECAY) {
      n_interactions++; // increment number of interactions, for plotting.

      rnd1=gRandom->Rndm();
      if (rnd1<0.65011)  // pick which type of decay it is.
	whichtype="hadrdecay";
      if (rnd1>=0.65011 && rnd1<0.8219)
	whichtype="mudecay";
      if (rnd1>=0.8219)
	whichtype="edecay";
           
      rnd1=1000.;
      rnd2=1000.;  // random numbers for throwing at dart board
      index_y=0;     
      
      if (whichtype=="hadrdecay") { // hadronic decay
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_hadrdecay[i],NPROB,rnd1,y);	
      }
      else if (whichtype=="edecay") { // e decay	
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_edecay[i],NPROB,rnd1,y);
      }
      else if (whichtype=="mudecay") { // mu decay
	rnd1=gRandom->Rndm();
	Picky(y_cumulative_tauon_mudecay[i],NPROB,rnd1,y);
      }
      
     
      if (y*plepton>max(em_secondaries_max, had_secondaries_max)) {  // if this is the biggest interaction yet,    
	if (whichtype=="edecay") // save it.
	  em_secondaries_max=y*plepton;
	if (whichtype=="hadrdecay")
	  had_secondaries_max=y*plepton;
      } //if     
    } //if (TAUDECAY)
  } //if (nutau)

} //GetSecondaries




// from icemc with tau mode (just removed hy and inu)
 int Secondaries::GetEMFrac(Settings *settings1,string nuflavor,
		     string current,
		     string taudecay,	      
		     double y,
		     //TH1F *hy,
                     double pnu,				  
                     //int inu,
		     double& emfrac,
		     double& hadfrac,
                     int& n_interactions, int taumodes1,double ptauf) {


  if (current=="cc")
    plepton=(1.-y)*pnu;
  else
    plepton=0.;
  
  if (nuflavor=="nue" && current=="cc") {
    emfrac=1.-y;
    hadfrac=y;
  }
  else if(nuflavor=="numu" && current=="cc") {
    emfrac=1.E-10;
    hadfrac=y;
  }
  else if(nuflavor=="nutau" && current=="cc") {
    // behaves like a muon
    if(taumodes1 ==1){//taumodes==1; tau created somewhere in rock and decays at posnu.
      emfrac = ptauf/(3.*pnu); //Neutrino creates tau (1-y), tau decays into electron, which has ~1/3 the energy.
      //cout<<"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ \n";
      //cout <<"pnu is "<<pnu<<"\n";
      //cout <<"ptauf inside get frac is "<<ptauf<<"\n";
      //cout <<"emfrac inside getfrac is "<<emfrac<<"\n";
      hadfrac = 1E-10;
    }
    else if (taumodes1 == 0){
    emfrac=1.E-10;
        hadfrac=y;
    }
    


  }
  else if (current=="nc") {
    emfrac=1.E-10;
    hadfrac=y;
  }


  em_secondaries_max =emfrac; // initialize search for maximum signal among primary, secondary interactions.
  had_secondaries_max=hadfrac;

  
  
  if (SECONDARIES==1 && current=="cc" && settings1->FORSECKEL!=1) {

    while (1) {

      GetSecondaries(settings1,nuflavor,plepton,em_secondaries_max,had_secondaries_max,n_interactions); // find how much em and hadronic energies comes from secondary interactions.  keep picking until you get a bunch of secondary interactions that conserve energy

      if (em_secondaries_max+had_secondaries_max<=plepton*(1.+1.E-5)) // if conserves energy, break.
	break;
      else {
	secondary_e_noncons++; //Record how many times we come up with something that doesn't conserve energy
	em_secondaries_max=emfrac;
	had_secondaries_max=hadfrac;
      } //else
    } //while(1)

    if ((em_secondaries_max+had_secondaries_max)>(emfrac+hadfrac)*pnu) { // if maximum signal from secondaries is larger than
                                                                         // signal from primary interaction
      emfrac=em_secondaries_max/pnu; // then use that one.
      hadfrac=had_secondaries_max/pnu;
      if (emfrac <= 1.E-10)
	emfrac=1.E-10;
      if (hadfrac<= 1.E-10)
	hadfrac=1.E-10;
    } //if
  } //if (charged current, secondaries on)

//--------------------------------------------------
//   if (nuflavor=="numu" && current=="cc" && n_interactions==0)
//     cout << "Look at this one.  inu is " << inu << "\n";
//-------------------------------------------------- 
  


  if ((y<0 || y>1) && y != -999.) 
    cout <<  "illegal y=" << y << "\n";
          
  if (emfrac+hadfrac>1.00001) {
    cout << "error emfrac,hadfrac=" << emfrac << " " << hadfrac << " " << emfrac+hadfrac << "\n";
    cout << "nuflavor,taudecay=" << nuflavor << " " << taudecay << "\n";
  } //if
  
  return 1;

} //GetEMFrac






// Eugene added version (no hy, inu)
 int Secondaries::GetEMFrac(Settings *settings1,string nuflavor,
		     string current,
		     string taudecay,	      
		     double y,
		     double pnu,				  
		     double& emfrac,
		     double& hadfrac,
		     int& n_interactions) {


  if (current=="cc")
    plepton=(1.-y)*pnu;
  else
    plepton=0.;
  
  if (nuflavor=="nue" && current=="cc") {
    emfrac=1.-y;
    hadfrac=y;
  }
  else if(nuflavor=="numu" && current=="cc") {
    emfrac=1.E-10;
    hadfrac=y;
  }
  else if(nuflavor=="nutau" && current=="cc") {
    // behaves like a muon
    

    emfrac=1.E-10;
    hadfrac=y;

    


  }
  else if (current=="nc") {
    emfrac=1.E-10;
    hadfrac=y;
  }


  em_secondaries_max =emfrac; // initialize search for maximum signal among primary, secondary interactions.
  had_secondaries_max=hadfrac;

  
  //cout<<"settings1->FORSECKEL = "<<settings1->FORSECKEL<<endl;
  
  if (SECONDARIES==1 && current=="cc" && settings1->FORSECKEL!=1) {

      cout<<"look for second interaction em, had frac"<<endl;

    while (1) {

      GetSecondaries(settings1,nuflavor,plepton,em_secondaries_max,had_secondaries_max,n_interactions ); // find how much em and hadronic energies comes from secondary interactions.  keep picking until you get a bunch of secondary interactions that conserve energy

      if (em_secondaries_max+had_secondaries_max<=plepton*(1.+1.E-5)) // if conserves energy, break.
	break;
      else {
	secondary_e_noncons++; //Record how many times we come up with something that doesn't conserve energy
	em_secondaries_max=emfrac;
	had_secondaries_max=hadfrac;
      } //else
    } //while(1)

    if ((em_secondaries_max+had_secondaries_max)>(emfrac+hadfrac)*pnu) { // if maximum signal from secondaries is larger than
                                                                         // signal from primary interaction
      emfrac=em_secondaries_max/pnu; // then use that one.
      hadfrac=had_secondaries_max/pnu;
      if (emfrac <= 1.E-10)
	emfrac=1.E-10;
      if (hadfrac<= 1.E-10)
	hadfrac=1.E-10;
    } //if
  } //if (charged current, secondaries on)

//--------------------------------------------------
//   if (nuflavor=="numu" && current=="cc" && n_interactions==0)
//     cout << "Look at this one.  inu is " << inu << "\n";
//-------------------------------------------------- 
  


  if ((y<0 || y>1) && y != -999.) 
    cout <<  "illegal y=" << y << "\n";
          
  if (emfrac+hadfrac>1.00001) {
    cout << "error emfrac,hadfrac=" << emfrac << " " << hadfrac << " " << emfrac+hadfrac << "\n";
    cout << "nuflavor,taudecay=" << nuflavor << " " << taudecay << "\n";
  } //if
  
  return 1;

} //GetEMFrac












//----------------------------------------------------------
//InitTauola()
//Initializes the tau decay information

 void Secondaries::InitTauola() {
  for(int k=0;k<5;k++)
    tauolainfile >> tauola[0][k];
  for(int i=1;i<N_TAUOLA;i++)
    for(int j=0;j<6;j++)
      tauolainfile >> tauola[i][j];

  return;
}//InitTauola

void Secondaries::GetTauDecay(string nuflavor,string current,string& taudecay,double& emfrac_db,double& hadfrac_db) {

  if (!(nuflavor=="nutau" || current=="cc" || interestedintaus))
    return;

  // if nu_tau choose tau decay type
  
  double rnd = gRandom->Rndm();
  int decay = static_cast<int>(rnd*(N_TAUOLA-2)+1);
  
  hadfrac_db=tauola[decay][3];
  emfrac_db=tauola[decay][4];
  
  if(tauola[decay][1]!=0)
    taudecay="m";
  else
    taudecay="e";
  

  if(taudecay=="m")
    secondbang=false; //for all muon decays, the interaction point chosen is the neutrino interaction since we don't detect the decay if
  //the tau decays into a muon.
  else {
    double rnd=gRandom->Rndm();
    if(rnd>TAUFRAC) {
      secondbang=true;
      count_nfb++;
    } else
      secondbang=false;
  }
  

} //GetTauDecay

//-----------------------------------------------------
//GetEMFracDB()
//Gets the emfrac_db and hadfrac_db for a doublebang

 void Secondaries::GetEMFracDB(double& emfrac_db, double& hadfrac_db) {


  double rnd = gRandom->Rndm();
  int decay = static_cast<int>(rnd*(N_TAUOLA-2)+1);

  hadfrac_db=tauola[decay][3];
  emfrac_db=tauola[decay][4];

  return;
}//GetEMFracDB






//------------------------------------------------------
//GetDBViewAngle()
//Gets the viewangle of the second bang

 double Secondaries::GetDBViewAngle(const Vector &refr, const Vector &nnu) {

  return ((nnu.ChangeCoord(refr)).Angle(z_axis));

}//GetDBViewAngle

//------------------------------------------------------
//GetFirstBang()
//Gets the position of the first bang when the interaction point is the tau decay point

//  void Secondaries::GetFirstBang(const Position &r_in, const Vector &nnu, Position &posnu, double len_int_kgm2, double chord, double &nuentrancelength) {
  
//   double weightbang;
//   double junk1;
//   double junk2;
//   double junk3;
//   int junk4,junk5,junk6;
//   double myair=0;

//   Vector r_out = r_in + chord*nnu;

//   antarctica->Getchord(len_int_kgm2,r_in,r_out,
// 		  junk1,weightbang,junk2,myair,junk3,junk4,junk5,junk6);
//   double r1,r2;
//   if(weightbang>.999)
//     r2=gRandom->Rndm()*chord;
//   else {
//     do {
//       r1=gRandom->Rndm();
//       r2=gRandom->Rndm()*chord;
//       r_out = r_in + r2*nnu;
//       antarctica->Getchord(len_int_kgm2,r_in,r_out,
// 		      junk1,weightbang,junk2,myair,junk3,junk4,junk5,junk6);
//     }
//     while(r1>1-weightbang);
//   }
//   posnu = r_in + r2*nnu;
//   nuentrancelength=r2;

//   return;
// }//GetFirstBang

//---------------------------------------------------------
//NFBWeight()
//Gets the weight of the tau decay for second bang events
  double Secondaries::NFBWeight(double ptau, double taulength) {
  
  double gamma=ptau/MTAU;
  double D=TAUDECAY_TIME*CLIGHT*gamma;

  return exp(-taulength/D);

}
void Secondaries::Picky(double *y_cumulative,int NPROB,double rnd,double& y) {
  for (int i=0;i<NPROB;i++) {
    if (y_cumulative[i]<=rnd && y_cumulative[i+1]>rnd) {
      y=(double)i/(double)NPROB;
      continue; // once you found the right bin, stop looping.
    } //if
  } //for
} //Picky

//-------------------------------------------------------------
double Secondaries::GetTauWeight(Primaries*primary1, Settings*settings1, double pnu, int nu_nubar, int currentint, double ptau_final, double Distance){
	double ptaui, ptaui_GeV,dptaui;//ptauFinal_GeV,
	double y, yweight;//d_dzPnu_surv;
	double Dcm, zcm, dzcm;//all in cm
	double  zdistance; //m
	Dcm=Distance*100.;//convert meters to cm.
	zdistance=0.; //m	
	double TauWeight;//total tau survival probability.
	double powerStep;
	double prob_at_zcm;
	
	dzcm=Distance;//will always have 10^5 steps.
	int numSteps_z=Dcm/dzcm; //use ints to keep precision
	
	//for(zcm=0.; zcm<=Dcm; zcm+=dzcm){
	for(int i=0; i<=numSteps_z; i++){
		//for loop integrates over all z's from 0 to Distance D.
		zcm=i*dzcm;
		zdistance=zcm/100.;//m	
		if(zcm==0.)
			powerStep=0.;
		else{
			//powerStep=fabs(log10(tauEnergyInitial(ptau, Dcm, zcm)-log10(ptaui)));
			powerStep=fabs(log10(TauEnergyInitial(ptau_final, Dcm, zcm)/ptaui));
		}
		ptaui=TauEnergyInitial(ptau_final, Distance, zdistance);
		ptaui_GeV=ptaui/1.E9;//eV
		// Solving for ptaui acts as the delta function d(f(ptaui, z')-ptau_final);
		//where f(ptaui, z') is from 
		//equation found using Equation 13 Case (III),
		//from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno.
		dptaui=powerStep*log(10.)*ptaui;//eV 
		y=1.-ptaui/pnu;
		
		yweight=primary1->Getyweight(pnu,y,nu_nubar,currentint);
		//(-1/pnu)=dy/dptaui.

		prob_at_zcm=d_dzPsurvNu(primary1,settings1,pnu, nu_nubar, currentint, zdistance)*probabilityTauSurv(ptaui, ptau_final)*yweight*(-1/pnu)*dzcm*dptaui;

		TauWeight+=prob_at_zcm;//sums all the steps together.
	}

	return TauWeight;//the total tau weight for a given Ev, Etau_f, and D.
} //GetTauWeight




// icemc new version (org version)
double Secondaries::GetTauWeight(Primaries *primary1, Settings *settings1,IceModel*antarctica1, double pnu, int nu_nubar, int currentint, double Etau_final, const Position posnu, const Position earth_in,
				 int& crust_entered, // 1 or 0
				 int& mantle_entered, // 1 or 0
				 int& core_entered, double *myxarray, double *myEarray, 
				 double *myyweightarray, double *mytausurvarray, double& tauchord, double *avgdensityarray, 
				 double *densityarray,int inu,double& TauWeight, double& weight_prob){
	
  // Settings *settings1=new Settings();
  //Primaries *primary1=new Primaries();
  Vector chord3;
  Vector nchord;
  
  int  N=1E3;
 

 //Find the chord, its length and its unit vector.
    chord3 = posnu - earth_in;
    double Distance=chord3.Mag();
    nchord = chord3 / Distance;
    tauchord=Distance;
    
   
   double Etaui_GeV,dEtau,Etau_finalGeV;
   double y, yweight;//d_dzPnu_surv;
   double yweight1=0;
   double  zdistance; //m
   
   zdistance=0.;	
   // double TauWeight=0;//total tau survival probability.
   double prob_at_z; //this will be used to get weight1
   double prob_at_z1;//this will be used for weight_prob
   double pnu_GeV;
   double Prob_Nu_Surv;
   double l;
   double tau_surv;
   double sigma = 0;
   double len_int_kgm2 =0;

   primary1->GetSigma(pnu,sigma,len_int_kgm2,settings1,nu_nubar,currentint);
	
	
   l = len_int_kgm2;
   double step=Tools::dMin(l/densities[1]/10,500.); //how big is the step size
  
   double dz=Distance/step;
   double avgdensity =0;//initilize average density.
   double avgdensity1=0;
   prob_at_z = 0;
   Position posnunow;
   double lat;
   double lon;
   cout <<"inu is "<<inu<<"\n";
   for(zdistance = 0; zdistance<Distance; zdistance +=step){
      //for loop integrates over all z's from 0 to Distance D.
     int i =(int)zdistance/step;
     double z = zdistance/step;
     Vector nchord1 =  zdistance*nchord;
     posnunow = earth_in + nchord1; //vector pointing to the step we are currently on.
     lat = posnunow.Lat();
     lon = posnunow.Lon();
     
     double geoid=antarctica1->Geoid(lat);
     double altitude=posnunow.Mag()-geoid; // what is the altitude of the point. 
     double surface_elevation = antarctica1->SurfaceAboveGeoid(lon,lat);
     double abovesurface=0;  
     
     if (zdistance ==0){
       avgdensity = antarctica1->GetDensity1(altitude,posnunow,posnu, crust_entered, mantle_entered, core_entered,abovesurface);
       if (abovesurface==1)
	 cout <<"altitude (tau) is "<<altitude<<"\n";
     }
     else{
        avgdensity1 = antarctica1->GetDensity1(altitude, posnunow,posnu, crust_entered, mantle_entered, core_entered,abovesurface);
      
	// if (avgdensity1 <920&&zdistance>0){
	
	// cout <<"////////////////////////////////////////////////// \n";
	// cout <<"inu is "<<inu<<"\n";
	// cout <<"avgdensity1 = "<<avgdensity1<<"\n";
	// cout <<"////////////////////////////////////////////////// \n";
	//  }
       avgdensity = avgdensity*z;
       avgdensity +=avgdensity1;  //add the new density to the old
       avgdensity  = avgdensity/(z+1); //divide by two to get the average
     }
      
     tau_surv =0;
     yweight =0;
     pnu_GeV = pnu/1.E9; //Convert E_nu to GeV.   
     Etau_finalGeV=Etau_final/1E9;//Convert Etau_final to GeV.
       
     Etaui_GeV=TauEnergyInitial(Etau_finalGeV, Distance, zdistance, avgdensity);//use Deltafunction to get initial tau energy from the final.
     if (Etaui_GeV>pnu_GeV){
       prob_at_z =0;
       yweight1=0;
     }
     else{
	Prob_Nu_Surv = avgdensity/l*exp(-zdistance*avgdensity*1./l);
	
	y=1.-Etaui_GeV/pnu_GeV;

	yweight=primary1->Getyweight(pnu,y,nu_nubar,currentint);
	
	yweight1 = yweight/.000967912;
	
	tau_surv = probabilityTauSurv(Etaui_GeV,Etau_finalGeV,avgdensity);
	
	double dEtauidEtauf = exp(B1*avgdensity*(Distance-zdistance))*Etaui_GeV/Etau_finalGeV;
	//cout <<"detauideTauf is "<<dEtauidEtauf<<"\n";
        dEtau = log(10)*Etaui_GeV*.1; //changing dEtau to log10 steps.
	
	//prob_at_z=Prob_Nu_Surv*tau_surv*yweight1*(1./pnu_GeV)*dz*dEtau; //old way
	prob_at_z=Prob_Nu_Surv*tau_surv*yweight1*dEtauidEtauf*dz;

	if (prob_at_z>1)
	  cout <<"detauideTauf is "<<dEtauidEtauf<<"\n";
	if(prob_at_z !=prob_at_z)
	  prob_at_z =0;
	
	if(Etaui_GeV < Etau_finalGeV){//THIS SHOULD NEVER HAPPEN. large z limit seems to make it do funky things.
	  
	  cout <<"NOT AGAIN! SOMETHING IS WRONG \n";
	  prob_at_z = 0;
	  }

	if (tau_surv == 0){
	  prob_at_z = 0;
	}
     }//Etaui>pnu
     
       if(i < N&&Etaui_GeV<=pnu_GeV){
	 
       	   myEarray[i]=Etaui_GeV;
	   myxarray[i]=zdistance;      
	   myyweightarray[i]= yweight1;
	   mytausurvarray[i]=tau_surv;
	   avgdensityarray[i]=avgdensity;
	   densityarray[i]=avgdensity1;
       }
       else if (i<N){
       	   myEarray[i]=0;
       	   myxarray[i]=zdistance;      
       	   myyweightarray[i]= 0;
       	   mytausurvarray[i]=0;
          avgdensityarray[i]=avgdensity;
          densityarray[i]=avgdensity1;
        }
       // cout <<"prob_at_z is "<<prob_at_z<<"\n";
	TauWeight+=prob_at_z;//sums all the steps together. units are dP/dE
	

   }//i loop
   
   return TauWeight;
} //GetTauWeight




// icemc new version (modified version)
double Secondaries::GetTauWeight(Primaries *primary1, Settings *settings1,IceModel*antarctica1, double pnu, int nu_nubar, int currentint, double Etau_final, const Position posnu, const Position earth_in,
				 int& crust_entered, // 1 or 0
				 int& mantle_entered, // 1 or 0
				 int& core_entered, 
				 double& TauWeight, double& weight_prob){
	
  // Settings *settings1=new Settings();
  //Primaries *primary1=new Primaries();
  Vector chord3;
  Vector nchord;
  
  int  N=1E3;
 

 //Find the chord, its length and its unit vector.
    chord3 = posnu - earth_in;
    double Distance=chord3.Mag();
    nchord = chord3 / Distance;
    //tauchord=Distance;
    
   
   double Etaui_GeV,dEtau,Etau_finalGeV;
   double y, yweight;//d_dzPnu_surv;
   double yweight1=0;
   double  zdistance; //m
   
   zdistance=0.;	
   // double TauWeight=0;//total tau survival probability.
   double prob_at_z; //this will be used to get weight1
   double prob_at_z1;//this will be used for weight_prob
   double pnu_GeV;
   double Prob_Nu_Surv;
   double l;
   double tau_surv;
   double sigma = 0;
   double len_int_kgm2 =0;

   primary1->GetSigma(pnu,sigma,len_int_kgm2,settings1,nu_nubar,currentint);
	
	
   l = len_int_kgm2;
   double step=Tools::dMin(l/densities[1]/10,500.); //how big is the step size
  
   double dz=Distance/step;
   double avgdensity =0;//initilize average density.
   double avgdensity1=0;
   prob_at_z = 0;
   Position posnunow;
   double lat;
   double lon;
   for(zdistance = 0; zdistance<Distance; zdistance +=step){
      //for loop integrates over all z's from 0 to Distance D.
     int i =(int)zdistance/step;
     double z = zdistance/step;
     Vector nchord1 =  zdistance*nchord;
     posnunow = earth_in + nchord1; //vector pointing to the step we are currently on.
     lat = posnunow.Lat();
     lon = posnunow.Lon();
     
     double geoid=antarctica1->Geoid(lat);
     double altitude=posnunow.Mag()-geoid; // what is the altitude of the point. 
     double surface_elevation = antarctica1->SurfaceAboveGeoid(lon,lat);
     double abovesurface=0;  
     
     if (zdistance ==0){
       avgdensity = antarctica1->GetDensity1(altitude,posnunow,posnu, crust_entered, mantle_entered, core_entered,abovesurface);
       if (abovesurface==1)
	 cout <<"altitude (tau) is "<<altitude<<"\n";
     }
     else{
        avgdensity1 = antarctica1->GetDensity1(altitude, posnunow,posnu, crust_entered, mantle_entered, core_entered,abovesurface);
      
	// if (avgdensity1 <920&&zdistance>0){
	
	// cout <<"////////////////////////////////////////////////// \n";
	// cout <<"inu is "<<inu<<"\n";
	// cout <<"avgdensity1 = "<<avgdensity1<<"\n";
	// cout <<"////////////////////////////////////////////////// \n";
	//  }
       avgdensity = avgdensity*z;
       avgdensity +=avgdensity1;  //add the new density to the old
       avgdensity  = avgdensity/(z+1); //divide by two to get the average
     }
      
     tau_surv =0;
     yweight =0;
     pnu_GeV = pnu/1.E9; //Convert E_nu to GeV.   
     Etau_finalGeV=Etau_final/1E9;//Convert Etau_final to GeV.
       
     Etaui_GeV=TauEnergyInitial(Etau_finalGeV, Distance, zdistance, avgdensity);//use Deltafunction to get initial tau energy from the final.
     if (Etaui_GeV>pnu_GeV){
       prob_at_z =0;
       yweight1=0;
     }
     else{
	Prob_Nu_Surv = avgdensity/l*exp(-zdistance*avgdensity*1./l);
	
	y=1.-Etaui_GeV/pnu_GeV;

	yweight=primary1->Getyweight(pnu,y,nu_nubar,currentint);
	
	yweight1 = yweight/.000967912;
	
	tau_surv = probabilityTauSurv(Etaui_GeV,Etau_finalGeV,avgdensity);
	
	double dEtauidEtauf = exp(B1*avgdensity*(Distance-zdistance))*Etaui_GeV/Etau_finalGeV;
	//cout <<"detauideTauf is "<<dEtauidEtauf<<"\n";
        dEtau = log(10)*Etaui_GeV*.1; //changing dEtau to log10 steps.
	
	//prob_at_z=Prob_Nu_Surv*tau_surv*yweight1*(1./pnu_GeV)*dz*dEtau; //old way
	prob_at_z=Prob_Nu_Surv*tau_surv*yweight1*dEtauidEtauf*dz;

	if (prob_at_z>1)
	  cout <<"detauideTauf is "<<dEtauidEtauf<<"\n";
	if(prob_at_z !=prob_at_z)
	  prob_at_z =0;
	
	if(Etaui_GeV < Etau_finalGeV){//THIS SHOULD NEVER HAPPEN. large z limit seems to make it do funky things.
	  
	  cout <<"NOT AGAIN! SOMETHING IS WRONG \n";
	  prob_at_z = 0;
	  }

	if (tau_surv == 0){
	  prob_at_z = 0;
	}
     }//Etaui>pnu
     
     /*
       if(i < N&&Etaui_GeV<=pnu_GeV){
	 
       	   myEarray[i]=Etaui_GeV;
	   myxarray[i]=zdistance;      
	   myyweightarray[i]= yweight1;
	   mytausurvarray[i]=tau_surv;
	   avgdensityarray[i]=avgdensity;
	   densityarray[i]=avgdensity1;
       }
       else if (i<N){
       	   myEarray[i]=0;
       	   myxarray[i]=zdistance;      
       	   myyweightarray[i]= 0;
       	   mytausurvarray[i]=0;
          avgdensityarray[i]=avgdensity;
          densityarray[i]=avgdensity1;
        }
    */
       // cout <<"prob_at_z is "<<prob_at_z<<"\n";
	TauWeight+=prob_at_z;//sums all the steps together. units are dP/dE
	

   }//i loop
   
   return TauWeight;
} //GetTauWeight




double Secondaries::d_dzPsurvNu(Primaries*primary1, Settings*settings1,double pnu, int nu_nubar, int currentint, double z_distance){//derivative with respect to z of Equation 28. -Connolly Calc 2011.
	//double A=1.;//constant that sets the total probability to unity
	double L=interactionLengthNu(primary1,settings1,pnu,nu_nubar,currentint);//meters
	double Lcm=L*100.;//cm
	double zcm=z_distance*100.;
	double dzPsurvNu=-A*exp(-zcm/Lcm)/Lcm; //units in (1/cm), PsurvNu returns dimensionless value for probability.
	//double dzPsurvNu=-A*exp(z_distance/L)/L; //units 1/m
	return dzPsurvNu; // cm^-1
}
double Secondaries::probabilityTauSurv(double ptaui, double ptau_final){
	//from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
	//Equation 16  & Equation 30.
 
/*	double B0,B1,E0;//parameterization using a logarithmic dependence on energy 
	//for B, the tau elecromagnetic energy loss parameter. 
	double p;//Density of Standard Rock. g/cm^3
	double mT;//mass of the Tau in Gev
	double cT;//Tau Decay length in cm
*/	
	double probSurv;

	probSurv=exp((mT*B1)/(cT*p*pow(B0,2))*(1/ptau_final*(1+log(ptau_final/E0))-1/ptaui*(1+log(ptaui/E0))))*exp(-mT/(cT*B0*p)*(1/ptau_final-1/ptaui));//Equation 30
	return probSurv;
}



double Secondaries::probabilityTauSurv(double ptaui, double ptau_final,double density){
	//from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
	//Equation 16  & Equation 30.
 
/*	double B0,B1,E0;//parameterization using a logarithmic dependence on energy 
	//for B, the tau elecromagnetic energy loss parameter. 
	ouble mT;//mass of the Tau in Gev
	double cT;//Tau Decay length in cm
*/	double p = density;
        double probSurv1;
	double mT = 1.77684; //mass of Tau in GeV
	double cT = 0.00008693; // m
	double A = mT*B1/(cT*p*pow(B0,2));//GeV
	double B = mT/(cT*B0*p);//GeV
	double C = (1/ptau_final)*(1+log(ptau_final/E0));//1/GeV
	double D = (1/ptaui)*(1+log(ptaui/E0));//1/GeV
	double F = (1/ptau_final)-(1/ptaui);//1/GeV
	//	cout<< "A,B,C are "<<A<<","<<B<<","<<C<<"\n";
	//	cout<< "D,F are "<<D<<","<<F<<"\n";
	probSurv1 = exp(A*(C-D)-(B*F));
	
	
	return probSurv1;//Correct
}


double Secondaries::TauEnergyInitial(double ptau_final, double Distance, double z_distance){
	//equation found using Equation 13 Case (III),
	//from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
	double ptaui;
	double ptau_finalGeV=ptau_final/1.E9;//converts eV to GeV.
	double zprime_cm=(Distance-z_distance)*100.;//(d-z)meters*100cm/1m=(d-z) cm.
	
	//all energies in GeV.
 /*	double B0;
	double B1;
	double E0;
	double p;//Density

	////////////////////// Units////////////////////
	p=2.65;            // g/cm^3
	B0=1.2*pow(10.,-6.); // cm^2/g
	B1=0.16*pow(10.,-6.);// cm^2/g
	E0=pow(10.,10.);     // GeV
 */
	double base=ptau_finalGeV/E0*exp(B0/B1*(1-exp(-B1*p*zprime_cm)));
	double power=1/exp(-B1*p*zprime_cm);
	ptaui=E0*pow(base,power);
	return ptaui;
}



double Secondaries::TauEnergyInitial(double ptau_final, double Distance, double z_distance, double density){
	//equation found using Equation 13 Case (III),
	//from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno. 
        double p = density;
        double ptaui;
	double zprime=(Distance-z_distance);
	
	//double base=ptau_finalGeV/E0*exp(B0/B1*(1-exp(-B1*p*zprime)));
	
	
	//double power=exp(B1*p*zprime);
	
	//ptaui=E0*pow(base,power); // Correct. Returns E tau initial in GeV.
	ptaui = E0*exp((log(ptau_final/E0)+B0/B1*(1-exp(-B1*p*zprime)))*exp(B1*p*zprime));
	
	return ptaui;
}



double Secondaries::interactionLengthNu(Primaries*primary1,Settings*settings1,double pnu,int nu_nubar,int currentint){
	//Equation 28 from Connolly Calc 2011.
	double Lcm;
	//double pnuGeV=pnu/1.E9;//Convert eV to GeV.
	double CrossSection;

	double sigma=0.;          //should this have to be declared here?
	double len_int_kgm2=0.;// & should this have to be declared here?
	
	CrossSection= primary1->GetSigma(pnu,sigma,len_int_kgm2,settings1,nu_nubar,currentint);
	//even though they are not usually    ^declared here^?
	//double Mn=1.672622*pow(10.,-24.);// nucleon/ proton mass in grams,also equal to 0.938 GeV.
	//double p=2.65; //density of standard rock g/cm^3

	Lcm=Mn/(CrossSection*1.E4*p);//Lcm is in cm. cross section is in m, must convert to cm.
	double L=Lcm/100.;//convert cm to m.
	return L;//meters.
}
