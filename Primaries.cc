
#include "TRandom3.h"
#include "Constants.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include "Primaries.h"
#include "Settings.h"
#include "counting.hh"
#include "Spectra.h"

#include "Vector.h"
#include "Position.h"
#include "EarthModel.h"
#include "IceModel.h"
#include "Detector.h"
#include "Ray.h"
#include "secondaries.hh"
#include "signal.hh"
#include "RaySolver.h"
#include "Report.h"

//#include "vector.hh"
//#include "position.hh"
//#include "earthmodel.hh"
//#include "icemodel.hh"

#include <cmath>



#include "TH2D.h"
#include "TCanvas.h"


ClassImp(Y);
ClassImp(Primaries);
ClassImp(Interaction);

const double Y::miny_low(0.00002);
const double Y::maxy_low(0.001);
const double Y::miny_high(0.001);
const double Y::maxy_high(1.);

/*
const double Interaction::nu_banana_theta_angle(-0.413 * 3.14159/180.);// don't let me use RADDEG which is annoying 
const double Interaction::altitude_nu_banana(-400.);//Depth of interaction of banana neutrino
const double Interaction::lat_nu_banana(0.); 
const double Interaction::lon_nu_banana(0.);
const double Interaction::banana_slopey(0.);//Turn slopyness off for banana plots (SLOPEY)
const double Interaction::nu_banana_phi_angle(0. * 3.14159/180.); 
const double Interaction::phi_nu_banana(3.14159/4); //Location in phi
const double Interaction::banana_observation_distance(600000.);//How far from the surface above the interaction are we when we measure the voltages? (meters) Note: Should be at least 100000 for best results.
const double Interaction::theta_nu_banana(170.*3.14159/180.);//Location of banana neutrino in theta
const double Interaction::pnu_banana(2.00E19);
const double Interaction::banana_y(0.2);//Elasticity.  0.2 is an average number.
const double Interaction::banana_signal_fluct(0.);//Turn off noise for banana plots (settings1->SIGNAL_FLUCT) (shouldn't matter)
const double Interaction::banana_sigma(0.);//NSIGMA in the case of a banana plot

*/








Primaries::Primaries(){//constructor

  // This is for parametrizations in Connolly et al. 2011  
  //in the form of [i][j] where i is neutrino type(nu_nubar) and j is current type, "nc" vs "cc".
  //[nu_nubar][currentint]
  //[0=nu, 1=nubar][0=neutral current, 1=charged current]
  //[0][0]->[nu][neutral current]
  //[0][1]->[nu][charged current]
  //[1][0]->[nubar][neutral current]
  //[1][1]->[nubar][charged current]
  
  //[nu][neutral current]
  c0[0][0]=-1.826;
  c1[0][0]=-17.31;
  c2[0][0]=-6.448; 
  c3[0][0]=1.431;
  c4[0][0]=-18.61;
  
  //[nu][charged current]
  c0[0][1]=-1.826;
  c1[0][1]=-17.31;
  c2[0][1]=-6.406; 
  c3[0][1]=1.431;
  c4[0][1]=-17.91;
  
  //[nubar][neutral current]	
  c0[1][0]=-1.033;
  c1[1][0]=-15.95;
  c2[1][0]= -7.296; 
  c3[1][0]=1.569;
  c4[1][0]=-18.30;
  
  //[nubar][charged current]
  c0[1][1]=-1.033;
  c1[1][1]=-15.95;
  c2[1][1]=-7.247; 
  c3[1][1]=1.569;
  c4[1][1]=-17.72;
  

  ////////////////////////////
  // upper bound
  //[nu][neutral current]
  c0_upper[0][0]= -1.456;
  c1_upper[0][0]= 32.23;
  c2_upper[0][0]= -32.32; 
  c3_upper[0][0]= 5.881;
  c4_upper[0][0]= -49.41;
  
  //[nu][charged current]
  c0_upper[0][1]= -1.456;
  c1_upper[0][1]= 33.47;
  c2_upper[0][1]= -33.02; 
  c3_upper[0][1]= 6.026;
  c4_upper[0][1]= -49.41;
  
  //[nubar][neutral current]	
  c0_upper[1][0]= -2.945;
  c1_upper[1][0]= 143.2;
  c2_upper[1][0]= -76.70; 
  c3_upper[1][0]= 11.75;
  c4_upper[1][0]= -142.8;
  
  //[nubar][charged current]
  c0_upper[1][1]= -2.945;
  c1_upper[1][1]= 144.5;
  c2_upper[1][1]= -77.44; 
  c3_upper[1][1]= 11.90;
  c4_upper[1][1]= -142.8;


  ////////////////////////////
  // lower bound
  //[nu][neutral current]
  c0_lower[0][0]= -15.35;
  c1_lower[0][0]= 16.16;
  c2_lower[0][0]= 37.71;
  c3_lower[0][0]= -8.801;
  c4_lower[0][0]= -253.1;
  
  //[nu][charged current]
  c0_lower[0][1]= -15.35;
  c1_lower[0][1]= 13.86;
  c2_lower[0][1]= 39.84; 
  c3_lower[0][1]= -9.205;
  c4_lower[0][1]= -253.1;
  
  //[nubar][neutral current]	
  c0_lower[1][0]= -13.08;
  c1_lower[1][0]= 15.17;
  c2_lower[1][0]= 31.19; 
  c3_lower[1][0]= -7.757;
  c4_lower[1][0]= -216.1;
  
  //[nubar][charged current]
  c0_lower[1][1]= -13.08;
  c1_lower[1][1]= 12.48;
  c2_lower[1][1]= 33.52; 
  c3_lower[1][1]= -8.191;
  c4_lower[1][1]= -216.1;
  //
  //




  char ch[50];
  string stmp;
  string sbase="fsigma";
  string bound_upper="upper";
  string bound_lower="lower";
  for(int i=0; i<=1;i++){ // nu, nubar
    for(int j=0; j<=1; j++){ // nc, cc
      sprintf(ch,"%d%d",i,j);
      stmp=ch;	
      //m_fsigma[i][j]=new TF1((sbase+stmp).c_str(),"pow(10, [1]+[2]*log(x-[0])+[3]*pow(log(x-[0]),2)+[4]/log(x-[0]))", 4., 21.);//check bounds. they're in log10 GeV.
      m_fsigma[i][j]=new TF1((sbase+stmp).c_str(),"pow(10, [1]+[2]*log(x-[0])+[3]*pow(log(x-[0]),2)+[4]/log(x-[0]))", 4., 12.);//check bounds. they're in log10 GeV.
      //x=log10(pnu/GeV).
      m_fsigma[i][j]->SetParameters(c0[i][j], c1[i][j], c2[i][j], c3[i][j], c4[i][j]);

      m_fsigma_upper[i][j]=new TF1((sbase+stmp+bound_upper).c_str(),"pow(10, [1]+[2]*log(x-[0])+[3]*pow(log(x-[0]),2)+[4]/log(x-[0]))", 4., 12.);//check bounds. they're in log10 GeV.
      //x=log10(pnu/GeV).
      m_fsigma_upper[i][j]->SetParameters(c0_upper[i][j], c1_upper[i][j], c2_upper[i][j], c3_upper[i][j], c4_upper[i][j]);

      m_fsigma_lower[i][j]=new TF1((sbase+stmp+bound_lower).c_str(),"pow(10, [1]+[2]*log(x-[0])+[3]*pow(log(x-[0]),2)+[4]/log(x-[0]))", 4., 12.);//check bounds. they're in log10 GeV.
      //x=log10(pnu/GeV).
      m_fsigma_lower[i][j]->SetParameters(c0_lower[i][j], c1_lower[i][j], c2_lower[i][j], c3_lower[i][j], c4_lower[i][j]);

      //"fsigma00"->[nu][neutral current]
      //"fsigma01"->[nu][charged current]
      //"fsigma10"->[nubar][neutral current]
      //"fsigma11"->[nubar][charged current]
    }		
  }




  m_csigma=new TCanvas("m_csigma","m_csigma title",1000, 700);
  m_hsigma=new TH2D("hsigma","title hsigma", 600, 7., 12., 600, -40., -30.);
  
  m_hsigma->SetTitle("log10 (pnu) vs.log10 Cross Section Sigma");
  m_hsigma->GetXaxis()->SetTitle("Log10(Ev/ GeV)");
  m_hsigma->GetYaxis()->SetTitle("log10(Cross Section/ m^2)");
  
  m_hsigma->Draw("scat");
  m_hsigma->SetMarkerStyle(7);
  m_hsigma->SetMarkerSize(3);
  
  // again y distributions from Connolly et al. 2011
  m_myY=new Y();
  
  //Low y/////////////////////////////From Table V. Connolly Calc 2011.
  //A_low[4];//same for any [i]nu_nubar and [j]currentint.
  A_low[0]=0.0;
  A_low[1]=0.0941;
  A_low[2]=4.72;
  A_low[3]=0.456;
  //high y////////////////////////////
  //[0=nu, 1=nubar][0=neutral current, 1=charged current]
  //
  //[nu][NC]
  A0_high[0][0]=-0.005;
  A1_high[0][0]=0.23;
  A2_high[0][0]=3.0;
  A3_high[0][0]=1.7;
  
  //[nu][CC]
  A0_high[0][1]=-0.008;
  A1_high[0][1]=0.26;
  A2_high[0][1]=3.0;
  A3_high[0][1]=1.7;
  
  //[nu_bar][NC]
  A0_high[1][0]=-0.005;
  A1_high[1][0]=0.23;
  A2_high[1][0]=3.0;
  A3_high[1][0]=1.7;
  
  //[nu_bar][CC]
  A0_high[1][1]=-0.0026;
  A1_high[1][1]=0.085;
  A2_high[1][1]=4.1;
  A3_high[1][1]=1.7;
  
  b0=2.55;
  b1=-0.0949; //C2_low=b0+b1*epsilon;
  
  ymin_low=2.E-5;
  ymax_low=1.E-3;
  ymin_high=ymax_low;
  ymax_high=1.;
  
  run_old_code=0;//for GetSigma() & Gety() runs of the old code if 1, else runs current code.

  // 0=Reno
  // 1=Connolly et al. 2011
  mine[0]=1.2E15;
  //mine[1]=1.E4;// minimum energy for cross section parametrizations
  mine[1]=1.E13;// minimum energy for cross section parametrizations
  maxe[0]=1.E21;
  maxe[1]=1.E21; // use the same upper limit for reno as for connolly et al.


  // added for air col
  GetAir(col1);





}

Primaries::~Primaries(){//default deconstructor
  m_hsigma->Draw("same");
  m_csigma->Print("sigmaCrossSection.pdf");
  delete m_hsigma;
  delete m_myY;
  for(int i=0; i<=1;i++){ // nu, nubar
    for(int j=0; j<=1; j++){ // nc, cc
      delete m_fsigma[i][j];

      delete m_fsigma_upper[i][j];
      delete m_fsigma_lower[i][j];
    }
  }
    delete m_csigma;


}//deconstructor



// old GetSigma function from icemc
// leave as a reference
int Primaries::GetSigma(double pnu,double& sigma,double &len_int_kgm2,Settings *settings1,int nu_nubar,int currentint){
  // calculate cross section
  if (pnu<mine[settings1->SIGMAPARAM] || pnu>maxe[settings1->SIGMAPARAM]) {
    cout <<  "Need a parameterization for this energy region.\n";
    return 0;
  } //if
  else {
   
    //nu_nubar=1;//default.
    //nu=0, nubar=1


    if(nu_nubar!=0 && nu_nubar!=1){   
      cout<<"nu_nubar is not defined correctly!\n";
      return 0;
    }
    if (currentint!=0 && currentint!=1){//default "cc"
      cout<<"Current is not cc or nc!\n";
      return 0;
    }
    
    if(settings1->SIGMAPARAM==0){ // Reno
      // fit to cross sections calculated by M.H. Reno using the same method as Gandhi et al, but with the CTEQ6-DIS parton distribution functions instead of the CTEQ4-DIS distribution functions
      sigma=(2.501E-39)*pow(pnu/1.E9,0.3076)*settings1->SIGMA_FACTOR; // 10^18 eV - 10^21 eV(use this one for ANITA)
      //sigma=(1.2873E-39)*pow(pnu/1.E9,0.33646)*SIGMA_FACTOR; // 10^17 eV - 10^20 eV (use this one for SalSA)
    }//old code
    else if (settings1->SIGMAPARAM==1) {//Connolly et al.
      double pnuGeV=pnu/1.E9;//Convert eV to GeV.
      double epsilon=log10(pnuGeV);
      sigma=settings1->SIGMA_FACTOR*(m_fsigma[nu_nubar][currentint]->Eval(epsilon))/1.E4;//convert cm to meters. multiply by (1m^2/10^4 cm^2).
      
      if(m_hsigma->GetEntries()<2000){
	m_hsigma->Fill(epsilon, log10(sigma));
      }
    }//else current code
  }//if
  // interaction length in kg/m^2
  
  len_int_kgm2=M_NUCL/sigma; // kg/m^2
  return 1;
} //GetSigma




int Primaries::GetSigma(double pnu,double& sigma,double &len_int_kgm2,Settings *settings1,int nu_nubar,int currentint, double &len_int_kgm2_total){
    double sigma_total;
  // calculate cross section
  if (pnu<mine[settings1->SIGMAPARAM] || pnu>maxe[settings1->SIGMAPARAM]) {
    cout <<  "Need a parameterization for this energy region.\n";
    return 0;
  } //if
  else {
   
    //nu_nubar=1;//default.
    //nu=0, nubar=1


    if(nu_nubar!=0 && nu_nubar!=1){   
      cout<<"nu_nubar is not defined correctly!\n";
      return 0;
    }
    if (currentint!=0 && currentint!=1){//default "cc"
      cout<<"Current is not cc or nc!\n";
      return 0;
    }

    
    if(settings1->SIGMAPARAM==0){ // Reno
      // fit to cross sections calculated by M.H. Reno using the same method as Gandhi et al, but with the CTEQ6-DIS parton distribution functions instead of the CTEQ4-DIS distribution functions
      sigma=(2.501E-39)*pow(pnu/1.E9,0.3076)*settings1->SIGMA_FACTOR; // 10^18 eV - 10^21 eV(use this one for ANITA)
      //sigma=(1.2873E-39)*pow(pnu/1.E9,0.33646)*SIGMA_FACTOR; // 10^17 eV - 10^20 eV (use this one for SalSA)

      sigma_total = sigma;

    }//old code
    else if (settings1->SIGMAPARAM==1) {//Connolly et al.
      double pnuGeV=pnu/1.E9;//Convert eV to GeV.
      double epsilon=log10(pnuGeV);


      if ( settings1->SIGMA_SELECT == 0 ) {// use mean value

          sigma=settings1->SIGMA_FACTOR*(m_fsigma[nu_nubar][currentint]->Eval(epsilon))/1.E4;//convert cm to meters. multiply by (1m^2/10^4 cm^2).
          sigma_total = (settings1->SIGMA_FACTOR*(m_fsigma[nu_nubar][0]->Eval(epsilon))/1.E4) + (settings1->SIGMA_FACTOR*(m_fsigma[nu_nubar][1]->Eval(epsilon))/1.E4);
      }
      else if ( settings1->SIGMA_SELECT == 1 ) {// use upper bound

          sigma=settings1->SIGMA_FACTOR*(m_fsigma_upper[nu_nubar][currentint]->Eval(epsilon))/1.E4;//convert cm to meters. multiply by (1m^2/10^4 cm^2).
          sigma_total = (settings1->SIGMA_FACTOR*(m_fsigma_upper[nu_nubar][0]->Eval(epsilon))/1.E4) + (settings1->SIGMA_FACTOR*(m_fsigma_upper[nu_nubar][1]->Eval(epsilon))/1.E4);
      }
      else if ( settings1->SIGMA_SELECT == 2 ) {// use lower bound

          sigma=settings1->SIGMA_FACTOR*(m_fsigma_lower[nu_nubar][currentint]->Eval(epsilon))/1.E4;//convert cm to meters. multiply by (1m^2/10^4 cm^2).
          sigma_total = (settings1->SIGMA_FACTOR*(m_fsigma_lower[nu_nubar][0]->Eval(epsilon))/1.E4) + (settings1->SIGMA_FACTOR*(m_fsigma_lower[nu_nubar][1]->Eval(epsilon))/1.E4);
      }
      else {// if other values are chosen, just use mean value

          sigma=settings1->SIGMA_FACTOR*(m_fsigma[nu_nubar][currentint]->Eval(epsilon))/1.E4;//convert cm to meters. multiply by (1m^2/10^4 cm^2).
          sigma_total = (settings1->SIGMA_FACTOR*(m_fsigma[nu_nubar][0]->Eval(epsilon))/1.E4) + (settings1->SIGMA_FACTOR*(m_fsigma[nu_nubar][1]->Eval(epsilon))/1.E4);
      }


      
      if(m_hsigma->GetEntries()<2000){
	m_hsigma->Fill(epsilon, log10(sigma));
      }
    }//else current code
  }//if
  // interaction length in kg/m^2
  
  len_int_kgm2=M_NUCL/sigma; // kg/m^2
  len_int_kgm2_total=M_NUCL/sigma_total; // kg/m^2
  return 1;
} //GetSigma




Vector Primaries::GetAnyDirection() {
  Vector output;
  double rndlist[2];
  gRandom->RndmArray(2,rndlist);
  
  costheta_nutraject=2*rndlist[0]-1;

 
  // pick a neutrino azimuthal angle
  phi_nutraject=2*PI*rndlist[1];
  
  // check that these give the right result
  double thetanu=acos(costheta_nutraject);
  
  double sinthetanu=sin(thetanu);
  
  // find direction vector of neutrino
  // **** are cosine and sine flipped?
  output.SetX(sinthetanu*cos(phi_nutraject));
  output.SetY(sinthetanu*sin(phi_nutraject));
  output.SetZ(costheta_nutraject);

  return output;


}



Vector Primaries::GetThatDirection( double theta, double d_theta ) {
  Vector output;
  double rndlist[2];
  gRandom->RndmArray(2,rndlist);
  
  costheta_nutraject=2*rndlist[0]-1;    // from -1 to 1
  costheta_nutraject= (costheta_nutraject * d_theta) + theta;
  double thetanu=costheta_nutraject;

  costheta_nutraject= cos( costheta_nutraject );


 
  // pick a neutrino azimuthal angle
  phi_nutraject=2*PI*rndlist[1];
  
  // check that these give the right result
  
  double sinthetanu=sin(thetanu);
  
  // find direction vector of neutrino
  // **** are cosine and sine flipped?
  output.SetX(sinthetanu*cos(phi_nutraject));
  output.SetY(sinthetanu*sin(phi_nutraject));
  output.SetZ(costheta_nutraject);

  return output;


}




// The interaction
double Primaries::Gety(Settings *settings1,double pnu,int nu_nubar,int currentint) {
  // THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM 
  //  Ghandhi,Reno,Quigg,Sarcevic  hep-ph/9512364
  //  (the curves are not in their later article.)
  //  There is also a slow energy dependence.
   //cout << "I'm here.\n";
  if(settings1->YPARAM==0){
  	double rnd;
  	double x = 0;
  	const double R1=0.36787944;  // 1/e
	const double R2=0.63212056;  // 1-r1
  
  	// generate according to Ghandi fig. 6 
  	// adjust exponent until looks like the curve
  	//  and has right mean.
  	//  (Note this is not the fcn, but the inverse of the integral...)
  


  	rnd = gRandom->Rndm(1); // (0,1)
  	//  cout << "R1, R2, rnd are " << R1 << " " << R2 << " " << rnd << "\n";

	x=pow(-log(R1+rnd*R2),2.5); 


  	return x;   
  }//old Gety

  else if (settings1->YPARAM==1) { //use prescription in Connolly et al.2011
  	//nu_nubar=0;
	double pnuGeV=pnu/1.E9;
        double epsilon=log10(pnuGeV);
	//double elast_y=m_myY->pickY(nu_nubar,currentint,pnuGeV);
	double elast_y=m_myY->pickY(nu_nubar,currentint,epsilon);
	return elast_y;   
  }//current Gety

} //Gety

double Primaries::Getyweight(double pnu, double y, int nu_nubar, int currentint){
	//from Connolly Calc 2011, Equations 9, 10, 11, 16, and 17.
	double dy=0.;//default
	//Ev, cc or nc, nu or nubar.
	
	double C0_highbar, C0_lowbar,C0_high, C0_low;//these C0's are normalization factors.
	double dNdy=0.;//default
	double U, W, B, T;//are added in to help with readability of equations.
	double C1_low, C2_low, C1_high, C1_highbar;
	double weighty;
	double epsilon=log10(pnu/1.E9);

	C2_low=b0+b1*epsilon;//Eq(17)
	C1_low=A_low[0]+A_low[1]*(-exp(-(epsilon-A_low[2])/A_low[3]));//Eq(16)
	int nu_nubar_tmp=0;
	C1_high=A0_high[nu_nubar_tmp][currentint]+A1_high[nu_nubar_tmp][currentint]*(-exp(-(epsilon-A2_high[nu_nubar_tmp][currentint])/A3_high[nu_nubar_tmp][currentint]));

	nu_nubar_tmp=1;
	C1_highbar=A0_high[nu_nubar_tmp][currentint]+A1_high[nu_nubar_tmp][currentint]*(-exp(-(epsilon-A2_high[nu_nubar_tmp][currentint])/A3_high[nu_nubar_tmp][currentint]));//Eq(16)

	if(nu_nubar==0) {
		U=1-1/C2_low;
		W=fabs( (ymax_high-C1_high)/(ymin_high-C1_high));
		B=(pow(ymax_low-C1_low, 1/C2_low)/(ymax_low-C1_high));
		T=B*((pow(ymax_low-C1_low, U)-pow(ymin_low-C1_low, U) )/U)+log(W);
		C0_high=1/T;	
		C0_low=C0_high*(pow(ymax_low-C1_low, 1/C2_low))/(ymax_low-C1_high);
		if(y<ymax_low){//Eq(9)
				dy=0.00002;
				dNdy=C0_low/pow(y-C1_low, 1/C2_low);//Eq(10)
		}
		else if(y>=ymax_low && y<1.){//Eq(9)
				dy=0.001;
				dNdy=C0_highbar/(y-C1_highbar);//Eq(10)
		}
		else{
			dNdy=0.;
      			cout<<"y value is outside of the domain of y.\n";
		}
	}
	else if(nu_nubar==1){
		U=1-1/C2_low;
		W=fabs( (ymax_high-C1_highbar)/(ymin_high-C1_highbar));
		B=(pow(ymax_low-C1_low, 1/C2_low)/(ymax_low-C1_highbar));
		T=B*((pow(ymax_low-C1_low, U)-pow(ymin_low-C1_low, U) )/U)+log(W);
		C0_highbar=1/T;	
		C0_lowbar=C0_highbar*(pow(ymax_low-C1_low, 1/C2_low))/(ymax_low-C1_highbar);
		if(y<ymax_low){
			dy=0.00002;
			dNdy=C0_lowbar/pow(y-C1_low, 1/C2_low);
		}
		else if(y>=ymax_low && y<1.){
			dy=0.001;
			dNdy=C0_highbar/(y-C1_highbar);
		}
		else{
			dNdy=0;
      cout<<"y value is outside of the domain of y.\n";
		}
	}
	else{
	    cout<<"Nu_nubar is not defined!\n";
	}		
	weighty=dNdy*dy;
	return weighty;
}//Getyweight

string Primaries::GetCurrent() {
  // choose CC or NC
  //  get from ratios in Ghandi etal paper
  // updated for the CTEQ6-DIS parton distribution functions
  string current;
  double rnd=gRandom->Rndm();
  if (rnd<=0.6865254) // 10^18 eV - 10^21 eV (use this one for ANITA)
//if (rnd<=0.6893498) // 10^17 eV - 10^20 eV (use this one for SalSA)
    current="cc";
  else
    current="nc";  
  return current;
} //GetCurrent

string Primaries::GetNuFlavor() {
    // pick a neutrino type, flavor ratio 1:1:1
  
  string nuflavor;

  double rnd=gRandom->Rndm();

  if (rnd<=(1./3.)) {  
    nuflavor="nue";
  } //if
  else if(rnd<=(2./3.)) { 
    nuflavor="numu";
  } //else if
  else if(rnd<=(1.)) { 
    nuflavor="nutau";
  } //else if
  else
    cout << "unable to pick nu flavor\n";
  return nuflavor;
} //GetNuFlavor




string Primaries::GetNuFlavor(Settings *settings1) {
    // pick a neutrino type, flavor ratio 1:1:1
  
  string nuflavor;

  double rnd=gRandom->Rndm();

  if ( settings1->SELECT_FLAVOR == 1 ) {
      nuflavor="nue";
  }
  else if ( settings1->SELECT_FLAVOR == 2 ) {
      nuflavor="numu";
  }
  else if ( settings1->SELECT_FLAVOR == 3 ) {
      nuflavor="nutau";
  }

  else { // not specifically select the flavor

      if (rnd<=(1./3.)) {  
        nuflavor="nue";
      } //if
      else if(rnd<=(2./3.)) { 
        nuflavor="numu";
      } //else if
      else if(rnd<=(1.)) { 
        nuflavor="nutau";
      } //else if
      else
        cout << "unable to pick nu flavor\n";
  }

  return nuflavor;
} //GetNuFlavor







int  Primaries::GetNuNuBar( string nuflavor ) {
     // depending on nu flavor, choose nu or nubar
     //
     // based on arXiv:1108.3163, section 3
      
    double rnd = gRandom->Rndm();

    int nu_nubar_out = 0; // start with nu
	
    if (nuflavor=="nue") {

        if ( rnd <= 0.78 )
            nu_nubar_out = 0;
        else 
            nu_nubar_out = 1;
    }
    else { // for both numu and nutau
            
        if ( rnd <= 0.61 )
            nu_nubar_out = 0;
        else 
            nu_nubar_out = 1;
    }

    return nu_nubar_out;
}






// copied from icemc.cc
double Primaries::GetThisAirColumn(Settings* settings1, Position r_in,Vector nnu,Position posnu, double& cosalpha,double& mytheta,
			double& cosbeta0,double& mybeta) {  
    double myair=0; // this is the output
    // it is the column of air in kg/m^2
    cosalpha=(r_in * nnu) / r_in.Mag(); // cosangle that the neutrino enters the earth wrt surface normal at its entrry point
    mytheta=(double)(acos(cosalpha)*DEGRAD)-90.; // turn this into an angle
    
    //------------------added on Dec 8------------------------
    if (settings1->ATMOSPHERE) {
	int index11=int(mytheta*10.); // which index this theta corresponds to
	int index12=index11+1;
	
	// find column of air at this theta
	myair=(col1[index11]+(col1[index12]-col1[index11])*(mytheta*10.-double(index11)))*10.;//unit is kg/m^2
    }
    else 
	myair=0.;//don't include effect of atmosphere
    
    //cout<<"mytheta="<<mytheta<<"; myair="<<myair<<endl;
    //------------------added on Dec 8------------------------
    
    
    cosbeta0= (posnu * nnu) / posnu.Mag(); // cos angle of neutrino wrt person standing over the interaction point
    mybeta=(double)(acos(cosbeta0)*DEGRAD)-90.; // turn that into a theta
    
    return myair;
    
}



// copied from icemc.cc
void Primaries::GetAir(double *col1) {
    double nothing;
    ifstream air1("data/atmosphere.dat"); // length of chord in air vs. theta (deg)
    //where theta is respect to "up"   
    // binned in 0.1 degrees
    for(int iii=0;iii<900;iii++) {
	air1>>nothing>>col1[iii];
    } // read in chord lengths
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Interaction::Interaction() {
    Initialize ();
    //default constructor
}

Interaction::~Interaction() {

    vmmhz1m.clear();
    vmmhz1m_em.clear();
    d_theta_em.clear();
    d_theta_had.clear();

}


void Interaction::Initialize() {
    // settings for GetSignal
    taudecay = "test_taudecay";
}

void Interaction::clear_useless(Settings *settings1){
 
  if(settings1->DATA_SAVE_MODE>0){
    
    shower_depth_m.clear();
    shower_Q_profile.clear();
    EM_shower_depth_m.clear();
    EM_shower_Q_profile.clear();
    HAD_shower_depth_m.clear();
    HAD_shower_Q_profile.clear();
  
  }
  
}

// remove ANITA option for ARASIM
//
/*
Interaction::Interaction(string inttype,Primaries *primary1,Settings *settings1,int whichray,Counting *count1) : banana_flavor("numu"), banana_current("nc"),  nu_banana(Position(theta_nu_banana,phi_nu_banana)) {

  noway=0;
  wheredoesitleave_err=0;
  neverseesice=0;

    wheredoesitenterice_err=0;
    toohigh=0;
    toolow=0;

  iceinteraction=0;
  dtryingdirection=0.;
  dnutries=0.;

    if (inttype=="banana") {
      nu_banana = (surface_over_banana_nu+altitude_nu_banana) * nu_banana;

      //Set neutrino direction
      nnu_banana = Vector(nu_banana_theta_angle + PI,nu_banana_phi_angle);
      nnu_banana = nnu_banana.ChangeCoord(nu_banana);
           
      current = banana_current;

      nuflavor = banana_flavor;


    }
    else {
    setNuFlavor(primary1,settings1,whichray,count1);
    setCurrent(primary1);
    //    setnu_nubar(primary1);//same function for inttype "banna" or otherwise.
    }
}

*/

Interaction::Interaction (double pnu, string nuflavor, int nu_nubar, int &n_interactions, IceModel *antarctica, Detector *detector, Settings *settings1, Primaries *primary1, Signal *signal, Secondaries *sec1 ) {

    Initialize ();


    if (settings1->NNU_THIS_THETA==1) {    // set specific theta angle for nnu
        nnu = primary1->GetThatDirection(settings1->NNU_THETA, settings1->NNU_D_THETA);
    }
    else { // nnu angle random
        nnu = primary1->GetAnyDirection();
    }
        
    cone_axis = nnu;



    setCurrent(primary1);   // set current of interaction (cc or nc) ! (this should be change if this is secondary interaction. if first interaction was cc, then there is no secondary cc)
    //cout<<"currentint : "<<currentint<<endl;

    // pick posnu (position where neutrino interact with ice
    // also they will calculate r_in (position nu enter the earth), r_enterice (position nu enter the ice), nuexitice (position nu exit the ice)
    //

    double L0 = 0.;


    if (settings1->CALPULSER_ON == 1)  // calpulser1 Hpol
//    {
//        PickExact(antarctica, detector, settings1, 50.27, 36.35*PI/180., 64.76*PI/180.); 
//        primary1->IsCalpulser = 1;
//    }
    {
        PickExact(antarctica, detector, settings1, 29.5, 2.81*PI/180., -93.6*PI/180.);
        primary1->IsCalpulser = 1;
    }
    else if(settings1->CALPULSER_ON == 2)  // calpulser2 Vpol
    {
        //PickExact(antarctica, detector, settings1, 47.18, 30.15*PI/180., -36.56*PI/180.); 
        PickExact(antarctica, detector, settings1, 47.18, -20.*PI/180., 34.*PI/180.); // calpulser 2 Vpol location
        primary1->IsCalpulser = 2;
    } 
    else if(settings1->CALPULSER_ON == 3)   // calpulser2 Hpol
    {
        //PickExact(antarctica, detector, settings1, 47.18, 30.15*PI/180., -36.56*PI/180.); 
        PickExact(antarctica, detector, settings1, 47.18, -28.*PI/180., 34.*PI/180.); // calpulser 2 Hpol location
        primary1->IsCalpulser = 3;
    } 
    else if(settings1->CALPULSER_ON == 4)   // calpulser2 middle of Vpol & Hpol
    {
        //PickExact(antarctica, detector, settings1, 47.18, 30.15*PI/180., -36.56*PI/180.); 
        PickExact(antarctica, detector, settings1, 47.18, -24.*PI/180., 34.*PI/180.); // calpulser 2 Hpol location
        primary1->IsCalpulser = 4;
    } 
    else if (settings1->CALPULSER_ON == 5)
    {
        // insert code to switch between the two calpulsers here
    }
    else if (settings1->CALPULSER_ON == 0){
        primary1->IsCalpulser = 0;

    if (settings1->INTERACTION_MODE == 0) {    // for pickunbiased. posnu will be selected the sphere around the stations
        //Interaction::PickUnbiased( antarctica );
        
        //pickposnu = 0;
        /*
        //  try PickUnbiased until it founds usable posnu
        while (pickposnu==0) {
            PickUnbiased( antarctica );
        }
        */
        //  try PickUnbiased just once, no matter pickposnu == 0 or 1
        //PickUnbiased( antarctica ); // PickUnbiased no more used

        L0 = PickNear_Sphere (antarctica, detector, settings1);

    }
    else if (settings1->INTERACTION_MODE == 1) {   // for picknear. posnu will be only near by ARA core with cylinderical volume
        
        //Interaction::PickNear (antarctica, detector, settings1);
        PickNear_Cylinder (antarctica, detector, settings1);

    } else if (settings1->INTERACTION_MODE == 2){       
//        PickExact(antarctica, detector, settings1, 400, 0, 0); 
          PickExact(antarctica, detector, settings1, 1000, -PI/4., -PI/3.);
    }


    }
    //cout<<" Finished Pick posnu, r_in, r_enterice, nuexitice!!"<<endl;


    // now set N at posnu
    signal->SetNDepth( antarctica->GetN( posnu ) );
    indexN = signal->N_DEPTH;
    changle = signal->changle;



    //sigma_err = primary1->GetSigma( pnu, sigma, len_int_kgm2, settings1, nu_nubar, currentint);
    sigma_err = primary1->GetSigma( pnu, sigma, len_int_kgm2, settings1, nu_nubar, currentint, len_int_kgm2_total );
//--------------------------------------------------
//     cout<<"len_int_kgm2 from GetSigma : "<<len_int_kgm2<<endl;
//     cout<<"sigma from GetSigma : "<<sigma<<endl;
//     cout<<" Finished SetSigma!!"<<endl;
//-------------------------------------------------- 

    double tmp; // for useless information

    myair = primary1->GetThisAirColumn( settings1, r_in, nnu, posnu, tmp, tmp, tmp, tmp);
    //cout<<" Finished GetThisAirColumn!!"<<endl;

    //cout<<"test EarthModel, radii[0] : "<<antarctica->radii[0]<<endl;


    // check tau mode
    int taumodes1 = 0.;

    if (settings1->taumodes ==1){
            double xrndm = gRandom->Rndm();
            if(xrndm <.5)
                    taumodes1 = 0;
            else
                    taumodes1= 1;
    }
    int trip =0;
    double lptau =0;
    double ptauf = 0.;
    
    if (settings1->taumodes ==1){
            ptauf =0;
            double Emin = 1E19*exp(-7.5);
            // while (ptauf <  1E19*exp(-7.5)){
            double ptaurndm= gRandom->Rndm();
            
            lptau =log10(Emin)+ptaurndm*(log10(pnu)-log10(Emin));
            ptauf = pow(10,lptau);
            
            if (ptauf < 1E16){
                    trip =1;
            }
    }


   
    elast_y = primary1->Gety(settings1, pnu, nu_nubar, currentint);  // set inelasticity
    //cout<<"set inelasticity : "<<elast_y<<endl;

   
    //sec1->GetEMFrac( settings1, nuflavor, current, taudecay, elast_y, pnu, emfrac, hadfrac, n_interactions);   // set em, had frac values.
    //cout<<"set emfrac : "<<emfrac<<" hadfrac : "<<hadfrac<<endl;
    sec1->GetEMFrac(settings1, nuflavor, current, taudecay, elast_y, pnu, emfrac, hadfrac, n_interactions, taumodes1, ptauf); // set em, had frac values



    // below this will be only done when we have chosen the usable posnu
    if ( pickposnu == 1 ) {



        //if (settings1->GETCHORD_MODE==0) {
        if (settings1->INTERACTION_MODE==1) {
            antarctica->Getchord(len_int_kgm2_total, r_in, posnu, 0, chord, weight, nearthlayers, myair, total_kgm2, crust_entered, mantle_entered, core_entered );
        }
        //else if (settings1->GETCHORD_MODE==1 && settings1->INTERACTION_MODE==0) {
        else if (settings1->INTERACTION_MODE==0) {
          antarctica->Getchord(primary1, settings1, antarctica, sec1, len_int_kgm2_total, r_in, r_enterice, nuexitice, posnu, 0, chord, probability, weight, nearthlayers, myair, total_kgm2, crust_entered, mantle_entered, core_entered, nuflavor, pnu, ptauf, nu_nubar, currentint, taumodes1 , L0);
        }

        /*
        else if (settings1->GETCHORD_MODE==1) {
            //antarctica->Getchord(len_int_kgm2_total, r_in, posnu, 0, chord, weight, nearthlayers, myair, total_kgm2, crust_entered, mantle_entered, core_entered );
            antarctica->Getchord(primary1, settings1, antarctica, sec1, len_int_kgm2_total, r_in, r_enterice, nuexitice, posnu, 0, chord, probability, weight, nearthlayers, myair, total_kgm2, crust_entered, mantle_entered, core_entered, nuflavor, pnu, ptauf, nu_nubar, currentint, taumodes1 );
        }
        */

        //cout<<" Finished Getchord!!"<<endl;




           if (settings1->SIMULATION_MODE == 0) { // freq domain simulation (old mode)

               // set vmmhz1m (which is generally used for all detector antennas)
               // vmmhz1m is calculated at 1m, cherenkov angle
               //


               for (int i=0; i<detector->GetFreqBin(); i++) {   // for detector freq bin numbers


                d_theta_em.push_back(0); // prepare d_theta_em and d_theta_had for GetSpread
                d_theta_had.push_back(0);
                vmmhz1m.push_back(0);
                vmmhz1m_em.push_back(0);


                   signal->GetSpread(pnu, emfrac, hadfrac, detector->GetFreq(i), d_theta_em[i], d_theta_had[i]);   // get max spread angle and save at d_theta_em[i] and d_theta_had[i]
                   //cout<<"Freq : "<<detector->GetFreq(i)<<endl;
                   //cout<<"GetSpread, theta_em : "<<d_theta_em[i]<<" theta_had : "<<d_theta_had[i]<<endl;

                   //vmmhz1m[i] = signal->GetVmMHz1m( pnu, detector->GetFreq(i) );   // get VmMHz at 1m at cherenkov angle at GetFreq(i)


                   if (primary1->IsCalpulser == 0){

                       if (settings1->AVZ_NORM_FACTOR_MODE == 0) { // use previous normalization factor ( with sqrt(2) )
                           vmmhz1m[i] = signal->GetVmMHz1m( pnu, detector->GetFreq(i) );   // get VmMHz at 1m at cherenkov angle at GetFreq(i)
                       }
                    
                       else if (settings1->AVZ_NORM_FACTOR_MODE == 1) { // use new normalization factor ( without sqrt(2) )
                           vmmhz1m[i] = signal->GetVmMHz1m( pnu, detector->GetFreq(i) ) / sqrt(2.);   // get VmMHz at 1m at cherenkov angle at GetFreq(i), cancel sqrt(2) factor
                       }

                   } else {
                       if (i <56){
                       vmmhz1m[i] = signal->GetVmMHz1mCalPulser( i );
                       } else {
                       vmmhz1m[i] =   signal->GetVmMHz1mCalPulser( 55 );
                       }
                       
                   }


                   
    //            std::cout << "Frequency bin:vmmhz1m:: " << i << " : " << vmmhz1m[i] << " : " << detector->GetFreqBin() << " : " << detector->GetFreq(i) << std::endl;
               }    // end detector freq bin numbers loop

           }// if SIMULATION_MODE = 0 (freq domain old method)


           else if (settings1->SIMULATION_MODE == 1) { // time domain simulation (new mode)

               // here, we just get the shower profile for EM, HAD showers 
               //
               // after we get raytrace solutions, and obtain the view angle, we calculate the signal at 1m and propagate them

               //cout<<"0";

               // only EM shower
               if ( settings1->SHOWER_MODE == 0 ) {
      
                   if ( emfrac > 1.e-10 ) {

                       //signal->GetShowerProfile( pnu*emfrac, 0, settings1->SHOWER_STEP, settings1->SHOWER_PARAM_MODEL, EM_shower_depth_m, EM_shower_Q_profile, EM_LQ );
                       signal->GetShowerProfile( pnu*emfrac, 0, settings1->SHOWER_STEP, settings1->SHOWER_PARAM_MODEL, shower_depth_m, shower_Q_profile, LQ );
                   }
                   else LQ = 0;

                   primary_shower = 0;
                   //cout<<"LQ : "<<LQ<<endl;
               }

               // only HAD shower
               else if ( settings1->SHOWER_MODE == 1 ) {

                   if ( hadfrac > 1.e-10 ) {
                       //signal->GetShowerProfile( pnu*hadfrac, 1, settings1->SHOWER_STEP, settings1->SHOWER_PARAM_MODEL, HAD_shower_depth_m, HAD_shower_Q_profile, HAD_LQ );
                       signal->GetShowerProfile( pnu*hadfrac, 1, settings1->SHOWER_STEP, settings1->SHOWER_PARAM_MODEL, shower_depth_m, shower_Q_profile, LQ );
                   }
                   else LQ = 0;


                   primary_shower = 1;
                   //cout<<"LQ : "<<LQ<<endl;
               }

               // use both EM and HAD 
               else if ( settings1->SHOWER_MODE == 2 ) {

                   if ( emfrac > 1.e-10 ) {

                       signal->GetShowerProfile( pnu*emfrac, 0, settings1->SHOWER_STEP, settings1->SHOWER_PARAM_MODEL, EM_shower_depth_m, EM_shower_Q_profile, EM_LQ );
                   }
                   else EM_LQ = 0;

                   if ( hadfrac > 1.e-10 ) {
               
                       signal->GetShowerProfile( pnu*hadfrac, 1, settings1->SHOWER_STEP, settings1->SHOWER_PARAM_MODEL, HAD_shower_depth_m, HAD_shower_Q_profile, HAD_LQ );
                   }
                   else HAD_LQ = 0;



                   if ( EM_LQ!=EM_LQ ) {// if nan
                       EM_LQ = 0.;
                   }
                   if ( HAD_LQ!=HAD_LQ ) {// if nan
                       HAD_LQ = 0.;
                   }

               
                   // get total LQ
                   LQ = EM_LQ + HAD_LQ;

                   /*
                   cout<<"EM_LQ : "<<EM_LQ<<endl;
                   cout<<"HAD_LQ : "<<HAD_LQ<<endl;
                   cout<<"LQ : "<<LQ<<endl;
                   */


                   //
                   // make shower profiles same length
                   if ( EM_LQ > 0 && HAD_LQ > 0 ) { // if both EM and HAD shower exist


                       double d_depth_EM = EM_shower_depth_m[1] - EM_shower_depth_m[0];
                       double d_depth_HAD = HAD_shower_depth_m[1] - HAD_shower_depth_m[0];

                       int lastbin;

                       while ( (int)EM_shower_Q_profile.size() != (int)HAD_shower_Q_profile.size() ) {

                           if ( (int)EM_shower_Q_profile.size() > (int)HAD_shower_Q_profile.size() ) {

                               lastbin = HAD_shower_Q_profile.size();

                               HAD_shower_Q_profile.push_back( 0 );

                               HAD_shower_depth_m.push_back( HAD_shower_depth_m[lastbin-1] + d_depth_HAD );

                           }

                           else if ( (int)EM_shower_Q_profile.size() < (int)HAD_shower_Q_profile.size() ) {

                               lastbin = EM_shower_Q_profile.size();

                               EM_shower_Q_profile.push_back( 0 );

                               EM_shower_depth_m.push_back( EM_shower_depth_m[lastbin-1] + d_depth_EM );

                           }

                       }

                   }

                   //cout<<"LQ : "<<LQ<<endl;

                   //signal->GetShowerProfile( pnu*emfrac, 0, settings1->SHOWER_STEP, settings1->SHOWER_PARAM_MODEL, EM_shower_depth_m, EM_shower_Q_profile, EM_LQ );
                   //signal->GetShowerProfile( pnu*hadfrac, 1, settings1->SHOWER_STEP, settings1->SHOWER_PARAM_MODEL, HAD_shower_depth_m, HAD_shower_Q_profile, HAD_LQ );
                   //cout<<"EM_LQ : "<<EM_LQ<<", HAD_LQ : "<<HAD_LQ<<endl;
                   //
               

               } // shower_mode = 2



               // only one dominant shower (either EM or HAD)
               else if ( settings1->SHOWER_MODE == 3 ) {


                   // find the dominant shower among EM and HAD
                   if ( emfrac > hadfrac ) {

                       signal->GetShowerProfile( pnu*emfrac, 0, settings1->SHOWER_STEP, settings1->SHOWER_PARAM_MODEL, shower_depth_m, shower_Q_profile, LQ );

                       primary_shower = 0;

                   }
                   else {

                       signal->GetShowerProfile( pnu*hadfrac, 1, settings1->SHOWER_STEP, settings1->SHOWER_PARAM_MODEL, shower_depth_m, shower_Q_profile, LQ );
                   
                       primary_shower = 1;
                   }


               } // shower_mode = 3




           } // tdomain mode


    }// if pickposnu


}


// GetSignal should move to Report class
/*
void Interaction::GetSignal (Settings *settings1, Primaries *primary1, Secondaries *sec1, Signal *signal, Detector *detector, RaySolver *raysolver, IceModel *icemodel, TH1F *hy, int inu, Report *report) {
//Interaction::GetSignal (Primaries *primary1, Settings *settings1, Secondaries *sec1, Signal *signal1, Detector *detector, RaySolver *raysolver, IceModel *icemodel, EvtReport *report1) {

    cout<<"test GetSignal started!"<<endl;


    // clear & initialize Report
    //
    report->Initialize (detector);


       //double elast_y = primary1->Gety(settings1, Evt->pnu, nu_nubar, currentint);  // set inelasticity
       elast_y = primary1->Gety(settings1, pnu, nu_nubar, currentint);  // set inelasticity
       cout<<"set inelasticity : "<<elast_y<<endl;

       sec1->GetEMFrac( settings1, nuflavor, current, taudecay, elast_y, hy, pnu, inu, emfrac, hadfrac, n_interactions);   // set em, had frac values.
       //sec1->GetEMFrac( settings1, Evt, taudecay, elast_y, hy, inu, emfrac, hadfrac, n_interactions);   // set em, had frac values.
       cout<<"set emfrac : "<<emfrac<<" hadfrac : "<<hadfrac<<endl;

       double vmmhz1m_sum;


       if (settings1->SIMULATION_MODE == 0) { // freq domain simulation (old mode)

           // set vmmhz1m (which is generally used for all detector antennas)
           // vmmhz1m is calculated for 1m, cherenkov angle
           //


           for (int i=0; i<detector->GetFreqBin(); i++) {   // for detector freq bin numbers

               signal->GetSpread(pnu, emfrac, hadfrac, detector->GetFreq(i), d_theta_em[i], d_theta_had[i]);   // get max spread angle and save at d_theta_em[i] and d_theta_had[i]
               cout<<"Freq : "<<detector->GetFreq(i)<<endl;
               cout<<"GetSpread, theta_em : "<<d_theta_em[i]<<" theta_had : "<<d_theta_had[i]<<endl;

               vmmhz1m[i] = signal->GetVmMHz1m( pnu, detector->GetFreq(i) );   // get VmMHz at 1m at cherenkov angle at GetFreq(i)
               cout<<"GetVmMHZ1m : "<<vmmhz1m[i]<<endl;



           }    // end detector freq bin numbers loop




           for (int i = 0; i< detector->params.number_of_stations; i++) {

               for (int j=0; j< detector->params.number_of_strings_station; j++) {

                   for (int k=0; k< detector->params.number_of_antennas_string; k++) {

                       // run ray solver, see if solution exist
                       // if not, skip (set something like Sol_No = 0;
                       // if solution exist, calculate view angle and calculate TaperVmMHz

                       if (pickposnu) {    // if posnu is selected inside the antarctic ic:"<<viewangle<<" th_em:"<<d_theta_em[l]<<" th_had:"<<d_theta_had[l]<<" emfrac:"<<emfrac<<" hadfrac:"<<hadfrac<<" vmmhz1m:"<<vmmhz1m[l]<<endl;e
                           
                           raysolver->Solve_Ray(posnu, detector->stations[i].strings[j].antennas[k], icemodel, ray_output);   // solve ray between source and antenna
                           //raysolver->Solve_Ray(Evt, detector->stations[i].strings[j].antennas[k], icemodel, ray_output);   // solve ray between source and antenna
                           
                          //cout<<"solution_toggle : "<<raysolver->solution_toggle<<endl;   
                            
                           ray_sol_cnt = 0;
                           
                           //if (ray_solver_toggle) {  // if there are solution from raysolver
                           if (raysolver->solution_toggle) {  // if there are solution from raysolver
                           //if (raysolver->solution_toggle && ray_solver_toggle) {  // if there are solution from raysolver
                              cout<<"ray_output size : "<<ray_output[0].size()<<endl;
                               
                               while ( ray_sol_cnt < ray_output[0].size() ) {   // for number of soultions (could be 1 or 2)
                                  cout<<"Path length : "<<ray_output[0][ray_sol_cnt]<<"\tView angle : "<<ray_output[1][ray_sol_cnt]<<"\tReceipt angle : "<<ray_output[2][ray_sol_cnt]<<endl;

                                   R1 = detector->stations[i].strings[j].antennas[k];
                                   R2 = posnu;
                                   viewangle = PI/2. - ray_output[1][ray_sol_cnt];

                                   launch_vector = (R1.Cross( R1.Cross(R2) )).Rotate(viewangle, R1.Cross(R2));
                                   viewangle = launch_vector.Angle(nnu);
                                   
                                   // store information to report
                                   report->stations[i].strings[j].antennas[k].view_ang.push_back(viewangle);
                                   report->stations[i].strings[j].antennas[k].rec_ang.push_back(ray_output[2][ray_sol_cnt]);
                                   report->stations[i].strings[j].antennas[k].Dist.push_back(ray_output[0][ray_sol_cnt]);
                                   
                                   report->stations[i].strings[j].antennas[k].vmmhz.resize(ray_sol_cnt+1);


                                   vmmhz1m_sum = 0;

                                   for (int l=0; l<detector->GetFreqBin(); l++) {   // for detector freq bin numbers


                                      cout<<"TaperVmMHz inputs VA:"<<viewangle<<" th_em:"<<d_theta_em[l]<<" th_had:"<<d_theta_had[l]<<" emfrac:"<<emfrac<<" hadfrac:"<<hadfrac<<" vmmhz1m:"<<vmmhz1m[l]<<endl;

                                       vmmhz1m_tmp = vmmhz1m[l];

                                       signal->TaperVmMHz( viewangle, d_theta_em[l], d_theta_had[l], emfrac, hadfrac, vmmhz1m_tmp, vmmhz1m_em[l]);
                                      cout<<"TaperVmMHz (1m at view angle) at "<<l<<"th bin : "<<vmmhz1m_tmp<<endl;

                                      vmmhz1m_tmp = vmmhz1m_tmp / ray_output[0][ray_sol_cnt] * exp(-ray_output[0][ray_sol_cnt]/icemodel->EffectiveAttenuationLength(settings1, posnu, 0));  // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
                                      cout<<"AttenLength : "<<icemodel->EffectiveAttenuationLength(settings1, posnu, 0)<<endl;

                                      vmmhz1m_sum += vmmhz1m_tmp;


                                       report->stations[i].strings[j].antennas[k].vmmhz[ray_sol_cnt].push_back( vmmhz1m_tmp );


                                   }// end for freq bin
                                   
                                   cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"].vmmhz1m["<<ray_sol_cnt<<"][0] : "<<report->stations[i].strings[j].antennas[k].vmmhz[ray_sol_cnt][0]<<endl;
                                       

                                   ray_sol_cnt++;
                               
                               }// end while number of solutions
                           
                           }// end if solution exist

                           else {
                               
                               //cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"].trg = "<<report->stations[i].strings[j].antennas[k].trg[ray_sol_cnt]<<"  No vmmhz1m data!"<<endl;

                           }
                       
                       }// end if posnu selected

                           else {
                               cout<<" No posnu!!!!!! No signals calculated at all!!"<<endl;
                           }


                           report->stations[i].strings[j].antennas[k].ray_sol_cnt = ray_sol_cnt;    // save number of RaySolver solutions
                       

                   }// for number_of_antennas_string

               }// for number_of_strings_station

           }// for number_of_stations




       }// if SIMULATION_MODE = 0 (freq domain old method)


       else if (settings1->SIMULATION_MODE == 1) { // time domain simulation (new mode)
           cout<<"Currently unavailable!!"<<endl;
           // we need to break
       }


}// end GetSignal

*/



int Interaction::PickUnbiased (IceModel *antarctica) {
    

    //Interaction::PickAnyDirection(); // first pick the neutrino direction

  double mincos=cos(antarctica->GetCOASTLINE()*RADDEG);
  double maxcos=cos(0.);
  double minphi=0.;
  double maxphi=2.*PI;
  double thisphi,thiscos,thissin;
  double theta=0.;
  double phi=0.;

  int ilon,ilat;    
  int e_coord,n_coord;
  double vol_thisbin=0.;
  double lon=0.;
  double lat=0.;
  
 
    thisphi=gRandom->Rndm()*(maxphi-minphi)+minphi;
    thiscos=gRandom->Rndm()*(maxcos-mincos)+mincos;
    thissin=sqrt(1.-thiscos*thiscos);
    Position thisr_in;// entrance point
    Position thisr_enterice;
Position thisr_enterice_tmp;
    Position thisnuexitearth;
    Position thisnuexitice;
    Position thisr_exitice;
    noway=0;
    wheredoesitleave_err=0;
    neverseesice=0;
    wheredoesitenterice_err=0;
    toohigh=0;
    toolow=0;

    thisr_in.SetXYZ(antarctica->R_EARTH*thissin*cos(thisphi),antarctica->R_EARTH*thissin*sin(thisphi),antarctica->R_EARTH*thiscos);



    if (thisr_in.Dot(nnu)>0)
      nnu=-1.*nnu;
    // does this intersect any ice
    //cout << "lat, coastline, cos are " << thisr_in.Lat() << " " << antarctica->GetCOASTLINE() << " " << cos(interaction1->nnu.Theta()) << "\n";
    if (thisr_in.Lat()>antarctica->GetCOASTLINE() && cos(nnu.Theta())<0) {
      noway=1;

      pickposnu=0;
      return 0; // there is no way it's going through the ice
    }

    int count1=0;
    int count2=0;

   
    if (Interaction::WhereDoesItLeave(thisr_in,nnu,antarctica,thisnuexitearth)) { // where does it leave Earth
        nuexit = thisnuexitearth;
      // really want to find where it leaves ice
      int err;
      // Does it leave in an ice bin
      if (antarctica->IceThickness(thisnuexitearth) && thisnuexitearth.Lat()<antarctica->GetCOASTLINE()) { // if this is an ice bin in the Antarctic
	//cout << "inu is " << inu << " it's in ice.\n";
	//cout << "this is an ice bin.\n";
	thisnuexitice=thisnuexitearth;
	thisr_exitice=thisnuexitearth;
	if (thisnuexitice.Mag()>antarctica->Surface(thisnuexitice)) { // if the exit point is above the surface
	//--------------------------------------------------
	// if (thisnuexitice.Mag()>Surface(thisnuexitice)) { // if the exit point is above the surface
	//-------------------------------------------------- 
	  if ((thisnuexitice.Mag()-antarctica->Surface(thisnuexitice))/cos(nnu.Theta())>5.E3) { 
	//--------------------------------------------------
	//   if ((thisnuexitice.Mag()-Surface(thisnuexitice))/cos(interaction1->nnu.Theta())>5.E3) { 
	//-------------------------------------------------- 
	    WhereDoesItExitIce(thisnuexitearth,nnu,5.E3, // then back up and find it more precisely
			       thisr_exitice, antarctica);
	    thisnuexitice=(5000.)*nnu;
	    thisnuexitice+=thisr_exitice;
	    count1++;
	  }
	  if ((thisnuexitice.Mag()-antarctica->Surface(thisnuexitice))/cos(nnu.Theta())>5.E2) {
	//--------------------------------------------------
	//   if ((thisnuexitice.Mag()-Surface(thisnuexitice))/cos(interaction1->nnu.Theta())>5.E2) {
	//-------------------------------------------------- 
	    
	    WhereDoesItExitIce(thisnuexitice,nnu,5.E2, // then back up and find it more precisely
			       thisr_exitice, antarctica);
	    thisnuexitice=5.E2*nnu;
	    thisnuexitice+=thisr_exitice;
	    count1++;
	  }
	  if ((thisnuexitice.Mag()-antarctica->Surface(thisnuexitice))/cos(nnu.Theta())>50.) {
	//--------------------------------------------------
	//   if ((thisnuexitice.Mag()-Surface(thisnuexitice))/cos(interaction1->nnu.Theta())>50.) {
	//-------------------------------------------------- 

	    WhereDoesItExitIce(thisnuexitice,nnu,50., // then back up and find it more precisely
			     thisr_exitice, antarctica);
	    count1++;
	  } // end third wheredoesitexit
	  thisnuexitice=thisr_exitice;
	} // if the exit point overshoots
	else
	  thisnuexitice=thisnuexitearth;

	// should also correct for undershooting
        if (count1>10) cout << "count1 is " << count1 << "\n";	  

      } // if it's an Antarctic ice bin


      else { // it leaves a rock bin so back up and find where it leaves ice
	//cout << "inu is " << inu << " it's in rock.\n";
	if (thisr_in.Distance(thisnuexitearth)>5.E4) {
	  count2++;
	  if (WhereDoesItExitIce(thisnuexitearth,nnu,5.E4, // then back up and find it more precisely
				 thisr_exitice, antarctica)) {
	    
	    thisnuexitice=(5.E4)*nnu;
	    thisnuexitice+=thisr_exitice;
	    //cout << "inu is " << inu << " I'm here 1.\n";

	  }
	  else {
	    neverseesice=1;
            pickposnu = 0;
	    return 0;
	  }
	}
	else
	  thisnuexitice=thisnuexitearth;
	//   WhereDoesItExitIce(inu,thisnuexit,interaction1->nnu,5.E4, // then back up and find it more precisely
// 			     thisr_exitice);
// 	  thisnuexit=5.E4*interaction1->nnu;
// 	  thisnuexit+=thisr_exitice;
	if (thisr_in.Distance(thisnuexitice)>5.E3) {

	  
	  if (WhereDoesItExitIce(thisnuexitice,nnu,5.E3, // then back up and find it more precisely
				  thisr_exitice, antarctica)) {
	    count2++;
	    //neverseesice=1;
	    thisnuexitice=5.E3*nnu;
	    thisnuexitice+=thisr_exitice;
	    //cout << "inu is " << inu << " I'm here 2\n";
	    //return 0;
	    
	  }
	}
	if (thisr_in.Distance(thisnuexitice)>5.E2) {


	  if (WhereDoesItExitIce(thisnuexitice,nnu,5.E2, // then back up and find it more precisely
				  thisr_exitice, antarctica)) {
	    count2++;
	    //interaction1->neverseesice=1;

	    thisnuexitice=5.E2*nnu;
	    thisnuexitice+=thisr_exitice;
	    //cout << "inu is " << inu << " I'm here 3\n";
	    //return 0;
	  }
	
	}
	if (thisr_in.Distance(thisnuexitice)>50.) {


	  if (WhereDoesItExitIce(thisnuexitice,nnu,50., // then back up and find it more precisely
				  thisr_exitice, antarctica)) {
	    //interaction1->neverseesice=1;
	    count2++;
	    //cout << "inu is " << inu << " I'm here 4\n";
	    //return 0;
	  }
	}
	thisnuexitice=thisr_exitice;
	if (count2>10) cout << "count1 is " << count2 << "\n";
	//	else return 0;  // never reaches any ice or is it because our step is too big

      } // if the nu leaves a rock bin


    } // end wheredoesitleave


    else {
      wheredoesitleave_err=1;
      pickposnu = 0;
      return 0;
    }

    // end finding where it leaves ice

// 	if (thisnuexit.Mag()<Surface(thisnuexit)) { // if the exit point is below the surface
// 	  WhereDoesItExitIceForward(thisnuexit,interaction1->nnu,20., // then find it more finely
// 			     thisr_exitice);
// 	  thisnuexit=thisr_enterice;
// 	  // then back up and find it more precisely
// 	}

    if (WhereDoesItEnterIce(thisnuexitearth,nnu,5.E3, // first pass with sort of course binning
			    thisr_enterice, antarctica)) {
      thisr_enterice_tmp=thisr_enterice+5.E3*nnu;
      //cout << "inu is " << inu << " thisr_enterice is ";thisr_enterice.Print();
      if (WhereDoesItEnterIce(thisr_enterice_tmp,nnu,20., // second pass with finer binning
			      thisr_enterice, antarctica)) {
	//cout << "inu is " << inu << " thisr_enterice is ";thisr_enterice.Print();
	//cout << "entersice is ";thisr_enterice.Print();
	//cout << "thisnuexitice is ";thisnuexitice.Print();
	pathlength_inice=thisr_enterice.Distance(thisnuexitice);
	//cout << "distance is " << distance << "\n";
	//cout << "inu " << inu << " thisr_enterice, thisnuexitice are ";thisr_enterice.Print();thisnuexitice.Print();
	posnu=pathlength_inice*gRandom->Rndm()*nnu;
	posnu=posnu+thisr_enterice;
	//cout << "inu" << inu << " thisr_enterice, thisnuexitice are ";thisr_enterice.Print();thisnuexitice.Print();
	//cout << "inu " << inu << " distance is " << distance << "\n";
      }
    }
    else {
      thisr_enterice=thisr_in;
      wheredoesitenterice_err=1;
      pickposnu = 0;
      return 0;
    }

    nuexitice=thisnuexitice;
    r_enterice=thisr_enterice;
    
    if (posnu.Mag()-antarctica->Surface(posnu)>0) {
//--------------------------------------------------
//     if (interaction1->posnu.Mag()-Surface(interaction1->posnu)>0) {
//-------------------------------------------------- 
      toohigh=1;
      //cout << "inu, toohigh is " << inu << " " << interaction1->toohigh << "\n";
      pickposnu = 0;
      return 0;
    }
    if (posnu.Mag()-antarctica->Surface(posnu)+antarctica->IceThickness(posnu)<0) {
//--------------------------------------------------
//     if (interaction1->posnu.Mag()-Surface(interaction1->posnu)+antarctica->IceThickness(interaction1->posnu)<0) {
//-------------------------------------------------- 
      toolow=1;
      //cout << "inu, toolow is " << inu << " " << interaction1->toolow << "\n";
      pickposnu = 0;
      return 0;
    }    
        
    // in case the code reachs here, (finally pick the usable location), now we can find r_in as posnu is selected
    r_in = antarctica->WhereDoesItEnter(posnu, nnu);

    pickposnu = 1;
    return 1;

}




//void Interaction::PickNear (IceModel *antarctica, Detector *detector, Settings *settings1) {
void Interaction::PickNear_Cylinder (IceModel *antarctica, Detector *detector, Settings *settings1) {


    double range = settings1->POSNU_RADIUS;   // test value, 2km radius. can be changed to read from Settings

    //thisphi=gRandom->Rndm()*(maxphi-minphi)+minphi;

    //Interaction::PickAnyDirection(); // first pick the neutrino direction

    //pick random posnu within boundary 2km radius
    double thisPhi = gRandom->Rndm() * (2*PI);
    //double thisR = gRandom->Rndm() * (range);
    double thisR = pow( gRandom->Rndm(), 0.5 ) * (range);   // for uniform distribution

    double X, Y, D;    // X,Y wrt detector core, and it's distance D
    
    //calculate posnu's X, Y wrt detector core
    if (detector->Get_mode() == 1 || detector->Get_mode() == 2 || detector->Get_mode() == 3) {   // detector mode is for ARA stations;
        X = detector->params.core_x + thisR*cos(thisPhi);
        Y = detector->params.core_y + thisR*sin(thisPhi);
        D = pow(X*X + Y*Y, 0.5);
        //interaction1->posnu.SetThetaPhi( D/antarctica->Surface(0., 0.), atan2(Y,X) ); 
    }
    //calculate posnu's X, Y wrt to (0,0)
    else {  // for mode = 0 (testbed)
        X = thisR*cos(thisPhi);
        Y = thisR*sin(thisPhi);
        D = pow(X*X + Y*Y, 0.5);
    }


    //Interaction::FlattoEarth(antarctica, X, Y, D);  //change to Earth shape and set depth (always in the ice)
    if (settings1->PICK_POSNU_DEPTH == 0) {
        FlattoEarth(antarctica, X, Y, D);  //change to Earth shape and set depth (always in the ice)
    }
    else if (settings1->PICK_POSNU_DEPTH == 1) {
        FlattoEarth_Near_Surface(antarctica, X, Y, D, settings1->MAX_POSNU_DEPTH);  //change to Earth shape and set depth (posnu depth has maximum value as MAX_POSNU_DEPTH)
    }



    pickposnu = 1;  // all PickNear sucess for pickposnu

    // set the position where nu enter the earth
    r_in = antarctica->WhereDoesItEnter(posnu, nnu);

    // set the position where nu exit the earth
    nuexit = antarctica->WhereDoesItLeave(posnu, nnu);

    // now set the position where nu enter the ice
    if (antarctica->IceThickness(r_in) && r_in.Lat()<antarctica->GetCOASTLINE()) { // if r_in (position where nu enter the earth) is antarctic ice
        r_enterice = r_in;  // nu enter the earth is same with nu enter the ice
    }
    else {  // nu enter the rock of earth. so we have to calculate the r_enterice
        Position thisnuenterice_tmp1;
        Position thisnuenterice_tmp2;
        // now first rough calculation with step size 5.E4.
        if (WhereDoesItEnterIce(posnu,nnu,5.E4,
			    thisnuenterice_tmp1, antarctica)) {
            thisnuenterice_tmp2=thisnuenterice_tmp1+5.E4*nnu;   // get one more step from 5.E4. calculation
            
            if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E3, // second pass with finer binning
			      thisnuenterice_tmp1, antarctica)) {
                thisnuenterice_tmp2=thisnuenterice_tmp1+5.E3*nnu;   // get one more step from 5.E3. calculation
        
                if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E2, // third pass with finer binning
			        thisnuenterice_tmp1, antarctica)) {
                    thisnuenterice_tmp2=thisnuenterice_tmp1+5.E2*nnu;   // get one more step from 5.E2. calculation

                    if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E1, // fourth pass with finer binning (final)
			            thisnuenterice_tmp1, antarctica)) {
                        thisnuenterice_tmp2=thisnuenterice_tmp1;   // max 50m step result
                    }
                }
            }
        }
        else {  // no result from the first step calculation
            cout<<"no nuenterice result from calculation!!!"<<endl;
            thisnuenterice_tmp2 = posnu;
        }
        r_enterice = thisnuenterice_tmp2;
    }// else; nu enter the rock of earth, so calculated the ice enter point


    // now we have to calcuate the nu ice exit position
    if (antarctica->IceThickness(nuexit) && nuexit.Lat()<antarctica->GetCOASTLINE()) { // if nuexit (position where nu exit the earth) is antarctic ice
        nuexitice = nuexit;  // nu exit the earth is same with nu exit the ice
    }
    else {  // nu exit the rock of earth. so we have to calculate the nuexitice
        Position thisnuexitice_tmp1;
        Position thisnuexitice_tmp2;
        // now first rough calculation with step size 5.E4.
        if (WhereDoesItExitIceForward(posnu,nnu,5.E4,
			    thisnuexitice_tmp1, antarctica)) {
            thisnuexitice_tmp2=thisnuexitice_tmp1-5.E4*nnu;   // get one more step from 5.E4. calculation
            
            if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E3, // second pass with finer binning
			      thisnuexitice_tmp1, antarctica)) {
                thisnuexitice_tmp2=thisnuexitice_tmp1-5.E3*nnu;   // get one more step from 5.E3. calculation
        
                if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E2, // third pass with finer binning
			        thisnuexitice_tmp1, antarctica)) {
                    thisnuexitice_tmp2=thisnuexitice_tmp1-5.E2*nnu;   // get one more step from 5.E2. calculation

                    if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E1, // fourth pass with finer binning (final)
			            thisnuexitice_tmp1, antarctica)) {
                        thisnuexitice_tmp2=thisnuexitice_tmp1;   // max 50m step result
                    }
                }
            }
        }
        else {  // no result from the first step calculation
            cout<<"no nuexitice result from calculation!!!"<<endl;
            thisnuexitice_tmp2 = posnu;
        }
        nuexitice = thisnuexitice_tmp2;
    }// else; nu enter the rock of earth, so calculated the ice enter point



}








//double Interaction::PickNear (IceModel *antarctica, Detector *detector, Settings *settings1) {
double Interaction::PickNear_Sphere (IceModel *antarctica, Detector *detector, Settings *settings1) {

  double L0=0.;
  const double range = settings1->POSNU_RADIUS;   // test value, 2km radius. can be changed to read from Settings

    //thisphi=gRandom->Rndm()*(maxphi-minphi)+minphi;

    //Interaction::PickAnyDirection(); // first pick the neutrino direction

    //pick random posnu within boundary 2km radius
    double thisPhi = gRandom->Rndm() * (2*PI);
    //double thisR = gRandom->Rndm() * (range);
    double thisR = pow( gRandom->Rndm(), 0.5 ) * (range);   // for uniform distribution




    // fristly pick random posnu within boundary 2km radius
    double rndX, rndY;
    do{
    rndX = (2.*gRandom->Rndm() -1.) * range;  // [-range:range]
    rndY = (2.*gRandom->Rndm() -1.) * range;
    }while(rndX*rndX + rndY*rndY > range*range);

    //printf("ttt %f %f ttt\n", rndX, rndY);

    Vector rndPosition(rndX, rndY, 0.);

    // now transform the plane to be vertical to the incident angle of a particle
    Vector normnnu = nnu.Unit();
    rndPosition.RotateUz(normnnu);
    double transX = rndPosition.GetX();
    double transY = rndPosition.GetY();
    double transZ = rndPosition.GetZ();


    double X, Y, Z;    // X,Y wrt detector core
    //calculate posnu's X, Y wrt detector core
    if (detector->Get_mode() == 1 || detector->Get_mode() == 2 || detector->Get_mode() == 3) {   // detector mode is for ARA stations;
      X = detector->params.core_x + transX;
      Y = detector->params.core_y + transY;
      Z = transZ;
    }
    //calculate posnu's X, Y wrt to (0,0)
    else {  // for mode = 0 (testbed)
      X = transX;
      Y = transY;
      Z = transZ;
    }
    //fprintf(stdout, "X: %f, Y: %f, Z: %f\n", 
    //X - detector->params.core_x, Y - detector->params.core_y, Z);

    // judge whether neutrino can interact or not
    double newx, newy, newz;

    bool interacted = Does_Interact(X - detector->params.core_x, Y - detector->params.core_y, Z, 
                                  nnu.Theta(), nnu.Phi(), range,
                                  newx, newy, newz, L0);
    if(!interacted){
      pickposnu = 0;
      return 0.;
    }

    newx += detector->params.core_x;
    newy += detector->params.core_y;

    FlattoEarth_Spherical(antarctica, newx, newy, newz);  //change to Earth shape and set depth (always in the ice)

    if(posnu.Mag()-antarctica->Surface(posnu)+antarctica->IceThickness(posnu)<0){
      toolow=1;
      //cout << "inu, toolow is " << inu << " " << interaction1->toolow << "\n";
      pickposnu = 0;
      return 0.;
    }

    pickposnu = 1;  // all PickNear sucess for pickposnu

    // set the position where nu enter the earth
    r_in = antarctica->WhereDoesItEnter(posnu, nnu);

    // set the position where nu exit the earth
    nuexit = antarctica->WhereDoesItLeave(posnu, nnu);

    // now set the position where nu enter the ice
    if (antarctica->IceThickness(r_in) && r_in.Lat()<antarctica->GetCOASTLINE()) { // if r_in (position where nu enter the earth) is antarctic ice
        r_enterice = r_in;  // nu enter the earth is same with nu enter the ice
    }
    else {  // nu enter the rock of earth. so we have to calculate the r_enterice
        Position thisnuenterice_tmp1;
        Position thisnuenterice_tmp2;
        // now first rough calculation with step size 5.E4.
        if (WhereDoesItEnterIce(posnu,nnu,5.E4,
			    thisnuenterice_tmp1, antarctica)) {
            thisnuenterice_tmp2=thisnuenterice_tmp1+5.E4*nnu;   // get one more step from 5.E4. calculation
            
            if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E3, // second pass with finer binning
			      thisnuenterice_tmp1, antarctica)) {
                thisnuenterice_tmp2=thisnuenterice_tmp1+5.E3*nnu;   // get one more step from 5.E3. calculation
        
                if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E2, // third pass with finer binning
			        thisnuenterice_tmp1, antarctica)) {
                    thisnuenterice_tmp2=thisnuenterice_tmp1+5.E2*nnu;   // get one more step from 5.E2. calculation

                    if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E1, // fourth pass with finer binning (final)
			            thisnuenterice_tmp1, antarctica)) {
                        thisnuenterice_tmp2=thisnuenterice_tmp1;   // max 50m step result
                    }
                }
            }
        }
        else {  // no result from the first step calculation
            cout<<"no nuenterice result from calculation!!!"<<endl;
            thisnuenterice_tmp2 = posnu;
        }
        r_enterice = thisnuenterice_tmp2;
    }// else; nu enter the rock of earth, so calculated the ice enter point


    // now we have to calcuate the nu ice exit position
    if (antarctica->IceThickness(nuexit) && nuexit.Lat()<antarctica->GetCOASTLINE()) { // if nuexit (position where nu exit the earth) is antarctic ice
        nuexitice = nuexit;  // nu exit the earth is same with nu exit the ice
    }
    else {  // nu exit the rock of earth. so we have to calculate the nuexitice
        Position thisnuexitice_tmp1;
        Position thisnuexitice_tmp2;
        // now first rough calculation with step size 5.E4.
        if (WhereDoesItExitIceForward(posnu,nnu,5.E4,
			    thisnuexitice_tmp1, antarctica)) {
            thisnuexitice_tmp2=thisnuexitice_tmp1-5.E4*nnu;   // get one more step from 5.E4. calculation
            
            if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E3, // second pass with finer binning
			      thisnuexitice_tmp1, antarctica)) {
                thisnuexitice_tmp2=thisnuexitice_tmp1-5.E3*nnu;   // get one more step from 5.E3. calculation
        
                if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E2, // third pass with finer binning
			        thisnuexitice_tmp1, antarctica)) {
                    thisnuexitice_tmp2=thisnuexitice_tmp1-5.E2*nnu;   // get one more step from 5.E2. calculation

                    if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E1, // fourth pass with finer binning (final)
			            thisnuexitice_tmp1, antarctica)) {
                        thisnuexitice_tmp2=thisnuexitice_tmp1;   // max 50m step result
                    }
                }
            }
        }
        else {  // no result from the first step calculation
            cout<<"no nuexitice result from calculation!!!"<<endl;
            thisnuexitice_tmp2 = posnu;
        }
        nuexitice = thisnuexitice_tmp2;
    }// else; nu enter the rock of earth, so calculated the ice enter point

    return L0;
}









void Interaction::PickExact (IceModel *antarctica, Detector *detector, Settings *settings1, double thisR, double thisTheta, double thisPhi) {
    
    double range = settings1->POSNU_RADIUS;   // test value, 2km radius. can be changed to read from Settings
    
    //thisphi=gRandom->Rndm()*(maxphi-minphi)+minphi;
    
    //Interaction::PickAnyDirection(); // first pick the neutrino direction
    
    //pick random posnu within boundary 2km radius
    //double thisR = gRandom->Rndm() * (range);
   // for uniform distribution
    
    double X, Y, Z, D;    // X,Y wrt detector core, and it's distance D
    double sumX = 0.;
    double sumY = 0.; 
    double sumZ = 0.;
    int count = 0;

    for (int i = 0; i < detector->stations[0].strings.size(); i++){
        for (int j = 0; j < detector->stations[0].strings[i].antennas.size(); j++){
//            std::cout << "AntennaPos:X:Y:Z:: " << detector->stations[0].strings[i].antennas[j].GetX() <<  " : " << detector->stations[0].strings[i].antennas[j].GetY() << " : " << detector->stations[0].strings[i].antennas[j].GetZ() <<std::endl;
            sumX = sumX + detector->stations[0].strings[i].antennas[j].GetX();
            sumY = sumY + detector->stations[0].strings[i].antennas[j].GetY();
            sumZ = sumZ + detector->stations[0].strings[i].antennas[j].GetZ();
            count++;
    }
    }
    double avgX = sumX/double(count);
    double avgY = sumY/double(count);
    double avgZ = sumZ/double(count);
    
    
//    std::cout << "DetectorStation:X:Y:: "  << detector->stations[0].GetX() << " : " << detector->stations[0].GetY() << std::endl;
    

    //calculate posnu's X, Y wrt detector core
    if (detector->Get_mode() == 1 || detector->Get_mode() == 2 ||detector->Get_mode() == 3) {   // detector mode is for ARA stations;
//        X = detector->params.core_x + thisR*cos(thisPhi)*cos(thisTheta);
//        Y = detector->params.core_y + thisR*sin(thisPhi)*cos(thisTheta);
        X = avgX + thisR*cos(thisPhi)*cos(thisTheta);
        Y = avgY + thisR*sin(thisPhi)*cos(thisTheta);
        D = pow(X*X + Y*Y, 0.5);
        //interaction1->posnu.SetThetaPhi( D/antarctica->Surface(0., 0.), atan2(Y,X) ); 
    }
    //calculate posnu's X, Y wrt to (0,0)
    else {  // for mode = 0 (testbed)
        X = thisR*cos(thisPhi)*cos(thisTheta);
        Y = thisR*sin(thisPhi)*cos(thisTheta);
        D = pow(X*X + Y*Y, 0.5);
    }

//    double centerZ = (detector->stations[0].strings[0].antennas[1].GetZ()+ detector->stations[0].strings[0].antennas[2].GetZ())/2.;
    Z = avgZ + thisR*sin(thisTheta);
//    std::cout << "CenterPosition:X:Y:Z:: "  << avgX << " : " << avgY << " : " << avgZ <<  std::endl;

    //std::cout << "Central position: " << centerZ << " : " << Z << " : " << X << " : " << Y << std::endl;

    
    posnu.SetXYZ(X,Y,Z);
//    std::cout << "Posnu:X:Y:Z:: "  << posnu.GetX() << " : " << posnu.GetY() << " : " << posnu.GetZ() <<  std::endl;

    //Interaction::FlattoEarth(antarctica, X, Y, D);  //change to Earth shape and set depth (always in the ice)
    //FlattoEarth(antarctica, X, Y, D);  //change to Earth shape and set depth (always in the ice)
    
    
    
    
    
    pickposnu = 1;  // all PickNear sucess for pickposnu
    
    // set the position where nu enter the earth
    r_in = antarctica->WhereDoesItEnter(posnu, nnu);
    
    // set the position where nu exit the earth
    nuexit = antarctica->WhereDoesItLeave(posnu, nnu);
    
    // now set the position where nu enter the ice
    if (antarctica->IceThickness(r_in) && r_in.Lat()<antarctica->GetCOASTLINE()) { // if r_in (position where nu enter the earth) is antarctic ice
        r_enterice = r_in;  // nu enter the earth is same with nu enter the ice
    }
    else {  // nu enter the rock of earth. so we have to calculate the r_enterice
        Position thisnuenterice_tmp1;
        Position thisnuenterice_tmp2;
        // now first rough calculation with step size 5.E4.
        if (WhereDoesItEnterIce(posnu,nnu,5.E4,
                                thisnuenterice_tmp1, antarctica)) {
            thisnuenterice_tmp2=thisnuenterice_tmp1+5.E4*nnu;   // get one more step from 5.E4. calculation
            
            if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E3, // second pass with finer binning
                                    thisnuenterice_tmp1, antarctica)) {
                thisnuenterice_tmp2=thisnuenterice_tmp1+5.E3*nnu;   // get one more step from 5.E3. calculation
                
                if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E2, // third pass with finer binning
                                        thisnuenterice_tmp1, antarctica)) {
                    thisnuenterice_tmp2=thisnuenterice_tmp1+5.E2*nnu;   // get one more step from 5.E2. calculation
                    
                    if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E1, // fourth pass with finer binning (final)
                                            thisnuenterice_tmp1, antarctica)) {
                        thisnuenterice_tmp2=thisnuenterice_tmp1;   // max 50m step result
                    }
                }
            }
        }
        else {  // no result from the first step calculation
            cout<<"no nuenterice result from calculation!!!"<<endl;
            thisnuenterice_tmp2 = posnu;
        }
        r_enterice = thisnuenterice_tmp2;
    }// else; nu enter the rock of earth, so calculated the ice enter point
    
    
    // now we have to calcuate the nu ice exit position
    if (antarctica->IceThickness(nuexit) && nuexit.Lat()<antarctica->GetCOASTLINE()) { // if nuexit (position where nu exit the earth) is antarctic ice
        nuexitice = nuexit;  // nu exit the earth is same with nu exit the ice
    }
    else {  // nu exit the rock of earth. so we have to calculate the nuexitice
        Position thisnuexitice_tmp1;
        Position thisnuexitice_tmp2;
        // now first rough calculation with step size 5.E4.
        if (WhereDoesItExitIceForward(posnu,nnu,5.E4,
                                      thisnuexitice_tmp1, antarctica)) {
            thisnuexitice_tmp2=thisnuexitice_tmp1-5.E4*nnu;   // get one more step from 5.E4. calculation
            
            if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E3, // second pass with finer binning
                                          thisnuexitice_tmp1, antarctica)) {
                thisnuexitice_tmp2=thisnuexitice_tmp1-5.E3*nnu;   // get one more step from 5.E3. calculation
                
                if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E2, // third pass with finer binning
                                              thisnuexitice_tmp1, antarctica)) {
                    thisnuexitice_tmp2=thisnuexitice_tmp1-5.E2*nnu;   // get one more step from 5.E2. calculation
                    
                    if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E1, // fourth pass with finer binning (final)
                                                  thisnuexitice_tmp1, antarctica)) {
                        thisnuexitice_tmp2=thisnuexitice_tmp1;   // max 50m step result
                    }
                }
            }
        }
        else {  // no result from the first step calculation
            cout<<"no nuexitice result from calculation!!!"<<endl;
            thisnuexitice_tmp2 = posnu;
        }
        nuexitice = thisnuexitice_tmp2;
    }// else; nu enter the rock of earth, so calculated the ice enter point
    
    
    
}
/*
void Interaction::PickExactXYZ (IceModel *antarctica, Detector *detector, Settings *settings1, double X, double Y, double Z) {
    
    double range = settings1->POSNU_RADIUS;   // test value, 2km radius. can be changed to read from Settings
    
    //thisphi=gRandom->Rndm()*(maxphi-minphi)+minphi;
    
    //Interaction::PickAnyDirection(); // first pick the neutrino direction
    
    //pick random posnu within boundary 2km radius
    //double thisR = gRandom->Rndm() * (range);
    // for uniform distribution
    
    double X, Y, Z, D;    // X,Y wrt detector core, and it's distance D
    double sumX = 0.;
    double sumY = 0.; 
    double sumZ = 0.;
    int count = 0;
    
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            sumX = sumX + detector->stations[0].strings[i].antennas[j].GetX();
            sumY = sumY + detector->stations[0].strings[i].antennas[j].GetY();
            sumZ = sumZ + detector->stations[0].strings[i].antennas[j].GetZ();
            count++;
        }
    }
    double avgX = sumX/double(count);
    double avgY = sumY/double(count);
    double avgZ = sumZ/double(count);
    
    
    //    std::cout << "DetectorStation:X:Y:: "  << detector->stations[0].GetX() << " : " << detector->stations[0].GetY() << std::endl;
    
    
    //calculate posnu's X, Y wrt detector core
    if (detector->Get_mode() == 1 || detector->Get_mode() == 2) {   // detector mode is for ARA stations;
        //        X = detector->params.core_x + thisR*cos(thisPhi)*cos(thisTheta);
        //        Y = detector->params.core_y + thisR*sin(thisPhi)*cos(thisTheta);
        X = avgX + thisR*cos(thisPhi)*cos(thisTheta);
        Y = avgY + thisR*sin(thisPhi)*cos(thisTheta);
        D = pow(X*X + Y*Y, 0.5);
        //interaction1->posnu.SetThetaPhi( D/antarctica->Surface(0., 0.), atan2(Y,X) ); 
    }
    //calculate posnu's X, Y wrt to (0,0)
    else {  // for mode = 0 (testbed)
        X = thisR*cos(thisPhi)*cos(thisTheta);
        Y = thisR*sin(thisPhi)*cos(thisTheta);
        D = pow(X*X + Y*Y, 0.5);
    }
    
    //    double centerZ = (detector->stations[0].strings[0].antennas[1].GetZ()+ detector->stations[0].strings[0].antennas[2].GetZ())/2.;
    Z = avgZ + thisR*sin(thisTheta);
    
    //std::cout << "Central position: " << centerZ << " : " << Z << " : " << X << " : " << Y << std::endl;
    
    
    posnu.SetXYZ(X,Y,Z);
    //    std::cout << "Posnu:X:Y:Z:: "  << posnu.GetX() << " : " << posnu.GetY() << " : " << posnu.GetZ() <<  std::endl;
    
    //Interaction::FlattoEarth(antarctica, X, Y, D);  //change to Earth shape and set depth (always in the ice)
    //FlattoEarth(antarctica, X, Y, D);  //change to Earth shape and set depth (always in the ice)
    
    
    
    
    
    pickposnu = 1;  // all PickNear sucess for pickposnu
    
    // set the position where nu enter the earth
    r_in = antarctica->WhereDoesItEnter(posnu, nnu);
    
    // set the position where nu exit the earth
    nuexit = antarctica->WhereDoesItLeave(posnu, nnu);
    
    // now set the position where nu enter the ice
    if (antarctica->IceThickness(r_in) && r_in.Lat()<antarctica->GetCOASTLINE()) { // if r_in (position where nu enter the earth) is antarctic ice
        r_enterice = r_in;  // nu enter the earth is same with nu enter the ice
    }
    else {  // nu enter the rock of earth. so we have to calculate the r_enterice
        Position thisnuenterice_tmp1;
        Position thisnuenterice_tmp2;
        // now first rough calculation with step size 5.E4.
        if (WhereDoesItEnterIce(posnu,nnu,5.E4,
                                thisnuenterice_tmp1, antarctica)) {
            thisnuenterice_tmp2=thisnuenterice_tmp1+5.E4*nnu;   // get one more step from 5.E4. calculation
            
            if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E3, // second pass with finer binning
                                    thisnuenterice_tmp1, antarctica)) {
                thisnuenterice_tmp2=thisnuenterice_tmp1+5.E3*nnu;   // get one more step from 5.E3. calculation
                
                if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E2, // third pass with finer binning
                                        thisnuenterice_tmp1, antarctica)) {
                    thisnuenterice_tmp2=thisnuenterice_tmp1+5.E2*nnu;   // get one more step from 5.E2. calculation
                    
                    if (WhereDoesItEnterIce(thisnuenterice_tmp2,nnu,5.E1, // fourth pass with finer binning (final)
                                            thisnuenterice_tmp1, antarctica)) {
                        thisnuenterice_tmp2=thisnuenterice_tmp1;   // max 50m step result
                    }
                }
            }
        }
        else {  // no result from the first step calculation
            cout<<"no nuenterice result from calculation!!!"<<endl;
            thisnuenterice_tmp2 = posnu;
        }
        r_enterice = thisnuenterice_tmp2;
    }// else; nu enter the rock of earth, so calculated the ice enter point
    
    
    // now we have to calcuate the nu ice exit position
    if (antarctica->IceThickness(nuexit) && nuexit.Lat()<antarctica->GetCOASTLINE()) { // if nuexit (position where nu exit the earth) is antarctic ice
        nuexitice = nuexit;  // nu exit the earth is same with nu exit the ice
    }
    else {  // nu exit the rock of earth. so we have to calculate the nuexitice
        Position thisnuexitice_tmp1;
        Position thisnuexitice_tmp2;
        // now first rough calculation with step size 5.E4.
        if (WhereDoesItExitIceForward(posnu,nnu,5.E4,
                                      thisnuexitice_tmp1, antarctica)) {
            thisnuexitice_tmp2=thisnuexitice_tmp1-5.E4*nnu;   // get one more step from 5.E4. calculation
            
            if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E3, // second pass with finer binning
                                          thisnuexitice_tmp1, antarctica)) {
                thisnuexitice_tmp2=thisnuexitice_tmp1-5.E3*nnu;   // get one more step from 5.E3. calculation
                
                if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E2, // third pass with finer binning
                                              thisnuexitice_tmp1, antarctica)) {
                    thisnuexitice_tmp2=thisnuexitice_tmp1-5.E2*nnu;   // get one more step from 5.E2. calculation
                    
                    if (WhereDoesItExitIceForward(thisnuexitice_tmp2,nnu,5.E1, // fourth pass with finer binning (final)
                                                  thisnuexitice_tmp1, antarctica)) {
                        thisnuexitice_tmp2=thisnuexitice_tmp1;   // max 50m step result
                    }
                }
            }
        }
        else {  // no result from the first step calculation
            cout<<"no nuexitice result from calculation!!!"<<endl;
            thisnuexitice_tmp2 = posnu;
        }
        nuexitice = thisnuexitice_tmp2;
    }// else; nu enter the rock of earth, so calculated the ice enter point
    
    
    
}
*/




int Interaction::WhereDoesItLeave( const Position &posnu, const Vector &ntemp, IceModel *antarctica, Position &r_out) {

    
    double distance=0;
    double posnu_length=posnu.Mag(); // distance from center of earth to interaction
    
    double lon,lat,lon_old,lat_old; //latitude, longitude indices for 1st and 2nd iteration
    lon = posnu.Lon(); // what latitude, longitude does interaction occur at
    lat = posnu.Lat();
    lon_old=lon; // save this longitude and latitude so we can refer to it later
    lat_old=lat;
    
    // use law of cosines to get distance from interaction to exit point for the ray
    // need to solve for that distance using the quadratic formula
    
    // angle between posnu and ntemp vector for law of cosines.
    double costheta=-1*(posnu*ntemp)/posnu_length;
    
    // a,b,c for quadratic formula, needed to solve for 
    double a=1;
    double b=-1*2*posnu_length*costheta;
    double c=posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2);

    
    if (b*b-4*a*c<0.) {
   
      return 0;
    }
    // else
//       cout << "positive.  c is " << c << "\n";
   

    // use the "+" solution because the other one is where the ray is headed downward toward the rock
    distance=(-1*b+sqrt(b*b-4*a*c))/2;
    
    
    // now here is the exit point for the ray
    r_out = posnu + distance*ntemp;
    
    lon = r_out.Lon(); // latitude and longitude of exit point
    lat = r_out.Lat();
    
    
    c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
    // sometimes though the new surface is lower than the one over posnu which causes a problem.
    if (b*b-4*a*c<0.) {
      //cout << "inu is " << inu << "\n";  
      // try halving the distance
      distance=distance/2.;
      //cout << "bad.  distance 1/2 is " << distance << "\n";
      r_out = posnu + distance*ntemp;
      lon = r_out.Lon(); // latitude and longitude of exit point
      lat = r_out.Lat();
      c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
      if (b*b-4*a*c<0.) { // if we still have the problem back up more
	distance=distance/2.; // now we are at 1/4 the distance
	//cout << "bad.  distance 1/4 is " << distance << "\n";
	r_out = posnu + distance*ntemp;
	lon = r_out.Lon(); // latitude and longitude of exit point
	lat = r_out.Lat();
	c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	if (b*b-4*a*c<0.) { // the problem is less then 1/4 of the way in
	  
	  distance=distance/2.; // now we are at 1/8 the distance
	  //cout << "bad.  distance 1/8 is " << distance << "\n";
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines

	  if (b*b-4*a*c<0.) {
	    // still have the problem so just make the distance 0
	    distance=0.; // now we are at 1/8 the distance
	    //cout << "bad.  distance is " << distance << "\n";
	    lon = posnu.Lon(); // latitude and longitude of exit point
	    lat = posnu.Lat();
	    r_out=antarctica->Surface(lon,lat)/posnu.Mag()*posnu;
	  }
// 	  else
// 	    cout << "good.\n";





	} // now we are at 1/8 the distance
	else {// if this surface is ok problem is between 1/4 and 1/2
	  distance=distance*1.5; // now we are at 3/8 the distance
	  //	cout << "good.  distance 3/8 is " << distance << "\n";
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines

	  if (b*b-4.*a*c<0.) {
	    distance=distance*2./3.; // go back to 1/4
	    r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	  //cout << "good at distance 1/4 is " << distance << "\n";
	  }
// 	  else
// 	    cout << "good.\n";

	} // now we are at 3/8 the distance
      

      } // now we are at 1/4 the distance
      else { // if this surface at 1/2 distance is ok see if we can go a little further
	distance=distance*1.5; // now we are at 3/4 the distance
	//cout << "good.  distance 3/4 is " << distance << "\n";
	r_out = posnu + distance*ntemp;
	lon = r_out.Lon(); // latitude and longitude of exit point
	lat = r_out.Lat();
	c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	if (b*b-4*a*c<0.) { // the problem is between 1/2 and 3/4 of the way in
	  
	  distance=distance*5./6.; // now we are at 5/8 the distance
	  //cout << "bad.  distance 5/8 is " << distance << "\n";
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines

	  if (b*b-4*a*c<0.) {
	    distance=distance*4./5.;
 	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	  //cout << "good at distance 3/4 is " << distance << "\n";
	  }
// 	  else 
// 	    cout << "good at distance 5/8.\n";



	} // now we are at 1/8 the distance
	else {// if this surface is ok problem is between 1/4 and 1/2
	  distance=distance*7./6.; // now we are at 7/8 the distance
	  //cout << "good.  distance 7/8 is " << distance << "\n";
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines

	  if (b*b-4*a*c<0) {
	    // now found the problem so go back to 3/4 distance
	    distance=distance*6./7.;
	    //cout << "good at  distance 3/4 is " << distance << "\n";
	    r_out = posnu + distance*ntemp;
	    lon = r_out.Lon(); // latitude and longitude of exit point
	    lat = r_out.Lat();
	    c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines


	  }
// 	  else
// 	    cout << "good.\n";

	} // now we are at 3/8 the distance







      } // now we are at 3/4 distance
    
    
    } // if exit point we initially found was not ok
    else {
      distance=(-1*b+sqrt(b*b-4*a*c))/2; // and quadratic formula
      r_out = posnu + distance*ntemp;
    }
//--------------------------------------------------
//     if (inu==1159426) {
//       cout << "a,b,c,b*b-4*a*c are " << a << " " << b << " " << c << " " << b*b-4*a*c << "\n";
//       cout << "posnu is ";posnu.Print();
//       cout << "distance is " << distance << "\n";
//       cout << "ntemp is ";ntemp.Print();
//     }
//-------------------------------------------------- 
    

    
    
    
    return 1;
}


int Interaction::WhereDoesItEnterIce ( const Position &posnu, const Vector &nnu, double stepsize, Position &r_enterice_output, IceModel *antarctica) {

  // now get exit point...
  //   see my geometry notes.
  // parameterize the neutrino trajectory and just see where it
  // crosses the earth radius.

  //  Position r_enterice;
  double distance=0;
  int left_edge=0;
  Position x = posnu;
  double x2;
  
  Position x_previous = posnu;

  double x_previous2= x_previous * x_previous;
  x2=x_previous2;
  
  double lon = x.Lon(),lat = x.Lat();
  double lon_old = lon,lat_old = lat;
  double local_surface = antarctica->Surface(lon,lat);
//--------------------------------------------------
//   double local_surface = Surface(lon,lat);
//-------------------------------------------------- 
  double rock_previous2= pow((local_surface - antarctica->IceThickness(lon,lat) - antarctica->WaterDepth(lon,lat)),2);
  double surface_previous2=pow(local_surface,2);

  double rock2=rock_previous2;
  double surface2=surface_previous2;
  int foundit=0;  // keeps track of whether you found an ice entrance point

  //  cout << "lon, lat are " << posnu.Lon() << " " << posnu.Lat() << "\n";
  //cout << "x2 at start is " << x2 << "\n";
  while (distance<2*local_surface+1000) {

    distance+=stepsize;

    x -= stepsize*nnu;
    x2=x*x;
    //cout << "x2 is " << x2 << "\n";
    lon = x.Lon();
    lat = x.Lat();

      double ice_thickness=antarctica->IceThickness(lon,lat);
    if (lon!=lon_old || lat!=lat_old) {
      local_surface = antarctica->Surface(lon,lat);
//--------------------------------------------------
//       local_surface = Surface(lon,lat);
//-------------------------------------------------- 

      //if (lat>antarctica->GetCOASTLINE()) 
      //left_edge=1;

      rock2=pow((local_surface - antarctica->IceThickness(lon,lat) - antarctica->WaterDepth(lon,lat)),2);
      surface2=pow(local_surface,2);    

      if (antarctica->Getice_model()==0) {
	if ((int)(lat)==antarctica->GetCOASTLINE() && rock_previous2 < x2 && surface2 > x2)
	  left_edge=1;
      } //if (Crust 2.0)
    } //if (neutrino has stepped into new lon/lat bin)

    if ((((x_previous2>rock_previous2 && x2<rock2) // crosses rock boundary from above
	 || (x_previous2<surface_previous2 && x2>surface2)) && ice_thickness>0 && lat<antarctica->GetCOASTLINE()) // crosses surface boundary from below
	|| left_edge) {
      //  cout << "lat, antarctica->GetCOASTLINE(), left_edge is " << lat << " " << antarctica->GetCOASTLINE()<< " " << left_edge << "\n";
      //cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
      r_enterice_output = x;
      // this gets you out of the loop.
      //continue;
      distance=3*antarctica->Geoid(lat);
      foundit=1;
      //cout << "foundit is " << foundit << "\n";
      //cout << "r_enterice is ";r_enterice.Print();
      //continue;
    } //if

    x_previous = x;
    x_previous2 = x2;
    //cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";

    if (lon!=lon_old || lat!=lat_old) {
      rock_previous2 = rock2;
      surface_previous2 = surface2;
      lat_old = lat;
      lon_old = lon;
    } //if

  } //while

  return foundit;
}//WhereDoesItEnterIce


int Interaction::WhereDoesItExitIce ( const Position &posnu, const Vector &nnu, double stepsize, Position &r_enterice_output, IceModel *antarctica) {

  // now get exit point...
  //   see my geometry notes.
  // parameterize the neutrino trajectory and just see where it
  // crosses the earth radius.

  //  Position r_enterice;
  double distance=0;
  int left_edge=0;
  Position x = posnu;
  double x2;
  
//--------------------------------------------------
//   if (inu==1491) {
//    
//     cout << "posnu is";posnu.Print();
//     cout << "nnu is ";nnu.Print();
//   }
//-------------------------------------------------- 
   

  Position x_previous = posnu;

  double x_previous2= x_previous * x_previous;
  x2=x_previous2;
  
  double lon = x.Lon(),lat = x.Lat();
  double lon_old = lon,lat_old = lat;
  double local_surface = antarctica->Surface(lon,lat);
//--------------------------------------------------
//   double local_surface = Surface(lon,lat);
//-------------------------------------------------- 
  double rock_previous2= pow((local_surface - antarctica->IceThickness(lon,lat) - antarctica->WaterDepth(lon,lat)),2);
  double surface_previous2=pow(local_surface,2);

  double rock2=rock_previous2;
  double surface2=surface_previous2;
  int foundit=0;  // keeps track of whether you found an ice entrance point

 

  //  cout << "lon, lat are " << posnu.Lon() << " " << posnu.Lat() << "\n";
  //cout << "x2 at start is " << x2 << "\n";
  int nsteps=0;
  while (distance<2*local_surface+1000) {
    //cout << "another step.\n";
    distance+=stepsize;
    nsteps++;
    //    cout << "inu, nsteps is " << inu << " " << nsteps << "\n";
    x -= stepsize*nnu;
    x2=x*x;
    //cout << "x2 is " << x2 << "\n";
    lon = x.Lon();
    lat = x.Lat();

      double ice_thickness=antarctica->IceThickness(lon,lat);
    if (lon!=lon_old || lat!=lat_old) {
      local_surface = antarctica->Surface(lon,lat);
//--------------------------------------------------
//       local_surface = Surface(lon,lat);
//-------------------------------------------------- 

      //if (lat>antarctica->GetCOASTLINE()) 
      //left_edge=1;

      rock2=pow((local_surface - antarctica->IceThickness(lon,lat) - antarctica->WaterDepth(lon,lat)),2);
      surface2=pow(local_surface,2);    

      if (antarctica->Getice_model()==0) {
	if ((int)(lat)==antarctica->GetCOASTLINE() && rock_previous2 < x2 && surface2 > x2)
	  left_edge=1;
      } //if (Crust 2.0)
    } //if (neutrino has stepped into new lon/lat bin)

//--------------------------------------------------
//     if (inu==1491 && nsteps<10)
//       cout << "inu, x_previous2, rock_previous2, x2, rock2 are " << inu << " " << x_previous2 << " " << rock_previous2 << " " << x2 << " " << rock2 << "\n";
//-------------------------------------------------- 

    if ((((x_previous2<rock_previous2 && x2>rock2) // crosses rock boundary from above
	 || (x_previous2>surface_previous2 && x2<surface2)) && ice_thickness>0 && lat<antarctica->GetCOASTLINE()) // crosses surface boundary from above
	|| left_edge) {
      //  cout << "lat, antarctica->GetCOASTLINE(), left_edge is " << lat << " " << antarctica->GetCOASTLINE()<< " " << left_edge << "\n";
      //cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
      r_enterice_output = x;
      // this gets you out of the loop.
      //continue;
      distance=3*antarctica->Geoid(lat);
      foundit=1;
      //cout << "foundit is " << foundit << "\n";
      //continue;
    } //if

    x_previous = x;
    x_previous2 = x2;
    //cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";

    if (lon!=lon_old || lat!=lat_old) {
      rock_previous2 = rock2;
      surface_previous2 = surface2;
      lat_old = lat;
      lon_old = lon;
    } //if

  } //while
//--------------------------------------------------
//   if (inu==0) {
//     cout << "r_enterice is ";r_enterice.Print();}
//-------------------------------------------------- 
  return foundit;
}//WhereDoesItExitIce



int Interaction::WhereDoesItExitIceForward ( const Position &posnu, const Vector &nnu, double stepsize, Position &r_enterice_output, IceModel *antarctica) { 

    // this is fixed version of icemc -> icemodel -> WhereDoesItExitIceForward

  // now get exit point...
  //   see my geometry notes.
  // parameterize the neutrino trajectory and just see where it
  // crosses the earth radius.

  //  Position r_enterice;
  double distance=0;
  int left_edge=0;
  Position x = posnu;
  double x2;
  
//--------------------------------------------------
//   if (inu==1491) {
//    
//     cout << "posnu is";posnu.Print();
//     cout << "nnu is ";nnu.Print();
//   }
//-------------------------------------------------- 
   

  Position x_previous = posnu;

  double x_previous2= x_previous * x_previous;
  x2=x_previous2;
  
  double lon = x.Lon(),lat = x.Lat();
  double lon_old = lon,lat_old = lat;
  double local_surface = antarctica->Surface(lon,lat);
//--------------------------------------------------
//   double local_surface = Surface(lon,lat);
//-------------------------------------------------- 
  double rock_previous2= pow((local_surface - antarctica->IceThickness(lon,lat) - antarctica->WaterDepth(lon,lat)),2);
  double surface_previous2=pow(local_surface,2);

  double rock2=rock_previous2;
  double surface2=surface_previous2;
  int foundit=0;  // keeps track of whether you found an ice entrance point

 

  //  cout << "lon, lat are " << posnu.Lon() << " " << posnu.Lat() << "\n";
  //cout << "x2 at start is " << x2 << "\n";
  int nsteps=0;
  while (distance<2*local_surface+1000) {
    //cout << "another step.\n";
    distance+=stepsize;
    nsteps++;
    //    cout << "inu, nsteps is " << inu << " " << nsteps << "\n";
    //x -= stepsize*nnu;
    x += stepsize*nnu;  // should step forward (not backward)
    x2=x*x;
    //cout << "x2 is " << x2 << "\n";
    lon = x.Lon();
    lat = x.Lat();

      double ice_thickness=antarctica->IceThickness(lon,lat);
    if (lon!=lon_old || lat!=lat_old) {
      local_surface = antarctica->Surface(lon,lat);
//--------------------------------------------------
//       local_surface = Surface(lon,lat);
//-------------------------------------------------- 

      //if (lat>antarctica->GetCOASTLINE()) 
      //left_edge=1;

      rock2=pow((local_surface - antarctica->IceThickness(lon,lat) - antarctica->WaterDepth(lon,lat)),2);
      surface2=pow(local_surface,2);    

      if (antarctica->Getice_model()==0) {
	if ((int)(lat)==antarctica->GetCOASTLINE() && rock_previous2 < x2 && surface2 > x2)
	  left_edge=1;
      } //if (Crust 2.0)
    } //if (neutrino has stepped into new lon/lat bin)

//--------------------------------------------------
//     if (inu==1491 && nsteps<10)
//       cout << "inu, x_previous2, rock_previous2, x2, rock2 are " << inu << " " << x_previous2 << " " << rock_previous2 << " " << x2 << " " << rock2 << "\n";
//-------------------------------------------------- 

    if ((((x_previous2>rock_previous2 && x2<rock2) // crosses rock boundary from above
	 || (x_previous2<surface_previous2 && x2>surface2)) && ice_thickness>0 && lat<antarctica->GetCOASTLINE()) // crosses surface boundary from above
	|| left_edge) {
      //  cout << "lat, antarctica->GetCOASTLINE(), left_edge is " << lat << " " << antarctica->GetCOASTLINE()<< " " << left_edge << "\n";
      //cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
      r_enterice_output = x;
      // this gets you out of the loop.
      //continue;
      distance=3*antarctica->Geoid(lat);
      foundit=1;
      //cout << "foundit is " << foundit << "\n";
      //continue;
    } //if

    x_previous = x;
    x_previous2 = x2;
    //cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";

    if (lon!=lon_old || lat!=lat_old) {
      rock_previous2 = rock2;
      surface_previous2 = surface2;
      lat_old = lat;
      lon_old = lon;
    } //if

  } //while
//--------------------------------------------------
//   if (inu==0) {
//     cout << "r_enterice is ";r_enterice.Print();}
//-------------------------------------------------- 
  return foundit;
}//WhereDoesItExitIceForward





void Interaction::FlattoEarth ( IceModel *antarctica, double X, double Y, double D) {
    posnu.SetThetaPhi( D/antarctica->Surface(0.,0.), atan2(Y,X) );
    posnu.SetR( gRandom->Rndm() * antarctica->IceThickness(posnu.Lon(), posnu.Lat()) + (antarctica->Surface(posnu.Lon(), posnu.Lat()) - antarctica->IceThickness(posnu.Lon(), posnu.Lat()) ) );
}



void Interaction::FlattoEarth_Near_Surface ( IceModel *antarctica, double X, double Y, double D, double max_depth) {
    posnu.SetThetaPhi( D/antarctica->Surface(0.,0.), atan2(Y,X) );
    posnu.SetR( antarctica->Surface(posnu.Lon(), posnu.Lat()) - (gRandom->Rndm() * max_depth) );
}


void Interaction::FlattoEarth_Spherical ( IceModel *antarctica, double X, double Y, double Z) {
  posnu.SetThetaPhi( sqrt(X*X+Y*Y)/antarctica->Surface(0.,0.), atan2(Y,X) );
  posnu.SetR( antarctica->Surface(posnu.Lon(), posnu.Lat()) + (Z) ); // note Z is negative
}


     
void Interaction::PickAnyDirection() {
  double rndlist[2];
  gRandom->RndmArray(2,rndlist);
  
  costheta_nutraject=2*rndlist[0]-1;

 
  // pick a neutrino azimuthal angle
  phi_nutraject=2*PI*rndlist[1];
  
  // check that these give the right result
  double thetanu=acos(costheta_nutraject);
  
  double sinthetanu=sin(thetanu);
  
  // find direction vector of neutrino
  // **** are cosine and sine flipped?
  nnu.SetX(sinthetanu*cos(phi_nutraject));
  nnu.SetY(sinthetanu*sin(phi_nutraject));
  nnu.SetZ(costheta_nutraject);


}

// removed this function as it's set in Event class
//
/*
void  Interaction::setNuFlavor(Primaries *primary1,Settings *settings1,int whichray,Counting *counting1) {
     // pick the neutrino flavor,  type of tau decay when relevant,
      //  lpm energy.
      nuflavor=primary1->GetNuFlavor();
      
      if (settings1->MINRAY==whichray) {// only increment neutrino flavor on first ray so you don't count twice
      
	if (nuflavor=="nue")
	  counting1->nnu_e++;      
	if (nuflavor=="numu")
	  counting1->nnu_mu++;      
	if (nuflavor=="nutau")
	  counting1->nnu_tau++;      
      }
      
      if (settings1->FORSECKEL==1) // For making array of signal vs. freq, viewangle, just use muon neutrinos
	nuflavor="nue";
      if (nuflavor=="nue")  //For outputting to file
	nuflavorint=1;
      else if (nuflavor=="numu")
	nuflavorint=2;
      else if (nuflavor=="nutau")
	nuflavorint=3;


}
*/


void  Interaction::setCurrent(Primaries *primary1) {
     // pick whether it is neutral current
      // or charged current
      current=primary1->GetCurrent();


      if (current=="cc")   //For outputting to file
	//currentint=1;
	currentint=kCC;
      else if(current=="nc")
	//currentint=2;   
	currentint=kNC;
  }//setCurrent




bool Interaction::Does_Interact(double x, double y, double z,
				double theta, double phi, double r,
				double &newx, double &newy, double &newz, double &L0)
{
  const double a = sin(theta)*cos(phi);
  const double b = sin(theta)*sin(phi);
  const double c = cos(theta);

  const double B = a*x + b*y + c*z;
  const double C = x*x + y*y + z*z - r*r;

  //fprintf(stdout, "B: %f  C: %f\n", B, C);

  const double t1 = -B + sqrt(B*B - C);
  const double t2 = -B - sqrt(B*B - C);

  //fprintf(stdout, "t1: %f  t2: %f\n", t1, t2);

  const double z1 = z + c*t1;
  const double z2 = z + c*t2;

  //fprintf(stdout, "z1: %f  z2: %f\n", z1, z2);

  double ratio=0.;
  if(z1 > 0. && z2 > 0.){
    return false;  // does not interact
  }else if(z1 < 0. && z2 < 0.){
    ratio = 1.;
  }else{
    if(z1 > z2)
      ratio = -z2/(z1-z2);
    else
      ratio = -z1/(z2-z1);
      
    if(ratio < 0.){
      fprintf(stderr, "Ratio must be positive; %f", ratio);
      exit(0);
    }
  }

  double rndratio = gRandom->Rndm() * ratio;

  double abst = fabs(t1-t2);
  double newt = rndratio*abst;

  //fprintf(stderr, "ratio: %f, newratio: %f\n", ratio, rndratio);

  double x1 = x + a*t1;
  double y1 = y + b*t1; 
  double x2 = x + a*t2;
  double y2 = y + b*t2;
  if(z1 > z2){
    if(cos(theta)>0.){
      newx = x2 + a*newt;
      newy = y2 + b*newt;
      newz = z2 + c*newt;
    }else{
      newx = x2 - a*newt;
      newy = y2 - b*newt;
      newz = z2 - c*newt;
    }
  }else{
    if(cos(theta)>0.){
      newx = x1 + a*newt;
      newy = y1 + b*newt;
      newz = z1 + c*newt;
    }else{
      newx = x1 - a*newt;
      newy = y1 - b*newt;
      newz = z1 - c*newt;
    }
  }

  // derive the length on which neutrinos are generated (L0)
  if(z1 < 0. && z2 < 0.){ // track is all inside the earth
    L0 = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
  }else{
    double t0 = -z/c;  // t at surface (z0=0) // note c is not speed of light
    double x0 = x + a*t0;
    double y0 = y + b*t0;
    const double z0 = 0.;
    if(z1 > z2){
      L0 = sqrt( (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0));
    }else{
      L0 = sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0));
    }
  }

  //fprintf(stdout, "##### new position: %f %f %f #####\n", newx, newy, newz);
  //fprintf(stdout, "### L0 %f ###", L0);

  return true;
}








///////////////////// Y ////////////////////
Y::Y() { // Constructor

  ffrac=new TF1("ffrac","[0]*sin([1]*(x-[2]))",7.,12.); // This is the fraction of the distribution in the low y region given by Equation 18. 
  
  ffrac->FixParameter(0,0.128); // These parameters are the same for all interaction types
  ffrac->FixParameter(1,-0.197);
  ffrac->FixParameter(2,21.8);

  string sbase="C1_high";
  char which[50];
  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      sprintf(which,"%d%d",i,j);
      string sname=sbase+which;
      fC1_high[i][j]=new TF1(sname.c_str(),"[0]+[1]*(-exp(-(x-[2])/[3]))",7.,12.); // parameterization of parameter C1 in the high y region according to Equation 16
    }
  }

  /*
// parameter A_0 in Table V for the high y region
  fC1_high[0][0]->FixParameter(0,-0.0026);//nubar, CC
  fC1_high[1][0]->FixParameter(0,-0.008); //nu,    CC
  fC1_high[0][1]->FixParameter(0,-0.005); //nubar, NC
  fC1_high[1][1]->FixParameter(0,-0.005); //nu,    NC

  // parameter A_1 in Table V for the high y region
  fC1_high[0][0]->FixParameter(1,0.085); // nubar, CC
  fC1_high[1][0]->FixParameter(1,0.26); // nu, CC
  fC1_high[0][1]->FixParameter(1,0.23); // nubar, NC
  fC1_high[1][1]->FixParameter(1,0.23); // nu, NC


  // parameter A_2 in Table V for the high y region
  fC1_high[0][0]->FixParameter(2,4.1); // nubar, CC
  fC1_high[1][0]->FixParameter(2,3.0); // nu, CC   
  fC1_high[0][1]->FixParameter(2,3.0); // nubar, NC
  fC1_high[1][1]->FixParameter(2,3.0); // nu, NC
  */
// parameter A_0 in Table V for the high y region (fixed order)
  fC1_high[0][0]->FixParameter(0,-0.005); //nu,    NC
  fC1_high[1][0]->FixParameter(0,-0.005); //nubar, NC
  fC1_high[0][1]->FixParameter(0,-0.008); //nu,    CC
  fC1_high[1][1]->FixParameter(0,-0.0026);//nubar, CC

  // parameter A_1 in Table V for the high y region
  fC1_high[0][0]->FixParameter(1,0.23); // nu, NC
  fC1_high[1][0]->FixParameter(1,0.23); // nubar, NC
  fC1_high[0][1]->FixParameter(1,0.26); // nu, CC
  fC1_high[1][1]->FixParameter(1,0.085); // nubar, CC


  // parameter A_2 in Table V for the high y region
  fC1_high[0][0]->FixParameter(2,3.0); // nu, NC
  fC1_high[1][0]->FixParameter(2,3.0); // nubar, NC
  fC1_high[0][1]->FixParameter(2,3.0); // nu, CC   
  fC1_high[1][1]->FixParameter(2,4.1); // nubar, CC

  // parameter A_3 in Table V for the high y region.  This parameter is the same for all four interaction types
  for (int i=0;i<2;i++) { // nu, nubar
    for (int j=0;j<2;j++) { // NC, CC
      fC1_high[i][j]->FixParameter(3,1.7);
    }
  }

  fC1_low=new TF1("C1_low","[0]+[1]*(-exp(-(x-[2])/[3]))",7.,12.); // parameterization of parameter C1 in the low y region according to Equation 16.
  // This parameterization is the same for all interaction types.
  
  fC1_low->FixParameter(0,0.);
  fC1_low->FixParameter(1,0.0941);
  fC1_low->FixParameter(2,4.72);
  fC1_low->FixParameter(3,0.456);

  fC2=new TF1("C2","[0]+[1]*x",7.,12.); // parameterization of parameter C2 in the low y region according to Equation 17.
  // This parameterization is the same for all interaction types.
  fC2->FixParameter(0,2.55);
  fC2->FixParameter(1,-9.49E-2);

  // For picking inelasticity in low y region according to Equation 14.
   fy0_low=new TF3("fy0_low","x+(z*([1]-x)^(-1./y+1)+(1-z)*([0]-x)^(-1./y+1))^(y/(y-1))"); // x=C_1, y=C_2, z=R
   fy0_low->SetParameter(0,0.00002);  // y_min
   fy0_low->SetParameter(1,0.001); // y_max
   
   // For picking inelasticity in high y region according to Equation 15.
   fy0_high=new TF2("fy0_high","([1]-x)^y/([0]-x)^(y-1.)+x"); // x=C_1, y=R
   fy0_high->SetParameter(0,0.001); // y_min
   fy0_high->SetParameter(1,1.); // y_max

}//Y Constructor


Y::~Y() {

    delete ffrac;
    delete fC1_low;
    delete fC2;
    delete fy0_low;
    delete fy0_high;
  for(int i=0; i<2;i++){ // nu, nubar
    for(int j=0; j<2; j++){ // nc, cc
      delete fC1_high[i][j];
    }
  }

}


double Y::pickY(int NU,int CURRENT,double e) {
	//e is in GeV.
  // Select a value of y that follows the appropriate distribution according to the prescription outlined in 
  // A. Connolly, R. Thorne and D. Waters, arXiv:1102.0691 [hep-ph]. 

  // pick a y region
 
  //double R1=Rand3Y.Rndm(); // choose our first random number
  double R1=gRandom->Rndm();

  int iyregion=0; // 0 for high y region 
  if (R1<ffrac->Eval(e)) // Is it going to be in low y region?
    iyregion=1; // 1 for low y region
 
  double C1_this;
  if (iyregion==0) // high y region
    C1_this=fC1_high[NU][CURRENT]->Eval(e); // C1 for this event
  else // low y region
    C1_this=fC1_low->Eval(e); // C1 for this event

  double C2_this=fC2->Eval(e); // C2 for this event
  
  // pick another random number
	double R2=gRandom->Rndm();//  double R2=Rand3Y.Rndm();
  double y0=0.;
  if (iyregion==0)  // high y region
    y0=fy0_high->Eval(C1_this,R2); // pick y0 according to Equation 15
  else if (iyregion==1)  // low y region
    y0=fy0_low->Eval(C1_this,C2_this,R2); // pick y0 according to Equation 14

  return y0;
}//pickY
