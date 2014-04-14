#include "signal.hh"
//#include "vector.hh"
#include "TF1.h"
#include "TRandom3.h"
//#include "vector.hh"
#include "Vector.h"
//#include "position.hh"
#include "Position.h"
#include "Settings.h"
#include <fstream>


#include "Constants.h"

using std::cout;


const double Signal::N_AIR(1.);                      // index of refraction of air
const double Signal::NICE(1.79);                      // index of refraction of ice
const double Signal::CHANGLE_ICE(acos(1/NICE));                      // index of refraction of ice
const double Signal::NSALT(2.45);                   // index of refracton for salt
//const double RHOSALT=2165;                 // density of salt (kg/m**3)
const double Signal::RHOSALT(2050.);                 // density of salt (kg/m**3)
const double Signal::RHOICE(917);                     // density of ice (kg/m**3)
const double Signal::RHOH20(1000);          // density of water (kg/m**3)
const double Signal::RHOAIR(1.25);          // density of air (kg/m**3)
const double Signal::RM_ICE(10.35); // moliere radius, in g/cm^2
const double Signal::RM_SALT(12.09); // moliere radius, in g/cm^2
const double Signal::KR_SALT(1.33); // constant in jaime's parameterization
const double Signal::KR_ICE(1.42); // constant in jaime's parameterization
//const double Signal::X0SALT=0.1024;                // radiation length of salt (meters)
const double Signal::X0SALT(0.1081);                // radiation length of salt (meters)
//const double Signal::ECSALT=40.;                   // critical energy in salt (MeV)
const double Signal::ECSALT(38.5);                   // critical energy in salt (MeV)
const double Signal::X0ICE(0.403); 
// //const double Signal::X0ICE=0.392; // radiation length of ice (meters)
const double Signal::ECICE(63.7);                     // critical energy in ice (MeV)
// //const double Signal::ECICE=73.0; // critical energy in ice (MeV)
const double Signal::AEX_ICE(1.);  //efficiency for producing charge asymmetry relative to ice.  1 by definition
 
const double Signal::ALPHAICE(1.32); // exponent that goes into cutting off the spectrum at high frequencies
const double Signal::AEX_SALT(0.684);  // efficiency for producing charge asymmetry relative to ice
const double Signal::ALPHASALT(1.27); // exponent that goes into cutting off the spectrum at high frequencies
const double Signal::KE_SALT(3.2E-16); // constant in jaime's parameterization, in V/cm/MHz
const double Signal::KL_SALT(21.12); //constant in jaime's parameterization
const double Signal::KDELTA_SALT(14.95); // constant in jaime's parameterization
const double Signal::KE_ICE(4.79E-16); // constant in jaime's parameterization, in V/cm/MHz
const double Signal::KL_ICE(23.80); //constant in jaime's parameterization
const double Signal::KDELTA_ICE(18.33); // constant in jaime's parameterization

const double Signal::KELVINS_ICE(250.+150.);          // temperature in Kelvin (ice+system)
const double Signal::KELVINS_SALT(500.);            // temperature in salt (350) + receiver temp (150)
const double Signal::BETAICE(2.25); // exponent, in jaime's parameterization
// double Signal::NU0_MODIFIED=0.; // nu_0 modified for a specific medium
// double Signal::NU_R;// parameter for signal parameterization
const double Signal::BETASALT(2.60); // exponent, in jaime's parameterization
const double Signal::VIEWANGLE_CUT(sqrt(5.)); // require viewangle is no more than 5 delta away from the cerenkov angle where


Signal::Signal() : N_DEPTH(1.79) {

  Initialize();
    
}



Signal::Signal(Settings *settings1) : N_DEPTH(1.79) {

  Initialize(settings1);

    
    
}

Signal::~Signal() {

}



 void Signal::InitializeMedium() {
  if (MEDIUM==1) {
    SetKelvins(KELVINS_SALT);
    
    SetRhoMedium(RHOSALT);
    SetNDepth(NSALT);
    //    changle=changle_salt;
    SetNMediumReceiver(NSALT);
    SetX0Medium(X0SALT);
    SetEcMedium(ECSALT);
    SetAexMedium(AEX_SALT);
    SetAlphaMedium(ALPHASALT);
    SetRmMedium(RM_SALT);
    SetKeMedium(KE_SALT); // constant in jaime's parameterization, in V/cm/MHz
    SetKlMedium(KL_SALT); //constant in jaime's parameterization
    SetKdelta_Medium(KDELTA_SALT); // constant in jaime's parameterization
    SetKrMedium(KR_SALT); // constant in jaime's parameterization
    SetBetaMedium(BETASALT); // exponent, in jaime's parameterization

  } //if (MEDIUM==1)
  else if (MEDIUM==0) {
    SetKelvins(KELVINS_ICE);
    SetNDepth(NICE);
    //changle=CHANGLE_ICE;
    SetRhoMedium(RHOICE);
    
    SetNMediumReceiver(N_AIR);
    SetX0Medium(X0ICE);
    SetEcMedium(ECICE);
    SetAexMedium(AEX_ICE);
    SetAlphaMedium(ALPHAICE);
    SetRmMedium(RM_ICE);
    SetKeMedium(KE_ICE); // constant in jaime's parameterization, in V/cm/MHz
    SetKlMedium(KL_ICE); //constant in jaime's parameterization
    SetKdelta_Medium(KDELTA_ICE); // constant in jaime's parameterization
    SetKrMedium(KR_ICE); // constant in jaime's parameterization
    SetBetaMedium(BETAICE); // exponent, in jaime's parameterization
  } //if (MEDIUM==0)
 
}

 void Signal::Initialize() {

  ReadCalPulserSpectrum();
     
  logscalefactor_taper=0.;
  JAIME_FACTOR=1.0; // factor to multiply Jaime's parameterization for error analysis

  x0ice=0.403; 
  //X0ICE=0.392; // radiation length of ice (meters)
  ecice=63.7;                     // critical energy in ice (MeV)
  //const static ECICE=73.0; // critical energy in ice (MeV)
  nice=1.79;                      // index of refraction of ice
  nfirn=1.3250;                   // index of refraction at the very surface - Peter
  invnfirn=1/nfirn; 
  invnice=1/nice;
  rhoice=917;                     // density of ice (kg/m**3)
  kelvins_ice=250.+150.;          // temperature in Kelvin (ice+system)
  changle_ice=acos(1./nice);
  aex_ice=1.;  //efficiency for producing charge asymmetry relative to ice.  1 by definition

  alphaice=1.32; // exponent that goes into cutting off the spectrum at high frequencies
  rm_ice=10.35; // moliere radius, in g/cm^2
  ke_ice=4.79E-16; // const staticant in jaime's parameterization, in V/cm/MHz
  kl_ice=23.80; //const staticant in jaime's parameterization
  kdelta_ice=18.33; // const staticant in jaime's parameterization
  kr_ice=1.42; // const staticant in jaime's parameterization
  betaice=2.25; // exponent, in jaime's parameterization
  nu0_modified=0.; // nu_0 modified for a specific medium

  freq_reference=1.E6; // reference frequency in MHz
  pnu_reference=1.E18; // reference energy in MHz



  if (WHICHPARAMETERIZATION==1) {
    nu_r=(RHOMEDIUM/1000.)
      //NU_R=(RHOMEDIUM/1000.) // density in g/cm^3
      /KR_MEDIUM/RM_MEDIUM*
      CLIGHT*100./N_DEPTH/sin(acos(1/N_DEPTH));
 
    vmmhz1m_reference=KE_MEDIUM/ECMEDIUM* // KE in V/cm/MHz^2, Ec in MeV
      (X0MEDIUM*100.) // radiation length in cm
      *freq_reference/1.E6 // frequency in MHz
      *sqrt(N_DEPTH*N_DEPTH-1)/N_DEPTH // sin(theta_c)
      *pnu_reference/1.E6 // energy in MeV
      *1./sin(changle); 

    //cout << "multiplying by 1/changle which is " << 1./sin(changle) << "\n";

    //    vmmhz1m*=1./(1.+pow(freq/NU_R,ALPHAMEDIUM));
    vmmhz1m_reference*=1./(1.+pow(freq_reference/nu_r,ALPHAMEDIUM));

  }
  else {
      //cout<<"whichparameterization : "<<WHICHPARAMETERIZATION<<"\n";
  }

  


}



 void Signal::Initialize(Settings *settings1) {

     ReadCalPulserSpectrum();

     
  SetParameterization(settings1->WHICHPARAMETERIZATION);

  logscalefactor_taper=0.;
  JAIME_FACTOR=1.0; // factor to multiply Jaime's parameterization for error analysis

  x0ice=0.403; 
  //X0ICE=0.392; // radiation length of ice (meters)
  ecice=63.7;                     // critical energy in ice (MeV)
  //const static ECICE=73.0; // critical energy in ice (MeV)
  nice=1.79;                      // index of refraction of ice
  nfirn=1.3250;                   // index of refraction at the very surface - Peter
  invnfirn=1/nfirn; 
  invnice=1/nice;
  rhoice=917;                     // density of ice (kg/m**3)
  kelvins_ice=250.+150.;          // temperature in Kelvin (ice+system)
  changle_ice=acos(1./nice);
  aex_ice=1.;  //efficiency for producing charge asymmetry relative to ice.  1 by definition

  alphaice=1.32; // exponent that goes into cutting off the spectrum at high frequencies
  rm_ice=10.35; // moliere radius, in g/cm^2
  ke_ice=4.79E-16; // const staticant in jaime's parameterization, in V/cm/MHz
  kl_ice=23.80; //const staticant in jaime's parameterization
  kdelta_ice=18.33; // const staticant in jaime's parameterization
  kr_ice=1.42; // const staticant in jaime's parameterization
  betaice=2.25; // exponent, in jaime's parameterization
  nu0_modified=0.; // nu_0 modified for a specific medium

  freq_reference=1.E6; // reference frequency in MHz
  pnu_reference=1.E18; // reference energy in MHz

  SetLPM(settings1->LPM);   // set LPM effect on/off value


  if (WHICHPARAMETERIZATION==1) {
    nu_r=(RHOMEDIUM/1000.)
      //NU_R=(RHOMEDIUM/1000.) // density in g/cm^3
      /KR_MEDIUM/RM_MEDIUM*
      CLIGHT*100./N_DEPTH/sin(acos(1/N_DEPTH));
 
    vmmhz1m_reference=KE_MEDIUM/ECMEDIUM* // KE in V/cm/MHz^2, Ec in MeV
      (X0MEDIUM*100.) // radiation length in cm
      *freq_reference/1.E6 // frequency in MHz
      *sqrt(N_DEPTH*N_DEPTH-1)/N_DEPTH // sin(theta_c)
      *pnu_reference/1.E6 // energy in MeV
      *1./sin(changle); 

    //cout << "multiplying by 1/changle which is " << 1./sin(changle) << "\n";

    //    vmmhz1m*=1./(1.+pow(freq/NU_R,ALPHAMEDIUM));
    vmmhz1m_reference*=1./(1.+pow(freq_reference/nu_r,ALPHAMEDIUM));

  }
  else {
      //cout<<"whichparameterization : "<<WHICHPARAMETERIZATION<<"\n";
  }

  


}


 void Signal::ReadCalPulserSpectrum(){
    int frequency;
    double vmmhz1m_temp;
    int count = 0;
    ifstream infile("calpulser_spectrum.txt");
    if (infile){
        while (1) {
            infile >> frequency >> vmmhz1m_temp;
            if (!infile.good()) break;
            CalpulserVmMHz1m[count] = vmmhz1m_temp;
            count++;
        }
    } else {
        std::cerr << "No calpulser spectrum file!" << std::endl;
    }
}



void Signal::GetVmMHz(double vmmhz_max,double vmmhz1m_max,double pnu,double *freq,double notch_min,double notch_max,double *vmmhz,int nfreq) {

  // parametrization from Jaime Alvarez Munhiz  
  //  here using astro-ph/0003315 
  
  for (int i=0;i<nfreq;i++) {
  
    vmmhz[i]=vmmhz_max
      //*1./FREQ_LOW*freq[i];

      
      *GetVmMHz1m(pnu,freq[i])/vmmhz1m_max;
    //if (WHICHPARAMETERIZATION==0)
    //vmmhz[i]*=(1./(1.+pow(freq[i]/NU0_MODIFIED,ALPHAMEDIUM)));
    //if (WHICHPARAMETERIZATION==1)
    //vmmhz[i]*=1./(1.+pow(freq[i]/NU_R,ALPHAMEDIUM));

    
    if (notch_min!=0 && notch_max!=0 && freq[i]>notch_min && freq[i]<notch_max)
      vmmhz[i]=0.;
  } //for

//   double sum[5]={0.};
//   if (WHICHPATH==4) {
//     for (int j=0;j<5;j++) {
//       for (int i=0;i<NFREQ;i++) {
// 	if (bwslice_min[j]<=freq[i] && bwslice_max[j]>freq[i]) { 
// 	  sum[j]+=GetVmMHz1m(pnu,freq[i],x0ice,ecice,N_DEPTH,AEXMEDIUM,WHICHPARAMETERIZATION)*(freq[i+1]-freq[i])/1.E6;
// 	}
//       }
//       cout << "j, sum is " << j << " " << sum[j] << "\n";
//     }
    
//   }
  
} //GetVmMHz

 double Signal::GetELPM() {

  // LPM
  // elpm =7.7 TeV/cm * rho * X0 in PDG, but our x0 is in meters
  // so use elpm =  7.7TeV/cm*X0 
  // X0 is radiation lengths in cm

  //double elpm=7.7E12*(X0ICE*100.);

  double elpm=2.E15*(X0MEDIUM/x0ice);  // this is what Jaime uses.  see caption under figure 4 of 0003315.
  return elpm;
} //GetELPM
 int Signal::GetLPM() {


  return LPM;
} //GetLPM
void Signal::GetSpread(double pnu,
	       double emfrac,
	       double hadfrac,
	       double freq,
		       //	       double n_depth,
		       //double X0DEPTH,

	       double& deltheta_em_max,
	       double& deltheta_had_max) {
  
  //  cout << KR_MEDIUM << " " << RM_MEDIUM << " " << KL_MEDIUM << " " << KE_MEDIUM << " " << ECMEDIUM << " " << X0MEDIUM << " " << ALPHAMEDIUM << " " << AEXMEDIUM << " " << KDELTA_MEDIUM << " " << BETAMEDIUM << " " << KELVINS << " " << JAIME_FACTOR << "\n";
  //cout << RHOSALT << " " << RHOICE << " " << RM_ICE << " " << RM_SALT << " " << KR_SALT << " " << KR_ICE << " " << X0SALT << " " << ECSALT << " " << X0ICE << " " << ECICE << " " << AEX_ICE << "\n";  
  // cout << ALPHAICE << " " << AEX_SALT << " " << ALPHASALT << " " << KE_SALT << " " << KL_SALT << " " << KDELTA_SALT << " " << KE_ICE << " " << KL_ICE << " " << KDELTA_ICE << " " << KELVINS_SALT << " " << BETAICE << " " << BETASALT << "\n";
  
  //  scale by how far off Cherenkov angle this viewing antenna is
  //  c.f. A-MZ  astro-ph/9706064 and astro-ph/0003315
  //  and for non-LPM (non-EM) showers from 
  // Phys.Lett.B434,396 (1998)  (astro-ph/9806098)
  //  The lengths are different hence the angular thickness of 
  //  the shower is different.  Get the angular thickness for
  //  both the EM and hadroic parts.

//   cout << KR_MEDIUM << " " << RM_MEDIUM << " " << KL_MEDIUM << " " << KE_MEDIUM << " " << ECMEDIUM << " " << X0MEDIUM << " " << ALPHAMEDIUM << " " << AEXMEDIUM << " " << KDELTA_MEDIUM << " " << BETAMEDIUM << " " << KELVINS << " " << JAIME_FACTOR << "\n";
//   cout << RHOSALT << " " << RHOICE << " " << RM_ICE << " " << RM_SALT << " " << KR_SALT << " " << KR_ICE << " " << X0SALT << " " << ECSALT << " " << X0ICE << " " << ECICE << " " << AEX_ICE << "\n";  
//   cout << ALPHAICE << " " << AEX_SALT << " " << ALPHASALT << " " << KE_SALT << " " << KL_SALT << " " << KDELTA_SALT << " " << KE_ICE << " " << KL_ICE << " " << KDELTA_ICE << " " << KELVINS_SALT << " " << BETAICE << " " << BETASALT << "\n";

//--------------------------------------------------
//   double elpm=GetLPM();
//-------------------------------------------------- 
  double elpm=GetELPM();

    //cout << "elpm is " << elpm << "\n";


  //  cout << "elpm is " << elpm << "\n";
  freq=freq/1.E6;  // frequency in MHz
  double showerlength=3.1;  //shower length in meters-gets a modification
                            //for em showers due to lpm effect.
  // this shower length is chosen somewhat arbitrarily, but is 
  // approximately the length of a shower in ice.
  // Then, the coefficient out front of the equations for
  // deltheta_em_max and deltheta_had_max are set so that
  // for ice, we get the equations in astro-ph/9706064
  // with the coefficient in front being 2.7 degrees.
  // I wanted to make the dependence on the shower length
  // and index of refraction explicit, so I pulled those variables
  // out of the equations for deltheta_em_max and deltheta_had_max.

  double em_eshower;  // em shower energy
  double had_eshower; // had shower energy
  double nu0; // reference frequency

  em_eshower=emfrac*pnu; // first, consider the electromagnetic shower.
  had_eshower=hadfrac*pnu;  // just the energy of the hadronic component of the shower

  // lengthen the shower to account for the lpm effect.
  // from astro-ph/9706064
  if (em_eshower<1.E15 || !LPM) 
    showerlength/=pow((em_eshower/1.e15),-0.03); 
  else 
    showerlength/=pow(elpm/(0.14*(em_eshower)+elpm),0.3);

  //  cout << "showerlength is " << showerlength << "\n";


    if (WHICHPARAMETERIZATION==0) {

      nu0=500.E6/1.E6*x0ice/X0MEDIUM; // for rego (astro-ph/9706064)

      // decoherence frequency scales with radiation length
      // when X0MEDIUM=X0ICE, nu0=500 MHz as in astro-ph/9706064


      // these equations are in astro-ph/9706064, but we have pulled
      // out the dependence on index of refraction and shower length. 
      // note that 12.32/sqrt(pow(n_depth,2)-1)*RADDEG/showerlength=2.7 degrees.
      // remember that Jaime has a factor of ln2 in the exponential here which we'll have to correct for further down
      deltheta_em_max=12.32/sqrt(pow(N_DEPTH,2)-1)*(nu0/freq)*RADDEG/showerlength;
      //cout<<"1) daltheta_em_max : "<<deltheta_em_max<<endl;

      if (hadfrac>0.00001) { // if there is a hadronic component
	
	
	// these equations are in astro-ph/9806098, but we have pulled
	// out the dependence on index of refraction and shower length.
	// remember that in this paper he includes a factor of ln2 in
	// the exponential, which we account for further down
	double epsilon=log10(had_eshower/1.E12);
	if (had_eshower>=1E12 && had_eshower<100.E12) 
	  deltheta_had_max=1.473/sqrt(pow(N_DEPTH,2)-1)*nu0/freq*RADDEG*(2.07-0.33*epsilon+(7.5e-2)*epsilon*epsilon);
	else if (had_eshower<100.E15) 
	  deltheta_had_max=1.473/sqrt(pow(N_DEPTH,2)-1)*nu0/freq*RADDEG*(1.744-(1.21e-2)*epsilon);
	else if (had_eshower<10.E18)   
	  deltheta_had_max=1.473/sqrt(pow(N_DEPTH,2)-1)*nu0/freq*RADDEG*(4.23-0.785*epsilon+(5.5e-2)*epsilon*epsilon);
	else {
	  //  beyond param, just use value at 10 EeV since slow variation
	  //  and parameterization might diverge
	  //  so scale from 10 EeV at 7.5% per decade (30/4=7.5)
	  //deltheta_had_max=1.473/sqrt(pow(N_DEPTH,2)-1)*nu0/freq*RADDEG*(4.23-0.785*7.+5.5e-2*49.);  // the last part in parenthesis if the previous equation evaluated at epsilon=7.
	  //deltheta_had_max=deltheta_had_max*(1.+(epsilon-7.)*0.075);
	  // It doesn't increase deltheta_had_max by 7.5% per decade anymore. Now it decreases the energy factor by 0.07 per decade.
	  deltheta_had_max=1.473/sqrt(pow(N_DEPTH,2)-1)*nu0/freq*RADDEG*(4.23-0.785*7.+5.5e-2*49. - (epsilon-7.)*0.07);
	} //else : beyond paramatrization
	deltheta_had_max/=sqrt(log(2.)); // in astro-ph/9706064, Jaime uses exp(-0.5* (theta-theta_c)^2/delta_had^2)

	// we adjust the delta's so that we can use the same form for both parameterizations: exp(-(theta-theta_c)^2/delta^2)

      }
      else
	deltheta_had_max=1.E-10;

      deltheta_em_max/=sqrt(log(2.)); // in astro-ph/9706064, Jaime uses exp(-0.5 (theta-theta_c)^2/delta_had^2)
      //cout<<"2) daltheta_em_max : "<<deltheta_em_max<<endl;

    }
    else if (WHICHPARAMETERIZATION==1) {

      //  cout << "I'm here inside GetSpread.\n";
      // we use the old parameterization for em showers
      nu0=500.E6/1.E6; // for rego (astro-ph/9706064)

      deltheta_em_max=12.32/sqrt(nice*nice-1)*(nu0/freq)*RADDEG/showerlength;


      // and then scale it according to astro-ph/0512337
      // Eq. 9
      deltheta_em_max*=RHOMEDIUM/rhoice
	/KDELTA_MEDIUM*kdelta_ice
	/X0MEDIUM*x0ice
	/sqrt(N_DEPTH*N_DEPTH-1)*sqrt(nice*nice-1);

      if (hadfrac>0.00001) { // if there is a hadronic component
      // for had showers, just use the one from astro-ph/0512337
      // Eq. 9
      // straight away
	deltheta_had_max=CLIGHT*100.// speed of light in cm/s
	  /(freq*1.E6)
	  *1/KDELTA_MEDIUM
	  /(X0MEDIUM*100.) // radiation length in cm
	  /sqrt(N_DEPTH*N_DEPTH-1.);   
	
      } //if (hadronic component)
      else 
	deltheta_had_max=1.E-10;
      
    } // if more recent parameterization


    


} //GetSpread


 double Signal::GetVmMHz1m(double pnu,double freq) { // constructor

  if (WHICHPARAMETERIZATION==0) {
    // parametrization from Jaime Alvarez Munhiz  
    //  here using astro-ph/0003315 
    double nu0=1150.E6/1.E6;
    //NU0_MODIFIED=nu0
    double nu0_modified=nu0
      *(x0ice/ecice)/(X0MEDIUM/ECMEDIUM)
      *(1/sqrt(N_DEPTH*N_DEPTH-1.))/(1/sqrt(nice*nice-1.));

    freq=freq/1.E6;  // frequency in MHz
    
    double factor=
      //1/sin(changle) // should be cerenkov angle for ice
      1/sqrt(1-1/(nice*nice)) // sin(changle) for ice
      *1/nu0 //
      *X0MEDIUM/x0ice  // track length *** use dE/dX rho instead
      *ecice/ECMEDIUM
      *AEXMEDIUM/aex_ice;  // to account for critical energy
    // to account for cerenkov threshold // maybe should be "a" instead

    vmmhz1m_max=factor*(2.53E-7)*(pnu/1.E12)*freq
      //      *(1./(1.+pow(freq/NU0_MODIFIED,ALPHAMEDIUM)))
      *(1./(1.+pow(freq/nu0_modified,1.44)));
  }
  else if (WHICHPARAMETERIZATION==1) {

 
    vmmhz1m_max=vmmhz1m_reference
      *freq/freq_reference
      *pnu/pnu_reference
      *1./(1.+pow(freq/nu_r,ALPHAMEDIUM))
      *(1.+pow(freq_reference/nu_r,ALPHAMEDIUM));

  }


  vmmhz1m_max=vmmhz1m_max/2.;  // This factor of 2 is to account for the 2 in the definition of the fourier transform in Equation 8 of the Halzen, Stanev and Zas paper.  The same factor of 2 seems to have propagated through subsequent Jaime papers.
  vmmhz1m_max=vmmhz1m_max*sqrt(2.);  // This is to account for the fact that the E fields quoted in the theory papers are double-sided in frequency (they extend from -F to F) whereas we are using it as a single-sided E-field (only from 0 to F).

  //  cout << "jaime_factor is " << JAIME_FACTOR << "\n";
  return vmmhz1m_max*JAIME_FACTOR;


  //  vmmhz1m=vmmhz1m/sqrt(1.E6/(BW/(double)NFREQ)); //THIS IS NEEDED TO CONSERVE ENERGY FOR DIFFERENT BIN WIDTHS.




//      // this is the old version
//      double factor=
//        X0MEDIUM/X0ICE  // track length
//        *(1-1/(N_DEPTH*N_DEPTH))/(1-1/(NICE*NICE)) // cerenkov index of refraction factor
//        *N_DEPTH/NICE // to account for cerenkov threshold
//        *ECICE/ECMEDIUM;  // to account for critical energy
    
//      double vmmhz1m=factor*(2.53E-7)*(pnu/1.E12)*(freq/nu0)*(1./(1.+pow(freq/nu0_modified,1.44)))*JAIME_FACTOR;

} //Signal constructor


double Signal::GetVmMHz1mCalPulser(int bin) { // constructor
    
    vmmhz1m_max = CalpulserVmMHz1m[bin];
//    std::cout << bin << " : " << vmmhz1m_max << std::endl;
    
//    vmmhz1m_max=vmmhz1m_max/2.;  // This factor of 2 is to account for the 2 in the definition of the fourier transform in Equation 8 of the Halzen, Stanev and Zas paper.  The same factor of 2 seems to have propagated through subsequent Jaime papers.
//    vmmhz1m_max=vmmhz1m_max*sqrt(2.);  // This is to account for the fact that the E fields quoted in the theory papers are double-sided in frequency (they extend from -F to F) whereas we are using it as a single-sided E-field (only from 0 to F).
    
    //  cout << "jaime_factor is " << JAIME_FACTOR << "\n";
//    return vmmhz1m_max*JAIME_FACTOR;
    return vmmhz1m_max;  
    
    //  vmmhz1m=vmmhz1m/sqrt(1.E6/(BW/(double)NFREQ)); //THIS IS NEEDED TO CONSERVE ENERGY FOR DIFFERENT BIN WIDTHS.
    
    
    
    
    //      // this is the old version
    //      double factor=
    //        X0MEDIUM/X0ICE  // track length
    //        *(1-1/(N_DEPTH*N_DEPTH))/(1-1/(NICE*NICE)) // cerenkov index of refraction factor
    //        *N_DEPTH/NICE // to account for cerenkov threshold
    //        *ECICE/ECMEDIUM;  // to account for critical energy
    
    //      double vmmhz1m=factor*(2.53E-7)*(pnu/1.E12)*(freq/nu0)*(1./(1.+pow(freq/nu0_modified,1.44)))*JAIME_FACTOR;
    
} //Signal constructor


 void Signal::SetParameterization(int whichparameterization) {

  WHICHPARAMETERIZATION=whichparameterization;
}


void Signal::TaperVmMHz(double viewangle,
		double deltheta_em,
		double deltheta_had,
		double emfrac,
		double hadfrac,

		double& vmmhz1m,
		double& vmmhz1m_em_obs) {

  //--EM
  
    bool calpulser = false;
  double vmmhz1m_em=0; // V/m/MHz at 1m due to EM component of shower
double vmmhz1m_had=0; // V/m/MHz at 1m due to HAD component of shower

  // this is the number that get exponentiated
  //  double rtemp=0.5*(viewangle-changle)*(viewangle-changle)/(deltheta_em*deltheta_em);
  double rtemp=(viewangle-changle)*(viewangle-changle)/(deltheta_em*deltheta_em);
  

  //cout << "dangle, deltheta_em is " << viewangle-changle << " " << deltheta_em << "\n";
  //cout << "rtemp (em) is " << rtemp << "\n";
  // the power goes like exp(-(theta_v-theta_c)^2/Delta^2)
    // so the e-field is the same with a 1/2 in the exponential
    
    if (emfrac>pow(10.,-10.)) { // if there is an em component
        if (calpulser == false){
            
            if (rtemp<=20) { // if the viewing angle is less than 20 sigma away from the cerankov angle
                // this is the effect of the em width on the signal
                vmmhz1m_em=vmmhz1m*exp(-rtemp);
                
            }
            else // if it's more than 20 sigma just set it to zero 
            {vmmhz1m_em=0.;}
            
        } else {
            vmmhz1m_em = vmmhz1m;
        }
    } else // if the em component is essentially zero than set this to zero
            vmmhz1m_em=0;

  //--HAD
  // this is the quantity that gets exponentiated

  rtemp=(viewangle-changle)*(viewangle-changle)/(deltheta_had*deltheta_had);

  //cout << "rtemp (had) is " << rtemp << "\n";

  if (hadfrac!=0) { // if there is a hadronic fraction
   if (calpulser == false){

    if (rtemp<20) { // if we are less than 20 sigma from cerenkov angle
      vmmhz1m_had=vmmhz1m*exp(-rtemp); // this is the effect of the hadronic width of the signal
   
    }
    else // if we're more than 20 sigma from cerenkov angle
      vmmhz1m_had=0.; // just set it to zero
   } else {
       vmmhz1m_had = vmmhz1m;
   }
  }
  else 
    vmmhz1m_had=0.;


  logscalefactor_taper=log10((emfrac*vmmhz1m_em+hadfrac*vmmhz1m_had)/vmmhz1m);


  //cout << "emfrac, vmmhz1m_em, hadfrac, vmmhz1m_had are " << emfrac << " " << vmmhz1m_em << " " << hadfrac << " " << vmmhz1m_had << "\n";
  vmmhz1m=sin(viewangle)*(emfrac*vmmhz1m_em+hadfrac*vmmhz1m_had);

  if (vmmhz1m==0)
    vmmhz1m_em_obs=0.;
  else
    vmmhz1m_em_obs=sin(viewangle)*emfrac*vmmhz1m_em/vmmhz1m;



} //TaperVmMHz
