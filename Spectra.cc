#include "Spectra.h"
#include "Tools.h"
#include "Settings.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>

ClassImp(Spectra);


Spectra::Spectra() {
    // default constructor
}


Spectra::~Spectra() {

}

Spectra::Spectra(Settings *settings1) {

  EXPONENT_model = settings1->EXPONENT;
  EXPONENT_min = settings1->EXPONENT_MIN;
  EXPONENT_max = settings1->EXPONENT_MAX;

  // initialize parameters!!

  E_bin = 12;
  const double dE_bin = (EXPONENT_max - EXPONENT_min) / (double)(E_bin - 1);

  double Emuons[E_bin]; // E dN/dE/dA/dt for neutrinos that are produced as muon neutrinos or muon antineutrinos.
  double Eelectrons[E_bin];// E dN/dE/dA/dt for neutrinos that are produced as electron neutrinos or muon antineutrinos.
  
  for (int i = 0; i < E_bin; i++) {
    energy.push_back(EXPONENT_min + ((double)i) * dE_bin);
    Emuons[i] = -30.;
    Eelectrons[i] = -30.;
  } //for

  // end of initialization!!


  if (EXPONENT_model <= 4.)  // dNdEdAdt ~ E^-(EXPONENT_model)
  {
      energy.clear();
      for (int i = 0; i < E_bin; i++) {   // set energy, EdNdEdAdt and E2dNdEdAdt
          energy.push_back(EXPONENT_min + ((double)i) * dE_bin);   // in log, in eV
          lgEdNdEdAdt.push_back(-(EXPONENT_model-1) * (energy[i] - 9.));    // in log, GeV
          lgE2dNdEdAdt.push_back(lgEdNdEdAdt[i] + (energy[i] - 9.));    // in log, in GeV
      }
  }

  else if (EXPONENT_model >= 10. && EXPONENT_model < 30.)
  {
      pnu_EXPONENT = EXPONENT_model;
      return;
  }
  
  else if (EXPONENT_model >= 510. && EXPONENT_model <= 650.) // g.n. added so we can have simulations at E^2=17.8 not just whole numbers
  {
      pnu_EXPONENT = (EXPONENT_model - 400) / 10;
      if (settings1->EVENT_GENERATION_MODE != 1){
	      cout<<"**************** energy is "<<pnu_EXPONENT<<" *******************"<<endl;
      }
      return;
  }
  
  else if (EXPONENT_model == 30.) // ESS baseline model. Used to be EXPONENT "0"
  {
    E_bin = 9;

    // electron component of Figure 4 of ES&S
    // astro-ph/0101216
    Eelectrons[0]=-17.2; // 16.
    Eelectrons[1]=-17.35; // 16.5
    Eelectrons[2]=-17.2; // 17.
    Eelectrons[3]=-17.1; // 17.5
    Eelectrons[4]=-17.2; // 18.
    Eelectrons[5]=-17.5; // 18.5
    Eelectrons[6]=-18.0; // 19
    Eelectrons[7]=-18.5; // 19.5
    Eelectrons[8]=-19.4; // 20.
    Eelectrons[9]=-30.; // 20.5 punt------ not using above here
    Eelectrons[10]=-30.; // 21.0 punt
    Eelectrons[11]=-30.; // 21.5 punt
  
    // lower curve of Figure 9 of ES&S
    // astro-ph/0101216
    Emuons[0]=-17.1;  //16.
    Emuons[1]=-16.6;  //16.5
    Emuons[2]=-16.3;  //17.
    Emuons[3]=-16.2; // 17.5
    Emuons[4]=-16.4; // 18.
    Emuons[5]=-16.7; // 18.5
    Emuons[6]=-17.3; // 19
    Emuons[7]=-17.95; // 19.5
    Emuons[8]=-18.85; // 20.
    Emuons[9]=-19.9; // 20.5 punt------not using above here
    Emuons[10]=-30.; // 21.0 punt
    Emuons[11]=-30.; // 21.5 punt
  
    energy.clear();
    for (int i=0;i<E_bin;i++) {
        energy.push_back(16.+((double)i)/2.);   // in log, in eV
        lgEdNdEdAdt.push_back(log10( pow(10., Eelectrons[i]) + pow(10., Emuons[i]) ));   // in log
        // I know that it have to change this to non-log later but I want to make them same way.
        lgE2dNdEdAdt.push_back(lgEdNdEdAdt[i] + (energy[i] - 9.));    // in log, in GeV
    }
  }

  else if (EXPONENT_model == 31.) // ESS-cosmological constant. Use Nu_el : Nu_mu ratio. used to be EXPONENT "5"
  {
    E_bin = 9;
    double emuratio[E_bin];

    // electron component of Figure 4 of ES&S
    // astro-ph/0101216
    Eelectrons[0]=-17.2; // 16.
    Eelectrons[1]=-17.35; // 16.5
    Eelectrons[2]=-17.2; // 17.
    Eelectrons[3]=-17.1; // 17.5
    Eelectrons[4]=-17.2; // 18.
    Eelectrons[5]=-17.5; // 18.5
    Eelectrons[6]=-18.0; // 19
    Eelectrons[7]=-18.5; // 19.5
    Eelectrons[8]=-19.4; // 20.
    Eelectrons[9]=-30.; // 20.5 punt------ not using above here
    Eelectrons[10]=-30.; // 21.0 punt
    Eelectrons[11]=-30.; // 21.5 punt
    
    // muon component of Figure 4 of ES&S
    // astro-ph/0101216
    Emuons[0]=-17.8; // 16.
    Emuons[1]=-17.4; // 16.5
    Emuons[2]=-17.; // 17.
    Emuons[3]=-16.75; // 17.5
    Emuons[4]=-16.9; // 18.
    Emuons[5]=-17.2; // 18.5
    Emuons[6]=-17.7; // 19
    Emuons[7]=-18.3; // 19.5
    Emuons[8]=-19.1; // 20.
    Emuons[9]=-30.; // 20.5 punt
    Emuons[10]=-30.; // 21.0 punt
    Emuons[11]=-30.; // 21.5 punt
    
    for(int i=0;i<E_bin;i++)
      emuratio[i]=Eelectrons[i]/Emuons[i];
    
    // upper curve in Figure 9 of ES&S
    // astro-ph/0101216
    Emuons[0]=-16.85;  //16.
    Emuons[1]=-16.4;  //16.5
    Emuons[2]=-16.05;  //17.
    Emuons[3]=-16.; // 17.5
    Emuons[4]=-16.15; // 18.
    Emuons[5]=-16.5; // 18.5
    Emuons[6]=-17.1; // 19
    Emuons[7]=-17.7; // 19.5
    Emuons[8]=-18.65; // 20.
    Emuons[9]=-19.75; // 20.5 punt
    Emuons[10]=-30.; // 21.0 punt
    Emuons[11]=-30.; // 21.5 punt
    
    
    for(int i=0;i<E_bin;i++) {
      Eelectrons[i]=Emuons[i]*emuratio[i];
      cout << "Eelectrons, Emuons are " << Eelectrons[i] << " " << Emuons[i] << "\n";
    }
  
    energy.clear();
    for (int i=0;i<E_bin;i++) {
      energy.push_back(16.+((double)i)/2.);   // in log, in eV
      lgEdNdEdAdt.push_back(log10( pow(10.,Eelectrons[i]) + pow(10.,Emuons[i]) ));  // in log.
      lgE2dNdEdAdt.push_back(lgEdNdEdAdt[i] + (energy[i] - 9.));    // in log, in GeV
      cout << "lgEdNdEdAdt are " << lgEdNdEdAdt[i] << "\n";
    }
    
  } // end if ESS-cosmological constant


  else if (EXPONENT_model > 31. && EXPONENT_model < 200.)  // use digitized flux from different models
  {
      switch ( (int)EXPONENT_model )
      {
          case 32:  // ESS3
              GetFlux("essfig9.dat");
              break;
          case 33:  // ESS4
              GetFlux("essbaseline.dat");
              break;
          case 34:  // ESS5
              GetFlux("ess_n0.dat");
              break;
          case 40:
              GetFlux("ahlers.dat");
              break;
          case 50:
              GetFlux("allard.dat");
              break;
          case 60:
              GetFlux("ave_max.dat");
              break;
          case 61:
              GetFlux("ave_min.dat");
              break;
          case 70:
              GetFlux("kkss_envo.dat");
              break;
          case 80:
              GetFlux("gzk_peter.dat");
              break;
          case 90:
              GetFlux("waxgzk.dat");
              break;
          case 100:
              GetFlux("e-2.dat");
              break;
          case 110:
              GetFlux("yuksel_grb.dat");
              break;
          case 111:
              GetFlux("yuksel_qso.dat");
              break;
          case 112:
              GetFlux("yuksel_sfh.dat");
              break;
          case 113: // equal number of triggers per lnE based on A23 2020 trigger-level effective area
              GetFlux("equalTriggersPerLnE.dat");
              break;
//          case 100:
//              GetFlux("berezinsky_saturate.dat");
//              break;
          default:
              cout<<"Error: Wrong input of EXPONENT!"<<endl;
      }
  }

  else if (EXPONENT_model == 200.)   // Iron model
  {
      GetFlux("Ave2005_Fe_Emax21.0.dat");
  }

  else if (EXPONENT_model == 201.)
  {
      GetFlux("Ave2005_Fe_Emax21.5.dat");
  }

  else if (EXPONENT_model == 202.)
  {
      GetFlux("Ave2005_Fe_Emax22.0.dat");
  }

  else if (EXPONENT_model == 203.)
  {
      GetFlux("Ave2005_Fe_hi_evo.dat");
  }

  else if (EXPONENT_model == 204.)
  {
      GetFlux("Ave2005_Fe_low_evo.dat");
  }

  else if (EXPONENT_model == 210.)
  {
      GetFlux("Stanev2008_heavy.dat");
  }

  else if (EXPONENT_model == 220.)
  {
      GetFlux("Kotera2010_Fe_only.dat");
  }

  else if (EXPONENT_model == 221.)
  {
      GetFlux("Kotera2010_Fe_rich.dat");
  }

  else if (EXPONENT_model == 222.)
  {
      GetFlux("Kotera2010_mix_max.dat");
  }

  else if (EXPONENT_model == 223.)
  {
      GetFlux("Kotera2010_mix_min.dat");
  }

  else if (EXPONENT_model == 224.)
  {
      GetFlux("Kotera2010_max.dat");
  }

  else if (EXPONENT_model == 225.)
  {
      GetFlux("CenA_Kachelriess.dat");
  }


  //
  // End of selecting the Model!!!
  //

  // From log to linear!!  
  for (int i = 0; i < E_bin; i++) {
      EdNdEdAdt.push_back(pow(10, lgEdNdEdAdt[i]));     //to linear 
      E2dNdEdAdt.push_back(pow(10, lgE2dNdEdAdt[i]));   //to linear
  }
  
  maxflux = Tools::dMax(&EdNdEdAdt[0], E_bin);

}
 

double  Spectra::GetNuEnergy() {
 
  if ((EXPONENT_model >= 10. && EXPONENT_model < 30.) || (EXPONENT_model >= 510. && EXPONENT_model <= 650.)) {
      return pow(10., pnu_EXPONENT);
  }
  
  double thisenergy = EXPONENT_min; // arbitrary initialisation
  double thisflux = 2.; // initialise higher than max
  double max = 1.;
  int energybin = 0; // arbitrary initialisation
  double maxenergy = Tools::dMax(&energy[0], E_bin);
  double minenergy = Tools::dMin(&energy[0], E_bin); 
  
  // this uses the dartboard approach
  while(thisflux > max) {
    // pick an energy  
    thisenergy = Rand3.Rndm() * (maxenergy - minenergy) + minenergy; // pick energy at random between the highest and lowest
    // the energy array is actually filled with energy exponents 
    
    // get interpolated flux at "thisenergy"
    max = pow(10, SimpleLinearInterpolation_value(E_bin, &energy[0], &lgEdNdEdAdt[0], thisenergy)) / maxflux; // normalize to 1
    thisflux = Rand3.Rndm(); // pick the flux at random between 0 and 1, if it's less than max it's a winner
  } //while
  return pow(10., thisenergy);
	
} //Pick Neutrino Energy

double Spectra::SimpleLinearInterpolation_value(int n1, double *x1, double *y1, double x2 ) {    // reads n1 array values x1, y1 and do simple linear interpolated value at x2 and return it
    //

    if ( x2 <= x1[0] ) 
      return y1[0];

    if ( x2 >= x1[n1-1] ) 
      return y1[n1-1];

    int cnt = 0;
    double output;

    // find the nearest bin for x2 value
    for ( int i=0; i<n1; i++) {

        if ( x1[i] >= x2 ) {

            cnt = i;
            output = y1[cnt-1] + (x2-x1[cnt-1])*(y1[cnt]-y1[cnt-1])/(x1[cnt]-x1[cnt-1]);

            break; // we found the value, get out from the loop
        }
    }

    return output;

}

inline void Spectra::GetFlux(string filename)
{
    ifstream influx(("./fluxes/"+filename).c_str());
    int NLINES;
    influx >> NLINES;   // Read how much lines in the file.
    cout<<"We are using "<<filename.c_str()<<" as the flux data."<<endl;
    cout<<"total lines in the file are "<<NLINES<<endl;
    E_bin = NLINES;
    
    energy.clear();
    for(int i=0;i<NLINES;i++) {
        double ebuff, fluxbuff;
        influx >> ebuff >> fluxbuff;
        energy.push_back(ebuff);
        lgE2dNdEdAdt.push_back(fluxbuff);
    }
    
    for(int i=0;i<E_bin;i++) {
        lgEdNdEdAdt.push_back(lgE2dNdEdAdt[i] + 9. - energy[i]);  // change from GeV to eV and E2dN -> EdN
    }
}

double *Spectra::Getenergy() {
    return &energy[0];
}

double *Spectra::GetEdNdEdAdt()  {
    return &EdNdEdAdt[0];
}

double *Spectra::GetE2dNdEdAdt() {
    return &E2dNdEdAdt[0];
}

double Spectra::GetEdNdEdAdt(double E_val) {
    
    double tmp_Get;
    if (E_val < energy[0]) {
        cout<<"Energy value is smaller than the energy boundary!\n";
        cout<<"Energy value is replaced to minimum value of energy bound : "<<energy[0]<<"\n";
        tmp_Get = EdNdEdAdt[0];
    }
    else if (E_val > energy[E_bin-1]) {
        cout<<"Energy value is bigger than the energy boundary!\n";
        cout<<"Energy value is replaced to maximum value of energy bound : "<<energy[E_bin-1]<<"\n";
        tmp_Get = EdNdEdAdt[E_bin-1];
    }
    else {
        tmp_Get = EdNdEdAdt[ Tools::Getifreq( E_val, energy[0], energy[E_bin-1], E_bin ) ];
    }
    return tmp_Get;
}

double Spectra::GetE2dNdEdAdt(double E_val) {
    double tmp_Get;
    if (E_val < energy[0]) {
        cout<<"Energy value is smaller than the energy boundary!\n";
        cout<<"Energy value is replaced to minimum value of energy bound : "<<energy[0]<<"\n";
        tmp_Get = E2dNdEdAdt[0];
    }
    else if (E_val > energy[E_bin-1]) {
        cout<<"Energy value is bigger than the energy boundary!\n";
        cout<<"Energy value is replaced to maximum value of energy bound : "<<energy[E_bin-1]<<"\n";
        tmp_Get = E2dNdEdAdt[E_bin-1];
    }
    else {
        //tmp_Get = sE2dNdEdAdt->Eval(E_val);
        tmp_Get = E2dNdEdAdt[ Tools::Getifreq( E_val, energy[0], energy[E_bin-1], E_bin ) ];
    }
    return tmp_Get;
}

double Spectra::Getmaxflux() {
    return maxflux;
}


int Spectra::GetE_bin() {
    return E_bin;
}


int Spectra::IsSpectrum() {
    int out;
    if (EXPONENT_model>=30. && EXPONENT_model <510) {
        out = 1;
    }
    else {
        out = 0;
    }
    return out;
}

int Spectra::IsMonoenergetic() {
    int out;
    if ((EXPONENT_model>0.&&EXPONENT_model<30.) || (EXPONENT_model>=510 && EXPONENT_model<=650) ) {
        out = 1;
    }
    else {
        out = 0;
    }
    return out;
}
