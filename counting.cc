//#include "vector.hh"
#include "Vector.h"
//#include "position.hh"
#include "Position.h"
#include "counting.hh"
#include "Tools.h"
#include "TMath.h"
#include "Constants.h"
#include "Event.h"


Counting::Counting() {
  Tools::Zero(npass,2);
  Tools::Zero(npassestrigger,2);
  Tools::Zero(nchanceinhell2,2);
  Tools::Zero(nviewanglecut,2);
  Tools::Zero(nchanceinhell,2);
  Tools::Zero(nchanceinhell_1overr,2);
  Tools::Zero(nchanceinhell_fresnel,2);
  Tools::Zero(nconverges,2);
  Tools::Zero(nacceptablerf,2);
  Tools::Zero(nraywithincontinent1,2); // reality check:  exiting ray is within 30 degrees of south pole
  Tools::Zero(nraywithincontinent2,2); // same, after next iteration.
  Tools::Zero(nraypointsup1,2); // same, after next iteration.
  Tools::Zero(nraypointsup2,2); // same, after next iteration.
  Tools::Zero(nnottoosmall,2); // same, after next iteration.
  Tools::Zero(nviewangle_lt_90,2); // same, after next iteration.
  Tools::Zero(ngoodfracs,2); // same, after next iteration.
  Tools::Zero(nbadfracs,2); // same, after next iteration.
  Tools::Zero(nnottir,2); // same, after next iteration.
  Tools::Zero(nentersice,2); // same, after next iteration.
  Tools::Zero(nabsorbed,2); // same, after next iteration.
  Tools::Zero(noway,2); // same, after next iteration.
  Tools::Zero(wheredoesitleave_err,2); // same, after next iteration.
  Tools::Zero(neverseesice,2); // same, after next iteration.
Tools::Zero(iceinteraction,2); // same, after next iteration.
Tools::Zero(inhorizon,2); // same, after next iteration.
 Tools::Zero(wheredoesitenterice_err,2); // same, after next iteration.
Tools::Zero(toohigh,2); // same, after next iteration.
Tools::Zero(toolow,2); // same, after next iteration.
  for (int i=0;i<NCOSTHETA;i++) {
    Tools::Zero(weights_rin[i],NPHI); // same, after next iteration.
  }
// variables for counting neutrinos and reporting results.
  nnu_e=0;  //counting the number of e,mu,tau neutrinos
  nnu_mu=0;  
  nnu_tau=0;



}

Counting::~Counting() {

}

void Counting::findCosthetaPhiBins(Position r_in,int &icostheta,int &iphi) {

  icostheta=(int)((cos(r_in.Theta())-COSTHETAMIN)/(COSTHETAMAX-COSTHETAMIN)*(double)NCOSTHETA);
  iphi=(int)((r_in.Phi()-PHIMIN)/(PHIMAX-PHIMIN)*(double)NPHI);


}
void Counting::IncrementWeights_r_in(Position r_in,double weight) {
  int iphi,icostheta;
  findCosthetaPhiBins(r_in,icostheta,iphi);
  weights_rin[icostheta][iphi]+=weight;

}


void Counting::findErrorOnSumWeights(double *eventsfound_binned,double &error_plus,double &error_minus) {
      for (int i=0;i<NBINS;i++) { // we are going to sum the error on the weight squared over the weight bins
	// in each bin, the error on the weight is the weight for that bin times the error on the number of events in that bin
	double thislogweight=((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT; // find log weight for this bin
	if (eventsfound_binned[i]<=20) {  // if the number of events in this bin <20, use poisson errors
	  error_plus+=pow(poissonerror_plus[(int)eventsfound_binned[i]]*pow(10.,thislogweight),2);
	  error_minus+=pow(poissonerror_minus[(int)eventsfound_binned[i]]*pow(10.,thislogweight),2);
	}
	else {// otherwise, use sqrt(n) errors
	  error_plus+=eventsfound_binned[i]*pow(pow(10.,thislogweight),2);
	  error_minus=error_plus;

	}

      }
      error_plus=sqrt(error_plus); // take the sqrt of the sum of the squares
      error_minus=sqrt(error_minus);
    }
//void Counting::incrementEventsFound(double weight,Interaction *interaction1) {
void Counting::incrementEventsFound(double weight, Event *event) {

  int index_weights=findWeightBin(log10(weight));


  // count number of events that pass, binned in weight
  if (index_weights<Counting::NBINS)
    eventsfound_binned[index_weights]++;
  
  // incrementing by flavor
  // also bin in weight for error calculation.
  if (event->nuflavor=="nue") { 
    sum[0]+=weight;
    eventsfound_binned_e[index_weights]++;
  }
  if (event->nuflavor=="numu") {
    sum[1]+=weight; // total sum of weights for this flavor
    eventsfound_binned_mu[index_weights]++;
  }
  if (event->nuflavor=="nutau") {
    sum[2]+=weight;
    eventsfound_binned_tau[index_weights]++;
  }
			     
			      
			      
}
int Counting::findWeightBin(double logweight) {
  // first, find which weight bin it is in
  int index_weights;
 if (logweight<MIN_LOGWEIGHT)  // underflows, set to 0th bin
    index_weights=0;
  else if (logweight>MAX_LOGWEIGHT) // overflows, set to last bin
    index_weights=NBINS-1;
  else
    // which index weight corresponds to.
    index_weights=(int)(((logweight-MIN_LOGWEIGHT)/(MAX_LOGWEIGHT-MIN_LOGWEIGHT))*(double)NBINS);  
 return index_weights;
}
