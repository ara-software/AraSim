#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include "TTreeIndex.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h" 
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
//#include "TObject.h"

//#include <fftw3.h>

using namespace std;

#include "Tools.h"
#include "Constants.h"
#include "Vector.h"
#include "Position.h"
#include "EarthModel.h"
#include "IceModel.h"
#include "Efficiencies.h"
#include "Spectra.h"
#include "Event.h"
#include "Trigger.h"
#include "Detector.h"
#include "Settings.h"
#include "counting.hh"
#include "Primaries.h"
#include "signal.hh"
#include "secondaries.hh"

#include "Ray.h"
#include "RaySolver.h"
#include "Report.h"







string outputdir="outputs";




//int main() {
int main(int argc, char **argv) {    // this is for manual power threshold value

RaySolver *raysolver = new RaySolver();


double src_x = atof(argv[1]);
double src_y = atof(argv[2]);
double src_z = atof(argv[3]);
cout<<"src_x : "<<src_x<<" src_y : "<<src_y<<" src_z : "<<src_z<<endl;

double trg_x = atof(argv[4]);
double trg_y = atof(argv[5]);
double trg_z = atof(argv[6]);
cout<<"trg_x : "<<trg_x<<" trg_y : "<<trg_y<<" trg_z : "<<trg_z<<endl;



raysolver->Solve_Ray_org ( src_x, src_y, src_z, trg_x, trg_y, trg_z, false );
//raysolver->Solve_Ray_org ( src_x, src_y, src_z, trg_x, trg_y, trg_z, true );

}



