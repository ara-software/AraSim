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
#include "Report.h"

#include "signal.hh"


#include "Ray.h"

//#include "FFTtools.h"

class EarthModel; //class


Detector *detector = 0; 
Settings *settings = 0;
Spectra *spectra = 0;
IceModel *icemodel = 0;
Event *event = 0;
Report *report = 0;
Trigger *trigger = 0;



string outputdir="outputs";




string getInfo_fromFilename (const char *runfile);

string getInfo_fromFilename_oldfiles (const char *runfile);

string getInfo_fromFilename_FdomainWoBug (const char *runfile);

string getInfo_fromFilename_TdomainWBug (const char *runfile);


//int main() {
int main(int argc, char **argv) {    // this is for manual power threshold value


    string readfile;
  if (argc<2) { // no setup file input, use default
      readfile = "outputs/AraOut.root";
  }
  else if (argc == 2) { // read file!!
      readfile = string( argv[1] );
  }
  else { // no mode for argc > 2!
      cout<<"reading many files!"<<endl;
      //cout<<"too many info! just use default AraOut.root file!"<<endl;
      //readfile = "outputs/AraOut.root";
  }



  // mean RMS values from Testbed data (which is calibrated for AraSim)
  double RMSV_ch[8];
  RMSV_ch[0] = 38.2552;
  RMSV_ch[1] = 49.7648;
  RMSV_ch[2] = 58.9413;
  RMSV_ch[3] = 44.0512;
  RMSV_ch[4] = 33.2493;
  RMSV_ch[5] = 48.6344;
  RMSV_ch[6] = 40.7179;
  RMSV_ch[7] = 36.3631;






//  Settings *settings = new Settings();

  //  Detector *detector=new Detector(settings->DETECTOR); // builds antenna array, 0 for testbed
//  Detector *detector=0; // builds antenna array, 0 for testbed
/*
  Detector *detector = 0; 
  Settings *settings = 0;
  Spectra *spectra = 0;
  IceModel *icemodel = 0;
  Event *event = 0;
  Report *report = 0;
  Trigger *trigger = 0;

*/

  //cout<<"construct detector"<<endl;

  Signal *signal = 0;

  Signal *signal_2 = 0;


  // global pass evt count
  int pass_evts = 0;

  //int string, ant;
  int ch_from_det;




  stringstream sstr;



  int total_evt;
  int passed_evt_run;
  double total_passed_evt = 0.;

  int passed_evt; // flag for passed evt or not


  double total_weight = 0.;

  double total_weight_el = 0.;
  double total_weight_mu = 0.;
  double total_weight_tau = 0.;



  double total_probability = 0.; // weight * probability

  double total_NNU = 0.;

  double Vice;
  double Vice_pre;

  double Aice;

  double Veff_we; // water equivalent
  double Veff;
  double Aeff;

  double len_int_kgm2_total; // interaction length in unit of km/m^2

  double len_int_m; // interaction length in the ice in unit of m

  int interaction_mode; // depending on interaction mode, we have to decide which way to calculate effective area Aeff
  // 1 : pick near = get Veff and then get Aeff from Veff
  // 0 : pick unbiased = get Aeff directly




  // test distribution map of PNU

  TH1D *hpnu = new TH1D ("hpnu", "", 50, 15, 22 );


  // get error 
  //
  Counting *count1 = new Counting();
  cout<<"called Counting"<<endl;




  // plot for V saturation
  TH1D *hPeakVpeak = new TH1D("hPeakVpeak", "", 100, 0. , 100. );





  // plots for POSNU RADIUS
  readfile = string( argv[1] ); // read first file
  TFile *AraFile_tmp=new TFile(( readfile ).c_str());
  TTree *AraTree_tmp=(TTree*)AraFile_tmp->Get("AraTree");
  AraTree_tmp->SetBranchAddress("settings",&settings);
  AraTree_tmp->SetBranchAddress("detector",&detector);
  AraTree_tmp->GetEvent(0);

  //double posnu_depth = 3000.;   // first estimated value, 3km depth
  double posnu_depth = 5000.;   // test 5km depth
  //double posnu_xmax = detector->params.core_x + settings->POSNU_RADIUS;
  //double posnu_xmin = detector->params.core_x - settings->POSNU_RADIUS;
  double posnu_xmax = detector->params.core_x + settings->POSNU_RADIUS + 1000.;
  double posnu_xmin = detector->params.core_x - settings->POSNU_RADIUS - 1000.;

  //double posnu_ymax = detector->params.core_y + settings->POSNU_RADIUS;
  //double posnu_ymin = detector->params.core_y - settings->POSNU_RADIUS;
  double posnu_ymax = detector->params.core_y + settings->POSNU_RADIUS + 1000.;
  double posnu_ymin = detector->params.core_y - settings->POSNU_RADIUS - 1000.;

  int nx, ny, nz;
  //nx = (int)( (posnu_xmax - posnu_xmin)/10 );   // 10m bin size
  //ny = (int)( (posnu_ymax - posnu_ymin)/10 );   // 10m bin size
  //nz = (int)( ( posnu_depth)/10 );   // 10m bin size
  nx = (int)( (posnu_xmax - posnu_xmin)/100 );   // 100m bin size
  ny = (int)( (posnu_ymax - posnu_ymin)/100 );   // 100m bin size
  nz = (int)( ( posnu_depth)/100 );   // 100m bin size

  cout<<"nx : "<<nx<<" ny : "<<ny<<" nz : "<<nz<<endl;

  cout<<"posnu_xmin : "<<posnu_xmin<<endl;
  cout<<"posnu_xmax : "<<posnu_xmax<<endl;
  cout<<"posnu_ymin : "<<posnu_ymin<<endl;
  cout<<"posnu_ymax : "<<posnu_ymax<<endl;
  cout<<"-posnu_depth : "<<-posnu_depth<<endl;

  /*
  TH2D *hposnu_xy = new TH2D ("hposnu_xy","",nx, posnu_xmin, posnu_xmax, ny, posnu_ymin, posnu_ymax);
  TH2D *hposnu_xz = new TH2D ("hposnu_xz","",nx, posnu_xmin, posnu_xmax, nz, -posnu_depth, 0.);
  TH2D *hposnu_rz = new TH2D ("hposnu_rz","",nx, 0, settings->POSNU_RADIUS + 1000., nz, -posnu_depth, 0.);

  TH1D *hposnu_ice = new TH1D ("hposnu_ice","", nz, 0., 5000.); // check ice depth
  */

  TH2D *hposnu_xy = new TH2D ("hposnu_xy","",100, posnu_xmin, posnu_xmax, 100, posnu_ymin, posnu_ymax);
  TH2D *hposnu_xz = new TH2D ("hposnu_xz","",100, posnu_xmin, posnu_xmax, 50, -posnu_depth, 0.);
  TH2D *hposnu_rz = new TH2D ("hposnu_rz","",100, 0, settings->POSNU_RADIUS + 1000., 50, -posnu_depth, 0.);

  TH1D *hposnu_ice = new TH1D ("hposnu_ice","", 50, 0., 5000.); // check ice depth


  // raysol r 
  double raysolR = settings->RAYSOL_RANGE;


  // plots for RAYSOL_RANGE
  TH1D *hposnuL = new TH1D ("hposnuL","", 100, 0., raysolR + 1000.); // check trig chs' distance between antenna and posnu (additional 1km)





  // plots for Offcone 
  TH1D *hoffcone = new TH1D ("hoffcone","", 60, 0., 30.); // check trig chs' strongest Vpeak chs' offcone angle



  // plot for weird ray tracing result channels
  TH1D *hweirdRayTrace = new TH1D ("hweirdRayTrace", "", 16, 0, 16 ); // from string0 to string3 (ant0 to ant4)





  int off011evts = 0;


  //string fileinfo = getInfo_fromFilename ( argv[1] );
  //string fileinfo = getInfo_fromFilename_oldfiles ( argv[1] );
  //string fileinfo = getInfo_fromFilename_FdomainWoBug ( argv[1] );
  string fileinfo = getInfo_fromFilename_TdomainWBug ( argv[1] );



  
  int string_no, ant_no;


  int openwell; // flag if file open well


  for (int run=1; run<argc; run++) {


      cout<<"run : "<<run<<" out of "<< argc<<endl;
      cout<<"run file : "<< argv[run]<<endl;

      openwell = 1;
  
  
      TFile *AraFile=TFile::Open( argv[run] );
      //TFile *AraFile=new TFile(( readfile ).c_str());
      //TFile *AraFile=new TFile((outputdir+"/AraOut.root").c_str());
      if (!AraFile) {
          cerr<<"Can't open file\n";
          return -1;
          openwell = 0;
      }
      cout<<"AraFile"<<endl;


      TTree *AraTree=(TTree*)AraFile->Get("AraTree");
      TTree *AraTree2=(TTree*)AraFile->Get("AraTree2");
      cout<<"AraTree"<<endl;
      AraTree->SetBranchAddress("detector",&detector);
      AraTree->SetBranchAddress("settings",&settings);
      AraTree->SetBranchAddress("spectra",&spectra);
      AraTree->SetBranchAddress("icemodel",&icemodel);
      AraTree->SetBranchAddress("trigger",&trigger);
      AraTree2->SetBranchAddress("event",&event);
      AraTree2->SetBranchAddress("report",&report);
      cout<<"branch detector"<<endl;
      
      AraTree->GetEvent(0);


      total_NNU += settings->NNU;

      total_evt = AraTree2->GetEntries();
      cout<<"AraTree2 entries : "<<total_evt<<endl;


      interaction_mode = settings->INTERACTION_MODE;
  
      passed_evt_run = 0;


      for (int inu=0;inu<total_evt;inu++) { // loop over neutrinos


      
          AraTree2->GetEvent(inu);

          passed_evt = 0;

  




           
          for (int i=0; i<detector->params.number_of_stations; i++) {
              
              if ( report->stations[i].Global_Pass ) {

                  passed_evt = 1;
                  break;
              }
              
          }




                      
          double repro_d_theta_em, repro_d_theta_had;

          double repro_d_theta_em_2, repro_d_theta_had_2;


          // collect weight values for passed events
          if ( passed_evt == 1 ) {



                  
              signal = new Signal (settings);
              signal->SetMedium(0);   // set medium as ice
              signal->SetNDepth( icemodel->GetN( event->Nu_Interaction[0].posnu) );


              double peakv_tmp = 0.;

              // loop over all stations, antennas and find the peak of Vpeak

              for (int i=0; i<(int)detector->params.number_of_stations; i++) {

                  //for (int j=0; j<(int)detector->params.number_of_strings_station; j++) {
                  for (int j=0; j<(int)detector->stations[i].strings.size(); j++) {
                        
                      //for (int k=0; k<(int)detector->params.number_of_antennas_string; k++) {
                      for (int k=0; k<(int)detector->stations[i].strings[j].antennas.size(); k++) {


                          if ( report->stations[i].strings[j].antennas[k].Trig_Pass > 0 ) {

                              double minoffcone = 100;

                              for (int l=0; l<(int)report->stations[i].strings[j].antennas[k].ray_sol_cnt; l++) {

                                  if ( report->stations[i].strings[j].antennas[k].PeakV[l] > peakv_tmp ) {

                                      peakv_tmp = report->stations[i].strings[j].antennas[k].PeakV[l];
                                  }


                                  // check the min offcone angle btw raysols
                                  if ( minoffcone > fabs(report->stations[i].strings[j].antennas[k].view_ang[l] - signal->changle)*DEGRAD ) {
                                  
                                      minoffcone = fabs(report->stations[i].strings[j].antennas[k].view_ang[l] - signal->changle)*DEGRAD;
                                      //cout<<"changle : "<<signal->changle*DEGRAD<<endl;
                                  }

                              }

                              // fill in offcone angle
                              hoffcone -> Fill( minoffcone, event->Nu_Interaction[0].weight );
                                      
                                      
                              // fill in posnuL
                              double dist = detector->stations[i].strings[j].antennas[k].Distance( event->Nu_Interaction[0].posnu );
                      
                              hposnuL -> Fill ( dist, event->Nu_Interaction[0].weight );

                              
                              bool RTbug = false;
                              // test if raysol dist is shorter than physical distance
                              for (int l=0; l<(int)report->stations[i].strings[j].antennas[k].ray_sol_cnt; l++) {


                                  if ( dist - report->stations[i].strings[j].antennas[k].Dist[l] >= 10. ){ // 10m or more difference!
                                      cout<<"string["<<j<<"].ant["<<k<<"].raysol["<<l<<"] ray sol dist ("<<report->stations[i].strings[j].antennas[k].Dist[l]<<") is shorter than the physical distance ("<<dist<<")!!"<<endl;
                                      RTbug = true;
                                  }
                              }


                              if ( RTbug == true ) {

                                  hweirdRayTrace->Fill( j*4 + k );
                                  //hweirdRayTrace = new TH1D ("hweirdRayTrace", "", 16, 0, 16 ); // from string0 to string3 (ant0 to ant4)
                              }


                          }
                      }
                  }
              }


              if ( peakv_tmp > 0. ) { // if peakv found

                  hPeakVpeak->Fill(peakv_tmp, event->Nu_Interaction[0].weight );
              }


















              hposnu_xy->Fill( event->Nu_Interaction[0].posnu.GetX(), event->Nu_Interaction[0].posnu.GetY(), event->Nu_Interaction[0].weight );
              hposnu_xz->Fill( event->Nu_Interaction[0].posnu.GetX(), event->Nu_Interaction[0].posnu.GetZ() - icemodel->Surface(0.,0.), event->Nu_Interaction[0].weight );

              double dx = detector->params.core_x - event->Nu_Interaction[0].posnu.GetX();
              double dy = detector->params.core_y - event->Nu_Interaction[0].posnu.GetY();

              hposnu_rz->Fill( sqrt(dx*dx + dy*dy)  , event->Nu_Interaction[0].posnu.GetZ() - icemodel->Surface(0.,0.), event->Nu_Interaction[0].weight );


              hposnu_ice->Fill( icemodel->IceThickness( event->Nu_Interaction[0].posnu ) );





              passed_evt_run++;

              total_passed_evt++;

          }

      } // loop over evts


      cout<<"Passed evts in this file : "<<passed_evt_run<<endl;




      // clear
      AraTree->ResetBranchAddresses();
      AraTree2->ResetBranchAddresses();
      AraFile->Close();
      delete AraFile;



  }// loop over files


      
  cout<<"Total passed evts : "<<total_passed_evt<<endl;
  if ( interaction_mode == 1 ) cout<<"Total weight : "<<total_weight<<endl;
  else if ( interaction_mode == 0 ) cout<<"Total probability : "<<total_probability<<endl;
  cout<<"Total NNU : "<<total_NNU<<endl;

  /*
  // fraction for each flavor
  cout<<"trig flavors"<<endl;
  cout<<"el. : "<<total_weight_el/total_weight<<endl;
  cout<<"mu. : "<<total_weight_mu/total_weight<<endl;
  cout<<"tau. : "<<total_weight_tau/total_weight<<endl;
  */






  // make integrated saturation study plots
  TCanvas *cSat = new TCanvas ("cSat","", 2400,1200);
  cSat->Divide(3,2);

  cSat->cd(1); // Vpeak
  cSat->cd(1)->SetLogx();

  hPeakVpeak->SetTitle("Peak of Vpeak from trig chs");
  hPeakVpeak->SetXTitle("Vpeak (V)");
  hPeakVpeak->SetYTitle("Weight");
  hPeakVpeak->SetMinimum(1.e-2);

  hPeakVpeak->Draw();


  cSat->cd(2); // raysol 

  hposnuL->SetTitle("Trig Chs' dist. btw Ant & posnu");
  hposnuL->SetXTitle("Distance (m)");
  hposnuL->SetYTitle("Weight");
  hposnuL->SetMinimum(1.e-2);

  hposnuL->Draw();


  cSat->cd(3); // offcone

  hoffcone->SetTitle("Offcone angle (abs)");
  hoffcone->SetXTitle("angle (deg)");
  hoffcone->SetYTitle("Weight");
  hoffcone->SetMinimum(1.e-2);

  hoffcone->Draw();


cSat->cd(4); // posnu xy
  hposnu_xy->SetTitle("Posnu XY plane");
  hposnu_xy->SetXTitle("X (m)");
  hposnu_xy->SetYTitle("Y (m)");
  hposnu_xy->Draw("COLZ");

cSat->cd(5); // xz
  hposnu_xz->SetTitle("Posnu XZ plane");
  hposnu_xz->SetXTitle("X (m)");
  hposnu_xz->SetYTitle("Z (m)");
  hposnu_xz->Draw("COLZ");

cSat->cd(6); // rho z
  hposnu_rz->SetTitle("Posnu RhoZ plane");
  hposnu_rz->SetXTitle("Rho (m)");
  hposnu_rz->SetYTitle("Z (m)");
  hposnu_rz->Draw("COLZ");


  sstr.str("");
  //sstr<<"./outputs/plots/hEMshower.pdf";
  //sstr<<"./outputs/plots/hPeakVpeak_"<<fileinfo<<".pdf";
  sstr<<"./outputs/plots/SatStudy/hSaturations_"<<fileinfo<<".pdf";
  cSat->SaveAs(sstr.str().c_str());







  /*
  // make Vpeak saturation study
  TCanvas *cPeakVpeak = new TCanvas ("cPeakVpeak","", 800,600);

  cPeakVpeak->cd();

  hPeakVpeak->SetTitle("Peak of Vpeak from trig chs");
  hPeakVpeak->SetXTitle("Vpeak (V)");
  hPeakVpeak->SetYTitle("Weight");
  hPeakVpeak->SetMinimum(1.e-2);

  hPeakVpeak->Draw();


  sstr.str("");
  //sstr<<"./outputs/plots/hEMshower.pdf";
  //sstr<<"./outputs/plots/hPeakVpeak_"<<fileinfo<<".pdf";
  sstr<<"./outputs/plots/SatStudy/hPeakVpeak_"<<fileinfo<<".pdf";
  cPeakVpeak->SaveAs(sstr.str().c_str());




  // make Raysol range saturation study
  TCanvas *cposnuL = new TCanvas ("cposnuL","", 800,600);

  cposnuL->cd();

  hposnuL->SetTitle("Trig Chs' dist. btw Ant & posnu");
  hposnuL->SetXTitle("Distance (m)");
  hposnuL->SetYTitle("Weight");
  hposnuL->SetMinimum(1.e-2);

  hposnuL->Draw();


  sstr.str("");
  //sstr<<"./outputs/plots/hEMshower.pdf";
  //sstr<<"./outputs/plots/hPeakVpeak_"<<fileinfo<<".pdf";
  sstr<<"./outputs/plots/SatStudy/hPOSNU_L_"<<fileinfo<<".pdf";
  cposnuL->SaveAs(sstr.str().c_str());




  // offcone angle saturation study
  TCanvas *coffcone = new TCanvas ("coffcone","", 800,600);

  coffcone->cd();

  hoffcone->SetTitle("Offcone angle (abs)");
  hoffcone->SetXTitle("angle (deg)");
  hoffcone->SetYTitle("Weight");
  hoffcone->SetMinimum(1.e-2);

  hoffcone->Draw();


  sstr.str("");
  //sstr<<"./outputs/plots/hEMshower.pdf";
  //sstr<<"./outputs/plots/hPeakVpeak_"<<fileinfo<<".pdf";
  sstr<<"./outputs/plots/SatStudy/hOffCone_"<<fileinfo<<".pdf";
  coffcone->SaveAs(sstr.str().c_str());





  // posnu r saturation study
TCanvas *cPosnu = new TCanvas("cPosnu","",4000,1000);
cPosnu->Divide(4,1);

cPosnu->cd(1);
  hposnu_xy->SetTitle("Posnu XY plane");
  hposnu_xy->SetXTitle("X (m)");
  hposnu_xy->SetYTitle("Y (m)");
  hposnu_xy->Draw("COLZ");

cPosnu->cd(2);
  hposnu_xz->SetTitle("Posnu XZ plane");
  hposnu_xz->SetXTitle("X (m)");
  hposnu_xz->SetYTitle("Z (m)");
  hposnu_xz->Draw("COLZ");

cPosnu->cd(3);
  hposnu_rz->SetTitle("Posnu RhoZ plane");
  hposnu_rz->SetXTitle("Rho (m)");
  hposnu_rz->SetYTitle("Z (m)");
  hposnu_rz->Draw("COLZ");

cPosnu->cd(4);
  hposnu_ice->SetTitle("Ice Thickness at Posnu");
  hposnu_ice->SetXTitle("Ice Thickness(m)");
  hposnu_ice->SetYTitle("evts");
  hposnu_ice->Draw();

//cPosnu->Print("POSNU_Hist.pdf");
    sstr.str("");
    //sstr<<"./outputs/plots/h2_E_SNR.pdf";
    //sstr<<"POSNU_Hist.pdf";
    sstr<<"./outputs/plots/SatStudy/POSNU_Hist_"<<fileinfo<<".pdf";
    cPosnu->SaveAs(sstr.str().c_str());

    */





  TCanvas *cRTbug = new TCanvas ("cRTbug","", 800,600);

  cRTbug->cd();
  cRTbug->cd()->SetLogy();

  hweirdRayTrace->SetTitle("Weird RayTracing Chs");
  hweirdRayTrace->SetXTitle("Ch number");
  hweirdRayTrace->SetYTitle("Evt");
  //hweirdRayTrace->SetMinimum(1.e-2);

  hweirdRayTrace->Draw();

  sstr.str("");
  sstr<<"./outputs/plots/SatStudy/hRTbug_"<<fileinfo<<".pdf";
  cRTbug->SaveAs(sstr.str().c_str());






  return 0;
}





string getInfo_fromFilename (const char *runfile) {
    
    string file = string( runfile );
    string chRun = "Fdomain_TestBed_";
    string chRun2 = ".run";
    size_t foundRun=file.find(chRun);
    size_t foundRun2=file.find(chRun2);
    if ( foundRun < 10 ) {
    
        chRun = "Tdomain_TestBed_";
        foundRun=file.find(chRun);
    }

    cout<<"file name first : "<<foundRun<<", second : "<<foundRun2<<endl;
    int diff = foundRun2 - foundRun;
    //string strRunNum = file.substr (foundRun + 16, foundRun + diff);
    string strRunNum = file.substr (foundRun + 16, diff);
    //int runNum = atoi(strRunNum.c_str());

    //return runNum;
    return strRunNum;
    
}



string getInfo_fromFilename_oldfiles (const char *runfile) {
    
    string file = string( runfile );
    string chRun = ".setup_";
    string chRun2 = ".run";
    size_t foundRun=file.find(chRun);
    size_t foundRun2=file.find(chRun2);
    if ( foundRun < 10 ) {
    
        chRun = "Tdomain_TestBed_";
        foundRun=file.find(chRun);
    }

    cout<<"file name first : "<<foundRun<<", second : "<<foundRun2<<endl;
    int diff = foundRun2 - foundRun;
    //string strRunNum = file.substr (foundRun + 16, foundRun + diff);
    string strRunNum = file.substr (foundRun + 7, diff);
    //int runNum = atoi(strRunNum.c_str());

    //return runNum;
    return strRunNum;
    
}



string getInfo_fromFilename_FdomainWoBug (const char *runfile) {
    
    string file = string( runfile );
    string chRun = "_ARA37_";
    string chRun2 = ".run";
    size_t foundRun=file.find(chRun);
    size_t foundRun2=file.find(chRun2);
    if ( foundRun < 10 ) {
    
        chRun = "Tdomain_TestBed_";
        foundRun=file.find(chRun);
    }

    cout<<"file name first : "<<foundRun<<", second : "<<foundRun2<<endl;
    int diff = foundRun2 - foundRun;
    //string strRunNum = file.substr (foundRun + 16, foundRun + diff);
    string strRunNum = file.substr (foundRun + 7, diff);
    //int runNum = atoi(strRunNum.c_str());

    //return runNum;
    return strRunNum;
    
}




string getInfo_fromFilename_TdomainWBug (const char *runfile) {
    
    string file = string( runfile );
    string chRun = ".setup_";
    string chRun2 = ".run";
    size_t foundRun=file.find(chRun);
    size_t foundRun2=file.find(chRun2);
    if ( foundRun < 10 ) {
    
        chRun = "Tdomain_TestBed_";
        foundRun=file.find(chRun);
    }

    cout<<"file name first : "<<foundRun<<", second : "<<foundRun2<<endl;
    int diff = foundRun2 - foundRun;
    //string strRunNum = file.substr (foundRun + 16, foundRun + diff);
    string strRunNum = file.substr (foundRun + 7, diff);
    //string strRunNum = file.substr (foundRun + 7, diff-4);
    //int runNum = atoi(strRunNum.c_str());

    //return runNum;
    return strRunNum;
    
}


