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

#include "Ray.h"

//#include "FFTtools.h"

class EarthModel; //class


string outputdir="outputs";


int getPeakBin(TGraph *gr);

double getPeak(TGraph *gr);


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
      cout<<"too many info! just use default AraOut.root file!"<<endl;
      readfile = "outputs/AraOut.root";
  }


//  Settings *settings = new Settings();

  //  Detector *detector=new Detector(settings->DETECTOR); // builds antenna array, 0 for testbed
//  Detector *detector=0; // builds antenna array, 0 for testbed
  Detector *detector = 0; 
  Settings *settings = 0;
  Spectra *spectra = 0;
  IceModel *icemodel = 0;
  Event *event = 0;
  Report *report = 0;
  Trigger *trigger = 0;
  cout<<"construct detector"<<endl;

  
  TFile *AraFile=new TFile(( readfile ).c_str());
  //TFile *AraFile=new TFile((outputdir+"/AraOut.root").c_str());
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
  cout<<"getevent"<<endl;
  cout << "I'm here.\n";
  cout << "GetGain(700,5,0,0) is " << detector->GetGain(700., 5., 0., 0) << "\n";
  cout << "GetGain(10,5,0,0) is " << detector->GetGain(10., 5., 0., 0) << "\n";

  cout<<"station x is "<<detector->stations[0].GetX()<<endl;
  cout<<"string x is "<<detector->stations[0].strings[0].GetX()<<endl;
  cout<<"antenna x is "<<detector->stations[0].strings[0].antennas[0].GetX()<<endl;
  cout<<"antenna Gain(700,5,0) is "<<detector->stations[0].strings[0].antennas[0].GetG(detector,700.,5.,0.)<<endl;

  cout<<"params.number_of_stations : "<<detector->params.number_of_stations<<endl;
  cout<<"params.station_spacing : "<<detector->params.station_spacing<<endl;

  cout<<"\n"<<endl;
  cout<<"Settings->NNU : "<<settings->NNU<<endl;
  cout<<"Settings->DETECTOR : "<<settings->DETECTOR<<endl;


cout<<"random energy from Spectra : "<<spectra->GetNuEnergy()<<endl;

cout<<"Detector static const double freq_init : "<<detector->Getfreq_init()<<endl;


cout<<"Detector -> freq_forfft[0] : "<<detector->freq_forfft[0]<<endl;
cout<<"Detector -> freq_forfft[9] : "<<detector->freq_forfft[9]<<endl;
cout<<"Detector -> freq_forfft[100] : "<<detector->freq_forfft[100]<<endl;

cout<<"icemodel surface : "<<icemodel->Surface(0.,0.)<<endl;

AraTree2->GetEvent(0);
cout<<"nnu x : "<<event->Nu_Interaction[0].nnu.GetX()<<endl;
AraTree2->GetEvent(1);
cout<<"nnu x : "<<event->Nu_Interaction[0].nnu.GetX()<<endl;

  int nnu_pass = 0; // number of nu events which passed PickUnbiased.
  double posnuX[settings->NNU];
  double posnuY[settings->NNU];
  double posnuR[settings->NNU];


  // global pass evt count
  int pass_evts = 0;

  int string, ant;
  int ch_from_det;



  // test if InstalledStations are working properly
//  cout<<"InstalledStations size : "<<detector->InstalledStations.size()<<endl;
//  cout<<"InstalledStations nStrings size : "<<detector->InstalledStations[0].nStrings<<endl;
//
//  cout<<"InstalledStations VHChannle size : "<<detector->InstalledStations[0].VHChannel.size()<<endl;
//  for (int i = 0; i < int(detector->InstalledStations[0].VHChannel.size()); i++){
//      cout<<"InstalledStations VHChannel["<<i<<"] size : "<<detector->InstalledStations[0].VHChannel[i].size()<<endl;
//  }



  double peak_plot_max = 1.e3;
  TH1D *hpeak[16];

  char title_test[100];

  for (int i=0; i<16; i++) {
      sprintf( title_test, "hpeakV_ch%d",i);
      hpeak[i] = new TH1D(title_test, "", 100, 0., peak_plot_max);
  }




    TFile *PeakVFile = new TFile("./outputs/plots/PeakVFile.root","RECREATE");
    TTree *Tree=new TTree("Tree","TestTree");

    Tree->Branch ("hpeakVch0", &hpeak[0]);
    Tree->Branch ("hpeakVch1", &hpeak[1]);
    Tree->Branch ("hpeakVch2", &hpeak[2]);
    Tree->Branch ("hpeakVch3", &hpeak[3]);
    Tree->Branch ("hpeakVch4", &hpeak[4]);
    Tree->Branch ("hpeakVch5", &hpeak[5]);
    Tree->Branch ("hpeakVch6", &hpeak[6]);
    Tree->Branch ("hpeakVch7", &hpeak[7]);






  int bin = settings->NFOUR/2;
  double wf_time_offset = -450.;// in ns


  stringstream sstr;

  //char title_test[100];

  TH1D *hist_peak[16];
  for (int i=0; i<16; i++) {
      sprintf( title_test, "hpeakV_ch%d",i);
      hist_peak[i] = new TH1D(title_test, "", 50, 1, bin);
  }

  int total_evt = settings->NNU;
  int num_plots = 10;


  for (int inu=0;inu<total_evt;inu++) { // loop over neutrinos


      AraTree2->GetEvent(inu);

      /*
      // save X, Y of posnus which passed PickUnbiased
      if ( event->Nu_Interaction[0].pickposnu ) {
          posnuX[nnu_pass] = event->Nu_Interaction[0].posnu.GetX();
          posnuY[nnu_pass] = event->Nu_Interaction[0].posnu.GetY();
          posnuR[nnu_pass] = event->Nu_Interaction[0].posnu.R();
          nnu_pass++;

          cout<<"evt no "<<inu<<"stations[0].strings[1].antennas[2].ray_sol_cnt : "<<report->stations[0].strings[1].antennas[2].ray_sol_cnt<<endl;
          */

/*
          if ( interaction->ray_solver_toggle ) {   // if ray_solver succeeded to get soutions
              cout<<"pass evt : "<<nnu_pass<<"\t";
              for (int i=0; i<interaction->ray_output[0].size(); i++) {
                  for (int j=0; j<interaction->ray_output.size(); j++) {
                      cout<<j<<"th ray_output : "<<interaction->ray_output[j][i]<<"\t";
                  }
                  cout<<"\n";
              }
          }// end if ray_solver_toggle
*/
      //}

      //if ( report->stations[0].Global_Pass ) {
      if ( report->stations[0].Global_Pass > 0) {

          pass_evts++;
          cout<<"passed evt : "<<pass_evts<<endl;

          //if (pass_evts == 1) {
          //if (pass_evts%50 == 0) {

              TCanvas *cTest = new TCanvas ("cTest","", 3200,3200);
              cTest->Divide(4,4);

              TGraph *gTest[16];

              //for ( int chID=0; chID<16; chID++) {
              for ( int chID=0; chID<detector->stations[0].number_of_antennas; chID++) {

                  detector->GetSSAfromChannel(0, chID+1, &ant, &string, settings);

                  //cout<<"save wf for ch "<<chID<<" ";
                  //cout<<"string:"<<string<<" ant:"<<ant<<" type:"<<detector->stations[0].strings[string].antennas[ant].type<<" ";
                  //cout<<"DAQchan:"<<detector->stations[0].strings[string].antennas[ant].DAQchan<<endl;

                  // plot waveform
                  double getx[bin];
                  double gety[bin];
                  for (int l=0; l<bin; l++) {

                      getx[l] = report->stations[0].strings[string].antennas[ant].time_mimic[l];// in ns
                      gety[l] = report->stations[0].strings[string].antennas[ant].V_mimic[l];
                      //cout<<"getx["<<l<<"] : "<<getx[l]<<"\tgety["<<l<<"] : "<<gety[l]<<endl;
                  }


                  gTest[chID] = new TGraph (bin, getx, gety);

                  cTest->cd(chID+1);
                  gTest[chID]->Draw("AL");


                  // fill peak hist
                  if (chID<16) { // only for BH chs
                      hist_peak[chID] -> Fill(getPeakBin(gTest[chID]) );
                      hpeak[chID] -> Fill( getPeak(gTest[chID]) );
                  }

              }// for chID


          //if (pass_evts%50 == 0) {
          //if (pass_evts% (total_evt/num_plots) == 0) {
          if (pass_evts < num_plots ) {

                  sstr.str("");
                  sstr<<"./outputs/plots/WF_AraOut.event"<<inu<<".pdf";
                  cTest->SaveAs(sstr.str().c_str());
          }



              delete cTest;
              //for ( int chID=0; chID<16; chID++) {
              for ( int chID=0; chID<detector->stations[0].number_of_antennas; chID++) {
                  delete gTest[chID];
              }



          //}// if pass_evt


      }// if global passed



  } // end loop over neutrinos


  // make canvas for hist peak and plot it
    TCanvas *cPeakbin = new TCanvas ("cPeakbin","", 3200,3200);
    cPeakbin->Divide(4,4);
    for (int i=0; i<16; i++) {
        cPeakbin->cd(i+1);
        hist_peak[i]->Draw();
    }
    sstr.str("");
    sstr<<"./outputs/plots/HistPeak_AraOut.pdf";
    cPeakbin->SaveAs(sstr.str().c_str());


    delete cPeakbin;
    for (int i=0; i<16; i++) {
        delete hist_peak[i];
    }

    int bin_filters = 60;

    double freq[bin_filters];
    double filterdb[bin_filters];
    double preampdb[bin_filters];
    double FOAMdb[bin_filters];

    for (int i=0; i<bin_filters; i++){
        //freq[i] = detector->GetFreq(i)*1.e6;// in MHz
        freq[i] = detector->GetFreq(i);
        filterdb[i] = detector->GetFilterGain(i);
        preampdb[i] = detector->GetPreampGain(i);
        FOAMdb[i] = detector->GetFOAMGain(i);
    }

    TGraph *gfilter = new TGraph (bin_filters, freq, filterdb);
    TGraph *gpreamp = new TGraph (bin_filters, freq, preampdb);
    TGraph *gFOAM = new TGraph (bin_filters, freq, FOAMdb);

    TCanvas *cFilters = new TCanvas ("cFilters","", 1600,1600);
    cFilters->Divide(2,2);

    cFilters->cd(1);
    gfilter->SetTitle("Filter gain");
    gfilter->GetHistogram()->SetXTitle("Freq (MHz)");
    gfilter->GetHistogram()->SetYTitle("Gain (dB)");
    gfilter->Draw("AL");
    cFilters->cd(2);
    gpreamp->SetTitle("Preamp gain");
    gpreamp->GetHistogram()->SetXTitle("Freq (MHz)");
    gpreamp->GetHistogram()->SetYTitle("Gain (dB)");
    gpreamp->Draw("AL");
    cFilters->cd(3);
    gFOAM->SetTitle("FOAM gain");
    gFOAM->GetHistogram()->SetXTitle("Freq (MHz)");
    gFOAM->GetHistogram()->SetYTitle("Gain (dB)");
    gFOAM->Draw("AL");

    sstr.str("");
    sstr<<"./outputs/plots/Filters_AraOut.pdf";
    cFilters->SaveAs(sstr.str().c_str());
    
    delete gfilter;
    delete gpreamp;
    delete gFOAM;
    delete cFilters;



    int bin_DB = settings->DATA_BIN_SIZE/2;

    double df_fft = 1./ ( (double)(settings->DATA_BIN_SIZE) * settings->TIMESTEP);


    double freq_DB[bin_DB];
    double filterdb_DB[bin_DB];
    double preampdb_DB[bin_DB];
    double FOAMdb_DB[bin_DB];

    for (int i=0; i<bin_DB; i++){
        //freq[i] = detector->GetFreq(i)*1.e6;// in MHz
        freq_DB[i] = (double)i * df_fft / (1.e6); // in MHz
        filterdb_DB[i] = detector->GetFilterGain_databin(i);
        preampdb_DB[i] = detector->GetPreampGain_databin(i);
        FOAMdb_DB[i] = detector->GetFOAMGain_databin(i);
    }

    TGraph *gfilter_DB = new TGraph (bin_DB, freq_DB, filterdb_DB);
    TGraph *gpreamp_DB = new TGraph (bin_DB, freq_DB, preampdb_DB);
    TGraph *gFOAM_DB = new TGraph (bin_DB, freq_DB, FOAMdb_DB);

    TCanvas *cFilters_DB = new TCanvas ("cFilters_DB","", 1600,1600);
    cFilters_DB->Divide(2,2);

    cFilters_DB->cd(1);
    gfilter_DB->SetTitle("Filter gain");
    gfilter_DB->GetHistogram()->SetXTitle("Freq (MHz)");
    gfilter_DB->GetHistogram()->SetYTitle("Gain (dB)");
    gfilter_DB->Draw("AL");
    cFilters_DB->cd(2);
    gpreamp_DB->SetTitle("Preamp gain");
    gpreamp_DB->GetHistogram()->SetXTitle("Freq (MHz)");
    gpreamp_DB->GetHistogram()->SetYTitle("Gain (dB)");
    gpreamp_DB->Draw("AL");
    cFilters_DB->cd(3);
    gFOAM_DB->SetTitle("FOAM gain");
    gFOAM_DB->GetHistogram()->SetXTitle("Freq (MHz)");
    gFOAM_DB->GetHistogram()->SetYTitle("Gain (dB)");
    gFOAM_DB->Draw("AL");

    sstr.str("");
    sstr<<"./outputs/plots/Filters_databinsize_AraOut.pdf";
    cFilters_DB->SaveAs(sstr.str().c_str());
    
    delete gfilter_DB;
    delete gpreamp_DB;
    delete gFOAM_DB;
    delete cFilters_DB;



  


    TCanvas *cPeakV = new TCanvas ("cPeakV","", 3200,3200);
    cPeakV->Divide(4,4);

    for (int i=0; i<16; i++) {

        cPeakV->cd(i+1);
        cPeakV->cd(i+1)->SetLogy();
        hpeak[i]->SetLineColor(2);
        hpeak[i]->Draw();

    }
    sstr.str("");
    sstr<<"./outputs/plots/histpeakV_AraOut.pdf";
    cPeakV->SaveAs(sstr.str().c_str());


    Tree->Fill();
    //hList->Write();
    PeakVFile->Write();
    PeakVFile->Close();


    delete cPeakV;
    for (int i=0; i<16; i++) {
        delete hpeak[i];
    }






    cout<<"passed evt : "<<pass_evts<<endl;


}




int getPeakBin(TGraph *gr) 
{
  double x,y;
  gr->GetPoint(0,x,y);
  double peakVal=y;
  int peakBin=0;
  for(int i=1;i<gr->GetN();i++) {
    gr->GetPoint(i,x,y);
    if(peakVal<y) {
      peakVal=y;
      peakBin=i;
    }      
  }
  return peakBin;
}

double getPeak(TGraph *gr)
{
  double x,y;
  gr->GetPoint(0,x,y);
  double peakVal=y*y;
  int peakBin=0;
  for(int i=1;i<gr->GetN();i++) {
    gr->GetPoint(i,x,y);
    if( peakVal<(y*y) ) {
      peakVal=(y*y);
      peakBin=i;
    }      
  }
  //return peakBin;
  return sqrt(peakVal);
}




