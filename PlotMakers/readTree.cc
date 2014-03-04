#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
//#include <stdlib.formh>
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

class EarthModel; //class
TStyle *RootStyle();
TStyle *color=RootStyle();

string outputdir="outputs";

int main() {
 gStyle=color;

//  Settings *settings = new Settings();

  //  Detector *detector=new Detector(settings->DETECTOR); // builds antenna array, 0 for testbed
//  Detector *detector=0; // builds antenna array, 0 for testbed
  Detector *detector = 0; 
  Settings *settings = 0;
  Spectra *spectra = 0;
  IceModel *icemodel = 0;
  Event *event = 0;
  Report *report = 0;
  //Declare number of trees to be read
  const int NUM_TREE=2;
  cout<<"construct detector"<<endl;

  TFile *AraFile[NUM_TREE];
  AraFile[0]=new TFile((outputdir+"/AraOut.root").c_str());
  AraFile[1]=new TFile((outputdir+"/AraOut2.root").c_str());
  cout<<"AraFile1"<<endl;
  cout<<"AraFile2"<<endl;
  TTree *AraTree[NUM_TREE];
TTree *AraTree2[NUM_TREE];
AraTree[0]=(TTree*)AraFile[0]->Get("AraTree");
AraTree[1]=(TTree*)AraFile[1]->Get("AraTree");
  AraTree2[0]=(TTree*)AraFile[0]->Get("AraTree2");
AraTree2[1]=(TTree*)AraFile[1]->Get("AraTree2");
  cout<<"AraTree"<<endl;
  AraTree[0]->SetBranchAddress("detector",&detector);
  AraTree[0]->SetBranchAddress("settings",&settings);
  AraTree[0]->SetBranchAddress("spectra",&spectra);
  AraTree[0]->SetBranchAddress("icemodel",&icemodel);
  AraTree2[0]->SetBranchAddress("event",&event);
  AraTree2[0]->SetBranchAddress("report",&report);
  
  //////////////////////////////////////
  //Redeclare everything for second AraTree
  AraTree[1]->SetBranchAddress("detector",&detector);
  AraTree[1]->SetBranchAddress("settings",&settings);
  AraTree[1]->SetBranchAddress("spectra",&spectra);
  AraTree[1]->SetBranchAddress("icemodel",&icemodel);
  AraTree2[1]->SetBranchAddress("event",&event);
  AraTree2[1]->SetBranchAddress("report",&report);
  ///////////////////////////////////////
  cout<<"branch detector"<<endl;

  AraTree[0]->GetEvent(0);
AraTree[1]->GetEvent(1);  
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

AraTree2[0]->GetEvent(0);
AraTree2[1]->GetEvent(0);

cout<<"nnu x : "<<event->nnu.GetX()<<endl;
AraTree2[0]->GetEvent(1);
AraTree2[1]->GetEvent(1);

cout<<"nnu x : "<<event->nnu.GetX()<<endl;

  int nnu_pass = 0; // number of nu events which passed PickUnbiased.
  double posnuX[settings->NNU];
  double posnuY[settings->NNU];
  double posnuR[settings->NNU];

  for (int inu=0;inu<settings->NNU;inu++) { // loop over neutrinos


      AraTree2[0]->GetEvent(inu);

      // save X, Y of posnus which passed PickUnbiased
      if ( event->Nu_Interaction[0].pickposnu ) {
          posnuX[nnu_pass] = event->Nu_Interaction[0].posnu.GetX();
          posnuY[nnu_pass] = event->Nu_Interaction[0].posnu.GetY();
          posnuR[nnu_pass] = event->Nu_Interaction[0].posnu.R();
          nnu_pass++;

          cout<<"evt no "<<inu<<"stations[0].strings[1].antennas[2].ray_sol_cnt : "<<report->stations[0].strings[1].antennas[2].ray_sol_cnt<<endl;

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
      }


  } // end loop over neutrinos



///////////////////////////////////////////
//  test Detector class
///////////////////////////////////////////




if ( settings->DETECTOR == 0 ) {

cout<<"\n\t Test reading antenna array infomation !"<<endl;
//cout<<"total number of strings : "<<(int)detector->Detector.params.number_of_strings<<endl;
cout<<"total number of strings : "<<(int)detector->params.number_of_strings<<endl;
cout<<"total number of antennas : "<<(int)detector->params.number_of_antennas<<endl;
 cout << "\n\n\n Testing \n\n\n";
cout<<"\nantenna0 position is"<<endl;
cout<<"x : "<<(double)detector->strings[0].GetX()<<endl;
cout<<"y : "<<(double)detector->strings[0].GetY()<<endl;
cout<<"z1 : "<<(double)detector->strings[0].antennas[0].GetZ()<<" type : "<<(int)detector->strings[0].antennas[0].type<<endl;
cout<<"z2 : "<<(double)detector->strings[0].antennas[1].GetZ()<<" type : "<<(int)detector->strings[0].antennas[1].type<<endl;

cout<<"\nantenna1 position is"<<endl;
cout<<"x : "<<(double)detector->strings[1].GetX()<<endl;
cout<<"y : "<<(double)detector->strings[1].GetY()<<endl;
cout<<"z1 : "<<(double)detector->strings[1].antennas[0].GetZ()<<" type : "<<(int)detector->strings[1].antennas[0].type<<endl;
cout<<"z2 : "<<(double)detector->strings[1].antennas[1].GetZ()<<" type : "<<(int)detector->strings[1].antennas[1].type<<endl;

cout<<"\nantenna2 position is"<<endl;
cout<<"x : "<<(double)detector->strings[2].GetX()<<endl;
cout<<"y : "<<(double)detector->strings[2].GetY()<<endl;
cout<<"z1 : "<<(double)detector->strings[2].antennas[0].GetZ()<<" type : "<<(int)detector->strings[2].antennas[0].type<<endl;
cout<<"z2 : "<<(double)detector->strings[2].antennas[1].GetZ()<<" type : "<<(int)detector->strings[2].antennas[1].type<<endl;

cout<<"\nantenna3 position is"<<endl;
cout<<"x : "<<(double)detector->strings[3].GetX()<<endl;
cout<<"y : "<<(double)detector->strings[3].GetY()<<endl;
cout<<"z1 : "<<(double)detector->strings[3].antennas[0].GetZ()<<" type : "<<(int)detector->strings[3].antennas[0].type<<endl;
cout<<"z2 : "<<(double)detector->strings[3].antennas[1].GetZ()<<" type : "<<(int)detector->strings[3].antennas[1].type<<endl;

}




/*

else if ( settings->DETECTOR == 1 ) {

cout<<"\n\t Test ARA-N array setting !"<<endl;
//cout<<"total number of strings : "<<(int)detector->Detector.params.number_of_strings<<endl;
cout<<"total number of strings : "<<(int)detector->params.number_of_strings<<endl;
cout<<"total number of antennas : "<<(int)detector->params.number_of_antennas<<endl;


TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,2400,700);
c1->Divide(3,1);

double x[(int)detector->params.number_of_stations], y[(int)detector->params.number_of_stations];

for (int i=0;i<(int)detector->params.number_of_stations;i++) {
    x[i] = (double)detector->stations[i].GetX();
    y[i] = (double)detector->stations[i].GetY();
}

TGraph *gr;
gr = new TGraph((int)detector->params.number_of_stations,x,y);

c1->cd(1);
gr->SetTitle("ARA station layout");
gr->GetHistogram()->SetMaximum(10000);
gr->GetHistogram()->SetMinimum(-10000);
gr->GetXaxis()->SetLimits(-10000,10000);
gr->GetHistogram()->SetXTitle("X (m)");
gr->GetHistogram()->SetYTitle("Y (m)");
gr->Draw("a*");


c1->cd(2);

int station_choice = 0;
//int station_choice = (int)detector->params.number_of_stations - 1;

double string_x[4], string_y[4];
double surface_x[4], surface_y[4];

for (int i=0;i<4;i++) {
    string_x[i] = (double)detector->stations[station_choice].strings[i].GetX();
    string_y[i] = (double)detector->stations[station_choice].strings[i].GetY();

    surface_x[i] = (double)detector->stations[station_choice].surfaces[i].GetX();
    surface_y[i] = (double)detector->stations[station_choice].surfaces[i].GetY();
}

TGraph *gr_string;
gr_string = new TGraph(4,string_x,string_y);

gr_string->SetTitle("Strings and surface antennas layout for each station");
//gr_string->GetHistogram()->SetMaximum(5300);
//gr_string->GetHistogram()->SetMinimum(5100);
//gr_string->GetXaxis()->SetLimits(-3100,-2900);
gr_string->GetHistogram()->SetMaximum( (int)detector->stations[station_choice].GetY() + 100);
gr_string->GetHistogram()->SetMinimum( (int)detector->stations[station_choice].GetY() - 100);
gr_string->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
gr_string->SetMarkerColor(4);
gr_string->SetMarkerSize(2);
gr_string->SetMarkerStyle(20);
gr_string->GetHistogram()->SetXTitle("X (m)");
gr_string->GetHistogram()->SetYTitle("Y (m)");
gr_string->Draw("ap");

TGraph *gr_surface;
gr_surface = new TGraph(4,surface_x,surface_y);
gr_surface->SetMarkerColor(2);
gr_surface->SetMarkerSize(2);
gr_surface->SetMarkerStyle(21);
gr_surface->Draw("p");


TLegend *Leg_string_surface = new TLegend(1., 0.95, 0.5,0.8);
Leg_string_surface -> AddEntry(gr_string, "Strings");
Leg_string_surface -> AddEntry(gr_surface, "Surface antennas");
Leg_string_surface -> Draw();


c1->cd(3);

int string_choice = 0;
//int station_choice = (int)detector->params.number_of_stations - 1;

double antenna_x[4], antenna_y[4];   // use x as x, y as z to see the depth layout

for (int i=0;i<4;i++) {
    antenna_x[i] = (double)detector->stations[station_choice].strings[string_choice].GetX();
    antenna_y[i] = (double)detector->stations[station_choice].strings[string_choice].antennas[i].GetZ();
}

TGraph *gr_antenna;
gr_antenna = new TGraph(4,antenna_x,antenna_y);

gr_antenna->SetTitle("Borehole antenna layout for each string");
//gr_string->GetHistogram()->SetMaximum(5300);
//gr_string->GetHistogram()->SetMinimum(5100);
//gr_string->GetXaxis()->SetLimits(-3100,-2900);
gr_antenna->GetHistogram()->SetMaximum( (int)detector->stations[station_choice].strings[string_choice].GetZ() );
gr_antenna->GetHistogram()->SetMinimum( (int)detector->stations[station_choice].strings[string_choice].antennas[0].GetZ() - 20);
gr_antenna->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
gr_antenna->GetYaxis()->SetTitle("Z (depth, m)");
gr_antenna->SetMarkerColor(4);
gr_antenna->SetMarkerSize(2);
gr_antenna->SetMarkerStyle(20);
gr_antenna->GetHistogram()->SetXTitle("X (m)");
gr_antenna->GetHistogram()->SetYTitle("Z (m)");
gr_antenna->Draw("ap");


c1->Print("ARA-37_station_layout.pdf");






}





*/






//else if ( settings->DETECTOR == 2 ) {
else {

cout<<"\n\t Test ARA-37 array setting !"<<endl;
//cout<<"total number of strings : "<<(int)detector->Detector.params.number_of_strings<<endl;
cout<<"total number of strings : "<<(int)detector->params.number_of_strings<<endl;
cout<<"total number of antennas : "<<(int)detector->params.number_of_antennas<<endl;


TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,4800,700);
c1->Divide(6,1);

double x[(int)detector->params.number_of_stations], y[(int)detector->params.number_of_stations];

for (int i=0;i<(int)detector->params.number_of_stations;i++) {
    x[i] = (double)detector->stations[i].GetX();
    y[i] = (double)detector->stations[i].GetY();
}

TGraph *gr;
gr = new TGraph((int)detector->params.number_of_stations,x,y);

c1->cd(1);
gr->SetTitle("ARA station layout");
gr->GetHistogram()->SetMaximum(detector->params.core_y + 10000);
gr->GetHistogram()->SetMinimum(detector->params.core_y - 10000);
gr->GetHistogram()->SetXTitle("X (m)");
gr->GetHistogram()->SetYTitle("Y (m)");
gr->GetYaxis()->SetTitleOffset(1.2);
gr->GetHistogram()->SetYTitle("Y (m)");
gr->GetXaxis()->SetLimits(detector->params.core_x-10000, detector->params.core_x+10000);
gr->Draw("a*");


c1->cd(2);

int station_choice = 0;

double string_x[4], string_y[4];
double surface_x[4], surface_y[4];

for (int i=0;i<4;i++) {
    string_x[i] = (double)detector->stations[station_choice].strings[i].GetX();
    string_y[i] = (double)detector->stations[station_choice].strings[i].GetY();

    surface_x[i] = (double)detector->stations[station_choice].surfaces[i].GetX();
    surface_y[i] = (double)detector->stations[station_choice].surfaces[i].GetY();
}

TGraph *gr_string;
gr_string = new TGraph(4,string_x,string_y);

gr_string->SetTitle("Strings and surface antennas layout for each station");
//gr_string->GetHistogram()->SetMaximum(5300);
//gr_string->GetHistogram()->SetMinimum(5100);
//gr_string->GetXaxis()->SetLimits(-3100,-2900);
gr_string->GetHistogram()->SetMaximum( (int)detector->stations[station_choice].GetY() + 100);
gr_string->GetHistogram()->SetMinimum( (int)detector->stations[station_choice].GetY() - 100);
gr_string->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
gr_string->SetMarkerColor(4);
gr_string->SetMarkerSize(2);
gr_string->SetMarkerStyle(20);
gr_string->GetHistogram()->SetXTitle("X (m)");
gr_string->GetHistogram()->SetYTitle("Y (m)");
gr_string->GetYaxis()->SetTitleOffset(1.2);
gr_string->Draw("ap");

TGraph *gr_surface;
gr_surface = new TGraph(4,surface_x,surface_y);
gr_surface->SetMarkerColor(2);
gr_surface->SetMarkerSize(2);
gr_surface->SetMarkerStyle(21);
gr_surface->Draw("p");


TLegend *Leg_string_surface = new TLegend(1., 0.95, 0.5,0.8);
Leg_string_surface -> AddEntry(gr_string, "Strings");
Leg_string_surface -> AddEntry(gr_surface, "Surface antennas");
Leg_string_surface -> Draw();



c1->cd(3);

int string_choice = 0;
//int station_choice = (int)detector->params.number_of_stations - 1;

double antenna_x[4], antenna_y[4];   // use x as x, y as z to see the depth layout

for (int i=0;i<4;i++) {
    antenna_x[i] = (double)detector->stations[station_choice].strings[string_choice].GetX();
    antenna_y[i] = (double)detector->stations[station_choice].strings[string_choice].antennas[i].GetZ();
}

TGraph *gr_antenna;
gr_antenna = new TGraph(4,antenna_x,antenna_y);

gr_antenna->SetTitle("Borehole antenna layout for each string");
//gr_string->GetHistogram()->SetMaximum(5300);
//gr_string->GetHistogram()->SetMinimum(5100);
//gr_string->GetXaxis()->SetLimits(-3100,-2900);
gr_antenna->GetHistogram()->SetMaximum( (int)detector->stations[station_choice].strings[string_choice].GetZ() );
gr_antenna->GetHistogram()->SetMinimum( (int)detector->stations[station_choice].strings[string_choice].antennas[0].GetZ() - 20);
gr_antenna->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
gr_antenna->GetYaxis()->SetTitle("Z (depth, m)");
gr_antenna->SetMarkerColor(4);
gr_antenna->SetMarkerSize(2);
gr_antenna->SetMarkerStyle(20);
gr_antenna->GetHistogram()->SetXTitle("X (m)");
gr_antenna->GetYaxis()->SetTitleOffset(1.2);
gr_antenna->Draw("ap");



c1->cd(4);


double station_x[(int)detector->params.number_of_stations], station_z[(int)detector->params.number_of_stations];

for (int i=0;i<(int)detector->params.number_of_stations;i++) {
    station_x[i] = (double)detector->stations[i].GetX();
    station_z[i] = (double)detector->stations[i].GetZ();
//    station_z[i] = (double)detector->stations[0].GetZ();
}

TGraph *gr_crosssection;
gr_crosssection = new TGraph((int)detector->params.number_of_stations,station_x,station_z);

gr_crosssection->SetTitle("Station layout CrossSection");
//--------------------------------------------------
// gr_crosssection->GetHistogram()->SetMaximum( (int)detector->stations[0].GetZ() + 10 );
// gr_crosssection->GetHistogram()->SetMinimum( (int)detector->stations[0].GetZ() - 10 );
//-------------------------------------------------- 
gr_crosssection->GetHistogram()->SetMaximum( detector->stations[detector->params.number_of_stations-1].GetZ() + 1000 );
gr_crosssection->GetHistogram()->SetMinimum( detector->stations[detector->params.number_of_stations-1].GetZ() - 1000 );
gr_crosssection->GetXaxis()->SetLimits(detector->params.core_x-10000, detector->params.core_x+10000);
gr_crosssection->GetHistogram()->SetXTitle("X (m)");
gr_crosssection->GetHistogram()->SetYTitle("Z (m)");
gr_crosssection->GetYaxis()->SetTitleOffset(1.2);
gr_crosssection->Draw("a*");


c1->cd(5);

TGraph *gr_posnu;
gr_posnu = new TGraph(nnu_pass,posnuX,posnuY);

gr_posnu->SetTitle("posnu");
//--------------------------------------------------
// gr_posnu->GetHistogram()->SetMaximum( icemodel->Surface(0.,0.)*sin(30.*RADDEG) );
// gr_posnu->GetHistogram()->SetMinimum( -icemodel->Surface(0.,0.)*sin(30.*RADDEG)  );
//-------------------------------------------------- 
gr_posnu->GetHistogram()->SetMaximum( detector->params.core_y+3000. );
gr_posnu->GetHistogram()->SetMinimum( detector->params.core_y-3000.  );
gr_posnu->GetHistogram()->SetXTitle("X (m)");
gr_posnu->GetHistogram()->SetYTitle("Y (m)");
gr_posnu->GetYaxis()->SetTitleOffset(1.2);
gr_posnu->GetHistogram()->SetYTitle("Y (m)");
gr_posnu->GetXaxis()->SetLimits( detector->params.core_x-3000., detector->params.core_x+3000.);
//--------------------------------------------------
// gr_posnu->GetXaxis()->SetLimits(-icemodel->Surface(0.,0.)*sin(30.*RADDEG),icemodel->Surface(0.,0.)*sin(30.*RADDEG));
//-------------------------------------------------- 
gr_posnu->Draw("a*");


c1->cd(6);

double cd6_x[settings->NNU];
double cd6_y_top[settings->NNU];
double cd6_y_bot[settings->NNU];

for(int i =0; i<settings->NNU; i++) {
      
    AraTree2[0]->GetEvent(i);

    cd6_x[i] = i;
    cd6_y_top[i] = icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() );
    cd6_y_bot[i] = icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() ) - icemodel->IceThickness( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() );
    
    if (icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() ) - posnuR[i] < 0) { 
        cout<<"Surface - posnuR : "<<icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() ) - posnuR[i]<<endl;
        cout<<"!offsurface"<<endl;
    }

    if (posnuR[i] - ( icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat()) - icemodel->IceThickness(event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat()) ) < 0)  {
        cout<<"posnuR - Icebottom : "<<posnuR[i] - ( icemodel->Surface( event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat() ) - icemodel->IceThickness(event->Nu_Interaction[0].posnu.Lon(), event->Nu_Interaction[0].posnu.Lat()) ) <<endl;
        cout<<"!offbottomice"<<endl;
    }

}


TGraph *gr_depth;
gr_depth = new TGraph(settings->NNU, cd6_x, posnuR);

TGraph *gr_top;
gr_top = new TGraph(settings->NNU, cd6_x, cd6_y_top);

TGraph *gr_bot;
gr_bot = new TGraph(settings->NNU, cd6_x, cd6_y_bot);


gr_depth->SetTitle("posnu_position");
gr_depth->GetHistogram()->SetMaximum( cd6_y_top[0] + 1000. );
gr_depth->GetHistogram()->SetMinimum( cd6_y_bot[0] - 1000.  );
gr_depth->GetHistogram()->SetXTitle("posnu evt number");
gr_depth->GetHistogram()->SetYTitle("depth (m)");
gr_depth->GetYaxis()->SetTitleOffset(1.2);
gr_depth->SetMarkerColor(2);
gr_depth->SetMarkerSize(1);
gr_depth->SetMarkerStyle(21);
gr_depth->Draw("ap");

gr_top->SetMarkerColor(3);
gr_top->SetMarkerSize(1);
gr_top->SetMarkerStyle(21);
gr_top->Draw("p");

gr_bot->SetMarkerColor(4);
gr_bot->SetMarkerSize(1);
gr_bot->SetMarkerStyle(21);
gr_bot->Draw("p");


TLegend *Leg_top_bot = new TLegend(1., 0.95, 0.5,0.8);
Leg_top_bot -> AddEntry(gr_top, "Ice Surface");
Leg_top_bot -> AddEntry(gr_depth, "posnu depth");
Leg_top_bot -> AddEntry(gr_bot, "Bedrock");
Leg_top_bot -> Draw();

c1->Print("ARA-37_station_layout.pdf");




}



// test
cout<<"station[0] x : "<<detector->stations[0].GetX()<<endl;
cout<<"string[0] x : "<<detector->stations[0].strings[0].GetX()<<endl;
cout<<"antenna[0] x : "<<detector->stations[0].strings[0].antennas[0].GetX()<<endl;
cout<<"antenna[0] Gain(700,5,0) : "<<detector->stations[0].strings[0].antennas[0].GetG(detector,700.,5.,0.)<<endl;
cout<<"GetGain(700,5,0,0) : "<<detector->GetGain(700.,5.,0.,0)<<endl;
cout<<"GetGain(10,5,0,0) : "<<detector->GetGain(10.,5.,0.,0)<<endl;

cout<<"params.number_of_stations : "<<detector->params.number_of_stations<<endl;
cout<<"params.station_spacing : "<<detector->params.station_spacing<<endl;

cout<<"Spectra energy values : "<<spectra->GetE_bin()<<endl;
for (int i=0;i<spectra->GetE_bin();i++) {
    cout<<"energy bin "<<i<<" : "<<spectra->energy[i]<<endl;
}

cout<<"IceModel R_EARTH : "<<icemodel->R_EARTH<<endl;

/////////////////////////////////////////////



double *energy = spectra->Getenergy();
int Ebin = spectra->GetE_bin();

cout<<"\n";
for (int i=0;i<Ebin;i++) {
    cout<<"energy["<<i<<"] : "<<energy[i]<<endl;
}

///////////////////////////////////////////////////

double E = 19.0;
cout<<"\n";
cout<<"Flux at 19 is : "<<spectra->GetEdNdEdAdt(E)<<endl;

///////////////////////////////////////////////////


//--------------------------------------------------
// TSpline3 *sp1;
// sp1 = spectra->GetSEdNdEdAdt();
//-------------------------------------------------- 

///////////////////////////////////////////////////





TGraph *GEdN;
//--------------------------------------------------
// GEdN = spectra->GetGEdNdEdAdt();
//-------------------------------------------------- 
GEdN = new TGraph( spectra->GetE_bin(), spectra->Getenergy(), spectra->GetEdNdEdAdt() );

TCanvas *c2 = new TCanvas("c2","A Simple Graph Example",200,10,1000,700);
c2 -> cd();
c2->SetLogy();
GEdN->SetTitle("Neutrino flux (ESS)");
GEdN->GetHistogram()->SetXTitle("log E");
GEdN->GetHistogram()->SetYTitle("Flux EdNdEdAdt (cm^2*str*s)");
GEdN->Draw("al");

//--------------------------------------------------
// sp1->SetLineColor(2);
// sp1->Draw("c same");
//-------------------------------------------------- 
c2 -> Print("GEdN.pdf");

//////////////////////////////////////////////////

// roughly 1 deg from the south pole, (approx 100km from south pole)
// ang resolution 0.1 deg.

int ang_step = 100;
double max_ang = 10.; //1 deg lat
//--------------------------------------------------
// double max_ang = 1.; //1 deg lat
//-------------------------------------------------- 

double lat[ang_step];
double Surf[ang_step];
double surf_abv_geo[ang_step];
double ice_bot[ang_step];
double ice_bot_ex[ang_step];

double ice_bot_min = 0.;
double ice_bot_min_ex = icemodel->Surface(0.,0.);;

for (int i=0;i<ang_step;i++) {
    lat[i] = (max_ang/(double)ang_step) * (double)i;
    Surf[i] = icemodel->Surface(0.,lat[i]);
    surf_abv_geo[i] = icemodel->SurfaceAboveGeoid(0.,lat[i]);
    ice_bot[i] = surf_abv_geo[i] - icemodel->IceThickness(0.,lat[i]);
    ice_bot_ex[i] = Surf[i] - icemodel->IceThickness(0.,lat[i]);
    if (ice_bot[i] < ice_bot_min) ice_bot_min = ice_bot[i];
    if (ice_bot_ex[i] < ice_bot_min_ex) ice_bot_min_ex = ice_bot_ex[i];
}

TGraph *G_surf_abv_geo;
G_surf_abv_geo = new TGraph(ang_step, lat, surf_abv_geo);
    
TGraph *G_ice_bot;
G_ice_bot = new TGraph(ang_step, lat, ice_bot);

TCanvas *cGeo = new TCanvas("cGeo","A Simple Graph Example",200,10,2000,700);
cGeo -> Divide(2,1);
cGeo -> cd(1);
G_surf_abv_geo->SetTitle("Ice surface and Ice thickness wrt Geoid");
G_surf_abv_geo->GetHistogram()->SetMaximum( surf_abv_geo[0] + 100);
G_surf_abv_geo->GetHistogram()->SetMinimum( ice_bot_min - 100);
//--------------------------------------------------
// G_surf_abv_geo->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
//-------------------------------------------------- 
G_surf_abv_geo->SetMarkerColor(4);
G_surf_abv_geo->SetMarkerSize(2);
G_surf_abv_geo->SetMarkerStyle(20);
G_surf_abv_geo->GetHistogram()->SetXTitle("Lat (deg)");
G_surf_abv_geo->GetHistogram()->SetYTitle("Z (m)");
G_surf_abv_geo->Draw("ap");

G_ice_bot->SetMarkerColor(2);
G_ice_bot->SetMarkerSize(2);
G_ice_bot->SetMarkerStyle(21);
G_ice_bot->Draw("p");

TLegend *Leg_Geo = new TLegend(1., 0.95, 0.5,0.8);
Leg_Geo -> AddEntry(G_surf_abv_geo, "Ice surface");
Leg_Geo -> AddEntry(G_ice_bot, "Ice bottom");
Leg_Geo -> Draw();


    
cGeo -> cd(2);

TGraph *G_surf_abv_geo2;
G_surf_abv_geo2 = new TGraph(ang_step, lat, Surf);
    
TGraph *G_ice_bot2;
G_ice_bot2 = new TGraph(ang_step, lat, ice_bot_ex);

G_surf_abv_geo2->SetTitle("Ice surface and Ice thickness");
G_surf_abv_geo2->GetHistogram()->SetMaximum( Surf[0] + 100);
G_surf_abv_geo2->GetHistogram()->SetMinimum( ice_bot_min_ex - 100);
//--------------------------------------------------
// G_surf_abv_geo->GetXaxis()->SetLimits( (int)detector->stations[station_choice].GetX() - 100, (int)detector->stations[station_choice].GetX() + 100 );
//-------------------------------------------------- 
G_surf_abv_geo2->SetMarkerColor(4);
G_surf_abv_geo2->SetMarkerSize(2);
G_surf_abv_geo2->SetMarkerStyle(20);
G_surf_abv_geo2->GetHistogram()->SetXTitle("Lat (deg)");
G_surf_abv_geo2->GetHistogram()->SetYTitle("Z (m)");
G_surf_abv_geo2->Draw("ap");

G_ice_bot2->SetMarkerColor(2);
G_ice_bot2->SetMarkerSize(2);
G_ice_bot2->SetMarkerStyle(21);
G_ice_bot2->Draw("p");

//--------------------------------------------------
// TLegend *Leg_Geo2 = new TLegend(1., 0.95, 0.5,0.8);
// Leg_Geo2 -> AddEntry(G_surf_abv_geo2, "Ice surface");
// Leg_Geo2 -> AddEntry(G_ice_bot2, "Ice bottom");
// Leg_Geo2 -> Draw();
//-------------------------------------------------- 




cGeo -> Print("GEOID1.pdf");



//////////////////////////////////////////////////


int evt_n, station_n, string_n, antenna_n, ray_sol_n;

evt_n = 0;
station_n = 0;
string_n = 1;
//string_n = 0;
antenna_n = 1;
//antenna_n = 0;
ray_sol_n = 0;

AraTree2[0]->GetEvent(evt_n);

double time[settings->NFOUR/2];
double V[settings->NFOUR/2];
double V_org[settings->NFOUR/2];

/*cout<<"view angle : "<<report->stations[station_n].strings[string_n].antennas[antenna_n].view_ang[ray_sol_n]*DEGRAD<<endl;

for (int i=0;i<settings->NFOUR/2;i++) {
    time[i] = report->stations[station_n].strings[string_n].antennas[antenna_n].time[ray_sol_n][i];
    V_org[i] = report->stations[station_n].strings[string_n].antennas[antenna_n].V[ray_sol_n][i];
    if (i<settings->NFOUR/4) {
        V[i+settings->NFOUR/4] = report->stations[station_n].strings[string_n].antennas[antenna_n].V[ray_sol_n][i];
    }
    else {
        V[i-settings->NFOUR/4] = report->stations[station_n].strings[string_n].antennas[antenna_n].V[ray_sol_n][i];
    }
}


TGraph *G_V_time;
G_V_time = new TGraph(settings->NFOUR/2, time, V);

TGraph *G_V_time_org;
G_V_time_org = new TGraph(settings->NFOUR/2, time, V_org);

TCanvas *cV_time = new TCanvas("cV_time","V(t)", 200,10,1400,700);
cV_time->Divide(2,1);
cV_time -> cd(1);
G_V_time_org->SetTitle("V(t) for evt %d, station[%d].string[%d].antenna[%d] (org)");
G_V_time_org->GetHistogram()->SetXTitle("time (s)");
G_V_time_org->GetHistogram()->SetYTitle("Voltage (V)");
G_V_time_org->Draw("al");

cV_time -> cd(2);
G_V_time->SetTitle("V(t) for evt 0, station[0].string[0].antenna[0] (fixed?)");
G_V_time->GetHistogram()->SetXTitle("time (s)");
G_V_time->GetHistogram()->SetYTitle("Voltage (V)");
G_V_time->Draw("al");

cV_time -> Print("V_time_evt0.pdf");


TGraph *G_V_time_zoom;
G_V_time_zoom = new TGraph(settings->NFOUR/2, time, V);

TCanvas *cV_time_zoom = new TCanvas("cV_time_zoom","V(t) for evt 0, station0.string0.antenna0",200,10,1400,700);
cV_time_zoom->Divide(2,1);
cV_time_zoom -> cd(1);
G_V_time->Draw("al");

cV_time_zoom -> cd(2);
G_V_time_zoom->SetTitle("V(t) for evt 0, station[0].string[0].antenna[0] (Zoomed)");
G_V_time_zoom->GetHistogram()->SetXTitle("time (s)");
G_V_time_zoom->GetHistogram()->SetYTitle("Voltage (V)");
G_V_time_zoom->GetXaxis()->SetLimits(8E-8,12E-8);
G_V_time_zoom->Draw("al");

cV_time_zoom -> Print("V_time_evt0_zoom.pdf");



int N_freq;
N_freq = detector->GetFreqBin();

double freq[N_freq];
double vmmhz_freq[N_freq];



for (int i=0;i<N_freq;i++) {
    freq[i] = detector->GetFreq(i); // freq in Hz
    vmmhz_freq[i] = report->stations[station_n].strings[string_n].antennas[antenna_n].vmmhz[ray_sol_n][i];
}


TGraph *G_vmmhz_freq;
G_vmmhz_freq = new TGraph(N_freq, freq, vmmhz_freq);


TCanvas *cVmMHz = new TCanvas("cVmMHz","V/m/MHz", 200,10,1000,700);

cVmMHz -> cd();
cVmMHz -> SetLogy();

G_vmmhz_freq->SetTitle("VmMHz for evt 0, station[0].string[0].antenna[0] (before antenna)");
G_vmmhz_freq->GetHistogram()->SetXTitle("freq (Hz)");
G_vmmhz_freq->GetHistogram()->SetYTitle("Signal (V/m/MHz)");
G_vmmhz_freq->Draw("al");

cVmMHz -> Print("VmMHz_evt0.pdf");


AraTree2->GetEvent(3);
for (int i=0; i<detector->params.number_of_stations; i++) {
    for (int j=0; j<detector->params.number_of_strings_station; j++) {
        for (int k=0; k<detector->params.number_of_antennas_string; k++) {
            cout<<"Rank of stations["<<i<<"].strings["<<j<<"].antennas["<<k<<"] = "<<report->stations[i].strings[j].antennas[k].Rank[0]<<endl;
            cout<<"PeakV of stations["<<i<<"].strings["<<j<<"].antennas["<<k<<"] = "<<report->stations[i].strings[j].antennas[k].PeakV[0]<<endl;
        }
    }
}
*/
////////////////////////////////////////
//Histogram of the neutrino flavor.
TCanvas *cflavor = new TCanvas("cflavor","cflavor",0,0,1000,700); 
 TH1F *hflavor=new TH1F("Flavor","",3,.5,3.5);

//Histogram of log10 of PNU
TCanvas *clogpnu = new TCanvas("logpnu","logpnu",0,0,1000,700); 
 TH1F *hlogpnu=new TH1F("logpnu","",50,15,21);
 
//Histogram of the number of interactions
TCanvas *cn_interactions = new TCanvas("cn_interactions","cn_interactions",0,0,1000,700); 
 cn_interactions->Divide(2,2);
 TH1F *hn_interactions=new TH1F("n_interactions","",2,1,3);
//Histogram of entry angle theta and histogram of cosine theta
TCanvas *ctheta = new TCanvas("cvectortheta","cvectortheta",0,0,1000,700); 
 ctheta->Divide(2,2);
 TH1F *htheta=new TH1F("vectortheta","",50,-PI/2.0-0.1,PI/2.0+0.1);
 double theta_array[settings->NNU];//used when making graphs of theta and costheta
 double costheta_array[settings->NNU];

//Histogram of entry angle phi
 TH1F *hphi=new TH1F("phi","",50,-.3,PI*2.0+0.1);
 double phi_array[settings->NNU];


//makes histogram and array of the viewing angle
TCanvas *cview_ang = new TCanvas("cview_ang","cview_ang",0,0,1000,700); 
 cview_ang->Divide(2,2);
 double launch_ang_array[settings->NNU];//array storing launch angle for all events
TCanvas *cview_ang_include = new TCanvas("view_ang_include","cview_ang_include",0,0,1000,700);//Canvas for graph of receiving angle with no sorting 
TCanvas *cview_ang_2plus = new TCanvas("cview_ang_2+","cview_ang_2+",0,0,1000,700);//Canvas for graph of receiving angle with no sorting 
TH1F *hview_ang=new TH1F("viewing angle","",50,-0.1,3.6);
TH1F *hview_ang_include=new TH1F("viewing angle including no solution","",50,-0.1,3.6);;//viewing angle histogram when there's no solution
 double view_ang_array[settings->NNU];
 double view_ang_array2plus[settings->NNU];//array for events where view_ang is two or more
 double event_depth[settings->NNU];

//Here, canvases are declared for histograms of receiving angle, distance from antenna, and launch angle 
 TCanvas *crec_ang = new TCanvas("crec_ang","crec_ang",0,0,1000,700);//Canvas for graph of receiving angle with no sorting 
 crec_ang->Divide(2,2);
 TCanvas *crec_ang_focus = new TCanvas("crec_ang_focus","crec_ang_focus",0,0,1000,700);//Canvas for graph of receiving angle with no sorting 
 crec_ang_focus->Divide(2,2); 
 TCanvas *crec_ang_super_focus = new TCanvas("crec_ang_focus","crec_ang_focus",0,0,1000,700);//Canvas for graph of receiving angle with no sorting 
crec_ang_super_focus->Divide(2,2); 
 TCanvas *crec_ang2 = new TCanvas("crec_ang2","crec_ang2",0,0,1000,700);
 crec_ang2->Divide (1,2); 
TCanvas *cdist = new TCanvas("cdist","cdist",0,0,1000,700);//Canvas for graph of distances 
 cdist->Divide(2,2); 
TCanvas *cdist_fern = new TCanvas("cdist_fern","cdist_fern",0,0,1000,700);//Canvas for graph of distances in relation to fern
 cdist_fern->Divide(2,2); 
 TCanvas *claunch= new TCanvas("claunch","claunch",0,0,1000,700);//Canvas for graph of distances 
 claunch->Divide(2,2);
TCanvas *cdist_sol = new TCanvas("cdist_sol","cdist_sol",0,0,1000,700);//Canvas for graph of distance, 1st and second solution
 cdist_sol->Divide(2,2);

 //canvases for histograms sorted by antenna strength

 TCanvas *cview_ang_multi = new TCanvas("cview_ang_multi","cview_ang_multi",0,0,1000,700); 
 TCanvas *crec_ang_multi = new TCanvas("crec_ang_multi","crec_ang_multi",0,0,1000,700); 
 TCanvas *creflect_ang_multi = new TCanvas("creflect_ang_multi","creflect_ang_multi",0,0,1000,700); 
 TCanvas *cdist_sorted = new TCanvas("cdist_sorted","cdist_sorted",0,0,1000,700); 

//Creates canvas for everything related to reflection angle
 TCanvas *creflect_ang = new TCanvas("creflect_ang","creflect_ang",0,0,1000,700); 
 creflect_ang->Divide(2,2);
TCanvas *crec_ang_not_hundred = new TCanvas("crec_ang_not_hundred","crec_ang_not_hundred",0,0,1000,700);  
TCanvas *cinteractions_1 = new TCanvas("cinteractions_1","cinteractions_1",0,0,1000,700);  
 cinteractions_1->Divide(2,2);
TCanvas *cinteractions_2 = new TCanvas("cinteractions_2","cinteractions_2",0,0,1000,700);  
 cinteractions_2->Divide(2,2);
TCanvas *cinteractions_3 = new TCanvas("cinteractions_3","cinteractions_3",0,0,1000,700);  
 cinteractions_3->Divide(2,2);
TCanvas *cinteractions_4 = new TCanvas("cinteractions_3","cinteractions_3",0,0,1000,700);  
 cinteractions_4->Divide(2,2);
TCanvas *cpol= new TCanvas("cpol","cpol",0,0,1000,700);  
 cpol->Divide(2,2);
TCanvas *cpol2= new TCanvas("cpol2","cpol2",0,0,1000,700);  
 cpol2->Divide(2,2);
TCanvas *cpol3= new TCanvas("cpol2","cpol2",0,0,1000,700);  
 cpol3->Divide(2,2);
 //Histogram declarations

TH1F *hrec_ang=new TH1F("receiving angle","",100,-6.5,6.5);//histogram for unsorted receiving angle
TH1F *hrec_ang_focus=new TH1F("receiving angle focus","",25,1.5,1.65);//histogram for unsorted receiving angle
TH1F *hrec_ang_super_focus_H=new TH1F("receiving angle super focus","",25,1.56,1.58);//histogram for unsorted receiving angle
TH1F *hrec_ang_super_focus_V=new TH1F("receiving angle super focus","",25,1.56,1.58);//histogram for unsorted receiving angle
TH1F *hrec_ang_super_focus_1st_H=new TH1F("receiving angle super focus 1st","",25,1.56,1.58);//histogram for unsorted receiving angle 
TH1F *hrec_ang_super_focus_1st_V=new TH1F("receiving angle super focus 1st","",25,1.56,1.58);//histogram for unsorted receiving angle 
TH1F *hrec_ang_second_antenna=new TH1F("receiving angle2nd antenna","",100,-.5,3.1);//wide histogram for second antenna of receiving angle
TH1F *hrec_ang_focus_second_antenna=new TH1F("receiving angle zoom 2nd antenna","",25,1.5,1.65);//histogram for second antenna receiving angle, zoomed in
TH1F *hdist=new TH1F("distance","",100,0.,3600);//histogram for unsorted distance
 TH1F *hdistweighted1=new TH1F("distance weighted1","",100,0.,3600);//histogram for unsorted distance modeled as circles of pie
 TH1F *hdistweighted2=new TH1F("distance weighted2","",30,0.,3600);//histogram for unsorted distance modeled as shells of sphere
 TH1F *hdist_fern=new TH1F("distance, in fern", " ", 100, 0., 3600);
 TH1F *hdist_not_fern=new TH1F("distance, not in fern", "", 100, 0., 3600);
 TH1F *hlaunch=new TH1F("launch angle", "", 100, 0., 2.0); //histogram for unsorted launch angle
//Reflection histograms declaration
 TH1F *hreflect_ang=new TH1F("reflect angle","",50,-2.0,-1.3);//histogram for only when the reflection angle exists
 TH1F *hrec_ang_not_hundred=new TH1F("receiving angle isn't 100","",100,2.5,2.6); //histogram showing receiving angle when reflection occurs
 TH1F *hdist_second_solution=new TH1F("dist_second_solution","",100,0,3600);

 //Some arrays we'll be using (some of them are probably irrelevant now)
 double rec_ang_array[settings->NNU];
 double rec_ang_twohalf_array[settings->NNU];
 int count_rec=0;
 double rec_ang_average=0;
 double reflect_ang_array[settings->NNU]; //array of reflection angles
 double rec_ang_not_hundred[settings->NNU]; //array of receiving angles when reflection occurs
 int count_reflect=0;

 double launch_ang_fern[settings->NNU];//array of launch angles when one is in fern
 double launch_ang_not_fern[settings->NNU];//array of launch angles when not in fern
 double rec_ang_fern[settings->NNU];//array of receiving angles when in fern
 double rec_ang_not_fern[settings->NNU];//array of receiving angles when not in fern
 //Make Zoomed versions of these last four arrays
 double launch_ang_zoom_fern[settings->NNU];//array of launch angles when one is in fern
 double launch_ang_zoom_not_fern[settings->NNU];//array of launch angles when not in fern
 double rec_ang_zoom_fern[settings->NNU];//array of receiving angles when in fern
 double rec_ang_zoom_not_fern[settings->NNU];//array of receiving angles when not in fern
 //more arrays
 double dist_array[settings->NNU];//array of distance for all events
 double  dist_fern_array[settings->NNU];//array of distance only in fern




 //Histograms to be filled as we loop through antennae, where they are filled according to antenna strength (i.e. the antenna with the strongest signal has all of its data go into the histograms labeled 1, second strongest in 2, etc.) 
 TH1F *hview_ang1=new TH1F("viewing angle","",30,-0.1,3.5);
 TH1F *hview_ang2=new TH1F("viewing angle","",30,-0.1,3.5);
 TH1F *hview_ang3=new TH1F("viewing angle","",30,-0.1,3.5);
 TH1F *hview_ang4=new TH1F("viewing angle","",30,-0.1,3.5);
 
 TH1F *hrec_ang1=new TH1F("receiving angle1","",25,-.3,3.5);
 TH1F *hrec_ang2=new TH1F("receiving angle2","",25,-.3,3.5);
 TH1F *hrec_ang3=new TH1F("receiving angle3","",25,-.3,3.5);
 TH1F *hrec_ang4=new TH1F("receiving angle4","",25,-.3,3.5);
 
 TH1F *hrec_ang_focus1=new TH1F("receiving angle, zoomed in1","", 25,1.5,1.65);
 TH1F *hrec_ang_focus2=new TH1F("receiving angle, zoomed in2","", 25,1.5,1.65);
 TH1F *hrec_ang_focus3=new TH1F("receiving angle, zoomed in3","", 25,1.5,1.65);
 TH1F *hrec_ang_focus4=new TH1F("receiving angle, zoomed in4","", 25,1.5,1.65);

TH1F *hreflect_ang1=new TH1F("rerflecting angle1","",100,-1.9,-1.5);
TH1F *hreflect_ang2=new TH1F("reflecting angle2","",100,-1.9,-1.5);
TH1F *hreflect_ang3=new TH1F("reflecting angle3","",100,-1.9,-1.5);
 TH1F *hreflect_ang4=new TH1F("reflecting angle4","",100,-1.9,-1.5);

 TH1F *hdist1=new TH1F("disance sorted1","", 25, 0,3500);
 TH1F *hdist2=new TH1F("disance sorted2","", 25, 0,3500);
 TH1F *hdist3=new TH1F("disance sorted3","", 25, 0,3500);
 TH1F *hdist4=new TH1F("disance sorted4","", 25, 0,3500);

 TH1F *hlaunch1=new TH1F("launch angle1", "", 25, 0., 2.0); //histogram for sorted launch angle, strongest signal
 TH1F *hlaunch2=new TH1F("launch angle2", "", 25, 0., 2.0); //histogram for sorted launch angle, 2nd strongest
 TH1F *hlaunch3=new TH1F("launch angle3", "", 25, 0., 2.0); //histogram for sorted launch angle, 3rd strongest
 TH1F *hlaunch4=new TH1F("launch angle4", "", 25, 0., 2.0); //histogram for sorted launch angle, 4th strongest

 //Interactions histograms
 TH1F  *hnu_nubar=new TH1F("nu_nubar","", 100, -2000000000, 2000000000);//histrogram for nu_nubar
 TH1F *hemfrac=new TH1F("emfrac","", 100, 0.1, 1.0);//histogram for emfrac
 TH1F *hemfrac_complement=new TH1F("emfrac_complement", "", 100, -7.0, .1);//histogram for log(1-emfrac)
 TH1F *hhadfrac=new TH1F("hadfrac","", 100, -6.0,0.1);//logarithmic histogram for hadfrac
 TH1F *hhad_em_sum=new TH1F("hhad_em_sum","", 100, 0, 1.1); //histogram of hadfrac+emfrac
 TH1F *hrst=new TH1F("rst","", 100, -3000000000, 3000000000);
 TH1F *hvmmhz1m_tmp=new TH1F("vmmhz1m_tmp","",100, -9.0, -6.5);
 TH1F *hvmmhz1m=new TH1F("vmmhz1m_tmp","", 100, 0,1.0);
 TH1F *hvmmhz1m_log=new TH1F("vmmhz1m_tmp_log","", 100, -7.0,0.1);
 TH1F *hd_theta_em=new TH1F("d_theta_em","", 100, 0.0, .15);
 TH1F *hd_theta_had=new TH1F("d_theta_had","", 100,0.045,0.0485);

 //Histogram for Polarization
 TH1F *hpol_17=new TH1F("pol_17","", 25, -1, 1);
TH1F *hpol_17_trig=new TH1F("pol_17_trig","", 25, -1, 1);
 TH1F *hpol_18=new TH1F("pol_18","", 25, -1, 1);
TH1F *hpol_18_trig=new TH1F("pol_18_trig","", 25, -1, 1);

 TH1F *hpol_factor=new TH1F("pol_factor","",25, -.5,1.5);
 TH1F *hpol_factor_trig=new TH1F("pol_factor_trig","",25, -.5,1.5);

 TH1F *hpolz_all=new TH1F("polz_all","",50,0,1);
 TH1F *hpolz_trig=new TH1F("polz_trig","",50,0,1);
 TH1F *hpol_flat_all=new TH1F("pol_flat_all","",50,0,1);
 TH1F *hpol_flat_trig=new TH1F("pol_flat_trig","",50,0,1);





 //We're interested in finding the maximum difference between viewing angles. Here we declare variables.
 double view; 
double maxview=-999;
 double minview=999;
 double maxdiff=0;
 double launch_debug;
 double maxlaunch=-999;//for debugging launch values
 double rec_ang;//just used to store long expression

//////////////////////////////////////////
 //Here we investigate the spike in the launch angle graph
   TCanvas *claunch_compare = new TCanvas("claunch_compare","claunch_compare",0,0,1000,700);//Canvas for all launch angle graphs
   claunch_compare->Divide(2,2);   
 TCanvas *crec_launch_fern = new TCanvas("crec_launch_fern", "crec_launch_fern",0,0,1000,700);//graph of receiving angle and launch angle in fern/not fern
 crec_launch_fern->Divide(2,2);
//Declare all the arrays we'll be using to make graphs
   double launch_ang_narrow[settings->NNU];
   double view_launch[settings->NNU];
   double rec_launch[settings->NNU];
   double dist_launch[settings->NNU];
   double reflect_launch[settings->NNU];
   double reflect_launch_not_hundred[settings->NNU];
   //These booleans are used to mark whether we've grabbed the interesting events in the dist-launch graph yet.   
bool peak_met=false;
   bool left_met=false;
   bool right_met=false;
   int progress_counter=0;
   double min_time=0;
   double max_time=0;
   double diff_time[settings->NNU];
   double diff_time_sub_hundred[settings->NNU];
   double diff_time_no_ref[settings->NNU];
   double cos_pol_ang[settings->NNU];
   //This is to note the size of the arrays that are divided up into different frequencies, such as vmmhz1m
   const int FREQ_SIZE=60;
   /////////////////////////////////////////
   //Get ready for file writing
   ofstream output;
   output.open("data.txt");
   ofstream vector_pair;
   vector_pair.open("math.txt");
   output << "Posnu          Rec_Ang       Launch_Ang View_Ang Dist      Dept   Rec_Ang\n";
   ofstream rec_pair_out;
   rec_pair_out.open("rec_pair.txt");
   ofstream launch2d;
   launch2d.open("vector2d.txt");
   launch2d << "Export[" << '"' << "2dvector.pdf" << '"' << ", {";
   /////////////////////////////////////////

cout<<"start loop"<<endl;
AraTree2[0]->GetEvent(0);
cout<<"nnu theta : "<<event->nnu.Theta()<<endl;

for(int i =0; i<settings->NNU; i++)
  {
  
  
  
    //if((100*i)/(settings->NNU)>progress_counter)
//	{
  cout  << "Progress is: " << (100*i)/(settings->NNU) << "%\n";
  progress_counter=(100*i)/(settings->NNU);
//	}
    AraTree2[0]->GetEvent(i);

    hflavor->Fill(event->nuflavorint);

    hlogpnu->Fill(log10(event->pnu));
    hn_interactions->Fill(event->n_interactions);

    theta_array[i]=event->nnu.Theta(); //makes array of theta for later use
    costheta_array[i]=cos(event->nnu.Theta()); //makes array of costheta for later use
    htheta->Fill(cos(event->nnu.Theta()));

    hphi->Fill(event->nnu.Phi());
    phi_array[i]=event->nnu.Phi(); //puts phi in array for later use
    if(report->stations[0].strings[0].antennas[0].ray_sol_cnt !=0)//This was commented out in earlier version for reasons (possibly good?) unknown
     {
     dist_array[i]=report->stations[0].strings[0].antennas[0].Dist[0];//fills distance array for all events
    launch_ang_array[i]=report->stations[0].strings[0].antennas[0].launch_ang[0]; //fills launch angle array
     }
   //view angle stuff
   if(report->stations[0].strings[0].antennas[0].ray_sol_cnt !=0)
     {
       hview_ang_include->Fill(report->stations[0].strings[0].antennas[0].view_ang[0]);
       hview_ang->Fill(report->stations[0].strings[0].antennas[0].view_ang[0]);
       view_ang_array[i]=report->stations[0].strings[0].antennas[0].view_ang[0];
       if(view_ang_array[i]>2)
	 {       
view_ang_array2plus[i]=report->stations[0].strings[0].antennas[0].view_ang[0];
 event_depth[i]=icemodel->Surface(event->Nu_Interaction[0].posnu)-event->Nu_Interaction[0].posnu.R();
//event_depth[i]=0;
 	 }//if 2+
       else
	 {
	   view_ang_array2plus[i]=0;
	   event_depth[i]=0;
	 }//else


       //receiving angle stuff
       rec_ang=report->stations[0].strings[0].antennas[0].rec_ang[0];
       hrec_ang->Fill(report->stations[0].strings[0].antennas[0].rec_ang[0]); //fill receving angle histogram
       hrec_ang_focus->Fill(report->stations[0].strings[0].antennas[0].rec_ang[0]); //I think this fills a histogram that doesn't actually get used...
       hrec_ang_super_focus_V->Fill(report->stations[0].strings[0].antennas[0].rec_ang[0]);//This fills the first numbered antennas rec ang into a histogram 1.56-1.58
       hrec_ang_super_focus_H->Fill(report->stations[0].strings[0].antennas[0].rec_ang[0]);//This fills the first numbered antennas rec ang into a histogram 1.56-1.58
       
rec_ang_array[i]=rec_ang; //Fills  receiving angle array that was used for debugging
       //   cout << rec_ang << "\n";
       if (rec_ang>2.5)//used to study spike at receiving angle 2.5
	 {
	   count_rec++;
	   rec_ang_twohalf_array[i]=rec_ang;
	   rec_ang_average+=rec_ang;
	 }
       else
	 {
	   rec_ang_twohalf_array[i]=0;
	 }
       //Fill histogram for second (numerical) antenna
   if(report->stations[0].strings[0].antennas[1].ray_sol_cnt !=0)
     {
       hrec_ang_second_antenna->Fill(report->stations[0].strings[0].antennas[1].rec_ang[0]);
      hrec_ang_focus_second_antenna->Fill(report->stations[0].strings[0].antennas[1].rec_ang[0]);
     } //if

   //Reflection histograms, filling
   if(report->stations[0].strings[0].antennas[0].ray_sol_cnt !=0 && report->stations[0].strings[0].antennas[0].reflect_ang[0] !=100) //if reflection occurs
     {
       hreflect_ang->Fill(report->stations[0].strings[0].antennas[0].reflect_ang[0]);
       reflect_ang_array[i]=report->stations[0].strings[0].antennas[0].reflect_ang[0];
       //cout << report->stations[0].strings[0].antennas[0].reflect_ang[0]<< "\n";
       count_reflect++;
       rec_ang_not_hundred[i]= report->stations[0].strings[0].antennas[0].rec_ang[0];
       hrec_ang_not_hundred->Fill(report->stations[0].strings[0].antennas[0].rec_ang[0]);
     }//if
  else
    {
      reflect_ang_array[i]=0;
      rec_ang_not_hundred[i]=0;
    }//else


 //fill distance histograms
       int dist=report->stations[0].strings[0].antennas[0].Dist[0];
       //cout<< " \n\n\n Distance is: " << dist << "\n\n\n";
   hdist->Fill(report->stations[0].strings[0].antennas[0].Dist[0]);
   hdistweighted1->Fill(report->stations[0].strings[0].antennas[0].Dist[0], 1/report->stations[0].strings[0].antennas[0].Dist[0]);
   hdistweighted2->Fill(report->stations[0].strings[0].antennas[0].Dist[0], 1/(report->stations[0].strings[0].antennas[0].Dist[0]*report->stations[0].strings[0].antennas[0].Dist[0]*report->stations[0].strings[0].antennas[0].Dist[0]*report->stations[0].strings[0].antennas[0].Dist[0]));
      hdist_second_solution->Fill(report->stations[0].strings[0].antennas[0].Dist[1]);
//if statements regarding whether events are in fern 
      if(icemodel->Surface(event->Nu_Interaction[0].posnu)-event->Nu_Interaction[0].posnu.R()<150)
     {
       //       cout << "\n filling in fern bin";
hdist_fern->Fill(report->stations[0].strings[0].antennas[0].Dist[0]);
 launch_ang_fern[i]=launch_ang_array[i];
 rec_ang_fern[i]=rec_ang_array[i];
 dist_fern_array[i]=dist_array[i];
     }
   else
     {
       //  cout << "\n filling out of fern bin";
hdist_not_fern->Fill(report->stations[0].strings[0].antennas[0].Dist[0]);
 launch_ang_not_fern[i]=launch_ang_array[i];
 rec_ang_not_fern[i]=rec_ang_array[i];
     }
   //reset variables for finding difference in viewing angle
   maxview=-999;
   minview=999;

  

  //add to launch angle histogram
   launch_debug=report->stations[0].strings[0].antennas[0].launch_ang[0];
   hlaunch->Fill(launch_debug);//fills launch angle histogram
   if (launch_debug>maxlaunch)
     {
       maxlaunch=launch_debug;
     }
   //   cout<<"Launch angle is: " << launch_debug << " and max so far is " << maxlaunch << "\n";

//For loop running through all the antennas of an event
// commented out for performance   
       for (int k=0; k< report->stations[0].strings.size();k++)
	 {
  
     for (int j=0; j<report->stations[0].strings[k].antennas.size(); j++)
       {
	 view=report->stations[0].strings[k].antennas[j].view_ang[0];
	 if(view>maxview)
	   {
	     maxview=view;
	   }//if
	 if(view<minview)
	   {
	     minview=view;
	   }//if
       //check to see if it's rank 1-4
	 //////////////////////////////
	 if (report->stations[0].strings[k].antennas[j].Rank[0]==1) //Fills each histogram corresponding to strongest antenna signal
	 {

	   hview_ang1->Fill(report->stations[0].strings[k].antennas[j].view_ang[0]);
	   hrec_ang1->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
	   hrec_ang_focus1->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
	   if (j % 2 ==0)
	     {
	       hrec_ang_super_focus_1st_V->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
	     }
	   else
	     {
	       hrec_ang_super_focus_1st_H->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
	     }
	   double reflect_ang=report->stations[0].strings[k].antennas[j].reflect_ang[0];
	   if(reflect_ang!=100.)
	     {
	       hreflect_ang1->Fill(reflect_ang);
	     }
	   hdist1->Fill(report->stations[0].strings[k].antennas[j].Dist[0]);
   hlaunch1->Fill(report->stations[0].strings[k].antennas[j].launch_ang[0]);
	 }
     ///////////////////////////////
	 if (report->stations[0].strings[k].antennas[j].Rank[0]==2)//fills histogram for 2nd strongest signal
	 {

	   hview_ang2->Fill(report->stations[0].strings[k].antennas[j].view_ang[0]);
	   hrec_ang2->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
hrec_ang_focus2->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
double reflect_ang=report->stations[0].strings[k].antennas[j].reflect_ang[0];
	   if(reflect_ang!=100.)
	     {
	       hreflect_ang2->Fill(reflect_ang);
	     }
	   hdist2->Fill(report->stations[0].strings[k].antennas[j].Dist[0]);	
   hlaunch2->Fill(report->stations[0].strings[k].antennas[j].launch_ang[0]);
 }

       ////////////////////////////////////
	 if (report->stations[0].strings[k].antennas[j].Rank[0]==3)//fills for third strongest signal
	 {

	   hview_ang3->Fill(report->stations[0].strings[k].antennas[j].view_ang[0]);
	   hrec_ang3->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
	   hrec_ang_focus3->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
double reflect_ang=report->stations[0].strings[k].antennas[j].reflect_ang[0];
	   if(reflect_ang!=100.)
	     {
	       hreflect_ang3->Fill(reflect_ang);
	     }
	   hdist3->Fill(report->stations[0].strings[k].antennas[j].Dist[0]);
   hlaunch3->Fill(report->stations[0].strings[k].antennas[j].launch_ang[0]);	 
}
       ///////////////////////////////////////
	 if (report->stations[0].strings[k].antennas[j].Rank[0]==4)//fills 4th strongest historam
	 {

	   hview_ang4->Fill(report->stations[0].strings[k].antennas[j].view_ang[0]);
	   hrec_ang4->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
	   hrec_ang_focus4->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
double reflect_ang=report->stations[0].strings[k].antennas[j].reflect_ang[0];
	   if(reflect_ang!=100.)
	     {
	       hreflect_ang4->Fill(reflect_ang);
	     }
	   hdist4->Fill(report->stations[0].strings[k].antennas[j].Dist[0]);
   hlaunch4->Fill(report->stations[0].strings[k].antennas[j].launch_ang[0]);
	 }//if rank=4
			    
       }       //for j
	 } // for k 
       //check to see  if  the view difference is maximum
       //       cout << "Max viewang is "  << maxview << "     Min view ang is " << minview << "     Difference is " << maxview-minview << "\n";	 
       if(maxview-minview>maxdiff)
	   {
	     maxdiff=maxview-minview;
	   }//if

    

       //Fill launch comparison junk
     if (report->stations[0].strings[0].antennas[0].launch_ang[0]>1.3 and report->stations[0].strings[0].antennas[0].launch_ang[0]<1.50)//This constricts it to the range where we see the bump in launch angle
       {
	 launch_ang_narrow[i]=launch_ang_array[i];
	 view_launch[i]=report->stations[0].strings[0].antennas[0].view_ang[0];
	 rec_launch[i]=report->stations[0].strings[0].antennas[0].rec_ang[0];
	
	 dist_launch[i]=report->stations[0].strings[0].antennas[0].Dist[0];
	 reflect_launch[i]=report->stations[0].strings[0].antennas[0].reflect_ang[0];
	 if(report->stations[0].strings[0].antennas[0].reflect_ang[0]!=100)
	   {
	     reflect_launch_not_hundred[i]=report->stations[0].strings[0].antennas[0].reflect_ang[0];
	   }//if reflect not hundred
	 else
	   {
reflect_launch_not_hundred[i]=0;
	   }//else 
	 launch_ang_zoom_fern[i]=launch_ang_fern[i];
	 if (launch_ang_zoom_fern[i]==0)
	   {
	     launch_ang_zoom_fern[i]=1.29;
	   }
	 rec_ang_zoom_fern[i]=rec_ang_fern[i];
	 launch_ang_zoom_not_fern[i]=launch_ang_not_fern[i];
	 if(launch_ang_zoom_not_fern[i]==0)
	   {
	     launch_ang_zoom_not_fern[i]=1.29;
	       }

	 rec_ang_zoom_not_fern[i]=rec_ang_not_fern[i];
       }//if launch is in range
     else
       {
	 launch_ang_narrow[i]= 1.29;
	 view_launch[i]=0;
	 rec_launch[i]=0;
	 dist_launch[i]=0;
	 reflect_launch[i]=0;
	 reflect_launch_not_hundred[i]=0;
	 launch_ang_zoom_fern[i]=1.29;
	 launch_ang_zoom_not_fern[i]=1.29;
       }//else

     }//if solution exists
  else
    {
      view_ang_array[i]=0;
      view_ang_array2plus[i]=0;
      event_depth[i]=0;
      rec_ang_array[i]=0;

	 launch_ang_narrow[i]= 1.29;
	 view_launch[i]=0;
	 rec_launch[i]=0;
	 dist_launch[i]=0;
	 reflect_launch[i]=0;
	 reflect_launch_not_hundred[i]=0;
    }//else
   /////////////////////////////////
   //Do Interaction class stuff

   int nu_nubar=event->Nu_Interaction[0].nu_nubar;
   //cout<< "\nNu_nubar is: " <<nu_nubar;
   hnu_nubar->Fill(nu_nubar);
   double emfrac=event->Nu_Interaction[0].emfrac;
   //cout << "\nemfrac is " << emfrac;
   hemfrac->Fill(emfrac);
   hemfrac_complement->Fill(log10(1-emfrac));
   double hadfrac=event->Nu_Interaction[0].hadfrac;
   //cout << "\nhadfrac is " << hadfrac;
   hhadfrac->Fill(log10(hadfrac));
   hhad_em_sum->Fill(emfrac+hadfrac);
   double elast_y=event->Nu_Interaction[0].elast_y;
   //cout << "\nelast_y is : " << elast_y;
   //helast_y->Fill(elast_y);
   int rst=event->Nu_Interaction[0].ray_sol_cnt;
   //cout << "\nrst is: " << rst;
   if(rst!=0)
     {
       hrst->Fill(rst);
     }
   double vmmhz1m_tmp=event->Nu_Interaction[0].vmmhz1m_tmp;
  
   hvmmhz1m_tmp->Fill(log10(vmmhz1m_tmp));

   ////////////////////////////////////////
   //for loops filling the vectors that have hundreds of components
double vmmhz1m_avg=0;   
 for(int i=0; i<event->Nu_Interaction[0].vmmhz1m.size()-1; i++)
     
     {
       vmmhz1m_avg+=event->Nu_Interaction[0].vmmhz1m[i];
     }
 vmmhz1m_avg=vmmhz1m_avg/(event->Nu_Interaction[0].vmmhz1m.size());

 hvmmhz1m->Fill(vmmhz1m_avg);
 hvmmhz1m_log->Fill(log10(vmmhz1m_avg));

 double vmmhz1m_em_avg=0;   
 for(int i=0; i<event->Nu_Interaction[0].vmmhz1m_em.size()-1; i++)
     
     {
       vmmhz1m_em_avg+=event->Nu_Interaction[0].vmmhz1m_em[i];
     }
 vmmhz1m_em_avg=vmmhz1m_em_avg/(event->Nu_Interaction[0].vmmhz1m_em.size());


double d_theta_em_avg=0;   
 for(int i=0; i<event->Nu_Interaction[0].d_theta_em.size()-1; i++)
     
     {
       d_theta_em_avg+=event->Nu_Interaction[0].d_theta_em[i];
     }
 d_theta_em_avg=d_theta_em_avg/(event->Nu_Interaction[0].d_theta_em.size());

 hd_theta_em->Fill(d_theta_em_avg);
 

double d_theta_had_avg=0;   
 for(int i=0; i<event->Nu_Interaction[0].d_theta_had.size()-1; i++)
     
     {
       d_theta_had_avg+=event->Nu_Interaction[0].d_theta_had[i];
     }//for loop
 d_theta_had_avg=d_theta_had_avg/(event->Nu_Interaction[0].d_theta_had.size());
 //cout << "\nd_theta_had_avg is: " <<d_theta_had_avg;
 hd_theta_had->Fill(d_theta_had_avg);
 //Output all the data for a few events
 if (i<30)
   {
     output <<  event->pnu <<  "   " << report->stations[0].strings[0].antennas[0].rec_ang[0] << "   " << report->stations[0].strings[0].antennas[0].launch_ang[0] << "   " << report->stations[0].strings[0].antennas[0].view_ang[0] << "   " << report->stations[0].strings[0].antennas[0].Dist[0] << "   " << icemodel->Surface(event->Nu_Interaction[0].posnu)-event->Nu_Interaction[0].posnu.R()<< "   " << report->stations[0].strings[0].antennas[0].reflect_ang[0] <<"\n";
   }
 if((launch_ang_array[i]< .7 && launch_ang_array[i]>.6) && peak_met==false)
   {
     output << "\nThis event is at peak of launch-dist graph\n";
  output <<  event->pnu <<  "   " << report->stations[0].strings[0].antennas[0].rec_ang[0] << "   " << report->stations[0].strings[0].antennas[0].launch_ang[0] << "   " << report->stations[0].strings[0].antennas[0].view_ang[0] << "   " << report->stations[0].strings[0].antennas[0].Dist[0] << "   " << icemodel->Surface(event->Nu_Interaction[0].posnu)-event->Nu_Interaction[0].posnu.R()<< "   " << report->stations[0].strings[0].antennas[0].reflect_ang[0] <<"\n\n";
  peak_met=true;
   }//if in range
 if((launch_ang_array[i]< 1.4 && launch_ang_array[i]>1.3) && left_met==false)
   {
     output << "\nThis event is at left part of gap in launch-dist graph\n";
     output <<  event->pnu <<  "   " << report->stations[0].strings[0].antennas[0].rec_ang[0] << "   " << report->stations[0].strings[0].antennas[0].launch_ang[0] << "   " << report->stations[0].strings[0].antennas[0].view_ang[0] << "   " << report->stations[0].strings[0].antennas[0].Dist[0] << "   " << icemodel->Surface(event->Nu_Interaction[0].posnu)-event->Nu_Interaction[0].posnu.R()<< "   " << report->stations[0].strings[0].antennas[0].reflect_ang[0] << "\n\n";
  left_met=true;
   }
 if((launch_ang_array[i]< 1.68 && launch_ang_array[i]>1.55) && right_met==false && dist_array[i]>1500)
   {
     output << "\nThis event is at right part of gap in launch-dist graph\n";
     output <<  event->pnu <<  "   " << report->stations[0].strings[0].antennas[0].rec_ang[0] << "   " << report->stations[0].strings[0].antennas[0].launch_ang[0] << "   " << report->stations[0].strings[0].antennas[0].view_ang[0] << "   " << report->stations[0].strings[0].antennas[0].Dist[0] << "   " << icemodel->Surface(event->Nu_Interaction[0].posnu)-event->Nu_Interaction[0].posnu.R()<< "   " << report->stations[0].strings[0].antennas[0].reflect_ang[0] << "\n\n";
  right_met=true;
   }
  if (report->stations[0].strings[0].antennas[0].ray_sol_cnt !=0)
 {
   //cout<<"raysolcnt : "<<report->stations[0].strings[0].antennas[0].ray_sol_cnt<<endl;
 min_time=0;
  max_time=0;

// cout << "Pol_vector theta is: " <<report->stations[0].strings[0].antennas[0].Pol_vector[0].Theta() << "\n\n";
// cout << "Energy is: " << event->pnu << "\n\n";
 hpol_17->Fill(cos(report->stations[0].strings[0].antennas[0].Pol_vector[0].Theta()), event->Nu_Interaction[0].weight);

 //do stuff with z polarization and xy polarization, triggered vs. nontriggered.
 double polz=abs(report->stations[0].strings[0].antennas[0].Pol_vector[0].GetZ()); 
 double pol_flat=sqrt(report->stations[0].strings[0].antennas[0].Pol_vector[0].GetX()*report->stations[0].strings[0].antennas[0].Pol_vector[0].GetX()+report->stations[0].strings[0].antennas[0].Pol_vector[0].GetY()*report->stations[0].strings[0].antennas[0].Pol_vector[0].GetY());
 hpol_factor->Fill(report->stations[0].strings[0].antennas[0].Pol_factor[0]);
 hpolz_all->Fill(polz,  event->Nu_Interaction[0].weight);
 hpol_flat_all->Fill(pol_flat,  event->Nu_Interaction[0].weight);
  if (report->stations[0].Global_Pass!=0)
 {
   hpol_factor_trig->Fill(report->stations[0].strings[0].antennas[0].Pol_factor[0],  event->Nu_Interaction[0].weight);
     hpol_17_trig-> Fill(cos(report->stations[0].strings[0].antennas[0].Pol_vector[0].Theta()), event->Nu_Interaction[0].weight);
     hpolz_trig->Fill(polz,  event->Nu_Interaction[0].weight);
     hpol_flat_trig->Fill(pol_flat,  event->Nu_Interaction[0].weight);
 cout << "Global pass";
 }
 
  //Here we find the difference between the maximum and minimum times of receiving the event, time and distance being equivalent. We then compare that to polarization angle
  for (int k=0; k< report->stations[0].strings.size();k++)
	 {
     for (int j=0; j<report->stations[0].strings[k].antennas.size(); j++)
       {
	 // cout << "\nWe are in loop, iteration k is: " << k << "j is: "<<j;
	 if (min_time> report->stations[0].strings[k].antennas[j].Dist[0] || min_time==0)
	   {
	     //cout << "\nLowering min.";
	     min_time=report->stations[0].strings[k].antennas[j].Dist[0];
	   }
	 if (max_time< report->stations[0].strings[k].antennas[j].Dist[0])
	   {
	     // cout << "\nRaising max.";
	     max_time=report->stations[0].strings[k].antennas[j].Dist[0];
	   }
       }//for k
	 }//for j
  diff_time[i]=max_time-min_time;
  if (diff_time[i]<100)
    {
      cout << "\ntriggered 100 if block, diff_time is: "<< diff_time[i];
      diff_time_sub_hundred[i]=diff_time[i];
      cout << "\ndiff<100 is: " << diff_time_sub_hundred[i];
    }
  else
    {
      cout << "triggered 100 else block";
      diff_time_sub_hundred[i]=0;
    }
  cout << "\nJust passed else block, diff<100 is: " << diff_time_sub_hundred[i];
  if(reflect_ang_array[i]>50)
    {
      cout << "triggered ref if block";
      diff_time_no_ref[i]=diff_time[i];
    }
  else
    {
      cout << "triggered ref else block";
      diff_time_no_ref[i]=0;
    }

  cos_pol_ang[i]=cos(report->stations[0].strings[0].antennas[0].Pol_vector[0].Theta());
  cout << "\nentering other else block, diff<100 is: " << diff_time_sub_hundred[i];
 }//if raysol!=0

  else
    {
 min_time=0;
  max_time=0;
diff_time[i]=0;
 diff_time_sub_hundred[i]=0;
 diff_time_no_ref[i]=0;
 cos_pol_ang[i]=0;
    }//else
  cout << "Just left other else block, diff<100 is: " << diff_time_sub_hundred[i];
  // cout << "\nMin time is: "<< min_time << "\nmax time is: " << max_time << "\nDiff is: " << diff_time[i]<< "\nDiff<100 is: "<< diff_time_sub_hundred[i] << "\nDiff_no_reflect is: " << diff_time_no_ref[i];

  //Here we output data on launch angle and receiving angle used to make vector plots
  if (i<100 && report->stations[0].strings[0].antennas[0].ray_sol_cnt !=0)
    {
  double posx=event->Nu_Interaction[0].posnu.GetX();
  cout << "\nposx is: " << posx;
  double posy=event->Nu_Interaction[0].posnu.GetY();
cout << "\nposy is: " << posy;
double posz=event->Nu_Interaction[0].posnu.GetZ();
cout << "\nposz is: " << posz; 
 posz=posz-6357500;
double antx=detector->stations[0].strings[0].antennas[0].GetX();
 cout << "\nantx is: " <<detector->stations[0].strings[0].antennas[0].GetX();
 double anty=detector->stations[0].strings[0].antennas[0].GetY();
 cout << "\nanty is: " <<detector->stations[0].strings[0].antennas[0].GetY();
 double antz=detector->stations[0].strings[0].antennas[0].GetZ();
 cout << "\nantz is: " <<detector->stations[0].strings[0].antennas[0].GetZ();
 double  diffx=antx-posx;
 cout << "\ndiff x is: " << antx-posx;
 double  diffy=anty-posy;
 cout << "\ndiff y is: " << anty-posy;
 double xylen=(double) sqrt(diffx*diffx + diffy*diffy);
 cout <<"\nxylen is: " << xylen;
 cout<< "\nlaunch ang is: " << launch_ang_array[i] << " sin is: " << sin(launch_ang_array[i]);
 double launchx=800*diffx/xylen*sin(launch_ang_array[i]);
 cout << "\ndiffx is: " << diffx;
 double launchy=800*diffy/xylen*sin(launch_ang_array[i]);
 cout << "\n diff y is: " << diffy;
 double launchz=800*cos(launch_ang_array[i]);
 cout << "launchz is: " << launchz;
 vector_pair << "{{" << posx << "," << posy << "," << posz << "},{" << launchx << "," << launchy << "," << launchz << "}}";
 if (i<(9))
   {
     vector_pair << ",";
   }
 double  recz=-600*cos(rec_ang_array[i]);
 double recx=-600*diffx/xylen *sin(rec_ang_array[i]);
 double recy=-600*diffy/xylen *sin(rec_ang_array[i]);
 rec_pair_out << "{{" << antx << "," << anty << "," << antz -6357500<< "},{" << recx << "," << recy << "," << recz << "}}";
 cout << "\nRec ang is: " << rec_ang_array[i] << "x, y, z are: " << recx << " " << recy<< " " << recz; 
 if (i<8)
   {
     rec_pair_out << ",";
   }
 //2d stuff
 double x2d=(double)sqrt((detector->stations[0].strings[0].antennas[0].GetX()-event->Nu_Interaction[0].posnu.GetX())*(detector->stations[0].strings[0].antennas[0].GetX()-event->Nu_Interaction[0].posnu.GetX())+(detector->stations[0].strings[0].antennas[0].GetY()-event->Nu_Interaction[0].posnu.GetY())*(detector->stations[0].strings[0].antennas[0].GetY()-event->Nu_Interaction[0].posnu.GetY()));
 double z2d= detector->stations[0].strings[0].antennas[0].GetZ()-event->Nu_Interaction[0].posnu.GetZ();
 double launchx2d=800*sin(launch_ang_array[i]); 
 double launchz2d=800*cos(launch_ang_array[i]);
 //double antx2d=sqrt(antx*antx+anty*anty);
 // double antz2d=antz-6357500;
 double recx2d=600*sin(rec_ang_array[i]);
 double recz2d=600*cos(rec_ang_array[i]);
 launch2d << "Show[Plot[-2000, {x,-3000,0}, Axes->False,PlotStyle->Brown],Plot[0, {x,-3000,0},Axes->False],ListVectorFieldPlot[{{{" <<-x2d << "," << -z2d << "},{" << launchx2d << "," << launchz2d << "}},{{0,0},{" << -recx2d << "," << -recz2d << "}}},Axes->False]],";
 
    }//if i<10 and report->stations[0].strings[0].antennas[0].ray_sol_cnt !=
  }
 
//for loop ends
//output antenna position for mathematica

double antx=detector->stations[0].strings[0].antennas[0].GetX();
 cout << "\nantx is: " <<detector->stations[0].strings[0].antennas[0].GetX();
 double anty=detector->stations[0].strings[0].antennas[0].GetY();
 cout << "\nanty is: " <<detector->stations[0].strings[0].antennas[0].GetY();
 double antz=detector->stations[0].strings[0].antennas[0].GetZ();
 cout << "\nantz is: " <<detector->stations[0].strings[0].antennas[0].GetZ();
vector_pair << "\n\n\n {" << antx << "," << anty << "," << antz-6357500 << "}";








///////////////////////////////////////////////////






for(int i=0;i<settings->NNU; i++)
{
AraTree2[1]->GetEvent(i);
 if (report->stations[0].strings[0].antennas[0].ray_sol_cnt !=0)
   {
hpol_18->Fill(cos(report->stations[0].strings[0].antennas[0].Pol_vector[0].Theta()), event->Nu_Interaction[0].weight);
if (report->stations[0].Global_Pass!=0)
 {
   cout << "Global pass 18";
hpol_18_trig->Fill(cos(report->stations[0].strings[0].antennas[0].Pol_vector[0].Theta()), event->Nu_Interaction[0].weight);
 }//if pass
   }//if ray_sol_cnt>0
 }//for loop

 launch2d << "}]";
 output.close();
vector_pair.close();
 rec_pair_out.close();
 launch2d.close();
/* 
//commented out for performance reasons
 //Here we average  every event's value at a given frequency, and throw that frequency's average into an array
 double vfa,vemfa,dtea,dtha;
 double vfa_array[event->Nu_Interaction[0].vmmhz1m.size()];
 double vemfa_array[event->Nu_Interaction[0].vmmhz1m_em.size()];
 double dtea_array[event->Nu_Interaction[0].d_theta_em.size()];
 double dtha_array[event->Nu_Interaction[0].d_theta_had.size()];
double generic[event->Nu_Interaction[0].vmmhz1m.size()];
 for (int j=0;j<FREQ_SIZE;j++)
   {
     cout << "\nFrequency: " << j;
     vfa=0;
     vemfa=0;
     dtea=0;
     dtha=0;
     for (int i=0;i<settings->NNU;i++)
       {
	 AraTree2[0]->GetEvent(i);
	 vfa+=event->Nu_Interaction[0].vmmhz1m[j];
	 vemfa+=event->Nu_Interaction[0].vmmhz1m[j];
	 dtea+=event->Nu_Interaction[0].d_theta_em[j];
	 dtha+=event->Nu_Interaction[0].d_theta_had[j];
       }
     vfa=vfa/event->Nu_Interaction[0].vmmhz1m.size();
     vemfa=vemfa/event->Nu_Interaction[0].vmmhz1m_em.size();
     dtea=dtea/event->Nu_Interaction[0].d_theta_em.size();
     dtha=dtha/event->Nu_Interaction[0].d_theta_had.size();
     //     cout << "\n" << vfa << "\n" << event->Nu_Interaction[0].vmmhz1m.size();
     vfa_array[j]=vfa;
     vemfa_array[j]=vemfa;
     dtea_array[j]=dtea;
     dtha_array[j]=dtha;
     generic[j]=j;
     }
*/		  
 //////////////////////////////////////
 hflavor->SetMinimum(0.);
 hflavor->SetXTitle("Integer flavor");
hflavor->SetYTitle("Number of Events");
cflavor->cd();
hflavor->Draw();
cflavor->Print("flavor1.pdf");


/////////////////////////////////////



 hlogpnu->SetMinimum(0.);
 hlogpnu->SetXTitle("Log10 of PNU");
 hlogpnu->SetYTitle("Number of Events");
clogpnu->cd();
hlogpnu->Draw();
clogpnu->Print("pnulog.pdf");

//////////////////////////////////////////////////

 hn_interactions->SetMinimum(0.);
 hn_interactions->SetXTitle("Number of Interactions");
 hn_interactions->SetYTitle("Number of Events");

 cn_interactions->cd(1);
 hlogpnu->Draw();
 cn_interactions->cd(2);
 hn_interactions->Draw();
 cn_interactions->cd(3);
 hflavor->Draw();
 //Here we print the last three plots
 cn_interactions->Print("flavor-pnu-ninteractions.pdf");

//////////////////////////////////////////////////


 htheta->SetMinimum(0.);
 htheta->SetXTitle("Cos(Theta)");
 htheta->SetYTitle("Number of Events");




//////////////////////////////////////////////

 hphi->SetMinimum(0.);
 hphi->SetXTitle("Phi-Component of Vector");
 hphi->SetYTitle("Number of Events");
 ctheta->cd(1);
 htheta->Draw();
 ctheta->cd(2);
 hphi->Draw();
 ctheta->Print("entry_angle.pdf");
///////////////////////////////////////////

 
   //Draw viewing angle histogram
 hview_ang->SetMinimum(0.);
 hview_ang->SetXTitle("Viewing Angle");
 hview_ang->SetYTitle("Number of Events");
 cview_ang->cd(1);
 hview_ang->Draw();

 //////////////////////////////////////////
 //graph of phi vs. viewing angle
 TGraph *Gphi_view=new TGraph(settings->NNU, phi_array, view_ang_array);
 Gphi_view->SetTitle("");
 Gphi_view->GetXaxis()->SetTitle("Phi");
 Gphi_view->GetYaxis()->SetTitle("Viewing Angle");
cview_ang->cd(2);
 Gphi_view->Draw("a*");


//////////////////////////////////////////
//graph of theta vs. viewing angle
TGraph *Gtheta_view=new TGraph(settings->NNU, theta_array, view_ang_array);
 Gtheta_view->SetTitle("");
 Gtheta_view->GetXaxis()->SetTitle("Theta");
 Gtheta_view->GetYaxis()->SetTitle("Viewing Angle");
 cview_ang->cd(3);
Gtheta_view->Draw("a*");

///////////////////////////////////////////
//graph of cos theta vs. viewing angle
 TGraph *Gcostheta_view=new TGraph(settings->NNU, costheta_array, view_ang_array);
 Gcostheta_view->SetTitle("");
 Gcostheta_view->GetXaxis()->SetTitle("Cosine Theta");
 Gcostheta_view->GetYaxis()->SetTitle("Viewing Angle");
 cview_ang->cd(4);
 Gcostheta_view->Draw("a*");
 cview_ang->Print("view_angle.pdf");
////////////////////////////////////////////////
 TGraph *Gview_ang_2plus=new TGraph(settings->NNU, event_depth, view_ang_array2plus);
 Gview_ang_2plus->SetTitle("");
 Gview_ang_2plus->GetXaxis()->SetTitle("Event Depth");
 Gview_ang_2plus->GetYaxis()->SetTitle("Viewing angle (2+)");
 cview_ang_2plus->cd();
 Gview_ang_2plus->Draw("a*");
 cview_ang_2plus->Print("view_angle_2plus.pdf");
/////////////////////////////////////////////////
//Draw histogram for viewing angle including when there's no solution
 hview_ang_include->SetMinimum(0.);
 hview_ang_include->SetXTitle("Viewing Angle, including no solution");
 hview_ang_include->SetYTitle("Number of Events");
 cview_ang_include->cd(1);
 hview_ang_include->Draw();
 cview_ang_include->Print("view_angle_include.pdf");


////////////////////////////////////////////////

///////////////////////////////////////////////
 

  
  //Viewing angle historgram, including when there's no solution
 //Receiving angle histogram (not sorted by signal strength
 hrec_ang->SetMinimum(0.);
 hrec_ang->SetXTitle("Receiving Angle");
 hrec_ang->SetYTitle("Number of Events");
 crec_ang->cd(1);
 hrec_ang->Draw();

 //Viewing angle (sorted by signal strength) histograms
 cview_ang_multi->cd();
 hview_ang2->SetMinimum(0.);
 hview_ang2->SetXTitle("Viewing Angle");
 hview_ang2->SetYTitle("Number of Events");
 hview_ang2->SetLineColor(6);
 hview_ang2->Draw();
 

 hview_ang1->SetMinimum(0.);
 hview_ang1->SetXTitle("Viewing Angle");
 hview_ang1->SetYTitle("Number of Events");
 hview_ang1->SetLineColor(2);
 hview_ang1->Draw("histsame");

 hview_ang3->SetMinimum(0.);
 hview_ang3->SetXTitle("Viewing Angle");
 hview_ang3->SetYTitle("Number of Events");
 hview_ang3->SetLineColor(8);
 hview_ang3->Draw("histsame");

 hview_ang4->SetMinimum(0.);
 hview_ang4->SetXTitle("Viewing Angle");
 hview_ang4->SetYTitle("Number of Events");
 hview_ang4->SetLineColor(4);
 hview_ang4->Draw("histsame");
 cview_ang_multi->Print("view_ang_multi.pdf");
 //////////////////////////////////////////
 //receiving angle, sorted by signal strength
 crec_ang_multi->cd();
 hrec_ang1->SetMinimum(0);
 hrec_ang1->SetXTitle("Receiving Angle");
 hrec_ang1->SetYTitle("Number of Events");
 hrec_ang1->SetLineColor(2);
 hrec_ang1->Draw();

 hrec_ang2->SetLineColor(4);
  hrec_ang2->Draw("histsame");

 hrec_ang3->SetLineColor(6);
 hrec_ang3->Draw("histsame");

 hrec_ang4->SetLineColor(8);
 hrec_ang4->Draw("histsame");

 crec_ang_multi->Print("rec_ang_multi.pdf");
 ///////////////////////////////////////////////
 //zoomed in rec angle drawings
 hrec_ang_focus->SetMinimum(0.);
 hrec_ang_focus->SetXTitle("Receiving Angle, focused");
 hrec_ang_focus->SetYTitle("Number of Events");
 crec_ang_focus->cd(1);
 hrec_ang_focus->Draw();

crec_ang_focus->cd(2);
 hrec_ang_focus1->SetMinimum(0);
 hrec_ang_focus1->SetXTitle("Receiving Angle, focused");
 hrec_ang_focus1->SetYTitle("Number of Events");
 hrec_ang_focus1->SetLineColor(2);
 hrec_ang_focus1->Draw();

 hrec_ang_focus2->SetLineColor(4);
  hrec_ang_focus2->Draw("histsame");

 hrec_ang_focus3->SetLineColor(6);
 hrec_ang_focus3->Draw("histsame");

 hrec_ang_focus4->SetLineColor(8);
 hrec_ang_focus4->Draw("histsame");

 crec_ang_focus->Print("rec_ang_focus.pdf");
 ///////////////////////////////////////////////
 //Draw super focused receiving angles
 crec_ang_super_focus->cd(1);
 hrec_ang_super_focus_V->SetMinimum(0);
 hrec_ang_super_focus_V->SetXTitle("Receiving Angle for first numbered antenna, focused");
 hrec_ang_super_focus_V->SetYTitle("Number of events");
 hrec_ang_super_focus_V->SetLineColor(2);
 hrec_ang_super_focus_V->Draw();

 hrec_ang_super_focus_H->SetLineColor(4);
 hrec_ang_super_focus_H->Draw("histsame");

crec_ang_super_focus->cd(2);
 hrec_ang_super_focus_1st_V->SetMinimum(0);
 hrec_ang_super_focus_1st_V->SetXTitle("Receiving Angle for first numbered antenna, focused");
 hrec_ang_super_focus_1st_V->SetYTitle("Number of events");
 hrec_ang_super_focus_1st_V->SetLineColor(2);
 hrec_ang_super_focus_1st_V->Draw();

 hrec_ang_super_focus_1st_H->SetLineColor(4);
 hrec_ang_super_focus_1st_H->Draw("histsame");

 crec_ang_super_focus->Print("rec_ang_super_focus.pdf");
 ///////////////////////////////////////////////
 //reflect_angle, sorted by signal strength
 creflect_ang_multi->cd();
 hreflect_ang1->SetMinimum(0);
 hreflect_ang1->SetXTitle("Reflection Angle");
 hreflect_ang1->SetYTitle("Number of Events");
 hreflect_ang1->SetLineColor(2);
 hreflect_ang1->Draw();

 hreflect_ang2->SetLineColor(4);
  hreflect_ang2->Draw("histsame");

 hreflect_ang3->SetLineColor(6);
 hreflect_ang3->Draw("histsame");

 hreflect_ang4->SetLineColor(8);
 hreflect_ang4->Draw("histsame");

 creflect_ang_multi->Print("reflect_ang_multi.pdf");
 //Create distance histograms

 cdist->cd(1);
 hdist->SetMinimum(0);
 hdist->SetXTitle("Distance between posnu and antenna");
 hdist->SetYTitle("Number of Events");
 hdist->Draw();
 

 cdist->cd(2);
 hdistweighted1->SetMinimum(0);
 hdistweighted1->SetXTitle("Distance between posnu and antenna, r weighted");
 hdistweighted1->SetYTitle("Number of Events");
 hdistweighted1->Draw();

 cdist->cd(3);
 hdistweighted2->SetMinimum(0);
 hdistweighted2->SetXTitle("Distance between posnu and antenna, rsquared weighted");
 hdistweighted2->SetYTitle("Number of Events");
 hdistweighted2->Draw();

 cdist->cd(4);
 hdist1->SetMinimum(0);
 hdist1->SetXTitle("Distance between posnu and antenna");
 hdist1->SetYTitle("Number of Events");
 hdist1->SetLineColor(2);
 hdist1->Draw();

 hdist2->SetLineColor(4);
  hdist2->Draw("histsame");

 hdist3->SetLineColor(6);
 hdist3->Draw("histsame");

 hdist4->SetLineColor(8);
 hdist4->Draw("histsame");

 cdist->Print("dist.pdf");
 //create histogram graphs taking fern into account
 cdist_fern->cd(1);
 hdist->SetMinimum(0);
 hdist->SetXTitle("Distance between posnu and antenna");
 hdist->SetYTitle("Number of Events");
 hdist->Draw();
 
 cdist_fern->cd(1);
 hdist_fern->SetMinimum(0);
 hdist_fern->SetXTitle("Distance between posnu and antenna, contained in fern");
 hdist_fern->SetYTitle("Number of Events");
 hdist_fern->Draw();


 cdist_fern->cd(2);
 hdist_not_fern->SetMinimum(0);
 hdist_not_fern->SetXTitle("Distance between posnu and antenna, not in fern");
 hdist_not_fern->SetYTitle("Number of Events");
 hdist_not_fern->Draw();

 /////////////////////////////////////////////

 TGraph *gdist_rec_fern=new TGraph(settings->NNU, dist_fern_array, rec_ang_fern);
 gdist_rec_fern->SetTitle("");
 gdist_rec_fern->GetXaxis()->SetTitle("Distance - fern");
 gdist_rec_fern->GetYaxis()->SetTitle("Receiving Angle - fern");
 cdist_fern->cd(3);
 gdist_rec_fern->Draw("a*");

cdist_fern->Print("dist_fern.pdf"); 
//Distance histograms, 1st solution vs. Second Solution
  cdist_sol->cd(1);
 hdist->SetMinimum(0);
 hdist->SetXTitle("Distance between posnu and antenna, solution 1");
 hdist->SetYTitle("Number of Events");
 hdist->Draw();

  cdist_sol->cd(2);
 hdist_second_solution->SetMinimum(0);
 hdist_second_solution->SetXTitle("Distance between posnu and antenna, solution 2");
 hdist_second_solution->SetYTitle("Number of Events");
 hdist_second_solution->Draw();
 


 //Drawing launch angle histograms
 claunch->cd(1);
 hlaunch->SetMinimum(0);
 hlaunch->SetXTitle("Launch Angle");
 hlaunch->SetYTitle("Number of Events");
 hlaunch->Draw();

 claunch->cd(2);
 hlaunch1->SetMinimum(0);
 hlaunch1->SetXTitle("Launch angle");
 hlaunch1->SetYTitle("Number of Events");
 hlaunch1->SetLineColor(2);
 hlaunch1->Draw();

 hlaunch2->SetLineColor(4);
 hlaunch2->Draw("histsame");

 hlaunch3->SetLineColor(6);
 hlaunch3->Draw("histsame");

 hlaunch4->SetLineColor(8);
 hlaunch4->Draw("histsame");
 claunch->Print("launch.pdf");

 //////////////////////////////////////////

 //Now we're going to make a receiving angle graph for each antenna
 /*
//This is commented out for performance reasons
   TCanvas *crec_ang_focus_all = new TCanvas("crec_ang_focus_all","crec_ang_focus_all",0,0,1000,700);//Canvas for graph of receiving angle with no sorting 
 crec_ang_focus_all->Divide(4,4); //we need a big canvas for it
	 TH1F *hrec_ang_focus_all[16];
      for (int k=0; k< report->stations[0].strings.size();k++)
	 {
     for (int j=0; j<report->stations[0].strings[k].antennas.size(); j++)
       {
	 cout<<"Doing receiving angle for string  "<< k << " antenna " << j<< "\n";
int num=k*4+j;
//Form("somestringstuff%u-for v", num
 hrec_ang_focus_all[num]=new TH1F(Form("rec_ang_focus_all%u",num),"", 25,1.5,1.65);

 for(int i=0; i<settings->NNU; i++)
   {
     AraTree2[0]->GetEvent(i);

   if(report->stations[0].strings[k].antennas[j].ray_sol_cnt !=0)//check if solution exists
     {
      hrec_ang_focus_all[num]->Fill(report->stations[0].strings[k].antennas[j].rec_ang[0]);
     } //if
   }//for i to NNU
 hrec_ang_focus_all[num]->SetMinimum(0.);
 hrec_ang_focus_all[num]->SetXTitle("Receiving Angle");
 hrec_ang_focus_all[num]->SetYTitle("Number of Events");
 crec_ang_focus_all->cd(num+1);
 hrec_ang_focus_all[num]->Draw();
  }//for j to number of antenna
	 }//for k to number of strings
	 crec_ang_focus_all->Print("rec_ang_all.pdf");
 */ 
 
     ////////////////////////////////////
     //Draw graphs
 TGraph *Glaunch_view=new TGraph(settings->NNU, launch_ang_narrow, view_launch);
 Glaunch_view->SetTitle("");
 Glaunch_view->GetXaxis()->SetTitle("Launch Angle");
 Glaunch_view->GetYaxis()->SetTitle("Viewing Angle");
 claunch_compare->cd(1);
 Glaunch_view->Draw("a*");

 ////////////////////////////////////////
 TGraph *Glaunch_rec=new TGraph(settings->NNU, launch_ang_narrow, rec_launch);
 Glaunch_rec->SetTitle("");

 Glaunch_rec->GetXaxis()->SetTitle("Launch Angle");
 Glaunch_rec->GetYaxis()->SetTitle("Receiving Angle");
 claunch_compare->cd(2);
 Glaunch_rec->Draw("a*");
 ////////////////////////////////////////
 TGraph *Glaunch_dist=new TGraph(settings->NNU, launch_ang_narrow, dist_launch);
 Glaunch_dist->SetTitle("");
 Glaunch_dist->GetXaxis()->SetTitle("Launch Angle");
 Glaunch_dist->GetYaxis()->SetTitle("Distance from Antenna to Posnu");
 claunch_compare->cd(3);
 Glaunch_dist->Draw("a*");
 ///////////////////////////////////////
 /*TGraph *Glaunch_reflect=new TGraph(settings->NNU, launch_ang_array, reflect_launch);
 Glaunch_reflect->SetTitle("");
 Glaunch_reflect->GetXaxis()->SetTitle("Launch Angle");
 Glaunch_reflect->GetYaxis()->SetTitle("Reflection Angle");
 claunch_compare->cd(4);
 Glaunch_reflect->Draw("a*");*/
 ///////////////////////////////////////
 TGraph *Glaunch_reflect_not_hundred=new TGraph(settings->NNU, launch_ang_narrow, reflect_launch_not_hundred);
 Glaunch_reflect_not_hundred->SetTitle("");
 Glaunch_reflect_not_hundred->GetXaxis()->SetTitle("Launch Angle");
 Glaunch_reflect_not_hundred->GetYaxis()->SetTitle("Reflection Angle that Exists");
 claunch_compare->cd(4);
 Glaunch_reflect_not_hundred->Draw("a*");
 //print the canvas
 claunch_compare->Print("launch_comparison.pdf");

 ////////////////////////////////////////
 //draw graphs of launch vs rec, fern vs. not
 TGraph *Glaunch_rec_fern=new TGraph(settings->NNU, launch_ang_fern, rec_ang_fern);
 Glaunch_rec_fern->SetTitle("");
 Glaunch_rec_fern->GetXaxis()->SetTitle("Launch Angle - fern");
 Glaunch_rec_fern->GetYaxis()->SetTitle("Receiving Angle - fern");
 crec_launch_fern->cd(1);
 Glaunch_rec_fern->Draw("a*");
 //////////////////////////////////////////
 TGraph *Glaunch_rec_not_fern=new TGraph(settings->NNU, launch_ang_not_fern, rec_ang_not_fern);
 Glaunch_rec_not_fern->SetTitle("");
 Glaunch_rec_not_fern->GetXaxis()->SetTitle("Launch Angle - not fern");
 Glaunch_rec_not_fern->GetYaxis()->SetTitle("Receiving Angle - not fern");
 crec_launch_fern->cd(2);
 Glaunch_rec_not_fern->Draw("a*");
 ///////////////////////////////////////////
 TGraph *Glaunch_rec_zoom_fern=new TGraph(settings->NNU, launch_ang_zoom_fern, rec_ang_zoom_fern);
 Glaunch_rec_zoom_fern->SetTitle("");
 Glaunch_rec_zoom_fern->GetXaxis()->SetTitle("Zoom Launch Angle - fern");
 Glaunch_rec_zoom_fern->GetYaxis()->SetTitle("Zoom Receiving Angle - fern");
 crec_launch_fern->cd(3);
 Glaunch_rec_zoom_fern->Draw("a*");
 //////////////////////////////////////////
 TGraph *Glaunch_rec_zoom_not_fern=new TGraph(settings->NNU, launch_ang_zoom_not_fern, rec_ang_zoom_not_fern);
 Glaunch_rec_zoom_not_fern->SetTitle("");
 Glaunch_rec_zoom_not_fern->GetXaxis()->SetTitle("Zoomed Launch Angle - not fern");
 Glaunch_rec_zoom_not_fern->GetYaxis()->SetTitle("Zoomed Receiving Angle - not fern");
 crec_launch_fern->cd(4);
 Glaunch_rec_zoom_not_fern->Draw("a*");
 crec_launch_fern->Print("rec_launch_fern.pdf");
 /////////////////////////////////////////////
 TGraph *glaunch_dist_wide = new TGraph(settings->NNU, launch_ang_array, dist_array);
 glaunch_dist_wide->SetTitle("");
 glaunch_dist_wide->GetXaxis()->SetTitle("Launch Angle");
 glaunch_dist_wide->GetYaxis()->SetTitle("Distance from Antenna to Posnu");
 cdist_sol->cd(3);
 glaunch_dist_wide->Draw("a*");
cdist_sol->Print("dist_sol.pdf");
/////////////////////////////////////////
//draw phi vs. receiving angle graph
 TGraph *Gphi_rec=new TGraph(settings->NNU, phi_array, rec_ang_array);
 Gphi_rec->SetTitle("");
 Gphi_rec->GetXaxis()->SetTitle("Phi");
 Gphi_rec->GetYaxis()->SetTitle("Receiving Angle");
 crec_ang->cd(2);
 Gphi_rec->Draw("a*");

//////////////////////////////////////////
//draw theta vs. receiving angle graph
TGraph *Gtheta_rec=new TGraph(settings->NNU, theta_array, rec_ang_array);
 Gtheta_rec->SetTitle("");
 Gtheta_rec->GetXaxis()->SetTitle("Theta");
 Gtheta_rec->GetYaxis()->SetTitle("Receiving Angle");
 crec_ang->cd(3);
Gtheta_rec->Draw("a*");

///////////////////////////////////////////
//Draw cos theta vs receiving angle graph
 TGraph *Gcostheta_rec=new TGraph(settings->NNU, costheta_array, rec_ang_array);
 Gcostheta_rec->SetTitle("");
 Gcostheta_rec->GetXaxis()->SetTitle("Cosine Theta");
 Gcostheta_rec->GetYaxis()->SetTitle("Receiving Angle");
 crec_ang->cd(4);
Gcostheta_rec->Draw("a*");
crec_ang->Print("rec_angle.pdf");
////////////////////////////////////////////
//Here we make a graph of the receiving angle for only the second antenna (which we think is an H antenna)

hrec_ang_second_antenna->SetMinimum(0.);
 hrec_ang_focus_second_antenna->SetMinimum(0.);
 hrec_ang_second_antenna->SetXTitle("Receiving Angle, 2nd antenna");
hrec_ang_focus_second_antenna->SetXTitle("Receiving Angle, 2nd antenna zoom");
 hrec_ang_focus_second_antenna->SetYTitle("Number of Events");
 hrec_ang_focus_second_antenna->SetYTitle("Number of Events");
 crec_ang2->cd(1);
 hrec_ang_second_antenna->Draw(); 
crec_ang2->cd(2);
 hrec_ang_focus_second_antenna->Draw();
 crec_ang2->Print("rec_ang_secondantenna.pdf");
////////////////////////////////////////////

//Draw reflection histogram
 hreflect_ang->SetMinimum(0.);
 hreflect_ang->SetXTitle("Reflect Angle");
 hreflect_ang->SetYTitle("Number of Events");
 creflect_ang->cd(1);
 hreflect_ang->Draw();
 
 //////////////////////////////////////////
 cout << "\nThe number of events where reflecting angle doesn't equal one hundred  is " << count_reflect << ".\n";

 //////////////////////////////////////////
 //graph phi vs reflection angle
 TGraph *Gphi_reflect=new TGraph(1000, phi_array, reflect_ang_array);
 Gphi_reflect->SetTitle("");
 Gphi_reflect->GetXaxis()->SetTitle("Phi");
 Gphi_reflect->GetYaxis()->SetTitle("Reflect Angle");
creflect_ang->cd(2);
 Gphi_reflect->Draw("a*");

//////////////////////////////////////////
//graph theta vs. reflection angle
TGraph *Gtheta_reflect=new TGraph(1000, theta_array, reflect_ang_array);
 Gtheta_reflect->SetTitle("");
 Gtheta_reflect->GetXaxis()->SetTitle("Theta");
 Gtheta_reflect->GetYaxis()->SetTitle("Reflect Angle");
 creflect_ang->cd(3);
Gtheta_reflect->Draw("a*");

///////////////////////////////////////////
//graph cos theta vs. reflection angle
 TGraph *Gcostheta_reflect=new TGraph(1000, costheta_array, reflect_ang_array);
 Gcostheta_reflect->SetTitle("");
 Gcostheta_reflect->GetXaxis()->SetTitle("Cosine Theta");
 Gcostheta_reflect->GetYaxis()->SetTitle("Reflect Angle");
 creflect_ang->cd(4);
Gcostheta_reflect->Draw("a*");
creflect_ang->Print("reflect_angle.pdf");
////////////////////////////////////////////////
//histogram for receiving angle when event is reflected
 hrec_ang_not_hundred->SetMinimum(0.);
 hrec_ang_not_hundred->SetXTitle("Receiving Angle when Reflect Angle isn't 100");
 hrec_ang_not_hundred->SetYTitle("Number of Events");
 crec_ang_not_hundred->cd();
 hrec_ang_not_hundred->Draw();
 crec_ang_not_hundred->Print("rec_ang_not_hundred.pdf");
////////////////////////////////////////////////
//Draw Interactions histograms
 hnu_nubar->SetMinimum(0.);
 hnu_nubar->SetXTitle("nu_nubar");
 hnu_nubar->SetYTitle("Number of Events");
 cinteractions_1->cd(1);
 hnu_nubar->Draw();


 /////////////////////////////////////////////////
 hemfrac->SetMinimum(0.);
 hemfrac->SetXTitle("emfrac");
 hemfrac->SetYTitle("Number of Events");
 cinteractions_1->cd(2);
 hemfrac->Draw();
 
 /////////////////////////////////////////////////
 hhadfrac->SetMinimum(0.);
 hhadfrac->SetXTitle("hadfrac, log10");
 hhadfrac->SetYTitle("Number of Events");
 cinteractions_1->cd(4);
 hhadfrac->Draw();
 ////////////////////////////////////////////////
 hrst->SetMinimum(0.);
 hrst->SetXTitle("ray_sol_count, not zero");
 hrst->SetYTitle("Number of Events");
 cinteractions_1->cd(3);
 hrst->Draw();
cinteractions_1->Print("interactions_1.pdf");
 /////////////////////////////////////////////////
 hvmmhz1m_tmp->SetMinimum(0.);
 hvmmhz1m_tmp->SetXTitle("vmhz1m_tmp, log10");
 hvmmhz1m_tmp->SetYTitle("Number of Events");
 cinteractions_2->cd(1);
 hvmmhz1m_tmp->Draw();
 /////////////////////////////////////////////////
 hvmmhz1m->SetMinimum(0.);
 hvmmhz1m->SetXTitle("vmhz1m average");
 hvmmhz1m->SetYTitle("Number of Events");
 cinteractions_2->cd(2);
 hvmmhz1m->Draw();
 
 /////////////////////////////////////////////////
 hvmmhz1m_log->SetXTitle("vmhz1m average, log10");
 hvmmhz1m_log->SetYTitle("Number of Events");
 cinteractions_2->cd(3);
 hvmmhz1m_log->Draw();
 /////////////////////////////////////////////////
 hd_theta_em->SetMinimum(0.);
 hd_theta_em->SetXTitle("d_theta_em average");
 hd_theta_em->SetYTitle("Number of Events");
 cinteractions_2->cd(4);
 hd_theta_em->Draw();
cinteractions_2->Print("interactions_2.pdf");
 /////////////////////////////////////////////////
 hd_theta_had->SetMinimum(0.);
 hd_theta_had->SetXTitle("d_theta_had average");
 hd_theta_had->SetYTitle("Number of Events");
 cinteractions_3->cd(1);
 hd_theta_had->Draw();
 /////////////////////////////////////////////////
 hemfrac_complement->SetMinimum(0.);
 hemfrac_complement->SetXTitle("log(1-emfrac)");
 hemfrac_complement->SetYTitle("Number of Events");
 cinteractions_3->cd(2);
 hemfrac_complement->Draw();
 //////////////////////////////////////////
 hhad_em_sum->SetMinimum(0.);
 hhad_em_sum->SetXTitle("hadron fraction + em fraction");
 hhad_em_sum->SetYTitle("Number of events");
 cinteractions_3->cd(3);
 hhad_em_sum->Draw(); 
 cinteractions_3->Print("interactions_3.pdf");
 ///////////////////////////////
		     //Commented out for performance
 /*TGraph *Gvfa=new TGraph(FREQ_SIZE, generic, vfa_array);
 Gvfa->SetTitle("");
 Gvfa->GetXaxis()->SetTitle("frequency");
 Gvfa->GetYaxis()->SetTitle("Average vmmhz1m");
 cinteractions_4->cd(1);
 Gvfa->Draw("a*");
 ///////////////////////////////////////////////
 TGraph *Gvemfa=new TGraph(FREQ_SIZE, generic, vemfa_array);
 Gvemfa->SetTitle("");
 Gvemfa->GetXaxis()->SetTitle("frequency");
 Gvemfa->GetYaxis()->SetTitle("Average vmmhz1m_em");
 cinteractions_4->cd(2);
 Gvemfa->Draw("a*");
 //////////////////////////////////////////////////
 TGraph *Gdtea=new TGraph(FREQ_SIZE, generic, dtea_array);
 Gdtea->SetTitle("");
 Gdtea->GetXaxis()->SetTitle("frequency");
 Gdtea->GetYaxis()->SetTitle("Average cone thickness em");
 cinteractions_4->cd(3);
 Gdtea->Draw("a*");
 //////////////////////////////////////////////////
 TGraph *Gdtha=new TGraph(FREQ_SIZE, generic, dtha_array);
 Gdtha->SetTitle("");
 Gdtha->GetXaxis()->SetTitle("frequency");
 Gdtha->GetYaxis()->SetTitle("Average cone thickness had");
 cinteractions_4->cd(4);
 Gdtha->Draw("a*");
 cinteractions_4->Print("interactions_4.pdf");
 */ 
/////////////////////////////////////
 hpol_17->SetMinimum(0.);
 hpol_17->SetXTitle("Cos(polarization-angle) 17&18");
 hpol_17->SetYTitle("Number of Events");
 cpol->cd(1);
 hpol_17->Draw();
	/////////////////////////////////////////////

hpol_18->SetLineColor(4);
 hpol_18->Draw("histsame");


 ////////////////////////////////////////////////

hpol_17_trig->SetMinimum(0.);
 hpol_17_trig->SetXTitle("Cos(polarization-angle) 17&18");
 hpol_17_trig->SetYTitle("Number of Events passing trigger");
 cpol->cd(2);
 hpol_17_trig->Draw();
	/////////////////////////////////////////////

hpol_18_trig->SetLineColor(4);
 hpol_18_trig->Draw("histsame");
 
 ////////////////////////////////////////////

 hpol_factor->SetMinimum(0.);
 hpol_factor->SetXTitle("Polarization factor, all events");
 hpol_factor->SetYTitle("Number of Events");
 cpol->cd(3);
 hpol_factor->Draw();
 /////////////////////////////////////////////
hpol_factor_trig->SetMinimum(0.);
 hpol_factor_trig->SetXTitle("Polarization Factor, triggered");
 hpol_factor_trig->SetYTitle("Number of Events");
 cpol->cd(4);
 hpol_factor_trig->Draw();
cpol->Print("pol.pdf");
 /////////////////////
////////////////////////////////////////////

 hpolz_all->SetMinimum(0.);
 hpolz_all->SetXTitle("Z-Polarization, all events");
 hpolz_all->SetYTitle("Number of Events");
 cpol2->cd(1);
 hpolz_all->Draw();
 ///////////////////////////////////////////////

 hpolz_trig->SetMinimum(0.);
 hpolz_trig->SetXTitle("Z-Polarization,triggered");
 hpolz_trig->SetYTitle("Number of Events");
 cpol2->cd(2);
 hpolz_trig->Draw();
 ///////////////////////////////////////////////

 hpol_flat_all->SetMinimum(0.);
 hpol_flat_all->SetXTitle("XY-Polarization, all events");
 hpol_flat_all->SetYTitle("Number of Events");
 cpol2->cd(3);
 hpol_flat_all->Draw();
 ///////////////////////////////////////////////

 hpol_flat_trig->SetMinimum(0.);
 hpol_flat_trig->SetXTitle("XY-Polarization,triggered");
 hpol_flat_trig->SetYTitle("Number of Events");
 cpol2->cd(4);
 hpol_flat_trig->Draw();
 cpol2->Print("pol-components.pdf");
 ///////////////////////////////////////////////
//draw phi vs. receiving angle graph
 TGraph *Gpol_time=new TGraph(settings->NNU,  cos_pol_ang, diff_time);
 Gpol_time->SetTitle("");
 Gpol_time->GetXaxis()->SetTitle("Cosine of polarization angle");
 Gpol_time->GetYaxis()->SetTitle("Difference in time between first and last station");
 cpol3->cd(1);
 Gpol_time->Draw("a*");
 ///////////////////////////////////////////////
//draw phi vs. receiving angle graph
 TGraph *Gpol_time2=new TGraph(settings->NNU,  cos_pol_ang, diff_time_sub_hundred);
 Gpol_time2->SetTitle("");
 Gpol_time2->GetXaxis()->SetTitle("Cosine of polarization angle");
 Gpol_time2->GetYaxis()->SetTitle("Sub hundred Difference in time between first and last station");
 cpol3->cd(2);
 Gpol_time2->Draw("a*");
 ///////////////////////////////////////////////
//draw phi vs. receiving angle graph
 TGraph *Gpol_time3=new TGraph(settings->NNU,  cos_pol_ang, diff_time_no_ref);
 Gpol_time3->SetTitle("");
 Gpol_time3->GetXaxis()->SetTitle("Cosine of polarization angle");
 Gpol_time3->GetYaxis()->SetTitle("No reflect Difference in time between first and last station");
 cpol3->cd(3);
 Gpol_time3->Draw("a*");
 cpol3->Print("pol3.pdf");
 /////////////////////////////////////////////////


 /////////////////////////////////////////////////
 cout << "The maximum difference in viewing angle is: " << maxdiff << "\n";
}


TStyle* RootStyle() {
  
  
    TStyle *RootStyle = new TStyle("Root-Style","The Perfect Style for Plots ;-)");
  
#ifdef __CINT__
    TStyle *GloStyle;
    GloStyle = gStyle;                          // save the global style reference
  
    gStyle = RootStyle;
#endif
    // otherwise you need to call TROOT::SetStyle("Root-Style")
  
    // Paper size
  
    RootStyle->SetPaperSize(TStyle::kUSLetter);
  
    // Canvas
  
    RootStyle->SetCanvasColor     (0);
    RootStyle->SetCanvasBorderSize(10);
    RootStyle->SetCanvasBorderMode(0);
    RootStyle->SetCanvasDefH      (600);
    RootStyle->SetCanvasDefW      (600);
    RootStyle->SetCanvasDefX      (10);
    RootStyle->SetCanvasDefY      (10);
  
    // Pads
  
    RootStyle->SetPadColor       (0);
    RootStyle->SetPadBorderSize  (10);
    RootStyle->SetPadBorderMode  (0);
    //  RootStyle->SetPadBottomMargin(0.13);
    RootStyle->SetPadBottomMargin(0.16);
    RootStyle->SetPadTopMargin   (0.08);
    RootStyle->SetPadLeftMargin  (0.18);
    RootStyle->SetPadRightMargin (0.05);
    RootStyle->SetPadGridX       (0);
    RootStyle->SetPadGridY       (0);
    RootStyle->SetPadTickX       (1);
    RootStyle->SetPadTickY       (1);
  
    // Frames
  
    RootStyle->SetFrameFillStyle ( 0);
    RootStyle->SetFrameFillColor ( 0);
    RootStyle->SetFrameLineColor ( 1);
    RootStyle->SetFrameLineStyle ( 0);
    RootStyle->SetFrameLineWidth ( 2);
    RootStyle->SetFrameBorderSize(10);
    RootStyle->SetFrameBorderMode( 0);
  
  
    // Histograms
  
    RootStyle->SetHistFillColor(0);
    RootStyle->SetHistFillStyle(1);
    RootStyle->SetHistLineColor(1);
    RootStyle->SetHistLineStyle(0);
    RootStyle->SetHistLineWidth(2);
  
    // Functions
  
    RootStyle->SetFuncColor(1);
    RootStyle->SetFuncStyle(0);
    RootStyle->SetFuncWidth(1);
  
    //Legends
  
    RootStyle->SetStatBorderSize(2);
    RootStyle->SetStatFont      (42);
    //  RootStyle->SetOptStat       (111111);
    RootStyle->SetOptStat       (0);
    RootStyle->SetStatColor     (0);
    RootStyle->SetStatX         (0.93);
    RootStyle->SetStatY         (0.90);
    RootStyle->SetStatFontSize  (0.07);
    //  RootStyle->SetStatW         (0.2);
    //  RootStyle->SetStatH         (0.15);
  
    // Labels, Ticks, and Titles
  
    RootStyle->SetTickLength ( 0.015,"X");
    RootStyle->SetTitleSize  ( 0.06,"X");
    RootStyle->SetTitleOffset( 1.30,"X");
    RootStyle->SetTitleBorderSize(0);
    //  RootStyle->SetTitleFontSize((double)3.);
    RootStyle->SetLabelOffset( 0.015,"X");
    RootStyle->SetLabelSize  ( 0.050,"X");
    RootStyle->SetLabelFont  ( 42   ,"X");
  
    RootStyle->SetTickLength ( 0.015,"Y");
    RootStyle->SetTitleSize  ( 0.06,"Y");
    RootStyle->SetTitleOffset( 1.100,"Y");
    RootStyle->SetLabelOffset( 0.015,"Y");
    RootStyle->SetLabelSize  ( 0.050,"Y");
    RootStyle->SetLabelFont  ( 42   ,"Y");
  
    RootStyle->SetTitleFont  (42,"XY");
    RootStyle->SetTitleColor  (1);
  
    // Options
  
    RootStyle->SetOptFit     (0);
  
    RootStyle->SetMarkerStyle(20);
    RootStyle->SetMarkerSize(0.4);
  
    //  cout << ">> Style initialized with the Root Style!" << endl;
    //  cout << ">> " << modified << endl << endl;
    return RootStyle;

}





