#include "Event.h"
#include <string>
#include "Settings.h"
#include "Spectra.h"
#include "Primaries.h"
#include "TMath.h"
#include "Report.h"

ClassImp(Event);

using namespace std;


Event::Event () {
    //default constructor
    Initialize ();
}

void Event::Initialize () {

    Nu_Interaction.clear();
    //test_report.clear();
    n_interactions = 1;  // total number of interaction (including secondaries)

}


//Event::Event (Settings *settings1, Spectra *spectra1, Primaries *primary1, IceModel *icemodel, Detector *detector, Signal *signal, Secondaries *sec1 ) {
Event::Event (Settings *settings1, Spectra *spectra1, Primaries *primary1, IceModel *icemodel, Detector *detector, Signal *signal, Secondaries *sec1, int event_num ) {

    Initialize ();

    inu_thrown = event_num;

    Choose_Evt_Type (settings1);

    if (Event_type == 0) { // if only neutrino events exist
       
        pnu = spectra1->GetNuEnergy();
        
        if (settings1->EVENT_GENERATION_MODE == 1){
            // pnu = pow(10., settings1->PNU[inu_thrown]);
            settings1->EXPONENT = settings1->PNU[inu_thrown];
            spectra1 = new Spectra(settings1->EXPONENT);
            settings1->SELECT_FLAVOR = settings1->NUFLAVORINT[inu_thrown];
            settings1->NU_NUBAR_SELECT_MODE = settings1->NUBAR[inu_thrown];
            settings1->SELECT_CURRENT = settings1->CURRENTINT[inu_thrown];
            settings1->INTERACTION_MODE = 2;
            settings1->POSNU_R = settings1->IND_POSNU_R[inu_thrown];
            settings1->POSNU_THETA = settings1->IND_POSNU_THETA[inu_thrown];
            settings1->POSNU_PHI = settings1->IND_POSNU_PHI[inu_thrown];
            settings1->NNU_THIS_THETA = 1;
            settings1->NNU_D_THETA = 0.0;
            settings1->NNU_THETA = settings1->IND_NNU_THETA[inu_thrown];
            settings1->NNU_THIS_PHI = 1;
            settings1->NNU_D_PHI = 0.0;
            settings1->NNU_PHI = settings1->IND_NNU_PHI[inu_thrown];
            settings1->YPARAM = 2;
            settings1->ELAST_Y = settings1->ELAST[inu_thrown];
        }
        pnu = spectra1->GetNuEnergy();
        // cout << pnu << endl;
        /*
        double hereTheta = 10.;
        Vector output;
        output.SetX(-1.*TMath::Sin(hereTheta*3.1415926535/180.));
        output.SetY(0.);
        output.SetZ(TMath::Cos(hereTheta*3.1415926535/180.));
        nnu = output;
        */
        //nuflavor = primary1->GetNuFlavor();
        nuflavor = primary1->GetNuFlavor(settings1);
        // nu_nubar = primary1->GetNuNuBar(nuflavor);
        nu_nubar = primary1->GetNuNuBar(nuflavor, settings1);


        /*
        if (settings1->NNU_THIS_THETA==1) {    // set specific theta angle for nnu
            nnu = primary1->GetThatDirection(settings1->NNU_THETA, settings1->NNU_D_THETA);
        }
        else { // nnu angle random
            nnu = primary1->GetAnyDirection();
        }
        */
        
        if (nuflavor=="nue"){
            nuflavorint=1;
        }
        else if (nuflavor=="numu"){
            nuflavorint=2;
        }
        else if (nuflavor=="nutau"){
            nuflavorint=3;
        }

        Interaction *Nu_temp;
        //Report *report_tmp;

        Nu_temp = new Interaction (pnu, nuflavor, nu_nubar, n_interactions, icemodel, detector, settings1, primary1, signal, sec1 );
        //report_tmp = new Report(detector ,settings1);
        
        Nu_Interaction.push_back(*Nu_temp);  // for the first interaction
        //test_report.push_back(*report_tmp);

        delete Nu_temp;

        // for multiple interactions...
        /*
        while (interaction_count < n_interactions) {    // not sure if this will work???
            Nu_tmp = new Interaction (...., n_interactions );
            Nu_Interaction.push_back( Nu_tmp );
        }
        */
        
    } 

    if (Event_type == 10) { // if only arbitrary events exist
        
  
        pnu = 0;
        // cout << pnu << endl;
        /*
        double hereTheta = 10.;
        Vector output;
        output.SetX(-1.*TMath::Sin(hereTheta*3.1415926535/180.));
        output.SetY(0.);
        output.SetZ(TMath::Cos(hereTheta*3.1415926535/180.));
        nnu = output;
        */
        // nuflavor = primary1->GetNuFlavor();
        nuflavor = "";
        nuflavorint = 0;
        nu_nubar = 0;


        /*
        if (settings1->NNU_THIS_THETA==1) {    // set specific theta angle for nnu
            nnu = primary1->GetThatDirection(settings1->NNU_THETA, settings1->NNU_D_THETA);
        }
        else { // nnu angle random
            nnu = primary1->GetAnyDirection();
        }
        */
        
        Interaction *Nu_temp;
        // Report *report_tmp;

        Nu_temp = new Interaction (settings1, detector, icemodel, primary1, signal );
        // report_tmp = new Report(detector ,settings1);
        
        Nu_Interaction.push_back(*Nu_temp);  // for the first interaction
        // test_report.push_back(*report_tmp);

        delete Nu_temp;

        // for multiple interactions...
        /*
        while (interaction_count < n_interactions) {    // not sure if this will work???
            Nu_tmp = new Interaction (...., n_interactions );
            Nu_Interaction.push_back( Nu_tmp );
        }
        */

    }

    IsCalpulser = primary1->IsCalpulser;

}

Event::~Event() {
    Nu_Interaction.clear();
}


void Event::delete_all() {
    Nu_Interaction.clear();
}


void Event::Choose_Evt_Type (Settings *settings1) {

    if (settings1->EVENT_TYPE==0) {
        //cout<<"Only Neutrino Events!"<<endl;
        Event_type = 0;
    }
    else if (settings1->EVENT_TYPE == 10){
        //            cout<<"Currently, only neutrino events possible!"<<endl;
        //            cout<<"Change Evt_type from "<<settings1->EVENT_TYPE<<" to 0"<<endl;
        Event_type = 10;
    }
    else {
        cout<<"Currently, only neutrino (EVET_TYPE=0) and arbitrary (EVENT_TYPE=10) events possible!"<<endl;
        cout<<"Change Evt_type from "<<settings1->EVENT_TYPE<<" to 0"<<endl;
        Event_type = 0;
    }
}


