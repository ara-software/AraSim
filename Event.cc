#include "Event.h"
#include <string>
#include "Settings.h"
#include "Spectra.h"
#include "Primaries.h"
#include "secondaries.hh"
#include "TMath.h"
#include "Report.h"
#include "Constants.h"

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

    // Set the default maximum number of interactions to n_interactions
    interaction_cnt_max = n_interactions;

    // If the event generation mode is set to 1 (event read-in mode), 
    // adjust interaction count to the number of interactions per neutrino primary
    // calculated in Settings.cc (ReadEvtFile) by looping over events with the same EVID
    if (settings1->EVENT_GENERATION_MODE==1){

        // Retrieve the number of interactions for the current neutrino primary
        interaction_cnt_max = settings1->INT_PER_NNU[event_num];

        // Initialize counter for the cumulative number of interactions
        // note that inu_thrown means "neutrinos thrown" typically but it also means
        // "interactions" now that lines in event lists can contain other interactions  
        // besides neutrino primaries
        inu_thrown = 0;

        // Loop over previous events to sum up the total number of interactions so far
        for (int i = 0; i < event_num; ++i) {
            inu_thrown += settings1->INT_PER_NNU[i];
        }  
    }

    for (int interaction_cnt = 1; interaction_cnt <= interaction_cnt_max; ++interaction_cnt){ 

        Choose_Evt_Type (settings1);

        if (Event_type == 0) { // if only neutrino events exist
       
            pnu = spectra1->GetNuEnergy();
        
            if (settings1->EVENT_GENERATION_MODE == 1){

                // Determine if PNU is provided as settings->EXPONENT 
                //   or actual energy and create spectra
                if (settings1->PNU[inu_thrown] < 1000) { // pnu is settings->EXPONENT
                    settings1->EXPONENT = settings1->PNU[inu_thrown];
                } 
                else if (settings1->PNU[inu_thrown] < 1e14) { // pnu is energy in GeV
                    pnu = settings1->PNU[inu_thrown] * 1e9; // convert to eV
                    settings1->EXPONENT = 400 + log10(pnu)*10.;
                }
                else { // pnu is energy in eV
                    pnu = settings1->PNU[inu_thrown];
                    settings1->EXPONENT = 400 + log10(pnu)*10.;
                }
                spectra1 = new Spectra(settings1);
                if (settings1->PNU[inu_thrown] < 1000) { 
                    // set the pnu that wasn't set manually bc 
                    // provided pnu was in exponent form
                    pnu = spectra1->GetNuEnergy();
                } 

                // Prepare the rest of the parameters
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

                //Add event ID
                event_ID.push_back(settings1->EVID[inu_thrown]);
            }

            //Calculate the interactions birth time for secondaries
            if (interaction_cnt == 1){
                // If this is the first interaction, clear the birth time storage
                interactions_birth_time.clear();

                // The first interaction is set to have a birth time of 0.0 (reference point)
                interactions_birth_time.push_back(0.0);

                // Store the index of the first neutrino primary for this event
                first_vertex_idx = inu_thrown;
            }
            else{
                // Compute the birth time for secondaries (inu_thrown) of a primary 
                // This is determined by dividing the interaction distance by the speed of light (CLIGHT),
                // giving the travel time from the previous interaction.
                // Assumes particles travel at ~CLIGHT.
                double tmp_birth_time = (sec1->get_interaction_distance(first_vertex_idx, inu_thrown, settings1)/CLIGHT);

                // Store the computed birth time for this secondary interaction
                interactions_birth_time.push_back(tmp_birth_time);
            }
            

            nuflavor = primary1->GetNuFlavor(settings1);
            nu_nubar = primary1->GetNuNuBar(nuflavor, settings1);

        
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

            Nu_temp = new Interaction (pnu, nuflavor, nu_nubar, n_interactions, icemodel, detector, settings1, primary1, signal, sec1 );        

            if(Nu_temp->sigma_err==1){
               // only if getting sigma was successful
               Nu_Interaction.push_back(*Nu_temp);  // for the first interaction
               delete Nu_temp;
            }
            else{
               //otherwise, just delete Nu_temp, but don't put it into Nu_Interaction
               delete Nu_temp;
               // and warn!
                cout<<"Warning! sigma_err is "<<Nu_temp->sigma_err<<endl;
                cout<<"This means the cross section was not found correctly!"<<endl;
                cout<<"Nu_Interaction will be empty!"<<endl;
            }

            inu_thrown++;
        } 
    }
    if (Event_type == 10) { // if only arbitrary events exist
        
  
        pnu = 0;
        nuflavor = "";
        nuflavorint = 0;
        nu_nubar = 0;
        
        Interaction *Nu_temp;
        // Report *report_tmp;

        Nu_temp = new Interaction (settings1, detector, icemodel, primary1, signal );
        // report_tmp = new Report(detector ,settings1);
        
        Nu_Interaction.push_back(*Nu_temp);  // for the first interaction
        // test_report.push_back(*report_tmp);

        // Initialize particle birth time to 0 (only relevant for multi-interaction events)
        interactions_birth_time.clear();
        interactions_birth_time.push_back(0.0);
        first_vertex_idx = inu_thrown;

        delete Nu_temp;

    }
    //Creating pulser event type that will be modelled after arbitrary event type (EVENT_TYPE=10). - JCF 4/6/2023
    if (Event_type == 11 or Event_type == 12) { // if only arbitrary events exist
        
  
        pnu = 0;
        nuflavor = "";
        nuflavorint = 0;
        nu_nubar = 0;
        
        Interaction *Nu_temp;

        Nu_temp = new Interaction (settings1, detector, icemodel, primary1, signal );
        
        Nu_Interaction.push_back(*Nu_temp);  // for the first interaction

        // Initialize particle birth time to 0 (only relevant for multi-interaction events)
        interactions_birth_time.clear();
        interactions_birth_time.push_back(0.0);
        first_vertex_idx = inu_thrown;

        delete Nu_temp;
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
        Event_type = 0;
    }
    else if (settings1->EVENT_TYPE == 10){
        Event_type = 10;
    }
    else if (settings1->EVENT_TYPE == 11){
        Event_type = 11;
    }
    else if (settings1->EVENT_TYPE == 12){
        Event_type = 12;
    }    
    else {
        cout<<"Currently, only neutrino (EVENT_TYPE=0), arbitrary (EVENT_TYPE=10), and pulser (EVENT_TYPE=11 or 12) events possible!"<<endl;
        cout<<"Change Evt_type from "<<settings1->EVENT_TYPE<<" to 0"<<endl;
        Event_type = 0;
    }
}
