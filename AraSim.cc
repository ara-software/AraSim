// ROOT includes
#include "TFile.h"
#include "TRandom3.h" 
#include "TTree.h"

// AraSim includes
//vector and position must be first
#include "Vector.h"
#include "Position.h"

#include "Constants.h"
#include "counting.hh"
#include "Detector.h"
#include "EarthModel.h"
#include "Efficiencies.h"
#include "Event.h"
#include "IceModel.h"
#include "Primaries.h"
#include "Ray.h"
#include "Report.h"
#include "RaySolver.h"
#include "secondaries.hh"
#include "Settings.h"
#include "signal.hh"
#include "Spectra.h"
#include "Tools.h"
#include "Trigger.h"
#include "Birefringence.h"

using namespace std;

#ifdef ARA_UTIL_EXISTS
    #include "UsefulIcrrStationEvent.h"
    ClassImp(UsefulIcrrStationEvent);
    #include "UsefulAtriStationEvent.h"
    ClassImp(UsefulAtriStationEvent);
#endif

class EarthModel; //class

void test();
void save_event_data(int, int, double*, int*, double*, double*, int*, Event*, Report*, TTree*, Counting*, Detector*, Settings*, Trigger*);
#ifdef ARA_UTIL_EXISTS
void save_useful_event(int, UsefulIcrrStationEvent*, UsefulAtriStationEvent*, double*, TTree*, Event*, Detector*, Report*, Settings*, Trigger*);
#endif

string outputdir="outputs";

int main(int argc, char **argv) {   // read setup.txt file
    
    Settings *settings1 = new Settings();

    cout<<"\n\tDefault values!"<<endl;
    cout<<"NNU : "<<settings1->NNU<<endl;
    cout<<"ICE_MODEL : "<<settings1->ICE_MODEL<<endl;
    cout<<"NOFZ : "<<settings1->NOFZ<<endl;
    cout<<"CONSTANTICETHICKNESS : "<<settings1->CONSTANTICETHICKNESS<<endl;
    cout<<"FIXEDELEVATION : "<<settings1->FIXEDELEVATION<<endl;
    cout<<"MOOREBAY : "<<settings1->MOOREBAY<<endl;
    cout<<"EXPONENT : "<<settings1->EXPONENT<<endl;
    cout<<"DETECTOR : "<<settings1->DETECTOR<<endl;

    string setupfile;
    string run_no;
    if (argc<2) { // no setup file input, use default
        setupfile = "setup.txt";
        cout<<"setupfile : "<<setupfile<<endl;
    }
    else if (argc > 1) { // read file!!
        setupfile = string( argv[1] );
        cout<<"setupfile : "<<setupfile<<endl;
    }
    if (argc > 2) { // read file!!
        run_no = string( argv[2] );
        cout<<"run number : "<<run_no<<endl;
    }
    if (argc > 3) { // read file!!
        outputdir = string( argv[3] );
        if(outputdir[outputdir.size()-1]=='/') outputdir=outputdir.substr(0,outputdir.size()-1); // make sure outputdir doesn't have a / at the end
        cout<<"outputdir : "<<outputdir<<endl;
    }

    settings1->ReadFile(setupfile);
    cout<<"Read "<<setupfile<<" file!"<<endl;

    int settings_compatibility_error = settings1->CheckCompatibilitiesSettings();
    if (settings_compatibility_error > 0) {
        cerr<<"There are "<< settings_compatibility_error<<" errors from settings. Check error messages."<<endl;
        return -1;
    }
 
    cout<<"\n\tNew values!"<<endl;
    cout<<"NNU : "<<settings1->NNU<<endl;
    cout<<"ICE_MODEL : "<<settings1->ICE_MODEL<<endl;
    cout<<"NOFZ : "<<settings1->NOFZ<<endl;
    cout<<"CONSTANTICETHICKNESS : "<<settings1->CONSTANTICETHICKNESS<<endl;
    cout<<"FIXEDELEVATION : "<<settings1->FIXEDELEVATION<<endl;
    cout<<"MOOREBAY : "<<settings1->MOOREBAY<<endl;
    cout<<"EXPONENT : "<<settings1->EXPONENT<<endl;
    cout<<"DETECTOR : "<<settings1->DETECTOR<<endl;
    cout<<"POSNU_RADIUS : "<<settings1->POSNU_RADIUS<<endl;
    cout << "EVENT_GENERATION_MODE: " << settings1->EVENT_GENERATION_MODE << endl;

    if (settings1->EVENT_GENERATION_MODE == 1){
        string evtfile = string(argv[argc - 1]);
        settings1->ReadEvtFile(evtfile);
        cout<<"Read "<< evtfile <<" file!"<<endl;
        cout << "EVID    NUFLAVORINT    NUBAR    PNU    CURRENTINT    IND_POSNU_R    IND_POSNU_THETA    IND_POSNU_PHI    IND_NNU_THETA    IND_NNU_PHI    ELAST" << endl;
        if (settings1->NNU == 0){
            // No events were read in from file, quit program
            cout<<"No events found in provided file. Exiting simulation."<<endl;
            return -1;
        }
    }

    // set gRandom as TRandom3 when settings1->RANDOM_MODE = 1
    if (settings1->RANDOM_MODE == 1) {
        // test TRandom3
        TRandom3 *test_randm3 = new TRandom3 (0);
        gRandom = test_randm3;
    } else {
        gRandom->SetSeed(settings1->SEED + atoi(run_no.c_str() ) );
            
    }
    cout<<"first random : "<<gRandom->Rndm()<<"\n";

    IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
    // IceModel inherits from EarthModel  

    cout<<endl;
    cout<<"Surface at (log:0, lat:0) : "<<icemodel->Surface(0., 0.)<<endl;
    cout<<"SurfaceAboveGeoid at (log:0, lat:0) : "<<icemodel->SurfaceAboveGeoid(0., 0.)<<endl;
    
    Detector *detector=new Detector(settings1, icemodel, setupfile ); // builds antenna array, 0 for testbed
    cout<<"end calling detector"<<endl;

    Birefringence *birefringence=new Birefringence(settings1);

    Trigger *trigger=new Trigger(detector, settings1); // builds the trigger  
    // Efficiencies *efficiencies=new Efficiencies(detector->getnRx(),outputdir); // keeps track of efficiencies at each stage of the simulation
    Efficiencies *efficiencies=new Efficiencies(100,outputdir); // keeps track of efficiencies at each stage of the simulation
    cout<<"called Efficiencies"<<endl;
    
    Spectra *spectra=new Spectra(settings1); // gets library (or whatever) of neutrino spectra
    cout<<"called Spectra"<<endl;

    Ray *ray = new Ray(); // construct Ray class
    cout<<"called Ray"<<endl;
    

    // test PickUnbiased in IceModel.
    Counting *count1 = new Counting();
    cout<<"called Counting"<<endl;

    Primaries *primary1 = new Primaries();
    cout<<"called Primaries"<<endl;

    int whichray = 0; // for test

    Event *event = new Event();
    cout<<"called Event"<<endl;

    Report *report = new Report();
    cout<<"called Evt"<<endl;

    // Build output files
    ofstream event_file;
    TFile *AraFile;
    if (argc > 2) {
        AraFile=new TFile((outputdir+"/AraOut."+setupfile.substr(setupfile.find_last_of("/")+1)+".run"+run_no+".root").c_str(),"RECREATE","ara");
    }
    else {
        AraFile=new TFile((outputdir+"/AraOut.root").c_str(),"RECREATE","ara");
    }
    TTree *AraTree=new TTree("AraTree","AraTree");    // for single entry
    TTree *AraTree2=new TTree("AraTree2","AraTree2"); //for many entries
    if (settings1->EVENT_GENERATION_MODE == 2){
        // Create file to save event list to 
        string output_file_name = (
            outputdir+"/AraOut."+
            setupfile.substr(setupfile.find_last_of("/")+1)+
            ".run"+run_no+".txt");
        event_file.open(output_file_name);
    }
    else{ 

        // Create tree for simulation outputs
        cout<<"assign AraFile, AraTrees"<<endl;

        AraTree->Branch("detector",&detector);
        cout<<"branch detector"<<endl;
        AraTree->Branch("icemodel",&icemodel);
        cout<<"branch icemodel"<<endl;
        AraTree->Branch("trigger",&trigger);
        cout<<"branch trigger"<<endl;
        AraTree->Branch("settings",&settings1);
        cout<<"branch settings"<<endl;
        AraTree->Branch("spectra",&spectra);
        cout<<"branch spectra"<<endl;
        AraTree2->Branch("event",&event);
        cout<<"branch Evt"<<endl;
        AraTree2->Branch("report",&report);
        cout<<"branch report"<<endl;

        cout<<"finish tree assign"<<endl;
    }

    RaySolver *raysolver = new RaySolver();
    cout<<"called RaySolver"<<endl;

    cout << "Make output file that is readable by AraRoot" << endl;

    #ifdef ARA_UTIL_EXISTS
    UsefulIcrrStationEvent *theIcrrEvent =0;
    UsefulAtriStationEvent *theAtriEvent =0;
    TTree *eventTree;
    double weight = 0.;
    if (settings1->EVENT_GENERATION_MODE != 2) {

        eventTree = new TTree("eventTree","Tree of ARA Events");
        eventTree->Branch("UsefulIcrrStationEvent", &theIcrrEvent);
        eventTree->Branch("UsefulAtriStationEvent",&theAtriEvent);
        eventTree->Branch("weight", &weight);

    }
    #endif


    cout<<"will call secondaries"<<endl;
    Secondaries *sec1 = new Secondaries (settings1);
    cout<<"will call signal"<<endl;
    Signal *signal = new Signal (settings1);
    signal->SetMedium(0);   // set medium as ice
    cout<<"finish calling secondaries and signal"<<endl;

    // before start looping events set noise values (this case, thermal)
    trigger->SetMeanRmsDiode(settings1, detector, report);

    if (settings1->NOISE_WAVEFORM_GENERATE_MODE == 0) {// noise waveforms will be generated for each evts
        trigger->ClearNoiseWaveforms();
    }

    // now in Trigger class, there will be meandiode, rmsdiode values for noise (we need this for trigger later)

    double max_dt = 0.; // max arrival time difference

    int Total_Global_Pass = 0;  // total global trigger passed number 
    double Total_Weight = 0.;
    double Total_Probability = 0.;

    double x_V[settings1->NFOUR/2];
    double y_V[settings1->NFOUR/2];

    double xbin[settings1->DATA_BIN_SIZE];
    for (int i=0; i<settings1->DATA_BIN_SIZE; i++) {
        xbin[i] = i;
    }

    cout<<"powerthreshold : "<<trigger->powerthreshold<<endl;

    int check_station_DC;

    ofstream TrigWind;
    TrigWind.open("outputs/TrigWindowStudy.txt");
                
    Total_Global_Pass = 0;
    cout<<"begin looping events!!"<<endl;

    double pre_posnu_x;
    double pre_posnu_y;
    double pre_posnu_z;

    double cur_posnu_x;
    double cur_posnu_y;
    double cur_posnu_z;

    cout << "Calpulser_on: " << settings1->CALPULSER_ON << endl;

    // test Detector set correctly
    cout<<"number of stations : "<<detector->params.number_of_stations << endl;
    cout<<"total number of antennas : "<<detector->params.number_of_antennas << endl;
    int ch_count = 0;
    for (int i=0; i<detector->params.number_of_stations; i++) {
        for (int j=0; j<detector->stations[i].strings.size(); j++) {
            for (int k=0; k<detector->stations[i].strings[j].antennas.size(); k++) {
                ch_count++;
                cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"] no_ch:"<<ch_count<<endl;
            }
        }
    }
         
    // check if settings have to compatibility problems
    // if there's any, stop AraSim
    settings_compatibility_error = settings1->CheckCompatibilitiesDetector(detector);
    if (settings_compatibility_error > 0) {
        cerr<<"There are "<< settings_compatibility_error<<" errors from settings after Detector class instance is initialized. Check error messages."<<endl;
        return -1;
    }
        
    #ifndef ARA_UTIL_EXISTS
        if (settings1->DETECTOR == 3 && settings1->READGEOM == 1){
            cerr << "ERROR::InstalledStation geometry not available without AraRoot installation!" << endl;
            return -1;
        }
    #endif

    // reset accumulative trig search bin info 
    settings1->ACCUM_TRIG_SEARCH_BINS_STATION0 = 0.;

    int nuLimit =0;
    if (settings1->EVENT_GENERATION_MODE == 1){ //event mode read in different single events
        nuLimit = settings1->NNU;
    }
    else if (settings1->ONLY_PASSED_EVENTS == 1){
        nuLimit = settings1->NNU_PASSED;
    }
    else {
        nuLimit = settings1->NNU;
    }
    int inu = 0;
    int Events_Thrown = 0;
    int Events_Passed = 0;

    while (inu < nuLimit){
        check_station_DC = 0;
        if ( settings1->DEBUG_MODE_ON==0 ) {
            std::cerr<<"*";
            if ( Events_Thrown%100 == 0 )
                cout<<"Thrown "<<Events_Thrown<<endl;
        }

        event = new Event ( settings1, spectra, primary1, icemodel, detector, signal, sec1, Events_Thrown );
        if(event->Nu_Interaction.size()<1){
            // If for some reason no interactions were placed into the event holder, continue.
            // This should be exceedingly rare, but could happen if something
            // goes wrong with interaction generation.
            // BAC added this specifically because if the neutrino is generated
            // with a higher energy than is available in the sigma parameterization
            // then no interaction will be added and Nu_Interaction will be empty.
            // But also, it's bad practice to assume the size of a vector anyway
            // so this is just better generically.
            cout<<"Warning! The interaction vector is empty!"<<endl;
            cout<<"Continuing on to the next event."<<endl;
            if (settings1->EVENT_GENERATION_MODE ==1){
                // If reading in events from a list, make sure you move to next event 
                Events_Thrown++;
                inu++;
            }
            continue;
        }
        event->inu_passed = -1;

        if ( settings1->EVENT_GENERATION_MODE == 2 ){
            // Write event lists to file, no simulation
            for (int interaction_i=0; interaction_i<event->Nu_Interaction.size(); interaction_i++){
                
                event_file << inu << " "; // EVID    
                event_file << event->nuflavorint << " "; // NUFLAVORINT       
                event_file << event->nu_nubar << " "; // NUBAR       
                event_file << event->pnu << " "; // PNU      
                event_file << event->Nu_Interaction[interaction_i].currentint    << " "; // CURRENTINT       
                event_file << event->Nu_Interaction[interaction_i].posnu.R()     << " "; // IND_POSNU_R      
                event_file << event->Nu_Interaction[interaction_i].posnu.Theta() << " "; // IND_POSNU_THETA       
                event_file << event->Nu_Interaction[interaction_i].posnu.Phi()   << " "; // IND_POSNU_PHI      
                event_file << event->Nu_Interaction[interaction_i].nnu.Theta()   << " "; // IND_NNU_THETA      
                event_file << event->Nu_Interaction[interaction_i].nnu.Phi()     << " "; // IND_NNU_PHI       
                event_file << event->Nu_Interaction[interaction_i].elast_y       << endl; // ELAST

                inu++;
                Events_Thrown++;
                event->delete_all();

            } // end add event to file

            continue; // Skip the simulation steps and move to next event
        }
                
        report = new Report(detector, settings1);
                        
        #ifdef ARA_UTIL_EXISTS
            theIcrrEvent = new UsefulIcrrStationEvent();
            theAtriEvent = new UsefulAtriStationEvent();
        #endif


        // go further only if we picked up usable posnu
        if (event->Nu_Interaction[0].pickposnu>0) {

            //--------------------------------------------------
            // cout<<"inu : "<<inu<<endl;
            // cout<<"event->pnu : "<<event->pnu<<endl;
            // cout<<"posnu : ";
            // event->Nu_Interaction[0].posnu.Print();
            // cout<<"nnu : ";
            // event->Nu_Interaction[0].nnu.Print();
            // cout<<"event->n_interactions : "<<event->n_interactions<<endl;
            // cout<<"nu_flavor : "<<event->nuflavor<<endl;
            // cout<<"event->Nu_Interaction[0].vmmhz1m[0] : "<<event->Nu_Interaction[0].vmmhz1m[0]<<endl;
            // cout<<"pickposnu : "<<event->Nu_Interaction[0].pickposnu<<endl;
            //-------------------------------------------------- 

            // Determine if the simulation is in debug mode or not
            int debugmode = 0;
            if (settings1->DEBUG_MODE_ON == 1 && Events_Thrown < settings1->DEBUG_SKIP_EVT) {
                debugmode = 1;  // skip most of computation intensive processes if debugmode == 1
            }
            else if (settings1->DEBUG_MODE_ON == 1 && Events_Thrown >= settings1->DEBUG_SKIP_EVT){
                cout << Events_Thrown << " " << endl;
            }

            // While there is at least 1 station with enough waveform to perform a trigger check,
            //   build waveforms for this event and check the trigger over it. 
            //   If there is a part of any waveform that wasn't analyzed and exists beyond the
            //   deadtime of the station, say there is still at least one station to trigger check. 
            int stations_to_trigger_check = report->stations.size();
            for (int station_i=0; station_i<report->stations.size(); station_i++){
                report->stations[station_i].next_trig_search_init = trigger->maxt_diode_bin + settings1->NFOUR;  
            }
            while (stations_to_trigger_check > 0) {
                
                // Ray trace and calculate signal for every ray from every interaction to every antenna
                report->CalculateSignals(debugmode, birefringence, detector, event, icemodel, raysolver, settings1, signal);

                // Combine all signal waveforms, add noise, perform trigger check on each station
                for (int station=0; station < detector->params.number_of_stations; station++) {
                    report->BuildAndTriggerOnWaveforms(
                        debugmode, station, Events_Thrown, 
                        report->stations[station].next_trig_search_init, detector, event, settings1, trigger);
                }

                // Save the event object to a temporary object in case we have more data to analyze
                Event *event_save = new Event(*event);

                // Count how many stations have waveform left to analyze (indicated by whether 
                //   `report->stations[station_i].next_trig_search_init` is equal to `-1` or not)
                stations_to_trigger_check = 0;
                int stations_next_trig_search_init[ detector->stations.size() ];
                for (int station_i=0; station_i<report->stations.size(); station_i++){
                    stations_next_trig_search_init[station_i] = report->stations[station_i].next_trig_search_init;
                    if (report->stations[station_i].next_trig_search_init != -1){
                        stations_to_trigger_check++;
                    }
                }

                // Save this event's event data, report data, and ARA_like data if ARA_UTIL_EXISTS
                save_event_data(
                    Events_Passed, inu, 
                    &max_dt, &Total_Global_Pass, &Total_Weight, &Total_Probability, &check_station_DC,
                    event, report, AraTree2,
                    count1, detector, settings1, trigger
                );
                #ifdef ARA_UTIL_EXISTS
                save_useful_event(
                    check_station_DC, 
                    theIcrrEvent, theAtriEvent, &weight, eventTree,
                    event,
                    detector, report, settings1, trigger);
                #endif

                // Ensure the calpulser didn't move randomly during the simulation
                if (settings1->CALPULSER_ON == 1) {
                    cur_posnu_x = event->Nu_Interaction[0].posnu.GetX();
                    cur_posnu_y = event->Nu_Interaction[0].posnu.GetY();
                    cur_posnu_z = event->Nu_Interaction[0].posnu.GetZ();
                    if (inu>0) {
                        if (pre_posnu_x==cur_posnu_x && pre_posnu_y==cur_posnu_y && pre_posnu_z==cur_posnu_z) {
                        }
                        else cout<<"posnu location changed!"<<endl;
                    }
                    pre_posnu_x = event->Nu_Interaction[0].posnu.GetX();
                    pre_posnu_y = event->Nu_Interaction[0].posnu.GetY();
                    pre_posnu_z = event->Nu_Interaction[0].posnu.GetZ();
                }

                // If there are stations to check for more trigger, overwrite 
                //   the probably, partitially cleared report and event objects
                //   with the saved full report and event objects
                if (stations_to_trigger_check >= 1){
                    delete report;
                    delete event;
                    
                    report = new Report(detector, settings1);
                    event = new Event(*event_save); 
                }

                for (int station_i=0; station_i<report->stations.size(); station_i++){
                    report->stations[station_i].next_trig_search_init = stations_next_trig_search_init[station_i];
                }

                // Delete the temporary report and event objects
                delete event_save;

            }

        } // if pickposnu > 0
        else {
            report->delete_all();
            event->delete_all();
        }

        if (settings1->EVENT_GENERATION_MODE == 1){
            inu++;
        }
        else if (settings1->ONLY_PASSED_EVENTS == 1){
            if (check_station_DC > 0){
                inu++;
            }
        }
        else {
            inu++;
        }
        if (check_station_DC > 0){
            Events_Passed++;
        }
        Events_Thrown++;

        delete event;
        delete report;
        #ifdef ARA_UTIL_EXISTS
            delete theIcrrEvent;
            delete theAtriEvent;
        #endif

    } // end loop over neutrinos

    settings1->NNU = Events_Thrown;
    settings1->NNU_PASSED = Total_Global_Pass;

    if (settings1->EVENT_GENERATION_MODE == 2){
        event_file.close();
        return 0;
    }

    // }// end trigger window loop
    TrigWind.close();

    ofstream weight_file;
    if (argc == 3) {
        weight_file.open(("./weight_output/weight_"+setupfile+".run"+run_no).c_str());
    }
    else if (argc >3 ){ // add the subdirectory for outputs
        weight_file.open((outputdir+"/weight_output/weight_"+setupfile+".run"+run_no).c_str());
    }
    else {
        weight_file.open(("./weight_output/weight_"+setupfile).c_str());
    }


    cout<<" end loop"<<endl;
    cout << "Total Events Thrown: " << Events_Thrown << endl;
    cout<<"Total_Global_Pass : "<<Total_Global_Pass<<endl;
    cout<<"Total_Weight : "<<Total_Weight<<endl;
    cout<<"Total_Probability : "<<Total_Probability<<endl;
                             
    if (settings1->INTERACTION_MODE==1) {
        weight_file << "Total_Weight="<<Total_Weight<<endl;
    }
    else if (settings1->INTERACTION_MODE==0) {
        weight_file << "Total_Probability="<<Total_Probability<<endl;
    }

    cout<<"weight bin values : ";
    for (int i=0; i<count1->NBINS-1; i++) {
        cout<<count1->eventsfound_binned[i]<<", ";
        weight_file << count1->eventsfound_binned[i]<<" ";
    }
    cout<<count1->eventsfound_binned[count1->NBINS-1];
    weight_file << count1->eventsfound_binned[count1->NBINS-1]<<"\n";
    weight_file.close();
    cout<<"\n\n";


    // if using picknear_cylinder method
    if (settings1->INTERACTION_MODE==1) {
        double IceVolume;
        IceVolume = PI * (settings1->POSNU_RADIUS) * (settings1->POSNU_RADIUS) * icemodel->IceThickness( detector->stations[0] );
        cout << "Radius: " << settings1->POSNU_RADIUS << " [m]" << endl;
        cout<<"IceVolume : "<<IceVolume<<endl;

        double Veff_test_we; // effective volume water equivalent
        double Veff_test; // effective volume ice

        // error bar for weight
        double error_plus = 0;
        double error_minus = 0;
        Counting::findErrorOnSumWeights( count1->eventsfound_binned, error_plus, error_minus );

        Veff_test_we = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * Total_Weight / (double)(settings1->NNU);
        Veff_test = IceVolume * 4. * PI * Total_Weight / (double)(settings1->NNU);
        error_plus = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * error_plus / (double)(settings1->NNU);
        error_minus = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * error_minus / (double)(settings1->NNU);

        cout<<"test Veff(ice) : "<<Veff_test<<" m3sr, "<<Veff_test*1.E-9<<" km3sr"<<endl;
        cout<<"test Veff(water eq.) : "<<Veff_test_we<<" m3sr, "<<Veff_test_we*1.E-9<<" km3sr"<<endl;
        cout<<"And Veff(water eq.) error plus : "<<error_plus*1.E-9<<" km3sr and error minus : "<<error_minus*1.E-9<<" km3sr"<<endl;
    }


     // if using picknear_sphere method
     else if (settings1->INTERACTION_MODE==0) {

        double IceArea;
        IceArea = PI * (settings1->POSNU_RADIUS) * (settings1->POSNU_RADIUS);
        cout << endl;
        cout<<"IceArea : "<< IceArea <<endl;

        double Aeff;
        Aeff = IceArea * Total_Probability / (double)(settings1->NNU);
        cout << "Aeff : " << Aeff << " [m^2]" << endl;

        // error bar for weight
        double error_plus = 0;
        double error_minus = 0;
        Counting::findErrorOnSumWeights( count1->eventsfound_binned, error_plus, error_minus );

        // account all factors to error
        error_plus = IceArea * error_plus / (double)(settings1->NNU);
        error_minus = IceArea * error_minus / (double)(settings1->NNU);

        cout<<"And Aeff error plus : "<<error_plus<<" and error minus : "<<error_minus<<" [m^2]"<<endl;

        // Aeff*sr
        cout<<"and Aeff*sr values are"<<endl;
        cout << "Aeff*sr : " << Aeff * 4.* PI << " [m^2sr]" <<", "<< Aeff * 4.* PI *1.e-6<<" [km^2sr]"<< endl;
    }

    // remove noisewaveform info if DATA_SAVE_MODE == 2
    // remove noisewaveform info if DATA_SAVE_MODE is not 0
    if (settings1->DATA_SAVE_MODE != 0) {
        trigger->v_noise_timedomain.clear();
        trigger->v_noise_timedomain_diode.clear();
    }
    if (settings1->DATA_SAVE_MODE == 2) {// in DATA_SAVE_MODE==2, remove noise spectrum before Rayleigh dist.
        trigger->Vfft_noise_before.clear();
    }
    
    AraTree->Fill();  // fill tree for one entry
    AraFile->Write();
    AraFile->Close();

    efficiencies->summarize(); // summarize the results in an output file  

    double freq[detector->GetFreqBin()], Filter[detector->GetFreqBin()];
    double Filter_E[detector->GetFreqBin()];

    for (int i=0; i<detector->GetFreqBin(); i++) {
        freq[i] = detector->GetFreq(i);    // in Hz
        Filter[i] = detector->GetFilterGain(i);    // in dB
        Filter_E[i] = pow(10., (detector->GetFilterGain(i))/20.);
    }

    cout<<"max_dt : "<<max_dt<<endl;
    cout<<"rmsdiode= "<<trigger->GetAntNoise_diodeRMS(0, settings1)<<endl;

    delete raysolver;
    delete icemodel;
    delete efficiencies;
    delete ray;
    delete detector;
    delete settings1;
    delete count1;
    delete primary1;
    delete trigger;
    delete spectra;
    delete sec1;
    delete signal;

    cout<<"outputdir= "<<outputdir<<endl;

    // Please do not delete this test line
    // Please leave 'test(); return 0;' as the last lines in 
    // AraSim.cc before '} //end main'
    // These test lines are used to verify that AraSim completed properly.

    test();
    return 0;   
} //end main


void test() {
    cout << "test is " << 0 << "\n";
}

void save_event_data(
    int Events_Passed, int inu, 
    double *max_dt, int *Total_Global_Pass, double *Total_Weight, double *Total_Probability, int *check_station_DC,
    Event *event, Report *report, TTree *AraTree2,
    Counting *count1, Detector *detector, Settings *settings1, Trigger *trigger
){
    // Save the report and event branches of the requested simulated event to the output TTree

    report->clear_useless(settings1);   // to reduce the size of output AraOut.root, remove some information
    report->ClearUselessfromConnect(detector, settings1, trigger);
    for(int i=0;i<event->Nu_Interaction.size(); i++){
        event->Nu_Interaction[i].clear_useless(settings1);
    }

    for (int i=0; i<detector->params.number_of_stations; i++) {
        
        // check the total global trigger passed
        if (report->stations[i].Global_Pass) {

            event->inu_passed = Events_Passed;  
                     
            if (*max_dt < report->stations[i].max_arrival_time - report->stations[i].min_arrival_time)
                *max_dt = report->stations[i].max_arrival_time - report->stations[i].min_arrival_time;

            cout<<"\nGlobal_Pass : "<<report->stations[i].Global_Pass<<" evt : "<<inu
                <<" added weight : "<<event->Nu_Interaction[0].weight<<endl;

            if ( *check_station_DC == 0) { // count trigger pass only once per event
                        
                *Total_Global_Pass += 1;
                *Total_Weight += event->Nu_Interaction[0].weight;
                *Total_Probability += event->Nu_Interaction[0].probability;
                        
                // test increment weight
                if (settings1->INTERACTION_MODE==1) {
                    count1->incrementEventsFound( event->Nu_Interaction[0].weight, event );
                }
                else if (settings1->INTERACTION_MODE==0) {
                    count1->incrementEventsFound( event->Nu_Interaction[0].probability, event );
                }
                else if (settings1->INTERACTION_MODE==3) {
                    count1->incrementEventsFound( event->Nu_Interaction[0].probability, event );
                }
                else if (settings1->INTERACTION_MODE==4) {
                    count1->incrementEventsFound( event->Nu_Interaction[0].weight, event );
                }    

            }
            *check_station_DC += 1;

        }
    }

    settings1->ACCUM_TRIG_SEARCH_BINS_STATION0 += report->stations[0].total_trig_search_bin;

    // test FILL_TREE_MODE
    if (settings1->FILL_TREE_MODE==0) { // fill event event  
        AraTree2->Fill();   //fill interaction every events
    }
    else if (settings1->FILL_TREE_MODE==1) { // fill only usable posnu event 
        if (event->Nu_Interaction[0].pickposnu>0) {
            AraTree2->Fill();   //fill interaction every events
        }
    }
    else if (settings1->FILL_TREE_MODE==2) { // fill only triggered event    
        if (*check_station_DC>0) {
            AraTree2->Fill();   //fill interaction every events
        }
    }

}


#ifdef ARA_UTIL_EXISTS
void save_useful_event(
    int check_station_DC, 
    UsefulIcrrStationEvent *theIcrrEvent, UsefulAtriStationEvent *theAtriEvent, double *weight, TTree *eventTree,
    Event *event, 
    Detector *detector, Report *report, Settings *settings1, Trigger *trigger){
    // Save simulated event data in a format similar to how ARA detectors save their events

    // Extract the number of channels in each detector
    for (int i=0; i<detector->params.number_of_stations; i++) {
        if (settings1->DATA_LIKE_OUTPUT != 0){
            if (settings1->DETECTOR == 3 && i == 0)
                // The testbed has 14 channels
                { theIcrrEvent->numRFChans = 14; }
            else if (settings1->DETECTOR == 4 && settings1->DETECTOR_STATION == 0)
                // The testbed has 14 channels
                { theIcrrEvent->numRFChans = 14; }
            else { 
                // All other stations have 16 channels
                theAtriEvent->fNumChannels = 20; // Includes 4 surface channels
                theIcrrEvent->numRFChans = 16; 
            }
        }
    }
    if (settings1->DATA_LIKE_OUTPUT !=0){
        int stationID;
        int stationIndex;
        if (settings1->DETECTOR == 4){
            // Traditional ARA1-5 Detectors
            stationID = settings1->DETECTOR_STATION;
            stationIndex = 0;
        }
        else if (settings1->DETECTOR == 5){
            // Phased Array
            stationID = 6;
            stationIndex = 0;
        } else {
            // Other detectors, e.g. the Testbed
            stationID = 0;
            stationIndex = 0;
        }

        if (report->stations[stationIndex].Global_Pass) {
            cout << endl << "Making useful event" << endl;
            report->MakeUsefulEvent(detector, settings1, trigger, stationID, stationIndex, theAtriEvent);
        }
        *weight = event->Nu_Interaction[0].weight;
    }   

    // Determine if the current event should be saved based on simulation settings
    if (settings1->FILL_TREE_MODE==0) { 
        // Save data for all events
        if (settings1->DATA_LIKE_OUTPUT==2) {
            // save all events whether they passed the trigger or not
            eventTree->Fill();
        }
        else if (settings1->DATA_LIKE_OUTPUT==1) {
            // only save events that passed the trigger
            if ( check_station_DC > 0 ) {
                eventTree->Fill();
            }
        }
    }
    else if (settings1->FILL_TREE_MODE==1) { 
        // Only save events if their location was "usable"
        if (event->Nu_Interaction[0].pickposnu>0) {
            if (settings1->DATA_LIKE_OUTPUT==2) {
                // save all events whether they passed the trigger or not
                eventTree->Fill();
            }
            else if (settings1->DATA_LIKE_OUTPUT==1) {
                // only save events that passed the trigger
                if ( check_station_DC > 0 ) {
                    eventTree->Fill();
                }
            }
        }
    }
    else if (settings1->FILL_TREE_MODE==2) { 
        // Only save events that triggered
        if (check_station_DC>0) {
            if (settings1->DATA_LIKE_OUTPUT==2) {
                // save all events whether they passed the trigger or not
                eventTree->Fill();
            }
            else if (settings1->DATA_LIKE_OUTPUT==1) {
                // only save events that passed the trigger
                if ( check_station_DC > 0 ) {
                    eventTree->Fill();
                }
            }
        }
    }

}
#endif
