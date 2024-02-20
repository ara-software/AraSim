#include "/cvmfs/ara.opensciencegrid.org/trunk/centos7/source/libRootFftwWrapper/include/FFTtools.h"

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

using namespace std;

#ifdef ARA_UTIL_EXISTS
    #include "UsefulIcrrStationEvent.h"
    ClassImp(UsefulIcrrStationEvent);
    #include "UsefulAtriStationEvent.h"
    ClassImp(UsefulAtriStationEvent);
#endif

// #include "/users/PAS0654/jflaherty13/source/AraRoot/AraRoot_build/include/AraGeomTool.h"
// #include "/users/PAS0654/jflaherty13/source/AraRoot/AraRoot_build/include/RayTraceCorrelator.h"
// #include "/users/PAS0654/jflaherty13/source/AraRoot/AraRoot_build/include/UsefulAtriStationEvent.h"

//TODO: Have outgoing pointer equal the incoming pointer by using a filler function to copy the original information, then I replace the voltage info with my deconvolved voltage.
UsefulAtriStationEvent *usefulAtriEvPtr;
UsefulAtriStationEvent *usefulAtriEvPtrOut;

int main(int argc, char **argv)
{
    if(argc<6) {
        std::cout << "Usage\n" << argv[0] << " <station> <config> <runnum> <input root file> <input reco file> <output_dir> \n";
        std::cout << "e.g.\n" << argv[0] << " 2 6 AraOut.root recangle_out_run<runnum>.root output/\n";
        return 0;
    }
    
    double interpV = 0.4;
    double interpH = 0.625;
    
    //Import AraRoot file
    printf("Opening root file...\n");
    TFile *fp = TFile::Open(argv[4]);
    if(!fp) { std::cerr << "Can't open file\n"; return -1; }
    printf("Root File opened!\n");
    
    //Import eventTree
    TTree *eventTree = (TTree*) fp->Get("eventTree");
    printf("Event tree opened!\n");
    if(!eventTree) { std::cerr << "Can't find eventTree\n"; return -1; }
    Long64_t numEntries=eventTree->GetEntries();
    cout << "eventTree has " << numEntries << " entries." << endl;
    RawAtriStationEvent *rawAtriEvPtr=0;
    
    // Check if sim or real data file by checking for existence of AraTree
    TTree *simSettingsTree;
    simSettingsTree=(TTree*) fp->Get("AraTree");
    bool dataLike = false;
    //data like
    if(!simSettingsTree) { 
        dataLike = true;
        std::cerr << "Can't find AraTree.  Importing as real data.\n";
        eventTree->SetBranchAddress("event",&rawAtriEvPtr);
        double weight = 1;
    }
    // sim like
    else {
        std::cerr << "AraTree exists.  Importing as simulated data.\n";
        eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
        double weight;
        eventTree->SetBranchAddress("weight", &weight);   
    }
        
    //Import vertex reco file
    printf("Opening reco file...\n");
    TFile *fp2 = TFile::Open(argv[5]);
    if(!fp2) { std::cerr << "Can't open file\n"; return -1; }
    printf("Reco File opened!\n");
    TTree *vertexReco = (TTree*) fp2->Get("vertexReco");
    double reco_arrivalThetas[16];
    double reco_arrivalPhis[16];
    double cutoffTime[16];
    
    // Testing using the true rf angles
    // vertexReco->SetBranchAddress("true_arrivalThetas", reco_arrivalThetas);
    // vertexReco->SetBranchAddress("true_arrivalPhis", reco_arrivalPhis);   
    // end testing
    vertexReco->SetBranchAddress("reco_arrivalThetas", reco_arrivalThetas);
    vertexReco->SetBranchAddress("reco_arrivalPhis", reco_arrivalPhis);
    vertexReco->SetBranchAddress("cutoffTime", cutoffTime);  
    
    printf("Vertex Reco tree opened!\n");
    
    printf("------------------\n");
    printf("Input files loaded.  Setting up detector stuff.\n");
    printf("------------------\n");
    
    string setupfile;
    setupfile = "SETUP/setup_variablePsi.txt";
    Settings *settings1 = new Settings();
    settings1->ReadFile(setupfile); 
    IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
    Detector *detector = new Detector(settings1, icemodel, setupfile);  
    Report *report = new Report(detector, settings1);
    
    settings1->NFOUR = 4096;
    
    cout << "Settings->TIMESTEP = " << settings1->TIMESTEP << endl;
    
    printf("------------------\n");
    printf("Make Output Files\n");
    printf("------------------\n");

    char outfile_name[400];
    sprintf(outfile_name, "%s/deconvolvedWaveforms_run_%s.root", argv[6], argv[3]);    
    
    std::cout<<"Output name is "<<outfile_name<<std::endl;
    TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
    if(!fpOut){ std::cerr<<"Cannot open output file "<<fpOut<<std::endl; return -1; }    
    
     TTree *outTree = new TTree("eventTree", "eventTree");
    outTree->Branch("UsefulAtriStationEvent", &usefulAtriEvPtrOut);
    
    //Need to grab lengths of voltage and time arrays from eventTree to initialize the branches in the outfile.
    Int_t fNumChannels; ///< The number of channels
    std::map< Int_t, std::vector <Double_t> > fTimesOut; ///< The times of samples
    std::map< Int_t, std::vector <Double_t> > fVoltsOut; ///< The voltages of samples    
    
    
    //Loop over events
    for(Long64_t event=0;event<numEntries;event++) {
        fp->cd();
        eventTree->GetEntry(event);
        vertexReco->GetEntry(event);
    
        std::cout<<"Looking at event number "<<event<<std::endl;
        
        if (dataLike) {
            cout << "Triggering datalike condition." << endl;
            delete usefulAtriEvPtr;         
            usefulAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);               
        }

        int vertexRecoElectToRFChan[] = {14,2,6,10,12,0,4,8,15,3,7,11,13,1,5,9};
        
        

        for(int i=0; i<16; i++){
            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);  //This is where the code breaks for real data.
            
            //Save initial and final time for truncating the padded arrays before output.
            double timeStart = gr->GetX()[0];
            double timeEnd = gr->GetX()[gr->GetN()-1];
            
            //Interpolate graph to 0.5 ns resolution
            gr = FFTtools::getInterpolatedGraph(gr,0.5);
            
            //Pad waveform to a factor of two. - JCF 9/27/2023
            if (gr->GetN() < settings1->NFOUR/2) {
                gr = FFTtools::padWaveToLength(gr, settings1->NFOUR/2);
            }
            // Padding 
            int waveform_bin = gr->GetN();
            
            
            double heff_lastbin;
            double freq_lastbin;
            double time[waveform_bin];
            double voltage[waveform_bin];
            double volts_forint[settings1->NFOUR / 2];
            double T_forint[settings1->NFOUR / 2];
            
            //TODO: This init_T isn't dynamic to the imported data.  Should make have it defined based on the input waveform.
            double init_T = settings1->TIMESTEP *-1.e9 *((double) settings1->NFOUR / 4);    // locate zero time at the middle and give random time shift       
            
            for(int k=0; k<waveform_bin; k++){
                time[k] = gr->GetX()[k];
                voltage[k] = gr->GetY()[k];
            }            
            delete gr;
            
            // for (int m = 0; m < settings1->NFOUR / 2; m++)
            for (int m = 0; m < 2048; m++)
            {
                T_forint[m] = -512 + m*0.5;   // in ns
            }
            
            //Importing the cutoff time between spicecore peaks
            double cutoffTimeChannel;

            if (!cutoffTime) {
                cutoffTimeChannel = time[waveform_bin-1];
            } else {
                cutoffTimeChannel = cutoffTime[i];
            }            
            //TODO: Add step that centers the waveform about zero in time, for purposes of the fourier transform.  Then save this shift and reapply it to restore the time-domain information after the InvFFT.            
            double timeshift = time[waveform_bin/2];
            cout << "timeshift = " << timeshift << endl;
    
            double freq_tmp, heff, antenna_theta, antenna_phi;  // values needed for apply antenna gain factor and prepare fft, trigger
            
            cout << "Importing angles." << endl;
            if (dataLike) {
                //Import RF angles use RF channel mapping
                antenna_theta = reco_arrivalThetas[i]*180/PI;
                antenna_phi = reco_arrivalPhis[i]*180/PI;
            }
            else {
                //Import RF angles using electric channel mapping
                antenna_theta = reco_arrivalThetas[vertexRecoElectToRFChan[i]]*180/PI;
                antenna_phi = reco_arrivalPhis[vertexRecoElectToRFChan[i]]*180/PI;
            }
            cout << "antenna_theta = " << antenna_theta << endl;
            cout << "antenna_phi = " << antenna_phi << endl;                    
            
            //Calculate polarization vector that inverts the polarization factor (makes dot products equal to one)
            double newPol_vectorX = -sin(antenna_phi*PI/180);
            double newPol_vectorY = cos(antenna_phi*PI/180);
            double newPol_vectorZ = -1/sin(antenna_theta*PI/180);
        
            //Testing using the actual polarization vector
            // double psi = argv[3]*PI/180;
            // double newPol_vectorX = -cos(psi)*cos(antenna_theta*PI/180)*cos(antenna_phi*PI/180) + sin(psi)*sin(antenna_phi*PI/180);
            // double newPol_vectorY = -cos(psi)*cos(antenna_theta*PI/180)*sin(antenna_phi*PI/180) - sin(psi)*cos(antenna_phi*PI/180);
            // double newPol_vectorZ = cos(psi)*sin(antenna_theta*PI/180);            

            Vector Pol_vector = Vector(newPol_vectorX, newPol_vectorY, newPol_vectorZ);

            double dF_Nnew;

            double nice = 1.79;


            int pol_ant;
            int gain_ch_no = i;          
            double Pol_factor;
            
            if (i < 8) {
                pol_ant=0;
            } else {
                pol_ant=1;
            }
            
            double dT_forfft = time[1] - time[0];
        
            int Ntmp = settings1->TIMESTEP *1.e9 / dT_forfft;
            
            int Nnew = 1;         
            while (Ntmp > 1)
            {
                Ntmp = Ntmp / 2;
                Nnew = Nnew *2;             
            }
            Nnew = Nnew * settings1->NFOUR / 2;            
            
            //Stealing antenna and electronic response steps from AraSim, but applying the inverse functions instead.
            
            double V_forfft[Nnew];
            double T_forfft[Nnew];
            
            for (int n = 0; n < Nnew; n++)
            {

                // make Tarray, Earray located at the center of Nnew array

                T_forfft[n] = time[waveform_bin / 2] - (dT_forfft *(double)(Nnew / 2 - n));
                // T_forfft[n] = time[waveform_bin / 2] - (dT_forfft *(double)(Nnew / 2 - n)) - timeshift;  //TODO: Apply time shift to center array in time domain.

                if ((n >= Nnew / 2 - waveform_bin / 2) &&
                    (n < Nnew / 2 + waveform_bin / 2))
                {
                    V_forfft[n] = voltage[n - (Nnew / 2 - waveform_bin / 2)];
                }
                else
                    V_forfft[n] = 0.;

            }            
            
            // get spectrum with zero padded WF
            Tools::realft(V_forfft, 1, Nnew); 
            
            dF_Nnew = 1. / ((double)(Nnew) *(dT_forfft) *1.e-9);    // in Hz

            freq_tmp = dF_Nnew *((double) Nnew / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

            freq_lastbin = freq_tmp;
         
            for (int n = 0; n < Nnew / 2; n++)
            // for (int n = 0; n < settings1->NFOUR / 2; n++)            
            {
                freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq
                heff_lastbin = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                        antenna_theta, antenna_phi, pol_ant),
                    freq_tmp, nice);             

                heff = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                        antenna_theta, antenna_phi, pol_ant),
                    freq_tmp, nice);
                // invert entire elect chain gain, phase
                if (n > 0)
                {                
                    report->InvertElect_Tdomain(freq_tmp *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no);
                }
                else
                {
                    report->InvertElect_Tdomain_FirstTwo(freq_tmp *1.e-6, freq_lastbin *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no);
                }
                if (n > 0)
                {
                    report->InvertAntFactors_Tdomain(detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, pol_ant),
                       heff, Pol_vector, pol_ant, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi);                   
                }
                else
                {
                    report->InvertAntFactors_Tdomain_FirstTwo(heff, heff_lastbin, Pol_vector, pol_ant, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi);

                }
                
                // Quick and dirty hack to filter out frequencies above 850 MHz and below 100 MHz.             
                if (freq_tmp > 850.*1.e6 or freq_tmp < 100.*1.e6) {
                    V_forfft[2*n] = 0;
                    V_forfft[2*n+1] = 0;
                }                   
                //Apply homemade butterworth filter of the fourth order
                // double freqMin = 150*1e6;
                // double freqMax = 300*1e6;
                
                //Trying user inputted butterworth filter
                double freqMin = atof(argv[7])*1e6;
                double freqMax = atof(argv[8])*1e6;
  
                
                double weight = 1;  // Setting initial weight to one, then applying bandpass.  Weight is then multiplied by signal in this bin.
                int order = 8;
                weight /= sqrt(1 + TMath::Power(freqMin/freq_tmp, 4*order));
                weight /= sqrt(1 + TMath::Power(freq_tmp/freqMax, 4*order));
                V_forfft[2*n] *= weight;
                V_forfft[2*n+1] *= weight; 
                //End Butterworth filter 
            }   // end for freq bin
                            
            // now get time domain waveform back by inv fft               
            
            Tools::realft(V_forfft, -1, Nnew);            
            
            //TODO: Make this 160 more data-driven.  Shouldn't be hard-coded, but constants that can be changed by the user.
            if (antenna_theta > 160 or antenna_theta < 90) {
                cout << "Event outside of theta range.  Setting voltage to zero and moving to next event." << endl;
                for (int i = 0; i < sizeof(V_forfft) / sizeof(V_forfft[0]); i++) {
                  V_forfft[i] = 1;
                }                
                // continue;
            }                 
            Tools::SincInterpolation(Nnew, T_forfft, V_forfft, settings1->NFOUR / 2, T_forint, volts_forint);
            // Tools::SincInterpolation(Nnew, T_forfft, V_forfft, waveform_bin, T_forint, volts_forint);      
            
            //TODO: Restore time shift in time domain
            // T_forint += timeshift;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   T_forint[i] += timeshift;
            // }                 
            
            //Now write deconvolved voltage data to file.
            for (int n = 0; n < settings1->NFOUR / 2; n++)
            // for (int n = 0; n < waveform_bin; n++)
            {
                int elecChan = AraGeomTool::Instance()->getElecChanFromRFChan(i, settings1->DETECTOR_STATION);
                // not pure noise mode (we need signal)
                usefulAtriEvPtrOut->fVolts[elecChan].push_back(settings1->ARBITRARY_EVENT_ATTENUATION *volts_forint[n] *2. / (double)(Nnew));  // 2/N for inverse FFT normalization factor
                usefulAtriEvPtrOut->fTimes[elecChan].push_back(T_forint[n]);                
            }
            usefulAtriEvPtrOut->stationId = settings1->DETECTOR_STATION;
            
        } //channel loop
        usefulAtriEvPtrOut->eventNumber = usefulAtriEvPtr->eventNumber;
        usefulAtriEvPtrOut->unixTime = usefulAtriEvPtr->unixTime;
        cout << "usefulAtriEvPtr->unixTime = " << usefulAtriEvPtr->unixTime << endl;
        cout << "usefulAtriEvPtrOut->unixTime = " << usefulAtriEvPtrOut->unixTime << endl;
        //Assign timestamp values to help identify calpulsers
        usefulAtriEvPtrOut->timeStamp = usefulAtriEvPtr->timeStamp;
        // Assign triggerInfo values to identify RF and software triggers.
        for (int bit = 0; bit < 4; bit++) {
            usefulAtriEvPtrOut->triggerInfo[bit] = usefulAtriEvPtr->triggerInfo[bit];
        }
        fpOut->cd();
        outTree->Fill();
        usefulAtriEvPtrOut->fVolts.clear();
        usefulAtriEvPtrOut->fTimes.clear();
    } //event loop

    fpOut->Write();
    fpOut->Close();
    fp->Close();
    fp2->Close();
    
} //main    