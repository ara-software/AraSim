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
        // UsefulAtriStationEvent *usefulAtriEvPtr;
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
    // IceModel *iceModel = new IceModel(0 + 1*10, 0, 0);
    Settings *settings1 = new Settings();
    // settings->ReadEvtFile(fp);
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
        
        // if (event != 19){
        //        continue;
        // }
    
        std::cout<<"Looking at event number "<<event<<std::endl;
        
        if (dataLike) {
            cout << "Triggering datalike condition." << endl;
            delete usefulAtriEvPtr;         
            usefulAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
            //Testing setting the outgoing pointer using the rawEvPtr to populate the calpulser/softtrigger parameters.
            // delete usefulAtriEvPtrOut;
            // usefulAtriEvPtrOut = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
//             if (rawAtriEvPtr->isCalpulserEvent()) {
//                 cout << "Is calpulser" << endl;
//                 usefulAtriEvPtrOut->isTrigType(1);
            
//             }
//             else if (rawAtriEvPtr->isSoftwareTrigger()) {
//                 cout << "Is soft trigger" << endl;
//                 usefulAtriEvPtrOut->isTrigType(2);
//             }
//             else {
//                 cout << "Is rf event" << endl;
//                 usefulAtriEvPtrOut->isTrigType(0);
//             }                    
        }

        int vertexRecoElectToRFChan[] = {14,2,6,10,12,0,4,8,15,3,7,11,13,1,5,9};
        
        

        for(int i=0; i<16; i++){
            // settings1->NFOUR = 4096;
            
            // if (i != 0) {
            //     continue;
            // }
            // cout << "aaaaa" << endl;
            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);  //This is where the code breaks for real data.
            
            //Save initial and final time for truncating the padded arrays before output.
            double timeStart = gr->GetX()[0];
            double timeEnd = gr->GetX()[gr->GetN()-1];
            
            //Interpolate graph to 0.5 ns resolution
            gr = FFTtools::getInterpolatedGraph(gr,0.5);
            
            //Pad waveform to a factor of two. - JCF 9/27/2023
            
            // cout << "gr->GetN() = "<< gr->GetN() << endl;
            if (gr->GetN() < settings1->NFOUR/2) {
                gr = FFTtools::padWaveToLength(gr, settings1->NFOUR/2);
            }
            // Padding 
            // cout << "bbbbbb" << endl;
            int waveform_bin = gr->GetN();
            // cout << "ccccccc" << endl;
            // cout << "gr->GetN() = "<< gr->GetN() << endl;
            //Setting NFOUR to 4096 for testing on SpiceCore events
            // settings1->NFOUR = 8192;
            
            
            double heff_lastbin;
            double freq_lastbin;
            double time[waveform_bin];
            double voltage[waveform_bin];
            // cout << "ddddddd" << endl;
            double volts_forint[settings1->NFOUR / 2];
            double T_forint[settings1->NFOUR / 2];
            // cout << "waveform_bin = "<< waveform_bin << endl;
            // cout << "settings1->TIMESTP = " << settings1->TIMESTEP << endl;
            // double volts_forint[ int(waveform_bin /(settings1->TIMESTEP*1e9))];
            // cout << "eeeeeee" << endl;
            // double T_forint[ int(waveform_bin /(settings1->TIMESTEP*1e9))]; 
            // cout << "fffffff" << endl;
            
            //TODO: This init_T isn't dynamic to the imported data.  Should make have it defined based on the input waveform.
            double init_T = settings1->TIMESTEP *-1.e9 *((double) settings1->NFOUR / 4);    // locate zero time at the middle and give random time shift
            // double init_T = ((double) waveform_bin / 2); 
            // cout << "ggggggg" << endl;
           


            // for (int n = 0; n < settings1->NFOUR / 2; n++)
            // {
            //     T_forint[n] = init_T + (double) n *settings1->TIMESTEP *1.e9;   // in ns
            // }            
            
            for(int k=0; k<waveform_bin; k++){
                time[k] = gr->GetX()[k];
                voltage[k] = gr->GetY()[k];
            }            
            delete gr;
            
            // cout << "init_T = " << init_T << endl;
            // cout << "Changing init_T based on input waveform." << endl;
            // init_T = time[0];
            // cout << "init_T = " << init_T << endl;
            
            // for (int m = 0; m < settings1->NFOUR / 2; m++)
            for (int m = 0; m < 2048; m++)
            {
                T_forint[m] = -512 + m*0.5;   // in ns
            }
            
            // std::cout << std::endl; 
            // cout << "T_forint (after definition) = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
            // std::cout << std::endl; 
            
            //Importing the cutoff time between spicecore peaks
            double cutoffTimeChannel;
            // if (dataLike) {
            //     cutoffTime = vertexReco->cutoffTime[i]; 
            // } else {
            //     cutoffTime = time[waveform_bin-1];
            // }
            if (!cutoffTime) {
                cutoffTimeChannel = time[waveform_bin-1];
            } else {
                cutoffTimeChannel = cutoffTime[i];
            }            
            //Apply cutoff time to waveform by creating mask array
            //First attempt is just set values past the timeCutoff equal to zero
            // voltage[time > cutoffTimeChannel] = 0;
            // time[time >cutoffTimeChannel] = 0;
            // for(int k=0; k<waveform_bin; k++){
            //     if(time[k] > cutoffTimeChannel){
            //         // time[k] = 0;
            //         voltage[k] = 0;                
            //     }
            // }               

            // std::cout << std::endl; 
            // cout << "T_forint (aaaaaaaaa) = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
            // std::cout << std::endl; 
            
            // cout << "CutoffTime = " << cutoffTimeChannel << endl;
            // cout << "Initial waveform v(t) length = " << waveform_bin << endl;
            // cout << "Initial waveform v(t) = " << endl;
            // for (int i = 0; i < sizeof(voltage) / sizeof(voltage[0]); i++) {
            //   std::cout << voltage[i] << ", ";
            // }            
            // std::cout << std::endl; 
            // cout << "T_forint (bbbbbbbbbb) = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
            // std::cout << std::endl; 
            //Add step that centers the waveform about zero in time, for purposes of the fourier transform.  Then save this shift and reapply it to restore the time-domain information after the InvFFT.
//             cout << "*****************************************************************" << endl;
//             cout << "time[0] = " << time[0] << endl;
//             cout << "time[(waveform_bin-1)/2] = " << time[(waveform_bin-1)/2] << endl;
//             cout << "time[waveform_bin/2] = " << time[waveform_bin/2] << endl;
//             cout << "time[waveform_bin-1] = " << time[waveform_bin-1] << endl;
//             cout << "*****************************************************************" << endl;
            
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
            
            // std::cout << std::endl; 
            // cout << "T_forint (cccccccccc) = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
            // std::cout << std::endl; 
            cout << "antenna_theta = " << antenna_theta << endl;
            cout << "antenna_phi = " << antenna_phi << endl;
            
            //Adding condition to filter out downgoing events and event near the theta poles.
            //I should write this on a per event basis rather than per channel.
            // //TODO: Make this 160 more data-driven.  Shouldn't be hard-coded, but constants that can be changed by the user.
            // if (antenna_theta > 160 or antenna_theta < 90) {
            //     cout << "Event outside of theta range.  Moving to next event." << endl;
            //     // for (int n = 0; n < settings1->NFOUR / 2; n++)
            //     // {
            //     //     // Push waveform of all zeros into outputfile.
            //     //     int elecChan = AraGeomTool::Instance()->getElecChanFromRFChan(i, settings1->DETECTOR_STATION);
            //     //     // not pure noise mode (we need signal)
            //     //     usefulAtriEvPtrOut->fVolts[elecChan].push_back(settings1->ARBITRARY_EVENT_ATTENUATION *volts_forint[n] *2. / (double)(Nnew));  // 2/N for inverse FFT normalization factor
            //     //     usefulAtriEvPtrOut->fTimes[elecChan].push_back(T_forint[n]);                
            //     // }                
            //     continue;
            // }            
            
            // cout << "theta = " << antenna_theta << endl;
            // cout << "phi = " << antenna_phi << endl;            
            
            //Calculate polarization vector that inverts the polarization factor (makes dot products equal to one)
            // polarization=np.array([-np.sin(phi),np.cos(phi),-1/np.sin(theta)]) via Jorge's pyrex application
            double newPol_vectorX = -sin(antenna_phi*PI/180);
            double newPol_vectorY = cos(antenna_phi*PI/180);
            double newPol_vectorZ = -1/sin(antenna_theta*PI/180);
            
            // std::cout << std::endl; 
            // cout << "T_forint (ddddddddd) = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
            // std::cout << std::endl;             
            //Testing using the actual polarization vector
            // double psi = argv[3]*PI/180;
            // double newPol_vectorX = -cos(psi)*cos(antenna_theta*PI/180)*cos(antenna_phi*PI/180) + sin(psi)*sin(antenna_phi*PI/180);
            // double newPol_vectorY = -cos(psi)*cos(antenna_theta*PI/180)*sin(antenna_phi*PI/180) - sin(psi)*cos(antenna_phi*PI/180);
            // double newPol_vectorZ = cos(psi)*sin(antenna_theta*PI/180);            

            Vector Pol_vector = Vector(newPol_vectorX, newPol_vectorY, newPol_vectorZ);

            double dF_Nnew;

            // int NFOUR = 1024;
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
            
            // cout << "dT_forfft = " << dT_forfft << endl;
        
            int Ntmp = settings1->TIMESTEP *1.e9 / dT_forfft;
            
            int Nnew = 1;
            // cout << "Ntmp = " << Ntmp << endl;
            // cout << "Nnew = " << Nnew << endl;            
            while (Ntmp > 1)
            {
                Ntmp = Ntmp / 2;
                Nnew = Nnew *2;
                // cout << "Ntmp = " << Ntmp << endl;
                // cout << "Nnew = " << Nnew << endl;                
            }
            Nnew = Nnew * settings1->NFOUR / 2;
            // Nnew = Nnew * waveform_bin;
            // cout << "Ntmp = " << Ntmp << endl;
            // cout << "Nnew = " << Nnew << endl;                

            // std::cout << std::endl; 
            // cout << "T_forint (eeeeeeee) = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
            // std::cout << std::endl;             
            
            //Stealing antenna and electronic response steps from AraSim, but applying the inverse functions instead.
            
            double V_forfft[Nnew];
            double T_forfft[Nnew];
            
            // cout << "waveform_bin = " << waveform_bin << endl;
            
            for (int n = 0; n < Nnew; n++)
            {

                // make Tarray, Earray located at the center of Nnew array

                T_forfft[n] = time[waveform_bin / 2] - (dT_forfft *(double)(Nnew / 2 - n));
                // T_forfft[n] = time[waveform_bin / 2] - (dT_forfft *(double)(Nnew / 2 - n)) - timeshift;  //Applying time shift to center array in time domain.

                if ((n >= Nnew / 2 - waveform_bin / 2) &&
                    (n < Nnew / 2 + waveform_bin / 2))
                {
                    V_forfft[n] = voltage[n - (Nnew / 2 - waveform_bin / 2)];
                }
                else
                    V_forfft[n] = 0.;

            }
//             std::cout << std::endl; 
//             cout << "T_forint (ffffffffff) = " << endl;
//             for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
//               std::cout << T_forint[i] << ", ";
//             }
//             std::cout << std::endl;                  
//             // cout << "Before FFT V_forfft[230] = " << V_forfft[230] << endl;
//             // cout << "Before FFT V_forfft[231] = " << V_forfft[231] << endl;
//             // cout << "gain_ch_no = " << gain_ch_no << endl;
//             cout << "Before FFT V_forfft = " << endl;
//             for (int i = 0; i < sizeof(V_forfft) / sizeof(V_forfft[0]); i++) {
//               std::cout << V_forfft[i] << ", ";
//             }
//             std::cout << std::endl;    
            
//             cout << "Before FFT T_forfft = " << endl;
//             for (int i = 0; i < sizeof(T_forfft) / sizeof(T_forfft[0]); i++) {
//               std::cout << T_forfft[i] << ", ";
//             }
//             std::cout << std::endl;
            
//             cout << "Before FFT time = " << endl;
//             for (int i = 0; i < sizeof(time) / sizeof(time[0]); i++) {
//               std::cout << time[i] << ", ";
//             }
//             std::cout << std::endl;             
            
            // get spectrum with zero padded WF
            Tools::realft(V_forfft, 1, Nnew); 
            // cout << "After FFT V_forfft = " << endl;
            // for (int i = 0; i < sizeof(V_forfft) / sizeof(V_forfft[0]); i++) {
            //   std::cout << V_forfft[i] << ", ";
            // }
            // std::cout << std::endl;            
            
            // cout << "After FFT V_forfft[230] = " << V_forfft[230] << endl;
            // cout << "After FFT V_forfft[231] = " << V_forfft[231] << endl;                  
            
            // cout << "Nnew = " << Nnew << endl;
            
            dF_Nnew = 1. / ((double)(Nnew) *(dT_forfft) *1.e-9);    // in Hz
            
            // cout << "dT_forfft = " << dT_forfft << endl;
            // cout << "Nnew = " << Nnew << endl;
            // cout << "dF_Nnew = " << dF_Nnew << endl;
            
            // cout << "dF_Nnew = " << dF_Nnew << endl;

            freq_tmp = dF_Nnew *((double) Nnew / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

            freq_lastbin = freq_tmp;
            
            // cout << "freq_tmp = " << freq_tmp << endl;
            
            // std::cout << std::endl; 
            // cout << "T_forint (ggggggggg) = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
            // std::cout << std::endl;      
            
            
            
            for (int n = 0; n < Nnew / 2; n++)
            // for (int n = 0; n < settings1->NFOUR / 2; n++)            
            {
                // if (n != 230){
                //    continue;
                // }                
                // cout << "**************************************************************" << endl;
                // cout << "n = " << n << endl;
                // cout << "aaa T_forint[0] = " << T_forint[0] << endl;
                freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq
                // cout << "bbb T_forint[0] = " << T_forint[0] << endl;
                heff_lastbin = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                        antenna_theta, antenna_phi, pol_ant),
                    freq_tmp, nice);                
                // cout << "ccc T_forint[0] = " << T_forint[0] << endl;
                // cout << "freq_tmp*1e-6 = " << freq_tmp*1e-6 << endl;
                // cout << "gain_ch_no = " << gain_ch_no << endl;
                // cout << "detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta, antenna_phi, pol_ant) = " << detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta, antenna_phi, pol_ant) << endl;
                // cout << "detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, pol_ant) = " << detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, pol_ant) << endl;                        
                // cout << "detector->GetElectGain_1D_OutZero( freq_tmp*1.e-6, gain_ch_no ) = " << detector->GetElectGain_1D_OutZero( freq_tmp*1.e-6, gain_ch_no ) << endl;
                // cout << "detector->GetElectPhase_1D(freq_tmp*1.e-6, gain_ch_no) = " << detector->GetElectPhase_1D(freq_tmp*1.e-6, gain_ch_no) << endl;                

                heff = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                        antenna_theta, antenna_phi, pol_ant),
                    freq_tmp, nice);
                // cout << "ddd T_forint[0] = " << T_forint[0] << endl;
                // cout << "heff = " << heff << endl;

                // cout << "n = " << n << endl;
                //
                // invert entire elect chain gain, phase
                //
                if (n > 0)
                {                                             
                    // cout << "Before invert elect (Real): "<< V_forfft[2 *n] << endl;
                    // cout << "Before invert elect (Imag): "<< V_forfft[2 *n + 1] << endl;                    
                    report->InvertElect_Tdomain(freq_tmp *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no);
                    // cout << "After invert elect (Real): "<< V_forfft[2 *n] << endl;
                    // cout << "After invert elect (Imag): "<< V_forfft[2 *n + 1] << endl;
                }
                else
                {
                    report->InvertElect_Tdomain_FirstTwo(freq_tmp *1.e-6, freq_lastbin *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no);
                }
                // cout << "eee T_forint[0] = " << T_forint[0] << endl;
                if (n > 0)
                {
                
                    // cout << "Before invert Ant (Real): "<< V_forfft[2 *n] << endl;
                    // cout << "Before invert Ant (Imag): "<< V_forfft[2 *n + 1] << endl;
                    report->InvertAntFactors_Tdomain(detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, pol_ant),
                       heff, Pol_vector, pol_ant, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi);
                    // cout << "After invert Ant (Real): "<< V_forfft[2 *n] << endl;
                    // cout << "After invert Ant (Imag): "<< V_forfft[2 *n + 1] << endl;                    
                }
                else
                {
                    report->InvertAntFactors_Tdomain_FirstTwo(heff, heff_lastbin, Pol_vector, pol_ant, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi);
                    // cout << "heff = " << heff << endl;
                    // cout << "heff_lastbin = " << heff_lastbin << endl;
                    // cout << "Pol_factor = " << Pol_factor << endl;

                }
                // cout << "fff T_forint[0] = " << T_forint[0] << endl;
                
                //Quick and dirty hack to filter out frequencies above 850 MHz and below 100 MHz.
                // cout << "******************************************" << endl;
                // cout << "Before Hard band-pass V_forfft = " << endl;
                // cout << "V_forfft[2 *n] = "<< V_forfft[2 *n] << endl;
                // cout << "V_forfft[2 *n + 1] = "<< V_forfft[2 *n + 1] << endl;                      
                // std::cout << std::endl;                
                if (freq_tmp > 850.*1.e6 or freq_tmp < 100.*1.e6) {
                    V_forfft[2*n] = 0;
                    V_forfft[2*n+1] = 0;
                }
                // cout << "ggg T_forint[0] = " << T_forint[0] << endl;
                // cout << "Before Butterworth V_forfft = " << endl;
                // cout << "V_forfft[2 *n] = "<< V_forfft[2 *n] << endl;
                // cout << "V_forfft[2 *n + 1] = "<< V_forfft[2 *n + 1] << endl;                     
                //Apply homemade butterworth filter of the fourth order
                // double freqMin = 150*1e6;
                // double freqMax = 300*1e6;
                
                //Trying user inputted butterworth filter
                double freqMin = atof(argv[7])*1e6;
                double freqMax = atof(argv[8])*1e6;
                cout << "freqMin = " << freqMin << endl;
                cout << "freqMax = " << freqMax << endl;
  
                
                double weight = 1;  //TODO:  Make this dynamic for simulations and real data. - JCF 10/4/2023
                int order = 8;
                weight /= sqrt(1 + TMath::Power(freqMin/freq_tmp, 4*order));
                weight /= sqrt(1 + TMath::Power(freq_tmp/freqMax, 4*order));
                // cout << "hhh T_forint[0] = " << T_forint[0] << endl;
//                 if (freq_tmp < freqMin) {
                        // weight /= sqrt(1 + TMath::Power(freqMin/freq_tmp, 4*order));
//                 }
//                 if (freq_tmp > freqMax) {
//                         weight /= sqrt(1 + TMath::Power(freq_tmp/freqMax, 4*order));
//                 }
                V_forfft[2*n] *= weight;
                V_forfft[2*n+1] *= weight; 
                // cout << "iii T_forint[0] = " << T_forint[0] << endl;
                //End Butterworth filter
                // cout << "After Butterworth V_forfft = " << endl;
                // cout << "V_forfft[2 *n] = "<< V_forfft[2 *n] << endl;
                // cout << "V_forfft[2 *n + 1] = "<< V_forfft[2 *n + 1] << endl;     
                // //
                // // invert entire elect chain gain, phase
                // //
                // if (n > 0)
                // {                                             
                //     // cout << "Before invert elect (Real): "<< V_forfft[2 *n] << endl;
                //     // cout << "Before invert elect (Imag): "<< V_forfft[2 *n + 1] << endl;                    
                //     report->InvertElect_Tdomain(freq_tmp *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no);
                //     // cout << "After invert elect (Real): "<< V_forfft[2 *n] << endl;
                //     // cout << "After invert elect (Imag): "<< V_forfft[2 *n + 1] << endl;                    
                // }
                // else
                // {
                //     report->InvertElect_Tdomain_FirstTwo(freq_tmp *1.e-6, freq_lastbin *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no);
                // }
            }   // end for freq bin
            
            // std::cout << std::endl; 
            // cout << "T_forint (hhhhhhhhh) = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
            // std::cout << std::endl;                  
            // now get time domain waveform back by inv fft
            // cout << "Nnew = " << Nnew << endl;            
            // cout << "Before InvFFT V_forfft[230] = " << V_forfft[230] << endl;
            // cout << "Before InvFFT V_forfft[231] = " << V_forfft[231] << endl;
            // cout << "Before InvFFT V_forfft = " << endl;
            // for (int i = 0; i < sizeof(V_forfft) / sizeof(V_forfft[0]); i++) {
            //   std::cout << V_forfft[i] << ", ";
            // }
            // std::cout << std::endl;
               
            
            Tools::realft(V_forfft, -1, Nnew);
            // cout << "V_forfft[230] = " << V_forfft[230] << endl;              
            // cout << "V_forfft[231] = " << V_forfft[231] << endl;
            // cout << "After InvFFT V_forfft = " << endl;
            // for (int i = 0; i < sizeof(V_forfft) / sizeof(V_forfft[0]); i++) {
            //   std::cout << V_forfft[i] << ", ";
            // }

            // std::cout << std::endl;            
            // cout << "Nnew = " << Nnew << endl;
            // cout << "sizeof(T_forfft) = " << sizeof(T_forfft) << endl;
            // cout << "sizeof(V_forfft) = " << sizeof(V_forfft) << endl;
            // cout << "settings1->NFOUR / 2 = " << settings1->NFOUR / 2 << endl;
            // cout << "sizeof(T_forint) = " << sizeof(T_forint) << endl;
            // cout << "sizeof(volts_forint) = " << sizeof(volts_forint) << endl;
            // cout << "T_forfft = " << endl;
            // for (int i = 0; i < sizeof(T_forfft) / sizeof(T_forfft[0]); i++) {
            //   std::cout << T_forfft[i] << ", ";
            // }
            // std::cout << std::endl; 
            // cout << "T_forint = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
//             std::cout << std::endl;
//             cout << "volts_forint = " << endl;
//             for (int i = 0; i < sizeof(volts_forint) / sizeof(volts_forint[0]); i++) {
//               std::cout << volts_forint[i] << ", ";
//             }
//             std::cout << std::endl; 
            
            
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
            // cout << "T_forint (after interpolation) = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
            // std::cout << std::endl;
            
            //Restore time shift in time domain
            // T_forint += timeshift;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   T_forint[i] += timeshift;
            // }            
            
            // cout << "volts_forint (after interpolation) = " << endl;
            // for (int i = 0; i < sizeof(volts_forint) / sizeof(volts_forint[0]); i++) {
            //   std::cout << volts_forint[i] << ", ";
            // }
            // std::cout << std::endl;
            // cout << "volts_forint (normalized) = " << endl;
            // for (int i = 0; i < sizeof(volts_forint) / sizeof(volts_forint[0]); i++) {
            //   std::cout << volts_forint[i]*2/Nnew << ", ";
            // }
            // std::cout << std::endl;            
            // usefulAtriEvPtrOut->fVolts[i].push_back(V_forfft);
            // usefulAtriEvPtrOut->fTimes[i].push_back(time);

            
            
            //Now write deconvolved voltage data to file.
            //    HOW DO??????
            //
            // int fVoltsMap[16] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27};
            // int fVoltsMap[16] = {8, 0, 24, 16, 9, 1, 25, 17, 10, 2, 26, 18, 11, 3, 27, 19};  //Troubleshooting volts map - JCF 7/11/2023
            // int fVoltsMap[16] = {8, 0, 16, 24, 1, 9, 17, 25, 2, 10, 18, 26, 3, 11, 19, 27};
            for (int n = 0; n < settings1->NFOUR / 2; n++)
            // for (int n = 0; n < waveform_bin; n++)
            {
                // int fVoltsMap[16] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27};
                int elecChan = AraGeomTool::Instance()->getElecChanFromRFChan(i, settings1->DETECTOR_STATION);
                // not pure noise mode (we need signal)
                usefulAtriEvPtrOut->fVolts[elecChan].push_back(settings1->ARBITRARY_EVENT_ATTENUATION *volts_forint[n] *2. / (double)(Nnew));  // 2/N for inverse FFT normalization factor
                usefulAtriEvPtrOut->fTimes[elecChan].push_back(T_forint[n]);                
                // // not pure noise mode (we need signal)
                // usefulAtriEvPtrOut->fVolts[fVoltsMap[i]].push_back(settings1->ARBITRARY_EVENT_ATTENUATION *volts_forint[n] *2. / (double)(Nnew));  // 2/N for inverse FFT normalization factor
                // usefulAtriEvPtrOut->fTimes[fVoltsMap[i]].push_back(T_forint[n]);
            }
            // cout << "sizeof(fVolts) = " << sizeof(usefulAtriEvPtrOut->fVolts) << endl;
            usefulAtriEvPtrOut->stationId = settings1->DETECTOR_STATION;
            
        } //channel loop
        usefulAtriEvPtrOut->eventNumber = usefulAtriEvPtr->eventNumber;
        usefulAtriEvPtrOut->unixTime = usefulAtriEvPtr->unixTime;
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