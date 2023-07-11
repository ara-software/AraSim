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

    // data like
    TTree *eventTree = (TTree*) fp->Get("eventTree");
    printf("Event tree opened!\n");
    if(!eventTree) { std::cerr << "Can't find eventTree\n"; return -1; }
    RawAtriStationEvent *rawAtriEvPtr=0;
    // eventTree->SetBranchAddress("event",&rawAtriEvPtr);    
    eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
    double weight;
    eventTree->SetBranchAddress("weight", &weight);
    Long64_t numEntries=eventTree->GetEntries();
    
    //Import vertex reco file
    printf("Opening reco file...\n");
    TFile *fp2 = TFile::Open(argv[5]);
    if(!fp2) { std::cerr << "Can't open file\n"; return -1; }
    printf("Reco File opened!\n");
    TTree *vertexReco = (TTree*) fp2->Get("vertexReco");
    double reco_arrivalThetas[16];
    double reco_arrivalPhis[16];
    
    // Testing using the true rf angles
    vertexReco->SetBranchAddress("true_arrivalThetas", reco_arrivalThetas);
    vertexReco->SetBranchAddress("true_arrivalPhis", reco_arrivalPhis);   
    // end testing
    // vertexReco->SetBranchAddress("reco_arrivalThetas", reco_arrivalThetas);
    // vertexReco->SetBranchAddress("reco_arrivalPhis", reco_arrivalPhis);
    
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
    // cout<<"number of stations : "<<detector->params.number_of_stations << endl;
    // cout<<"total number of antennas : "<<detector->params.number_of_antennas << endl;    
    Report *report = new Report(detector, settings1);
    // configure settings    
    // settings->Z_THIS_TOLERANCE = 1; // higher tolerance
    // settings->Z_TOLERANCE = 0.05;
    // settings->NOFZ=1; // make sure n(z) is turned on
    // settings->RAY_TRACE_ICE_MODEL_PARAMS = 0; // set the ice model as user requested
    
    cout << "Settings->TIMESTEP = " << settings1->TIMESTEP << endl;
    
    printf("------------------\n");
    printf("Make Output Files\n");
    printf("------------------\n");

    char outfile_name[400];
    sprintf(outfile_name, "%s/deconvolvedWaveforms_run_%s.root", argv[6], argv[3]);    
    
    std::cout<<"Output name is "<<outfile_name<<std::endl;
    TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
    if(!fpOut){ std::cerr<<"Cannot open output file "<<fpOut<<std::endl; return -1; }    
    
    // TTree *outTree = new TTree("deconvolvedWaveform", "deconvolvedWaveform");
     TTree *outTree = new TTree("eventTree", "eventTree");
    // UsefulAtriStationEvent *usefulAtriEvPtrOut = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
    outTree->Branch("UsefulAtriStationEvent", &usefulAtriEvPtrOut);
    
    //Need to grab lengths of voltage and time arrays from eventTree to initialize the branches in the outfile.
    Int_t fNumChannels; ///< The number of channels
    std::map< Int_t, std::vector <Double_t> > fTimesOut; ///< The times of samples
    std::map< Int_t, std::vector <Double_t> > fVoltsOut; ///< The voltages of samples    
    // eventTree->GetEntry(0);
    // TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(0);
    // int waveform_bin = gr->GetN();
    // outTree->Branch("deconvolvedVoltage", fVoltsOut, 
    
    
    //Loop over events
    for(Long64_t event=0;event<numEntries;event++) {
        fp->cd();
        eventTree->GetEntry(event);
        vertexReco->GetEntry(event);
        
        // if (event != 0){
        //        continue;
        // }
    
        std::cout<<"Looking at event number "<<event<<std::endl;
        
        // UsefulAtriStationEvent *usefulAtriEvPtrOut = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
        // outTree->Branch("UsefulAtriStationEvent", &usefulAtriEvPtrOut);        
        
        //TODO: Loop over channels  MAP THESE TO ARASIM CHANNELS YOU DOLT - Love, Justin
        //simLikeToDataLikeMask = [14,2,6,10,12,0,4,8,15,3,7,11,13,1,5,9] **use this for vertex reco angles
        int vertexRecoElectToRFChan[] = {14,2,6,10,12,0,4,8,15,3,7,11,13,1,5,9};
        for(int i=0; i<16; i++){
            
            // if (i != 5){
            //    continue;
            // }
            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);
            int waveform_bin = gr->GetN();
            
            double heff_lastbin;
            double freq_lastbin;
            double time[waveform_bin];
            double voltage[waveform_bin];
            double volts_forint[settings1->NFOUR / 2];
            double T_forint[settings1->NFOUR / 2];
            double init_T = settings1->TIMESTEP *-1.e9 *((double) settings1->NFOUR / 4);    // locate zero time at the middle and give random time shift

            for (int n = 0; n < settings1->NFOUR / 2; n++)
            {
                T_forint[n] = init_T + (double) n *settings1->TIMESTEP *1.e9;   // in ns
            }            
            
            for(int k=0; k<waveform_bin; k++){
                time[k] = gr->GetX()[k];
                voltage[k] = gr->GetY()[k];
            }            
            delete gr;
            cout << "Initial waveform v(t) = " << endl;
            for (int i = 0; i < sizeof(voltage) / sizeof(voltage[0]); i++) {
              std::cout << voltage[i] << ", ";
            }            

    
            double freq_tmp, heff, antenna_theta, antenna_phi;  // values needed for apply antenna gain factor and prepare fft, trigger
            
            antenna_theta = reco_arrivalThetas[vertexRecoElectToRFChan[i]]*180/PI;
            antenna_phi = reco_arrivalPhis[vertexRecoElectToRFChan[i]]*180/PI;
            

            cout << "antenna_theta = " << antenna_theta << endl;
            cout << "antenna_phi = " << antenna_phi << endl;                
            
            // cout << "theta = " << antenna_theta << endl;
            // cout << "phi = " << antenna_phi << endl;            
            
            //Calculate polarization vector that inverts the polarization factor (makes dot products equal to one)
            // polarization=np.array([-np.sin(phi),np.cos(phi),-1/np.sin(theta)]) via Jorge's pyrex application
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

            // int NFOUR = 1024;
            double nice = 1.79;


            int pol_ant;
            int gain_ch_no = i;
            // int gain_ch_no = vertexRecoElectToRFChan[i];  //Testing channel mapping for electronics gain.            
            double Pol_factor;
            
            if (i < 8) {
                pol_ant=0;
            } else {
                pol_ant=1;
            }
            
            // cout << "time[0] = " << time[0] << endl;
            // cout << "time[1] = " << time[1] << endl;
            
            double dT_forfft = time[1] - time[0];
            
            cout << "dT_forfft = " << dT_forfft << endl;
        
            int Ntmp = settings1->TIMESTEP *1.e9 / dT_forfft;
            
            int Nnew = 1;
            cout << "Ntmp = " << Ntmp << endl;
            cout << "Nnew = " << Nnew << endl;            
            while (Ntmp > 1)
            {
                Ntmp = Ntmp / 2;
                Nnew = Nnew *2;
                cout << "Ntmp = " << Ntmp << endl;
                cout << "Nnew = " << Nnew << endl;                
            }
            Nnew = Nnew * settings1->NFOUR / 2;
            cout << "Ntmp = " << Ntmp << endl;
            cout << "Nnew = " << Nnew << endl;                

            //Stealing antenna and electronic response steps from AraSim, but applying the inverse functions instead.
            
            double V_forfft[Nnew];
            double T_forfft[Nnew];
            
            cout << "waveform_bin = " << waveform_bin << endl;
//             waveform_bin = 100;
            
//             cout << "redefining waveform_bin = " << waveform_bin << endl;
            
            for (int n = 0; n < Nnew; n++)
            {

                // make Tarray, Earray located at the center of Nnew array

                T_forfft[n] = time[waveform_bin / 2] - (dT_forfft *(double)(Nnew / 2 - n));

                if ((n >= Nnew / 2 - waveform_bin / 2) &&
                    (n < Nnew / 2 + waveform_bin / 2))
                {
                    V_forfft[n] = voltage[n - (Nnew / 2 - waveform_bin / 2)];
                }
                else
                    V_forfft[n] = 0.;

            }
            // cout << "Before FFT V_forfft[230] = " << V_forfft[230] << endl;
            // cout << "Before FFT V_forfft[231] = " << V_forfft[231] << endl;
            cout << "gain_ch_no = " << gain_ch_no << endl;
            cout << "Before FFT V_forfft = " << endl;
            for (int i = 0; i < sizeof(V_forfft) / sizeof(V_forfft[0]); i++) {
              std::cout << V_forfft[i] << ", ";
            }
            std::cout << std::endl;            
            
            // get spectrum with zero padded WF
            Tools::realft(V_forfft, 1, Nnew); 
            cout << "After FFT V_forfft = " << endl;
            for (int i = 0; i < sizeof(V_forfft) / sizeof(V_forfft[0]); i++) {
              std::cout << V_forfft[i] << ", ";
            }
            std::cout << std::endl;            
            
            // cout << "After FFT V_forfft[230] = " << V_forfft[230] << endl;
            // cout << "After FFT V_forfft[231] = " << V_forfft[231] << endl;                  
            
            // cout << "Nnew = " << Nnew << endl;
            
            dF_Nnew = 1. / ((double)(Nnew) *(dT_forfft) *1.e-9);    // in Hz
            
            cout << "dT_forfft = " << dT_forfft << endl;
            cout << "Nnew = " << Nnew << endl;
            cout << "dF_Nnew = " << dF_Nnew << endl;
            
            // cout << "dF_Nnew = " << dF_Nnew << endl;

            freq_tmp = dF_Nnew *((double) Nnew / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

            freq_lastbin = freq_tmp;
            
            // cout << "freq_tmp = " << freq_tmp << endl;
            

            // for (int n = 0; n < settings1->NFOUR / 2; n++)
            
            for (int n = 0; n < Nnew / 2; n++)
            {
                // if (n != 230){
                //    continue;
                // }                
                // cout << "**************************************************************" << endl;
                freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq
                
                heff_lastbin = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                        antenna_theta, antenna_phi, pol_ant),
                    freq_tmp, nice);                
                
                // cout << "freq_tmp*1e-6 = " << freq_tmp*1e-6 << endl;
                // cout << "gain_ch_no = " << gain_ch_no << endl;
                // cout << "detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta, antenna_phi, pol_ant) = " << detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta, antenna_phi, pol_ant) << endl;
                // cout << "detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, pol_ant) = " << detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, pol_ant) << endl;                        
                // cout << "detector->GetElectGain_1D_OutZero( freq_tmp*1.e-6, gain_ch_no ) = " << detector->GetElectGain_1D_OutZero( freq_tmp*1.e-6, gain_ch_no ) << endl;
                // cout << "detector->GetElectPhase_1D(freq_tmp*1.e-6, gain_ch_no) = " << detector->GetElectPhase_1D(freq_tmp*1.e-6, gain_ch_no) << endl;                

                heff = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                        antenna_theta, antenna_phi, pol_ant),
                    freq_tmp, nice);
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
                
                //Quick and dirty hack to filter out frequencies above 850 MHz.
                if (freq_tmp > 825.*1.e6 or freq_tmp < 140.*1.e6) {
                    V_forfft[2*n] = 0;
                    V_forfft[2*n+1] = 0;
                }

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
            // now get time domain waveform back by inv fft
            // cout << "Nnew = " << Nnew << endl;            
            // cout << "Before InvFFT V_forfft[230] = " << V_forfft[230] << endl;
            // cout << "Before InvFFT V_forfft[231] = " << V_forfft[231] << endl;
            cout << "Before InvFFT V_forfft = " << endl;
            for (int i = 0; i < sizeof(V_forfft) / sizeof(V_forfft[0]); i++) {
              std::cout << V_forfft[i] << ", ";
            }

            std::cout << std::endl;
               
            
            Tools::realft(V_forfft, -1, Nnew);
            // cout << "V_forfft[230] = " << V_forfft[230] << endl;              
            // cout << "V_forfft[231] = " << V_forfft[231] << endl;
//             cout << "After InvFFT V_forfft = " << endl;
//             for (int i = 0; i < sizeof(V_forfft) / sizeof(V_forfft[0]); i++) {
//               std::cout << V_forfft[i] << ", ";
//             }

//             std::cout << std::endl;            
//             cout << "Nnew = " << Nnew << endl;
//             cout << "sizeof(T_forfft) = " << sizeof(T_forfft) << endl;
//             cout << "sizeof(V_forfft) = " << sizeof(V_forfft) << endl;
//             cout << "settings1->NFOUR / 2 = " << settings1->NFOUR / 2 << endl;
//             cout << "sizeof(T_forint) = " << sizeof(T_forint) << endl;
//             cout << "sizeof(volts_forint) = " << sizeof(volts_forint) << endl;
//             cout << "T_forfft = " << endl;
//             for (int i = 0; i < sizeof(T_forfft) / sizeof(T_forfft[0]); i++) {
//               std::cout << T_forfft[i] << ", ";
//             }
//             std::cout << std::endl; 
//             cout << "T_forint = " << endl;
//             for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
//               std::cout << T_forint[i] << ", ";
//             }
//             std::cout << std::endl;
//             cout << "volts_forint = " << endl;
//             for (int i = 0; i < sizeof(volts_forint) / sizeof(volts_forint[0]); i++) {
//               std::cout << volts_forint[i] << ", ";
//             }
//             std::cout << std::endl;              
            Tools::SincInterpolation(Nnew, T_forfft, V_forfft, settings1->NFOUR / 2, T_forint, volts_forint);
            // cout << "T_forint (after interpolation) = " << endl;
            // for (int i = 0; i < sizeof(T_forint) / sizeof(T_forint[0]); i++) {
            //   std::cout << T_forint[i] << ", ";
            // }
            // std::cout << std::endl;
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
            cout << "sizeof(fVolts) = " << sizeof(usefulAtriEvPtrOut->fVolts) << endl;
            usefulAtriEvPtrOut->stationId = settings1->DETECTOR_STATION;
            
        } //channel loop
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