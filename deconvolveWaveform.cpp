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
        std::cout << "e.g.\n" << argv[0] << " 2 6 AraOut.root output/\n";
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
    IceModel *iceModel = new IceModel(0 + 1*10, 0, 0);
    Settings *settings = new Settings();
    // settings->ReadEvtFile(fp);
    settings->ReadFile(setupfile);
    Detector *detector = new Detector(settings, iceModel, setupfile);
    Report *report = new Report(detector, settings);
    // configure settings    
    // settings->Z_THIS_TOLERANCE = 1; // higher tolerance
    // settings->Z_TOLERANCE = 0.05;
    // settings->NOFZ=1; // make sure n(z) is turned on
    // settings->RAY_TRACE_ICE_MODEL_PARAMS = 0; // set the ice model as user requested
    
    cout << "Settings->TIMESTEP = " << settings->TIMESTEP << endl;
    
    printf("------------------\n");
    printf("Make Output Files\n");
    printf("------------------\n");

    char outfile_name[400];
    sprintf(outfile_name, "%s/deconvolvedWaveforms_run_%s.root", argv[6], argv[3]);    
    
    std::cout<<"Output name is "<<outfile_name<<std::endl;
    TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
    if(!fpOut){ std::cerr<<"Cannot open output file "<<fpOut<<std::endl; return -1; }    
    
    TTree *outTree = new TTree("deconvolvedWaveform", "deconvolvedWaveform");
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
        
        if (event != 0){
               continue;
        }
    
        std::cout<<"Looking at event number "<<event<<std::endl;
        
        // UsefulAtriStationEvent *usefulAtriEvPtrOut = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
        // outTree->Branch("UsefulAtriStationEvent", &usefulAtriEvPtrOut);        
        
        //Loop over channels  MAP THESE TO ARASIM CHANNELS YOU DOLT - Love, Justin
        for(int i=0; i<16; i++){
            
            if (i != 0){
               continue;
            }
            TGraph *gr = usefulAtriEvPtr->getGraphFromRFChan(i);
            int waveform_bin = gr->GetN();
            
            double heff_lastbin;
            double freq_lastbin;
            double time[waveform_bin];
            double voltage[waveform_bin];
            double volts_forint[settings->NFOUR / 2];
            double T_forint[settings->NFOUR / 2];
            
            for(int k=0; k<waveform_bin; k++){
                time[k] = gr->GetX()[k];
                voltage[k] = gr->GetY()[k];
            }            
            delete gr;
            

    
            double freq_tmp, heff, antenna_theta, antenna_phi;  // values needed for apply antenna gain factor and prepare fft, trigger
            
            antenna_theta = reco_arrivalThetas[i]*180/PI;
            antenna_phi = reco_arrivalPhis[i]*180/PI;
            
            cout << "antenna_theta = " << antenna_theta << endl;
            cout << "antenna_phi = " << antenna_phi << endl;                
            
            // cout << "theta = " << antenna_theta << endl;
            // cout << "phi = " << antenna_phi << endl;            
            
            //Calculate polarization vector that inverts the polarization factor (makes dot products equal to one)
            // polarization=np.array([-np.sin(phi),np.cos(phi),-1/np.sin(theta)]) via Jorge's pyrex application
            // double newPol_vectorX = -sin(antenna_phi);
            // double newPol_vectorY = cos(antenna_phi);
            // double newPol_vectorZ = -1/sin(antenna_theta);
            
            //Testing using the actual polarization vector
            double psi = 45*PI/180;
            double newPol_vectorX = -cos(psi)*cos(antenna_theta)*cos(antenna_phi) + sin(psi)*sin(antenna_phi);
            double newPol_vectorY = -cos(psi)*cos(antenna_theta)*sin(antenna_phi) - sin(psi)*cos(antenna_phi);
            double newPol_vectorZ = cos(psi)*sin(antenna_theta);            

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
            
            // cout << "time[0] = " << time[0] << endl;
            // cout << "time[1] = " << time[1] << endl;
            
            double dT_forfft = time[1] - time[0];
            
            cout << "dT_forfft = " << dT_forfft << endl;
        
            int Ntmp = settings->TIMESTEP *1.e9 / dT_forfft;
            
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
            Nnew = Nnew * settings->NFOUR / 2;
            cout << "Ntmp = " << Ntmp << endl;
            cout << "Nnew = " << Nnew << endl;                

            //Stealing antenna and electronic response steps from AraSim, but applying the inverse functions instead.
            
            double V_forfft[Nnew];
            double T_forfft[Nnew];
            
            
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
            cout << "Before FFT V_forfft[230] = " << V_forfft[230] << endl;
            cout << "Before FFT V_forfft[231] = " << V_forfft[231] << endl;            
            
            // get spectrum with zero padded WF
            Tools::realft(V_forfft, 1, Nnew); 
            
            cout << "After FFT V_forfft[230] = " << V_forfft[230] << endl;
            cout << "After FFT V_forfft[231] = " << V_forfft[231] << endl;                  
            
            // cout << "Nnew = " << Nnew << endl;
            
            dF_Nnew = 1. / ((double)(Nnew) *(dT_forfft) *1.e-9);    // in Hz
            
            // cout << "dF_Nnew = " << dF_Nnew << endl;

            freq_tmp = dF_Nnew *((double) Nnew / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

            freq_lastbin = freq_tmp;
            
            // cout << "freq_tmp = " << freq_tmp << endl;
            

            for (int n = 0; n < settings->NFOUR / 2; n++)
            {
                // if (n != 230){
                //    continue;
                // }                

                freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq
                
                // cout << "freq_tmp = " << freq_tmp << endl;

                heff = report->GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                        antenna_theta, antenna_phi, pol_ant),
                    freq_tmp, nice);
                // cout << "heff = " << heff << endl;

                // cout << "n = " << n << endl;
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
                    cout << "heff = " << heff << endl;
                    cout << "heff_lastbin = " << heff_lastbin << endl;
                    cout << "Pol_factor = " << Pol_factor << endl;

                }

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
            }   // end for freq bin
            // now get time domain waveform back by inv fft
            // cout << "Nnew = " << Nnew << endl;            
            cout << "Before InvFFT V_forfft[230] = " << V_forfft[230] << endl;
            cout << "Before InvFFT V_forfft[231] = " << V_forfft[231] << endl;
            Tools::realft(V_forfft, -1, Nnew);
            cout << "V_forfft[230] = " << V_forfft[230] << endl;              
            cout << "V_forfft[231] = " << V_forfft[231] << endl;           
            
            Tools::SincInterpolation(Nnew, T_forfft, V_forfft, settings->NFOUR / 2, T_forint, volts_forint);
            // usefulAtriEvPtrOut->fVolts[i].push_back(V_forfft);
            // usefulAtriEvPtrOut->fTimes[i].push_back(time);

            
            
            //Now write deconvolved voltage data to file.
            //    HOW DO??????
            //
            
            for (int n = 0; n < settings->NFOUR / 2; n++)
            {

                // not pure noise mode (we need signal)
                usefulAtriEvPtrOut->fVolts[i].push_back(settings->ARBITRARY_EVENT_ATTENUATION *volts_forint[n] *2. / (double)(Nnew));  // 2/N for inverse FFT normalization factor
                usefulAtriEvPtrOut->fTimes[i].push_back(T_forint[n]);
                // cout << "volts_forint[" << n << "] = " << volts_forint[n] << endl;
                // cout << "T_forint[" << n << "] = " << T_forint[n] << endl;
            }
            
        } //channel loop
        fpOut->cd();
        outTree->Fill();
    } //event loop

    fpOut->Write();
    fpOut->Close();
    fp->Close();
    fp2->Close();
    
} //main    