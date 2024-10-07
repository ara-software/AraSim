//--------------------------------------------------
// class Antenna_Response
//-------------------------------------------------- 
//

#ifndef REPORT_H
#define REPORT_H

#include <vector>
#include "TGraph.h"

#include "Position.h"

#ifndef __CINT__
//Include output format to enable reading by analysis software AraRoot
#ifdef ARA_UTIL_EXISTS
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraGeomTool.h"
#endif
#endif

#ifdef ARA_UTIL_EXISTS
class UsefulIcrrStationEvent;
class UsefulAtriStationEvent;
class AraGeomTool;
#endif

class Birefringence;
class Detector;
class Antenna;
class Event;
class RaySolver;
class Signal;
class IceModel;
class Interaction;
class Settings;
class Vector;
class Trigger;

class FFTWComplex;

using namespace std;

class Surface_antenna_r {
    public:

        ClassDef(Surface_antenna_r,1);
};

class Antenna_r {
    public:

        // one dimention for number of solutions, (if there,) another dimention for array of information
        //

        int ray_sol_cnt;    // number of RaySolver solutions

        //vector <int> trg;    // if antenna recieved any signal or not. 0 : no signal,  1 : yes signal

        vector < vector <double> > view_ang;    //viewing angle
        vector < vector <double> > launch_ang;  //launch angle
        vector < vector <double> > rec_ang;     //receiving angle phi (in radians)
        vector < vector <double> > phi_rec;     // receiving phi angle,in antenna's coord system.
        vector < vector <double> > theta_rec;     // receiving theta angle, in antenna's coord system.
        vector < vector <double> > phi_launch;     // launch phi angle,in antenna's coord system.
        vector < vector <double> > theta_launch;     // launch theta angle, in antenna's coord system.    
        vector < vector <double> > reflect_ang; // surface reflection angle (if 100 : no reflection case)
        vector < vector <double> > Dist;        //Distance between posnu and antenna
        vector < vector <double> > L_att;        //Attenuation factor
        vector < vector <double> > arrival_time;        //time from posnu to antenna (t:0 at posnu)
        vector < vector <int> > reflection;     // non-reflected : 0,  reflected : 1
        vector < vector < Position > > Pol_vector;   // polarization vector at the antenna
        //vector <Position> n_H;  // normalized vector for H pol
        //vector <Position> n_V;  // normalized vector for V pol

        //! Save every ray steps between the vertex (source) and an antenna (target), unless DATA_SAVE_MODE is 2. 02-12-2021 -MK-
        //! These xz coordinates were calculated after we convert the earth coordinates to flat coordinates by the RaySolver::Earth_to_Flat_same_angle()
        vector < vector < vector < vector <double> > > > ray_step;

        // below freq domain simulation output
        vector < vector < vector <double> > > vmmhz;  // signal V/m/MHz for each freq bin
        //
        vector < vector < vector <double> > > Heff;  // effective height for each freq bin
        vector < vector <double> > Mag;  // magnification factor
        vector < vector <double> > Fresnel;  // Fresnel factor
        vector < vector <double> > Pol_factor;  // Polarization factor
        vector < vector <double> > Pol_factorH;  //Hpol Polarization factor
        vector < vector <double> > Pol_factorV;  //Vpol Polarization factor
        
        vector < vector < vector <double> > > Vm_zoom;  // E field before ant T-domain
        vector < vector < vector <double> > > Vm_zoom_T;  // E field before ant T-domain time
        //vector < vector <double> > Vm_wo_antfactor;  // before applying ApplyAntFactors
        //vector < vector <double> > VHz_antfactor;  // after applying ApplyAntFactors to vmmhz above ( 1/sqrt2 * 1/dt * 0.5 * heff * pol_factor )
        //vector < vector <double> > VHz_filter;  // after applying ApplyAntFactors above and then apply filter gain from detector->GetFilterGain
        //
        int skip_bins[2]; // for two ray sols

        int Nnew[2]; // new number of bins for V_fotfft array

        vector < vector < vector <double> > > Vfft;  // signal V preparing for FFT
        vector < vector < vector <double> > > Vfft_noise;  // noise V preparing for FFT


        // below time domain simulation output
        vector <double> time;   // time of time domain Askaryan radiation
        vector <double> time_mimic;   // time of time domain Askaryan radiation (same time range with data)
        vector <double> V_mimic;    // signal + noise waveform which mimics the data (size : NFOUR/2 bin)

        int global_trig_bin; // from V_mimic [0, NFOUR/2] bins, where global trigger occured

        vector < vector < vector <double> > > V;   // For each ray individually, volt signal with all factors applied (from fft, excludes gain offse)
        vector <double> V_convolved;    // After convolution of all rays, volt signal with all factors applied (from Convolve_Signal, excludes gain offset)
        vector <double> V_noise;        // noise voltage waveform with all factors applied (from Convolve_Signal, excludes gain offset)

        vector < vector <int> > SignalExt; // flag if actual signal exist for the ray trace solution

        vector < vector <int> > SignalBin; // the bin number where the center of signal located. we can compare this value to Trig_Pass value to have likely triggered ray trace solution

        vector < vector <double> > SignalBinTime; ///< the time of center of bin where signal should locate after sim decided the readout window. MK added -2023-05-18-

        vector <int> noise_ID;      // information about which pure noise waveform is used for trigger analysis

        vector < vector <double> > PeakV;  // peak voltage in time domain
        vector <int> Rank;      // rank of peak voltage between antennas (Rank = 0 for 0 signal)

        //	vector <double> PeakV_fromFullWaveform; // peak voltage in time domain taken from full waveform, including noise at the time of signal insertion
        //	vector <int> Rank_fromFullWaveform;      // rank of peak voltage between antennas (Rank = 0 for 0 signal)

        int Trig_Pass; // 0 if not passed the trigger, 1 if passed the trigger
        //vector <int> Trig_Pass; // 0 if not passed the trigger, 1 if passed the trigger

        int Likely_Sol[2]; // comparing Trig_Pass and SignalBin value, this value returns which ray trace solution has been triggered (not perfect but most likely)

        int SingleChannelTriggers; // how many bins passed the threshold in this channel (should be equal to size of SCT_threshold_pass).
        vector <double> SCT_threshold_pass; // for each bin that passed, what was the threshold value at which it passed (for TRIG_SCAN_MODE only). 
	
        long TotalBinsScannedPerChannel;
	
        vector <int> TooMuch_Tdelay;    // 0 is PeakV is located inside the DATA_BIN_SIZE array,  1 is when PeakV is located outside the DATA_BIN_SIZE so that we can't correctly check if it is triggered or not

        void Prepare_Outputs(int n_interactions);
        void Find_Likely_Sol();
        int Get_Max_SignalBin();
        void Get_Brightest_Interaction(int (*brightest_event)[2]);
        void clear ();  // clear all vector format information for next event
        void clear_useless ( Settings *settings1 );  // clear all vector information which are useless

        ClassDef(Antenna_r,3);
};

class String_r {
    public:
        //int trg;    // if any antenna trigg in the event. 0 : no antenna in the string trg
                    //                                    1 : 1 or more antenna trg

        vector <Antenna_r> antennas;

        ClassDef(String_r,1);
};

class Station_r {
    public:
        //int trg;    // if any antenna trigg in the event. 0 : no antenna trg
                    //                                    1: 1 or more antenna trg 
        vector <String_r> strings;
        vector <Surface_antenna_r> surfaces;

        double min_arrival_time;    // for each station, minimum arrival time (include all ray_solves). this will be used for time delay between antennas.
        double max_arrival_time;    // for each station, maximum arrival time (include all ray_solves). this will be used for time delay between antennas.
        double max_PeakV;           // for each station, maximum PeakV value (include all ray_solves). this will also be used for time delay plot (to set same vertical scale)
        int Total_ray_sol;          // total number of ray_sols in the stations. If there is 0 Total_ray_sol, we don't need to do trigger check while there is any Total_ray_sol, we do trigger check.
        int Global_Pass;    // if global trigger passed or not: 0 = not passed, >0 passed, number indicates the first bin in the triggered window of the waveform at which the global trigger passed

        int total_trig_search_bin;  // total number of bins for searching trigger. 

//         int numChan;
// 	int numChanVpol;
// 	int numChanHpol;
        
        // TDR is for Tunnel Diode Response i.e. the value on which the trigger happened
        vector <double> TDR_all;
	vector <double> TDR_all_sorted;
	vector <double> TDR_Hpol_sorted;
	vector <double> TDR_Vpol_sorted;
        
        ClassDef(Station_r,3);
};

class CircularBuffer{

    public:
        int i;
        int mode; // in mode 1 only check number of values above threshold, in >1 check what the best value is too
        int changelog;// if the best value changed, this is returned by add and fill
        int N;// size of buffer
        double *buffer;
        double pthresh; // general run's pthresh
        double best_value;// best value
        double temp_value;// copy of best value, updates each call (so we can zero it when sorting)
        double epsilon; // best value needs to be this close to current last value
        double last_value; // value leaving buffer
        int addToNPass; // number of values above pthresh inside buffer

        CircularBuffer(int size, double threshold, int scan_mode) : mode(scan_mode), N(size), pthresh(threshold) 
          { i=0; best_value=0; temp_value=0; last_value=0; addToNPass=0; epsilon=1e-6; buffer=new double[N]; for(int j=0;j<N;j++) buffer[j]=0; }
        ~CircularBuffer(){ delete [] buffer; }

        int add(double input_value);
        int fill(double input_value);
        double findBestValue();
        int numBinsToOldestTrigger();
        int numBinsToLatestTrigger();

};


class Report {
    private:
        vector <double> noise_phase;    // random noise phase generated in GetNoisePhase()


        // variables we need for trigger
           // test selecting noise waveform

           int noise_pass_nogo; // index for checking if any same noise_ID is used in different chs.
           int N_noise;     // needed number of noise waveforms (most cases, we will need only 1)
           int noise_ID[5];    // selected noise waveform ID (we should not need 5 noise waveforms, but just in case)
           int ch_ID;   // channel ID
           //double Full_window[detector->params.number_of_strings_station * detector->params.number_of_antennas_string][settings1->DATA_BIN_SIZE];    // entire window for trigger check (diode convlv results for all antennas in a station)
           //vector < vector <double> > Full_window;  // entire window for trigger check (diode convlv results for all antennas in a station)
           int max_total_bin;   // to save time, use only necessary number of bins
           int remain_bin;      // the bin number for not using entire DATA_BIN_SIZE array
           vector < vector <int> > signal_bin;      // the center of bin where signal should locate

           int triggerCheckLoop(Settings *settings1, Detector *detector, Event *event, Trigger *trigger, int stationID, int trig_search_init, int max_total_bin, int trig_window_bin, int scan_mode=1);
// 	   int triggerCheckLoopScan();
// 	   int triggerCheckLoopScanNumbers();
	   
           int saveTriggeredEvent(Settings *settings1, Detector *detector, Event *event, Trigger *trigger, int stationID, int trig_search_init, int max_total_bin, int trig_window_bin, int last_trig_bin);

           vector < vector < vector <double> > > RayStep;


    public:
           /*
           double Full_window[16][16384];  // test with array, not vector, diode response
           double Full_window_V[16][16384];  // test with array, not vector, voltage waveform
           */
           vector <int> Passed_chs;
        //int trg;    // if any antenna in entire detectors trg. 0 : no antenna trg
                    //                                         1 : 1 or more antenna trg

        Report ();
        Report (Detector *detector, Settings *settings1);
        ~Report ();
    //make the UsefulIcrrStationEvent for use with AraRoot
    //UsefulIcrrStationEvent theUsefulEvent;

    
        void Initialize (Detector *detector, Settings *settings1);

        //void Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger);
        //void Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger, UsefulIcrrStationEvent *theUsefulEvent);


//        void Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger, UsefulIcrrStationEvent *theUsefulEvent);
    
//    void Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger);
    
    void Connect_Interaction_Detector_V2 (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Birefringence *birefringence, Settings *settings1, Trigger *trigger, int evt);     
    void rerun_event(Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, Birefringence *birefringence, IceModel *icemodel, Settings *settings, int which_solution,
        vector<int> &numSolutions, vector<vector<vector<double> > > &traceTimes, vector<vector<vector<double> > > &traceVoltages
        );
    void ModelRay(
        int ray_idx, vector< vector< double > > ray_output, int interaction_idx, double *T_forint, 
        Antenna_r *antenna_r, Antenna *antenna_d,  int i, int j, int k, 
        int debugmode,  Birefringence *birefringence, Detector *detector, 
        Event *event, IceModel *icemodel, Settings *settings, Signal *signal);
    void InitializeNNew(Antenna_r *antenna, int interaction_idx, int ray_idx, double dT, Settings *settings1);
    void GetRayParameters(
        Antenna_r *antenna_r, Antenna *antenna_d, Interaction interaction, int interaction_idx,
        int i, int j, int k, int ray_idx, vector<vector< double > > ray_output,
        Vector *n_trg_pokey, Vector *n_trg_slappy, Vector *Pol_vector_src,
        Position *launch_vector, Position *receive_vector,
        IceModel *icemodel, Settings *settings1 );
    void PropagateSignal(
        double dT_forfft, int efield_length, vector< double > efield_time, vector< double > efield, double *T_forint,
        int interaction_idx, int ray_idx, vector<vector< double > > ray_output, Position launch_vector, double time_diff_birefringence, 
        Vector Pol_vector_src, Vector Pol_vector, Vector n_trg_slappy, Vector n_trg_pokey, 
        Antenna_r *antenna_r, Antenna *antenna_d, int gain_ch_no, int j, int k,
        Birefringence *birefringence, Detector *detector, Event *event, IceModel *icemodel, Settings *settings);
    
    // Phased Array functions
    bool isTrigger(double eff);
    void checkPATrigger(
        int i, Detector *detector, Event *event, int evt, Trigger *trigger, Settings *settings1, 
        int trig_search_init, int max_total_bin);    
    double interpolate(double *xdata,double *ydata, double xi, int numData);
    
#ifdef ARA_UTIL_EXISTS

    void MakeUsefulEvent(Detector *detector, Settings *settings1, Trigger *trigger, int stationID, int stationIndex, UsefulIcrrStationEvent *theUsefulEvent);
    void MakeUsefulEvent(Detector *detector, Settings *settings1, Trigger *trigger, int stationID, int stationIndex, UsefulAtriStationEvent *theUsefulEvent);
#endif
    
    void ClearUselessfromConnect(Detector *detector, Settings *settings1, Trigger *trigger);

        // Signal+noise convolution functions
        void Convolve_Signals(    
            Antenna_r *antenna, int channel_number, int station_number,
            Event *event, Settings *settings1, Trigger *trigger, Detector *detector);
        void Select_Wave_Convlv_Exchange( // Convolve the signal from 1 ray
            vector <double> &V, 
            int BINSIZE, vector <double> *V_signal); 
        void Select_Wave_Convlv_Exchange( // Convolve the signal from 2 rays
            int signalbin_1, int signalbin_2, 
            vector <double> &V1, vector <double> &V2, 
            int BINSIZE, vector <double> *V_signal); 
        void Select_Wave_Convlv_Exchange( // Convolve the signal from 3 rays
            int signalbin_0, int signalbin_1, int signalbin_2, 
            vector <double> &V0, vector <double> &V1, vector <double> &V2, 
            int BINSIZE, vector <double> *V_signal);
        void Combine_Waveforms( // combine the signal from two vectors into single vector
            int signalbin_0, int signalbin_1,
            vector<double> V0, vector<double> V1,
            int* signalbin_combined, vector<double>* V_combined);
        void GetNoiseThenConvolve(
            Antenna_r *antenna, vector <double> V_signal,
            int BINSIZE, int this_signalbin, int n_connected_rays, 
            int channel_index, int station_number, 
            Settings *settings1, Trigger *trigger, Detector *detector);
        void GetAntennaNoiseWF(
            int signalbin, 
            int wf_length, int BINSIZE, int ID, int StationIndex, vector <double> *V_noise_only,
            Settings *settings1, Trigger *trigger, Detector *detector);


        void Apply_Gain_Offset(Settings *settings1, Trigger *trigger, Detector *detector, int ID, int StationIndex); // we need to apply a gain offset to the basic waveforms.



        int GetChNumFromArbChID( Detector *detector, int ID, int StationIndex, Settings *settings1);// get actual ch number from arb chID

        Vector GetPolarization (Vector &nnu, Vector &launch_vector);

        void GetParameters (Position &src, Position &trg, Vector &nnu, double &viewangle, double receive_angle, Vector &launch_vector, Vector &receive_vector, Vector &n_trg_slappy, Vector &n_trg_pokey );    // get viewangle, launch, receive vectors  (it reads launch angle as a viewangle and returns actual viewangle)

        double GaintoHeight(double gain, double freq, double n_medium, double Z_A=50);
        
        double calculatePolFactor(Vector &Pol_vector, int ant_type, double antenna_theta, double antenna_phi);

        void ApplyAntFactors(double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &pol_factor, double &vmmhz, double antenna_theta, double antenna_phi);

        void ApplyAntFactors_Tdomain(double AntPhase, double heff, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img, Settings *settings1, double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode=false, bool applyInverse=false);

        void ApplyAntFactors_Tdomain_FirstTwo ( double heff, double heff_lastbin, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_bin0, double &vm_bin1, double antenna_theta, double antenna_phi,  double freq, bool useInTransmitterMode=false, bool applyInverse=false);
    
        void InvertAntFactors_Tdomain(double AntPhase, double heff, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img, Settings *settings1, double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode=false);

        void InvertAntFactors_Tdomain_FirstTwo ( double heff, double heff_lastbin, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_bin0, double &vm_bin1, double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode=false);


        void ApplyElect_Tdomain(double freq, Detector *detector, double &vm_real, double &vm_img, int gain_ch_no, Settings *settings1, bool applyInverse=false);

        void ApplyElect_Tdomain_FirstTwo(double freq0, double freq1, Detector *detector, double &vm_bin0, double &vm_bin1, int gain_ch_no, Settings *settings1, bool applyInverse=false);
    
        void InvertElect_Tdomain(double freq, Detector *detector, double &vm_real, double &vm_img, int gain_ch_no, Settings *settings1);

        void InvertElect_Tdomain_FirstTwo(double freq0, double freq1, Detector *detector, double &vm_bin0, double &vm_bin1, int gain_ch_no, Settings *settings1);
    
        void ApplySplitterFactor(double &vm_real, double &vm_img, Detector *detector, Settings *settings1, bool applyInverse=false);



        void ApplyFilter(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFilter_databin(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFilter_NFOUR(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFilter_OutZero (double freq, Detector *detector, double &vmmhz);


        // apply gain in Preamp
        void ApplyPreamp(int bin_n, Detector *detector, double &vmmhz);
        void ApplyPreamp_databin(int bin_n, Detector *detector, double &vmmhz);
        void ApplyPreamp_NFOUR(int bin_n, Detector *detector, double &vmmhz);
        void ApplyPreamp_OutZero (double freq, Detector *detector, double &vmmhz);

	void ApplyNoiseFig_databin(int ch, int bin_n, Detector *detector, double &vmmhz, Settings *settings1);

        // apply gain in FOAM
        void ApplyFOAM(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFOAM_databin(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFOAM_NFOUR(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFOAM_OutZero (double freq, Detector *detector, double &vmmhz);


        // apply RFCM gain
        void ApplyRFCM(int ch, int bin_n, Detector *detector, double &vmmhz, double RFCM_OFFSET);
        void ApplyRFCM_databin(int ch, int bin_n, Detector *detector, double &vmmhz, double RFCM_OFFSET);


        void GetAngleAnt(Vector &rec_vector, Position &antenna, double &ant_theta, double &ant_phi);
        void GetAngleLaunch(Vector &launch_vector, double &launch_theta, double &launch_phi);

        // Noise Functions
        void GetNoiseWaveforms(Settings *settings1, Detector *detector, double vhz_noise, double *vnoise);
        void GetNoiseWaveforms_ch(Settings *settings1, Detector *detector, double vhz_noise, double *vnoise, int ch);
        void GetNoisePhase(Settings *settings1);
        void Prepare_Antenna_Noise(    
            int debugmode, int ch_ID, 
            int station_number, int string_number, int antenna_number,
            Settings *settings1, Trigger *trigger, Detector *detector
        );

        void MakeArraysforFFT(Settings *settings1, Detector *detector, int StationIndex, vector <double> &vsignal_array, double *vsignal_forfft);
        void MakeArraysforFFT(Settings *settings1, Detector *detector, int StationIndex, double *vsignal_array, double *vsignal_forfft);

        void MakeArraysforFFT_noise(Settings *settings1, Detector *detector,  int StationIndex, vector <double> &vsignal_array, double *vsignal_forfft);


        double FindPeak (double *waveform, int n);  // same with icemc; trigger->AntTrigger::FindPeak
        double FindPeak(vector< double > waveform, int n);

        void SetRank(Detector *detector); // set rank (rank of strength of signal at each antenna)



        int GetChannelNum8_LowAnt(int string_num, int antenna_num); // just return ch numbers 1-8 for antenna 0-1 (bottom antennas) and higher ch numbers for antenna 2-3 (top antennas) this is used for only TRIG_ONLY_LOW_CH_ON=1 mode with 


	TGraph *getWaveform(Detector *detector, int ch, int station_i=0, int event_num=0, int run_num=0);

	vector<TGraph*> getWaveformVector(Detector *detector, int station_i=0, int event_num=0, int run_num=0);
	vector<TGraph*> getWaveformVectorVpol(Detector *detector, int station_i=0, int event_num=0, int run_num=0);
	vector<TGraph*> getWaveformVectorHpol(Detector *detector, int station_i=0, int event_num=0, int run_num=0);

       int getNumOfSignalledAnts(Station_r station);

        double get_SNR(vector<double> signal_array, vector<double> noise_array);
	
        vector <double> Vfft_noise_after;   // noise Vfft after get_random_rician
        vector <double> Vfft_noise_before;   // noise Vfft before get_random_rician
        //vector <double> V_noise_timedomain;   // noise V timedomain after get_random_rician and inverse fft
        double Vfft_noise_org;              // V/Hz for thermal noise from Johnson-Nyquist

        vector <double> V_total_forconvlv; // vector array for pure signal diode convlv result


        void clear_useless(Settings *settings1);   // to reduce the size of output AraOut.root, remove some information

        void delete_all(); // delete all informations test

        vector <Station_r> stations;
        vector <String_r> strings;

        double RandomTshift; // for t-domain signal, a factor for random init time shift

        double init_T; // locate zero time at the middle and give random time shift (for interpolated waveforms)

        // Phased Array variables
        double pa_force_trigger_snr = 3.5; 
            // SNR that should always trigger, 
            // used (eg) when triggering on noise only events
        double pa_snr_cap = 25.;
            // KAH thinks this is the max SNR she had efficiencies calculated for
            // KAH says the PA has a SNR cap of 25 in practice and may have a link showing this.

        ClassDef(Report,1);

};

#endif  //REPORT_H
