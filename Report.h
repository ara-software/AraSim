//--------------------------------------------------
// class Antenna_Response
//-------------------------------------------------- 
//

#ifndef REPORT_H
#define REPORT_H

#include <vector>

#ifndef __CINT__
//Include output format to enable reading by analysis software AraRoot
#ifdef ARA_UTIL_EXISTS
#include "UsefulIcrrStationEvent.h"
#endif
#endif

#ifdef ARA_UTIL_EXISTS
class UsefulIcrrStationEvent;
#endif

class Detector;
class Event;
class RaySolver;
class Signal;
class IceModel;
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

        vector <double> view_ang;    //viewing angle
        vector <double> launch_ang;  //launch angle
        vector <double> rec_ang;     //receiving angle
        vector <double> reflect_ang; // surface reflection angle (if 100 : no reflection case)
        vector <double> Dist;        //Distance between posnu and antenna
        vector <double> L_att;        //Attenuation factor
        vector <double> arrival_time;        //time from posnu to antenna (t:0 at posnu)
        vector <int> reflection;     // non-reflected : 0,  reflected : 1
        vector <Position> Pol_vector;   // polarization vector at the antenna
        //vector <Position> n_H;  // normalized vector for H pol
        //vector <Position> n_V;  // normalized vector for V pol

        // below freq domain simulation output
        vector < vector <double> > vmmhz;  // signal V/m/MHz for each freq bin
        //
        vector < vector <double> > Heff;  // effective height for each freq bin
        vector <double> Mag;  // magnification factor
        vector <double> Fresnel;  // Fresnel factor
        vector <double> Pol_factor;  // Polarization factor
        
        vector < vector <double> > Vm_zoom;  // E field before ant T-domain
        vector < vector <double> > Vm_zoom_T;  // E field before ant T-domain time
        //vector < vector <double> > Vm_wo_antfactor;  // before applying ApplyAntFactors
        //vector < vector <double> > VHz_antfactor;  // after applying ApplyAntFactors to vmmhz above ( 1/sqrt2 * 1/dt * 0.5 * heff * pol_factor )
        //vector < vector <double> > VHz_filter;  // after applying ApplyAntFactors above and then apply filter gain from detector->GetFilterGain
        //
        int skip_bins[2]; // for two ray sols

        int Nnew[2]; // new number of bins for V_fotfft array

        vector < vector <double> > Vfft;  // signal V preparing for FFT
        vector < vector <double> > Vfft_noise;  // noise V preparing for FFT


        // below time domain simulation output
        vector <double> time;   // time of time domain Askaryan radiation
        vector <double> time_mimic;   // time of time domain Askaryan radiation (same time range with data)
        vector <double> V_mimic;    // signal + noise waveform which mimics the data (size : NFOUR/2 bin)

        int global_trig_bin; // from V_mimic [0, NFOUR/2] bins, where global trigger occured

        vector < vector <double> > Ax;     // vector potential x component
        vector < vector <double> > Ay;
        vector < vector <double> > Az;
        vector < vector <double> > V;   // volt signal with all factors applied (as far as we can) (from fft)

        vector <int> SignalExt; // flag if actual signal exist for the ray trace solution

        vector <int> SignalBin; // the bin number where the center of signal located. we can compare this value to Trig_Pass value to have likely triggered ray trace solution


        vector <int> noise_ID;      // information about which pure noise waveform is used for trigger analysis

        //vector < vector <double> > V_noise; // volt noise signal (with all factors applied as far as we can) (from thermal noise + fft)

        //vector < vector <double> > V_total; // volt signal + noise with all factors applied as far as we can

        //vector < vector <double> > V_total_diode;   // volt signal + noise with all factors (as far as we can) and convlution with diode (time domain)

        //vector < vector <double> > V_total_timedelay;   // volt signal + noise with all factors applied and time delay between antennas
        //
        //
        vector <double> PeakV;  // peak voltage in time domain
        vector <int> Rank;      // rank of peak voltage between antennas (Rank = 0 for 0 signal)
        int Trig_Pass; // 0 if not passed the trigger, 1 if passed the trigger
        //vector <int> Trig_Pass; // 0 if not passed the trigger, 1 if passed the trigger

        int Likely_Sol; // comparing Trig_Pass and SignalBin value, this value returns which ray trace solution has been triggered (not perfect but most likely)


        vector <int> TooMuch_Tdelay;    // 0 is PeakV is located inside the DATA_BIN_SIZE array,  1 is when PeakV is located outside the DATA_BIN_SIZE so that we can't correctly check if it is triggered or not

        
        void clear ();  // clear all vector format information for next event
        void clear_useless ( Settings *settings1 );  // clear all vector information which are useless

        ClassDef(Antenna_r,1);
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
           vector <int> signal_bin;      // the center of bin where signal should locate
           vector <int> signal_dbin;     // the bin difference between signal bins
           vector <int> connect_signals;    // if ray_sol time delay is small enough to connect each other

           int triggerCheckLoop(Settings *settings1, Detector *detector, Event *event, Trigger *trigger, int stationID, int trig_search_init, int max_total_bin, int trig_window_bin, int scan_mode=1);
// 	   int triggerCheckLoopScan();
// 	   int triggerCheckLoopScanNumbers();
	   
           int saveTriggeredEvent(Settings *settings1, Detector *detector, Event *event, Trigger *trigger, int stationID, int trig_search_init, int max_total_bin, int trig_window_bin, int last_trig_bin);


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
    
    void Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger, int evt);    
    
    
#ifdef ARA_UTIL_EXISTS

    void MakeUsefulEvent(Detector *detector, Settings *settings1, Trigger *trigger, int stationID, UsefulIcrrStationEvent *theUsefulEvent);
#endif
    
    void ClearUselessfromConnect(Detector *detector, Settings *settings1, Trigger *trigger);

    
        void Select_Wave_Convlv_Exchange(Settings *settings1, Trigger *trigger, Detector *detector, int signalbin, vector <double> &V, int *noise_ID, int ID, int StationIndex);   // literally get noise waveform from trigger class and add signal voltage "V" and do convlv. convlv result will replace the value in Full_window array
        
        void Select_Wave_Convlv_Exchange(Settings *settings1, Trigger *trigger, Detector *detector, int signalbin_1, int signalbin_2, vector <double> &V1, vector <double> &V2, int *noise_ID, int ID, int StationIndex);   // literally get noise waveform from trigger class and add signal voltage "V" and do convlv. convlv result will replace the value in Full_window array

        void Select_Wave_Convlv_Exchange(Settings *settings1, Trigger *trigger, Detector *detector, int signalbin_0, int signalbin_1, int signalbin_2, vector <double> &V0, vector <double> &V1, vector <double> &V2, int *noise_ID, int ID, int StationIndex);   // literally get noise waveform from trigger class and add signal voltage "V" and do convlv. convlv result will replace the value in Full_window array


        void Apply_Gain_Offset(Settings *settings1, Trigger *trigger, Detector *detector, int ID, int StationIndex); // we need to apply a gain offset to the basic waveforms.



        int GetChNumFromArbChID( Detector *detector, int ID, int StationIndex, Settings *settings1);// get actual ch number from arb chID

        Vector GetPolarization (Vector &nnu, Vector &launch_vector);

        void GetParameters (Position &src, Position &trg, Vector &nnu, double &viewangle, double receive_angle, Vector &launch_vector, Vector &receive_vector, Vector &n_trg_slappy, Vector &n_trg_pokey );    // get viewangle, launch, receive vectors  (it reads launch angle as a viewangle and returns actual viewangle)

        double GaintoHeight(double gain, double freq, double n_medium);

        void ApplyAntFactors(double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &pol_factor, double &vmmhz);

        void ApplyAntFactors_Tdomain(double AntPhase, double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img, Settings *settings1);

        void ApplyAntFactors_Tdomain_Transmitter(double AntPhase, double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img, Settings *settings1);

        void ApplyAntFactors_Tdomain_FirstTwo ( double heff, double heff_lastbin, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_bin0, double &vm_bin1);


        void ApplyElect_Tdomain(double freq, Detector *detector, double &vm_real, double &vm_img, Settings *settings1);

        void ApplyElect_Tdomain_FirstTwo(double freq0, double freq1, Detector *detector, double &vm_bin0, double &vm_bin1);



        void ApplyFilter(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFilter_databin(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFilter_NFOUR(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFilter_OutZero (double freq, Detector *detector, double &vmmhz);


        // apply gain in Preamp
        void ApplyPreamp(int bin_n, Detector *detector, double &vmmhz);
        void ApplyPreamp_databin(int bin_n, Detector *detector, double &vmmhz);
        void ApplyPreamp_NFOUR(int bin_n, Detector *detector, double &vmmhz);
        void ApplyPreamp_OutZero (double freq, Detector *detector, double &vmmhz);

        // apply gain in FOAM
        void ApplyFOAM(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFOAM_databin(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFOAM_NFOUR(int bin_n, Detector *detector, double &vmmhz);
        void ApplyFOAM_OutZero (double freq, Detector *detector, double &vmmhz);


        // apply RFCM gain
        void ApplyRFCM(int ch, int bin_n, Detector *detector, double &vmmhz, double RFCM_OFFSET);
        void ApplyRFCM_databin(int ch, int bin_n, Detector *detector, double &vmmhz, double RFCM_OFFSET);


        void GetAngleAnt(Vector &rec_vector, Position &antenna, double &ant_theta, double &ant_phi);

        void GetNoiseWaveforms(Settings *settings1, Detector *detector, double vhz_noise, double *vnoise);
        void GetNoiseWaveforms_ch(Settings *settings1, Detector *detector, double vhz_noise, double *vnoise, int ch);

        void GetNoisePhase(Settings *settings1);

        void MakeArraysforFFT(Settings *settings1, Detector *detector, int StationIndex, vector <double> &vsignal_array, double *vsignal_forfft);
        void MakeArraysforFFT(Settings *settings1, Detector *detector, int StationIndex, double *vsignal_array, double *vsignal_forfft);

        void MakeArraysforFFT_noise(Settings *settings1, Detector *detector,  int StationIndex, vector <double> &vsignal_array, double *vsignal_forfft);


        double FindPeak (double *waveform, int n);  // same with icemc; trigger->AntTrigger::FindPeak

        void SetRank(Detector *detector); // set rank (rank of strength of signal at each antenna)



        int GetChannelNum8_LowAnt(int string_num, int antenna_num); // just return ch numbers 1-8 for antenna 0-1 (bottom antennas) and higher ch numbers for antenna 2-3 (top antennas) this is used for only TRIG_ONLY_LOW_CH_ON=1 mode with 




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


        // test T domain waveform
        //static const int outbin = 50;
        static const int outbin = 64;
        double Tarray[outbin];
        double Earray[outbin];

        double init_T; // locate zero time at the middle and give random time shift (for interpolated waveforms)


        ClassDef(Report,1);

};

#endif  //REPORT_H
