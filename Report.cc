#include "Detector.h"
#include "Report.h"
#include "Event.h"
#include "RaySolver.h"
#include "signal.hh"
#include "IceModel.h"
#include "Settings.h"
#include "Vector.h"
#include "Tools.h"
#include "Trigger.h"
#include "Constants.h"
#include "Birefringence.h"

#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "TRandom3.h"
#include "TMath.h"

#include <cstdlib>

ClassImp(Report);
ClassImp(Antenna_r);
ClassImp(Surface_antenna_r);
ClassImp(String_r);
ClassImp(Station_r);
ClassImp(CircularBuffer);

Report::Report() {
}

Report::Report(Detector *detector, Settings *settings1) {
    // Default constructor

    Initialize(detector, settings1);
    
}



Report::~Report() {
    // default destructor

    stations.clear();
    strings.clear();

    Passed_chs.clear();
    Vfft_noise_after.clear();
    Vfft_noise_before.clear();
    //V_noise_timedomain.clear();

}

void Report::delete_all() {

    stations.clear();
    strings.clear();

    Passed_chs.clear();
    Vfft_noise_after.clear();
    Vfft_noise_before.clear();

    noise_phase.clear();    // random noise phase generated in GetNoisePhase()
    signal_bin.clear();      // the center of bin where signal should locate
    signal_dbin.clear();     // the bin difference between signal bins
    connect_signals.clear();    // if ray_sol time delay is small enough to connect each other
    Passed_chs.clear();


}

void Report::Initialize(Detector *detector, Settings *settings1) {
    
    // clear information stored in (but there shouldn't be. just to make sure)
    //
    //stations.clear();
    delete_all();


    // tmp for push_back vector structure
    Antenna_r tmp_antenna;
    String_r tmp_string;
    Station_r tmp_station;
    Surface_antenna_r tmp_surface;

    for (int i=0; i<detector->params.number_of_stations; i++) {
        // vector stations
        stations.push_back(tmp_station);
        for (int j=0; j<detector->stations[i].strings.size(); j++) {
            // vector strings
            stations[i].strings.push_back(tmp_string);
            for (int k=0; k<detector->stations[i].strings[j].antennas.size(); k++) {
                // vector antennas
                stations[i].strings[j].antennas.push_back(tmp_antenna);
            }
        }
        for (int j=0; j<detector->stations[i].surfaces.size(); j++) {
            // vector surface antennas
            stations[i].surfaces.push_back(tmp_surface);
        }

	if(settings1->TRIG_SCAN_MODE>0){// scan Pthresh mode

	  int numChan=0; 
	  int numChanVpol=0;
	  int numChanHpol=0;
	  
	  for(int j=0;j<detector->stations[i].strings.size(); j++){
	    
	    for (int k=0;k<detector->stations[i].strings[j].antennas.size();k++) {
	      
	      int string_i = detector->getStringfromArbAntID( i, numChan);
              int antenna_i = detector->getAntennafromArbAntID( i, numChan);
	      if (detector->stations[i].strings[string_i].antennas[antenna_i].type == 0) numChanVpol++;
	      if (detector->stations[i].strings[string_i].antennas[antenna_i].type == 1) numChanHpol++;
	      numChan++;
	    }// for k
	    
	  }// for j
	  
	  stations[i].TDR_all.clear();
	  for(int ch=0;ch<numChan; ch++ ) stations[i].TDR_all.push_back(0);
	  stations[i].TDR_all_sorted.clear();
	  if(settings1->TRIG_MODE==0) for(int ch=0;ch<numChan; ch++ ) stations[i].TDR_all_sorted.push_back(0);

	  stations[i].TDR_Vpol_sorted.clear();
	  if(settings1->TRIG_MODE==1) for(int ch=0;ch<numChanVpol; ch++) stations[i].TDR_Vpol_sorted.push_back(0);

	  stations[i].TDR_Hpol_sorted.clear();
	  if(settings1->TRIG_MODE==1) for(int ch=0;ch<numChanHpol; ch++ ) stations[i].TDR_Hpol_sorted.push_back(0);

	  
	}// if TRIG_SCAN_MODE 
	
    }// for i (number of stations)




}



void Antenna_r::clear() {   // if any vector variable added in Antenna_r, need to be added here!
    
    view_ang.clear();
    launch_ang.clear();
    rec_ang.clear();
    reflect_ang.clear();
    Dist.clear();
    L_att.clear();
    arrival_time.clear();
    reflection.clear();
    Pol_vector.clear();
    vmmhz.clear();
    Heff.clear();
    Mag.clear();
    Fresnel.clear();
    Pol_factor.clear();
    Pol_factorH.clear();
    Pol_factorV.clear();
    phi_rec.clear();
    theta_rec.clear();
    //VHz_antfactor.clear();
    //VHz_filter.clear();
    Vfft.clear();
    Vfft_noise.clear();

    ray_step.clear();

    time.clear();
    time_mimic.clear();
    V_noise.clear();
    V_convolved.clear();
    V_mimic.clear();
    Ax.clear();
    Ay.clear();
    Az.clear();

    V.clear();

    noise_ID.clear();

    PeakV.clear();
    Rank.clear();
    TooMuch_Tdelay.clear();

    Trig_Pass = 0;


    // additional for before ant waveform
    //Vm_wo_antfactor.clear();
    Vm_zoom.clear();
    Vm_zoom_T.clear();

    SignalBin.clear();
    SignalBinTime.clear();
    SignalExt.clear(); 
    
    SCT_threshold_pass.clear();
    
}



void Antenna_r::clear_useless(Settings *settings1) {   // to reduce the size of output AraOut.root, remove some information


    if (settings1->DATA_SAVE_MODE == 1) {
    Heff.clear();
    
    //VHz_antfactor.clear();
    //VHz_filter.clear();
    Vfft.clear();
    Vfft_noise.clear();

    Ax.clear();
    Ay.clear();
    Az.clear();

    V.clear();

    //Trig_Pass.clear();
    TooMuch_Tdelay.clear();

    // need or not?
    //Pol_vector.clear();
    vmmhz.clear();
    //Mag.clear();
    //Fresnel.clear();
    //Pol_factor.clear();
    //
    // additional for before ant waveform
    //Vm_wo_antfactor.clear();
    Vm_zoom.clear();
    Vm_zoom_T.clear();

    }
    else if (settings1->DATA_SAVE_MODE == 2) {
    
    //! clear the ray step to reduce the size of output AraOut.root
    ray_step.clear();  
    
    Heff.clear();
    //VHz_antfactor.clear();
    //VHz_filter.clear();
    Vfft.clear();
    Vfft_noise.clear();

    Ax.clear();
    Ay.clear();
    Az.clear();

    V.clear();

    //Trig_Pass.clear();
    TooMuch_Tdelay.clear();

    // need or not?
    //Pol_vector.clear();
    vmmhz.clear();
    //Mag.clear();
    //Fresnel.clear();
    //Pol_factor.clear();
    
    // clear global trigger waveform info also
    time.clear();
    time_mimic.clear();
    V_convolved.clear();
    V_noise.clear();
    V_mimic.clear();


    // additional for before ant waveform
    //Vm_wo_antfactor.clear();
    Vm_zoom.clear();
    Vm_zoom_T.clear();
    
    }


}

    
int CircularBuffer::add(double input_value){

    changelog=0;
    temp_value=best_value;
    
    if(buffer[i]<pthresh) addToNPass--; // if the value leaving the buffer is over threshold, we reduce the counter.
    if(mode>1) last_value=buffer[i]; // in mode 1 we don't care about the values, just about the addToNPass
        
    if(input_value<pthresh){
        addToNPass++; // if value entering buffer is over threshold we increase the counter.
        buffer[i]=input_value;
    }
    else buffer[i]=0; // if under threshold just insert zero
      
    // improve best threshold value in the buffer
    if ( mode>1 && buffer[i]<best_value ) { 
        best_value=buffer[i]; 
        temp_value=best_value; 
        changelog=1;
    }

    if( mode>1 && std::fabs(last_value-best_value)<epsilon ) { 
        // i.e. if last_value==best_value. this happens when the best value leaves buffer
        best_value=findBestValue(); // rescan whole buffer. this should be rather rare...
        temp_value=best_value;
        changelog=1;
    }
       
    i++;
    if(i==N) i=0;
        
    return changelog; 
    // in mode>1: return value is zero if no changes, 
    // and non-zero if there are changes to best_value (for mode==1 it's always zero);
      
}// add 

int CircularBuffer::fill(double input_value){// buffer is just filling, don't use last_value
    
    changelog=0;
    temp_value=best_value;
    
    if(input_value<pthresh){
        addToNPass++; // if value entering buffer is over threshold we increase the counter.
        buffer[i]=input_value;
        // cout<<"addToNPass++\n";
    }
    else buffer[i]=0; // if under threshold add just zero
    
    if(mode>1&&buffer[i]<best_value){ // improve best threshold value in the buffer
        best_value=buffer[i]; 
        temp_value=best_value; changelog=1;
    }

    i++;
    if(i==N) i=0;

    return changelog;

}// fill
    
double CircularBuffer::findBestValue(){
    double temp_best=0;
    for(int ii=0;ii<N;ii++) {
        if ( buffer[ii]<temp_best ) {
            temp_best=buffer[ii];
        }
    }
    return temp_best;
}
    
int CircularBuffer::numBinsToOldestTrigger(){
    
    if(addToNPass<1){ 
        cerr<<"ERROR: there are no triggers in this buffer right now!\n"; 
        return -1;
    }

    int j; // the number of bins after "i" where the first trigger is found 
    for(j=1; j<N; j++){
        // i-1 is b/c we did i++ in the last call to add/fill.
        if(buffer[(i-1+j)%N]<0) break; // if there's any value here, its because it passed the threshold, so take it
    } // for j

    return N-j; // the backward count of how many bins between i and the earliest trigger... 
    
} // numBinsToOldestTrigger
    
int CircularBuffer::numBinsToLatestTrigger(){
    
    if(addToNPass<1){ 
        cerr<<"ERROR: there are no triggers in this buffer right now!\n"; 
        return -1;
    }

    int j; // the number of bins after "i" where the first trigger is found 
    for(j=0; j<N; j++){

        int bin=i-1-j;

        if(bin<0) bin=N+(i-j);

        if(buffer[bin]<0){// if there's any value here, its because it passed the threshold, so take it
            break;
        }
    }// for j
    
    return j; // the count of how many bins between i and the latest trigger... 
    
}// numBinsToLatestTrigger



void Report::clear_useless(Settings *settings1) {   // to reduce the size of output AraOut.root, remove some information



    if (settings1->DATA_SAVE_MODE != 0) {

        // also clear all vector info to reduce output root file size
        noise_phase.clear();
        signal_bin.clear();
        signal_dbin.clear();
        connect_signals.clear();
        Passed_chs.clear();
        Vfft_noise_after.clear();
        Vfft_noise_before.clear();
        //V_noise_timedomain.clear();
        // done clear vector info in report head
        //
    
        V_total_forconvlv.clear();
	RayStep.clear();

    }

}

void Report::Connect_Interaction_Detector_V2(Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Birefringence *birefringence, Settings *settings1, Trigger *trigger, int evt) {


    int ray_sol_cnt;
    double viewangle;
    Position launch_vector; // direction of ray at the source
    Position receive_vector;    // direction of ray at the target antenna
    Vector n_trg_pokey; // unit pokey vector at the target
    Vector n_trg_slappy;    // unit slappy vector at the target
    vector<vector < double>> ray_output;
    double noise_rms;
    double all_receive_ang[2];

    double vmmhz1m_tmp, vmmhz1m_sum, vmmhz1m_em;    // currently not using vmmhz1m_em
    Position Pol_vector;    // polarization vector at the source
    double mag; // magnification factor. it can vary in case of plane / spherical wave
    double fresnel; // fresnel factor
    double Pol_factor;  // polarization factor
    double tmp; // for non use return values
    double Pol_factorV;
    double Pol_factorH;
    double phi_rec;
    double theta_rec;

    double freq_tmp, heff, antenna_theta, antenna_phi, launch_theta, launch_phi;  // values needed for apply antenna gain factor and prepare fft, trigger
    double heff_Tx, Tx_theta, Tx_phi, heff_Tx_lastbin; // Values for transmitting antenna mode.
    double volts_forfft[settings1->NFOUR / 2];  // array for fft
    double dT_forfft;
    double volts_forint[settings1->NFOUR / 2];  // array for interpolation
    double T_forint[settings1->NFOUR / 2];  // array for interpolation

    double dF_NFOUR = 1. / ((double)(settings1->NFOUR / 2) *settings1->TIMESTEP);   // in Hz

    int waveformLength = settings1->WAVEFORM_LENGTH;
    int waveformCenter = settings1->WAVEFORM_CENTER;

    double dF_Nnew;

    double heff_lastbin;
    double freq_lastbin;

    int check_toomuch_Tdelay;   // return value from MixSignalNoise_Tdelay

    double min_arrival_time_tmp;    // min arrival time between all antennas, raysolves
    double max_arrival_time_tmp;    // max arrival time between all antennas, raysolves
    double max_PeakV_tmp;   // max PeakV of all antennas in the station

    int trig_window_bin = (int)(settings1->TRIG_WINDOW / settings1->TIMESTEP);  // coincidence window bin for trigger

    RandomTshift = gRandom->Rndm();

    init_T = settings1->TIMESTEP *-1.e9 *((double) settings1->NFOUR / 4 + RandomTshift);    // locate zero time at the middle and give random time shift

    // decide whether debug mode or not
    int debugmode = 0;
    if (settings1->DEBUG_MODE_ON == 1 && evt < settings1->DEBUG_SKIP_EVT) debugmode = 1;
    else if (settings1->DEBUG_MODE_ON == 1 && evt >= settings1->DEBUG_SKIP_EVT) cout << evt << " " << endl;
    // skip most of computation intensive processes if debugmode == 1

    int N_pass; // number of trigger passed channels (antennas)
    int N_pass_V;   // number of trigger passed channels (Vpol antennas)
    int N_pass_H;   // number of trigger passed channels (Hpol antennas)

    for (int i = 0; i < detector->params.number_of_stations; i++)
    {

        min_arrival_time_tmp = 10.; // first min_arrival_time is unreasonably big value
        max_arrival_time_tmp = 0.;  // first max_arrival_time is unreasonably small value
        max_PeakV_tmp = 0.; // first max_PeakV_tmp is 0.

        stations[i].Total_ray_sol = 0;  // initial Total_ray_sol value

        for (int j = 0; j < detector->stations[i].strings.size(); j++)
        {
            
            init_T = settings1->TIMESTEP *-1.e9 *((double) settings1->NFOUR / 4 + RandomTshift);    // locate zero time at the middle and give random time shift
            for (int n = 0; n < settings1->NFOUR / 2; n++)
            {
                T_forint[n] = init_T + (double) n *settings1->TIMESTEP *1.e9;   // in ns
            }

            for (int k = 0; k < detector->stations[i].strings[j].antennas.size(); k++)
            {
                // cout << i << " : " << j << " : " << k << endl;

                stations[i].strings[j].antennas[k].clear(); // clear data in antenna which stored in previous event

		// This (gain_ch_no) is used for per-channel gain implementation. 
		// It is used in all instances of ApplyElect_Tdomain() and ApplyElect_Tdomain_FirstTwo(), to indicate channel number
		// Note that channel numbering is different for DETECTOR==4 than for the other modes (1-3). See that in the definition of GetChannelfromStringAntenna() 
		int gain_ch_no;
		if (settings1->DETECTOR==4 || settings1->DETECTOR==5){
			gain_ch_no = detector->GetChannelfromStringAntenna (i, j, k, settings1)-1;
		}
		else{
			gain_ch_no = detector->GetChannelfromStringAntenna (i, j, k, settings1);
		}
		
		// run ray solver, see if solution exist
                // if not, skip (set something like Sol_No = 0;
                // if solution exist, calculate view angle and calculate TaperVmMHz

                // added one more condition to run raysolver (direct distance is less than 3km)
                //

                // cout << event->Nu_Interaction[0].posnu.GetX() << " : " << event->Nu_Interaction[0].posnu.GetY() << " : " << event->Nu_Interaction[0].posnu.GetZ() << endl;
                // cout << event->Nu_Interaction[0].pickposnu << " : " << event->Nu_Interaction[0].posnu.Distance(detector->stations[i].strings[j].antennas[k]) << " : " << settings1->RAYSOL_RANGE << endl; 
                if (event->Nu_Interaction[0].pickposnu && event->Nu_Interaction[0].posnu.Distance(detector->stations[i].strings[j].antennas[k]) <= settings1->RAYSOL_RANGE)
                {
                    // if posnu is selected inside the antarctic ice
                    // cout << i << " : " << j << " : " << k << endl;

                    RayStep.clear();    // remove previous values
                    raysolver->Solve_Ray(event->Nu_Interaction[0].posnu, detector->stations[i].strings[j].antennas[k], icemodel, ray_output, settings1, RayStep);   // solve ray between source and antenna

                    ray_sol_cnt = 0;

                    if (raysolver->solution_toggle)
                    {
                        // if there are solution from raysolver

                        while (ray_sol_cnt < ray_output[0].size())
                        {
                            // for number of soultions (could be 1 or 2)


			    double time_diff_birefringence = birefringence->Time_Diff_TwoRays(RayStep[ray_sol_cnt][0], RayStep[ray_sol_cnt][1], ray_output[3][ray_sol_cnt], event->Nu_Interaction[0].posnu_from_antcen, settings1); // calculate time differences for birefringence 


                            stations[i].strings[j].antennas[k].arrival_time.push_back(ray_output[4][ray_sol_cnt]);

                            //! Save every ray steps between the vertex (source) and an antenna (target), unless DATA_SAVE_MODE is 2. 02-12-2021 -MK-
                            //! These xz coordinates were calculated after we convert the earth coordinates to flat coordinates by the RaySolver::Earth_to_Flat_same_angle()
                            stations[i].strings[j].antennas[k].ray_step.resize(ray_sol_cnt + 1);    ///< resize by number of ray solutions
                            stations[i].strings[j].antennas[k].ray_step[ray_sol_cnt].resize(2); ///< resize by xz values
                            for (int steps = 0; steps < (int) RayStep[ray_sol_cnt][0].size(); steps++)
                            {
                                ///< push back each ray step coordinates
                                stations[i].strings[j].antennas[k].ray_step[ray_sol_cnt][0].push_back(RayStep[ray_sol_cnt][0][steps]);
                                stations[i].strings[j].antennas[k].ray_step[ray_sol_cnt][1].push_back(RayStep[ray_sol_cnt][1][steps]);
                            }

                            // get ice attenuation factor

                            double IceAttenFactor = 1.;
                            if (settings1->USE_ARA_ICEATTENU == 1)
                            {
                                // use new ARA measured ice attenuation values

                                double dx, dz, dl;
                                for (int steps = 1; steps < (int) RayStep[ray_sol_cnt][0].size(); steps++)
                                {

                                    dx = RayStep[ray_sol_cnt][0][steps - 1] - RayStep[ray_sol_cnt][0][steps];
                                    dz = RayStep[ray_sol_cnt][1][steps - 1] - RayStep[ray_sol_cnt][1][steps];
                                    dl = sqrt((dx *dx) + (dz *dz));

                                    // Skipping attenuation calculation when the distance between two RaySteps is 0. Prevening adds -nan into the IceAttenFactor. (MK 2021)
                                    if (dl > 0)
                                    {
                                        // use new ice model
                                        // use the midpoint of the array to calculate the attenuation length, instead of the end of the ray (BAC 2020)
                                        IceAttenFactor *= (exp(-dl / icemodel->GetARAIceAttenuLength(-RayStep[ray_sol_cnt][1][steps])) + exp(-dl / icemodel->GetARAIceAttenuLength(-RayStep[ray_sol_cnt][1][steps - 1]))) / 2;
                                    }
                                }
                            }
                            else if (settings1->USE_ARA_ICEATTENU == 0)
                            {
                                // use old method
                                IceAttenFactor = exp(-ray_output[0][ray_sol_cnt] / icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0));
                            }

                            if (debugmode == 0)
                            {

                                // set viewangle, launch_vector, receive vectors
                                viewangle = ray_output[1][ray_sol_cnt];
                                GetParameters(event->Nu_Interaction[0].posnu,   // posnu
                                    detector->stations[i].strings[j].antennas[k],   // trg antenna
                                    event->Nu_Interaction[0].nnu,   // nnu
                                    viewangle,  // inputs launch_angle, returns viewangle
                                    ray_output[2][ray_sol_cnt], // receive_angle
                                    launch_vector, receive_vector,
                                    n_trg_slappy, n_trg_pokey);

                                // check viewangle that if ray in near Cherenkov cone
                                if (viewangle * DEGRAD > 55. && viewangle * DEGRAD < 57.)
                                {
                                    // if viewangle is 56 deg +- 1 deg
                                    //cout<<"near cone! view angle : "<<viewangle * DEGRAD <<"  station["<<i<<"].string["<<j<<"].antenna["<<k<<"] with  ray_sol_cnt : "<<ray_sol_cnt<<endl;
                                }

                                // store information to report
                                stations[i].strings[j].antennas[k].view_ang.push_back(viewangle);
                                stations[i].strings[j].antennas[k].launch_ang.push_back(ray_output[1][ray_sol_cnt]);
                                stations[i].strings[j].antennas[k].rec_ang.push_back(ray_output[2][ray_sol_cnt]);
                                stations[i].strings[j].antennas[k].Dist.push_back(ray_output[0][ray_sol_cnt]);
                                stations[i].strings[j].antennas[k].L_att.push_back(icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0));
                                stations[i].strings[j].antennas[k].reflect_ang.push_back(ray_output[3][ray_sol_cnt]);
                                stations[i].strings[j].antennas[k].vmmhz.resize(ray_sol_cnt + 1);
                                stations[i].strings[j].antennas[k].Heff.resize(ray_sol_cnt + 1);
                                stations[i].strings[j].antennas[k].Vm_zoom.resize(ray_sol_cnt + 1);
                                stations[i].strings[j].antennas[k].Vm_zoom_T.resize(ray_sol_cnt + 1);
                                stations[i].strings[j].antennas[k].Vfft.resize(ray_sol_cnt + 1);
                                stations[i].strings[j].antennas[k].Vfft_noise.resize(ray_sol_cnt + 1);
                                stations[i].strings[j].antennas[k].V.resize(ray_sol_cnt + 1);
                                stations[i].strings[j].antennas[k].SignalExt.resize(ray_sol_cnt + 1);

                                // calculate the polarization vector at the source
                                Pol_vector = GetPolarization(event->Nu_Interaction[0].nnu, launch_vector);

				Vector Pol_vector_src = Pol_vector; //store the src Pol			

                                icemodel->GetFresnel(ray_output[1][ray_sol_cnt],    // launch_angle
                                    ray_output[2][ray_sol_cnt], // rec_angle
                                    ray_output[3][ray_sol_cnt], // reflect_angle
                                    event->Nu_Interaction[0].posnu,
                                    launch_vector,
                                    receive_vector,
                                    settings1,
                                    fresnel,
                                    Pol_vector);    // input src Pol and return Pol at trg
                                icemodel->GetMag(
                                    mag, 
                                    ray_output[0][ray_sol_cnt], // ray path length
                                    ray_output[1][ray_sol_cnt], // zenith angle of ray at launch
                                    ray_output[2][ray_sol_cnt], // zenith angle of ray upon receipt
                                    ray_sol_cnt,
                                    event->Nu_Interaction[0].posnu, // Neutrino
                                    detector->stations[i].strings[j].antennas[k], // Antenna
                                    -0.01, // 1cm antenna shift, inspired from NuRadioMC
                                    icemodel, settings1
                                );

                                if (ray_output[3][ray_sol_cnt] < PI / 2.)
                                {
                                    // when not reflected at the surface, angle = 100
                                    stations[i].strings[j].antennas[k].reflection.push_back(1); // say this is reflected ray
                                }
                                else
                                {
                                    stations[i].strings[j].antennas[k].reflection.push_back(0); // say this is not reflected ray
                                }

                                stations[i].strings[j].antennas[k].Pol_vector.push_back(Pol_vector);    // this Pol_vector is for the target antenna
                                stations[i].strings[j].antennas[k].Mag.push_back(mag);  // magnification factor
                                stations[i].strings[j].antennas[k].Fresnel.push_back(fresnel);  // Fresnel factor

                                vmmhz1m_sum = 0;

                                // get the arrival angle at the antenna, and store the relevant polarization factors
                                GetAngleAnt(receive_vector, detector->stations[i].strings[j].antennas[k], antenna_theta, antenna_phi);  // get theta, phi for signal ray arrived at antenna
                                GetAngleLaunch(launch_vector, launch_theta, launch_phi);

                                Vector thetaHat = Vector(cos(antenna_theta *(PI / 180)) *cos(antenna_phi *(PI / 180)),
                                    cos(antenna_theta *(PI / 180)) *sin(antenna_phi *(PI / 180)),
                                    -sin(antenna_theta *(PI / 180)));

                                Vector phiHat = Vector(-sin(antenna_phi *(PI / 180)),
                                    cos(antenna_phi *(PI / 180)),
                                    0);
                                stations[i].strings[j].antennas[k].Pol_factorH.push_back(abs(phiHat *Pol_vector));
                                stations[i].strings[j].antennas[k].Pol_factorV.push_back(abs(thetaHat *Pol_vector));
                                stations[i].strings[j].antennas[k].phi_rec.push_back(antenna_phi *(PI / 180));
                                stations[i].strings[j].antennas[k].theta_rec.push_back(antenna_theta *(PI / 180));
                                stations[i].strings[j].antennas[k].phi_launch.push_back(launch_phi *(PI / 180));
                                stations[i].strings[j].antennas[k].theta_launch.push_back(launch_theta *(PI / 180));                                

                                // old freq domain signal mode (AVZ model)
                                if (settings1->SIMULATION_MODE == 0)
                                {

                                    // initially give raysol has actual signal
                                    stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] = 1;

                                    double vmmhz_filter[(int)(detector->GetFreqBin())];

                                    for (int l = 0; l < detector->GetFreqBin(); l++)
                                    {
                                        // for detector freq bin numbers

                                        //cout<<"TaperVmMHz inputs VA:"<<viewangle<<" th_em:"<<event->Nu_Interaction[0].d_theta_em[l]<<" th_had:"<<event->Nu_Interaction[0].d_theta_had[l]<<" emfrac:"<<event->Nu_Interaction[0].emfrac<<" hadfrac:"<<event->Nu_Interaction[0].hadfrac<<" vmmhz1m:"<<event->Nu_Interaction[0].vmmhz1m[l]<<endl;
                                        //                                       switch (event->IsCalpulser){
                                        //                                           case 0:
                                        if (event->IsCalpulser > 0)
                                        {
                                            vmmhz1m_tmp = event->Nu_Interaction[0].vmmhz1m[l] *settings1->CALPUL_AMP;
                                            //vmmhz1m_tmp = event->Nu_Interaction[0].vmmhz1m[l];// calpulser -> let's use slight offset from cone
                                            //signal->TaperVmMHz(settings1->CALPUL_OFFCONE_ANGLE*RADDEG, event->Nu_Interaction[0].d_theta_em[l], event->Nu_Interaction[0].d_theta_had[l], event->Nu_Interaction[0].emfrac, event->Nu_Interaction[0].hadfrac, vmmhz1m_tmp, vmmhz1m_em);
                                        }
                                        else
                                        {
                                            vmmhz1m_tmp = event->Nu_Interaction[0].vmmhz1m[l];
                                            signal->TaperVmMHz(viewangle, event->Nu_Interaction[0].d_theta_em[l], event->Nu_Interaction[0].d_theta_had[l], event->Nu_Interaction[0].emfrac, event->Nu_Interaction[0].hadfrac, vmmhz1m_tmp, vmmhz1m_em);
                                        }

                                        //                                               break;
                                        //                                           case 1:
                                        //                                               vmmhz1m_tmp = 0;
                                        //                                               break;
                                        //                                           case 2:
                                        //                                               vmmhz1m_tmp = 0;
                                        //                                               break;
                                        //                                           default:
                                        //                                               vmmhz1m_tmp = event->Nu_Interaction[0].vmmhz1m[l];
                                        //                                               break;
                                        //                                       }

                                        //                                       signal->TaperVmMHz(viewangle, event->Nu_Interaction[0].d_theta_em[l], event->Nu_Interaction[0].d_theta_had[l], event->Nu_Interaction[0].emfrac, event->Nu_Interaction[0].hadfrac, vmmhz1m_tmp, vmmhz1m_em);
                                        //cout<<"TaperVmMHz (1m at view angle) at "<<l<<"th bin : "<<vmmhz1m_tmp<<endl;

                                        // multiply all factors for traveling ice
                                        //vmmhz1m_tmp = vmmhz1m_tmp / ray_output[0][ray_sol_cnt] *exp(-ray_output[0][ray_sol_cnt]/icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0)) *mag * fresnel;   // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
                                        if (settings1->USE_ARA_ICEATTENU == 1 || settings1->USE_ARA_ICEATTENU == 0)
                                        {
                                            vmmhz1m_tmp = vmmhz1m_tmp / ray_output[0][ray_sol_cnt] *IceAttenFactor *mag * fresnel;  // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
                                        }
                                        else if (settings1->USE_ARA_ICEATTENU == 2)
                                        {
                                            double IceAttenFactor = 1.;
                                            double dx, dz, dl;
                                            for (int steps = 1; steps < (int) RayStep[ray_sol_cnt][0].size(); steps++)
                                            {
                                                dx = RayStep[ray_sol_cnt][0][steps - 1] - RayStep[ray_sol_cnt][0][steps];
                                                dz = RayStep[ray_sol_cnt][1][steps - 1] - RayStep[ray_sol_cnt][1][steps];
                                                dl = sqrt((dx *dx) + (dz *dz));

                                                // Skipping attenuation calculation when the distance between two RaySteps is 0. Prevening adds -nan into the IceAttenFactor. (MK 2021)
                                                if (dl > 0)
                                                {
                                                    // IceAttenFactor *=  exp(-dl / icemodel->GetFreqDepIceAttenuLength(-RayStep[ray_sol_cnt][1][steps], detector->GetFreq(l) / 1e9));
                                                    // use ray midpoint for attenuation calculation
                                                    IceAttenFactor *= (exp(-dl / icemodel->GetFreqDepIceAttenuLength(-RayStep[ray_sol_cnt][1][steps], detector->GetFreq(l) / 1e9)) +
                                                        exp(-dl / icemodel->GetFreqDepIceAttenuLength(-RayStep[ray_sol_cnt][1][steps - 1], detector->GetFreq(l) / 1e9))
                                                   ) / 2.;  // 1e9 to convert to GHz
                                                }
                                            }
                                            vmmhz1m_tmp = vmmhz1m_tmp / ray_output[0][ray_sol_cnt] *IceAttenFactor *mag * fresnel;  // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
                                        }
                                        //cout<<"AttenLength : "<<icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0)<<endl;

                                        vmmhz1m_sum += vmmhz1m_tmp;

                                        stations[i].strings[j].antennas[k].vmmhz[ray_sol_cnt].push_back(vmmhz1m_tmp);

                                        freq_tmp = detector->GetFreq(l);    // freq in Hz

                                        //cout << "Check 1" << endl;
                                        /*
                                        // Get ant gain with 2-D interpolation (may have bug?) 
                                        //
                                        heff = GaintoHeight(detector->stations[i].strings[j].antennas[k].GetG(detector, freq_tmp*1.E-6, // to MHz
                                                    antenna_theta, antenna_phi), 
                                                    freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]));
                                         */
                                        heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta,  antenna_phi,  detector->stations[i].strings[j].antennas[k].type, j, k),
                                                            freq_tmp, 
                                                            icemodel->GetN(detector->stations[i].strings[j].antennas[k]),
                                                            detector->GetImpedance(freq_tmp*1.E-6, detector->stations[i].strings[j].antennas[k].type, j, k));                                        

                                        //cout<<"n_medium : "<<icemodel->GetN(detector->stations[i].strings[j].antennas[k])<<endl;
                                        //cout<<"gain : "<<detector->stations[i].strings[j].antennas[k].GetG(detector, freq_tmp*1.E-6, antenna_theta, antenna_phi)<<endl;
                                        //cout<<"heff : "<<heff<<endl;
                                        stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(heff);

                                        // apply pol factor, heff
                                        if (event->IsCalpulser == 1)
                                        {
                                            //cout<<"set signal pol as Hpol for Calpulser1 evts"<<endl;
                                            Pol_vector = n_trg_slappy;
                                        }
                                        else if (event->IsCalpulser == 2)
                                        {
                                            //cout<<"set signal pol as Vpol for Calpulser2 evts"<<endl;
                                            Pol_vector = n_trg_pokey;
                                        }
                                        else if (event->IsCalpulser == 3)
                                        {
                                            //cout<<"set signal pol as Hpol for Calpulser2 evts"<<endl;
                                            Pol_vector = n_trg_slappy;
                                        }
                                        else if (event->IsCalpulser == 4)
                                        {
                                            //cout<<"set signal pol as Vpol + Hpol for Calpulser2 evts"<<endl;
                                            Pol_vector = n_trg_slappy + n_trg_pokey;
                                        }

                                        ApplyAntFactors(heff, n_trg_pokey, n_trg_slappy, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, vmmhz1m_tmp, antenna_theta, antenna_phi);

                                        //cout << "Check 2" << endl;
                                        //stations[i].strings[j].antennas[k].VHz_antfactor[ray_sol_cnt].push_back(vmmhz1m_tmp);

                                        // apply filter
                                        ApplyFilter(l, detector, vmmhz1m_tmp);

                                        // apply Preamp gain
                                        ApplyPreamp(l, detector, vmmhz1m_tmp);

                                        // apply FOAM gain
                                        ApplyFOAM(l, detector, vmmhz1m_tmp);

                                        //stations[i].strings[j].antennas[k].VHz_filter[ray_sol_cnt].push_back(vmmhz1m_tmp);
                                        vmmhz_filter[l] = vmmhz1m_tmp;
                                    }   // end for freq bin

                                    stations[i].strings[j].antennas[k].Pol_factor.push_back(Pol_factor);

                                    //cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"].vmmhz1m["<<ray_sol_cnt<<"][0] : "<<stations[i].strings[j].antennas[k].vmmhz[ray_sol_cnt][0]<<endl;

                                    //MakeArraysforFFT(settings1, detector, i, stations[i].strings[j].antennas[k].VHz_filter[ray_sol_cnt], volts_forfft);
                                    MakeArraysforFFT(settings1, detector, i, vmmhz_filter, volts_forfft);

                                    // save freq domain array which is prepaired for realft
                                    for (int n = 0; n < settings1->NFOUR / 2; n++)
                                    {
                                        stations[i].strings[j].antennas[k].Vfft[ray_sol_cnt].push_back(volts_forfft[n]);
                                    }

                                    // now, after realft, volts_forfft is time domain signal at backend of antenna
                                    Tools::realft(volts_forfft, -1, settings1->NFOUR / 2);
                                    //Tools::realft(volts_forfft,1,settings1->NFOUR/2);

                                    //cout<<"Finished getting V signal part!!"<<endl;

                                    stations[i].strings[j].antennas[k].PeakV.push_back(FindPeak(volts_forfft, settings1->NFOUR / 2));

                                    // Vfft_noise_org is in fft freq bin!!
                                    // same unit with Vfft[V] but filter not applied

                                    Tools::NormalTimeOrdering(settings1->NFOUR / 2, volts_forfft);
                                    //cout<<"finished NormalTimeOrdering!!"<<endl;

                                    for (int n = 0; n < settings1->NFOUR / 2; n++)
                                    {

                                        if (settings1->TRIG_ANALYSIS_MODE != 2)
                                        {
                                            // not pure noise mode (we need signal)
                                            stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(volts_forfft[n]);
                                        }
                                        else if (settings1->TRIG_ANALYSIS_MODE == 2)
                                        {
                                            // pure noise mode (set signal to 0)
                                            stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(0.);
                                        }

                                        //stations[i].strings[j].antennas[k].time[ray_sol_cnt].push_back(stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt] + (double)(n - settings1->NFOUR/4)*settings1->TIMESTEP);      // time at 0 s is when ray started at the posnu
                                    }
                                    //cout<<"finished push_back V, V_noise V_total, and time!!"<<endl;
                                }   // if SIMULATION_MODE = 0

                                else if (settings1->SIMULATION_MODE == 1)
                                {

                                    // if event is not calpulser
                                    if (event->IsCalpulser == 0)
                                    {

                                        if (settings1->EVENT_TYPE == 0)
                                        {
                                            // see if integrated shower profile LQ is non-zero
                                            // and near the cone viewangle
                                            if (event->Nu_Interaction[0].LQ > 0 && (fabs(viewangle - signal->CHANGLE_ICE) <= settings1->OFFCONE_LIMIT *RADDEG))
                                            {
                                                // initially give raysol has actual signal
                                                stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] = 1;

                                                // let's make NFOUR/2 bin of time domain pure signal part for now
                                                // later once we understand how to apply antenna phase, total electronics with phase, apply those

                                                double atten_factor = 0.;
                                                if (settings1->USE_ARA_ICEATTENU == 1 || settings1->USE_ARA_ICEATTENU == 0)
                                                {
                                                    atten_factor = 1. / ray_output[0][ray_sol_cnt] *IceAttenFactor *mag * fresnel;  // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
                                                }
                                                else if (settings1->USE_ARA_ICEATTENU == 2)
                                                {
                                                    atten_factor = 1. / ray_output[0][ray_sol_cnt] *mag * fresnel;  //apply freq dependent IceAttenFactor later
                                                }

                                                // signal before the antenna (get signal at 1m and apply atten factor)
                                                signal->GetVm_FarField_Tarray(event, settings1, viewangle, atten_factor, outbin, Tarray, Earray, stations[i].strings[j].antennas[k].skip_bins[ray_sol_cnt]);

                                                dT_forfft = Tarray[1] - Tarray[0];  // step in ns

                                                int Ntmp = settings1->TIMESTEP *1.e9 / dT_forfft;
                                                stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = 1;
                                                while (Ntmp > 1)
                                                {
                                                    Ntmp = Ntmp / 2;
                                                    stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] *2;
                                                }
                                                stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] *settings1->NFOUR / 2;
                                                // now new NFOUR for zero padding

                                                // now we have to make NFOUR/2 number of bins with random init time
                                                //
                                                // as a test, make first as it is and zero pad

                                                double V_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];
                                                double T_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];

						int max_bire_ray_cnt = settings1->BIREFRINGENCE + 1; // rays in birefringence per ray solution

                                                double V_forfft_bire[max_bire_ray_cnt][stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]]; // for the waveforms of the rays in birefringence
	
						for ( int bire_ray_cnt = 0; bire_ray_cnt < max_bire_ray_cnt; bire_ray_cnt++ )
						{
		
							max_bire_ray_cnt = birefringence->Reflected_ray_remove_bire(ray_output[3][ray_sol_cnt], max_bire_ray_cnt); //change to 1 if the ray solution is reflected

                                                	for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]; n++)
                                                	{

                                                    		if (n < outbin)
                                                    		{
                                                        		stations[i].strings[j].antennas[k].Vm_zoom[ray_sol_cnt].push_back(Earray[n]);
                                                        		stations[i].strings[j].antennas[k].Vm_zoom_T[ray_sol_cnt].push_back(Tarray[n]);
                                                    		}

                                                    		// make Tarray, Earray located at the center of Nnew array

                                                    		T_forfft[n] = Tarray[outbin / 2] - (dT_forfft *(double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - n));

                                                    		if ((n >= stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - outbin / 2) &&
                                                        	(n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 + outbin / 2))
                                                    		{
                                                        		V_forfft[n] = Earray[n - (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - outbin / 2)];
                                                    		}
                                                    		else
                                                        		V_forfft[n] = 0.;
                                                	}

                                                	// just get peak from the array
                                                	stations[i].strings[j].antennas[k].PeakV.push_back(FindPeak(Earray, outbin));

							int T_shift_bire = int(time_diff_birefringence/dT_forfft); //time shift for birefringence
							double split_factor_bire = birefringence->Power_split_factor(Pol_vector_src, bire_ray_cnt, ray_output[3][ray_sol_cnt], settings1); //split power factor for birefringence
							birefringence->Time_shift_and_power_split(V_forfft, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], T_shift_bire, split_factor_bire, bire_ray_cnt, max_bire_ray_cnt, settings1); // apply time differences and power split

                                                	// this forward fft volts_forfft is now in unit of V at each freq we can just apply each bin's gain factor to each freq bins
                                                	// without any phase consideration,
                                                	// apply same factor to both real, img parts

                                                	// get spectrum with zero padded WF
                                                	Tools::realft(V_forfft, 1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

                                                	dF_Nnew = 1. / ((double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]) *(dT_forfft) *1.e-9);    // in Hz

                                                	freq_tmp = dF_Nnew *((double) stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

                                                	freq_lastbin = freq_tmp;

							birefringence->Principal_axes_polarization(Pol_vector, bire_ray_cnt, max_bire_ray_cnt, settings1); //For birefringence, modify the polarization at the antennas

                                                	/*
                                                	// Get ant gain with 2-D interpolation 
                                                 	*/
                            
                                                    heff_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type, j, k),
                                                                                freq_tmp, 
                                                                                icemodel->GetN(detector->stations[i].strings[j].antennas[k]),
                                                                                detector->GetImpedance(freq_tmp*1.E-6, detector->stations[i].strings[j].antennas[k].type, k));                              


                                                	for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2; n++)
                                                	{
                                                    		freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq
                                                        
                                                            heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta,  antenna_phi,  detector->stations[i].strings[j].antennas[k].type, j, k),
                                                                                freq_tmp, 
                                                                                icemodel->GetN(detector->stations[i].strings[j].antennas[k]),
                                                                                detector->GetImpedance(freq_tmp*1.E-6, detector->stations[i].strings[j].antennas[k].type, k));                                                        

                                                    		stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(heff);

                                                    		//apply freq dependent attenuation model
                                                    		if (settings1->USE_ARA_ICEATTENU == 2)
                                                    		{
                                                        		double IceAttenFactor = 1.;
                                                        		double dx, dz, dl;
                                                        		for (int steps = 1; steps < (int) RayStep[ray_sol_cnt][0].size(); steps++)
                                                        		{
                                                            			dx = RayStep[ray_sol_cnt][0][steps - 1] - RayStep[ray_sol_cnt][0][steps];
                                                            			dz = RayStep[ray_sol_cnt][1][steps - 1] - RayStep[ray_sol_cnt][1][steps];
                                                            			dl = sqrt((dx *dx) + (dz *dz));

                                                            			// Skipping attenuation calculation when the distance between two RaySteps is 0. Prevening adds -nan into the IceAttenFactor. (MK 2021)
                                                            			if (dl > 0)
                                                            			{
                                                                			// use ray midpoint for attenuation calculation
                                                                			IceAttenFactor *= (exp(-dl / icemodel->GetFreqDepIceAttenuLength(-RayStep[ray_sol_cnt][1][steps], freq_tmp *1.E-9)) +
                                                                			exp(-dl / icemodel->GetFreqDepIceAttenuLength(-RayStep[ray_sol_cnt][1][steps - 1], freq_tmp *1.E-9))
                                                               				) / 2.;  // 1e9 for conversion to GHz
                                                            			}
                                                        		}
                                                        		//cout << "apply IceAttenFactor to the real part of fft. V_forfft[2 *n] = " << V_forfft[2 *n] << " *" << IceAttenFactor << endl;
                                                        		V_forfft[2 *n] *= IceAttenFactor;   // apply IceAttenFactor to the real part of fft
                                                        		V_forfft[2 *n + 1] *= IceAttenFactor;   // apply IceAttenFactor to the imag part of fft
                                                    		}

                                                    		// apply ant factors
                                                    		if (n > 0)
                                                    		{
                                                                
                                                                ApplyAntFactors_Tdomain(detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type),
                                                                heff, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], settings1, antenna_theta, antenna_phi, freq_tmp);                                                                
                                                    		}
                                                    		else
                                                    		{
                                                        		ApplyAntFactors_Tdomain_FirstTwo(heff, heff_lastbin, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi, freq_tmp);
                                                    		}

                                                    		// apply entire elect chain gain, phase
                                                    		if (n > 0)
                                                    		{
                                                        		ApplyElect_Tdomain(freq_tmp *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings1);
						    		}
                                                    		else
                                                    		{
                                                        		ApplyElect_Tdomain_FirstTwo(freq_tmp *1.e-6, freq_lastbin *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings1);
                                                    		}
                                                	}   // end for freq bin

                                                	// now get time domain waveform back by inv fft
                                                	Tools::realft(V_forfft, -1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);
			
							birefringence->Store_V_forfft_for_interference(V_forfft, V_forfft_bire[bire_ray_cnt], stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]); //Store waveforms from birefringence for interference

						} // end for bire_ray_cnt 

						birefringence->Two_rays_interference(V_forfft, V_forfft_bire[0], V_forfft_bire[1], stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], max_bire_ray_cnt, settings1); //Apply interference of two rays from birefringence

                                                // do linear interpolation
                                                // changed to sinc interpolation Dec 2020 by BAC
                                                Tools::SincInterpolation(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], T_forfft, V_forfft, settings1->NFOUR / 2, T_forint, volts_forint);

                                                stations[i].strings[j].antennas[k].Pol_factor.push_back(Pol_factor);
                                                for (int n = 0; n < settings1->NFOUR / 2; n++)
                                                {

                                                    if (settings1->TRIG_ANALYSIS_MODE != 2)
                                                    {
                                                        // not pure noise mode (we need signal)
                                                        stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(volts_forint[n] *2. / (double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]));  // 2/N for inverse FFT normalization factor
                                                    }
                                                    else if (settings1->TRIG_ANALYSIS_MODE == 2)
                                                    {
                                                        // pure noise mode (set signal to 0)
                                                        stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(0.);
                                                    }
                                                }
                                            }

                                            else
                                            {
                                                // no signal generating

                                                // initially give raysol has actual signal
                                                stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] = 0;

                                                // if no signal, push_back 0 values (otherwise the value inside will remain as old value)
                                                for (int n = 0; n < settings1->NFOUR / 2; n++)
                                                {

                                                    if (n < outbin)
                                                    {
                                                        stations[i].strings[j].antennas[k].Vm_zoom[ray_sol_cnt].push_back(0.);
                                                        stations[i].strings[j].antennas[k].Vm_zoom_T[ray_sol_cnt].push_back(n);
                                                    }
                                                    stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(0.);
                                                }

                                                stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = settings1->NFOUR / 2;
                                                // just get peak from the array
                                                stations[i].strings[j].antennas[k].PeakV.push_back(0.);
                                            }
                                            
//                                             // @Justin:  Make this dynamic for user to set polarization and check if it does ray-tracing.
                                            
//                                             double psi = TMath::DegToRad()*settings1->CLOCK_ANGLE;
//                                             double theta = acos(receive_vector[2]); //receive_vector is a unit vector
//                                             double phi = atan2(receive_vector[1],receive_vector[0]);
                                                                             
//                                             //Justin's method
//                                             double newPol_vectorX = -cos(psi)*cos(theta)*cos(phi) + sin(psi)*sin(phi);
//                                             double newPol_vectorY = -cos(psi)*cos(theta)*sin(phi) - sin(psi)*cos(phi);
//                                             double newPol_vectorZ = cos(psi)*sin(theta);
                                            
//                                             Vector Pol_vector = Vector(newPol_vectorX, newPol_vectorY, newPol_vectorZ);
//                                             //Justin's Method
                                            
                                        } // neutrino events
                                        else if (settings1->EVENT_TYPE == 10)
                                        {

                                            /*
                                            NB: This has not be "cleaned up" in the way that the neutrino mode (EVENT_TYPE==10) has above
                                            (BAC June 2022)
                                            */

                                            // initially give raysol has actual signal
                                            stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] = 1;

                                            int waveform_bin = (int) signal->ArbitraryWaveform_V.size();
                                            //                  cout << waveform_bin << endl;

                                            //dT_forfft = Tarray[1] - Tarray[0];    // step in ns
                                            //dT_forfft = detector->CalPulserWF_ns[1] - detector->CalPulserWF_ns[0];    // step in ns
                                            dT_forfft = signal->ArbitraryWaveform_T[1] - signal->ArbitraryWaveform_T[0];    // step in ns

                                            //                  cout << "dT_forfft: " << dT_forfft << endl;
                                            int Ntmp = settings1->TIMESTEP *1.e9 / dT_forfft;
                                            //                  cout << Ntmp << endl;

                                            stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = 1;
                                            while (Ntmp > 1)
                                            {
                                                Ntmp = Ntmp / 2;
                                                stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] *2;
                                            }
                                            stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] *settings1->NFOUR / 2;
                                            // now new NFOUR for zero padding

                                            // now we have to make NFOUR/2 number of bins with random init time
                                            //
                                            // as a test, make first as it is and zero pad

                                            double V_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];
                                            double T_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];

                                            for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]; n++)
                                            {
                                                //cout << n << endl;
                                                if (n < waveform_bin)
                                                {
                                                    stations[i].strings[j].antennas[k].Vm_zoom[ray_sol_cnt].push_back(signal->ArbitraryWaveform_V[n]);
                                                    stations[i].strings[j].antennas[k].Vm_zoom_T[ray_sol_cnt].push_back(signal->ArbitraryWaveform_T[n]);
                                                }

                                                // make Tarray, Earray located at the center of Nnew array

                                                T_forfft[n] = signal->ArbitraryWaveform_T[waveform_bin / 2] - (dT_forfft *(double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - n));

                                                if ((n >= stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - waveform_bin / 2) &&
                                                    (n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 + waveform_bin / 2))
                                                {
                                                    V_forfft[n] = signal->ArbitraryWaveform_V[n - (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - waveform_bin / 2)];
                                                }
                                                else
                                                    V_forfft[n] = 0.;
                                            }

                                            // just get peak from the array
                                            stations[i].strings[j].antennas[k].PeakV.push_back(FindPeak(V_forfft, waveform_bin));

                                            // get spectrum with zero padded WF
                                            Tools::realft(V_forfft, 1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

                                            dF_Nnew = 1. / ((double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]) *(dT_forfft) *1.e-9);    // in Hz

                                            freq_tmp = dF_Nnew *((double) stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

                                            freq_lastbin = freq_tmp;

                                            Pol_vector = n_trg_pokey;

                                            //
                                            //
                                            for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2; n++)
                                            {

                                                freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

                                                heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta,  antenna_phi,  detector->stations[i].strings[j].antennas[k].type, j, k),
                                                                    freq_tmp, 
                                                                    icemodel->GetN(detector->stations[i].strings[j].antennas[k]),
                                                                    detector->GetImpedance(freq_tmp*1.E-6, detector->stations[i].strings[j].antennas[k].type, k));

                                                stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(heff);

                                                if (n > 0)
                                                {
                                                    ApplyAntFactors_Tdomain(detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type),
                                                        heff, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], settings1, antenna_theta, antenna_phi, freq_tmp);
                                                }
                                                else
                                                {
                                                    ApplyAntFactors_Tdomain_FirstTwo(heff, heff_lastbin, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi, freq_tmp);
                                                }

                                                //
                                                // apply entire elect chain gain, phase
                                                //
                                                if (n > 0)
                                                {
                                                    ApplyElect_Tdomain(freq_tmp *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings1);
						}
                                                else
                                                {
                                                    ApplyElect_Tdomain_FirstTwo(freq_tmp *1.e-6, freq_lastbin *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings1);
                                                }
                                            }   // end for freq bin

                                            // now get time domain waveform back by inv fft
                                            Tools::realft(V_forfft, -1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

                                            // we need to do normal time ordering as we did zero padding(?)
                                            // If signal is located at the center, we don't need to do NormalTimeOrdering???
                                            //Tools::NormalTimeOrdering(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], V_forfft);

                                            // skip linear interpolation for now
                                            // changed to sinc interpolation Dec 2020 by BAC
                                            Tools::SincInterpolation(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], T_forfft, V_forfft, settings1->NFOUR / 2, T_forint, volts_forint);

                                            // check what we save as V[], volts_forint? or volts_forfft

                                            for (int n = 0; n < settings1->NFOUR / 2; n++)
                                            {

                                                if (settings1->TRIG_ANALYSIS_MODE != 2)
                                                {
                                                    // not pure noise mode (we need signal)
                                                    stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(settings1->ARBITRARY_EVENT_ATTENUATION *volts_forint[n] *2. / (double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]));  // 2/N for inverse FFT normalization factor
                                                }
                                                else if (settings1->TRIG_ANALYSIS_MODE == 2)
                                                {
                                                    // pure noise mode (set signal to 0)
                                                    stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(0.);
                                                }
                                            }
                                        }   // Arbitrary Events
                                        
                                        //Attempting to create pulser events using arbitrary events as a framework. - JCF 4/6/2023
                                        else if (settings1->EVENT_TYPE == 11)
                                        {

                                            /*
                                            NB: This has not be "cleaned up" in the way that the neutrino mode (EVENT_TYPE==10) has above
                                            (BAC June 2022)
                                            */

                                            // initially give raysol has actual signal
                                            stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] = 1;

                                            int waveform_bin = (int) signal->PulserWaveform_V.size();
                                            //                  cout << waveform_bin << endl;

                                            //dT_forfft = Tarray[1] - Tarray[0];    // step in ns
                                            //dT_forfft = detector->CalPulserWF_ns[1] - detector->CalPulserWF_ns[0];    // step in ns
                                            dT_forfft = signal->PulserWaveform_T[1] - signal->PulserWaveform_T[0];    // step in ns

                                            // cout << "dT_forfft: " << dT_forfft << endl;
                                            int Ntmp = settings1->TIMESTEP *1.e9 / dT_forfft;
                                            // cout << "Ntmp = " << Ntmp << endl;
                                            
                                            stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = 1;
                                            // cout << "Nnew = " << stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] << endl;
                                            while (Ntmp > 1)
                                            {
                                                Ntmp = Ntmp / 2;
                                                stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] *2;
                                            }
                                            // cout << "Ntmp = " << Ntmp << endl;
                                            stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] *settings1->NFOUR / 2;
                                            // cout << "Nnew = " << stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] << endl;
                                            // now new NFOUR for zero padding

                                            // now we have to make NFOUR/2 number of bins with random init time
                                            //
                                            // as a test, make first as it is and zero pad

                                            double V_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];
                                            double T_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];
                                            
					    int max_bire_ray_cnt = settings1->BIREFRINGENCE + 1;  // rays in birefringence per ray solution
                                            double V_forfft_bire[max_bire_ray_cnt][stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]]; // for the waveforms of the rays in birefringence

                                            for ( int bire_ray_cnt = 0; bire_ray_cnt < max_bire_ray_cnt; bire_ray_cnt++ )
					    {

						max_bire_ray_cnt = birefringence->Reflected_ray_remove_bire(ray_output[3][ray_sol_cnt], max_bire_ray_cnt); //change to 1 if the ray solution is reflected					
                                            	for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]; n++)
                                            	{
                                                	//cout << n << endl;
                                                	if (n < waveform_bin)
                                                	{
                                                    		stations[i].strings[j].antennas[k].Vm_zoom[ray_sol_cnt].push_back(signal->PulserWaveform_V[n]);
                                                    		stations[i].strings[j].antennas[k].Vm_zoom_T[ray_sol_cnt].push_back(signal->PulserWaveform_T[n]);
                                                	}

                                                	// make Tarray, Earray located at the center of Nnew array

                                                	T_forfft[n] = signal->PulserWaveform_T[waveform_bin / 2] - (dT_forfft *(double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - n));

                                                	if ((n >= stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - waveform_bin / 2) &&
                                                    	   (n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 + waveform_bin / 2))
                                                	{
                                                    		V_forfft[n] = signal->PulserWaveform_V[n - (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - waveform_bin / 2)];
                                                	}
                                                	else
                                                    		V_forfft[n] = 0.;
                                            	}                                         

                                            	// just get peak from the array
                                            	// stations[i].strings[j].antennas[k].PeakV.push_back(-1.);    // just let -1.
                                            	stations[i].strings[j].antennas[k].PeakV.push_back(FindPeak(V_forfft, waveform_bin));

                                            	//Defining polarization at the source (using launch_vector)
                                            
                                            	double psi = TMath::DegToRad()*settings1->CLOCK_ANGLE;
                                            	double theta = acos(launch_vector[2]); //launch_vector is a unit vector
                                            	double phi = atan2(launch_vector[1],launch_vector[0]);
                                                                             
                                            	//Justin's method
                                            	double newPol_vectorX = -cos(psi)*cos(theta)*cos(phi) + sin(psi)*sin(phi);
                                            	double newPol_vectorY = -cos(psi)*cos(theta)*sin(phi) - sin(psi)*cos(phi);
                                            	double newPol_vectorZ = cos(psi)*sin(theta);
                                            
                                            	Vector Pol_vector = Vector(newPol_vectorX, newPol_vectorY, newPol_vectorZ);                                            


						int T_shift_bire = int(time_diff_birefringence/dT_forfft); //time shift for birefringence
                                                double split_factor_bire = birefringence->Power_split_factor(Pol_vector, bire_ray_cnt, ray_output[3][ray_sol_cnt], settings1); //split power factor for birefringence
                                                birefringence->Time_shift_and_power_split(V_forfft, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], T_shift_bire, split_factor_bire, bire_ray_cnt, max_bire_ray_cnt, settings1); // apply time differences and power split

                                            	// get spectrum with zero padded WF
                                            	Tools::realft(V_forfft, 1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);                                            

                                            	dF_Nnew = 1. / ((double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]) *(dT_forfft) *1.e-9);    // in Hz

                                            	freq_tmp = dF_Nnew *((double) stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

                                            	freq_lastbin = freq_tmp;
                                                
                                                heff_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type, j, k),
                                                                            freq_tmp, 
                                                                            icemodel->GetN(detector->stations[i].strings[j].antennas[k]),
                                                                            detector->GetImpedance(freq_tmp*1.E-6, detector->stations[i].strings[j].antennas[k].type, k));                                                                                 
                                            
                                            	icemodel->GetFresnel(
                                                    ray_output[1][ray_sol_cnt],    // launch_angle
                                                    ray_output[2][ray_sol_cnt], // rec_angle
                                                    ray_output[3][ray_sol_cnt], // reflect_angle
                                                    event->Nu_Interaction[0].posnu,
                                                    launch_vector,
                                                    receive_vector,
                                                    settings1,
                                                    fresnel,
                                                    Pol_vector);    // input src Pol and return Pol at trg
                                                icemodel->GetMag(
                                                    mag, 
                                                    ray_output[0][ray_sol_cnt], // ray path length
                                                    ray_output[1][ray_sol_cnt], // zenith angle of ray at launch
                                                    ray_output[2][ray_sol_cnt], // zenith angle of ray upon receipt
                                                    ray_sol_cnt,
                                                    event->Nu_Interaction[0].posnu, // Neutrino
                                                    detector->stations[i].strings[j].antennas[k], // Antenna
                                                    -0.01, // 1cm antenna shift, inspired from NuRadioMC
                                                    icemodel, settings1
                                                );

						birefringence->Principal_axes_polarization(Pol_vector, bire_ray_cnt, max_bire_ray_cnt, settings1); //For birefringence, modify the polarization at the antennas                                            

                                            	for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2; n++)
                                            	{

                                                	freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

                                                    heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta,  antenna_phi,  detector->stations[i].strings[j].antennas[k].type, j, k),
                                                                        freq_tmp, 
                                                                        icemodel->GetN(detector->stations[i].strings[j].antennas[k]),
                                                                        detector->GetImpedance(freq_tmp*1.E-6, detector->stations[i].strings[j].antennas[k].type, k));                                                    

                                                	stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(heff);
                                           		
                                                    if (n > 0)
                                                	{
                                                        ApplyAntFactors_Tdomain(detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type),
                                                            heff, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], settings1, antenna_theta, antenna_phi, freq_tmp);
                                                    
                                                	}
                                                	else
                                                	{
                                                    		ApplyAntFactors_Tdomain_FirstTwo(heff, heff_lastbin, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi, freq_tmp);
                                                    
                                                	}

                                                	//
                                                	// apply entire elect chain gain, phase
                                                	//
                                                	if (n > 0)
                                                	{                                                  
                                                    		ApplyElect_Tdomain(freq_tmp *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings1);                                              
                                                	}
                                                	else
                                                	{
                                                    		ApplyElect_Tdomain_FirstTwo(freq_tmp *1.e-6, freq_lastbin *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings1);
                                                	}
                                            	}   // end for freq bin

                                            	Tools::realft(V_forfft, -1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);                                            
                        
						birefringence->Store_V_forfft_for_interference(V_forfft, V_forfft_bire[bire_ray_cnt], stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]); //Store waveforms from birefringence for interference                	
				
					    } //end for bire_ray_cnt    

					    birefringence->Two_rays_interference(V_forfft, V_forfft_bire[0], V_forfft_bire[1], stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], max_bire_ray_cnt, settings1); //Apply interference of two rays from birefringence						

					    Tools::SincInterpolation(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], T_forfft, V_forfft, settings1->NFOUR / 2, T_forint, volts_forint);

                                            for (int n = 0; n < settings1->NFOUR / 2; n++)
                                            {

                                                if (settings1->TRIG_ANALYSIS_MODE != 2)
                                                {
                                                    // not pure noise mode (we need signal)
                                                    stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(settings1->ARBITRARY_EVENT_ATTENUATION *volts_forint[n] *2. / (double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]));  // 2/N for inverse FFT normalization factor

                                                    
                                                    
                                                }
                                                else if (settings1->TRIG_ANALYSIS_MODE == 2)
                                                {
                                                    // pure noise mode (set signal to 0)
                                                    stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(0.);
                                                }
                                            }
                                        } // Simple Pulser Simulation                                        
                                        
                                        
                                        //Attempting to simulate PVA pulser.  Starting separate from previous pulser event type to avoid breaking things. - JCF 1/9/2024
                                        else if (settings1->EVENT_TYPE == 12)
                                        {

                                            // Import Voltage versus time fed into antenna (via data from Alisa in IDL2_InputVoltageVersusTime.txt)
                                            stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] = 1;
                                            int waveform_bin = (int) signal->InputVoltage_V.size();
                                            dT_forfft = signal->InputVoltage_T[1] - signal->InputVoltage_T[0];    // step in ns
                                            int Ntmp = settings1->TIMESTEP *1.e9 / dT_forfft;
                                            stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = 1;
                                            while (Ntmp > 1)
                                            {
                                                Ntmp = Ntmp / 2;
                                                stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] *2;
                                            }
                                            stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] *settings1->NFOUR / 2;                                    

                                            // Pad the input voltage and take the FFT
                                            double V_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];
                                            double T_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];

					    int max_bire_ray_cnt = settings1->BIREFRINGENCE + 1;  // rays in birefringence per ray solution
                                            double V_forfft_bire[max_bire_ray_cnt][stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]]; // for the waveforms of the rays in birefringence

                                            for ( int bire_ray_cnt = 0; bire_ray_cnt < max_bire_ray_cnt; bire_ray_cnt++ )
					    {

						max_bire_ray_cnt = birefringence->Reflected_ray_remove_bire(ray_output[3][ray_sol_cnt], max_bire_ray_cnt); //change to 1 if the ray solution is reflected					
                                            	for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]; n++)
                                            	{
                                                	if (n < waveform_bin)
                                                	{
                                                    		stations[i].strings[j].antennas[k].Vm_zoom[ray_sol_cnt].push_back(signal->InputVoltage_V[n]*1e-3);  //Input voltage needs to be converted to millivolts
                                                    		stations[i].strings[j].antennas[k].Vm_zoom_T[ray_sol_cnt].push_back(signal->InputVoltage_T[n]);
                                                	}

                                                	// make Tarray, Earray located at the center of Nnew array

                                                	T_forfft[n] = signal->InputVoltage_T[waveform_bin / 2] - (dT_forfft *(double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - n));

                                                	if ((n >= stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - waveform_bin / 2) &&
                                                    	   (n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 + waveform_bin / 2))
                                                	{
                                                    		V_forfft[n] = signal->InputVoltage_V[n - (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - waveform_bin / 2)]*1e-3;  //Converting Input voltage to millivolts
                                                	}
                                                	else
                                                    		V_forfft[n] = 0.;

                                            	} 
                                                //End padding

                                            	// just get peak from the array
                                            	// stations[i].strings[j].antennas[k].PeakV.push_back(-1.);    // just let -1.
                                            	stations[i].strings[j].antennas[k].PeakV.push_back(FindPeak(V_forfft, waveform_bin));  //TODO: This needs to be the peak of the electric field at the source.

                                            	//Defining polarization at the source (using launch_vector)
                                            
                                            	// double psi = TMath::DegToRad()*settings1->CLOCK_ANGLE;
                                                double psi = 0;  //In the absence of cross-pol, the polarization angle is nominally zero for an antenna in the ice that is azimuthally symmetric. Cross pol is the next step in implementation. - JCF 2/9/2024
                                            	double theta = acos(launch_vector[2]); //launch_vector is a unit vector
                                            	double phi = atan2(launch_vector[1],launch_vector[0]);
                                                                             
                                            	//Justin's method
                                            	double newPol_vectorX = -cos(psi)*cos(theta)*cos(phi) + sin(psi)*sin(phi);
                                            	double newPol_vectorY = -cos(psi)*cos(theta)*sin(phi) - sin(psi)*cos(phi);
                                            	double newPol_vectorZ = cos(psi)*sin(theta);
                                            
                                            	Vector Pol_vector = Vector(newPol_vectorX, newPol_vectorY, newPol_vectorZ);                                            


						int T_shift_bire = int(time_diff_birefringence/dT_forfft); //time shift for birefringence
                                                double split_factor_bire = birefringence->Power_split_factor(Pol_vector, bire_ray_cnt, ray_output[3][ray_sol_cnt], settings1); //split power factor for birefringence
                                                birefringence->Time_shift_and_power_split(V_forfft, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], T_shift_bire, split_factor_bire, bire_ray_cnt, max_bire_ray_cnt, settings1); // apply time differences and power split

                                            	// Get FFT of input voltage
                                            	Tools::realft(V_forfft, 1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);                                            

                                            	dF_Nnew = 1. / ((double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]) *(dT_forfft) *1.e-9);    // in Hz

                                            	freq_tmp = dF_Nnew *((double) stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

                                            	freq_lastbin = freq_tmp;
                                                
                                                //Apply transmitting antenna effective height and impedance to calculate electric field at 1m from Tx.
                                                
                                                //Defining the launch angles using the polarization vector calculations above, converted to degrees.
                                                Tx_theta = theta*180/PI;
                                                Tx_phi = phi*180/PI;
                                                
                                                //Tx effective height for last bin.  Currently locked to standard ARA Vpol and HPol antennas.  Need to add selection mode.
                                                
                                                heff_Tx_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, Tx_theta, Tx_phi, 0, 0, 0, true),
                                                                               freq_tmp, 
                                                                               icemodel->GetN(detector->stations[i].strings[j].antennas[k]),
                                                                               detector->GetImpedance(freq_tmp*1.E-6, 0, 0, true));                                                                                      
                                                //End Tx effective height for last bin
                                                
                                                //Apply effective height of last bin for receiving antenna
                                                heff_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type, j, k),
                                                                            freq_tmp, 
                                                                            icemodel->GetN(detector->stations[i].strings[j].antennas[k]),
                                                                            detector->GetImpedance(freq_tmp*1.E-6, detector->stations[i].strings[j].antennas[k].type, k));                                                
                                                //end effective height of last bin for receiving antenna.
                                            
                                                //Apply Fresnel factors for magnification and 1/r dependence
                                            	icemodel->GetFresnel(ray_output[1][ray_sol_cnt],    // launch_angle
                                                	ray_output[2][ray_sol_cnt], // rec_angle
                                                	ray_output[3][ray_sol_cnt], // reflect_angle
                                                	event->Nu_Interaction[0].posnu,
                                                	launch_vector,
                                                	receive_vector,
                                                	settings1,
                                                	fresnel,
                                                	Pol_vector);    // input src Pol and return Pol at trg
                                                //End Fresnel application

                                                //Apply birefringence.
                                                birefringence->Principal_axes_polarization(Pol_vector, bire_ray_cnt, max_bire_ray_cnt, settings1); //For birefringence, modify the polarization at the antennas       

                                                //Begin frequency binning
                                            	for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2; n++)
                                            	{

                                                    // Calculate effective height for transmitting antenna in current frequency bin
                                                    heff_Tx = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, Tx_theta, Tx_phi, 0, 0, 0, true),
                                                                           freq_tmp, 
                                                                           icemodel->GetN(detector->stations[i].strings[j].antennas[k]),
                                                                           detector->GetImpedance(freq_tmp*1.E-6, 0, 0, true));                                                    
                                                    // End Tx effective height calculation                                                    
                                            
                                                    // Calculate effective height for receiving antenna in current frequency bin
                                                    heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type, j, k),
                                                                        freq_tmp, 
                                                                        icemodel->GetN(detector->stations[i].strings[j].antennas[k]),
                                                                        detector->GetImpedance(freq_tmp*1.E-6, detector->stations[i].strings[j].antennas[k].type, k));                                                    
                                                    // End Rx effective height calculation

                                                	stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(heff);
                                                    
                                                	freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq                                                      
                                                    
                                                    //Apply Tx antenna factors
                                                    if (n > 0)
                                                	{
                                                            ApplyAntFactors_Tdomain(detector->GetAntPhase_1D(freq_tmp *1.e-6, Tx_theta, Tx_phi, 0),
                                                                heff_Tx, Pol_vector, 0, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], settings1, Tx_theta, Tx_phi, freq_tmp, detector->GetImpedance(freq_tmp*1.E-6, 0, 0, true), true);
                                                	}
                                                	else
                                                	{
                                                    		ApplyAntFactors_Tdomain_FirstTwo(heff_Tx, heff_Tx_lastbin, Pol_vector, 0, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], Tx_theta, Tx_phi, freq_tmp);
                                                    
                                                	}
                                                    //End Tx antenna factors
                                                    
                                                	freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq                                                 
                                           		
                                                    //Apply Rx antenna factors
                                                    if (n > 0)
                                                	{
                                                            ApplyAntFactors_Tdomain(detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type),
                                                                heff, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], settings1, antenna_theta, antenna_phi, freq_tmp);
                                                	}
                                                	else
                                                	{
                                                    		ApplyAntFactors_Tdomain_FirstTwo(heff, heff_lastbin, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi, freq_tmp);
                                                    
                                                	}
                                                    //End Rx antenna factors                                                        

                                                	//
                                                	// apply entire elect chain gain, phase
                                                	//
                                                	if (n > 0)
                                                	{                                                  
                                                    		ApplyElect_Tdomain(freq_tmp *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings1);                                              
                                                	}
                                                	else
                                                	{
                                                    		ApplyElect_Tdomain_FirstTwo(freq_tmp *1.e-6, freq_lastbin *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings1);
                                                	}                                                   

                                            	}   // end for freq bin
                                                //End amplification at receiving antenna

                                            	Tools::realft(V_forfft, -1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);                                                 
                        
						birefringence->Store_V_forfft_for_interference(V_forfft, V_forfft_bire[bire_ray_cnt], stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]); //Store waveforms from birefringence for interference                	
				
					    } //end for bire_ray_cnt    

					    birefringence->Two_rays_interference(V_forfft, V_forfft_bire[0], V_forfft_bire[1], stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], max_bire_ray_cnt, settings1); //Apply interference of two rays from birefringence						

					    Tools::SincInterpolation(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], T_forfft, V_forfft, settings1->NFOUR / 2, T_forint, volts_forint);                                              

                                            for (int n = 0; n < settings1->NFOUR / 2; n++)
                                            {

                                                if (settings1->TRIG_ANALYSIS_MODE != 2)
                                                {
                                                    // not pure noise mode (we need signal)
                                                    stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(settings1->ARBITRARY_EVENT_ATTENUATION *volts_forint[n] *2. / (double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]));  // 2/N for inverse FFT normalization factor
                                                    
                                                    
                                                }
                                                else if (settings1->TRIG_ANALYSIS_MODE == 2)
                                                {
                                                    // pure noise mode (set signal to 0)
                                                    stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(0.);
                                                }
                                            }
                                        } // PVA Pulser Events
                                        
                                        
                                    }   // if not calpulser event

                                    // if calpulser event
                                    else if (event->IsCalpulser > 0)
                                    {

                                        // initially give raysol has actual signal
                                        stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] = 1;

                                        int CP_bin = (int) detector->CalPulserWF_ns.size();

                                        //dT_forfft = Tarray[1] - Tarray[0];    // step in ns
                                        dT_forfft = detector->CalPulserWF_ns[1] - detector->CalPulserWF_ns[0];  // step in ns

                                        int Ntmp = settings1->TIMESTEP *1.e9 / dT_forfft;
                                        stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = 1;
                                        while (Ntmp > 1)
                                        {
                                            Ntmp = Ntmp / 2;
                                            stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] *2;
                                        }
                                        stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] *settings1->NFOUR / 2;
                                        // now new NFOUR for zero padding

                                        // now we have to make NFOUR/2 number of bins with random init time
                                        //
                                        // as a test, make first as it is and zero pad

                                        double V_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];
                                        double T_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];

                                        for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]; n++)
                                        {

                                            if (n < CP_bin)
                                            {
                                                stations[i].strings[j].antennas[k].Vm_zoom[ray_sol_cnt].push_back(detector->CalPulserWF_V[n]);
                                                stations[i].strings[j].antennas[k].Vm_zoom_T[ray_sol_cnt].push_back(detector->CalPulserWF_ns[n]);
                                            }

                                            // make Tarray, Earray located at the center of Nnew array
                                            //T_forfft[n] = Tarray[outbin/2] - (dT_forfft*(double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]/2)) + (double)n*dT_forfft;

                                            T_forfft[n] = detector->CalPulserWF_ns[CP_bin / 2] - (dT_forfft *(double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - n));

                                            if ((n >= stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - CP_bin / 2) &&
                                                (n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 + CP_bin / 2))
                                            {
                                                V_forfft[n] = detector->CalPulserWF_V[n - (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2 - CP_bin / 2)];
                                            }
                                            else
                                                V_forfft[n] = 0.;

                                            //stations[i].strings[j].antennas[k].Vm_wo_antfactor[ray_sol_cnt].push_back(V_forfft[n]);
                                        }

                                        // just get peak from the array
                                        //stations[i].strings[j].antennas[k].PeakV.push_back(FindPeak(detector->CalPulserWF_V, CP_bin));
                                        stations[i].strings[j].antennas[k].PeakV.push_back(-1.);    // just let -1.

                                        // get spectrum with zero padded WF
                                        //Tools::realft(volts_forfft,1,settings1->NFOUR/2);
                                        Tools::realft(V_forfft, 1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

                                        //dF_outbin = 1./((double)(outbin) *(Tarray[1]-Tarray[0])*1.e-9);   // in Hz
                                        dF_Nnew = 1. / ((double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]) *(dT_forfft) *1.e-9);    // in Hz

                                        //freq_tmp = dF_Nnew*((double)stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]/2.+1.+0.5);// in Hz 0.5 to place the middle of the bin and avoid zero freq
                                        freq_tmp = dF_Nnew *((double) stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

                                        freq_lastbin = freq_tmp;

                                        // heff last bin for transmitter ant
                                        double heff_lastbin_trans;
                                        double ant_theta_trans = ray_output[1][ray_sol_cnt] *DEGRAD;    // from 0 to 180
                                        //cout<<"ant theta trans : "<<ant_theta_trans<<"deg"<<endl;
                                        heff_lastbin_trans = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6, // to MHz
                                                ant_theta_trans, antenna_phi, detector->stations[i].strings[j].antennas[k].type, j, k),
                                            freq_tmp, icemodel->GetN(event->Nu_Interaction[0].posnu));                                        

                                        // heff last bin for receiver ant
                                        /*
                                        // Get ant gain with 2-D interpolation (may have bug?) 
                                        //
                                        heff_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6,    // to MHz
                                        antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
                                        freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]));
                                        */
                                        heff_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                                                antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type, j, k),
                                            freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]));                                        

                                        // apply calpulser waveform
                                        // apply pol factor, heff
                                        if (event->IsCalpulser == 1)
                                        {
                                            //cout<<"set signal pol as Hpol for Calpulser1 evts"<<endl;
                                            Pol_vector = n_trg_slappy;
                                        }
                                        else if (event->IsCalpulser == 2)
                                        {
                                            //cout<<"set signal pol as Vpol for Calpulser2 evts"<<endl;
                                            Pol_vector = n_trg_pokey;
                                        }
                                        else if (event->IsCalpulser == 3)
                                        {
                                            //cout<<"set signal pol as Hpol for Calpulser2 evts"<<endl;
                                            Pol_vector = n_trg_slappy;
                                        }
                                        else if (event->IsCalpulser == 4)
                                        {
                                            //cout<<"set signal pol as Vpol + Hpol for Calpulser2 evts"<<endl;
                                            Pol_vector = n_trg_slappy + n_trg_pokey;
                                        }

                                        //
                                        //
                                        for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] / 2; n++)
                                        {

                                            freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

                                            //
                                            // apply ant factors (transmitter ant)
                                            //
                                            heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                                                    ant_theta_trans, antenna_phi, detector->stations[i].strings[j].antennas[k].type, j, k),
                                                freq_tmp, icemodel->GetN(event->Nu_Interaction[0].posnu));                                            
                                            //
                                            if (n > 0)
                                            {
                                                ApplyAntFactors_Tdomain(detector->GetAntPhase_1D(freq_tmp *1.e-6, ant_theta_trans, antenna_phi, detector->stations[i].strings[j].antennas[k].type),
                                                    heff, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], settings1, true, antenna_theta, antenna_phi, freq_tmp);
                                            }
                                            else
                                            {
                                                ApplyAntFactors_Tdomain_FirstTwo(heff, heff_lastbin_trans, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi, freq_tmp);
                                            }

                                            //
                                            // apply ant factors (receiver ant)
                                            //
                                            /*
                                              heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6,  // to MHz
                                              antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
                                              freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]));
                                            */
                                            heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp *1.E-6,   // to MHz
                                                    antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type, j, k),
                                                freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]));                                            

                                            stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(heff);

                                            if (n > 0)
                                            {
                                                ApplyAntFactors_Tdomain(detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type),
                                                    heff, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], settings1, antenna_theta, antenna_phi, freq_tmp);
                                            }
                                            else
                                            {
                                                ApplyAntFactors_Tdomain_FirstTwo(heff, heff_lastbin, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1], antenna_theta, antenna_phi, freq_tmp);
                                            }

                                            //
                                            // apply entire elect chain gain, phase
                                            //
                                            if (n > 0)
                                            {
                                                ApplyElect_Tdomain(freq_tmp *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings1);
                                            }
                                            else
                                            {
                                                ApplyElect_Tdomain_FirstTwo(freq_tmp *1.e-6, freq_lastbin *1.e-6, detector, V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings1);
                                            }
                                        }   // end for freq bin

                                        // now get time domain waveform back by inv fft
                                        Tools::realft(V_forfft, -1, stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

                                        // we need to do normal time ordering as we did zero padding(?)
                                        // If signal is located at the center, we don't need to do NormalTimeOrdering???
                                        //Tools::NormalTimeOrdering(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], V_forfft);

                                        // skip linear interpolation for now
                                        // changed to sinc interpolation Dec 2020 by BAC
                                        Tools::SincInterpolation(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], T_forfft, V_forfft, settings1->NFOUR / 2, T_forint, volts_forint);

                                        // check what we save as V[], volts_forint? or volts_forfft

                                        for (int n = 0; n < settings1->NFOUR / 2; n++)
                                        {

                                            if (settings1->TRIG_ANALYSIS_MODE != 2)
                                            {
                                                // not pure noise mode (we need signal)
                                                //stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(volts_forfft[n]);
                                                //stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(V_forfft[n] *2./(double)(settings1->NFOUR/2));    // 2/N for inverse FFT normalization factor
                                                //stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(volts_forint[n] *2./(double)(settings1->NFOUR/2));    // 2/N for inverse FFT normalization factor
                                                stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(volts_forint[n] *2. / (double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]));  // 2/N for inverse FFT normalization factor
                                            }
                                            else if (settings1->TRIG_ANALYSIS_MODE == 2)
                                            {
                                                // pure noise mode (set signal to 0)
                                                stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(0.);
                                            }
                                        }
                                    }   // Calpulser events

                                }   // if SIMULATION_MODE = 1

                                // check max_PeakV
                                if (max_PeakV_tmp < stations[i].strings[j].antennas[k].PeakV[ray_sol_cnt])
                                {
                                    max_PeakV_tmp = stations[i].strings[j].antennas[k].PeakV[ray_sol_cnt];
                                }
                            }   // if not debug mode

                            // check min_arrival_time
                            if (min_arrival_time_tmp > stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt])
                            {
                                min_arrival_time_tmp = stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt];
                            }

                            // check max_arrival_time
                            if (max_arrival_time_tmp < stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt])
                            {
                                max_arrival_time_tmp = stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt];
                            }

                            ray_sol_cnt++;
                        }   // end while number of solutions

                    }   // end if solution exist

                    else
                    {
                        //cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"].trg = "<<stations[i].strings[j].antennas[k].trg[ray_sol_cnt]<<"  No vmmhz1m data!"<<endl;
                    }
                }   // end if posnu selected

                else
                {
                    //cout<<" No posnu!!!!!! No signals calculated at all!!"<<endl;
                    ray_sol_cnt = 0;
                }

                stations[i].strings[j].antennas[k].ray_sol_cnt = ray_sol_cnt;   // save number of RaySolver solutions

                stations[i].Total_ray_sol += ray_sol_cnt;   // add ray_sol_cnt to Total_ray_sol

            }   // for number_of_antennas_string

        }   // for number_of_strings_station

        // set each station's min/max arrival time
        stations[i].min_arrival_time = min_arrival_time_tmp;
        stations[i].max_arrival_time = max_arrival_time_tmp;

        // set each station's max PeakV
        stations[i].max_PeakV = max_PeakV_tmp;
    }   // for number_of_stations

    // clear ray_output info
    ray_output.clear();

    // do only if it's not in debugmode
    if (debugmode == 0)
    {
        // after all values are stored in Report, set ranking of signal between antennas
        SetRank(detector);
    }

    // now loop over all antennas again to make DATA_BIN_SIZE array for signal + noise. (with time delay)
    // with that signal + noise array, we'll do convolution with diode response.
    // with the convolution result, we'll do trigger check
    for (int i = 0; i < detector->params.number_of_stations; i++)
    {

        N_pass = 0;
        N_pass_V = 0;
        N_pass_H = 0;

        stations[i].Global_Pass  = 0;
        int check_passed_global_trigger = 0;    // this switch determines if station globally triggers (in all TRIG_SCAN_MODEs)

        // If we're looking to trigger on signal, check that there 
        //   is enough signal that reached the station. 
        // Otherwise, don't worry about checking for a trigger
        int ants_with_nonzero_signal = 0;
        if (settings1->TRIG_ANALYSIS_MODE == 2) { 
            // TRIG_ANALYSIS_MODE=2 is the noise-only mode, 
            //   say that all antennas have good signal
            ants_with_nonzero_signal = 16; 
        }
        else {
            // If rays connected to the station, count how many antennas 
            //   have sufficient signal from the signal-only waveform
            if (stations[i].Total_ray_sol) { 
                ants_with_nonzero_signal = getNumOfSignalledAnts(stations[i]);
            } 
        }

        if (stations[i].Total_ray_sol && ants_with_nonzero_signal)
        {
            // if there is any ray_sol (don't check trigger if there is no ray_sol at all)

            // calculate total number of bins we need to do trigger check
            max_total_bin = (stations[i].max_arrival_time - stations[i].min_arrival_time) / settings1->TIMESTEP + settings1->NFOUR *3 + trigger->maxt_diode_bin;    // make more time

            //stations[i].max_total_bin = max_total_bin;

            // test generating new noise waveform for only stations if there's any ray trace solutions
            if (settings1->NOISE_WAVEFORM_GENERATE_MODE == 0)
            {
                // noise waveforms will be generated for each evts

                // redefine DATA_BIN_SIZE
                int DATA_BIN_SIZE_tmp;
                for (int DBS = 10; DBS < 16; DBS++)
                {
                    DATA_BIN_SIZE_tmp = (int) pow(2., (double) DBS);
                    if (DATA_BIN_SIZE_tmp > max_total_bin) DBS = 16;    // come out
                }
                settings1->DATA_BIN_SIZE = DATA_BIN_SIZE_tmp;
                // cout<<"new DATA_BIN_SIZE : "<<DATA_BIN_SIZE_tmp<<endl;
                // cout<<"max_total_bin : "<<max_total_bin<<endl;

                // reset all filter values in Detector class
                detector->get_NewDiodeModel(settings1);
                detector->ReadFilter_New(settings1);
                detector->ReadPreamp_New(settings1);
                detector->ReadFOAM_New(settings1);

                if (settings1->USE_TESTBED_RFCM_ON == 1)
                {
                    detector->ReadRFCM_New(settings1);
                }
                if (settings1->NOISE == 1 && settings1->DETECTOR == 3)
                {
                    detector->ReadRayleigh_New(settings1);
                }

                // TODO: I think this is where the Rayleigh reading will go for this next version of the code
                // if (settings1->NOISE==1 && settings1->DETECTOR==4 || settings1->DETECTOR==5) {
                //     detector->ReadRayleigh_Station(settings1);
                // }

                // reset Trigger class noise temp values
                trigger->Reset_V_noise_freqbin(settings1, detector);

                // now call new noise waveforms with new DATA_BIN_SIZE
                trigger->GetNewNoiseWaveforms(settings1, detector, this);
                // cout << "New noise waveforms gotten" << endl;
            }

            // do only if it's not in debugmode
            if (debugmode == 1) N_noise = 1;
            // now, check if DATA_BIN_SIZE is enough for total time delay between antennas
            else
            {
                N_noise = (int)(max_total_bin / settings1->DATA_BIN_SIZE) + 1;
            }
            // cout<<"N_noise : "<<N_noise<<endl;

            if (N_noise > 1) cout << "N_noise : " << N_noise << " max_total_bin : " << max_total_bin << " might cause error!!" << endl;
            // mostly N_noise should be "1"

            // now, check the number of bins we need for portion of noise waveforms
            remain_bin = max_total_bin % settings1->DATA_BIN_SIZE;
            ch_ID = 0;
            int ants_with_sufficient_SNR = 0;
            for (int j = 0; j < detector->stations[i].strings.size(); j++)
            {
                // cout << j << endl;
                for (int k = 0; k < detector->stations[i].strings[j].antennas.size(); k++)
                {

                    // Fill trigger->Full_window and trigger->Full_window_V with noise waveforms
                    Prepare_Antenna_Noise(debugmode, ch_ID, i, j, k, settings1, trigger, detector);

                    // currently there is a initial spoiled bins (maxt_diode_bin) 
                    // at the initial Full_window "AND" at the initial of connected noisewaveform 
                    // (we can fix this by adding values but not accomplished yet)

                    // Combine signals from all rays, add noise, and convolve them through the tunnel diode

                    if (debugmode == 0) {
                        Convolve_Signals(
                            &stations[i].strings[j].antennas[k], ch_ID, i,
                            settings1, trigger, detector
                        );
                    }

                    // Get the SNR from this antenna
                    double ant_SNR = 0.;
                    if (settings1->TRIG_ANALYSIS_MODE==2){
                        // TRIG_ANALYSIS_MODE=2 is the noise-only mode
                        // Set the SNR to an arbitarily high number so it passes
                        //   the coming insufficient signal check
                        ant_SNR = 100.;
                    }
                    else {
                        // For all other simulations, use previously calculated noise RMS
                        
                        // Steal the noise RMS from the trigger class and pass 
                        //   it as the noise WF to get_SNR() (since the RMS of an 
                        //   array with one element is the absolute value of that element)
                        vector <double> tmp_noise_RMS;
                        int trigger_ch_ID = GetChNumFromArbChID(detector, ch_ID, i, settings1) - 1;
                        double ant_noise_voltage_RMS = trigger->GetAntNoise_voltageRMS(trigger_ch_ID, settings1);
                        tmp_noise_RMS.push_back( ant_noise_voltage_RMS );

                        // Calculate SNR in this antenna
                        ant_SNR = get_SNR( 
                            stations[i].strings[j].antennas[k].V_convolved, 
                            tmp_noise_RMS);

                    }

                    // Log if this antenna has a strong SNR or not
                    if ( ant_SNR > 0.01 ) ants_with_sufficient_SNR++;

                    // Apply gain factors
                    if ((debugmode == 0) && (settings1->USE_MANUAL_GAINOFFSET == 1) )
                    {
                        Apply_Gain_Offset(settings1, trigger, detector, ch_ID, i);

                    } 

                    ch_ID++;    // now to next channel

                }   // for antennas

            }   // for strings

            // Check for a trigger on this event
            // do only if it's not in debugmode
            //   and if at least 1 antenna has sufficient signal
            if ( (debugmode == 0) && (ants_with_sufficient_SNR) )
            {

                // before we move to next station, do trigger check here!!!

                int trig_i, trig_j, trig_bin;
                int trig_search_init;

                // parts that are added for fixed non-trigger passed chs' V_mimic (fixed V_mimic)
                int last_trig_bin;  // stores last trigger passed bin number
                // mode select for non-trigger passed chs' V_mimic
                int V_mimic_mode = settings1->V_MIMIC_MODE;
                // 0 for orginal style (save the middle of trig_window)
                // 1 for saving waveform starting from last_trig_bin
                // 2 for saving waveform where last_trig_bin located on the middle of the waveform

                int trig_mode = settings1->TRIG_MODE;
                // global trigger mode
                // 0 for orginal N_TRIG out of all channels
                // 1 for new stations, N_TRIG_V out of Vpol channels or N_TRIG_H out of Hpol channels

                int check_ch;

                trig_search_init = trigger->maxt_diode_bin + settings1->NFOUR;  // give some time shift for mimicing force trig events
                trig_i = trig_search_init;

                // save trig search bin info default (if no global trig, this value saved)
                stations[i].total_trig_search_bin = max_total_bin - trig_search_init;

                if (settings1->TRIG_SCAN_MODE==5){ // Trigger mode for phased array
                    checkPATrigger(
                        i, all_receive_ang, viewangle, ray_sol_cnt, 
                        detector, event, evt, trigger, settings1, 
                        trig_search_init, max_total_bin
                    );
                }
                else if (settings1->TRIG_SCAN_MODE == 0)
                {
                    // ********************old mode left as-is ********************

                    // avoid really long trig_window_bin case (change trig_window to check upto max_total_bin)
                    if (max_total_bin - trig_window_bin <= trig_i) trig_window_bin = max_total_bin - trig_i - 1;

                    //while (trig_i < settings1->DATA_BIN_SIZE - trig_window_bin) {
                    while (trig_i < max_total_bin - trig_window_bin)
                    {

                        N_pass = 0;
                        N_pass_V = 0;
                        N_pass_H = 0;
                        last_trig_bin = 0;
                        Passed_chs.clear();

                        //for (trig_j=0; trig_j < ch_ID; trig_j++) {        // loop over all channels
                        trig_j = 0;
                        while (trig_j < ch_ID)
                        {

                            int string_i = detector->getStringfromArbAntID(i, trig_j);
                            int antenna_i = detector->getAntennafromArbAntID(i, trig_j);
			    int channel_num = detector->GetChannelfromStringAntenna(i, string_i, antenna_i, settings1);

			    if (!(settings1->DETECTOR==4 || settings1->DETECTOR==5)){
 			    	channel_num = channel_num+1; // Channel numbering is different for DETECTOR=(1,2,3) than for DETECTOR = 4 in GetChannelfromStringAntenna(), it needs that shift 
 			    }


			    if( detector->GetTrigMasking(channel_num-1)==0){ //Antenna Masking (masked_ant=0 means this antenna should be ignored from trigger)
				trig_j++;
				continue;	
			    }

			    int offset = detector->GetTrigOffset(channel_num-1, settings1);

                            // If Testbed simulation, check if we want to use BH chs only for trigger analysis
                            if ((settings1->TRIG_ONLY_BH_ON == 1) && (settings1->DETECTOR == 3))
                            {

                                // check if this channel is BH ch (DAQchan)
                                if (detector->stations[i].strings[string_i].antennas[antenna_i].DAQchan == 0)
                                {

                                    trig_bin = 0;
                                    while (trig_bin < trig_window_bin)
                                    {

                                        //cout<<"trig_bin : "<<trig_bin<<endl;

                                        double diode_noise_RMS = trigger->GetAntNoise_diodeRMS(channel_num-1, settings1);
                                        if (trigger->Full_window[trig_j][trig_i + trig_bin] < (detector->GetThres(i, channel_num - 1, settings1) * diode_noise_RMS *detector->GetThresOffset(i, channel_num - 1, settings1)))
                                        {
                                            // if this channel passed the trigger!
                                            //cout<<"trigger passed at bin "<<trig_i+trig_bin<<" ch : "<<trig_j<<endl;
                                            //stations[i].strings[(int)((trig_j)/4)].antennas[(int)((trig_j)%4)].Trig_Pass = trig_i+trig_bin;
                                            stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = trig_i + trig_bin;
                                            N_pass++;
                                            if (detector->stations[i].strings[string_i].antennas[antenna_i].type == 0)
                                            {
                                                // Vpol
                                                N_pass_V++;
                                            }
                                            if (detector->stations[i].strings[string_i].antennas[antenna_i].type == 1)
                                            {
                                                // Hpol
                                                N_pass_H++;
                                            }
                                            if (last_trig_bin < trig_i + trig_bin) last_trig_bin = trig_i + trig_bin;   // added for fixed V_mimic
                                            trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
                                            Passed_chs.push_back(trig_j);
                                        }

                                        trig_bin++;
                                    }
                                }
                            }

                            // For non-Testbed simulations, check if we just want to use first/lower 8 chs' thres values
                            else if ((settings1->TRIG_ONLY_LOW_CH_ON == 1) && (settings1->DETECTOR != 3))
                            {

                                // reset channel numbers so that bottom antennas have ch 1-8
                                channel_num = GetChannelNum8_LowAnt(string_i, antenna_i);

                                if (antenna_i < 2)
                                {
                                    // only antenna 0, 1 which are bottom 2 antennas

                                    // set channel_num as new value (antenna 0, 1 are only possible antennas for channel_num 1 - 8)

                                    trig_bin = 0;
                                    while (trig_bin < trig_window_bin)
                                    {

                                        double diode_noise_RMS = trigger->GetAntNoise_diodeRMS(channel_num-1, settings1);
                                        // with threshold offset by chs
                                        if (trigger->Full_window[trig_j][trig_i + trig_bin] < (detector->GetThres(i, channel_num - 1, settings1) * diode_noise_RMS *detector->GetThresOffset(i, channel_num - 1, settings1)))
                                        {
                                            // if this channel passed the trigger!
                                            //cout<<"trigger passed at bin "<<trig_i+trig_bin<<" ch : "<<trig_j<<endl;
                                            //stations[i].strings[(int)((trig_j)/4)].antennas[(int)((trig_j)%4)].Trig_Pass = trig_i+trig_bin;
                                            stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = trig_i + trig_bin;
                                            N_pass++;
                                            if (detector->stations[i].strings[string_i].antennas[antenna_i].type == 0)
                                            {
                                                // Vpol
                                                N_pass_V++;
                                            }
                                            if (detector->stations[i].strings[string_i].antennas[antenna_i].type == 1)
                                            {
                                                // Hpol
                                                N_pass_H++;
                                            }
                                            if (last_trig_bin < trig_i + trig_bin) last_trig_bin = trig_i + trig_bin;   // added for fixed V_mimic
                                            trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
                                            Passed_chs.push_back(trig_j);
                                        }

                                        trig_bin++;
                                    }
                                }
                            }

                            //else if (settings1->TRIG_ONLY_BH_ON == 0) {
                            // other cases: use all possible chs for trigger analysis
                            else
                            {
                                
                                trig_bin = 0;
                                while (trig_bin < trig_window_bin)
                                {

                                    double diode_noise_RMS = trigger->GetAntNoise_diodeRMS(channel_num-1, settings1);
                                    if( trig_i+offset+trig_bin >= settings1->DATA_BIN_SIZE ) break; //if trigger window hits wf end, cannot scan this channel further with this trig_i
                                    if (trigger->Full_window[trig_j][trig_i + trig_bin + offset] < (detector->GetThres(i, channel_num - 1, settings1) * diode_noise_RMS *detector->GetThresOffset(i, channel_num - 1, settings1)))
                                    {
                                        // if this channel passed the trigger!
                                        //cout<<"trigger passed at bin "<<trig_i+trig_bin<<" ch : "<<trig_j<<endl;
                                        //stations[i].strings[(int)((trig_j)/4)].antennas[(int)((trig_j)%4)].Trig_Pass = trig_i+trig_bin;
                                        stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = trig_i + trig_bin + offset;
                                        N_pass++;
                                        if (detector->stations[i].strings[string_i].antennas[antenna_i].type == 0)
                                        {
                                            // Vpol
                                            N_pass_V++;
                                        }
                                        if (detector->stations[i].strings[string_i].antennas[antenna_i].type == 1)
                                        {
                                            // Hpol
                                            N_pass_H++;
                                        }
                                        if (last_trig_bin < trig_i + trig_bin) last_trig_bin = trig_i + trig_bin;   // added for fixed V_mimic
                                        trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
                                        Passed_chs.push_back(trig_j);
                                    }

                                    trig_bin++;
                                }
                            }

                            // check all triggered channels not just 3
                            //if (N_pass > 2) trig_j += ch_ID;      // if the number of passed channels is 3 or more, no need to check other remaining channels as this station is trigged!
                            //else trig_j++;    // if station not passed the trigger, just go to next channel
                            trig_j++;   // if station not passed the trigger, just go to next channel

                        }   // while trig_j < ch_ID

                        //if (N_pass > settings1->N_TRIG-1) {   // now as global trigged!! = more or eq to N_TRIG triggered
                        if (((trig_mode == 0) && (N_pass > settings1->N_TRIG - 1))  // trig_mode = 0 case!
                            ||  // or
                            ((trig_mode == 1) && ((N_pass_V > settings1->N_TRIG_V - 1) || (N_pass_H > settings1->N_TRIG_H - 1)))    // trig_mode = 1 case!
                       )
                        {

                            check_ch = 0;
                            //stations[i].Global_Pass = trig_i;
                            stations[i].Global_Pass = last_trig_bin;    // where actually global trigger occured

                            //trig_i = settings1->DATA_BIN_SIZE;    // also if we know this station is trigged, don't need to check rest of time window
                            trig_i = max_total_bin; // also if we know this station is trigged, don't need to check rest of time window
                            for (int ch_loop = 0; ch_loop < ch_ID; ch_loop++)
                            {
                                //cout << ch_loop << "/" << ch_ID << endl;
                                int string_i = detector->getStringfromArbAntID(i, ch_loop);
                                int antenna_i = detector->getAntennafromArbAntID(i, ch_loop);
                                //          cout << "string:antenna: " << string_i << " : " << antenna_i << endl;
                                stations[i].strings[string_i].antennas[antenna_i].Likely_Sol = -1;  // no likely init
                                if (ch_loop == Passed_chs[check_ch] && check_ch < N_pass)
                                {
                                    // added one more condition (check_ch < N_Pass) for bug in vector Passed_chs.clear()???

                                    // store which ray sol is triggered based on trig time
                                    //

                                    int mindBin = 1.e9; // init big value
                                    int dBin = 0;

                                    for (int m = 0; m < stations[i].strings[string_i].antennas[antenna_i].ray_sol_cnt; m++)
                                    {
                                        // loop over raysol numbers

                                        if (stations[i].strings[string_i].antennas[antenna_i].SignalExt[m])
                                        {

                                            dBin = abs(stations[i].strings[string_i].antennas[antenna_i].SignalBin[m] - stations[i].strings[string_i].antennas[antenna_i].Trig_Pass);

                                            if (dBin < mindBin)
                                            {
                                                stations[i].strings[string_i].antennas[antenna_i].Likely_Sol = m;   // store the ray sol number which is minimum difference between Trig_Pass bin
                                                mindBin = dBin;
                                            }
                                        }
                                    }

                                    //skip this passed ch as it already has bin info
                                    //cout<<"trigger passed at bin "<<stations[i].strings[(int)((ch_loop)/4)].antennas[(int)((ch_loop)%4)].Trig_Pass<<"  passed ch : "<<ch_loop<<" Direct dist btw posnu : "<<event->Nu_Interaction[0].posnu.Distance(detector->stations[i].strings[(int)((ch_loop)/4)].antennas[(int)((ch_loop)%4)])<<endl;
                                    //cout<<endl<<"trigger passed at bin "<<stations[i].strings[string_i].antennas[antenna_i].Trig_Pass<<"  passed ch : "<<ch_loop<<" ("<<detector->stations[i].strings[string_i].antennas[antenna_i].type<<"type) Direct dist btw posnu : "<<event->Nu_Interaction[0].posnu.Distance(detector->stations[i].strings[string_i].antennas[antenna_i])<<" noiseID : "<<stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];
                                    //cout<<endl<<"trigger passed at bin "<<stations[i].strings[string_i].antennas[antenna_i].Trig_Pass<<"  passed ch : "<<ch_loop<<" ("<<detector->stations[i].strings[string_i].antennas[antenna_i].type<<"type) Direct dist btw posnu : "<<event->Nu_Interaction[0].posnu.Distance(detector->stations[i].strings[string_i].antennas[antenna_i])<<" noiseID : "<<stations[i].strings[string_i].antennas[antenna_i].noise_ID[0]<<" ViewAngle : "<<stations[i].strings[string_i].antennas[antenna_i].view_ang[0]*DEGRAD;

                                    if (settings1->TRIG_ONLY_LOW_CH_ON == 0)
                                    {
                                        cout << endl << "trigger passed at bin " << stations[i].strings[string_i].antennas[antenna_i].Trig_Pass << "  passed ch : " << ch_loop << " (" << detector->stations[i].strings[string_i].antennas[antenna_i].type << "type) Direct dist btw posnu : " << event->Nu_Interaction[0].posnu.Distance(detector->stations[i].strings[string_i].antennas[antenna_i]) << " noiseID : " << stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];
                                        if (stations[i].strings[string_i].antennas[antenna_i].Likely_Sol != -1)
                                        {
                                            cout << " ViewAngle : " << stations[i].strings[string_i].antennas[antenna_i].view_ang[0] *DEGRAD << " LikelyTrigSignal : " << stations[i].strings[string_i].antennas[antenna_i].Likely_Sol;
                                        }
                                    }
                                    else if (settings1->TRIG_ONLY_LOW_CH_ON == 1)
                                    {
                                        cout << endl << "trigger passed at bin " << stations[i].strings[string_i].antennas[antenna_i].Trig_Pass << "  passed ant: str[" << string_i << "].ant[" << antenna_i << "] (" << detector->stations[i].strings[string_i].antennas[antenna_i].type << "type) Direct dist btw posnu : " << event->Nu_Interaction[0].posnu.Distance(detector->stations[i].strings[string_i].antennas[antenna_i]) << " noiseID : " << stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];
                                        if (stations[i].strings[string_i].antennas[antenna_i].Likely_Sol != -1)
                                        {
                                            cout << " ViewAngle : " << stations[i].strings[string_i].antennas[antenna_i].view_ang[0] *DEGRAD << " LikelyTrigSignal : " << stations[i].strings[string_i].antennas[antenna_i].Likely_Sol;
                                        }
                                    }

                                    //cout << "True station number: " << detector->InstalledStations[0].VHChannel[string_i][antenna_i] << endl;
                                    //cout << event->Nu_Interaction[0].posnu[0] << " : " <<  event->Nu_Interaction[0].posnu[1] << " : " << event->Nu_Interaction[0].posnu[2] << endl;
                                    check_ch++;

                                    // now save the voltage waveform to V_mimic
                                    //

                                    for (int mimicbin = 0; mimicbin < waveformLength; mimicbin++)
                                    {

                                        // new DAQ waveform writing mechanism test
                                        if (V_mimic_mode == 0)
                                        {
                                            // Global passed bin is the center of the window
                                            stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][last_trig_bin + waveformCenter - waveformLength / 2 + mimicbin]) *1.e3);   // save in mV
                                            stations[i].strings[string_i].antennas[antenna_i].time.push_back(last_trig_bin + waveformCenter - waveformLength / 2 + mimicbin);
                                            stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back((waveformCenter - waveformLength / 2 + mimicbin) *settings1->TIMESTEP *1.e9);    // save in ns
                                        }
                                        else if (V_mimic_mode == 1)
                                        {
                                            // Global passed bin is the center of the window + delay to each chs from araGeom
                                            stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin]) *1.e3);   // save in mV
                                            stations[i].strings[string_i].antennas[antenna_i].time.push_back(last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin);
                                            stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back((-(detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin) *settings1->TIMESTEP *1.e9);   // save in ns
                                        }
                                        else if (V_mimic_mode == 2)
                                        {
                                            // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
                                            stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin]) *1.e3);    // save in mV
                                            stations[i].strings[string_i].antennas[antenna_i].time.push_back(last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin);
                                            stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back((-(detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin) *settings1->TIMESTEP *1.e9 + detector->params.TestBed_WFtime_offset_ns);    // save in ns
                                        }
                                        if (mimicbin == 0) {
                                            for (int m = 0; m < stations[i].strings[string_i].antennas[antenna_i].ray_sol_cnt; m++) { ///< calculates time of center of each rays signal based on readout window time config
                                                double signal_center_offset = (double)(stations[i].strings[string_i].antennas[antenna_i].SignalBin[m] - stations[i].strings[string_i].antennas[antenna_i].time[0]) * settings1->TIMESTEP * 1.e9;
                                                double signal_center_time = signal_center_offset + stations[i].strings[string_i].antennas[antenna_i].time_mimic[0];
                                                //! signal_center_offset: time offset between beginning of readout window and center of signal
                                                //! signal_center_time: time of center of signal based on readout window time config
                                                stations[i].strings[string_i].antennas[antenna_i].SignalBinTime.push_back(signal_center_time);
                                            }
                                        }
                                    }

                                    // set global_trig_bin values
                                    if (V_mimic_mode == 0)
                                    {
                                        // Global passed bin is the center of the window
                                        stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (waveformLength / 2 - waveformCenter);
                                    }
                                    else if (V_mimic_mode == 1)
                                    {
                                        // Global passed bin is the center of the window + delay to each chs from araGeom
                                        stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin) - waveformCenter + waveformLength / 2;
                                    }
                                    else if (V_mimic_mode == 2)
                                    {
                                        // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
                                        stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) - waveformCenter + waveformLength / 2;
                                    }
                                }
                                else
                                {
                                    //stations[i].strings[(int)((ch_loop)/4)].antennas[(int)((ch_loop)%4)].Trig_Pass = stations[i].Global_Pass + settings1->NFOUR/4;    // so that global trig is 
                                    //stations[i].strings[(int)((ch_loop)/4)].antennas[(int)((ch_loop)%4)].Trig_Pass = 0.;
                                    stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = 0.;

                                    // new DAQ waveform writing mechanism test
                                    for (int mimicbin = 0; mimicbin < waveformLength; mimicbin++)
                                    {
                                        if (V_mimic_mode == 0)
                                        {
                                            // Global passed bin is the center of the window
                                            stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][last_trig_bin + waveformCenter - waveformLength / 2 + mimicbin]) *1.e3);   // save in mV
                                            stations[i].strings[string_i].antennas[antenna_i].time.push_back(last_trig_bin + waveformCenter - waveformLength / 2 + mimicbin);
                                            stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back((waveformCenter - waveformLength / 2 + mimicbin) *settings1->TIMESTEP *1.e9);    // save in ns
                                        }
                                        else if (V_mimic_mode == 1)
                                        {
                                            // Global passed bin is the center of the window + delay to each chs from araGeom
                                            stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin]) *1.e3);   // save in mV
                                            stations[i].strings[string_i].antennas[antenna_i].time.push_back(last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin);
                                            stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back((-(detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin) *settings1->TIMESTEP *1.e9);   // save in ns
                                        }
                                        if (V_mimic_mode == 2)
                                        {
                                            // Global passed bin is the center of the window + delay to each chs from araGeom + fitted delay
                                            //stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back(trigger->Full_window_V[ch_loop][ last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop]-detector->params.TestBed_BH_Mean_delay_bin) - settings1->NFOUR/4 + mimicbin ]);
                                            stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin]) *1.e3);    // save in mV
                                            stations[i].strings[string_i].antennas[antenna_i].time.push_back(last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin);
                                            stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back((-(detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin) *settings1->TIMESTEP *1.e9 + detector->params.TestBed_WFtime_offset_ns);    // save in ns
                                        }
                                        if (mimicbin == 0) {
                                            for (int m = 0; m < stations[i].strings[string_i].antennas[antenna_i].ray_sol_cnt; m++) { ///< calculates time of center of each rays signal based on readout window time config
                                                double signal_center_offset = (double)(stations[i].strings[string_i].antennas[antenna_i].SignalBin[m] - stations[i].strings[string_i].antennas[antenna_i].time[0]) * settings1->TIMESTEP * 1.e9;
                                                double signal_center_time = signal_center_offset + stations[i].strings[string_i].antennas[antenna_i].time_mimic[0];
                                                //! signal_center_offset: time offset between beginning of readout window and center of signal
                                                //! signal_center_time: time of center of signal based on readout window time config
                                                stations[i].strings[string_i].antennas[antenna_i].SignalBinTime.push_back(signal_center_time);
                                            }
                                        }
                                    }

                                    // set global_trig_bin values
                                    if (V_mimic_mode == 0)
                                    {
                                        // Global passed bin is the center of the window
                                        stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = waveformLength / 2 - waveformCenter;
                                    }
                                    else if (V_mimic_mode == 1)
                                    {
                                        // Global passed bin is the center of the window + delay to each chs from araGeom
                                        stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin) + waveformLength / 2 - waveformCenter;
                                    }
                                    else if (V_mimic_mode == 2)
                                    {
                                        // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
                                        stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) + waveformLength / 2 - waveformCenter;
                                    }

                                    // done V_mimic for non-triggered chs (done fixed V_mimic)
                                }

                                double arrivtime = stations[i].strings[string_i].antennas[antenna_i].arrival_time[0];
                                double X = detector->stations[i].strings[string_i].antennas[antenna_i].GetX();
                                double Y = detector->stations[i].strings[string_i].antennas[antenna_i].GetY();
                                double Z = detector->stations[i].strings[string_i].antennas[antenna_i].GetZ();
                                //std::cout << "Arrival time:X:Y:Z " << arrivtime << " : " << X << " : " << Y << " : " << Z << std::endl;
                                /*
                                int AraRootChannel = 0;
                                AraRootChannel = detector->GetChannelfromStringAntenna (i, string_i, antenna_i, settings1);

                                int UsefulEventBin;
                                if (settings1->NFOUR/2<EFFECTIVE_LAB3_SAMPLES*2) UsefulEventBin = settings1->NFOUR/2;
                                else UsefulEventBin = EFFECTIVE_LAB3_SAMPLES*2;

                                //for (int mimicbin=0; mimicbin<settings1->NFOUR/2; mimicbin++) {
                                for (int mimicbin=0; mimicbin < UsefulEventBin; mimicbin++) {
                                    theUsefulEvent->fVoltsRF[AraRootChannel-1][mimicbin] = stations[i].strings[string_i].antennas[antenna_i].V_mimic[mimicbin];
                                    //theUsefulEvent->fTimesRF[AraRootChannel-1][mimicbin] = double(stations[i].strings[string_i].antennas[antenna_i].time[mimicbin])*settings1->TIMESTEP*1.0E9;
                                    theUsefulEvent->fTimesRF[AraRootChannel-1][mimicbin] = stations[i].strings[string_i].antennas[antenna_i].time_mimic[mimicbin];
                                    //cout << theUsefulEvent->fVoltsRF[ch_loop][mimicbin] << endl;
                                    //cout << theUsefulEvent->fTimesRF[ch_loop][mimicbin] <<endl;
                                }
                                //theUsefulEvent->fNumPointsRF[ch_loop] = EFFECTIVE_SAMPLES * 2;
                                theUsefulEvent->fNumPointsRF[ch_loop] = UsefulEventBin;
                                //cout << " : " << theUsefulEvent->fNumPointsRF[ch_loop] << endl;
                                 */
                            }

                            //cout<<"Global trigger passed!!, N_pass : "<<N_pass<<endl;
                            Passed_chs.clear();

                            // save trig search bin info
                            stations[i].total_trig_search_bin = stations[i].Global_Pass + trig_window_bin - trig_search_init;
                        }   // if global trig!
                        else
                        {
                            trig_i++;   // also if station not passed the trigger, just go to next bin
                            /*
                            for (int ch_loop=0; ch_loop < ch_ID; ch_loop++) {
                                int string_i = detector->getStringfromArbAntID(i, ch_loop);
                                int antenna_i = detector->getAntennafromArbAntID(i, ch_loop);
                                int AraRootChannel = 0;
                                AraRootChannel = detector->GetChannelfromStringAntenna (i, string_i, antenna_i, settings1);

                                int UsefulEventBin;
                                if (settings1->NFOUR/2<EFFECTIVE_LAB3_SAMPLES*2) UsefulEventBin = settings1->NFOUR/2;
                                else UsefulEventBin = EFFECTIVE_LAB3_SAMPLES*2;

                                //for (int mimicbin=0; mimicbin<settings1->NFOUR/2; mimicbin++) {
                                for (int mimicbin=0; mimicbin < UsefulEventBin; mimicbin++) {
                                    theUsefulEvent->fVoltsRF[AraRootChannel-1][mimicbin] = 0;
                                    theUsefulEvent->fTimesRF[AraRootChannel-1][mimicbin] = 0;
                                    //cout << theUsefulEvent->fVoltsRF[ch_loop][mimicbin] << endl;
                                    //cout << theUsefulEvent->fTimesRF[ch_loop][mimicbin] <<endl;
                                }
                                //theUsefulEvent->fNumPointsRF[ch_loop] = EFFECTIVE_SAMPLES * 2;
                                theUsefulEvent->fNumPointsRF[ch_loop] = UsefulEventBin;
                            }
                            */
                        }
                    }   // while trig_i

                }   // if TRIG_SCAN_MODE==0

                else if (settings1->TRIG_SCAN_MODE > 0) triggerCheckLoop(settings1, detector, event, trigger, i, trig_search_init, max_total_bin, trig_window_bin, settings1->TRIG_SCAN_MODE);

                //         else if(settings1->TRIG_SCAN_MODE==2) triggerCheckLoopScan();

                //         else if(settings1->TRIG_SCAN_MODE==3) triggerCheckLoopScanNumbers();
            }   // if it's not debugmode

            // delete noise waveforms 
            if (settings1->NOISE_WAVEFORM_GENERATE_MODE == 0)
            {
                // noise waveforms will be generated for each evts
                // remove noise waveforms for next evt
                trigger->ClearNoiseWaveforms();
            }
        }   // if there is any ray_sol in the station

    }   // for stations

    // also clear all vector info to reduce output root file size
    clear_useless(settings1);   // to reduce the size of output AraOut.root, remove some information

}   // end Connect_Interaction_Detector


void Report::rerun_event(Event *event, Detector *detector, 
    RaySolver *raysolver, Signal *signal, 
    IceModel *icemodel, Settings *settings, int which_solution,
    vector<int> &numSolutions, vector<vector<vector<double> > > &traceTimes,
    vector<vector<vector<double> > > &traceVoltages
    ){

    // which solution = 0 (direct only), 1 (refracted/reflected only), 2 (both)
    std::vector<int> which_sols_to_search;
    if (which_solution==0){
        which_sols_to_search.push_back(0);
    }
    else if(which_solution==1){
        which_sols_to_search.push_back(1);
    }
    else if(which_solution==2){
        which_sols_to_search.push_back(0);
        which_sols_to_search.push_back(1);
    }

    // we create T_forint, which is the final time sampling of the traces before they are combined
    double RandomTshift = 0. ; // for replication, we have no choice by to set this
    double init_T;
    double T_forint[settings->NFOUR/2];

    int num_strings = detector->stations[0].strings.size();
    int num_antennas = detector->stations[0].strings[0].antennas.size();

    for(int j=0; j<num_strings; j++){

        // Normal way of calculating init_T and T_forint - Using user defined TIMESTEP in settings
        init_T = settings->TIMESTEP *-1.e9 *((double) settings->NFOUR / 4 + RandomTshift);    // locate zero time at the middle and give random time shift
        for (int n = 0; n < settings->NFOUR / 2; n++)
        {
            T_forint[n] = init_T + (double) n *settings->TIMESTEP *1.e9;   // in ns
        }

        for(int k=0; k<num_antennas; k++){
	   
            // This (gain_ch_no) is used for per-channel gain implementation. 
            // It is used in all instances of ApplyElect_Tdomain() and ApplyElect_Tdomain_FirstTwo(), to indicate channel number
            // Note that channel numbering is different for DETECTOR==4 than for the rest (1-3). See that in the definition of GetChannelfromStringAntenna()  
            int gain_ch_no;
            if (settings->DETECTOR==4 || settings->DETECTOR==5){
                    gain_ch_no = detector->GetChannelfromStringAntenna (0, j, k, settings)-1;
            }       
            else{
                    gain_ch_no = detector->GetChannelfromStringAntenna (0, j, k, settings);
            }
  

            int idx = ((j*4)+k);

            // redo the ray tracing
            vector<vector<double> > ray_output;
            vector<vector<vector<double> > > Ray_Step;
            raysolver->Solve_Ray(event->Nu_Interaction[0].posnu,
                detector->stations[0].strings[j].antennas[k],
                icemodel, ray_output, settings, Ray_Step);

            // if there is a solution
            if (raysolver->solution_toggle){

                // numSolutions[idx]=ray_output[0].size();
                numSolutions[idx]=ray_output[0].size();

                int ray_sol_cnt=0;
                while (ray_sol_cnt < ray_output[0].size()){

                    if(std::find(which_sols_to_search.begin(), which_sols_to_search.end(), ray_sol_cnt) == which_sols_to_search.end()){
                        // if this solution is not meant to be searched, skip it
                        ray_sol_cnt++;
                        continue;
                    }

                    double viewangle = ray_output[1][ray_sol_cnt];
                    Position launch_vector;
                    Position receive_vector;
                    Vector n_trg_pokey;
                    Vector n_trg_slappy;

                    GetParameters(
                        event->Nu_Interaction[0].posnu,
                        detector->stations[0].strings[j].antennas[k],
                        event->Nu_Interaction[0].nnu,
                        viewangle,
                        ray_output[2][ray_sol_cnt],
                        launch_vector, receive_vector,
                        n_trg_slappy, n_trg_pokey
                        );
                    // get polarization
                    Vector Pol_vector = GetPolarization(
                        event->Nu_Interaction[0].nnu, launch_vector);

                    double fresnel, mag;
                    icemodel->GetFresnel(
                        ray_output[1][ray_sol_cnt],
                        ray_output[2][ray_sol_cnt],
                        ray_output[3][ray_sol_cnt],
                        event->Nu_Interaction[0].posnu,
                        launch_vector, receive_vector,
                        settings, fresnel, Pol_vector
                        );
                    icemodel->GetMag(
                        mag, 
                        ray_output[0][ray_sol_cnt], // ray path length
                        ray_output[1][ray_sol_cnt], // zenith angle of ray at launch
                        ray_output[2][ray_sol_cnt], // zenith angle of ray upon receipt
                        ray_sol_cnt,
                        event->Nu_Interaction[0].posnu, // Neutrino
                        detector->stations[0].strings[j].antennas[k], // Antenna
                        -0.01, // 1cm antenna shift, inspired from NuRadioMC
                        icemodel, settings
                    );

                    // get arrival angle at the antenna
                    double antenna_theta, antenna_phi;
                    GetAngleAnt(
                        receive_vector, 
                        detector->stations[0].strings[j].antennas[k],
                        antenna_theta, antenna_phi
                        );
                        
                    // get launch angle of signal
                    double launch_theta, launch_phi;
                    GetAngleLaunch(
                        launch_vector, 
                        launch_theta, launch_phi
                        );                        
                    
                    // this is the 1/R and fresnel and focusing effect
                    double atten_factor = 1. / ray_output[0][ray_sol_cnt] * mag * fresnel;

                    // get the electric field
                    int local_outbin=64;
                    double local_Earray[local_outbin];
                    double local_Tarray[local_outbin];
                    int local_skipbins;
                    signal->GetVm_FarField_Tarray(event, settings, viewangle,
                        atten_factor, local_outbin, local_Tarray, local_Earray, local_skipbins
                        );

                    double dT_forfft = local_Tarray[1] - local_Tarray[0];
                    int Ntmp = settings->TIMESTEP*1.e9 / dT_forfft;
                    int Nnew=1;
                    while(Ntmp > 1){
                        Ntmp = Ntmp/2;
                        Nnew = Nnew*2;
                    }
                    Nnew = Nnew * settings->NFOUR/2;

                    double V_forfft[Nnew];
                    double T_forfft[Nnew];

                    for(int n=0; n<Nnew; n++){
                        T_forfft[n] = local_Tarray[local_outbin/2] - (dT_forfft*(double)(Nnew/2 - n));
                        if ( (n >= Nnew/2 - local_outbin/2) && (n < Nnew/2 + local_outbin/2) ) {
                            V_forfft[n] = local_Earray[ n - (Nnew/2 - local_outbin/2) ];
                        }
                        else{
                            V_forfft[n] = 0.;
                        }
                    }

                    // transform to the frequency domain
                    Tools::realft(V_forfft, 1, Nnew);

                    double dF_Nnew = 1./( (double)(Nnew) * (dT_forfft)*1.e-9 ); // in Hz
                    double freq_tmp = dF_Nnew*((double)Nnew/2.+0.5);// in Hz 0.5 to place the middle of the bin and avoid zero freq
                    double freq_lastbin = freq_tmp;

                    // save the gain and heff of the last bin, which must be treated specially
                    double gain_lastbin = detector->GetGain_1D_OutZero(
                        freq_tmp*1.E-6, // Hz
                        antenna_theta, antenna_phi,
                        detector->stations[0].strings[j].antennas[k].type,
                        j, k, false
                        );
                    double heff_lastbin = GaintoHeight(gain_lastbin, freq_tmp,
                        icemodel->GetN(detector->stations[0].strings[j].antennas[k])
                        );

                    // integrate along the path and save the frequency dependent attenuation
                    double attenuations[Nnew/2];
                    std::fill_n(attenuations, Nnew/2, 1.); // all of the attenuations initially start at 1

                    for (int steps = 1; steps < (int) Ray_Step[ray_sol_cnt][0].size(); steps++) {
                        double dx = Ray_Step[ray_sol_cnt][0][steps - 1] - Ray_Step[ray_sol_cnt][0][steps];
                        double dz = Ray_Step[ray_sol_cnt][1][steps - 1] - Ray_Step[ray_sol_cnt][1][steps];
                        double dl = sqrt((dx * dx) + (dz * dz));

                        // Skipping attenuation calculation when the distance between two RaySteps is 0. Prevening adds -nan into the IceAttenFactor. (MK 2021)
                        if (dl > 0){
                            for(int n=0; n<Nnew/2; n++){
                                freq_tmp = dF_Nnew*((double)n+0.5);
                                // use ray midpoint for attenuation calculation
                                double IceAttenFactor = (  exp(-dl / icemodel->GetFreqDepIceAttenuLength(-Ray_Step[ray_sol_cnt][1][steps], freq_tmp * 1.E-9))
                                                         + exp(-dl / icemodel->GetFreqDepIceAttenuLength(-Ray_Step[ray_sol_cnt][1][steps-1], freq_tmp * 1.E-9))
                                                        )/2.; // 1e9 for conversion to GHz
                                // increase the attenuation
                                attenuations[n]*=IceAttenFactor;
                            }
                        }
                    }

                    // loop over frequency bins
                    for(int n=0; n<Nnew/2; n++){
                        freq_tmp = dF_Nnew*((double)n+0.5);

                        // apply the attenuation for this frequency bin
                        V_forfft[2*n]*=attenuations[n];
                        V_forfft[2*n + 1]*=attenuations[n];

                        // get the antenna gain (and then heff) and phase
                        double gain = detector->GetGain_1D_OutZero(
                            freq_tmp*1.E-6, // Hz
                            antenna_theta, antenna_phi,
                            detector->stations[0].strings[j].antennas[k].type,
                            j, k, false
                            );
                        double heff = GaintoHeight(gain, freq_tmp,
                            icemodel->GetN(detector->stations[0].strings[j].antennas[k])
                            );
                        double phase = detector->GetAntPhase_1D(
                            freq_tmp*1.E-6, // Hz
                            antenna_theta, antenna_phi, 
                            detector->stations[0].strings[j].antennas[k].type
                            );

                        // apply the antenna factors (this will also apply the polarization)
                        // and the electronics response
                        double Pol_factor;
                        if(n>0){
                            ApplyAntFactors_Tdomain(
                                phase, heff, Pol_vector, 
                                detector->stations[0].strings[j].antennas[k].type,
                                Pol_factor, V_forfft[2*n], V_forfft[2*n + 1],
                                settings, antenna_theta, antenna_phi, freq_tmp
                                );
                            ApplyElect_Tdomain(freq_tmp*1.e-6, detector,
                                V_forfft[2*n], V_forfft[2*n + 1], gain_ch_no, settings
                                );
                        }
                        else{
                            ApplyAntFactors_Tdomain_FirstTwo(
                                heff, heff_lastbin, Pol_vector, 
                                detector->stations[0].strings[j].antennas[k].type,
                                Pol_factor, V_forfft[2*n], V_forfft[2*n + 1],
                                antenna_theta, antenna_phi, freq_tmp
                                );
                            ApplyElect_Tdomain_FirstTwo(freq_tmp*1.e-6,
                                freq_lastbin*1.e-6, detector,
                                V_forfft[2*n], V_forfft[2*n + 1], gain_ch_no, settings
                                );
                        }
                    }

                    // fft back to time domain
                    Tools::realft(V_forfft, -1, Nnew);

                    // interpolate to final time sampling
                    double volts_forint[settings->NFOUR/2];
                    Tools::SincInterpolation(Nnew, T_forfft, V_forfft,
                        settings->NFOUR/2, T_forint, volts_forint
                    );

                    // fix the normalization on the inverse fft (2/N)
                    for(int n=0; n<settings->NFOUR/2; n++){
                        volts_forint[n] *= 2./Nnew;
                    }

                    // store the result
                    for(int n=0; n<settings->NFOUR/2; n++){
                        traceTimes[idx][ray_sol_cnt][n] = T_forint[n];
                        traceVoltages[idx][ray_sol_cnt][n] = volts_forint[n];
                    }

                    ray_sol_cnt++;
                }
            }
        }
    }
}

int Report::triggerCheckLoop(Settings *settings1, Detector *detector, Event *event, Trigger *trigger, int stationID, int trig_search_init, int max_total_bin, int trig_window_bin, int scan_mode ){
 
  int i=stationID;
  
  int numChan=stations[i].TDR_all.size();
  int numChanVpol=stations[i].TDR_Vpol_sorted.size();
  int numChanHpol=stations[i].TDR_Hpol_sorted.size();
  

  double powerthreshold=settings1->POWERTHRESHOLD;
  // this value is compared to POWERTHRESHOLD for local trigger.  
  
  int first_trigger=0;
  
  double Pthresh_value[numChan];
  CircularBuffer **buffer=new CircularBuffer*[numChan];
  for(int trig_j=0;trig_j<numChan; trig_j++){// initialize Trig_Pass and buffers
        
    int string_i = detector->getStringfromArbAntID( i, trig_j);
    int antenna_i = detector->getAntennafromArbAntID( i, trig_j);

    trig_window_bin = (int)(settings1->TRIG_WINDOW / settings1->TIMESTEP);  // coincidence window bin for trigger
    
    Pthresh_value[trig_j]=0;
    buffer[trig_j]=new CircularBuffer( trig_window_bin, powerthreshold, scan_mode);
    
    stations[i].strings[string_i].antennas[antenna_i].SingleChannelTriggers=0;
    stations[i].strings[string_i].antennas[antenna_i].TotalBinsScannedPerChannel=0;
    
//     int string_i = detector->getStringfromArbAntID( i, trig_j);
//     int antenna_i = detector->getAntennafromArbAntID( i, trig_j);
//     stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = 0.;
    
  }// for trig_j
  
  double *TDR_all_sorted_temp;
  double *TDR_Vpol_sorted_temp;
  double *TDR_Hpol_sorted_temp;
  
  if(scan_mode>1){
    
    TDR_all_sorted_temp=new double[numChan];
    TDR_Vpol_sorted_temp=new double[numChanVpol];
    TDR_Hpol_sorted_temp=new double[numChanHpol];
    
    // only need to initialize for scan_mode>1
    for(int trig_j=0;trig_j<numChan;trig_j++) TDR_all_sorted_temp[trig_j]=0;
    for(int trig_j=0;trig_j<numChanVpol;trig_j++) TDR_Vpol_sorted_temp[trig_j]=0;
    for(int trig_j=0;trig_j<numChanHpol;trig_j++) TDR_Hpol_sorted_temp[trig_j]=0;
          
  }// if scan_mode>1

  int global_pass_bit=0; // whether this event passes (in all windows)
  int check_TDR_configuration=0; // check if we need to reorder our TDR arrays
  int SCTR_cluster_bit[numChan];   
  
  for(int trig_j=0;trig_j<numChan;trig_j++) SCTR_cluster_bit[trig_j]=0;
  
  for(int trig_i = trig_search_init; trig_i < max_total_bin; trig_i++) { // scan the different window positions

    // for trigger check:
    int N_pass = 0;
    int N_pass_V = 0;
    int N_pass_H = 0;
    int window_pass_bit = 0; // Whether this trig_i window passes

    check_TDR_configuration=0;
      
    for (int trig_j=0; trig_j<numChan; trig_j++){
      
      int string_i = detector->getStringfromArbAntID( i, trig_j);
      int antenna_i = detector->getAntennafromArbAntID( i, trig_j);

        trig_window_bin = (int)(settings1->TRIG_WINDOW / settings1->TIMESTEP);  // coincidence window bin for trigger
      
      if (settings1->TRIG_ONLY_BH_ON==0
	               || 
         (settings1->TRIG_ONLY_BH_ON==1 && settings1->DETECTOR==3 && detector->stations[i].strings[string_i].antennas[antenna_i].DAQchan==0)
                       ||
	 (settings1->TRIG_ONLY_LOW_CH_ON==1 && settings1->DETECTOR!=3 && antenna_i<2 ) ){ // channel filter: choose if to use lower/borehole channels or not
      
      int channel_num = detector->GetChannelfromStringAntenna ( i, string_i, antenna_i, settings1 );

      // assign Pthresh a value 
      double diode_noise_RMS = trigger->GetAntNoise_diodeRMS(channel_num-1, settings1);
      Pthresh_value[trig_j]=trigger->Full_window[trig_j][trig_i]/(diode_noise_RMS * detector->GetThresOffset( i, channel_num-1,settings1) );

	// this is to count how many local trigger clusters there are 
      if(Pthresh_value[trig_j]<powerthreshold){
	
	if(SCTR_cluster_bit[trig_j]==0) stations[i].strings[string_i].antennas[antenna_i].SingleChannelTriggers++;
		
	// records all the different Pthresh values that caused local trigger.
        if(settings1->TRIG_SCAN_MODE>2){ 
	  
	  if(SCTR_cluster_bit[trig_j]==0){// if first trigger in cluster
	    
	    stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.push_back(Pthresh_value[trig_j]);
	    
	    
	  }
	  else{// choose the highest trigger value (most negative) in cluster
	   
	    if(Pthresh_value[trig_j]<stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.back()) stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.back()=Pthresh_value[trig_j];
	    
	  }
	  
	}// trig scan mode > 2
	
	SCTR_cluster_bit[trig_j]=1;
	
      }// if local trigger
      else SCTR_cluster_bit[trig_j]=0;// if no local trigger, set zero to start a new cluster at next local trigger
      
      // and how many bins scanned
      stations[i].strings[string_i].antennas[antenna_i].TotalBinsScannedPerChannel++;
      


      
      						    
      // fill the buffers (if any changes occur mark check_TDR_configuration as non-zero)
      if(trig_i<trig_search_init+trig_window_bin) check_TDR_configuration+=buffer[trig_j]->fill(Pthresh_value[trig_j]);
      else check_TDR_configuration+=buffer[trig_j]->add(Pthresh_value[trig_j]);
	
      if(buffer[trig_j]->addToNPass>0){// if there is at least one value above threshold in the buffer, this is ++
	  
	N_pass++;
	if(detector->stations[i].strings[string_i].antennas[antenna_i].type == 0) N_pass_V++;
	if(detector->stations[i].strings[string_i].antennas[antenna_i].type == 1) N_pass_H++;

	  
// 	  if(last_trig_bin<trig_i) last_trig_bin=trig_i;// not sure we need this variable in this mode... 
	  
       }// if addToNPass>0
	
      }// non-testbed case (i.e. use all channels)
      
    }// for trig_j (channel scan)

    // check if global trigger...
    
    if( (settings1->TRIG_MODE==0&&( N_pass >= settings1->N_TRIG )) || 
	(settings1->TRIG_MODE==1&&( N_pass_V >= settings1->N_TRIG_V || 
				    N_pass_H >= settings1->N_TRIG_H )
	 ) 
	//|| (settings1->TRIG_MODE==2&&( N_pass_0 >= settings1->N_TRIG_0 || N_pass_1 >= settings1->N_TRIG_1 ))
	){ // if there's a trigger !
      

      global_pass_bit=1;  
      window_pass_bit=1;
      if(first_trigger==0){ // if this is the first trigger, mark this position and save event

	first_trigger=1;

	for(int trig_j=0;trig_j<numChan;trig_j++){
	 
	  int string_i = detector->getStringfromArbAntID( i, trig_j);
	  int antenna_i = detector->getAntennafromArbAntID( i, trig_j);

// 	  if(buffer[trig_j]->addToNPass>0) stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = trig_i-buffer[trig_j]->numBinsToLatestTrigger(); // mark the bin on which we triggered...
	  if(buffer[trig_j]->addToNPass>0) stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = trig_i-buffer[trig_j]->numBinsToOldestTrigger(); // mark the bin on which we triggered...
	  else stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = 0.;
	  
	}// for trig_j
	
	saveTriggeredEvent(settings1, detector, event, trigger, stationID, trig_search_init, max_total_bin, trig_window_bin, trig_i);
// 	cout<<"\nPthresh value=";
// 	if(scan_mode==1) for(int trig_j=0;trig_j<numChan; trig_j++) cout<<" "<<Pthresh_value[trig_j];
// 	if(scan_mode>1)  for(int trig_j=0;trig_j<numChan; trig_j++) cout<<" "<<buffer[trig_j]->best_value;
// 	cout<<"\n";
	
      }// first trigger
      
//       if(scan_mode==1) return last_trig_bin; //  if we aren't going to scan all the Pthresh values, just return
      if(scan_mode==1) return trig_i; //  if we aren't going to scan all the Pthresh values, just return
    }
    
    // if there's a trigger and anything changes in the buffers, restock the TDR arrays
    if( scan_mode>1 && check_TDR_configuration && window_pass_bit ){
      	    
      for(int trig_j=0;trig_j<numChan;trig_j++) TDR_all_sorted_temp[trig_j]=0;
      for(int trig_j=0;trig_j<numChanVpol;trig_j++) TDR_Vpol_sorted_temp[trig_j]=0;
      for(int trig_j=0;trig_j<numChanHpol;trig_j++) TDR_Hpol_sorted_temp[trig_j]=0;

      for(int trig_j=0;trig_j<numChan;trig_j++){// fill the TDR (unsorted) arrays if they improved... 
	
	if(buffer[trig_j]->best_value<stations[i].TDR_all[trig_j]) stations[i].TDR_all[trig_j]=buffer[trig_j]->best_value;

      }// for trig_j

      if(settings1->TRIG_MODE==0){ // for N out of all mode
      
	if(N_pass>=settings1->N_TRIG) for(int ii=0;ii<N_pass; ii++){// find the N_pass best channel's TDR and store them.
	 
	  double best_thresh=0;
	  int best_chan=0;
	  
	  for(int trig_j=0;trig_j<numChan;trig_j++) if(buffer[trig_j]->temp_value<best_thresh){ best_thresh=buffer[trig_j]->temp_value; best_chan=trig_j;}
	  
	  buffer[best_chan]->temp_value=0;
	  
	  TDR_all_sorted_temp[ii]=best_thresh;
	  
	}// for ii
	
	
	// debug output:
	if(TDR_all_sorted_temp[0]>TDR_all_sorted_temp[1]||TDR_all_sorted_temp[1]>TDR_all_sorted_temp[2]){
	   
	  cout<<"\n";
	  for(int p=0;p<80;p++) cout<<"*";
	  cout<<"\n  ordering problem: "<<TDR_all_sorted_temp[0]<<" "<<TDR_all_sorted_temp[1]<<" "<<TDR_all_sorted_temp[2]<<"\n";
	  for(int p=0;p<80;p++) cout<<"*";
	  cout<<"\n";
	    
	}
	
      }// if trig_mode==0
      if(settings1->TRIG_MODE==1){ // for N out of either polarization

	// for Vpol only:
	if(N_pass_V>=settings1->N_TRIG_V) for(int ii=0;ii<N_pass_V; ii++){// find the N_pass best channel's TDR and store them.

	  double best_thresh=0;
	  int best_chan=0;
	  
	  for(int trig_j=0;trig_j<numChan;trig_j++){
	    
	    int string_i = detector->getStringfromArbAntID( i, trig_j);
	    int antenna_i = detector->getAntennafromArbAntID( i, trig_j);  

	    if(detector->stations[i].strings[string_i].antennas[antenna_i].type == 0&&buffer[trig_j]->temp_value<best_thresh){ 
	      
	      best_thresh=buffer[trig_j]->temp_value; 
	      best_chan=trig_j;
	      
	    }// if best
	  }// for trig_j
	  buffer[best_chan]->temp_value=0;
	  
	  TDR_Vpol_sorted_temp[ii]=best_thresh;
	  
	}// for ii
	  
	  // debug output:
// 	  if(TDR_Vpol_sorted_temp[0]>TDR_Vpol_sorted_temp[1]||TDR_Vpol_sorted_temp[1]>TDR_Vpol_sorted_temp[2]){
// 	   
// 	    cout<<"\n";
// 	    for(int p=0;p<80;p++) cout<<"*";
// 	    cout<<"\n  ordering problem, Vpol: "<<TDR_Vpol_sorted_temp[0]<<" "<<TDR_Vpol_sorted_temp[1]<<" "<<TDR_Vpol_sorted_temp[2]<<"\n";
// 	    for(int p=0;p<80;p++) cout<<"*";
// 	    cout<<"\n";
// 	    
// 	  }
	// for Hpol only
	if(N_pass_H>=settings1->N_TRIG_H) for(int ii=0;ii<N_pass_H; ii++){// find the N_pass best channel's TDR and store them.

	  double best_thresh=0;
	  int best_chan=0;
	  
	  for(int trig_j=0;trig_j<numChan;trig_j++){ 
	    
	    int string_i = detector->getStringfromArbAntID( i, trig_j);
	    int antenna_i = detector->getAntennafromArbAntID( i, trig_j);
	 
	    if(detector->stations[i].strings[string_i].antennas[antenna_i].type==1&&buffer[trig_j]->temp_value<best_thresh){ 
	      
	      best_thresh=buffer[trig_j]->temp_value; 
	      best_chan=trig_j;
	      
	    }// if best
	  }// for trig_j
	
	  buffer[best_chan]->temp_value=0;
	  
	  TDR_Hpol_sorted_temp[ii]=best_thresh;
	  

	  	  
	}// for ii
	
//           // debug output:
// 	  if(TDR_Hpol_sorted_temp[0]>TDR_Hpol_sorted_temp[1]||TDR_Hpol_sorted_temp[1]>TDR_Hpol_sorted_temp[2]){
// 	   
// 	    cout<<"\n";
// 	    for(int p=0;p<80;p++) cout<<"*";
// 	    cout<<"\n  ordering problem, Hpol: "<<TDR_Hpol_sorted_temp[0]<<" "<<TDR_Hpol_sorted_temp[1]<<" "<<TDR_Hpol_sorted_temp[2]<<"\n";
// 	    for(int p=0;p<80;p++) cout<<"*";
// 	    cout<<"\n";
// 	    
// 	  }
      }// if trig_mode==1

      /*          
      if(settings1->TRIG_MODE==2){ // for N out of arbitrary sets of antennas

	// for Set 0 only:
	if(N_pass_0>=settings1->N_TRIG_0) for(int ii=0;ii<N_pass_0; ii++){// find the N_pass best channel's TDR and store them.

	  double best_thresh=0;
	  int best_chan=0;
	  
	  for(int trig_j=0;trig_j<numChan;trig_j++){
	    
	    int string_i = detector->getStringfromArbAntID( i, trig_j);
	    int antenna_i = detector->getAntennafromArbAntID( i, trig_j);	 

	    if((antenna_i == 0 || antenna_i == 2 ) && buffer[trig_j]->temp_value<best_thresh){ 
	      
	      best_thresh=buffer[trig_j]->temp_value; 
	      best_chan=trig_j;
	      
	    }// if best
	  }// for trig_j
	  buffer[best_chan]->temp_value=0;
	  
	  TDR_Vpol_sorted_temp[ii]=best_thresh;
	  
	}// for ii
	  
	  // debug output:
// 	  if(TDR_Vpol_sorted_temp[0]>TDR_Vpol_sorted_temp[1]||TDR_Vpol_sorted_temp[1]>TDR_Vpol_sorted_temp[2]){
// 	   
// 	    cout<<"\n";
// 	    for(int p=0;p<80;p++) cout<<"*";
// 	    cout<<"\n  ordering problem, Vpol: "<<TDR_Vpol_sorted_temp[0]<<" "<<TDR_Vpol_sorted_temp[1]<<" "<<TDR_Vpol_sorted_temp[2]<<"\n";
// 	    for(int p=0;p<80;p++) cout<<"*";
// 	    cout<<"\n";
// 	    
// 	  }
	// for Hpol only
	if(N_pass_1>=settings1->N_TRIG_1) for(int ii=0;ii<N_pass_1; ii++){// find the N_pass best channel's TDR and store them.

	  double best_thresh=0;
	  int best_chan=0;
	  
	  for(int trig_j=0;trig_j<numChan;trig_j++){ 
	    
	    int string_i = detector->getStringfromArbAntID( i, trig_j);
	    int antenna_i = detector->getAntennafromArbAntID( i, trig_j);
	 
	    if((antenna_i==1 || antenna_i == 3) && buffer[trig_j]->temp_value<best_thresh){ 
	      
	      best_thresh=buffer[trig_j]->temp_value; 
	      best_chan=trig_j;
	      
	    }// if best
	  }// for trig_j
	
	  buffer[best_chan]->temp_value=0;
	  
	  TDR_Hpol_sorted_temp[ii]=best_thresh;
	  

	  	  
	}// for ii
*/

    
    // check if temp TDR arrays improved, if so update TDR arrays:
    if(settings1->TRIG_MODE==0){
      
      if(N_pass>=settings1->N_TRIG) for(int ii=0;ii<N_pass;ii++) if(TDR_all_sorted_temp[ii]<stations[i].TDR_all_sorted[ii]) stations[i].TDR_all_sorted[ii]=TDR_all_sorted_temp[ii];
      
    }
    if(settings1->TRIG_MODE==1){
     
      if(N_pass_V>=settings1->N_TRIG_V) for(int ii=0;ii<N_pass_V;ii++) if(TDR_Vpol_sorted_temp[ii]<stations[i].TDR_Vpol_sorted[ii]) stations[i].TDR_Vpol_sorted[ii]=TDR_Vpol_sorted_temp[ii];
      if(N_pass_H>=settings1->N_TRIG_H) for(int ii=0;ii<N_pass_H;ii++) if(TDR_Hpol_sorted_temp[ii]<stations[i].TDR_Hpol_sorted[ii]) stations[i].TDR_Hpol_sorted[ii]=TDR_Hpol_sorted_temp[ii];
      
      // for this mode, can get TDR_all_sorted from these two arrays:
    }

    /*
    if(settings1->TRIG_MODE==2){
     
      if(N_pass_0>=settings1->N_TRIG_0) for(int ii=0;ii<N_pass_0;ii++) if(TDR_Vpol_sorted_temp[ii]<stations[i].TDR_Vpol_sorted[ii]) stations[i].TDR_Vpol_sorted[ii]=TDR_Vpol_sorted_temp[ii];
      if(N_pass_1>=settings1->N_TRIG_1) for(int ii=0;ii<N_pass_1;ii++) if(TDR_Hpol_sorted_temp[ii]<stations[i].TDR_Hpol_sorted[ii]) stations[i].TDR_Hpol_sorted[ii]=TDR_Hpol_sorted_temp[ii];
      
      // for this mode, can get TDR_all_sorted from these two arrays:
    }
    */

    
    }// if trigger and buffer changed
    
  }// while trig_i
  
  if(scan_mode>1&&global_pass_bit){
    
    if(settings1->TRIG_MODE==0){
     
      cout<<"\nPthresh best: ";
      for(int ii=0;ii<3;ii++) cout<<" "<<stations[i].TDR_all_sorted[ii];
      cout<<"\n";
        
      // debug output:
      if(stations[i].TDR_all_sorted[0]>stations[i].TDR_all_sorted[1]||stations[i].TDR_all_sorted[1]>stations[i].TDR_all_sorted[2]){
	   
	cout<<"\n";
	for(int p=0;p<80;p++) cout<<"*";
	cout<<"\n  ordering problem: "<<stations[i].TDR_all_sorted[0]<<" "<<stations[i].TDR_all_sorted[1]<<" "<<stations[i].TDR_all_sorted[2]<<"\n";
	for(int p=0;p<80;p++) cout<<"*";
	cout<<"\n";
		
      }// ordering problem
      
   
    }// trig mode 0
    
    
    if(settings1->TRIG_MODE==1){
      cout<<"\nPthresh best: ";
      cout<<"  Vpol: "; for(int ii=0;ii<stations[i].TDR_Vpol_sorted.size();ii++) cout<<" "<<stations[i].TDR_Vpol_sorted[ii];
      cout<<"  Hpol: "; for(int ii=0;ii<stations[i].TDR_Hpol_sorted.size();ii++) cout<<" "<<stations[i].TDR_Hpol_sorted[ii];
      cout<<"\n";
      
        
      // debug output:
      if(stations[i].TDR_Vpol_sorted[0]>stations[i].TDR_Vpol_sorted[1]||stations[i].TDR_Vpol_sorted[1]>stations[i].TDR_Vpol_sorted[2]){
	   
	cout<<"\n";
	for(int p=0;p<80;p++) cout<<"*";
	cout<<"\n  ordering problem (final) Vpol: "<<stations[i].TDR_Vpol_sorted[0]<<" "<<stations[i].TDR_Vpol_sorted[1]<<" "<<stations[i].TDR_Vpol_sorted[2]<<"\n";
	for(int p=0;p<80;p++) cout<<"*";
	cout<<"\n";
	
      }// ordering problem
	  
      // debug output:
      if(stations[i].TDR_Hpol_sorted[0]>stations[i].TDR_Hpol_sorted[1]||stations[i].TDR_Hpol_sorted[1]>stations[i].TDR_Hpol_sorted[2]){
	   
	cout<<"\n";
	for(int p=0;p<80;p++) cout<<"*";
	cout<<"\n  ordering problem (final), Hpol: "<<stations[i].TDR_Hpol_sorted[0]<<" "<<stations[i].TDR_Hpol_sorted[1]<<" "<<stations[i].TDR_Hpol_sorted[2]<<"\n";
	for(int p=0;p<80;p++) cout<<"*";
	cout<<"\n";
		
      }// ordering problem
         
    }// trig mode 1
    
  }// if scan_mode>1 and global pass
  

  for(int trig_j=0;trig_j<numChan; trig_j++) delete buffer[trig_j];
  delete [] buffer;
  
  if(scan_mode>1){
    
    delete [] TDR_all_sorted_temp;
    delete [] TDR_Vpol_sorted_temp;
    delete [] TDR_Hpol_sorted_temp;
  }
  
  
  return stations[i].Global_Pass;
  
}

int Report::saveTriggeredEvent(Settings *settings1, Detector *detector, Event *event, Trigger *trigger, int stationID, int trig_search_init, int max_total_bin, int trig_window_bin, int last_trig_bin){
 
  int i=stationID;
  int numChan=stations[i].TDR_all.size();
  cout<<"saving event"<<endl;
  stations[i].Global_Pass = last_trig_bin;

  int waveformLength = settings1->WAVEFORM_LENGTH;
  int waveformCenter = settings1->WAVEFORM_CENTER;



  for(int trig_j=0; trig_j<numChan;trig_j++){
    
    int string_i = detector->getStringfromArbAntID( i, trig_j);
    int antenna_i = detector->getAntennafromArbAntID( i, trig_j);
    
    if(stations[i].strings[string_i].antennas[antenna_i].Trig_Pass){// if this channel triggered
    
    stations[i].strings[string_i].antennas[antenna_i].Likely_Sol = -1; // no likely init
    
    int mindBin = 1.e9; // init big values
    int dBin = 0;
    
    for (int m=0; m<stations[i].strings[string_i].antennas[antenna_i].ray_sol_cnt; m++) {   // loop over raysol numbers

      if ( stations[i].strings[string_i].antennas[antenna_i].SignalExt[m] ) {

	dBin = abs( stations[i].strings[string_i].antennas[antenna_i].SignalBin[m] - stations[i].strings[string_i].antennas[antenna_i].Trig_Pass );
               
	if ( dBin < mindBin ) {
    
	  stations[i].strings[string_i].antennas[antenna_i].Likely_Sol = m; // store the ray sol number which is minimum difference between Trig_Pass bin
	  mindBin = dBin;
        	
	}
      }
      
      
    }// for m (ray sol numbers)
    
    
    
    
    
      if ( settings1->TRIG_ONLY_LOW_CH_ON==0 ) {

	cout<<endl<<"trigger passed at bin "<<stations[i].strings[string_i].antennas[antenna_i].Trig_Pass<<"  passed ch : "<<trig_j<<" ("<<detector->stations[i].strings[string_i].antennas[antenna_i].type<<"type) Direct dist btw posnu : "<<event->Nu_Interaction[0].posnu.Distance( detector->stations[i].strings[string_i].antennas[antenna_i] )<<" noiseID : "<<stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];
        
	if (stations[i].strings[string_i].antennas[antenna_i].Likely_Sol != -1) {
        
	  cout<<" ViewAngle : "<<stations[i].strings[string_i].antennas[antenna_i].view_ang[0]*DEGRAD<<" LikelyTrigSignal : "<<stations[i].strings[string_i].antennas[antenna_i].Likely_Sol;
        	  
	}
      }
      else if ( settings1->TRIG_ONLY_LOW_CH_ON==1 ) {
      
	cout<<endl<<"trigger passed at bin "<<stations[i].strings[string_i].antennas[antenna_i].Trig_Pass<<"  passed ant: str["<<string_i<<"].ant["<<antenna_i<<"] ("<<detector->stations[i].strings[string_i].antennas[antenna_i].type<<"type) Direct dist btw posnu : "<<event->Nu_Interaction[0].posnu.Distance( detector->stations[i].strings[string_i].antennas[antenna_i] )<<" noiseID : "<<stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];
        
	if (stations[i].strings[string_i].antennas[antenna_i].Likely_Sol != -1) {
        
	  cout<<" ViewAngle : "<<stations[i].strings[string_i].antennas[antenna_i].view_ang[0]*DEGRAD<<" LikelyTrigSignal : "<<stations[i].strings[string_i].antennas[antenna_i].Likely_Sol;
          	  
	}
        
      }
      
      

       
    }// if Trig_Pass
      
// now save the voltage waveform to V_mimic
      for (int mimicbin=0; mimicbin < waveformLength/2; mimicbin++) {

	// new DAQ waveform writing mechanism test
        if (settings1->V_MIMIC_MODE == 0) { // Global passed bin is the center of the window
        
	  stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back( ( trigger->Full_window_V[trig_j][ last_trig_bin - waveformLength/2 + waveformCenter + mimicbin ] )*1.e3 );// save in mV
          stations[i].strings[string_i].antennas[antenna_i].time.push_back(  last_trig_bin - waveformLength/2 + waveformCenter + mimicbin );
          stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back( ( waveformLength/2 + waveformCenter + mimicbin) * settings1->TIMESTEP*1.e9  );// save in ns
        }
        else if (settings1->V_MIMIC_MODE == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
          stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back( ( trigger->Full_window_V[trig_j][  last_trig_bin - (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin) - waveformLength/2 + waveformCenter + mimicbin ] )*1.e3 );// save in mV
          stations[i].strings[string_i].antennas[antenna_i].time.push_back(  last_trig_bin - (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin) - waveformLength/2 + waveformCenter + mimicbin );
          stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back( ( -(detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin) - waveformLength/2 + waveformCenter + mimicbin) * settings1->TIMESTEP*1.e9  );// save in ns
          	  
	}
          
          else if (settings1->V_MIMIC_MODE == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
          
           stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back( ( trigger->Full_window_V[trig_j][  last_trig_bin - (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) - waveformLength/2 + waveformCenter + mimicbin ] )*1.e3 );// save in mV
	   stations[i].strings[string_i].antennas[antenna_i].time.push_back(  last_trig_bin - (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) - waveformLength/2 + waveformCenter + mimicbin );
           stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back( ( -(detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) - waveformLength/2 + waveformCenter + mimicbin) * settings1->TIMESTEP*1.e9 + detector->params.TestBed_WFtime_offset_ns );// save in ns
           	    
	  }
      if (mimicbin == 0) {
        for (int m = 0; m < stations[i].strings[string_i].antennas[antenna_i].ray_sol_cnt; m++) { ///< calculates time of center of each rays signal based on readout window time config
            double signal_center_offset = (double)(stations[i].strings[string_i].antennas[antenna_i].SignalBin[m] - stations[i].strings[string_i].antennas[antenna_i].time[0]) * settings1->TIMESTEP * 1.e9;
            double signal_center_time = signal_center_offset + stations[i].strings[string_i].antennas[antenna_i].time_mimic[0];
            //! signal_center_offset: time offset between beginning of readout window and center of signal
            //! signal_center_time: time of center of signal based on readout window time config
            stations[i].strings[string_i].antennas[antenna_i].SignalBinTime.push_back(signal_center_time);
        }
      }      
  	
      }
      
            
      // set global_trig_bin values
       if (settings1->V_MIMIC_MODE == 0) { // Global passed bin is the center of the window
         stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = waveformLength/2 + waveformCenter ;
       }
       else if (settings1->V_MIMIC_MODE == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
         stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin) + waveformLength/2 + waveformCenter ;
       }
       else if (settings1->V_MIMIC_MODE == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
          stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) + waveformLength/2 + waveformCenter ;
       }
      
      
      double arrivtime = stations[i].strings[string_i].antennas[antenna_i].arrival_time[0];
      double X = detector->stations[i].strings[string_i].antennas[antenna_i].GetX();
      double Y = detector->stations[i].strings[string_i].antennas[antenna_i].GetY();
      double Z = detector->stations[i].strings[string_i].antennas[antenna_i].GetZ();
      
      
      stations[i].total_trig_search_bin = stations[i].Global_Pass + trig_window_bin - trig_search_init;
      
      
      
      
      
      
      
    
  }// for trig_j
  
  
   if(settings1->OUTPUT_TDR_GRAPH>0){
     
     settings1->OUTPUT_TDR_GRAPH--;
     
     TGraph **gr=new TGraph*[numChan];
     
     for(int trig_j=0;trig_j<numChan;trig_j++){
      
       int string_i = detector->getStringfromArbAntID( i, trig_j);
       int antenna_i = detector->getAntennafromArbAntID( i, trig_j);
       int channel_num = detector->GetChannelfromStringAntenna ( i, string_i, antenna_i, settings1 );
       double thresh_value=0;
       
      // assign Pthresh a value 
      double diode_noise_RMS = trigger->GetAntNoise_diodeRMS(channel_num-1, settings1);
      thresh_value=detector->GetThres(i, channel_num-1, settings1) * diode_noise_RMS * detector->GetThresOffset( i, channel_num-1,settings1);
       
       gr[trig_j]=new TGraph();
       for(int trig_i=0;trig_i<settings1->DATA_BIN_SIZE/2;trig_i++) gr[trig_j]->SetPoint(trig_i, (trig_search_init+trig_i), trigger->Full_window[trig_j][trig_i]);
       gr[trig_j]->SetNameTitle(Form("TDR_waveform%dC%02d",  settings1->OUTPUT_TDR_GRAPH, trig_j), Form("Tunnel diode response waveform %d, channel %02d, trig_pass= %d, P_{th}= %le; time bins; power after convolution with tunnel diode", settings1->OUTPUT_TDR_GRAPH, trig_j, stations[i].strings[string_i].antennas[antenna_i].Trig_Pass, thresh_value));
       gr[trig_j]->Write();
       delete gr[trig_j];
     }
     
     delete [] gr;
     
	
   }
  
  return 1;
  
  
  
}// saveTriggeredEvent

#ifdef ARA_UTIL_EXISTS
void Report::MakeUsefulEvent(Detector *detector, Settings *settings1, Trigger *trigger, int stationID, int stationIndex, UsefulIcrrStationEvent *theUsefulEvent) {
    if (stationID < detector->params.number_of_stations){
        int i = stationID;
        cout << stationID << endl;
	int ch_limit;
	if (stationID == 0){
	  ch_limit = 14;
	} else {
	  ch_limit = 16;
	}

        for (int ch_loop=0; ch_loop<ch_limit; ch_loop++) {
            int string_i = detector->getStringfromArbAntID( stationIndex, ch_loop);
            int antenna_i = detector->getAntennafromArbAntID( stationIndex, ch_loop);
            int AraRootChannel = 0;
            AraRootChannel = detector->GetChannelfromStringAntenna (i, string_i, antenna_i, settings1);
            
            int UsefulEventBin;
            if ( settings1->NFOUR/2 < EFFECTIVE_LAB3_SAMPLES*2) UsefulEventBin = settings1->NFOUR/2;
            else UsefulEventBin = EFFECTIVE_LAB3_SAMPLES*2;
            
            //for (int mimicbin=0; mimicbin<settings1->NFOUR/2; mimicbin++) {
            for (int mimicbin=0; mimicbin<UsefulEventBin; mimicbin++) {
                if (stations[stationIndex].Global_Pass > 0){
                    theUsefulEvent->fVoltsRF[AraRootChannel-1][mimicbin] = stations[i].strings[string_i].antennas[antenna_i].V_mimic[mimicbin];
                    theUsefulEvent->fTimesRF[AraRootChannel-1][mimicbin] = stations[i].strings[string_i].antennas[antenna_i].time_mimic[mimicbin];
                }
                else {
                    theUsefulEvent->fVoltsRF[AraRootChannel-1][mimicbin] = 0.;
                    theUsefulEvent->fTimesRF[AraRootChannel-1][mimicbin] = 0.;
                }
            }
            //theUsefulEvent->fNumPointsRF[ch_loop] = EFFECTIVE_SAMPLES * 2;
            theUsefulEvent->fNumPointsRF[ch_loop] = UsefulEventBin;
        }
    }
}
#endif

#ifdef ARA_UTIL_EXISTS
void Report::MakeUsefulEvent(Detector *detector, Settings *settings1, Trigger *trigger, int stationID, int stationIndex, UsefulAtriStationEvent *theUsefulEvent) {
  //    if (stationID < detector->params.number_of_stations){

        int i = stationID;
	int stationID_AraRoot = settings1->DETECTOR_STATION_ARAROOT;
	cout << "StationID: " << stationID << endl;	
	cout << "StationID_AraRoot: " << stationID_AraRoot << endl;
	theUsefulEvent->fNumChannels = 32;
	theUsefulEvent->stationId = stationID_AraRoot;
	
	//	cout << endl << stationID << endl;
	
	int ch_limit;
	if (stationID == 0){
	  ch_limit = 14;
	} else {
	  ch_limit = 16;
	}

	int maxElecChans = 32;
	
	for (int ch_loop=0; ch_loop < ch_limit; ch_loop++) {
	  //	  int elecChan = AraGeom->getElecChanFromRFChan(ch_loop, stationID);
	  int elecChan = AraGeomTool::Instance()->getElecChanFromRFChan(ch_loop, stationID_AraRoot);
	  int string_i = 0;
	  int antenna_i = 0;
	  detector->GetSSAfromChannel(stationID, ch_loop, &antenna_i, &string_i, settings1);

	  //	  cout << ch_loop << " : " << elecChan << " : " << string_i << " : " << antenna_i << endl;

	  //	  int string_i = detector->getStringfromArbAntID( stationIndex, ch_loop);
	  //	  int antenna_i = detector->getAntennafromArbAntID( stationIndex, ch_loop);
	  int AraRootChannel = 0;
	  AraRootChannel = detector->GetChannelfromStringAntenna (stationID, string_i, antenna_i, settings1);
	  
	  int UsefulEventBin;
	  //	  if ( settings1->NFOUR/2 < EFFECTIVE_LAB3_SAMPLES*2) UsefulEventBin = settings1->NFOUR/2;
	  //	  else UsefulEventBin = EFFECTIVE_LAB3_SAMPLES*2;

	  UsefulEventBin = settings1->WAVEFORM_LENGTH;

	  vector < double > volts;
	  volts.resize(UsefulEventBin);
	  vector < double > times;
	  times.resize(UsefulEventBin);

	  //	  theUsefulEvent->fVolts[AraRootChannel-1].resize(UsefulEventBin);
	  //	  theUsefulEvent->fTimes[AraRootChannel-1].resize(UsefulEventBin);
	  //	  cout << string_i << " : " << antenna_i << endl;
	  
	  //for (int mimicbin=0; mimicbin<settings1->NFOUR/2; mimicbin++) {
	  for (int mimicbin=0; mimicbin<UsefulEventBin; mimicbin++) {
	    if (stations[stationIndex].Global_Pass > 0){
	      //	      cout << "Test 1" << endl;
	      volts[mimicbin] = stations[stationIndex].strings[string_i].antennas[antenna_i].V_mimic[mimicbin];
	      times[mimicbin] = stations[stationIndex].strings[string_i].antennas[antenna_i].time_mimic[mimicbin];
	      
	      
	    }
	    else {
	      volts[mimicbin] = 0.;
	      times[mimicbin] = 0.;
	    }
	  }
	  theUsefulEvent->fVolts.insert( std::pair < int, std::vector < double > > (elecChan, volts ));
	  theUsefulEvent->fTimes.insert( std::pair < int, std::vector < double > > (elecChan, times ));

	  //	  cout << elecChan << " : " << theUsefulEvent->fTimes[elecChan][1] <<  " : " << stations[stationIndex].strings[string_i].antennas[antenna_i].time_mimic[1] << endl;
	  //	  cout << elecChan << " : " << theUsefulEvent->fVolts[elecChan][1] <<  " : " << stations[stationIndex].strings[string_i].antennas[antenna_i].V_mimic[1] << endl;


	  volts.clear();
	  times.clear();

        }
	//  }
}
#endif


void Report::ClearUselessfromConnect(Detector *detector, Settings *settings1, Trigger *trigger){
    
    for (int i = 0; i< detector->params.number_of_stations; i++) {
        // now remove all information which are useless
        for (int c_j=0; c_j< detector->stations[i].strings.size(); c_j++) {
            for (int c_k=0; c_k< detector->stations[i].strings[c_j].antennas.size(); c_k++) {
                stations[i].strings[c_j].antennas[c_k].clear_useless(settings1);  // clear data in antenna which stored in previous event
                //stations[i].strings[c_j].antennas[c_k].clear();  // clear data in antenna which stored in previous event
            }
        }
        //clear_useless(settings1);
    }
}



void Report::Convolve_Signals(
    Antenna_r *antenna, int channel_index, int station_number, 
    Settings *settings1, Trigger *trigger, Detector *detector
){
    // For the provided antenna: 
    // Calculate bin values for the signal and determine if how many rays fit in the readout window
    // Combine signals from rays and load noise waveforms
    // Pass the noise+signal waveform through the tunnel diode
    // Enforce voltage saturation

    // init : bin + maxt_diode_bin, fin : bin + NFOUR/2

    // Clear old signals
    signal_bin.clear();
    signal_dbin.clear();
    connect_signals.clear();

    // Determine the bin (index) where the signal will arrive in output waveforms
    //   and flag which signals will fit into the same time window
    for (int m = 0; m < antenna->ray_sol_cnt; m++)
    {
        // loop over raysol numbers

        // Store the bin where the singal is located
        signal_bin.push_back(
            (antenna->arrival_time[m] - stations[station_number].min_arrival_time) / (settings1->TIMESTEP) 
            + settings1->NFOUR *2 + trigger->maxt_diode_bin);
        antenna->SignalBin.push_back(signal_bin[m]);

        // Determine if any rays are connected to each other
        // If connect_signals[m] == 1, signal from the ray at index `m` is 
        //   readout in the same waveform as the signal from the ray at index `m-1`.
        if (m > 0)
        {
            signal_dbin.push_back(signal_bin[m] - signal_bin[m - 1]);
            if (signal_dbin[m - 1] < settings1->NFOUR / 2)
            {
                // if two ray_sol time delay is smaller than time window
                connect_signals.push_back(1);
            }
            else
            {
                connect_signals.push_back(0);
            }
        }
        else if (antenna->ray_sol_cnt == 1 && m == 0)
        {
            // if there's only one solution
            connect_signals.push_back(0);
        }

    }

    // Set default signal wf length
    int BINSIZE = settings1->NFOUR/2;

    // If there are no ray solutions to this antenna, fill all the arrays with 0s
    //   then get noise and convolve
    if ( antenna->ray_sol_cnt == 0 ){
        
        // Save signal wf as array of 0s
        vector <double> V_signal;
        for (int bin=0; bin<BINSIZE; bin++) V_signal.push_back(0.);
        for (int bin=0; bin<BINSIZE; bin++) antenna->V_convolved.push_back(0.);
        for (int bin=0; bin<BINSIZE; bin++) antenna->V_noise.push_back(0.);

        // Convolve the noise signal and add to the array used for triggering
        GetNoiseThenConvolve(
            antenna, V_signal,
            BINSIZE, BINSIZE/2, 0, // waveform_length, this_signalbin, n_connected_rays
            channel_index, station_number, 
            settings1, trigger, detector);
            
    }
    // If there are ray solutions, prepare signal wf array with 0s 
    //   and save an array with the full noise-only waveform to the antenna
    else {

        // Initialize waveform length as 2 BINSIZES longer than the last ray's signal bin
        int array_length = signal_bin[antenna->ray_sol_cnt-1] + 2*BINSIZE;

        // Make array of 0s for array we're saving all ray signals to
        for (int i=0; i<array_length; i++) antenna->V_convolved.push_back(0.);
        for (int i=0; i<array_length; i++) antenna->V_noise.push_back(0.);

    }
    
    // Loop over ray solutions and get signals from each
    // Loop does not run if there are no ray solutions
    for (int m = 0; m < antenna->ray_sol_cnt; m++)
    {

        // Grab signal-only waveform, what the signal bin of the main ray is, 
        //   and the number of rays in this tunnel diode window (connected rays)
        // V_signal will have a length of BINSIZE if only 1 ray exists in a waveform
        //   array with length NFOUR/2 and BINSIZE*2 if 2 or more fit into one tunnel diode window.
        int this_signalbin = 0;
        int n_connected_rays = 0;
        vector <double> V_signal;
        GetAntennaSignalWF(
            m, &n_connected_rays, &this_signalbin, BINSIZE, antenna, &V_signal,
            settings1, trigger, detector);

        // Get the noise-only wf from just this signal-window then convolve 
        //   through the tunnel diode
        if (n_connected_rays>0) {
            GetNoiseThenConvolve(
                antenna, V_signal,
                BINSIZE, this_signalbin, n_connected_rays, 
                channel_index, station_number, 
                settings1, trigger, detector);
        }

    }   // end loop over ray solutions

}


void Report::GetAntennaSignalWF(
    int raysol, int *n_connected_rays, int *this_signalbin,
    int BINSIZE, Antenna_r *antenna, vector <double> *V_signal,
    Settings *settings1, Trigger *trigger, Detector *detector
){
    // Determine which rays fall within the same trigger bin (BINSIZE) 
    //   and combine them if so

    // loop over raysol numbers
    // when ray_sol_cnt == 0, this loop inside codes will not run

    if (raysol == 0)
    {
        // if it's first sol

        if (connect_signals[raysol] == 1)
        {
            *n_connected_rays = 2;
            *this_signalbin = signal_bin[raysol];

            // Combine signals from m and m+1 ray in double sized (NFOUR) signal array
            Select_Wave_Convlv_Exchange( 
                signal_bin[raysol], signal_bin[raysol + 1], 
                antenna->V[raysol], antenna->V[raysol + 1], 
                BINSIZE, V_signal);
        }
        else if (connect_signals[raysol] == 0)
        {
            *n_connected_rays = 1;
            *this_signalbin = signal_bin[raysol];

            // Make a normally sized (NFOUR/2) signal array with just this ray 
            Select_Wave_Convlv_Exchange(
                antenna->V[raysol], 
                BINSIZE, V_signal);
        }
    }
    else
    {
        // if it's not the first sol

        if (raysol + 1 < antenna->ray_sol_cnt)
        {
            // if there is next raysol

            if (connect_signals[raysol] == 1)
            {
                // next raysol is connected

                if (connect_signals[raysol - 1] == 1)
                {
                    // and previous raysol also connected

                    *n_connected_rays = 3;
                    *this_signalbin = signal_bin[raysol];

                    // Combine signals from all 3 rays in a double sized (NFOUR) signal array
                    Select_Wave_Convlv_Exchange(
                        signal_bin[raysol - 1], signal_bin[raysol], signal_bin[raysol + 1], 
                        antenna->V[raysol - 1], antenna->V[raysol], antenna->V[raysol + 1], 
                        BINSIZE, V_signal);
                }
                else if (connect_signals[raysol - 1] == 0)
                {
                    // and previous raysol not connected

                    *n_connected_rays = 2;
                    *this_signalbin = signal_bin[raysol];

                    // Combine signals from m and m+1 ray in a double sized (NFOUR) signal array
                    Select_Wave_Convlv_Exchange(
                        signal_bin[raysol], signal_bin[raysol + 1], 
                        antenna->V[raysol], antenna->V[raysol + 1], 
                        BINSIZE, V_signal);
                }
            }
            else if (connect_signals[raysol] == 0)
            {
                // next raysol is not connected

                if (connect_signals[raysol - 1] == 1)
                {
                    // and previous raysol is connected

                    // skip the process as this should have done before

                }
                else if (connect_signals[raysol - 1] == 0)
                {
                    // and previous raysol not connected

                    *n_connected_rays = 1;
                    *this_signalbin = signal_bin[raysol];

                    // Make a normally sized (NFOUR/2) signal array with just this ray 
                    Select_Wave_Convlv_Exchange(
                        antenna->V[raysol], 
                        BINSIZE, V_signal);
                }
            }
        }
        else
        {
            // there is no next raysol (this "m" is the last raysol)

            if (connect_signals[raysol - 1] == 1)
            {
                // and previous raysol is connected

                // skip the process as this should have done before

            }
            else if (connect_signals[raysol - 1] == 0)
            {
                // and previous raysol is not connected

                *n_connected_rays = 1;
                *this_signalbin = signal_bin[raysol];

                // Make a normally sized (NFOUR/2) signal array with just this ray 
                Select_Wave_Convlv_Exchange(
                    antenna->V[raysol], 
                    BINSIZE, V_signal);
            }
        }

    }   // if not the first raysol (all other raysols)

}


void Report::Select_Wave_Convlv_Exchange(
    vector <double> &V, 
    int BINSIZE, vector <double> *V_signal
) {
    // Convolve a single signal into the signal array
    
    // Clear previous waveform data
    V_signal->clear();
    
    // Save the signal voltage waveform
    for (int bin=0; bin<BINSIZE; bin++) {
        V_signal->push_back( V[bin] );
    }

}


void Report::Select_Wave_Convlv_Exchange(
    int signalbin_1, int signalbin_2, 
    vector <double> &V1, vector <double> &V2, 
    int BINSIZE, vector <double> *V_signal
) {
    // Convolve 2 signals into one signal array
    
    // Clear previous waveform data
    V_signal->clear();
    
    // Save the signal voltage waveform
    int signal_dbin = signalbin_2 - signalbin_1;
    for (int bin=0; bin<BINSIZE*2; bin++) {

        if (bin < signal_dbin) {  // bins where only first signal is shown
            V_signal->push_back( V1[bin] );
        }
        else if (bin < BINSIZE) { // bins where first + second signal is shown
            V_signal->push_back( V1[bin] + V2[bin - signal_dbin] );
        }
        else if (bin < BINSIZE + signal_dbin) { // bins where only second signal is shown
            V_signal->push_back( V2[bin - signal_dbin] );
        }
        else {
            V_signal->push_back(0.);
        }

    }

}


void Report::Select_Wave_Convlv_Exchange(
    int signalbin_0, int signalbin_1, int signalbin_2, 
    vector <double> &V0, vector <double> &V1, vector <double> &V2, 
    int BINSIZE, vector <double> *V_signal
) {
    // Convolve 3 signals into one signal array
    
    // Clear previous waveform data  
    V_signal->clear();
    
    // Save the signal voltage waveform
    int signal_dbin = signalbin_2 - signalbin_1;
    int signal_dbin0 = signalbin_1 - signalbin_0;
    for (int bin=0; bin<BINSIZE*2; bin++) {

        if (bin < signal_dbin) {  // bins where no second signal is shown
            if ( signal_dbin0 + bin < BINSIZE ) {   // previous signal is also here!
                V_signal->push_back( V0[ signal_dbin0 + bin] );
            }
            else {  // no previous signal, and next signal
                V_signal->push_back( V1[bin] );
            }
        }
        else if (bin < BINSIZE) { // bins where first + second signal is shown
            if ( signal_dbin0 + bin < BINSIZE ) {   // previous signal is also here!
                V_signal->push_back( V1[bin] + V0[ signal_dbin0 + bin] + V2[bin - signal_dbin] );
            }
            else {  // no previous signal, and next signal
                V_signal->push_back( V1[bin] + V2[bin - signal_dbin] );
            }
        }
        else if (bin < BINSIZE + signal_dbin) { // bins where only second signal is shown
            V_signal->push_back( V2[bin - signal_dbin] );
        }
        else {
            V_signal->push_back(0.);
        }

    }

}


void Report::GetNoiseThenConvolve(
    Antenna_r *antenna, vector <double> V_signal,
    int BINSIZE, int this_signalbin, int n_connected_rays, 
    int channel_index, int station_number, 
    Settings *settings1, Trigger *trigger, Detector *detector
){
    // Get noise waveform, signal waveform, combine them, 
    //   convolve them through the tunnel diode, apply voltage saturation,
    //   then save the noise and signals to the 
    //   `Antenna_r` object and the `trigger` class
    
    // Extend the length of this waveform we're constructing if more than 1 ray connected
    int wf_length = 0;
    int min_wf_bin = this_signalbin-BINSIZE/2+(trigger->maxt_diode_bin);
    int max_wf_bin = 0;
    vector <double> diode_response;
    if ( n_connected_rays > 1 ) { // multiple ray solutions in one window
        wf_length = BINSIZE * 2;
        max_wf_bin = this_signalbin + BINSIZE/2 + BINSIZE;
        diode_response = detector->fdiode_real_double;
    }
    else if ( antenna->ray_sol_cnt == 0 ){ // No rays connected to this antenna
        this_signalbin = BINSIZE/2;
        wf_length = BINSIZE;
        min_wf_bin = 0;
        max_wf_bin = BINSIZE;
        diode_response = detector->fdiode_real;
    }
    else { // Only one ray signal in the window
        wf_length = BINSIZE;
        max_wf_bin = this_signalbin + BINSIZE/2;
        diode_response = detector->fdiode_real;
    }

    // Get noise-only waveform
    vector <double> V_noise;
    GetAntennaNoiseWF(
        this_signalbin, wf_length, BINSIZE, channel_index, station_number, &V_noise, 
        settings1, trigger, detector);

    // Create noise+signal waveforms
    V_total_forconvlv.clear();
    for (int bin=0; bin<wf_length; bin++){
        V_total_forconvlv.push_back( V_signal.at(bin) + V_noise[bin]);
    }

    // Push noise+signal waveform through the tunnel diode
    trigger->myconvlv( V_total_forconvlv, wf_length, diode_response, V_total_forconvlv);

    // Export our convolved waveforms to trigger->Full_window and trigger->Full_window_V
    for (int bin=min_wf_bin; bin<max_wf_bin; bin++) {

        // Export WFs with newly-processed signals
        trigger->Full_window[channel_index][bin] = V_total_forconvlv[bin - this_signalbin + BINSIZE/2];
        trigger->Full_window_V[channel_index][bin] += V_signal[bin - this_signalbin + BINSIZE/2];
        antenna->V_convolved[bin] += V_signal[bin - this_signalbin + BINSIZE/2];
        antenna->V_noise[bin] += V_noise[bin - this_signalbin + BINSIZE/2];

        // Add electronics saturation effect to WF we'll perform trigger check on
        if ( trigger->Full_window_V[channel_index][bin] > settings1->V_SATURATION ) {
            trigger->Full_window_V[channel_index][bin] = settings1->V_SATURATION;
        }
        else if ( trigger->Full_window_V[channel_index][bin] < -1.*settings1->V_SATURATION ) {
            trigger->Full_window_V[channel_index][bin] = -1.*settings1->V_SATURATION;
        }

    } 

}


void Report::GetAntennaNoiseWF(
    int signalbin, 
    int wf_length, int BINSIZE, int channel_index, int StationIndex, vector <double> *V_noise_only,
    Settings *settings1, Trigger *trigger, Detector *detector
){
    // Save a `wf_length` long noise waveform to the provided `V_noise_only` array
    //   with the signal bin located at index `BINSIZE/2`

    // Clear old noise waveforms
    V_noise_only->clear();

    // Pick noise waveform
    vector< vector <double> > *noise_wf;
    int channel_number = GetChNumFromArbChID(detector, channel_index, StationIndex, settings1) - 1;
    if ( settings1->NOISE_CHANNEL_MODE==0) {
        noise_wf = &trigger->v_noise_timedomain;
    }
    // Use channel-by-channel temperatures to generate noise
    else if ( settings1->NOISE_CHANNEL_MODE==1) {
        noise_wf = &trigger->v_noise_timedomain_ch[ channel_number ];
    }
    // Only use channel-by-channel temperatures to generate noise for the first 8 channels. 
    // Use the same temperature for the remaining 8.
    else if ( settings1->NOISE_CHANNEL_MODE==2) {
        // If this channel is one of the first 8, use channel specific noise
        if ( channel_number < 8) {
            noise_wf = &trigger->v_noise_timedomain_ch[ channel_number ];
        }
        // This channel is NOT one of the first 8, use same noise for these channels
        else {
            noise_wf = &trigger->v_noise_timedomain_ch[ 8 ];
        }
    }
    else{
        cout<<"Cannot find noise waveform for channel "<<channel_index;
        cout<<" on station "<<StationIndex;
        cout<<" for requested NOISE_CHANNEL_MODE "<<settings1->NOISE_CHANNEL_MODE<<endl;
        throw std::runtime_error("");
    }

    // Loop over bins and get noise voltage value for each
    int bin_value;
    int noise_ID_index;
    int noise_wf_index;
    for (int bin=0; bin<wf_length; bin++) {
        bin_value = signalbin - BINSIZE/2 + bin;
        noise_ID_index = bin_value / settings1->DATA_BIN_SIZE;
        noise_wf_index = bin_value % settings1->DATA_BIN_SIZE;
        V_noise_only->push_back(
            noise_wf->at( noise_ID[noise_ID_index] ).at( noise_wf_index )
        );
    } // end loop over bins

}


void Report::Apply_Gain_Offset(Settings *settings1, Trigger *trigger, Detector *detector, int ID, int StationIndex ) {

    int string_num = detector->getStringfromArbAntID( StationIndex, ID );
    int ant_num = detector->getAntennafromArbAntID( StationIndex, ID );

    //int channel_num = detector->GetChannelfromStringAntenna ( StationIndex, string_num, ant_num );
    int channel_num = detector->GetChannelfromStringAntenna ( StationIndex, string_num, ant_num, settings1 );
    //cout<<"station "<<StationIndex<<" ch"<<channel_num<<" applying gain offset "<<detector->GetGainOffset( StationIndex, channel_num, settings1 )<<" applying"<<endl;

    for (int bin=0; bin<settings1->DATA_BIN_SIZE; bin++) {   // test for full window
        
        trigger->Full_window[ID][bin] = ( trigger->Full_window[ID][bin] * detector->GetGainOffset( StationIndex, channel_num-1, settings1 ) * detector->GetGainOffset( StationIndex, channel_num-1, settings1 ) ); // offset in voltage factor so we need power (V^2 factor to diode response)
        trigger->Full_window_V[ID][bin] = ( trigger->Full_window_V[ID][bin] * detector->GetGainOffset( StationIndex, channel_num-1, settings1 ) ); // gain in voltage factor
    }

}



int Report::GetChNumFromArbChID( Detector *detector, int ID, int StationIndex, Settings *settings1 ) {

    int string_num = detector->getStringfromArbAntID( StationIndex, ID );
    int ant_num = detector->getAntennafromArbAntID( StationIndex, ID );
    
    int StationID = detector->stations[StationIndex].StationID;
    //    cout << "Station ID: " << StationID << endl;
    //    cout << "string_num: " << string_num << endl;
    //    cout << "ant_num: " << ant_num << endl;
    //    cout << StationID << " : " << string_num << " : " << ant_num << endl;
    
    //int channel_num = detector->GetChannelfromStringAntenna ( StationIndex, string_num, ant_num );
    int channel_num = detector->GetChannelfromStringAntenna ( StationID, string_num, ant_num, settings1 );

    return channel_num;
}



Vector Report::GetPolarization (Vector &nnu, Vector &launch_vector) {
    // copy from icemc GetPolarization

    // Want to find a unit vector in the same plane as
    // nnu and launch_vector, but perpendicular to launch_vector, pointing away
    // from nnu.
    
    // cross nnu with launch_vector to get the direction of the B field.
    Vector n_bfield = nnu.Cross(launch_vector);
    
    // cross b-field with nrf2_iceside to get the polarization vector.
    Vector n_pol = n_bfield.Cross(launch_vector);
    
    n_pol = n_pol.Unit();
    
    // check and make sure E-field is pointing in the right direction.
    if (nnu*launch_vector>0 && n_pol*nnu>0)
	cout << "error in GetPolarization.  \n";
    
     /*
      * BAC, AC, and JT 2020/11/3
      * 
      * We confirmed that we think this calculation is correct up to a minus sign
      * e.g. outward or inward on the cone depending on whether
      * one is outside or inside of the cherenkov cone
     */
    
    return n_pol;
} //GetPolarization


void Report::GetParameters( Position &src, Position &trg, Vector &nnu, double &viewangle, double receive_angle, Vector &launch_vector, Vector &receive_vector, Vector &n_trg_slappy, Vector &n_trg_pokey) {

    viewangle = PI/2. - viewangle;  // viewangle was actually launch angle

    launch_vector = (trg.Cross( trg.Cross(src) )).Rotate(viewangle, trg.Cross(src));
    launch_vector = launch_vector.Unit();
    viewangle = launch_vector.Angle(nnu);
    
     /*
      * BAC, AC, and JT 2020/11/3
      * 
      * n_trg_pokey points from the center of the earth to the antenna location
      * so for any antenna somewhere on the earth, this points to approximately
      * the center axis of the antenna
      *
      * n_trg_slappy is parallel to the surface and perpendicular to the plane
      * defined by the vector pointing to ant and vertex
     */
    
    //cout<<"launch_vector angle between R1 (trg) : "<<launch_vector.Angle(trg)<<"\n";

    receive_vector = trg.Rotate( receive_angle, src.Cross(trg) );
    receive_vector = receive_vector.Unit();

    n_trg_pokey = trg.Unit();
    n_trg_slappy = (trg.Cross(src)).Unit();


}


double Report::GaintoHeight(double gain, double freq, double n_medium, double Z_A) {
    
    return sqrt((gain*CLIGHT*CLIGHT*Zr) / (4*PI*n_medium*freq*freq*Z0));

}

double Report::calculatePolFactor(Vector &Pol_vector, int ant_type, double antenna_theta, double antenna_phi){
    // convert to radians
     antenna_phi*=(PI/180);
     antenna_theta*=(PI/180);

     // calculate the local thetaHat and phiHat vectors
     Vector thetaHat = Vector(cos(antenna_theta)*cos(antenna_phi),
                              cos(antenna_theta)*sin(antenna_phi),
                              -sin(antenna_theta));

     Vector phiHat = Vector(-sin(antenna_phi),
                            cos(antenna_phi),
                            0);
     double pol_factor=0.;
     if (ant_type == 0) {    // if v pol
         pol_factor = Pol_vector *thetaHat;
     }
     else if (ant_type == 1) {   // if h pol
         pol_factor = Pol_vector *phiHat;
     }
     pol_factor = pol_factor;
     return pol_factor;
 }


void Report::ApplyAntFactors(double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &pol_factor, double &vmmhz, double antenna_theta, double antenna_phi) {  // vmmhz is input and output. output will have some antenna factors on it

    //double pol_factor;
    pol_factor = calculatePolFactor(Pol_vector, ant_type, antenna_theta, antenna_phi);

    // apply 3dB spliter, d nu to prepare FFT
    // now actually vmmhz is not V/m/MHz but V/m/Hz unit
    //vmmhz = vmmhz/sqrt(2.)/(settings1->TIMESTEP*1.E6); //sqrt(2) for 3dB spliter for TURF, SURF
    vmmhz = vmmhz/sqrt(2.)/(1.E6); //sqrt(2) for 3dB spliter for TURF, SURF. 1/TIMESTEP moved to MakeArraysforFFT
    // 1/(1.E6) for V/MHz to V/Hz


    // apply antenna effective height and 0.5 (to calculate power with heff), and polarization factor
    // not vmmhz is actually V/Hz unit
    vmmhz = vmmhz * 0.5 * heff * pol_factor;

    // now we have to use MakeArraysforFFT to change signal arrays for FFT
}



void Report::ApplyAntFactors_Tdomain (double AntPhase, double heff, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img, Settings *settings1, double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode, bool applyInverse) {  
    /*  Report::ApplyAntFactors_Tdomain()
    Purpose: 
        Multiply (or divide) voltage in Fourier space by the antenna gain and phase.
    
    Inputs:
        AntPhase:  Phase of the antenna, usually calculated with Detector::GetAntPhase_1D()
        heff:  Effective height of the antenna, usually calculated with Report::GainToHeight
        Pol_vector:  polarization vector of electric field as it reaches the antenna.  Used to calculate pol_factor
        ant_type:  integer indicates Vpol (0) or HPol (1) antenna.  Use to calculate pol_factor
        pol_factor:  The projection of the E-field onto the antenna.  Essentially the cos(theta) term in a dot product between the E-field and effective height vector.
        vm_real:  Real component of the voltage in this frequency bin of the fourier space.
        vm_img:  Imaginary component of the voltage in this frequency bin of the fourier space.
        settings1:  The ARA settings object containing the setup of the simulation or detector.
        antenna_theta:  Incoming angle of the k-vector of the E-field when incident on the antenna.
        antenna_phi:  Incoming angle of the k-vector of the E-field when incident on the antenna.
        freq:  Frequency of the current frequency bin.
        useInTransmitterMode:  Boolean that dictates whether this antenna is in Transmitter mode (Tx) or Receiver mode (Rx)
        applyInverse:  Boolean that dictates if we're applying the antenna factors (projecting E-field onto h to get voltage) or inverting the factors (dividing h out of voltage to get E-field).
        
    Outputs:
        New values of the following via pass-by-reference:
        pol_factor, vm_real, vm_img
    
    */
    double phaseSign = 1.;
    double amplitudeSign = 1.;
    //If using in transmitter mode, the phase gets a minus sign since the signal is out-going from the antenna rather than incoming.
    if(useInTransmitterMode==true){ phaseSign*=-1.;};
    //If using this function to invert the antenna response, we apply a minus sign the phase in order to undo the phase applied by the antenna.  We also divide the amplitude by the factor rather than multiply.  To do this, we simply apply a -1 to the power of the factor applied to the amplitude so that it is divided out.
    if(applyInverse==true){ phaseSign*=-1.; amplitudeSign*=-1;};

    //Calculate the polarization factor, which is essentially the vector component of the Electric field that is projected onto the antenna.
    pol_factor = calculatePolFactor(Pol_vector, ant_type, antenna_theta, antenna_phi);
    if ( settings1->PHASE_SKIP_MODE != 1 ) {
        double phase_current;
        if ( vm_real != 0. ) {
            phase_current = atan( vm_img / vm_real );
            // phase in +-PI range
            if (vm_real<0.) {
                if (vm_img>0.) phase_current += PI;
                else if (vm_img<0.) phase_current -= PI;
            }
        }
        else {
            if ( vm_img>0. ) phase_current = PI;
            else if (vm_img<0.) phase_current = -PI;
            else phase_current = 0.;
        }
        //Calculate amplitude via the real and imaginary components.
        double v_amp  = sqrt(vm_real*vm_real + vm_img*vm_img);
        //Apply effective height and polarization factors to ampltitude.
        v_amp *= pow(heff*pol_factor, amplitudeSign);
        //If in transmitter mode, we must apply additional frequency and impedance terms to the amplitude.
        if (useInTransmitterMode==true){ 
            phase_current += PI/2;
            // The factors of two in this line are currently up for debate.  Will be cleaned up in future push - JCF 3/2/2024 
            v_amp *= pow(freq/CLIGHT*(Z0/Zr)/4/sqrt(2.), amplitudeSign);          
        }
        //Calculate the real and imaginary terms using the new ampltitude and phase.
        vm_real = v_amp * cos( phase_current + (phaseSign * AntPhase*RADDEG) );
        vm_img =  v_amp * sin( phase_current + (phaseSign * AntPhase*RADDEG) );
    }

    else { // only amplitude
        vm_real *= pow(heff * pol_factor, amplitudeSign);
        vm_img *= pow(heff * pol_factor, amplitudeSign); 
    }
}



void Report::ApplyAntFactors_Tdomain_FirstTwo (double heff, double heff_lastbin, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_bin0, double &vm_bin1, double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode, bool applyInverse) {
    /* Report::ApplyAntFactors_Tdomain_FirstTwo()
    Purpose: 
        Multiply (or divide) first and last bin of voltage in Fourier space by the antenna gain and phase.
        
    For a more verbose explanation, see Report::ApplyAntFactors_Tdomain.  This simply applies the operation to the first and last bin, which much be done separately as a consequence of FFT's.  Does not calculate phase.
    */

    double amplitudeSign = 1.;
    //If using this function to invert the antenna response, we apply a minus sign the phase in order to undo the phase applied by the antenna.  We also divide the amplitude by the factor rather than multiply.  To do this, we simply apply a -1 to the power of the factor applied to the amplitude so that it is divided out.
    if(applyInverse==true){amplitudeSign*=-1;};

    //double pol_factor;
    pol_factor = calculatePolFactor(Pol_vector, ant_type, antenna_theta, antenna_phi);

    //First step in the voltage calculation.  Both Tx and Rx mode use this.
    vm_bin0 *= pow(heff * pol_factor, amplitudeSign);
    vm_bin1 *= pow(heff_lastbin * pol_factor, amplitudeSign);
    
    if (useInTransmitterMode) {
        // The factors of two in these two lines are currently up for debate.  Will be cleaned up in future push - JCF 3/2/2024  
        vm_bin0 *= pow(freq/CLIGHT*(Z0/(Zr))/4/sqrt(2.), amplitudeSign);
        vm_bin1 *= pow(freq/CLIGHT*(Z0/(Zr))/4/sqrt(2.), amplitudeSign);
    }
    else {
        vm_bin0 *= pow(heff * pol_factor, amplitudeSign);
        vm_bin1 *= pow(heff_lastbin * pol_factor, amplitudeSign);
    }

}

void Report::InvertAntFactors_Tdomain (double AntPhase, double heff, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img, Settings *settings1, double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode) {  
    /* Report::InvertAntFactors_Tdomain()
    Purpose:  Inverts the antenna factors in a convenient way by simply calling ApplyAntFactors with the boolean applyInverse enabled.
    */
    ApplyAntFactors_Tdomain (AntPhase, heff, Pol_vector, ant_type, pol_factor, vm_real, vm_img, settings1, antenna_theta, antenna_phi, freq, useInTransmitterMode, true);

}



void Report::InvertAntFactors_Tdomain_FirstTwo (double heff, double heff_lastbin, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_bin0, double &vm_bin1, double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode) { 
    /* Report::InvertAntFactors_Tdomain_FirstTwo()
    Purpose:  Inverts the antenna factors in a convenient way by simply calling ApplyAntFactors with the boolean applyInverse enabled.
    */

    ApplyAntFactors_Tdomain_FirstTwo (heff, heff_lastbin, Pol_vector, ant_type, pol_factor, vm_bin0, vm_bin1, antenna_theta, antenna_phi, freq, useInTransmitterMode, true);

}



void Report::ApplyFilter(int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetFilterGain(bin_n) )/20.);   // from dB to unitless gain for voltage

}


void Report::ApplyPreamp(int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetPreampGain(bin_n) )/20.);   // from dB to unitless gain for voltage

}


void Report::ApplyFOAM(int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetFOAMGain(bin_n) )/20.);   // from dB to unitless gain for voltage

}



void Report::ApplyFilter_NFOUR (int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetFilterGain_NFOUR(bin_n) )/20.);   // from dB to unitless gain for voltage

}

void Report::ApplyFilter_OutZero (double freq, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetFilterGain_1D_OutZero(freq) )/20.);   // from dB to unitless gain for voltage

}


void Report::ApplyElect_Tdomain(double freq, Detector *detector, double &vm_real, double &vm_img, int gain_ch_no, Settings *settings1, bool applyInverse) { 
    /*  Report::ApplyElect_Tdomain()
    Purpose: 
        Multiply (or divide) voltage in Fourier space by the electronics gain and phase.
    
    Inputs:
        freq:  Frequency of the current frequency bin.
        detector: ARA detector object, which contains the config and antenna info to look up the gain and phase.
        vm_real:  Real component of the voltage in this frequency bin of the fourier space.
        vm_img:  Imaginary component of the voltage in this frequency bin of the fourier space.
        gain_ch_no:  Integer that represents the RF channel of the antenna, so that it can look up the gain and phase info.
        settings1:  The ARA settings object containing the setup of the simulation or detector.
        applyInverse:  Boolean that dictates if we're applying the electronics factors or inverting the factors.
        
    Outputs:
        New values of the following via pass-by-reference:
        vm_real, vm_img
    
    */
    double phaseSign = 1.;
    double amplitudeSign = 1.;
    //If using this function to invert the electronics response, we apply a minus sign the phase in order to undo the phase applied by the electronics.  We also divide the amplitude by the factor rather than multiply.  To do this, we simply apply a -1 to the power of the factor applied to the amplitude so that it is divided out.
    if(applyInverse==true){ phaseSign*=-1.; amplitudeSign*=-1;};

    if ( settings1->PHASE_SKIP_MODE == 0 ) {

        double phase_current;

        if ( vm_real != 0. ) {

            phase_current = atan( vm_img / vm_real );

            // phase in +-PI range
            if (vm_real<0.) {
                if (vm_img>0.) phase_current += PI;
                else if (vm_img<0.) phase_current -= PI;
            }
        }
        else {

            if ( vm_img>0. ) phase_current = PI;
            else if (vm_img<0.) phase_current = -PI;
            else phase_current = 0.;
        }

        // V amplitude
        double v_amp  = sqrt(vm_real*vm_real + vm_img*vm_img) * pow(detector->GetElectGain_1D_OutZero( freq, gain_ch_no ),amplitudeSign); // apply gain (unitless) to amplitude

        vm_real = v_amp * cos( phase_current - phaseSign*detector->GetElectPhase_1D(freq, gain_ch_no) );
        vm_img = v_amp * sin( phase_current - phaseSign*detector->GetElectPhase_1D(freq, gain_ch_no ) );
    }

    else {

        vm_real = vm_real * pow(detector->GetElectGain_1D_OutZero( freq, gain_ch_no),amplitudeSign); // only amplitude

        vm_img = vm_img * pow(detector->GetElectGain_1D_OutZero( freq, gain_ch_no),amplitudeSign); // only amplitude
    }
    
    //Apply power splitter/attenuator based on station.
    ApplySplitterFactor(vm_real, vm_img, detector, settings1, applyInverse);

}




void Report::ApplyElect_Tdomain_FirstTwo(double freq0, double freq1, Detector *detector, double &vm_bin0, double &vm_bin1, int gain_ch_no, Settings *settings1, bool applyInverse) {  // read elect chain gain (unitless), phase (rad) and apply to V/m
    /*  Report::ApplyElect_Tdomain_FirstTwo()
    Purpose: 
        Multiply (or divide) first two bins of voltage in Fourier space by the electronics gain and phase.
    
    For a more verbose explanation, see Report::ApplyElect_Tdomain.  This simply applies the operation to the first and last bin, which much be done separately as a consequence of FFT's.  Does not calculate phase.

    
    */
    double amplitudeSign = 1.;
    //If using this function to invert the electronics response, we apply a minus sign the phase in order to undo the phase applied by the electronics.  We also divide the amplitude by the factor rather than multiply.  To do this, we simply apply a -1 to the power of the factor applied to the amplitude so that it is divided out.
    if(applyInverse==true){amplitudeSign*=-1;};

    vm_bin0 = vm_bin0 * pow(detector->GetElectGain_1D_OutZero( freq0 , gain_ch_no),amplitudeSign);
    vm_bin1 = vm_bin1 * pow(detector->GetElectGain_1D_OutZero( freq1 , gain_ch_no),amplitudeSign);

    //Apply power splitter/attenuator based on station.
    ApplySplitterFactor(vm_bin0, vm_bin1, detector, settings1, applyInverse);
}

void Report::ApplySplitterFactor(double &vm_real, double &vm_img, Detector* detector, Settings *settings1, bool applyInverse) {
    // Apply splitter/attenuation factor in the digitizer path based on station.
    
    //If we wish to invert the splitter factor, we must divide the factor out rather than multiply.  To accomplish this, we simply raise the factor to the -1 power in the multiplication step.
    double amplitudeSign=1;
    if (applyInverse==true){amplitudeSign*=-1;};    
  
    const double factor = detector->GetSplitterFactor(settings1);

    vm_real *= pow(factor, amplitudeSign);
    vm_img *= pow(factor, amplitudeSign);

    return;
}

void Report::InvertElect_Tdomain(double freq, Detector *detector, double &vm_real, double &vm_img, int gain_ch_no, Settings *settings1) {  // read elect chain gain (unitless), phase (rad) and apply to V/m
    /* Report::InvertElect_Tdomain()
    Purpose:  Inverts the electronics factors in a convenient way by simply calling ApplyElecFactors with the boolean applyInverse enabled.
    */
    ApplyElect_Tdomain(freq, detector, vm_real, vm_img, gain_ch_no, settings1, true);
  
}

void Report::InvertElect_Tdomain_FirstTwo(double freq0, double freq1, Detector *detector, double &vm_bin0, double &vm_bin1, int gain_ch_no, Settings *settings1) {  // read elect chain gain (unitless), phase (rad) and apply to V/m
    /* Report::InvertElect_Tdomain_FirstTwo()
    Purpose:  Inverts the electronics factors in a convenient way by simply calling ApplyElecFactors with the boolean applyInverse enabled.
    */
    ApplyElect_Tdomain_FirstTwo(freq0, freq1, detector, vm_bin0, vm_bin1, gain_ch_no, settings1, true);
    
}





void Report::ApplyPreamp_NFOUR (int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetPreampGain_NFOUR(bin_n) )/20.);   // from dB to unitless gain for voltage

}


void Report::ApplyPreamp_OutZero (double freq, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetPreampGain_1D_OutZero(freq) )/20.);   // from dB to unitless gain for voltage

}


void Report::ApplyFOAM_NFOUR (int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetFOAMGain_NFOUR(bin_n) )/20.);   // from dB to unitless gain for voltage

}


void Report::ApplyFOAM_OutZero (double freq, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetFOAMGain_1D_OutZero(freq) )/20.);   // from dB to unitless gain for voltage

}



void Report::ApplyRFCM(int ch, int bin_n, Detector *detector, double &vmmhz, double RFCM_OFFSET) {  // read RFCM gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetRFCMGain(ch,bin_n) + RFCM_OFFSET )/20.);   // from dB to unitless gain for voltage

}


void Report::ApplyFilter_databin(int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetFilterGain_databin(bin_n) )/20.);   // from dB to unitless gain for voltage

}

void Report::ApplyPreamp_databin(int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetPreampGain_databin(bin_n) )/20.);   // from dB to unitless gain for voltage

}

void Report::ApplyFOAM_databin(int bin_n, Detector *detector, double &vmmhz) {  // read filter gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetFOAMGain_databin(bin_n) )/20.);   // from dB to unitless gain for voltage

}


void Report::ApplyRFCM_databin(int ch, int bin_n, Detector *detector, double &vmmhz, double RFCM_OFFSET) {  // read RFCM gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetRFCMGain_databin(ch,bin_n) + RFCM_OFFSET )/20.);   // from dB to unitless gain for voltage

}


void Report::ApplyNoiseFig_databin(int ch, int bin_n, Detector *detector, double &vmmhz, Settings *settings1) {  // read noise figure and apply unitless gain to vmmhz
  //cout<<"Entered ApplyNoiseFig_databin    ";
  //cout<<"vmmhz: "<<vmmhz;
  
  double tempNoise = vmmhz*vmmhz;
  //  cout<<"    1: "<<detector->GetNoiseFig_databin(ch,bin_n);
  //  cout<<"    2: "<<detector->GetTransm_databin(ch, bin_n);
  //  cout<<"    3: "<<settings1->NOISE_TEMP;
  //FILE *fp = fopen("noiseTemp_NOISE_CHANNEL_MODE_0.txt","a+");
  if(detector->GetNoiseFig_databin(ch,bin_n)>1.0){vmmhz = TMath::Sqrt( tempNoise*(detector->GetTransm_databin(ch, bin_n)) + tempNoise/(settings1->NOISE_TEMP)*220.0*(detector->GetNoiseFig_databin(ch,bin_n) - 1.0) ) ;  
    //cout<<"   vmmhz   "<<vmmhz;
    //cout<<"   bin: "<<bin_n<<" freq: "<<detector->GetFreq(bin_n) << " noise temp: " << 246.0*((detector->GetTransm_databin(ch, bin_n)) + 1.0/(settings1->NOISE_TEMP)*220.0*(detector->GetNoiseFig_databin(ch,bin_n) - 1.0) ) << endl;  
    //fprintf(fp, "%f   ",246.0*((detector->GetTransm_databin(ch, bin_n)) + 1.0/(settings1->NOISE_TEMP)*220.0*(detector->GetNoiseFig_databin(ch,bin_n) - 1.0) ));
}
  else{
    vmmhz = vmmhz;
    //fprintf(fp, "%f   ",246.0);
  }
  
  //fclose(fp);
  //    if(ch==0 && bin_n==2050) cout << "Noise ch: "<< ch <<"   "<< detector->GetNoiseFig_databin(ch,bin_n) << "   " << detector->GetTransm_databin(ch, bin_n)
//		<< "   " << TMath::Sqrt( tempNoise ) << "   " << vmmhz << endl;
  //  if(ch==8 && bin_n==2050) cout << "Noise ch: "<< ch <<"   "<< detector->GetNoiseFig_databin(ch,bin_n) << "   " << detector->GetTransm_databin(ch, bin_n)
  //		<< "   " << TMath::Sqrt( tempNoise ) << "   " << vmmhz << endl;
  //    if(ch==24 && bin_n==2000) cout << "Noise ch: "<< ch << "   " << detector->GetNoiseFig_databin(ch,bin_n)*220.0/246.0 << "   " << TMath::Sqrt(detector->GetTransm_databin(ch, bin_n) ) << endl;
  
}




void Report::GetAngleAnt(Vector &rec_vector, Position &antenna, double &ant_theta, double &ant_phi) {   //ant_theta and ant_phi is in degree 

    /*
     * 2020-12-07 BAC
     * We take the receive vector, which currently specificies where the signal is *going*
     * And flip it to specify where the signal is *coming from*
     * And then get the theta and phi, and return them in *degrees*.
     *
     * This is necessary because AraSim requires antenna_theta and antenna_phi
     * to indicate "where the signal came from" in order to correctly
     * calculate the polarization factors and the gain.
    */

    Vector flip_receive_vector = -1. * rec_vector;
    ant_theta = flip_receive_vector.Theta() * TMath::RadToDeg(); // return in degrees
    ant_phi = flip_receive_vector.Phi() * TMath::RadToDeg();

}

void Report::GetAngleLaunch(Vector &launch_vector, double &launch_theta, double &launch_phi) {
    
    /*
    2024-07-01 JCF
    Takes the launch vector of the signal and calculates the theta and phi of the launch vector in the 
    station-centric coordinates and returns them in degrees; where theta=0 is along the upward vertical axis,
    theta=90 is along the horizon, and theta=180 is along the downward vertical axis.
    
    This is necessary because for source reconstruction, we need the launch vector of the RF signal in 
    conjunction with the polarization to find the neutrino trajectory.
    */
    
    launch_theta = launch_vector.Theta() * TMath::RadToDeg();
    launch_phi = launch_vector.Phi() * TMath::RadToDeg();

}


// generate DATA_BIN_SIZE sized noise waveform array in time domain
void Report::GetNoiseWaveforms(Settings *settings1, Detector *detector, double v_noise, double *vnoise) {

    if (settings1->NOISE == 0) {    // NOISE == 0 : flat thermal noise with Johnson-Nyquist noise
        //V_noise_fft_bin = sqrt( (double)(NFOUR/2) * 50. * KBOLTZ * T_noise / (2. * TIMESTEP) ); 

        Vfft_noise_after.clear();  // remove previous Vfft_noise values
        Vfft_noise_before.clear();  // remove previous Vfft_noise values
        //V_noise_timedomain.clear(); // remove previous V_noise_timedomain values


        double V_tmp; // copy original flat H_n [V] value
        double current_amplitude, current_phase;

        GetNoisePhase(settings1); // get random phase for noise

        //MakeArraysforFFT_noise(settings1, detector, vhz_noise_tmp , vnoise);
        // returned array vnoise currently have real value = imag value as no phase term applied

        //for (int k=0; k<settings1->NFOUR/4; k++) {
        for (int k=0; k<settings1->DATA_BIN_SIZE/2; k++) {

            V_tmp = v_noise / sqrt(2.) ; // copy original flat H_n [V] value, and apply 1/sqrt2 for SURF/TURF divide same as signal

            if (settings1->USE_TESTBED_RFCM_ON == 0) {
                ApplyFilter_databin(k, detector, V_tmp);
                ApplyPreamp_databin(k, detector, V_tmp);
                ApplyFOAM_databin(k, detector, V_tmp);
		if (settings1->APPLY_NOISE_FIGURE==1){
		  ApplyNoiseFig_databin(0,k,detector, V_tmp, settings1);
		}
            }
            else if (settings1->USE_TESTBED_RFCM_ON == 1) {
                // apply RFCM gain
                // as this mode don't have different chs' noise waveform separately, just use ch0 RFCM gain...
                cerr<<"Trying not allowed mode : NOISE_CHANNEL_MODE=0 and USE_TESTBED_RFCM_ON=1 same time!"<<endl;
                break;
                //ApplyRFCM_databin(0, detector, V_tmp);
            }


            Vfft_noise_before.push_back( V_tmp );


            current_phase = noise_phase[k];

            //Tools::get_random_rician( 0., 0., sqrt(2./ M_PI) * V_tmp, current_amplitude, current_phase);    // use real value array value
            Tools::get_random_rician( 0., 0., sqrt(2./M_PI)/1.177 * V_tmp, current_amplitude, current_phase);    // use real value array value, extra 1/1.177 to make total power same with "before random_rician".

            // vnoise is currently noise spectrum (before fft, unit : V)
           //vnoise[2 * k] = sqrt(current_amplitude) * cos(noise_phase[k]);
           //vnoise[2 * k + 1] = sqrt(current_amplitude) * sin(noise_phase[k]);
           vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
           vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);



           //vnoise[2 * k] = (V_tmp) * cos(noise_phase[k]);
           //vnoise[2 * k + 1] = (V_tmp) * sin(noise_phase[k]);
           

            Vfft_noise_after.push_back( vnoise[2*k] );
            Vfft_noise_after.push_back( vnoise[2*k+1] );

            // inverse FFT normalization factor!
            vnoise[2 * k] *= 2./((double)settings1->DATA_BIN_SIZE);
            vnoise[2 * k + 1] *= 2./((double)settings1->DATA_BIN_SIZE);


        }


        // now vnoise is time domain waveform
        Tools::realft( vnoise, -1, settings1->DATA_BIN_SIZE);

        
        // save timedomain noise to Report class
        /*
        for (int k=0; k<settings1->DATA_BIN_SIZE; k++) {
            V_noise_timedomain.push_back( vnoise[k] );
        }
        */

    }
    else {  // currently there are no more options!
        cout<<"no noise option for NOISE = "<<settings1->NOISE<<endl;
    }


}



// generate DATA_BIN_SIZE sized noise waveform array in time domain
void Report::GetNoiseWaveforms_ch(Settings * settings1, Detector * detector, double v_noise, double * vnoise, int ch) {
    //  cout << "channel: " << ch << endl;
    if (settings1 -> NOISE == 0) { // NOISE == 0 : flat thermal noise with Johnson-Nyquist noise

        Vfft_noise_after.clear(); // remove previous Vfft_noise values
        Vfft_noise_before.clear(); // remove previous Vfft_noise values
        //V_noise_timedomain.clear(); // remove previous V_noise_timedomain values

        double V_tmp; // copy original flat H_n [V] value
        double current_amplitude, current_phase;

        GetNoisePhase(settings1); // get random phase for noise

        for (int k = 0; k < settings1 -> DATA_BIN_SIZE / 2; k++) {

            V_tmp = v_noise / sqrt(2.); // copy original flat H_n [V] value, and apply 1/sqrt2 for SURF/TURF divide same as signal

            if (settings1 -> USE_TESTBED_RFCM_ON == 0) {
                ApplyFilter_databin(k, detector, V_tmp);
                ApplyPreamp_databin(k, detector, V_tmp);
                ApplyFOAM_databin(k, detector, V_tmp);
                if (settings1 -> APPLY_NOISE_FIGURE == 1) {
                    ApplyNoiseFig_databin(ch % 16, k, detector, V_tmp, settings1);
                }
            } else if (settings1 -> USE_TESTBED_RFCM_ON == 1) {
                // apply RFCM gain
                ApplyRFCM_databin(ch, k, detector, V_tmp, settings1 -> RFCM_OFFSET);
            }

            Vfft_noise_before.push_back(V_tmp);

            current_phase = noise_phase[k];

            //Tools::get_random_rician( 0., 0., sqrt(2./ M_PI) * V_tmp, current_amplitude, current_phase);    // use real value array value
            Tools::get_random_rician(0., 0., sqrt(2. / M_PI) / 1.177 * V_tmp, current_amplitude, current_phase); // use real value array value, extra 1/1.177 to make total power same with "before random_rician".

            // vnoise is currently noise spectrum (before fft, unit : V)
            //vnoise[2 * k] = sqrt(current_amplitude) * cos(noise_phase[k]);
            //vnoise[2 * k + 1] = sqrt(current_amplitude) * sin(noise_phase[k]);
            vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
            vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);

            Vfft_noise_after.push_back(vnoise[2 * k]);
            Vfft_noise_after.push_back(vnoise[2 * k + 1]);

            // inverse FFT normalization factor!
            vnoise[2 * k] *= 2. / ((double) settings1 -> DATA_BIN_SIZE);
            vnoise[2 * k + 1] *= 2. / ((double) settings1 -> DATA_BIN_SIZE);

        }

        //	cout << "After applying noise figure" << endl;

        // now vnoise is time domain waveform
        Tools::realft(vnoise, -1, settings1 -> DATA_BIN_SIZE);

        // save timedomain noise to Report class
        /*
        for (int k=0; k<settings1->DATA_BIN_SIZE; k++) {
            V_noise_timedomain.push_back( vnoise[k] );
        }
        */
        //	cout << "After fft" << endl;

    } else if (settings1 -> NOISE == 1) { // NOISE == 1 : use Rayleigh dist. fits

        Vfft_noise_after.clear(); // remove previous Vfft_noise values
        Vfft_noise_before.clear(); // remove previous Vfft_noise values

        double V_tmp; // copy original flat H_n [V] value
        double current_amplitude, current_phase;

        GetNoisePhase(settings1); // get random phase for noise

        // BAC, June 17 2022
        // Do something different for TestBed and Deep stations.
        // I'm not thrilled about this solution, but it's the way forward that requires
        // the least re-factorization of the Testbed code

        if(settings1->DETECTOR==3 && settings1->DETECTOR_STATION==0){

            for (int k = 0; k < settings1 -> DATA_BIN_SIZE / 2; k++) {

                if (ch < detector -> RayleighFit_ch) {

                    Vfft_noise_before.push_back(detector -> GetRayleighFit_databin(ch, k));

                    current_phase = noise_phase[k];

                    V_tmp = detector -> GetRayleighFit_databin(ch, k) * sqrt((double) settings1 -> DATA_BIN_SIZE / (double)(settings1 -> NFOUR / 2.));

                    //Tools::get_random_rician( 0., 0., sqrt(2./ M_PI) * V_tmp, current_amplitude, current_phase);    // use real value array value
                    //Tools::get_random_rician( 0., 0., sqrt(2./M_PI)/1.177 * V_tmp, current_amplitude, current_phase);    // use real value array value, extra 1/1.177 to make total power same with "before random_rician".
                    Tools::get_random_rician(0., 0., V_tmp, current_amplitude, current_phase); // use real value array value, extra 1/1.177 to make total power same with "before random_rician".

                } else {

                    V_tmp = v_noise / sqrt(2.); // copy original flat H_n [V] value, and apply 1/sqrt2 for SURF/TURF divide same as signal

                    if (settings1 -> USE_TESTBED_RFCM_ON == 0) {
                        ApplyFilter_databin(k, detector, V_tmp);
                        ApplyPreamp_databin(k, detector, V_tmp);
                        ApplyFOAM_databin(k, detector, V_tmp);

                    } else if (settings1 -> USE_TESTBED_RFCM_ON == 1) {
                        // apply RFCM gain
                        ApplyRFCM_databin(ch, k, detector, V_tmp, settings1 -> RFCM_OFFSET);
                    }

                    Vfft_noise_before.push_back(V_tmp);

                    current_phase = noise_phase[k];

                    //Tools::get_random_rician( 0., 0., sqrt(2./ M_PI) * V_tmp, current_amplitude, current_phase);    // use real value array value
                    Tools::get_random_rician(0., 0., sqrt(2. / M_PI) / 1.177 * V_tmp, current_amplitude, current_phase); // use real value array value, extra 1/1.177 to make total power same with "before random_rician".

                }

                // vnoise is currently noise spectrum (before fft, unit : V)
                //vnoise[2 * k] = sqrt(current_amplitude) * cos(noise_phase[k]);
                //vnoise[2 * k + 1] = sqrt(current_amplitude) * sin(noise_phase[k]);
                vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
                vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);

                Vfft_noise_after.push_back(vnoise[2 * k]);
                Vfft_noise_after.push_back(vnoise[2 * k + 1]);

                // inverse FFT normalization factor!
                vnoise[2 * k] *= 2. / ((double) settings1 -> DATA_BIN_SIZE);
                vnoise[2 * k + 1] *= 2. / ((double) settings1 -> DATA_BIN_SIZE);

            }
        }
        else if(settings1->DETECTOR==4 && settings1->DETECTOR_STATION>0){

            // check to make sure we have this channel available
            if(ch < 0 || ch >= detector->RayleighFit_ch){
                char errorMessage[400];
                sprintf(errorMessage,
                    "L%d: You have asked for Rayleigh Fits for ch %d, which is not supported",
                    __LINE__, ch
                );
                throw std::runtime_error(errorMessage);
            }

            // get the fits for this specific station
            // this function will throw exceptions if the station doesn't exist
            // so we can call this safely
            auto fits_for_this_station = detector->GetRayleighFitVector_databin(settings1->DETECTOR_STATION, settings1);

            // to normalize the bin content, we also need to keep delta f
            double this_delta_f = 1./(settings1->DATA_BIN_SIZE * settings1->TIMESTEP);

            // loop over frequency bins
            for(int k=0; k<settings1->DATA_BIN_SIZE/2; k++){

                Vfft_noise_before.push_back(fits_for_this_station[ch][k]);
                current_phase = noise_phase[k];
                V_tmp = fits_for_this_station[ch][k];
                
                // the right normalization factor is N * sqrt(deltaF)
                V_tmp *= double(settings1->DATA_BIN_SIZE);
                V_tmp *= sqrt(this_delta_f);

                Tools::get_random_rician(0., 0., V_tmp, current_amplitude, current_phase); // draw a random number from the distribution with this rayleigh fit parameter

                // set the real and imaginary components of the FFT
                vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
                vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);
                
                // stash those values
                Vfft_noise_after.push_back(vnoise[2 * k]);
                Vfft_noise_after.push_back(vnoise[2 * k + 1]);

                // and apply the inverse FFT normalization factor
                vnoise[2 * k] *= 2. / ((double) settings1 -> DATA_BIN_SIZE);
                vnoise[2 * k + 1] *= 2. / ((double) settings1 -> DATA_BIN_SIZE);

            }

            // real FT back to get vnoise in time domain waveform; 
            Tools::realft(vnoise, -1, settings1 -> DATA_BIN_SIZE);

        }
        else if(settings1->DETECTOR==5 && settings1->DETECTOR_STATION>0){

            // check to make sure we have this channel available
            if(ch < 0 || ch >= detector->RayleighFit_ch){
                char errorMessage[400];
                sprintf(errorMessage,
                    "L%d: You have asked for Rayleigh Fits for ch %d, which is not supported",
                    __LINE__, ch
                );
                throw std::runtime_error(errorMessage);
            }

            // get the fits for this specific station
            // this function will throw exceptions if the station doesn't exist
            // so we can call this safely
            auto fits_for_this_station = detector->GetRayleighFitVector_databin(settings1->DETECTOR_STATION, settings1);

            // to normalize the bin content, we also need to keep delta f
            double this_delta_f = 1./(settings1->DATA_BIN_SIZE * settings1->TIMESTEP);

            // loop over frequency bins
            for(int k=0; k<settings1->DATA_BIN_SIZE/2; k++){

                Vfft_noise_before.push_back(fits_for_this_station[ch][k]);
                current_phase = noise_phase[k];
                V_tmp = fits_for_this_station[ch][k];
                
                // the right normalization factor is N * sqrt(deltaF)
                V_tmp *= double(settings1->DATA_BIN_SIZE);
                V_tmp *= sqrt(this_delta_f);

                Tools::get_random_rician(0., 0., V_tmp, current_amplitude, current_phase); // draw a random number from the distribution with this rayleigh fit parameter

                // set the real and imaginary components of the FFT
                vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
                vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);
                
                // stash those values
                Vfft_noise_after.push_back(vnoise[2 * k]);
                Vfft_noise_after.push_back(vnoise[2 * k + 1]);

                // and apply the inverse FFT normalization factor
                vnoise[2 * k] *= 2. / ((double) settings1 -> DATA_BIN_SIZE);
                vnoise[2 * k + 1] *= 2. / ((double) settings1 -> DATA_BIN_SIZE);

            }

            // real FT back to get vnoise in time domain waveform; 
            Tools::realft(vnoise, -1, settings1 -> DATA_BIN_SIZE);

        }

    } else { // currently there are no more options!
        cout << "no noise option for NOISE = " << settings1 -> NOISE << endl;
    }

}




void Report::GetNoisePhase(Settings *settings1) {

    noise_phase.clear();    // remove all previous noise_phase vector values

    //for (int i=0; i<settings1->NFOUR/4; i++) {
    for (int i=0; i<settings1->DATA_BIN_SIZE/2; i++) {  // noise with DATA_BIN_SIZE bin array
        noise_phase.push_back(2*PI*gRandom->Rndm());  // get random phase for flat thermal noise
    }
}


void Report::Prepare_Antenna_Noise(
    int debugmode, int channel_index, 
    int station_number, int string_number, int antenna_number,
    Settings *settings1, Trigger *trigger, Detector *detector
){

    // Determine which noise waveforms from the trigger class we'd like to use
    //   for each antenna then export said waveforms into trigger->Full_window
    //   and trigger->Full_window_V 

    for (int l = 0; l < N_noise; l++)
    {

        // select noise waveform/ID from trigger class
        // if we are sharing same noise waveform for all chs, make sure diff chs use diff noise waveforms
        if (settings1->NOISE_CHANNEL_MODE == 0)
        {

            // get random noise_ID and check if there are same noise_ID in different ch.
            // if there's same noise_ID, get noise_ID again until no noise_ID are same between chs
            noise_pass_nogo = 1;
            while (noise_pass_nogo)
            {
                noise_ID[l] = (int)(settings1->NOISE_EVENTS *gRandom->Rndm());
                noise_pass_nogo = 0;
                for (int j_sub = 0; j_sub < string_number + 1; j_sub++)
                {
                    for (int k_sub = 0; k_sub < detector->stations[station_number].strings[string_number].antennas.size(); k_sub++)
                    {
                        if (j_sub == string_number)
                        {
                            // if we are checking current string
                            if (k_sub < antenna_number)
                            {
                                // check antennas before current antenna
                                if (noise_ID[l] == stations[station_number].strings[j_sub].antennas[k_sub].noise_ID[0])
                                {
                                    // check only first one for now;;;
                                    noise_pass_nogo = 1;
                                }
                            }
                        }
                        else // if we are checking previous string, check upto entire antennas
                        {
                            // Avoids segfault from string 2 having 1 antenna for PA DETECTOR_STATION 3
                            if (
                                settings1->DETECTOR==5 && // Phased Array detector mode
                                settings1->DETECTOR_STATION==3 && // second PA detector configuration (PA + 7 ARA VPols)
                                j_sub==2 && // String 2
                                k_sub==1 // 2nd antenna that doesn't exist but this block of code expects it to
                            ) {
                                continue;
                            } 

                            // check only first one for now;;;
                            if (noise_ID[l] == stations[station_number].strings[j_sub].antennas[k_sub].noise_ID[0])
                            {
                                noise_pass_nogo = 1;
                            }
                        }
                    } // end loop over k_sub (antennas)
                } // end loop over j_sub (strings)
            }   // end while noise_pass_nogo

        }   // end if NOISE_CHANNEL_MODE = 0

        // if we are using diff noise waveform for diff chs, just pick any noise waveform
        else if (settings1->NOISE_CHANNEL_MODE == 1)
        {
            noise_ID[l] = (int)(settings1->NOISE_EVENTS *gRandom->Rndm());
        }   // end if NOISE_CHANNEL_MODE = 1

        // if we are sharing same noise waveform for first 8 chs and share same noisewaveforms for others,  just pick any noise waveform
        else if (settings1->NOISE_CHANNEL_MODE == 2)
        {
            noise_ID[l] = (int)(settings1->NOISE_EVENTS *gRandom->Rndm());
        }   // end if NOISE_CHANNEL_MODE = 2

        // save noise ID
        stations[station_number].strings[string_number].antennas[antenna_number].noise_ID.push_back(noise_ID[l]);

        // Export noise waveforms to trigger->Full_window and Full_window_V, only if it's not in debugmode
        if (debugmode == 0)
        {

            // Export noise generated from a single temperature
            if (settings1->NOISE_CHANNEL_MODE == 0)
            {
                if (l == N_noise - 1)
                {
                    // when it's final noise waveform
                    for (int bin = 0; bin < settings1->DATA_BIN_SIZE; bin++)
                    {
                        // test for full window
                        trigger->Full_window[channel_index][bin] = (trigger->v_noise_timedomain_diode[noise_ID[l]][bin]);
                        trigger->Full_window_V[channel_index][bin] = (trigger->v_noise_timedomain[noise_ID[l]][bin]);
                    }
                    // cout<<"last noise filled in Full_window!"<<endl;
                }
                else
                {
                    // when it's not final noise waveform
                    // cout<<"full noise will fill in Full_window!"<<endl;
                    for (int bin = 0; bin < settings1->DATA_BIN_SIZE; bin++)
                    {
                        trigger->Full_window[channel_index][bin] = (trigger->v_noise_timedomain_diode[noise_ID[l]][bin]);
                        trigger->Full_window_V[channel_index][bin] = (trigger->v_noise_timedomain[noise_ID[l]][bin]);
                    }
                }
            } // end NOISE_CHANNEL_MODE 0

            // Export noise generated from different temperatures for each antenna
            else if (settings1->NOISE_CHANNEL_MODE == 1)
            {
                if (l == N_noise - 1)
                {
                    // when it's final noise waveform
                    for (int bin = 0; bin < settings1->DATA_BIN_SIZE; bin++)
                    {
                        // test for full window
                        trigger->Full_window[channel_index][bin] = (
                            trigger->v_noise_timedomain_diode_ch[
                                GetChNumFromArbChID(detector, channel_index, station_number, settings1) - 1
                            ][ noise_ID[l] ][ bin ] );
                        trigger->Full_window_V[channel_index][bin] = (
                            trigger->v_noise_timedomain_ch[
                                GetChNumFromArbChID(detector, channel_index, station_number, settings1) - 1
                            ][ noise_ID[l] ][ bin ] );
                    }
                    //cout<<"last noise filled in Full_window!"<<endl;
                }
                else
                {
                    // when it's not final noise waveform
                    //cout<<"full noise will fill in Full_window!"<<endl;
                    for (int bin = 0; bin < settings1->DATA_BIN_SIZE; bin++)
                    {
                        trigger->Full_window[channel_index][bin] = (
                            trigger->v_noise_timedomain_diode_ch[
                                GetChNumFromArbChID(detector, channel_index, station_number, settings1) - 1
                            ][ noise_ID[l] ][ bin ] );
                        trigger->Full_window_V[channel_index][bin] = (
                            trigger->v_noise_timedomain_ch[
                                GetChNumFromArbChID(detector, channel_index, station_number, settings1) - 1
                            ][ noise_ID[l] ][ bin ] );
                    }
                }
                
            } // end NOISE_CHANNEL_MODE 1

            // Export noise generated where first 8 channels have their own temperatures and the remaining 8 share one temperature
            else if (settings1->NOISE_CHANNEL_MODE == 2)
            {

                if ((GetChNumFromArbChID(detector, channel_index, station_number, settings1) - 1) < 8)
                {

                    if (l == N_noise - 1)
                    {
                        // when it's final noise waveform
                        for (int bin = 0; bin < settings1->DATA_BIN_SIZE; bin++)
                        {
                            // test for full window
                            trigger->Full_window[channel_index][bin] = (
                                trigger->v_noise_timedomain_diode_ch[
                                    GetChNumFromArbChID(detector, channel_index, station_number, settings1) - 1
                                ][ noise_ID[l] ][ bin ]);
                            trigger->Full_window_V[channel_index][bin] = (
                                trigger->v_noise_timedomain_ch[
                                    GetChNumFromArbChID(detector, channel_index, station_number, settings1) - 1
                                ][ noise_ID[l] ][ bin ] );
                        }
                        // cout<<"last noise filled in Full_window!"<<endl;
                    }
                    else
                    {
                        // when it's not final noise waveform
                        // cout<<"full noise will fill in Full_window!"<<endl;
                        for (int bin = 0; bin < settings1->DATA_BIN_SIZE; bin++)
                        {
                            trigger->Full_window[channel_index][bin] = (
                                trigger->v_noise_timedomain_diode_ch[
                                    GetChNumFromArbChID(detector, channel_index, station_number, settings1) - 1
                                ][ noise_ID[l] ][ bin ] );
                            trigger->Full_window_V[channel_index][bin] = (
                                trigger->v_noise_timedomain_ch[
                                    GetChNumFromArbChID(detector, channel_index, station_number, settings1) - 1
                                ][ noise_ID[l] ][ bin ] );
                        }
                    }

                } // end NOISE_CHANNEL_MODE 2

                // Export noise generated from the same temperature for all antennas
                else
                {

                    if (l == N_noise - 1)
                    {
                        // when it's final noise waveform
                        for (int bin = 0; bin < settings1->DATA_BIN_SIZE; bin++)
                        {
                            // test for full window
                            trigger->Full_window[channel_index][bin] = (trigger->v_noise_timedomain_diode_ch[8][noise_ID[l]][bin]);
                            trigger->Full_window_V[channel_index][bin] = (trigger->v_noise_timedomain_ch[8][noise_ID[l]][bin]);
                        }
                        // cout<<"last noise filled in Full_window!"<<endl;
                    }
                    else
                    {
                        // when it's not final noise waveform
                        // cout<<"full noise will fill in Full_window!"<<endl;
                        for (int bin = 0; bin < settings1->DATA_BIN_SIZE; bin++)
                        {
                            trigger->Full_window[channel_index][bin] = (trigger->v_noise_timedomain_diode_ch[8][noise_ID[l]][bin]);
                            trigger->Full_window_V[channel_index][bin] = (trigger->v_noise_timedomain_ch[8][noise_ID[l]][bin]);
                        }
                    }

                } 

            } // end consideration for all other NOISE_CHANNEL_MODE options

        }   // if we are not in debug mode

    } // end loop for l over N_noise

}


    

void Report::MakeArraysforFFT(Settings *settings1, Detector *detector, int StationIndex, vector <double> &vsignal_array, double *vsignal_forfft) {



    // from icemc, anita class MakeArraysforFFT
    int NFOUR = settings1->NFOUR;
    double TIMESTEP = settings1->TIMESTEP;

//    int NFOUR = detector->stations[StationIndex].NFOUR;
//    int TIMESTEP = detector->stations[StationIndex].TIMESTEP;
    Tools::Zero(vsignal_forfft,NFOUR/2);
    


    if (settings1->AVZ_NORM_FACTOR_MODE == 0) { // use previous normalization factors


        double previous_value_e_even=0.;
        double previous_value_e_odd=0.;
        int count_nonzero=0;
        int iprevious=0;
        int ifirstnonzero=-1;
        int ilastnonzero=2000;
        //for (int i=0;i<NFREQ;i++) {
        for (int i=0;i<detector->GetFreqBin();i++) {
            
            // freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
            // but there are only NFOUR/4 different values
            // it's the index among the NFOUR/4 that we're interested in
            int ifour=Tools::Getifreq(detector->GetFreq(i),detector->freq_forfft[0],detector->freq_forfft[NFOUR/2-1],NFOUR/4);
            
            if (ifour!=-1 && 2*ifour+1<NFOUR/2) {
                count_nonzero++;
                if (ifirstnonzero==-1)
                    ifirstnonzero=ifour;
                
                vsignal_forfft[2*ifour]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP); // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft
                
                //      cout << "ifour, vsignal is " << ifour << " " << vsignal_e_forfft[2*ifour] << "\n";
                
                vsignal_forfft[2*ifour+1]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP); // phase is 90 deg.
                // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting
                
                // how about we interpolate instead of doing a box average
                
                for (int j=iprevious+1;j<ifour;j++) {
                    vsignal_forfft[2*j]=previous_value_e_even+(vsignal_forfft[2*ifour]-previous_value_e_even)*(double)(j-iprevious)/(double)(ifour-iprevious);
                    //	cout << "j, vsignal is " << j << " " << vsignal_e_forfft[2*j] << "\n";
                    
                    vsignal_forfft[2*j+1]=previous_value_e_odd+(vsignal_forfft[2*ifour+1]-previous_value_e_odd)*(double)(j-iprevious)/(double)(ifour-iprevious);
                }
                
                ilastnonzero=ifour;
                iprevious=ifour;
                previous_value_e_even=vsignal_forfft[2*ifour];
                previous_value_e_odd=vsignal_forfft[2*ifour+1];
            }
            
        } // end loop over nfreq
        


        // don't applying any extra factor for the change in array (change of bin size)
        // as change in the bin size doesn't matter for the total energy
        // total energy is just same as integral over frequency range and this frequency range will not change
        //
          //cout << "ifirstnonzero, ilastnonzero are " << ifirstnonzero << " " << ilastnonzero << "\n";
          //cout << "non zero count is " << count_nonzero << "\n";
          //cout << "ratio is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
        for (int j=0;j<NFOUR/4;j++) {
            vsignal_forfft[2*j]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
            vsignal_forfft[2*j+1]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
        }



        
        //  Tools::InterpolateComplex(vsignal_e_forfft,NFOUR/4);
        //Tools::InterpolateComplex(vsignal_h_forfft,NFOUR/4);
        for (int ifour=0;ifour<NFOUR/4;ifour++) {



            vsignal_forfft[2*ifour]*=cos(settings1->PHASE*PI/180.);
            vsignal_forfft[2*ifour+1]*=sin(settings1->PHASE*PI/180.);
            
            //--------------------------------------------------
            // if (!PULSER) {
            //     
            //     vsignal_forfft[2*ifour]*=cos(phase*PI/180.);
            //     vsignal_forfft[2*ifour+1]*=sin(phase*PI/180.);
            //     
            //     
            // }
            // else {
            //     vsignal_forfft[2*ifour]*=cos(v_phases[ifour]*PI/180.);
            //     vsignal_forfft[2*ifour+1]*=sin(v_phases[ifour]*PI/180.);
            //     
            // }	  	  	  
            //-------------------------------------------------- 
            
            
        }
    }
    else if (settings1->AVZ_NORM_FACTOR_MODE == 1) { // use new (fixed) normalization factors
    
        double dF = 1. / ((double)(NFOUR/2) * TIMESTEP); // in Hz
        //cout<<"dF1 : "<<dF<<endl;
            
        double dF_org = detector->GetFreq(1) - detector->GetFreq(0); // in Hz
            
        //double dF_factor = sqrt( dF_org / dF );
         double dF_factor = 1.;
         //double dF_factor = sqrt( dF / dF_org );
            
         int Norg = detector->GetFreqBin();
                                    
         double FreqOrg[Norg+1];
         double VmMHzOrg[Norg+1];
                                            
         double FreqNFOUR[NFOUR/4+1]; // one more bin
         double VmMHzNFOUR[NFOUR/4+1]; // one more bin
         
         for (int i=0;i<Norg+1;i++) {
                                                        
             if ( i==0 ) {
                                                                       
                 VmMHzOrg[i] = 0.;
                 FreqOrg[i] = 0.;
             }
             else {
                                                                                                              
                 VmMHzOrg[i] = vsignal_array[i-1];
                 FreqOrg[i] = detector->GetFreq(i-1);
             }
                                                                                                                                      
         }
                                                                                                                                          
         for (int ifour=0;ifour<NFOUR/4+1;ifour++) {
                                                                                                                                                   
             FreqNFOUR[ifour] = dF * (double)ifour;
             VmMHzNFOUR[ifour] = 0.;
                                                                                                                                                                
         }
         
           
         //Tools::SimpleLinearInterpolation_OutZero( Norg, FreqOrg, vsignal_array, NFOUR/4+1, FreqNFOUR, VmMHzNFOUR );
         Tools::SimpleLinearInterpolation_OutZero( Norg+1, FreqOrg, VmMHzOrg, NFOUR/4+1, FreqNFOUR, VmMHzNFOUR );
          
           
         for (int ifour=1;ifour<NFOUR/4;ifour++) {
           
             // same amplitude, 2/(NFOUR/2) for inverse FFT normalization factor
             vsignal_forfft[2*ifour] = VmMHzNFOUR[ifour] * 2/((double)NFOUR/2) * dF_factor / TIMESTEP; // change to V/Hz, apply norm factor, change to Hn
             vsignal_forfft[2*ifour+1] = VmMHzNFOUR[ifour] * 2/((double)NFOUR/2) * dF_factor / TIMESTEP;
          
             // apply phase
             vsignal_forfft[2*ifour]*=cos(settings1->PHASE*PI/180.);
             vsignal_forfft[2*ifour+1]*=sin(settings1->PHASE*PI/180.);
         }
           
         // first and last freq bin read values to 0, 1 bin
         vsignal_forfft[0] = VmMHzNFOUR[0] * 2/((double)NFOUR/2) * dF_factor / TIMESTEP;
         vsignal_forfft[1] = VmMHzNFOUR[NFOUR/4] * 2/((double)NFOUR/2) * dF_factor / TIMESTEP;
    }



    
}




void Report::MakeArraysforFFT(Settings *settings1, Detector *detector, int StationIndex, double *vsignal_array, double *vsignal_forfft) {



    // from icemc, anita class MakeArraysforFFT
    int NFOUR = settings1->NFOUR;
    double TIMESTEP = settings1->TIMESTEP;

//    int NFOUR = detector->stations[StationIndex].NFOUR;
//    int TIMESTEP = detector->stations[StationIndex].TIMESTEP;
    Tools::Zero(vsignal_forfft,NFOUR/2);


    if (settings1->AVZ_NORM_FACTOR_MODE == 0) { // use previous normalization factors

    
        double previous_value_e_even=0.;
        double previous_value_e_odd=0.;
        int count_nonzero=0;
        int iprevious=0;
        int ifirstnonzero=-1;
        int ilastnonzero=2000;
        //for (int i=0;i<NFREQ;i++) {
        for (int i=0;i<detector->GetFreqBin();i++) {
            
            // freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
            // but there are only NFOUR/4 different values
            // it's the index among the NFOUR/4 that we're interested in
            int ifour=Tools::Getifreq(detector->GetFreq(i),detector->freq_forfft[0],detector->freq_forfft[NFOUR/2-1],NFOUR/4);
            
            if (ifour!=-1 && 2*ifour+1<NFOUR/2) {
                count_nonzero++;
                if (ifirstnonzero==-1)
                    ifirstnonzero=ifour;
                
                vsignal_forfft[2*ifour]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP); // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft
                
                //      cout << "ifour, vsignal is " << ifour << " " << vsignal_e_forfft[2*ifour] << "\n";
                
                vsignal_forfft[2*ifour+1]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP); // phase is 90 deg.
                // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting
                
                // how about we interpolate instead of doing a box average
                
                for (int j=iprevious+1;j<ifour;j++) {
                    vsignal_forfft[2*j]=previous_value_e_even+(vsignal_forfft[2*ifour]-previous_value_e_even)*(double)(j-iprevious)/(double)(ifour-iprevious);
                    //	cout << "j, vsignal is " << j << " " << vsignal_e_forfft[2*j] << "\n";
                    
                    vsignal_forfft[2*j+1]=previous_value_e_odd+(vsignal_forfft[2*ifour+1]-previous_value_e_odd)*(double)(j-iprevious)/(double)(ifour-iprevious);
                }
                
                ilastnonzero=ifour;
                iprevious=ifour;
                previous_value_e_even=vsignal_forfft[2*ifour];
                previous_value_e_odd=vsignal_forfft[2*ifour+1];
            }
            
        } // end loop over nfreq
        


        // don't applying any extra factor for the change in array (change of bin size)
        // as change in the bin size doesn't matter for the total energy
        // total energy is just same as integral over frequency range and this frequency range will not change
        //
          //cout << "ifirstnonzero, ilastnonzero are " << ifirstnonzero << " " << ilastnonzero << "\n";
          //cout << "non zero count is " << count_nonzero << "\n";
          //cout << "ratio is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
        for (int j=0;j<NFOUR/4;j++) {
            vsignal_forfft[2*j]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
            vsignal_forfft[2*j+1]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
        }



        
        //  Tools::InterpolateComplex(vsignal_e_forfft,NFOUR/4);
        //Tools::InterpolateComplex(vsignal_h_forfft,NFOUR/4);
        for (int ifour=0;ifour<NFOUR/4;ifour++) {



            vsignal_forfft[2*ifour]*=cos(settings1->PHASE*PI/180.);
            vsignal_forfft[2*ifour+1]*=sin(settings1->PHASE*PI/180.);
            
            //--------------------------------------------------
            // if (!PULSER) {
            //     
            //     vsignal_forfft[2*ifour]*=cos(phase*PI/180.);
            //     vsignal_forfft[2*ifour+1]*=sin(phase*PI/180.);
            //     
            //     
            // }
            // else {
            //     vsignal_forfft[2*ifour]*=cos(v_phases[ifour]*PI/180.);
            //     vsignal_forfft[2*ifour+1]*=sin(v_phases[ifour]*PI/180.);
            //     
            // }	  	  	  
            //-------------------------------------------------- 
            
            
        }
    }
    else if (settings1->AVZ_NORM_FACTOR_MODE == 1) { // use new (fixed) normalization factors
    
        double dF = 1. / ((double)(NFOUR/2) * TIMESTEP); // in Hz
        //cout<<"dF1 : "<<dF<<endl;
            
        double dF_org = detector->GetFreq(1) - detector->GetFreq(0); // in Hz
            
        //double dF_factor = sqrt( dF_org / dF );
         double dF_factor = 1.;
         //double dF_factor = sqrt( dF / dF_org );
            
         int Norg = detector->GetFreqBin();
                                    
         double FreqOrg[Norg+1];
         double VmMHzOrg[Norg+1];
                                            
         double FreqNFOUR[NFOUR/4+1]; // one more bin
         double VmMHzNFOUR[NFOUR/4+1]; // one more bin
         
         for (int i=0;i<Norg+1;i++) {
                                                        
             if ( i==0 ) {
                                                                       
                 VmMHzOrg[i] = 0.;
                 FreqOrg[i] = 0.;
             }
             else {
                                                                                                              
                 VmMHzOrg[i] = vsignal_array[i-1];
                 FreqOrg[i] = detector->GetFreq(i-1);
             }
                                                                                                                                      
         }
                                                                                                                                          
         for (int ifour=0;ifour<NFOUR/4+1;ifour++) {
                                                                                                                                                   
             FreqNFOUR[ifour] = dF * (double)ifour;
             VmMHzNFOUR[ifour] = 0.;
                                                                                                                                                                
         }
         
           
         //Tools::SimpleLinearInterpolation_OutZero( Norg, FreqOrg, vsignal_array, NFOUR/4+1, FreqNFOUR, VmMHzNFOUR );
         Tools::SimpleLinearInterpolation_OutZero( Norg+1, FreqOrg, VmMHzOrg, NFOUR/4+1, FreqNFOUR, VmMHzNFOUR );
          
           
         for (int ifour=1;ifour<NFOUR/4;ifour++) {
           
             // same amplitude, 2/(NFOUR/2) for inverse FFT normalization factor
             vsignal_forfft[2*ifour] = VmMHzNFOUR[ifour] * 2/((double)NFOUR/2) * dF_factor / TIMESTEP; // change to V/Hz, apply norm factor, change to Hn
             vsignal_forfft[2*ifour+1] = VmMHzNFOUR[ifour] * 2/((double)NFOUR/2) * dF_factor / TIMESTEP;
          
             // apply phase
             vsignal_forfft[2*ifour]*=cos(settings1->PHASE*PI/180.);
             vsignal_forfft[2*ifour+1]*=sin(settings1->PHASE*PI/180.);
         }
           
         // first and last freq bin read values to 0, 1 bin
         vsignal_forfft[0] = VmMHzNFOUR[0] * 2/((double)NFOUR/2) * dF_factor / TIMESTEP;
         vsignal_forfft[1] = VmMHzNFOUR[NFOUR/4] * 2/((double)NFOUR/2) * dF_factor / TIMESTEP;
    }



}






void Report::MakeArraysforFFT_noise(Settings *settings1, Detector *detector, int StationIndex, vector <double> &vsignal_array, double *vsignal_forfft) {



    // from icemc, anita class MakeArraysforFFT
    int NFOUR = settings1->NFOUR;
    double TIMESTEP = settings1->TIMESTEP;

//    int NFOUR = detector->stations[StationIndex].NFOUR;
//    int TIMESTEP = detector->stations[StationIndex].TIMESTEP;
    
    Tools::Zero(vsignal_forfft,NFOUR/2);
    
    double previous_value_e_even=0.;
    double previous_value_e_odd=0.;
    int count_nonzero=0;
    int iprevious=0;
    int ifirstnonzero=-1;
    int ilastnonzero=2000;
    //for (int i=0;i<NFREQ;i++) {
    for (int i=0;i<detector->GetFreqBin();i++) {
	
	// freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
	// but there are only NFOUR/4 different values
	// it's the index among the NFOUR/4 that we're interested in
	int ifour=Tools::Getifreq(detector->GetFreq(i),detector->freq_forfft[0],detector->freq_forfft[NFOUR/2-1],NFOUR/4);
	
	if (ifour!=-1 && 2*ifour+1<NFOUR/2) {
	    count_nonzero++;
	    if (ifirstnonzero==-1)
		ifirstnonzero=ifour;
	    
	    vsignal_forfft[2*ifour]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP); // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft
	    
	    //      cout << "ifour, vsignal is " << ifour << " " << vsignal_e_forfft[2*ifour] << "\n";
	    
	    vsignal_forfft[2*ifour+1]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP); // phase is 90 deg.
	    // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting
	    
	    // how about we interpolate instead of doing a box average
	    
	    for (int j=iprevious+1;j<ifour;j++) {
		vsignal_forfft[2*j]=previous_value_e_even+(vsignal_forfft[2*ifour]-previous_value_e_even)*(double)(j-iprevious)/(double)(ifour-iprevious);
		//	cout << "j, vsignal is " << j << " " << vsignal_e_forfft[2*j] << "\n";
		
		vsignal_forfft[2*j+1]=previous_value_e_odd+(vsignal_forfft[2*ifour+1]-previous_value_e_odd)*(double)(j-iprevious)/(double)(ifour-iprevious);
	    }
	    
	    ilastnonzero=ifour;
	    iprevious=ifour;
	    previous_value_e_even=vsignal_forfft[2*ifour];
	    previous_value_e_odd=vsignal_forfft[2*ifour+1];
	}
	
    } // end loop over nfreq
    


    // don't applying any extra factor for the change in array (change of bin size)
    // as change in the bin size doesn't matter for the total energy
    // total energy is just same as integral over frequency range and this frequency range will not change
    //
      //cout << "ifirstnonzero, ilastnonzero are " << ifirstnonzero << " " << ilastnonzero << "\n";
      //cout << "non zero count is " << count_nonzero << "\n";
      //cout << "ratio is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
    for (int j=0;j<NFOUR/4;j++) {
	vsignal_forfft[2*j]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
	vsignal_forfft[2*j+1]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
    }


    // phase for noise will be applied in GetNoiseWaveforms

    
}




                                       



double Report::FindPeak(double *waveform, int n) {  // same with icemc, trigger-> AntTrigger::FindPeak

    double peak=0.;

    for (int i=0;i<n;i++) {
        if (fabs(waveform[i])>peak)
            peak = fabs(waveform[i]);
    }
    return peak;
}

void Report::SetRank(Detector *detector) {
    // think about this!
    // test 1

    int current = 0;    // current checking rank
    int check;  // check is any change in rank
    double maxpeak; // maxpeak value
    double pre_maxpeak; // maxpeak at previous round


    //cout<<"SetRank step1) : find all PeakV == 0"<<endl;
        for (int i = 0; i< detector->params.number_of_stations; i++) {

            for (int j=0; j< detector->stations[i].strings.size(); j++) {

                for (int k=0; k< detector->stations[i].strings[j].antennas.size(); k++) {

                    if (stations[i].strings[j].antennas[k].ray_sol_cnt) {
                    if (stations[i].strings[j].antennas[k].PeakV[0] == 0.) {
                        stations[i].strings[j].antennas[k].Rank.push_back(0);  // rank 0, PeakV is 0, non-countable rank.
                    }
                    else {
                        stations[i].strings[j].antennas[k].Rank.push_back(current+1);  // elses, if PeakV is not 0, set Rank as 1 (for now)
                    }
                    }   // if ray_sol_cnt != 0
                }
            }
        }


    //cout<<"Start while loop for Ranking!!"<<endl;

    check = 1;
    pre_maxpeak = 1.E5; // unreasonably big value which real PeakV will never reach
    while (check!=0) {
        check=0;
        maxpeak = 0.;
        for (int i = 0; i< detector->params.number_of_stations; i++) {

            for (int j=0; j< detector->stations[i].strings.size(); j++) {

                for (int k=0; k< detector->stations[i].strings[j].antennas.size(); k++) {

                    //for (int l=0; l< detector->stations[i].strings[j].antennas[k].ray_sol_cnt; l++) {
                    //
                    //lets start with 1st ray_sol only
                    //

                    if (stations[i].strings[j].antennas[k].ray_sol_cnt) {
                    if (stations[i].strings[j].antennas[k].Rank[0] != 0) {  // there is non-zero value! and ray_sol_cnt also non-zero!
                    if (stations[i].strings[j].antennas[k].PeakV[0] < pre_maxpeak) {

                        if (maxpeak < stations[i].strings[j].antennas[k].PeakV[0] ) {
                            maxpeak = stations[i].strings[j].antennas[k].PeakV[0];

                            stations[i].strings[j].antennas[k].Rank[0] = current+1;
                            check++;
                        } // is maxpeak < PeakV
                        else {
                            stations[i].strings[j].antennas[k].Rank[0] = current + 2;
                        } // else

                    } // if rank > current
                    }
                    }

                } // for antennas
            } // for strings
        } // for stations

        current++;
        pre_maxpeak = maxpeak;


    }   // while


    //cout<<"END Ranking!!"<<endl;



}



int Report::GetChannelNum8_LowAnt(int string_num, int antenna_num) {

    int outputch = 16;

    // just give ch numbers 1-8 for antenna 0 - 1 while higher ch numbers for antenna 2 - 3
    //
    if ( string_num == 0 && antenna_num == 0 ) {
        outputch = 1;
    }
    else if ( string_num == 0 && antenna_num == 1 ) {
        outputch = 2;
    }
    else if ( string_num == 1 && antenna_num == 0 ) {
        outputch = 3;
    }
    else if ( string_num == 1 && antenna_num == 1 ) {
        outputch = 4;
    }
    else if ( string_num == 2 && antenna_num == 0 ) {
        outputch = 5;
    }
    else if ( string_num == 2 && antenna_num == 1 ) {
        outputch = 6;
    }
    else if ( string_num == 3 && antenna_num == 0 ) {
        outputch = 7;
    }
    else if ( string_num == 3 && antenna_num == 1 ) {
        outputch = 8;
    }

    // higher antennas
    else if ( string_num == 0 && antenna_num == 2 ) {
        outputch = 9;
    }
    else if ( string_num == 0 && antenna_num == 3 ) {
        outputch = 10;
    }
    else if ( string_num == 1 && antenna_num == 2 ) {
        outputch = 11;
    }
    else if ( string_num == 1 && antenna_num == 3 ) {
        outputch = 12;
    }
    else if ( string_num == 2 && antenna_num == 2 ) {
        outputch = 13;
    }
    else if ( string_num == 2 && antenna_num == 3 ) {
        outputch = 14;
    }
    else if ( string_num == 3 && antenna_num == 2 ) {
        outputch = 15;
    }
    else if ( string_num == 3 && antenna_num == 3 ) {
        outputch = 16;
    }

    return outputch;

}

                        


                        

TGraph *Report::getWaveform(Detector *detector, int ch, int station_i, int event_num, int run_num){
 
  int string_i = detector->getStringfromArbAntID( station_i, ch);
  int antenna_i = detector->getAntennafromArbAntID( station_i, ch);
  
  TGraph *gr=new TGraph();
  
  int N=stations[station_i].strings[string_i].antennas[antenna_i].V_mimic.size();
  
  for(int i=0;i<N;i++){
   
    gr->SetPoint(i,stations[station_i].strings[string_i].antennas[antenna_i].time_mimic[i], stations[station_i].strings[string_i].antennas[antenna_i].V_mimic[i]);
    
  }// for i
  
  gr->SetNameTitle(Form("WF%dS%02dC%02d", event_num, station_i, ch), Form("Simulated waveform %d, station %d, channel %d;time (ns); voltage (mV)", event_num, station_i, ch));
  
  return gr;
  
}


vector<TGraph*> Report::getWaveformVector(Detector *detector, int station_i, int event_num, int run_num){

  vector<TGraph*> waveforms;
  
  for(int ch=0;ch<detector->stations[station_i].number_of_antennas; ch++){
      
    int string_i = detector->getStringfromArbAntID( station_i, ch);
    int antenna_i = detector->getAntennafromArbAntID( station_i, ch);
    
    TGraph *gr=new TGraph();
      
    int N=stations[station_i].strings[string_i].antennas[antenna_i].V_mimic.size();
    
    for(int i=0;i<N;i++){
   
      gr->SetPoint(i,stations[station_i].strings[string_i].antennas[antenna_i].time_mimic[i], stations[station_i].strings[string_i].antennas[antenna_i].V_mimic[i]);
    
    }// for i
  
    gr->SetNameTitle(Form("WF%dS%02dC%02d", event_num, station_i, ch), Form("Simulated waveform %d, station %d, channel %d;time (ns); voltage (mV)", event_num, station_i, ch));
    
    waveforms.push_back(gr);
    
  }// for ch
  
  return waveforms;
  
}

vector<TGraph*> Report::getWaveformVectorVpol(Detector *detector, int station_i, int event_num, int run_num){
 
  vector<TGraph*> waveforms;
  
  for(int ch=0;ch<detector->stations[station_i].number_of_antennas; ch++){
      
    int string_i = detector->getStringfromArbAntID( station_i, ch);
    int antenna_i = detector->getAntennafromArbAntID( station_i, ch);
    
    if(detector->stations[station_i].strings[string_i].antennas[antenna_i].type!=0) continue; // jump to next channel if this isn't Vpol
    
    TGraph *gr=new TGraph();
      
    int N=stations[station_i].strings[string_i].antennas[antenna_i].V_mimic.size();
    
    for(int i=0;i<N;i++){
   
      gr->SetPoint(i,stations[station_i].strings[string_i].antennas[antenna_i].time_mimic[i], stations[station_i].strings[string_i].antennas[antenna_i].V_mimic[i]);
    
    }// for i
  
    gr->SetNameTitle(Form("WF%dS%02dC%02d", event_num, station_i, ch), Form("Simulated waveform %d, station %d, channel %d;time (ns); voltage (mV)", event_num, station_i, ch));
    
    waveforms.push_back(gr);
    
  }// for ch
  
  return waveforms;
  
}

vector<TGraph*> Report::getWaveformVectorHpol(Detector *detector, int station_i, int event_num, int run_num){
 
  vector<TGraph*> waveforms;
  
  for(int ch=0;ch<detector->stations[station_i].number_of_antennas; ch++){
      
    int string_i = detector->getStringfromArbAntID( station_i, ch);
    int antenna_i = detector->getAntennafromArbAntID( station_i, ch);
    
    if(detector->stations[station_i].strings[string_i].antennas[antenna_i].type!=1) continue; // jump to next channel if this isn't Hpol
    
    TGraph *gr=new TGraph();
      
    int N=stations[station_i].strings[string_i].antennas[antenna_i].V_mimic.size();
    
    for(int i=0;i<N;i++){
   
      gr->SetPoint(i,stations[station_i].strings[string_i].antennas[antenna_i].time_mimic[i], stations[station_i].strings[string_i].antennas[antenna_i].V_mimic[i]);
    
    }// for i
  
    gr->SetNameTitle(Form("WF%dS%02dC%02d", event_num, station_i, ch), Form("Simulated waveform %d, station %d, channel %d;time (ns); voltage (mV)", event_num, station_i, ch));
    
    waveforms.push_back(gr);
    
  }// for ch
  
  return waveforms;
  
}


vector<double> Report::getHitTimesVector(Detector *detector, int station_i, int polarization){

  vector<double> times;
    
  for(int ch=0;ch<detector->stations[station_i].number_of_antennas;ch++){
    
    int string_i = detector->getStringfromArbAntID( station_i, ch);
    int antenna_i = detector->getAntennafromArbAntID( station_i, ch);
    
    if(polarization>=0 && detector->stations[station_i].strings[string_i].antennas[antenna_i].type!=polarization) continue; // jump to next channel if this isn't Hpol/Vpol

    
    if(stations[station_i].strings[string_i].antennas[antenna_i].arrival_time.size()){
   
      times.push_back(1e9*stations[station_i].strings[string_i].antennas[antenna_i].arrival_time[0]);// get the direct beam arrival time
    
    }
  
    else times.push_back(-1000); // if there's no ray-trace solution just put something in.
  
  }
  
  return times;
  
}

int Report::getNumOfSignalledAnts(Station_r station){
    // Count and return the number of antennas in the provided station that 
    //   have waveforms, from signal only, that exceed 0.0

    int ants_with_good_signal = 0;
    vector <double> summary;

    for ( // loop over each string in the station
        int s = 0; s < station.strings.size(); s++
    ) { 
        for ( // loop over each antenna on the string
            int a = 0; a < station.strings[s].antennas.size(); a++
        ) { 

            double max_signal = 0.0; 

            // If there are ray solutions to this antenna, get the max 
            //   (of the absolute) value of their signal-only waveforms
            if ( station.strings[s].antennas[a].ray_sol_cnt > 0 ){

                for ( // Loop over ray solutions
                    int ray=0; ray<station.strings[s].antennas[a].ray_sol_cnt; ray++
                ){
                    // Get max (of the absolute) value in waveform
                    double ray_max_signal = Tools::getMaxAbsoluteMagnitude(
                        station.strings[s].antennas[a].V[ray]
                    );
                    // Check if the max (abs) value is greater than 
                    //   the signals from the other rays
                    if ( ray_max_signal > max_signal ){
                        max_signal = ray_max_signal;
                    }
                } 

            } // end if (ray solutions exist)

            // If this antenna had a non-zero signal-only signal, 
            //   increment the good antenna tracker
            if ( max_signal > 0.0){ 
                ants_with_good_signal++;
                summary.push_back(max_signal);
            }

        } // end loop over antennas on string
    } // end loop over strings on the antenna

    return ants_with_good_signal;
  
}

double Report::get_SNR(vector<double> signal_array, vector<double> noise_array){
    // Returns the SNR of the provided signal compared to the provided noise

    // if signal_array or noise_array have a size of 0, tell stderr
    if ( signal_array.size() == 0 ){
        cerr<<"Signal array provided to Report::get_SNR has size 0. ";
        cerr<<"Function will return inaccurate values."<<endl;
        throw("Cancelling Simulation");
    }
    if ( noise_array.size() == 0 ){
        cerr<<"Noise array provided to Report::get_SNR has size 0. ";
        cerr<<"Function will return inaccurate values."<<endl;
        throw("Cancelling Simulation");
    }

    // Get peak-to-peak of signal
    double min_value =  10000000.;
    double max_value = -10000000.;
    for (int i=0; i<signal_array.size(); i++){
        if ( signal_array[i] < min_value ) min_value = signal_array[i];
        else if ( signal_array[i] > max_value ) max_value = signal_array[i];
    }
    double p2p = max_value - min_value;

    // Get RMS
    double rms = 0;
    for (int i=0; i<noise_array.size(); i++){
        rms += noise_array[i] * noise_array[i];
    }
    rms = sqrt( rms / noise_array.size() );

    // Catch potential errors
    if ( (p2p==-20000000.) || (p2p!=p2p) ){ 
        // if p2p is original value or nan
        cerr<<"In Report::get_SNR, peak-to-peak not calculated correctly. ";
        cerr<<"Will set peak-to-peak to 0."<<endl;
        p2p = 0;
    }
    if ( (rms==0.) || (rms!=rms) ){ 
        // if RMS is 0 or nan
        cerr<<"In Report::get_SNR, RMS was not calculated properly. ";
        cerr<<"Will set RMS to 1 to avoid divide-by-zero."<<endl;
        rms = 1;
    }

    // Calculate and return SNR
    double snr = p2p / 2. / rms;
    return snr;

}

vector<double> Report::getHitTimesVectorVpol(Detector *detector, int station_i){

  return getHitTimesVector(detector, station_i, 0);
  
}

vector<double> Report::getHitTimesVectorHpol(Detector *detector, int station_i){
   
  return getHitTimesVector(detector, station_i, 1);
  
}


bool Report::isTrigger(double eff){
  if (eff >= 1.0) return true;
  float randNum = gRandom->Rndm();
  if (randNum < eff) return true;
  return false;
}

void Report::checkPATrigger(
    int i, double all_receive_ang[2], double &viewangle, int ray_sol_cnt,
    Detector *detector, Event *event, int evt, Trigger *trigger, Settings *settings1, 
    int trig_search_init, int max_total_bin
){
    // Calculates max SNR in topmost PA Vpol 
    //   and multplies it by viewing angle factors. 
    // Then determines signal efficiency by interpolating oberved SNR from 
    //   efficiency vs SNR data.
    // Triggers if signal efficiency is above certain treshold.
    // If triggered, saves relevant information.

    //cout <<"successfully made it to PA Trigger!" << endl;
    
    // For phased array, waveform length is 680 ns, but 
    // for PA trigger check only 20 ns around the signal bin.
    // This is to avoid getting the second ray
    // KAH and ARB are not sure where the 1200 number comes from
    int BINSIZE = 1200/(settings1->TIMESTEP*1.e9);  // Number of bins (aka datapoints) of data to save

    int waveformLength = settings1->WAVEFORM_LENGTH;
    int waveformCenter = settings1->WAVEFORM_CENTER;
    int raySolNum = 0;
    int dsignalBin = 0;
    viewangle=viewangle*180.0/PI;
    bool searchSecondRay = true;
    if (ray_sol_cnt == 2){
        dsignalBin = abs(signal_bin[0] - signal_bin[1]); //original kah
        dsignalBin = abs(stations[i].strings[0].antennas[8].SignalBin[0]-stations[i].strings[0].antennas[8].SignalBin[1]);
    }
    bool hasTriggered = false;

    // If the antenna we use to trigger the PA doesn't have any ray solutions, 
    //   do not perform the trigger check.
    if (stations[i].strings[0].antennas[8].ray_sol_cnt == 0){
        return;
    }

    //cout << "time to trigger " << endl;
    while(raySolNum < stations[i].strings[0].antennas[8].SignalBin.size()){

        int signalbinPA = stations[i].strings[0].antennas[8].SignalBin[raySolNum]; //new kah
        int bin_value;
        double ant_SNR;
        if(settings1->TRIG_ANALYSIS_MODE == 2) { // Noise only triggers
            ant_SNR = pa_force_trigger_snr; 
        }
        // // ARB 7/7/23: 
        // //    I don't trust the signal+noise trigger estimator right now 
        // //    so it's commented out below. This is because the trigger 
        // //    is determined based on PA trigger efficiency that was
        // //    calculated with signal-only triggers IIRC. This means
        // //    using signal+noise will artificially inflate the trigger rate.
        // //    I think we'll need to develop trigger efficiency as a function
        // //    of signal+noise data or actually build the actual physical
        // //    triggering pipeline (forming beams, etc) before we can 
        // //    run PA simulations over signal+noise. 
        // //    This is a project for the future. 
        // else if (settings1->TRIG_ANALYSIS_MODE==1) // Noise + signal triggers
        //     // Estimate average SNR in topmost vpol
        //     if(stations[i].strings[0].antennas[8].V.size()>raySolNum) {   
        //
        //         // Steal the noise RMS from the trigger class and pass 
        //         //   it as the noise WF to get_SNR() (since the RMS of an 
        //         //   array with one element is the absolute value of that element)
        //         vector <double> tmp_noise_RMS;
        //         int trigger_ch_ID = GetChNumFromArbChID(detector, 8, i, settings1) - 1;
        //         double ant_noise_voltage_RMS = trigger->GetAntNoise_voltageRMS(trigger_ch_ID, settings1);
        //         tmp_noise_RMS.push_back( ant_noise_voltage_RMS );
        //
        //         // Calculate SNR in this antenna
        //         ant_SNR = get_SNR( 
        //             stations[i].strings[0].antennas[8].V_convolved, 
        //             tmp_noise_RMS);
        //
        //     }
        //     else {
        //         ant_SNR = 0.0;
        //     }
        else { // signal only triggers
            // Estimate average SNR in topmost vpol
            if(stations[i].strings[0].antennas[8].V.size()>raySolNum) {   

                // Steal the noise RMS from the trigger class and pass 
                //   it as the noise WF to get_SNR() (since the RMS of an 
                //   array with one element is the absolute value of that element)
                vector <double> tmp_noise_RMS;
                int trigger_ch_ID = GetChNumFromArbChID(detector, 8, i, settings1) - 1;
                double ant_noise_voltage_RMS = trigger->GetAntNoise_voltageRMS(trigger_ch_ID, settings1);
                tmp_noise_RMS.push_back( ant_noise_voltage_RMS );

                // Calculate SNR in this antenna
                ant_SNR = get_SNR( 
                    stations[i].strings[0].antennas[8].V_convolved, 
                    tmp_noise_RMS);
                
            }
            else {
                ant_SNR = 0.0;
            }
        }

        // Cap SNR
        if(ant_SNR>pa_snr_cap) ant_SNR = pa_snr_cap;

        //scale snr to reflect the angle
        all_receive_ang[raySolNum] = all_receive_ang[raySolNum]*180.0/PI-90.0;
        double snr_50 = interpolate(
            trigger->angle_PA, trigger->aSNR_PA, // x and y coordinates of curve to interpolate
            all_receive_ang[raySolNum], // x value to interpolate y value for
            (*(&trigger->angle_PA+1) - trigger->angle_PA) - 1 // len(ang_data) - 1
        );
        ant_SNR = ant_SNR*2.0/snr_50;

        // Estimate the PA signal efficiency of for this SNR from curve of efficiency vs data
        double eff = interpolate(
            trigger->snr_PA, trigger->eff_PA, // x and y coordinates of curve to interpolate
            ant_SNR, // x value to interpolate y value for
            (*(&trigger->snr_PA+1) - trigger->snr_PA) - 1 // len(snr_PA) - 1
        ); 
        
        if(ant_SNR > 0.5){
                   
            if(isTrigger(eff)){ // if a randomly selected value is greater than the PA efficiency we calculated, this event triggers
                cout<<endl<<"PA trigger ~~~ raySolNum: "<< raySolNum;
                cout<<"  SNR: "<<ant_SNR<<"  Event Number : "<<evt;
                cout<<"  PA efficiency : "<<eff<<endl;

                if (hasTriggered) {
                    cout<<"Weight for Second Ray trigger is: "<<event->Nu_Interaction[0].weight<<endl;
                    break;
                }
                viewAngle = viewangle;
                my_averageSNR = ant_SNR;
                my_raysol = raySolNum;
                my_receive_ang = all_receive_ang[raySolNum];
                int last_trig_bin = signalbinPA;
                int my_ch_id = 0;
                stations[i].Global_Pass = last_trig_bin;
                for (size_t str = 0; str < detector->stations[i].strings.size(); str++) {
                    for (size_t ant = 0; ant < detector->stations[i].strings[str].antennas.size(); ant++) {
                        double peakvalue = 0;
                        for (int bin=0; bin<BINSIZE; bin++) {

                            bin_value = signalbinPA - BINSIZE/2 + bin;
                            stations[i].strings[str].antennas[ant].V_mimic.push_back(trigger->Full_window_V[my_ch_id][bin_value]);// save in V (kah)
                            stations[i].strings[str].antennas[ant].time.push_back( bin_value );
                            stations[i].strings[str].antennas[ant].time_mimic.push_back( ( bin) * settings1->TIMESTEP*1.e9 );// save in ns
                            if (TMath::Abs(trigger->Full_window_V[ant][bin_value]) > peakvalue) {
                                peakvalue = TMath::Abs(trigger->Full_window_V[my_ch_id][bin_value]);
                            }
                            
                        }//end bin
                        my_ch_id ++;
                        //cout<<" Peak Value for ant "<<ant<<" is "<<peakvalue<<endl;
                    }//end ant
                }//end detector

                hasTriggered = true;

            }//end efficiency if

        }//end avgsnr if

        // Save information like in triggerCheckLoop() and saveTriggeredEvent()
        if (hasTriggered==true){

            int numChan=stations[i].TDR_all.size();
            int numChanVpol=stations[i].TDR_Vpol_sorted.size();
            int numChanHpol=stations[i].TDR_Hpol_sorted.size();

            double powerthreshold = 2.0;
            double Pthresh_value[numChan];
            CircularBuffer **buffer=new CircularBuffer*[numChan];

            int trig_window_bin = (int)(settings1->TRIG_WINDOW / settings1->TIMESTEP);  // coincidence window bin for trigger

            int SCTR_cluster_bit[numChan];
            for(int trig_j=0;trig_j<numChan;trig_j++) SCTR_cluster_bit[trig_j]=0;

            double *TDR_all_sorted_temp;
            double *TDR_Vpol_sorted_temp;
            double *TDR_Hpol_sorted_temp;
            
            if(settings1->TRIG_SCAN_MODE>1){ // prepare TDR storage arrays and initialize all values to 0
                TDR_all_sorted_temp=new double[numChan];
                TDR_Vpol_sorted_temp=new double[numChanVpol];
                TDR_Hpol_sorted_temp=new double[numChanHpol];
                for(int trig_j=0; trig_j<numChan; trig_j++) {
                    TDR_all_sorted_temp[trig_j] =0;
                }
                for(int trig_j=0; trig_j<numChanVpol; trig_j++) {
                    TDR_Vpol_sorted_temp[trig_j]=0;
                }
                for(int trig_j=0; trig_j<numChanHpol; trig_j++) {
                    TDR_Hpol_sorted_temp[trig_j]=0;  
                }
            } // if scan_mode>1

            int check_TDR_configuration=0; // check if we need to reorder our TDR arrays
            int first_trigger=0;

            // Calculate PThresh information (taken from triggerCheckLoop)
            for(int trig_j=0;trig_j<numChan; trig_j++){// initialize Trig_Pass and buffers for each channel
                int string_i = detector->getStringfromArbAntID( i, trig_j);
                int antenna_i = detector->getAntennafromArbAntID( i, trig_j);

                Pthresh_value[trig_j]=0;
                buffer[trig_j]=new CircularBuffer(trig_window_bin, powerthreshold, settings1->TRIG_SCAN_MODE);

                stations[i].strings[string_i].antennas[antenna_i].SingleChannelTriggers=0;
                stations[i].strings[string_i].antennas[antenna_i].TotalBinsScannedPerChannel=0;

            } // end channel loop 

            // Loop over window bins, get PThresh and decide if we need to update sorted arrays
            int N_pass_V = 0;
            int window_pass_bit = 0; // Whether this trig_i window passes
            int bin_to_save_on = -1; 
            for(int trig_i = trig_search_init; trig_i < max_total_bin; trig_i++) { // scan the different window positions
                
                // FOR EACH CHANNEL
                for (int trig_j=0; trig_j<numChan; trig_j++){
                
                    int string_i = detector->getStringfromArbAntID( i, trig_j);
                    int antenna_i = detector->getAntennafromArbAntID( i, trig_j);

                    int channel_num = detector->GetChannelfromStringAntenna( 5, string_i, antenna_i, settings1 );

                    // assign Pthresh a value 
                    double diode_noise_RMS = trigger->GetAntNoise_diodeRMS(channel_num-1, settings1);
                    Pthresh_value[trig_j] = (
                        trigger->Full_window[trig_j][trig_i] / 
                        ( diode_noise_RMS  * detector->GetThresOffset(i, channel_num-1,settings1) )
                    );

                    // this is to count how many local trigger clusters there are 
                    if(Pthresh_value[trig_j]<powerthreshold){
            
                        if(SCTR_cluster_bit[trig_j]==0) stations[i].strings[string_i].antennas[antenna_i].SingleChannelTriggers++;
                    
                        // records all the different Pthresh values that caused local trigger.
                        if(settings1->TRIG_SCAN_MODE>2){  // save PThresh to ant SCT_threshold_pass if its the best or first
                            if(SCTR_cluster_bit[trig_j]==0){// if first trigger in cluster
                                stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.push_back(Pthresh_value[trig_j]);
                            }
                            else{// choose the highest trigger value (most negative) in cluster
                                if(Pthresh_value[trig_j]<stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.back()) stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.back()=Pthresh_value[trig_j];
                            }
                        }// trig scan mode > 2
                        
                        SCTR_cluster_bit[trig_j]=1;
            
                    }// if local trigger
                    else SCTR_cluster_bit[trig_j]=0;// if no local trigger, set zero to start a new cluster at next local trigger
                    
                    // and how many bins scanned
                    stations[i].strings[string_i].antennas[antenna_i].TotalBinsScannedPerChannel++;

                    // fill the buffers (if any changes occur mark check_TDR_configuration as non-zero)
                    if(trig_i<trig_search_init+trig_window_bin) check_TDR_configuration+=buffer[trig_j]->fill(Pthresh_value[trig_j]);
                    else check_TDR_configuration+=buffer[trig_j]->add(Pthresh_value[trig_j]);

                    // if there is at least one value above threshold in the buffer, this is ++
                    if ( 
                        buffer[trig_j]->addToNPass>0 && 
                        string_i == 0 && 
                        detector->stations[i].strings[string_i].antennas[antenna_i].type == 0
                    ) {
                        N_pass_V++;
                    }
                
                }// end iteration over channels to calculate PThresh

                // Add in PThresh calculators inspired by triggerCheckLoop()
                if ( N_pass_V > 3 ) {

                    window_pass_bit=1;

                    // if this is the first trigger, mark this position and save event
                    if(first_trigger==0){ 
                        first_trigger=1;
                        // FOR EACH CHANNEL
                        for(int trig_j=0;trig_j<numChan;trig_j++){ // Set Trig_Pass for each ant
                            int string_i = detector->getStringfromArbAntID( i, trig_j);
                            int antenna_i = detector->getAntennafromArbAntID( i, trig_j);
                            if(buffer[trig_j]->addToNPass>0) stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = trig_i-buffer[trig_j]->numBinsToOldestTrigger(); // mark the bin on which we triggered...
                            else stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = 0.;
                        }// for trig_j
                        bin_to_save_on = trig_i;
                    }// first trigger

                } // end if N_Pass_Vpol > some value

                if(settings1->TRIG_SCAN_MODE>1&&check_TDR_configuration&&window_pass_bit){// if there's a trigger and anything changes in the buffers, restock the TDR arrays
                        
                    for(int trig_j=0;trig_j<numChan;trig_j++) TDR_all_sorted_temp[trig_j]=0;
                    for(int trig_j=0;trig_j<numChanVpol;trig_j++) TDR_Vpol_sorted_temp[trig_j]=0;
                    for(int trig_j=0;trig_j<numChanHpol;trig_j++) TDR_Hpol_sorted_temp[trig_j]=0;

                    for(int trig_j=0;trig_j<numChan;trig_j++){// fill the TDR (unsorted) arrays if they improved... 
                        if(buffer[trig_j]->best_value<stations[i].TDR_all[trig_j]) stations[i].TDR_all[trig_j]=buffer[trig_j]->best_value;
                    }// for trig_j

                    // Changes TDR sorting and buffer[best_chan] for Vpol only:
                    for(int ii=0;ii<numChanVpol; ii++){// find the best channel's TDR and store them.

                        double best_thresh=0;
                        int best_chan=0;
                        
                        // Pick a new best_thresh and best_chan
                        for(int trig_j=0;trig_j<numChan;trig_j++){
                            int string_i = detector->getStringfromArbAntID( i, trig_j);
                            int antenna_i = detector->getAntennafromArbAntID( i, trig_j);  
                            if(detector->stations[i].strings[string_i].antennas[antenna_i].type == 0&&buffer[trig_j]->temp_value<best_thresh){ 
                                best_thresh=buffer[trig_j]->temp_value; 
                                best_chan=trig_j;
                            }// if best
                        }// for trig_j

                        buffer[best_chan]->temp_value=0;
                        TDR_Vpol_sorted_temp[ii]=best_thresh;
                    
                    }// end sort Vpol tdr values
                    
                    // Changes TDR sorting and buffer[best_chan] for Hpol only
                    for(int ii=0;ii<numChanHpol; ii++){// find the best channel's TDR and store them.

                        double best_thresh=0;
                        int best_chan=0;
                        
                        // Pick a new best_chan and best_thresh
                        for(int trig_j=0;trig_j<numChan;trig_j++){ 
                            int string_i = detector->getStringfromArbAntID( i, trig_j);
                            int antenna_i = detector->getAntennafromArbAntID( i, trig_j);
                            if(detector->stations[i].strings[string_i].antennas[antenna_i].type==1&&buffer[trig_j]->temp_value<best_thresh){ 
                                best_thresh=buffer[trig_j]->temp_value; 
                                best_chan=trig_j;
                            }// if best
                        }// for trig_j
                        
                        buffer[best_chan]->temp_value=0;
                        TDR_Hpol_sorted_temp[ii]=best_thresh;
                
                    }// end sort Hpol tdr values

                    // Update TDR arrays (no matter what)
                    if(settings1->TRIG_MODE==0){
                        for(int ii=0;ii<stations[i].TDR_all_sorted.size();ii++) {
                            if(TDR_all_sorted_temp[ii]<stations[i].TDR_all_sorted[ii]) {
                                stations[i].TDR_all_sorted[ii]=TDR_all_sorted_temp[ii];
                        } }  
                    }
                    if(settings1->TRIG_MODE==1){
                        for(int ii=0;ii<stations[i].TDR_Vpol_sorted.size();ii++) {
                            if(TDR_Vpol_sorted_temp[ii]<stations[i].TDR_Vpol_sorted[ii]) {
                                stations[i].TDR_Vpol_sorted[ii]=TDR_Vpol_sorted_temp[ii];
                        } } 
                        for(int ii=0;ii<stations[i].TDR_Hpol_sorted.size();ii++) {
                            if(TDR_Hpol_sorted_temp[ii]<stations[i].TDR_Hpol_sorted[ii]) {
                                stations[i].TDR_Hpol_sorted[ii]=TDR_Hpol_sorted_temp[ii];
                        } } 
                    }

                }// if trigger and buffer changed

            } // end iteration over windows/bins

            if ( bin_to_save_on == -1 ) bin_to_save_on = signalbinPA;
            
            // Do what saveTriggeredEvent() does
            for(int trig_j=0; trig_j<numChan;trig_j++){

                int string_i = detector->getStringfromArbAntID( i, trig_j);
                int antenna_i = detector->getAntennafromArbAntID( i, trig_j);
                    
                stations[i].strings[string_i].antennas[antenna_i].Likely_Sol = -1; // no likely init
                    
                int mindBin = 1.e9; // init big values
                int dBin = 0;
                    
                for (int m=0; m<stations[i].strings[string_i].antennas[antenna_i].ray_sol_cnt; m++) {   // loop over raysol numbers
                    if ( stations[i].strings[string_i].antennas[antenna_i].SignalExt[m] ) {
                        dBin = abs( 
                            stations[i].strings[string_i].antennas[antenna_i].SignalBin[m] 
                            - stations[i].strings[string_i].antennas[antenna_i].Trig_Pass 
                        );
                        if ( dBin < mindBin ) {
                            stations[i].strings[string_i].antennas[antenna_i].Likely_Sol = m; // store the ray sol number which is minimum difference between Trig_Pass bin
                            mindBin = dBin;    
                        }
                    }     
                } // for m (ray sol numbers)
                        
                // set global_trig_bin values
                if (settings1->V_MIMIC_MODE == 0) { // Global passed bin is the center of the window
                stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = waveformLength/2 + waveformCenter ;
                }
                else if (settings1->V_MIMIC_MODE == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
                stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin) + waveformLength/2 + waveformCenter ;
                }
                else if (settings1->V_MIMIC_MODE == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
                    stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) + waveformLength/2 + waveformCenter ;
                }
                
                double arrivtime = stations[i].strings[string_i].antennas[antenna_i].arrival_time[0];
                double X = detector->stations[i].strings[string_i].antennas[antenna_i].GetX();
                double Y = detector->stations[i].strings[string_i].antennas[antenna_i].GetY();
                double Z = detector->stations[i].strings[string_i].antennas[antenna_i].GetZ();
                
                stations[i].total_trig_search_bin = stations[i].Global_Pass + trig_window_bin - trig_search_init;


            } // end for trig_j in numchans

        } // end if hasTriggered==true
        
        raySolNum++;
        if(hasTriggered==true) break;
        if(searchSecondRay == false) break;

        //if(stations[i].Global_Pass > 0) break;
    }//while ray solve

}

double Report::interpolate(double *xdata,double *ydata, double xi, int numData)
{
    double result = 0.0; // Initialize result
    double x1, y1, x2, y2; // Data points about which linear interpolation will take place
    double c, m; //slope and constant in y = mx + c
    //cout<<"Last Value "<<xdata[numData]<<endl;
    if (xi >= xdata[numData]) return ydata[numData];
    if (xi <= xdata[0]) return ydata[0];
    for (int i=0; i<numData; i++)
    {
        if (i == numData - 1){
          result = ydata[numData];
        }
        if (xi > xdata[i] && xi < xdata[i+1]){
          x1 = xdata[i]; x2 = xdata[i+1]; y1 = ydata[i]; y2 = ydata[i+1];
          c = (y2*x1 - x2*y1)/(x1-x2);
          m = (y1-y2)/(x1-x2);
          result = m*xi + c; //linear interpolation
          break;
        }
        if (xdata[i] == x1) result = ydata[i];
    }
    //cout<<"x1 "<<x1<<" x2 "<<x2<<" xi "<<xi<<endl;

    return result;
}

//Adding function for padding waveforms to take FFT


