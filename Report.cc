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

}

void Report::delete_all() {

  stations.clear();
  strings.clear();
  Passed_chs.clear();
  Vfft_noise_after.clear();
  Vfft_noise_before.clear();
  noise_phase.clear();
  Passed_chs.clear();

}

void Report::Initialize(Detector *detector, Settings *settings1) {

  // Clear information stored in this object (but there shouldn't be. just to make sure)
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

      // vector antennas
      for (int k=0; k<detector->stations[i].strings[j].antennas.size(); k++)
        stations[i].strings[j].antennas.push_back(tmp_antenna);
    }

    // vector surface antennas
    for (int j=0; j<detector->stations[i].surfaces.size(); j++)
      stations[i].surfaces.push_back(tmp_surface);

    if(settings1->TRIG_SCAN_MODE>0){// scan Pthresh mode

      int numChan=0;
      int numChanVpol=0;
      int numChanHpol=0;

      for(int j=0;j<detector->stations[i].strings.size(); j++){

        for (int k=0;k<detector->stations[i].strings[j].antennas.size();k++) {

          int string_i = detector->getStringfromArbAntID( i, numChan);
          int antenna_i = detector->getAntennafromArbAntID( i, numChan);

          if (detector->stations[i].strings[string_i].antennas[antenna_i].type == 0)
            numChanVpol++;

          if (detector->stations[i].strings[string_i].antennas[antenna_i].type == 1)
            numChanHpol++;

          numChan++;
        }// for k

      }// for j

      stations[i].TDR_all.clear();
      for(int ch=0;ch<numChan; ch++ )
        stations[i].TDR_all.push_back(0);

      stations[i].TDR_all_sorted.clear();
      if(settings1->TRIG_MODE==0) {
        for(int ch=0;ch<numChan; ch++ )
          stations[i].TDR_all_sorted.push_back(0);
      }

      stations[i].TDR_Vpol_sorted.clear();
      if(settings1->TRIG_MODE==1) {
        for(int ch=0;ch<numChanVpol; ch++)
          stations[i].TDR_Vpol_sorted.push_back(0);
      }

      stations[i].TDR_Hpol_sorted.clear();
      if(settings1->TRIG_MODE==1) {
        for(int ch=0;ch<numChanHpol; ch++ )
          stations[i].TDR_Hpol_sorted.push_back(0);
      }

    }// if TRIG_SCAN_MODE

  }// for i (number of stations)

}



void Antenna_r::Prepare_Outputs(int n_interactions) {
  // Initialize all vectors used to save event information so that they have one
  //   index per interaction in this event.
  // If any objects are added to Antenna_r, they may need to be added here!

  view_ang.resize(n_interactions);
  launch_ang.resize(n_interactions);
  phi_launch.resize(n_interactions);
  theta_launch.resize(n_interactions);
  rec_ang.resize(n_interactions);
  reflect_ang.resize(n_interactions);
  Dist.resize(n_interactions);
  L_att.resize(n_interactions);
  arrival_time.resize(n_interactions);
  reflection.resize(n_interactions);
  Pol_vector.resize(n_interactions);
  vmmhz.resize(n_interactions);
  Heff.resize(n_interactions);
  Heff_copol.resize(n_interactions);
  Heff_crosspol.resize(n_interactions);
  Mag.resize(n_interactions);
  Fresnel.resize(n_interactions);
  Pol_factor.resize(n_interactions);
  Pol_factorH.resize(n_interactions);
  Pol_factorV.resize(n_interactions);
  phi_rec.resize(n_interactions);
  theta_rec.resize(n_interactions);
  Vfft.resize(n_interactions);
  Vfft_noise.resize(n_interactions);

  ray_step.resize(n_interactions);

  V.resize(n_interactions);

  PeakV.resize(n_interactions);

  Trig_Pass = 0;

  Vm_zoom.resize(n_interactions);
  Vm_zoom_T.resize(n_interactions);

  SignalBin.resize(n_interactions);
  SignalBinTime.resize(n_interactions);
  SignalExt.resize(n_interactions);

}



void Antenna_r::clear() {
  // Empty all vectors used to save event information.
  // If any objects are added to Antenna_r, they may need to be added here!

  view_ang.clear();
  launch_ang.clear();
  phi_launch.clear();
  theta_launch.clear();
  rec_ang.clear();
  reflect_ang.clear();
  Dist.clear();
  L_att.clear();
  arrival_time.clear();
  reflection.clear();
  Pol_vector.clear();
  vmmhz.clear();
  Heff.clear();
  Heff_copol.clear();
  Heff_crosspol.clear();
  Mag.clear();
  Fresnel.clear();
  Pol_factor.clear();
  Pol_factorH.clear();
  Pol_factorV.clear();
  phi_rec.clear();
  theta_rec.clear();
  Vfft.clear();
  Vfft_noise.clear();

  ray_step.clear();

  time.clear();
  time_mimic.clear();
  V_noise.clear();
  V_convolved.clear();
  V_mimic.clear();

  V.clear();

  noise_ID.clear();

  PeakV.clear();
  Rank.clear();
  TooMuch_Tdelay.clear();

  Vm_zoom.clear();
  Vm_zoom_T.clear();

  SignalBin.clear();
  SignalBinTime.clear();
  SignalExt.clear();

  SCT_threshold_pass.clear();

  // Reset variables to values we've designated as "null" values
  global_trig_bin = -1;
  Likely_Sol[0] = -1;
  Likely_Sol[1] = -1;
  Nnew[0] = 1;
  Nnew[1] = 1;
  ray_sol_cnt = 0;
  SingleChannelTriggers = 0;
  skip_bins[0] = 1;
  skip_bins[1] = 1;
  TotalBinsScannedPerChannel = 0;
  Trig_Pass = 0;

}



void Antenna_r::clear_useless(Settings *settings1) {
  // If the user requests that less information be stored in the output file,
  //   clear some information here from the Antenna_r object

  // DATA_SAVE_MODE meanings from Settings.h
  // 0 : save all information generated by the Report class
  // 1 : light mode, remove most generated data except for physics event info
  //     (vertex, direction, energy, etc) and the final waveform
  // 2 : light mode 2, only save physics information (energy, geometric info),
  //     Remove all waveform information (include noise waveforms in trigger class)

  // Remove data from intermediate waveform objects
  if (settings1->DATA_SAVE_MODE >= 1) {

    // Clear EM signal and waveform data for each interaction and ray
    vmmhz.clear();
    Heff.clear();
    Heff_copol.clear();
    Heff_crosspol.clear();
    Vfft.clear();
    Vfft_noise.clear();
    V.clear();

    TooMuch_Tdelay.clear();

    // additional for before ant waveform
    Vm_zoom.clear();
    Vm_zoom_T.clear();
    V_convolved.clear();
    V_noise.clear();

  }

  // Remove final waveform data and ray path data
  if (settings1->DATA_SAVE_MODE >= 2) {

    //! clear the ray step to reduce the size of output AraOut.root
    ray_step.clear();

    // clear global trigger waveform info also
    time.clear();
    time_mimic.clear();
    V_mimic.clear();

  }

  return;
}


int CircularBuffer::add(double input_value){

  changelog=0;
  temp_value=best_value;

  if(buffer[i]<pthresh)
    addToNPass--; // if the value leaving the buffer is over threshold, we reduce the counter.
  if(mode>1)
    last_value=buffer[i]; // in mode 1 we don't care about the values, just about the addToNPass

  if(input_value<pthresh){
    addToNPass++; // if value entering buffer is over threshold we increase the counter.
    buffer[i]=input_value;
  }
  else
    buffer[i]=0; // if under threshold just insert zero

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
  if(i==N)
    i=0;

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
  }
  else
    buffer[i]=0; // if under threshold add just zero

  if(mode>1&&buffer[i]<best_value){ // improve best threshold value in the buffer
    best_value=buffer[i];
    temp_value=best_value; changelog=1;
  }

  i++;
  if(i==N)
    i=0;

  return changelog;

}// fill

double CircularBuffer::findBestValue(){

  double temp_best=0;
  for(int ii=0;ii<N;ii++) {
    if ( buffer[ii]<temp_best )
      temp_best=buffer[ii];
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
      if(buffer[(i-1+j)%N]<0)
        break; // if there's any value here, its because it passed the threshold, so take it
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

    if(bin<0)
      bin=N+(i-j);

    if(buffer[bin]<0)// if there's any value here, its because it passed the threshold, so take it
      break;
  }// for j

  return j; // the count of how many bins between i and the latest trigger...

}// numBinsToLatestTrigger



void Report::clear_useless(Settings *settings1) {
  // If the user requests that less information be stored in the output file,
  //   clear some information here from the Report object

  if (settings1->DATA_SAVE_MODE != 0) {

      noise_phase.clear();
      Passed_chs.clear();
      Vfft_noise_after.clear();
      Vfft_noise_before.clear();
      V_total_forconvlv.clear();
      RayStep.clear();

  }

  return;
}


void Report::CalculateSignals(
    int debugmode,
    Birefringence *birefringence, Detector *detector, Event *event,
    IceModel *icemodel, RaySolver *raysolver, Settings *settings1, Signal *signal
){
  // For each each station, perform the ray tracing from each cascade to each antenna
  //   and determine the voltage readout in the antenna for each connected ray

  double RandomTshift = gRandom->Rndm(); // for t-domain signal, a factor for random init time shift

  // Loop over stations
  for (int i = 0; i < detector->stations.size(); i++)
  {

    double min_arrival_time_tmp = 10.; // min arrival time between all antennas, raysolves, first min_arrival_time is unreasonably big value
    double max_arrival_time_tmp = 0.;  // max arrival time between all antennas, raysolves, first max_arrival_time is unreasonably small value
    double max_PeakV_tmp = 0.;         // max PeakV of all antennas in the station, first max_PeakV_tmp is 0.

    stations[i].Total_ray_sol = 0;  // initial Total_ray_sol value

    // Loop over strings on this station
    for (int j = 0; j < detector->stations[i].strings.size(); j++)
    {

      // Create the array used for interpolation
      double init_T = settings1->TIMESTEP *-1.e9 *((double) settings1->NFOUR / 4 + RandomTshift); // locate zero time at the middle and give random time shift
      double T_forint[settings1->NFOUR / 2];
      for (int n = 0; n < settings1->NFOUR / 2; n++)
        T_forint[n] = init_T + (double) n *settings1->TIMESTEP *1.e9; // in ns

      // Loop over antennas on this string
      for (int k = 0; k < detector->stations[i].strings[j].antennas.size(); k++) {

        stations[i].strings[j].antennas[k].clear(); // clear data in antenna which stored in previous event
        stations[i].strings[j].antennas[k].Prepare_Outputs(event->Nu_Interaction.size()); // Resize output arrays to the number of interactions

        // Loop over interactions in this event
        for (int interaction_idx=0; interaction_idx<event->Nu_Interaction.size(); interaction_idx++) {

          // Perform the ray tracing between this interaction (interaction_idx)
          //   and the antenna (k) then calculate the voltage response in the antenna

          int ray_sol_cnt = 0;
          vector<vector < double>> ray_output;

          // Check that the event is in a valid location and close enough to the antenna before raytracing
          if (  event->Nu_Interaction[interaction_idx].pickposnu && // If the event location is in Antarctic Ice
                event->Nu_Interaction[interaction_idx].posnu.Distance(detector->stations[i].strings[j].antennas[k]) <= settings1->RAYSOL_RANGE
                  // If the event is closer to the antenna than the user-defined maximum distance for ray tracing (settings1->RAYSOL_RANGE)
          ) {

            // Clear and calculate all ray solutions between the cascade and antenna
            RayStep.clear();
            raysolver->Solve_Ray(
                event->Nu_Interaction[interaction_idx].posnu, detector->stations[i].strings[j].antennas[k],
                icemodel, ray_output, settings1, RayStep
            );

            // If there are ray tracing solutions from the cascade to the antenna
            if (raysolver->solution_toggle) {

              // Loop over the number of ray solutions (could be 1 or 2)
              while (ray_sol_cnt < ray_output[0].size()) {

                // Calculate parameters of this ray tracing solution such as
                //   attenuation factors, antenna viewing angles, polarization, etc
                ModelRay(
                    ray_sol_cnt, ray_output, interaction_idx, T_forint,
                    &stations[i].strings[j].antennas[k],
                    &detector->stations[i].strings[j].antennas[k],
                    i, j, k, debugmode,
                    birefringence, detector, event, icemodel, settings1, signal
                );

                // check max_PeakV
                if ( debugmode == 0){
                  if (max_PeakV_tmp < stations[i].strings[j].antennas[k].PeakV[interaction_idx][ray_sol_cnt])
                    max_PeakV_tmp = stations[i].strings[j].antennas[k].PeakV[interaction_idx][ray_sol_cnt];
                }
                // check min_arrival_time
                if (min_arrival_time_tmp > stations[i].strings[j].antennas[k].arrival_time[interaction_idx][ray_sol_cnt])
                  min_arrival_time_tmp = stations[i].strings[j].antennas[k].arrival_time[interaction_idx][ray_sol_cnt];

                // check max_arrival_time
                if (max_arrival_time_tmp < stations[i].strings[j].antennas[k].arrival_time[interaction_idx][ray_sol_cnt])
                  max_arrival_time_tmp = stations[i].strings[j].antennas[k].arrival_time[interaction_idx][ray_sol_cnt];

                ray_sol_cnt++;

              }   // end while number of solutions

            }   // end if solution exist

          }   // end if event position is acceptable

          else {
            ray_sol_cnt = 0;
          }

          stations[i].strings[j].antennas[k].ray_sol_cnt += ray_sol_cnt;   // save number of RaySolver solutions
          stations[i].Total_ray_sol += ray_sol_cnt;   // add ray_sol_cnt to Total_ray_sol

        }

      }   // for number_of_antennas_string

    }   // for number_of_strings_station

    // set each station's min/max arrival time
    stations[i].min_arrival_time = min_arrival_time_tmp;
    stations[i].max_arrival_time = max_arrival_time_tmp;

    // set each station's max PeakV
    stations[i].max_PeakV = max_PeakV_tmp;
  }   // for number_of_stations

  // If the simulation is not in debug mode, set ranking of signal between
  //   antennas after all values are stored in Report
  if (debugmode == 0){
    SetRank(detector);
  }

  return;
}

void Report::BuildAndTriggerOnWaveforms(
    int debugmode, int station_index, int evt, int trig_search_init,
    Detector *detector, Event *event, Settings *settings1, Trigger *trigger
){
  // Loop over all antennas to make DATA_BIN_SIZE array for signal + noise. (with time delay)
  // With that signal + noise array, we'll do the convolution with the diode response.
  // With the convolution result, we'll do a trigger check.

  // When this variable is non-zero, it indicates that the station triggered
  // Reinitialize it back to 0
  stations[station_index].Global_Pass  = 0;

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
    if (stations[station_index].Total_ray_sol)
        ants_with_nonzero_signal = getNumOfSignalledAnts(stations[station_index]);
  }

  if (stations[station_index].Total_ray_sol && ants_with_nonzero_signal)
  {
    // if there is any ray_sol (don't check trigger if there is no ray_sol at all)

    // calculate total number of bins we need to do trigger check
    // to save time, use only necessary number of bins
    // the max number of bins should be (max_arrival_bin - min_arrival_bin) + size of a signal
    // we define the size of a signal as BINSIZE = NFOUR/2 later in the code
    // there could be two ray solutions in an interaction. We add 2*BINSIZE = NFOUR for safety
    int max_total_bin = (
        (stations[station_index].max_arrival_time - stations[station_index].min_arrival_time)
        / settings1->TIMESTEP) + settings1->NFOUR;

    // Compute next power of two larger than max_total_bin
    max_total_bin = (int)std::pow(2, std::ceil(std::log2(max_total_bin)));

    //Analyze more bins for historically unspecified reasons.
    //Could be to adjust for an offset in indexing and signals that come after the last event temporally
    //However AraSim may break without this extra padding, so adjust carefully
    //This offset is probably to account for offset between indices of V_convolved and V_signal in GetNoiseThenConvolve()
    max_total_bin += settings1->NFOUR *3 + trigger->maxt_diode_bin;

    // Decide if new noise waveforms will be generated for every new event.
    // Otherwise, the same noise waveforms will be used for the entire simulation run.
    if (settings1->NOISE_WAVEFORM_GENERATE_MODE == 0) {
      // Generate brand new noise waveforms for each event

      // redefine DATA_BIN_SIZE
      int DATA_BIN_SIZE_tmp;
      int DBS = 10;  // Starting power of 2
      while (true) {
        DATA_BIN_SIZE_tmp = (int) pow(2., (double) DBS);
        if (DATA_BIN_SIZE_tmp > max_total_bin)
          break; // Exit the loop if the size exceeds max_total_bin

        DBS++; // Increment DBS for the next power of 2
      }
      settings1->DATA_BIN_SIZE = DATA_BIN_SIZE_tmp;

      // Resize Full_window and Full_window_V to size 16 for 16 channels
      trigger->Full_window.resize(16);
      trigger->Full_window_V.resize(16);
      for (int i = 0; i < 16; i++) {
        trigger->Full_window[i].clear();
        trigger->Full_window[i].resize(DATA_BIN_SIZE_tmp, 0.0);

        trigger->Full_window_V[i].clear();
        trigger->Full_window_V[i].resize(DATA_BIN_SIZE_tmp, 0.0);
      }

      // reset all filter values in Detector class
      detector->get_NewDiodeModel(settings1);
      detector->ReadFilter_New(settings1);
      detector->ReadPreamp_New(settings1);
      detector->ReadFOAM_New(settings1);

      if (settings1->USE_TESTBED_RFCM_ON == 1)
        detector->ReadRFCM_New(settings1);
      if (settings1->NOISE == 1 && settings1->DETECTOR == 3)
        detector->ReadRayleigh_New(settings1);

      // reset Trigger class noise temp values
      trigger->Reset_V_noise_freqbin(settings1, detector);

      // now call new noise waveforms with new DATA_BIN_SIZE
      trigger->GetNewNoiseWaveforms(settings1, detector, this);

    }

    // If the simulation is not in debug mode, check if DATA_BIN_SIZE is large
    //   enough for the total time delay between antennas via this N_noise
    //   variable. It should evaluate to 1 in most cases.
    if (debugmode == 1) {
        N_noise = 1;
    }
    else {
        N_noise = (int)(max_total_bin / settings1->DATA_BIN_SIZE) + 1;
    }
    if (N_noise > 1) {
        cout << "N_noise : " << N_noise << " max_total_bin : " << max_total_bin << " might cause error!!" << endl;
    }

    // For all antennas on all strings, get the noise waveform and combine
    //   it with the signal waveforms from all connected rays then convolve the
    //   full waveform through the tunnel diode.
    int ch_ID = 0;
    int ants_with_sufficient_SNR = 0;
    for (int j = 0; j < detector->stations[station_index].strings.size(); j++)
    {
      // Loop over strings

      for (int k = 0; k < detector->stations[station_index].strings[j].antennas.size(); k++)
      {
        // Loop over antennas

        // Fill trigger->Full_window and trigger->Full_window_V with noise waveforms
        Prepare_Antenna_Noise(debugmode, ch_ID, station_index, j, k, settings1, trigger, detector);

        // currently there is a initial spoiled bins (maxt_diode_bin)
        // at the initial Full_window "AND" at the initial of connected noisewaveform
        // (we can fix this by adding values but not accomplished yet)

        // If the simulation is not in debug mode, combine the voltage response
        //   from all rays into a single vector, add the noise waveform,
        //   and convolve this complete waveform with the tunnel diode.
        if (debugmode == 0) {
            Convolve_Signals(
                &stations[station_index].strings[j].antennas[k], ch_ID, station_index,
                event, settings1, trigger, detector
            );
        }

        // Calculate the SNR in this antenna from the event signal only
        double ant_SNR = 0.;
        if (settings1->TRIG_ANALYSIS_MODE==2)
          // TRIG_ANALYSIS_MODE=2 is the noise-only mode
          // Set the SNR to an arbitarily high number so it passes
          //   the coming insufficient signal check
          ant_SNR = 100.;
        else {

          // Use the noise RMS from the trigger class and pass
          //   it as the noise WF to get_SNR() (since the RMS of an
          //   array with one element is the absolute value of that element)
          vector <double> tmp_noise_RMS;
          int trigger_ch_ID = GetChNumFromArbChID(detector, ch_ID, station_index, settings1) - 1;
          double ant_noise_voltage_RMS = trigger->GetAntNoise_voltageRMS(trigger_ch_ID, settings1);
          tmp_noise_RMS.push_back( ant_noise_voltage_RMS );

          // Calculate the SNR in this antenna from the signal only
          ant_SNR = get_SNR(
              stations[station_index].strings[j].antennas[k].V_convolved,
              tmp_noise_RMS);

        }

        // Log if this antenna has a decent signal-only SNR or not.
        if ( ant_SNR > 0.01 )
          ants_with_sufficient_SNR++;

        // Apply gain factors
        if ((debugmode == 0) && (settings1->USE_MANUAL_GAINOFFSET == 1) )
          Apply_Gain_Offset(settings1, trigger, detector, ch_ID, station_index);

        ch_ID++; // now to next channel

      } // for antennas

    } // for strings

    // Check for a trigger on this event
    // Do only if it's not in debugmode and if at least 1 antenna has sufficient signal
    // Previously, events with very low neutrino signals were triggering on
    //   noise and inflating effective volume calculations. By skipping
    //   events with insufficient signal in all antennas, we reduce
    //   contributions from noise-only triggered events in our effective
    //   volume calculations which require a pure sample of signal-only
    //   triggered events.
    if ( (debugmode == 0) && (ants_with_sufficient_SNR) ) {

      // coincidence window bin for trigger
      int trig_window_bin = (int)(settings1->TRIG_WINDOW / settings1->TIMESTEP);

      // save trig search bin info default (if no global trig, this value saved)
      stations[station_index].total_trig_search_bin = max_total_bin - trig_search_init;

      if (settings1->TRIG_SCAN_MODE==5) {
        // Trigger mode for phased array
        checkPATrigger(
            station_index, detector, event, evt, trigger, settings1,
            trig_search_init, max_total_bin
        );
      }
      else if (settings1->TRIG_SCAN_MODE == 0) {
        // Original and default AraSim trigger mode
        triggerCheck_ScanMode0(
            trig_search_init, max_total_bin, trig_window_bin,
            ch_ID, station_index, detector, event, settings1, trigger
        );
      }
      else if (settings1->TRIG_SCAN_MODE > 0) {
        // Newer trigger check loops that use the Circular Buffer and
        //   save some extra information. More information in log.txt.
        triggerCheckLoop(
            settings1, detector, event, trigger,
            station_index, trig_search_init, max_total_bin, trig_window_bin,
            settings1->TRIG_SCAN_MODE
        );
      }

    }   // if it's not debugmode

    // delete noise waveforms
    if (settings1->NOISE_WAVEFORM_GENERATE_MODE == 0)
      // noise waveforms will be generated for each evts
      // remove noise waveforms for next evt
      trigger->ClearNoiseWaveforms();

  }   // if there is any ray_sol in the station

  // Check if there's additional waveform to trigger on. Save this event and rerun trigger if so.
  stations[station_index].next_trig_search_init = -1;
  if (stations[station_index].Global_Pass) {

    bool analyze_more_waveform = false;

    for (int string=0; string<stations[station_index].strings.size(); string++){
      for (int antenna=0; antenna<stations[station_index].strings[string].antennas.size(); antenna++){

        // Calculate the next bin that could be checked for a trigger based on this antenna's trigger bin
        int this_next_trig_search_init = (
            stations[station_index].Global_Pass // trigger bin in V_convolved
            + (settings1->WAVEFORM_LENGTH - stations[station_index].strings[string].antennas[antenna].global_trig_bin)
                // number of bins between trigger bin and end of readout window.
                // Sometimes includes channel-specific delays
            + settings1->DEADTIME/settings1->TIMESTEP // deadtime in units of bins
        );

        // Save this antenna's "next initial trigger bin" as the whole station's
        //   "next initial trigger bin" if it's larger than what was previously saved
        if (this_next_trig_search_init > stations[station_index].next_trig_search_init)
          stations[station_index].next_trig_search_init = this_next_trig_search_init;
      }
    }

    for (int string=0; string<stations[station_index].strings.size(); string++){

      for (int antenna=0; antenna<stations[station_index].strings[string].antennas.size(); antenna++){

        // Find the last bin in V_convolved that has nonzero signal
        int signal_length = stations[station_index].strings[string].antennas[antenna].V_convolved.size();
        int last_nonzero_bin = 0;
        for (int bin=signal_length-1; bin>-1; bin--) {
          if ( stations[station_index].strings[string].antennas[antenna].V_convolved[bin] != 0. ) {
            last_nonzero_bin = bin;
            break;
          }
        }

        // If there is signal that exists later than what would be the next trig_search_init,
        //   indicate that we need to analyze more signal.
        if ( last_nonzero_bin > stations[station_index].next_trig_search_init )
          analyze_more_waveform = true;

      }
    }

    // If we have no waveform left to analyze, reset the `next_trig_search_init` to `-1`
    //   to indicate that there's nothing to trigger check after reading out this event so far.
    if ( !analyze_more_waveform )
      stations[station_index].next_trig_search_init = -1;
  }

  return;
}

void Antenna_r::Find_Likely_Sol(){
  // Estimate which interaction and which ray from said interaction triggered
  //   the detector based on the timing between the bin that triggered and the
  //   SignalBin of each ray from each simulated cascade.

  int mindBin = 1.e9; // init big value
  int dBin = 0;

  Likely_Sol[0] = -1;
  Likely_Sol[1] = -1;

  for (int interaction_index=0; interaction_index<SignalBin.size(); interaction_index++) { // Loop over interactions
    for (int ray_index = 0; ray_index < SignalBin[interaction_index].size(); ray_index++) { // loop over raysol numbers
      if (SignalExt[interaction_index][ray_index]) {
        dBin = abs(SignalBin[interaction_index][ray_index] - Trig_Pass);
        if (dBin < mindBin) {
          // store the ray sol number which is minimum difference between Trig_Pass bin
          Likely_Sol[0] = interaction_index;
          Likely_Sol[1] = ray_index;
          mindBin = dBin;
        }
      }
    }
  }

  return;
}

void Antenna_r::Get_Brightest_Interaction(int (*brightest_event)[2]){
  // Determine which interaction and which ray from that interaction yielded
  //   the highest voltage response in this antenna.

  brightest_event[0][0] = -1;
  brightest_event[0][1] = -1;
  double maxV = -1;
  for (int interaction_idx=0; interaction_idx<V.size(); interaction_idx++){
    for (int ray=0; ray<V[interaction_idx].size(); ray++){
      double this_maxV = Tools::getMaxAbsoluteMagnitude(V[interaction_idx][ray]);
      if (this_maxV > maxV){
        maxV = this_maxV;
        brightest_event[0][0] = interaction_idx;
        brightest_event[0][1] = ray;
      }
    }
  }

  return;
}


void Report::ModelRay(
    int ray_idx, vector< vector< double > > ray_output, int interaction_idx, double *T_forint,
    Antenna_r *antenna_r, Antenna *antenna_d,  int i, int j, int k,
    int debugmode,  Birefringence *birefringence, Detector *detector,
    Event *event, IceModel *icemodel, Settings *settings, Signal *signal
){
    // Get the voltage response in an antenna along one ray path from a cascade by
    //   - Calculating the electric field signal at the cascade
    //   - Transforming that into the signal at the antenna by modifying the
    //       field according to attenuation, fresnel, etc along the ray path
    //   - Convolving the signal at the antenna with the antenna's effective height,
    //       gain, and other parameters to obtain the voltage response to the cascade

    // This (gain_ch_no) is used for per-channel gain implementation.
    // It is used in all instances of ApplyElect_Tdomain() and ApplyElect_Tdomain_FirstTwo(), to indicate channel number
    // Note that channel numbering is different for DETECTOR==4 than for the other modes (1-3). See that in the definition of GetChannelfromStringAntenna()
    int gain_ch_no;
    if (settings->DETECTOR==4 || settings->DETECTOR==5) {
        gain_ch_no = detector->GetChannelfromStringAntenna (i, j, k, settings)-1;
    }
    else {
        gain_ch_no = detector->GetChannelfromStringAntenna (i, j, k, settings);
    }

    double time_diff_birefringence = birefringence->Time_Diff_TwoRays(
        RayStep[ray_idx][0], RayStep[ray_idx][1], ray_output[3][ray_idx],
        event->Nu_Interaction[interaction_idx].posnu_from_antcen, settings
    ); // calculate time differences for birefringence

    antenna_r->arrival_time[interaction_idx].push_back(ray_output[4][ray_idx] + event->interactions_birth_time[interaction_idx]);

    //! Save every ray steps between the vertex (source) and an antenna (target), unless DATA_SAVE_MODE is 2.
    //! These xz coordinates were calculated after we convert the earth coordinates to flat coordinates by the RaySolver::Earth_to_Flat_same_angle()
    antenna_r->ray_step[interaction_idx].resize(ray_idx + 1); ///< resize by number of ray solutions
    antenna_r->ray_step[interaction_idx][ray_idx].resize(2); ///< resize by xz values
    for (int steps = 0; steps < (int) RayStep[ray_idx][0].size(); steps++) {
        ///< push back each ray step coordinates
        antenna_r->ray_step[interaction_idx][ray_idx][0].push_back(RayStep[ray_idx][0][steps]);
        antenna_r->ray_step[interaction_idx][ray_idx][1].push_back(RayStep[ray_idx][1][steps]);
    }

    // get ice attenuation factor
    double IceAttenFactor = 1.;
    if (settings->USE_ARA_ICEATTENU == 1) { // use new ARA measured ice attenuation values
        double dx, dz, dl;
        for (int steps = 1; steps < (int) RayStep[ray_idx][0].size(); steps++) {

            dx = RayStep[ray_idx][0][steps - 1] - RayStep[ray_idx][0][steps];
            dz = RayStep[ray_idx][1][steps - 1] - RayStep[ray_idx][1][steps];
            dl = sqrt((dx *dx) + (dz *dz));

            // Skipping attenuation calculation when the distance between two RaySteps is 0.
            // PrevenTing adds -nan into the IceAttenFactor.
            if (dl > 0) {
                // use new ice model
                // use the midpoint of the array to calculate the attenuation length, instead of the end of the ray
                IceAttenFactor *= (
                    exp(-dl / icemodel->GetARAIceAttenuLength(-RayStep[ray_idx][1][steps])) +
                    exp(-dl / icemodel->GetARAIceAttenuLength(-RayStep[ray_idx][1][steps - 1]))
                ) / 2;
            }

        }
    }
    else if (settings->USE_ARA_ICEATTENU == 0) { // use old method
        IceAttenFactor = exp(-ray_output[0][ray_idx] / icemodel->EffectiveAttenuationLength(settings, event->Nu_Interaction[interaction_idx].posnu, 0));
    }

    // If the simulation is not in debug mode, calculate the field from the
    //   cascade and transform it into a voltage read out by the antenna
    if (debugmode == 0) {

        Vector n_trg_pokey; // unit pokey vector at the target
        Vector n_trg_slappy; // unit slappy vector at the target
        Vector Pol_vector_src; // Polarization at the source since Pol_vector is polarization vector at antenna
        Position launch_vector; // direction of ray at the source
        Position receive_vector; // direction of ray at the target antenna
        GetRayParameters(
            antenna_r, antenna_d, event->Nu_Interaction[interaction_idx], interaction_idx,
            i, j, k, ray_idx, ray_output,
            &n_trg_pokey, &n_trg_slappy, &Pol_vector_src,
            &launch_vector, &receive_vector,
            icemodel, settings
        );

        double viewangle = antenna_r->view_ang[interaction_idx][ray_idx];
        double mag = antenna_r->Mag[interaction_idx][ray_idx];
        double fresnel = antenna_r->Fresnel[interaction_idx][ray_idx];
        double antenna_theta = antenna_r->theta_rec[interaction_idx][ray_idx] * 180 / PI;
        double antenna_phi = antenna_r->phi_rec[interaction_idx][ray_idx] * 180 / PI;
        Vector Pol_vector = antenna_r->Pol_vector[interaction_idx][ray_idx];

        double vmmhz1m_tmp = 0;
        double vmmhz1m_sum = 0;
        double vmmhz1m_em  = 0;

        // old freq domain signal mode (AVZ model)
        if (settings->SIMULATION_MODE == 0) {

            // initially give raysol has actual signal
            antenna_r->SignalExt[interaction_idx][ray_idx] = 1;

            double vmmhz_filter[(int)(detector->GetFreqBin())];

            double Pol_factor = 0;

            // In frequency space, modify the signal from the cascade according
            //   to the ice attenuation, the antennas effective height,
            //   the antenna's Filter, the Preamplifiers, and the FOAM
            for (int l = 0; l < detector->GetFreqBin(); l++) { // for detector freq bin numbers

                if (event->IsCalpulser > 0) {
                    vmmhz1m_tmp = event->Nu_Interaction[interaction_idx].vmmhz1m[l] *settings->CALPUL_AMP;
                }
                else {
                    vmmhz1m_tmp = event->Nu_Interaction[interaction_idx].vmmhz1m[l];
                    signal->TaperVmMHz(
                        viewangle,
                        event->Nu_Interaction[interaction_idx].d_theta_em[l],
                        event->Nu_Interaction[interaction_idx].d_theta_had[l],
                        event->Nu_Interaction[interaction_idx].emfrac,
                        event->Nu_Interaction[interaction_idx].hadfrac,
                        vmmhz1m_tmp, vmmhz1m_em);
                }

                // multiply all factors for traveling through ice
                if (settings->USE_ARA_ICEATTENU == 1 || settings->USE_ARA_ICEATTENU == 0) {
                    // now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
                    vmmhz1m_tmp = vmmhz1m_tmp / ray_output[0][ray_idx] *IceAttenFactor *mag * fresnel;
                }
                else if (settings->USE_ARA_ICEATTENU == 2) {

                    double IceAttenFactor = 1.;
                    double dx, dz, dl;
                    for (int steps = 1; steps < (int) RayStep[ray_idx][0].size(); steps++) {
                        dx = RayStep[ray_idx][0][steps - 1] - RayStep[ray_idx][0][steps];
                        dz = RayStep[ray_idx][1][steps - 1] - RayStep[ray_idx][1][steps];
                        dl = sqrt((dx *dx) + (dz *dz));

                        // Skipping attenuation calculation when the distance between two RaySteps is 0.
                        // Prevening adds -nan into the IceAttenFactor.
                        if (dl > 0)
                        {
                            // use ray midpoint for attenuation calculation
                            IceAttenFactor *= (exp(-dl / icemodel->GetFreqDepIceAttenuLength(-RayStep[ray_idx][1][steps], detector->GetFreq(l) / 1e9)) +
                                exp(-dl / icemodel->GetFreqDepIceAttenuLength(-RayStep[ray_idx][1][steps - 1], detector->GetFreq(l) / 1e9))
                            ) / 2.;  // 1e9 to convert to GHz
                        }
                    }

                    // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
                    vmmhz1m_tmp = vmmhz1m_tmp / ray_output[0][ray_idx] *IceAttenFactor *mag * fresnel;

                }

                vmmhz1m_sum += vmmhz1m_tmp;

                antenna_r->vmmhz[interaction_idx][ray_idx].push_back(vmmhz1m_tmp);

                double freq_tmp = detector->GetFreq(l);    // freq in Hz

                // Get ant gain with 2-D interpolation
                double heff = GaintoHeight(
                    detector->GetGain_1D_OutZero(
                        freq_tmp *1.E-6, antenna_theta,  antenna_phi,
                        antenna_d->type, j, k),
                    freq_tmp,
                    icemodel->GetN(*antenna_d),
                    detector->GetImpedance(freq_tmp*1.E-6, antenna_d->type, j, k));

                antenna_r->Heff[interaction_idx][ray_idx].push_back(heff);

                // apply pol factor, heff
                if (event->IsCalpulser == 1) {
                    Pol_vector = n_trg_slappy;
                }
                else if (event->IsCalpulser == 2) {
                    Pol_vector = n_trg_pokey;
                }
                else if (event->IsCalpulser == 3) {
                    Pol_vector = n_trg_slappy;
                }
                else if (event->IsCalpulser == 4) {
                    Pol_vector = n_trg_slappy + n_trg_pokey;
                }

                ApplyAntFactors(
                    heff, n_trg_pokey, n_trg_slappy, Pol_vector, antenna_d->type,
                    Pol_factor, vmmhz1m_tmp, antenna_theta, antenna_phi);
                ApplyFilter_OutZero(freq_tmp, detector, vmmhz1m_tmp);
                ApplyPreamp_OutZero(freq_tmp, detector, vmmhz1m_tmp);
                ApplyFOAM_OutZero(freq_tmp, detector, vmmhz1m_tmp);

                vmmhz_filter[l] = vmmhz1m_tmp;

            }   // end for freq bin

            antenna_r->Pol_factor[interaction_idx].push_back(Pol_factor);

            double volts_forfft[settings->NFOUR / 2];  // array for fft
            MakeArraysforFFT(settings, detector, i, vmmhz_filter, volts_forfft);

            // save freq domain array which is prepaired for realft
            for (int n = 0; n < settings->NFOUR / 2; n++) {
                antenna_r->Vfft[interaction_idx][ray_idx].push_back(volts_forfft[n]);
            }

            // now, after realft, volts_forfft is time domain signal at backend of antenna
            Tools::realft(volts_forfft, -1, settings->NFOUR / 2);

            antenna_r->PeakV[interaction_idx].push_back(FindPeak(volts_forfft, settings->NFOUR / 2));

            Tools::NormalTimeOrdering(settings->NFOUR / 2, volts_forfft);

            // Store the signal from this interaction and this ray to this antenna
            for (int n = 0; n < settings->NFOUR / 2; n++) {
                if (settings->TRIG_ANALYSIS_MODE != 2) { // not pure noise mode (we need signal)
                    antenna_r->V[interaction_idx][ray_idx].push_back(volts_forfft[n]);
                }
                else if (settings->TRIG_ANALYSIS_MODE == 2) { // pure noise mode (set signal to 0)
                    antenna_r->V[interaction_idx][ray_idx].push_back(0.);
                }
            }

        }   // if SIMULATION_MODE = 0

        else if (settings->SIMULATION_MODE == 1) { // Time domain simulation mode

            // if event is not calpulser
            if (event->IsCalpulser == 0) {

                if (settings->EVENT_TYPE == 0) {
                    // Get the signal from Neutrino events

                    // see if integrated shower profile LQ is non-zero and near the cone viewangle
                    static const int outbin = 64;
                    if (event->Nu_Interaction[interaction_idx].LQ > 0 && (fabs(viewangle - signal->CHANGLE_ICE) <= settings->OFFCONE_LIMIT *RADDEG)) {

                        double atten_factor = 0.;
                        if (settings->USE_ARA_ICEATTENU == 1 || settings->USE_ARA_ICEATTENU == 0) {
                            atten_factor = 1. / ray_output[0][ray_idx] *IceAttenFactor *mag * fresnel;
                        }
                        else if (settings->USE_ARA_ICEATTENU == 2) {
                            atten_factor = 1. / ray_output[0][ray_idx] *mag * fresnel;  //apply freq dependent IceAttenFactor later
                        }

                        // signal before the antenna (get signal at 1m and apply atten factor)
                        double Tarray[outbin];
                        double Earray[outbin];
                        signal->GetVm_FarField_Tarray(
                            event, settings, viewangle, atten_factor, outbin, Tarray, Earray,
                            antenna_r->skip_bins[ray_idx], interaction_idx);

                        double dT_forfft = Tarray[1] - Tarray[0];  // step in ns

                        // Determine the number of bins that will be used to build the waveform that will be fourier transformed
                        DetermineWFBins(antenna_r, interaction_idx, ray_idx, dT_forfft, settings);

                        // Convert the time array so it works with signal calculator
                        vector< double > Tarray_vector;
                        for (int bin=0; bin<outbin; bin++){
                            Tarray_vector.push_back( Tarray[bin] );
                        }
                        vector< double > Earray_vector;
                        for (int bin=0; bin<outbin; bin++){
                            Earray_vector.push_back( Earray[bin] );
                        }

                        PropagateSignal(
                            dT_forfft, outbin, Tarray_vector, Earray_vector, T_forint,
                            interaction_idx, ray_idx, ray_output, launch_vector, time_diff_birefringence,
                            Pol_vector_src, Pol_vector, n_trg_slappy, n_trg_pokey,
                            antenna_r, antenna_d,
                            gain_ch_no, j, k, birefringence, detector, event, icemodel, settings);

                    }
                    else {
                        // no signal generating

                        // initially give raysol has actual signal
                        antenna_r->SignalExt[interaction_idx][ray_idx] = 0;

                        // if no signal, push_back 0 values (otherwise the value inside will remain as old value)
                        for (int n = 0; n < settings->NFOUR / 2; n++) {
                            if (n < outbin) {
                                antenna_r->Vm_zoom[interaction_idx][ray_idx].push_back(0.);
                                antenna_r->Vm_zoom_T[interaction_idx][ray_idx].push_back(n);
                            }
                            antenna_r->V[interaction_idx][ray_idx].push_back(0.);
                        }

                        antenna_r->Nnew[ray_idx] = settings->NFOUR / 2;

                        antenna_r->PeakV[interaction_idx].push_back(0.);

                    }

                } // neutrino events
                else if (settings->EVENT_TYPE == 10) {
                    // Get the signal from Arbitrary Events

                    int waveform_bin = (int) signal->ArbitraryWaveform_V.size();
                    double dT_forfft = signal->ArbitraryWaveform_T[1] - signal->ArbitraryWaveform_T[0]; // step in ns

                    // Determine the number of bins that will be used to build the waveform that will be fourier transformed
                    DetermineWFBins(antenna_r, interaction_idx, ray_idx, dT_forfft, settings);

                    // Convert time array to vector double so it'll work with the signal calculator
                    vector< double > ArbitraryWaveform_T_vector;
                    for (int bin=0; bin < waveform_bin; bin++){
                        ArbitraryWaveform_T_vector.push_back(signal->ArbitraryWaveform_T[bin]);
                    }
                    vector< double > ArbitraryWaveform_V_vector;
                    for (int bin=0; bin < waveform_bin; bin++){
                        ArbitraryWaveform_V_vector.push_back(signal->ArbitraryWaveform_V[bin]);
                    }

                    PropagateSignal(
                        dT_forfft, waveform_bin, ArbitraryWaveform_T_vector, ArbitraryWaveform_V_vector, T_forint,
                        interaction_idx, ray_idx, ray_output, launch_vector, time_diff_birefringence,
                        Pol_vector_src, Pol_vector, n_trg_slappy, n_trg_pokey,
                        antenna_r, antenna_d,
                        gain_ch_no, j, k, birefringence, detector, event, icemodel, settings);

                }   // Arbitrary Events
                else if (settings->EVENT_TYPE == 11) {
                    // Get the signal from simple pulser events

                    int waveform_bin = (int) signal->PulserWaveform_V.size();
                    double dT_forfft = signal->PulserWaveform_T[1] - signal->PulserWaveform_T[0]; // step in ns

                    // Determine the number of bins that will be used to build the waveform that will be fourier transformed
                    DetermineWFBins(antenna_r, interaction_idx, ray_idx, dT_forfft, settings);

                    // Define polarization at the source (using launch_vector (a unit vector))
                    double psi = TMath::DegToRad()*settings->CLOCK_ANGLE;
                    double theta = acos(launch_vector[2]);
                    double phi = atan2(launch_vector[1],launch_vector[0]);

                    //Justin's method
                    double newPol_vectorX = -cos(psi)*cos(theta)*cos(phi) + sin(psi)*sin(phi);
                    double newPol_vectorY = -cos(psi)*cos(theta)*sin(phi) - sin(psi)*cos(phi);
                    double newPol_vectorZ = cos(psi)*sin(theta);
                    Vector Pol_vector = Vector(newPol_vectorX, newPol_vectorY, newPol_vectorZ);

                    icemodel->GetFresnel(
                        ray_output[1][ray_idx], // launch_angle
                        ray_output[2][ray_idx], // rec_angle
                        ray_output[3][ray_idx], // reflect_angle
                        event->Nu_Interaction[interaction_idx].posnu,
                        launch_vector,
                        receive_vector,
                        settings,
                        fresnel,
                        Pol_vector);    // input src Pol and return Pol at trg
                    icemodel->GetMag(
                        mag,
                        ray_output[0][ray_idx], // ray path length
                        ray_output[1][ray_idx], // zenith angle of ray at launch
                        ray_output[2][ray_idx], // zenith angle of ray upon receipt
                        ray_idx,
                        event->Nu_Interaction[interaction_idx].posnu, // Neutrino
                        *antenna_d, // Antenna
                        -0.01, // 1cm antenna shift, inspired from NuRadioMC
                        icemodel, settings
                    );

                    PropagateSignal(
                        dT_forfft, waveform_bin, signal->PulserWaveform_T, signal->PulserWaveform_V, T_forint,
                        interaction_idx, ray_idx, ray_output, launch_vector, time_diff_birefringence,
                        Pol_vector_src, Pol_vector, n_trg_slappy, n_trg_pokey,
                        antenna_r, antenna_d,
                        gain_ch_no, j, k, birefringence, detector, event, icemodel, settings);

                } // Simple Pulser Simulation

                // PVA pulser simulation.  Separate from previous pulser event type to avoid breaking things.
                else if (settings->EVENT_TYPE == 12) {
                    // Get the signal from PVA Pulser events

                    // Import Voltage versus time fed into antenna (via data from Alisa in IDL2_InputVoltageVersusTime.txt)

                    int waveform_bin = (int) signal->InputVoltage_V.size();
                    double dT_forfft = signal->InputVoltage_T[1] - signal->InputVoltage_T[0];    // step in ns

                    // Determine the number of bins that will be used to build the waveform that will be fourier transformed
                    DetermineWFBins(antenna_r, interaction_idx, ray_idx, dT_forfft, settings);

                    //Defining polarization at the source (using launch_vector (a unit vector))
                    // double psi = TMath::DegToRad()*settings->CLOCK_ANGLE;
                    double psi = 0;  // In the absence of cross-pol, the polarization angle is nominally zero for
                                     //   an antenna in the ice that is azimuthally symmetric.
                                     // Cross pol psi is updated after emission from a Tx
                    double theta = acos(launch_vector[2]);
                    double phi = atan2(launch_vector[1],launch_vector[0]);

                    //Justin's method
                    double newPol_vectorX = -cos(psi)*cos(theta)*cos(phi) + sin(psi)*sin(phi);
                    double newPol_vectorY = -cos(psi)*cos(theta)*sin(phi) - sin(psi)*cos(phi);
                    double newPol_vectorZ = cos(psi)*sin(theta);
                    Vector Pol_vector = Vector(newPol_vectorX, newPol_vectorY, newPol_vectorZ);

                    //Apply Fresnel factors for magnification and 1/r dependence
                    icemodel->GetFresnel(
                        ray_output[1][ray_idx], // launch_angle
                        ray_output[2][ray_idx], // rec_angle
                        ray_output[3][ray_idx], // reflect_angle
                        event->Nu_Interaction[interaction_idx].posnu,
                        launch_vector,
                        receive_vector,
                        settings,
                        fresnel,
                        Pol_vector);    // input src Pol and return Pol at trg

                    PropagateSignal(
                        dT_forfft, waveform_bin, signal->InputVoltage_T, signal->InputVoltage_V, T_forint,
                        interaction_idx, ray_idx, ray_output, launch_vector, time_diff_birefringence,
                        Pol_vector_src, Pol_vector, n_trg_slappy, n_trg_pokey,
                        antenna_r, antenna_d,
                        gain_ch_no, j, k, birefringence, detector, event, icemodel, settings);


                } // PVA Pulser Events

            }   // if not calpulser event

            // if calpulser event
            else if (event->IsCalpulser > 0) {
                // Get the signal from Calpulser events

                int CP_bin = (int) detector->CalPulserWF_ns.size();

                double dT_forfft = detector->CalPulserWF_ns[1] - detector->CalPulserWF_ns[0];  // step in ns

                // Determine the number of bins that will be used to build the waveform that will be fourier transformed
                DetermineWFBins(antenna_r, interaction_idx, ray_idx, dT_forfft, settings);

                PropagateSignal(
                    dT_forfft, CP_bin, detector->CalPulserWF_ns, detector->CalPulserWF_V, T_forint,
                    interaction_idx, ray_idx, ray_output, launch_vector, time_diff_birefringence,
                    Pol_vector_src, Pol_vector, n_trg_slappy, n_trg_pokey,
                    antenna_r, antenna_d,
                    gain_ch_no, j, k, birefringence, detector, event, icemodel, settings);

            }   // Calpulser events

        }   // if SIMULATION_MODE = 1

    }   // if not debug mode

}

void Report::GetRayParameters(
    Antenna_r *antenna_r, Antenna *antenna_d,
    Interaction interaction, int interaction_idx,
    int i, int j, int k,
    int ray_idx, vector<vector< double >> ray_output,
    Vector *n_trg_pokey, Vector *n_trg_slappy, Vector *Pol_vector_src,
    Position *launch_vector, Position *receive_vector,
    IceModel *icemodel, Settings *settings1
){
    // For a given ray path connecting a cascade to an antenna, calculate and
    //   save parameters including the antenna's viewing angle of the ray,
    //   attenuation length, magnification factors, etc.
    // Also initializes the objects that store electric field and voltage
    //   signals from the cascade by resizing those objects to the number of rays.

    // set viewangle, launch_vector, receive vectors
    double viewangle = ray_output[1][ray_idx];
    GetParameters(
        interaction.posnu, // posnu
        *antenna_d, // trg antenna
        interaction.nnu, // nnu
        viewangle, // inputs launch_angle, returns viewangle
        ray_output[2][ray_idx], // receive_angle
        *launch_vector, *receive_vector,
        *n_trg_slappy, *n_trg_pokey);

    // store information to report
    antenna_r->view_ang[interaction_idx].push_back(viewangle);
    antenna_r->launch_ang[interaction_idx].push_back(ray_output[1][ray_idx]);
    antenna_r->rec_ang[interaction_idx].push_back(ray_output[2][ray_idx]);
    antenna_r->Dist[interaction_idx].push_back(ray_output[0][ray_idx]);
    antenna_r->L_att[interaction_idx].push_back(icemodel->EffectiveAttenuationLength(settings1, interaction.posnu, 0));
    antenna_r->reflect_ang[interaction_idx].push_back(ray_output[3][ray_idx]);
    antenna_r->vmmhz[interaction_idx].resize(ray_idx + 1);
    antenna_r->Heff[interaction_idx].resize(ray_idx + 1);
    antenna_r->Heff_copol[interaction_idx].resize(ray_idx + 1);
    antenna_r->Heff_crosspol[interaction_idx].resize(ray_idx + 1);
    antenna_r->Vm_zoom[interaction_idx].resize(ray_idx + 1);
    antenna_r->Vm_zoom_T[interaction_idx].resize(ray_idx + 1);
    antenna_r->Vfft[interaction_idx].resize(ray_idx + 1);
    antenna_r->Vfft_noise[interaction_idx].resize(ray_idx + 1);
    antenna_r->V[interaction_idx].resize(ray_idx + 1);
    antenna_r->SignalExt[interaction_idx].resize(ray_idx + 1);

    // calculate the polarization vector at the source
    Position Pol_vector = GetPolarization(interaction.nnu, *launch_vector); // polarization vector at the source
    *Pol_vector_src = Pol_vector; //store the src Pol

    double fresnel = 0;
    icemodel->GetFresnel(
        ray_output[1][ray_idx], // launch_angle
        ray_output[2][ray_idx], // rec_angle
        ray_output[3][ray_idx], // reflect_angle
        interaction.posnu,
        *launch_vector,
        *receive_vector,
        settings1,
        fresnel,
        Pol_vector); // input src Pol and return Pol at trg
    double mag = 0;
    icemodel->GetMag(
        mag,
        ray_output[0][ray_idx], // ray path length
        ray_output[1][ray_idx], // zenith angle of ray at launch
        ray_output[2][ray_idx], // zenith angle of ray upon receipt
        ray_idx,
        interaction.posnu, // Neutrino
        *antenna_d, // Antenna
        -0.01, // 1cm antenna shift, inspired from NuRadioMC
        icemodel, settings1
    );

    if (ray_output[3][ray_idx] < PI / 2.) {
        // when not reflected at the surface, angle = 100
        antenna_r->reflection[interaction_idx].push_back(1); // say this is reflected ray
    }
    else {
        antenna_r->reflection[interaction_idx].push_back(0); // say this is not reflected ray
    }

    antenna_r->Pol_vector[interaction_idx].push_back(Pol_vector);    // this Pol_vector is for the target antenna
    antenna_r->Mag[interaction_idx].push_back(mag);  // magnification factor
    antenna_r->Fresnel[interaction_idx].push_back(fresnel);  // Fresnel factor

    // get the arrival angle at the antenna, and store the relevant polarization factors
    double antenna_theta = 0;
    double antenna_phi = 0;
    GetAngleAnt(*receive_vector, *antenna_d, antenna_theta, antenna_phi);  // get theta, phi for signal ray arrived at antenna
    double launch_theta = 0;
    double launch_phi = 0;
    GetAngleLaunch(*launch_vector, launch_theta, launch_phi);
    Vector thetaHat = Vector(cos(antenna_theta *(PI / 180)) *cos(antenna_phi *(PI / 180)),
        cos(antenna_theta *(PI / 180)) *sin(antenna_phi *(PI / 180)),
        -sin(antenna_theta *(PI / 180)));
    Vector phiHat = Vector(-sin(antenna_phi *(PI / 180)),
        cos(antenna_phi *(PI / 180)),
        0);
    antenna_r->Pol_factorH[interaction_idx].push_back(abs(phiHat *Pol_vector));
    antenna_r->Pol_factorV[interaction_idx].push_back(abs(thetaHat *Pol_vector));
    antenna_r->phi_rec[interaction_idx].push_back(antenna_phi *(PI / 180));
    antenna_r->theta_rec[interaction_idx].push_back(antenna_theta *(PI / 180));
    antenna_r->phi_launch[interaction_idx].push_back(launch_phi *(PI / 180));
    antenna_r->theta_launch[interaction_idx].push_back(launch_theta *(PI / 180));

}

void Report::DetermineWFBins(
    Antenna_r *antenna, int interaction_idx, int ray_idx, double dT, Settings *settings1
){
    // Chooses the number of bins for the waveforms that will be fourier transformed
    // The Fourier transform performs more quickly if the waveform has a length
    //   that is a power of two.
    antenna->SignalExt[interaction_idx][ray_idx] = 1;
    int Ntmp = settings1->TIMESTEP *1.e9 / dT;
    antenna->Nnew[ray_idx] = 1;
    while (Ntmp > 1) {
        Ntmp = Ntmp / 2;
        antenna->Nnew[ray_idx] = antenna->Nnew[ray_idx] *2;
    }
    antenna->Nnew[ray_idx] = antenna->Nnew[ray_idx] *settings1->NFOUR / 2;
}


void Report::PropagateSignal(
    double dT_forfft, int efield_length, vector< double > efield_time, vector< double > efield, double *T_forint,
    int interaction_idx, int ray_idx, vector<vector< double > > ray_output, Position launch_vector, double time_diff_birefringence,
    Vector Pol_vector_src, Vector Pol_vector, Vector n_trg_slappy, Vector n_trg_pokey,
    Antenna_r *antenna_r, Antenna *antenna_d, int gain_ch_no, int j, int k,
    Birefringence *birefringence, Detector *detector, Event *event, IceModel *icemodel, Settings *settings
){
    // For the provided `antenna_r`, interaction, and ray:
    //   take the EM signal at the source and calculate how it's read out
    //   as a voltage waveform by the antenna by applying all attenuation,
    //   birefringence, antenna, and other factors.

    // now we have to make NFOUR/2 number of bins with random init time
    // as a test, make first as it is and zero pad
    double V_forfft[antenna_r->Nnew[ray_idx]];
    double T_forfft[antenna_r->Nnew[ray_idx]];

    double antenna_theta = antenna_r->theta_rec[interaction_idx][ray_idx] * 180 / PI;
    double antenna_phi = antenna_r->phi_rec[interaction_idx][ray_idx] * 180 / PI;

    double Pol_factor;  // polarization factor

    int max_bire_ray_cnt = settings->BIREFRINGENCE + 1; // rays in birefringence per ray solution

    double V_forfft_bire[max_bire_ray_cnt][antenna_r->Nnew[ray_idx]]; // for the waveforms of the rays in birefringence

    for ( int bire_ray_cnt = 0; bire_ray_cnt < max_bire_ray_cnt; bire_ray_cnt++ ) {

        max_bire_ray_cnt = birefringence->Reflected_ray_remove_bire(
            ray_output[3][ray_idx], max_bire_ray_cnt); //change to 1 if the ray solution is reflected

        for (int n = 0; n < antenna_r->Nnew[ray_idx]; n++) {

            if (n < efield_length) {
                antenna_r->Vm_zoom[interaction_idx][ray_idx].push_back(efield[n]);
                antenna_r->Vm_zoom_T[interaction_idx][ray_idx].push_back(efield_time[n]);
            }

            // make efield_time, efield located at the center of Nnew array
            T_forfft[n] = efield_time[efield_length / 2] - (dT_forfft *(double)(antenna_r->Nnew[ray_idx] / 2 - n));

            if ((n >= antenna_r->Nnew[ray_idx] / 2 - efield_length / 2) &&
                (n <  antenna_r->Nnew[ray_idx] / 2 + efield_length / 2)   ) {
                V_forfft[n] = efield[n - (antenna_r->Nnew[ray_idx] / 2 - efield_length / 2)];
            }
            else{
                V_forfft[n] = 0.;
            }

        }

        // just get peak from the array
        antenna_r->PeakV[interaction_idx].push_back(FindPeak(efield, (int)efield_length));

        int T_shift_bire = int(time_diff_birefringence/dT_forfft); //time shift for birefringence
        double split_factor_bire = birefringence->Power_split_factor(
            Pol_vector_src, bire_ray_cnt, ray_output[3][ray_idx], settings
        ); //split power factor for birefringence
        birefringence->Time_shift_and_power_split(
            V_forfft, antenna_r->Nnew[ray_idx],
            T_shift_bire, split_factor_bire, bire_ray_cnt, max_bire_ray_cnt, settings
        ); // apply time differences and power split

        // this forward fft volts_forfft is now in unit of V at each freq we can just apply each bin's gain factor to each freq bins
        // without any phase consideration,
        // apply same factor to both real, img parts

        // get spectrum with zero padded WF
        Tools::realft(V_forfft, 1, antenna_r->Nnew[ray_idx]);

        double dF_Nnew = 1. / ((double)(antenna_r->Nnew[ray_idx]) *(dT_forfft) *1.e-9);    // in Hz

        double freq_tmp = dF_Nnew *((double) antenna_r->Nnew[ray_idx] / 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

        double freq_lastbin = freq_tmp;

        //For birefringence, modify the polarization at the antennas
        birefringence->Principal_axes_polarization(Pol_vector, bire_ray_cnt, max_bire_ray_cnt, settings);

        // Get ant gain with 2-D interpolation
        double heff_lastbin = GaintoHeight(
            detector->GetGain_1D_OutZero(
                freq_tmp *1.E-6, antenna_theta, antenna_phi,
                antenna_d->type, j, k),
            freq_tmp,
            icemodel->GetN(*antenna_d),
            detector->GetImpedance(freq_tmp*1.E-6, antenna_d->type, k));

        // Tx effective height for last bin.  Currently locked to standard ARA Vpol and HPol antennas.
        // Need to add selection mode.
        // Only used in EVENT_MODE == 12
        double Tx_theta = 0;
        double Tx_phi = 0;
        double heff_Tx_lastbin = 0;
        if (settings->EVENT_TYPE == 11 || settings->EVENT_TYPE == 12){ // Pulser simulations

            //Defining polarization at the source (using launch_vector (a unit vector))
            double theta = acos(launch_vector[2]);
            double phi = atan2(launch_vector[1],launch_vector[0]);

            double Tx_theta = theta*180/PI;
            double Tx_phi = phi*180/PI;
            double heff_Tx_lastbin = GaintoHeight(
                detector->GetGain_1D_OutZero(freq_tmp *1.E-6, Tx_theta, Tx_phi, 0, 0, 0, true),
                freq_tmp,
                icemodel->GetN(*antenna_d),
                detector->GetImpedance(freq_tmp*1.E-6, 0, 0, true));
        }
        else if (event->IsCalpulser > 0){ // Calibration Pulser simulations
            Tx_theta = ray_output[1][ray_idx] *DEGRAD;    // from 0 to 180
            heff_Tx_lastbin = GaintoHeight(
                detector->GetGain_1D_OutZero(freq_tmp *1.E-6, Tx_theta, antenna_phi, antenna_d->type, j, k),
                freq_tmp, icemodel->GetN(event->Nu_Interaction[interaction_idx].posnu));

            if (event->IsCalpulser == 1) {
                Pol_vector = n_trg_slappy;
            }
            else if (event->IsCalpulser == 2) {
                Pol_vector = n_trg_pokey;
            }
            else if (event->IsCalpulser == 3) {
                Pol_vector = n_trg_slappy;
            }
            else if (event->IsCalpulser == 4) {
                Pol_vector = n_trg_slappy + n_trg_pokey;
            }

        }

        // Apply frequency-dependent scaling factors to the waveform like
        //   ice attenuation, antenna effective height, and electronics gain
        for (int n = 0; n < antenna_r->Nnew[ray_idx] / 2; n++) { // Loop over frequency bins

            freq_tmp = dF_Nnew *((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

            double heff = GaintoHeight(
                detector->GetGain_1D_OutZero(freq_tmp *1.E-6, antenna_theta,  antenna_phi, antenna_d->type, j, k),
                freq_tmp,
                icemodel->GetN(*antenna_d),
                detector->GetImpedance(freq_tmp*1.E-6, antenna_d->type, k));
            antenna_r->Heff[interaction_idx][ray_idx].push_back(heff);

            //apply freq dependent attenuation model if in neutrino mode
            if (settings->EVENT_TYPE == 0 && settings->USE_ARA_ICEATTENU == 2) {

                double IceAttenFactor = 1.;
                double dx, dz, dl;
                for (int steps = 1; steps < (int) RayStep[ray_idx][0].size(); steps++) {
                        dx = RayStep[ray_idx][0][steps - 1] - RayStep[ray_idx][0][steps];
                        dz = RayStep[ray_idx][1][steps - 1] - RayStep[ray_idx][1][steps];
                        dl = sqrt((dx *dx) + (dz *dz));

                        // Skipping attenuation calculation when the distance between two RaySteps is 0.
                        // Prevening adds -nan into the IceAttenFactor.
                        if (dl > 0) {
                            // use ray midpoint for attenuation calculation
                            IceAttenFactor *= (
                                exp(-dl / icemodel->GetFreqDepIceAttenuLength(-RayStep[ray_idx][1][steps],     freq_tmp *1.E-9)) +
                                exp(-dl / icemodel->GetFreqDepIceAttenuLength(-RayStep[ray_idx][1][steps - 1], freq_tmp *1.E-9))
                            ) / 2.;  // 1e9 for conversion to GHz
                        }

                }

                V_forfft[2 *n] *= IceAttenFactor;   // apply IceAttenFactor to the real part of fft
                V_forfft[2 *n + 1] *= IceAttenFactor;   // apply IceAttenFactor to the imag part of fft

            }

            // Apply transmitter effective height if in PVA Pulser mode
            if (settings->EVENT_TYPE == 12){
                double heff_Tx = GaintoHeight(
                    detector->GetGain_1D_OutZero(freq_tmp *1.E-6, Tx_theta, Tx_phi, 0, 0, 0, true),
                    freq_tmp,
                    icemodel->GetN(*antenna_d),
                    detector->GetImpedance(freq_tmp*1.E-6, 0, 0, true));
                if (n > 0) {
                    ApplyAntFactors_Tdomain(
                        detector->GetAntPhase_1D(freq_tmp *1.e-6, Tx_theta, Tx_phi, 0),
                        heff_Tx, Pol_vector, 0, Pol_factor, V_forfft[2 *n], V_forfft[2 *n + 1],
                        settings, Tx_theta, Tx_phi, freq_tmp, detector->GetImpedance(freq_tmp*1.E-6, 0, 0, true), true);
                }
                else {
                    ApplyAntFactors_Tdomain_FirstTwo(
                        heff_Tx, heff_Tx_lastbin, Pol_vector, 0, Pol_factor,
                        V_forfft[2 *n], V_forfft[2 *n + 1], settings, Tx_theta, Tx_phi, freq_tmp);
                }
            }
            else if (event->IsCalpulser > 0){
                // apply ant factors (transmitter ant)
                heff = GaintoHeight(
                    detector->GetGain_1D_OutZero(
                        freq_tmp *1.E-6,   // to MHz
                        Tx_theta, antenna_phi, antenna_d->type, j, k),
                    freq_tmp, icemodel->GetN(event->Nu_Interaction[interaction_idx].posnu));

                if (n > 0) {
                    ApplyAntFactors_Tdomain(
                        detector->GetAntPhase_1D(freq_tmp *1.e-6, Tx_theta, antenna_phi, antenna_d->type),
                        heff, Pol_vector, antenna_d->type, Pol_factor,
                        V_forfft[2 *n], V_forfft[2 *n + 1], settings, true, antenna_theta, antenna_phi, freq_tmp);
                }
                else {
                    ApplyAntFactors_Tdomain_FirstTwo(
                        heff, heff_Tx_lastbin, Pol_vector, antenna_d->type, Pol_factor,
                        V_forfft[2 *n], V_forfft[2 *n + 1], settings, antenna_theta, antenna_phi, freq_tmp);
                }
            }

            // apply ant factors
            if (n > 0) {
                ApplyAntFactors_Tdomain(
                    detector->GetAntPhase_1D(freq_tmp *1.e-6, antenna_theta, antenna_phi, antenna_d->type),
                    heff, Pol_vector, antenna_d->type, Pol_factor,
                    V_forfft[2 *n], V_forfft[2 *n + 1], settings, antenna_theta, antenna_phi, freq_tmp);
            }
            else {
                ApplyAntFactors_Tdomain_FirstTwo(
                    heff, heff_lastbin, Pol_vector, antenna_d->type, Pol_factor,
                    V_forfft[2 *n], V_forfft[2 *n + 1], settings, antenna_theta, antenna_phi, freq_tmp);
            }

            // apply entire elect chain gain, phase
            if (n > 0) {
                ApplyElect_Tdomain(
                    freq_tmp *1.e-6, detector,
                    V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings);
            }
            else {
                ApplyElect_Tdomain_FirstTwo(
                    freq_tmp *1.e-6, freq_lastbin *1.e-6, detector,
                    V_forfft[2 *n], V_forfft[2 *n + 1], gain_ch_no, settings);
            }

        }   // end for freq bin

        // now get time domain waveform back by inv fft
        Tools::realft(V_forfft, -1, antenna_r->Nnew[ray_idx]);

        birefringence->Store_V_forfft_for_interference(
            V_forfft, V_forfft_bire[bire_ray_cnt], antenna_r->Nnew[ray_idx]
        ); //Store waveforms from birefringence for interference

    } // end for bire_ray_cnt

    birefringence->Two_rays_interference(
        V_forfft, V_forfft_bire[0], V_forfft_bire[1], antenna_r->Nnew[ray_idx], max_bire_ray_cnt, settings
    ); //Apply interference of two rays from birefringence

    antenna_r->Pol_factor[interaction_idx].push_back(Pol_factor);

    // do linear interpolation
    // changed to sinc interpolation Dec 2020 by BAC
    double volts_forint[settings->NFOUR / 2];  // array for interpolation
    Tools::SincInterpolation(antenna_r->Nnew[ray_idx], T_forfft, V_forfft, settings->NFOUR / 2, T_forint, volts_forint);

    for (int n = 0; n < settings->NFOUR / 2; n++) {
        if (settings->TRIG_ANALYSIS_MODE != 2) {
            // not pure noise mode (we need signal)
            antenna_r->V[interaction_idx][ray_idx].push_back(
                volts_forint[n] *2. / (double)(antenna_r->Nnew[ray_idx])
            );  // 2/N for inverse FFT normalization factor
        }
        else if (settings->TRIG_ANALYSIS_MODE == 2) {
            // pure noise mode (set signal to 0)
            antenna_r->V[interaction_idx][ray_idx].push_back(0.);
        }
    }

}

void Report::triggerCheck_ScanMode0(
    int trig_search_init, int max_total_bin, int trig_window_bin,
    int n_scanned_channels, int station_index,
    Detector *detector, Event *event, Settings *settings1, Trigger *trigger
){
    // Check the tunnel diode for each antenna for sufficient signal to trigger
    //   then perform station-wide trigger check. Transforms waveforms into
    //   the data-like V_mimic objects
    // Original AraSim trigger checking mode

    int this_bin = trig_search_init;

    // avoid really long trig_window_bin case (change trig_window to check upto max_total_bin)
    if (max_total_bin - trig_window_bin <= this_bin) {
        trig_window_bin = max_total_bin - this_bin - 1;
    }

    int trig_mode = settings1->TRIG_MODE;
    // global trigger mode
    // 0 for orginal N_TRIG out of all channels
    // 1 for stations 1-5, N_TRIG_V out of Vpol channels or N_TRIG_H out of Hpol channels

    // mode select for non-trigger passed chs' V_mimic
    int V_mimic_mode = settings1->V_MIMIC_MODE;
    // 0 for orginal style (save the middle of trig_window)
    // 1 for saving waveform starting from last_trig_bin
    // 2 for saving waveform where last_trig_bin located on the middle of the waveform

    // Initialize staiton objects
    Station_r *station_r = &stations[station_index];
    ARA_station *station_d = &detector->stations[station_index];

    // Loop over global time bins and identify the first, if any, time bin
    //    with sufficient signal to trigger in the appropriate amount of antennas
    while (this_bin < max_total_bin - trig_window_bin) { // loop over waveform bins

        int N_pass = 0;
        int N_pass_V = 0;
        int N_pass_H = 0;
        int last_trig_bin = 0; // stores last trigger passed bin number
        Passed_chs.clear();

        int antenna_counter = 0;
        while (antenna_counter < n_scanned_channels) { // loop over antennas

            int string_i = detector->getStringfromArbAntID(station_index, antenna_counter);
            int antenna_i = detector->getAntennafromArbAntID(station_index, antenna_counter);
            int channel_num = detector->GetChannelfromStringAntenna(station_index, string_i, antenna_i, settings1);

            // Channel numbering is different for DETECTOR=(1,2,3) than for
            // DETECTOR = 4 and 5 in GetChannelfromStringAntenna(), it needs that shift
            if (!(settings1->DETECTOR==4 || settings1->DETECTOR==5)){
                channel_num = channel_num+1;
            }

            // Antenna Masking (masked_ant=0 means this antenna should be ignored from trigger)
            if( detector->GetTrigMasking(channel_num-1)==0){
                antenna_counter++;
                continue;
            }

            int offset = detector->GetTrigOffset(channel_num-1, settings1);

            // If Testbed simulation, check if we want to use BH chs only for trigger analysis
            if ((settings1->TRIG_ONLY_BH_ON == 1) && (settings1->DETECTOR == 3)) {

                // check if this channel is BH ch (DAQchan)
                if (station_d->strings[string_i].antennas[antenna_i].DAQchan == 0) {

                    int trig_bin = 0;
                    while (trig_bin < trig_window_bin) {

                        double diode_noise_RMS = trigger->GetAntNoise_diodeRMS(channel_num-1, settings1);
                        if (
                            trigger->Full_window[antenna_counter][this_bin + trig_bin] <
                            (   detector->GetThres(station_index, channel_num - 1, settings1)
                                * diode_noise_RMS
                                * detector->GetThresOffset(station_index, channel_num - 1, settings1))
                        ) {
                            // if this channel passed the trigger!

                            station_r->strings[string_i].antennas[antenna_i].Trig_Pass = this_bin + trig_bin;
                            N_pass++;
                            if (station_d->strings[string_i].antennas[antenna_i].type == 0) { // Vpol
                                N_pass_V++;
                            }
                            if (station_d->strings[string_i].antennas[antenna_i].type == 1) { // Hpol
                                N_pass_H++;
                            }
                            if (last_trig_bin < this_bin + trig_bin)
                                last_trig_bin = this_bin + trig_bin;   // added for fixed V_mimic
                            trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
                            Passed_chs.push_back(antenna_counter);

                        }

                        trig_bin++;
                    }
                }
            }

            // For non-Testbed simulations, check if we just want to use first/lower 8 chs' thres values
            else if ((settings1->TRIG_ONLY_LOW_CH_ON == 1) && (settings1->DETECTOR != 3)) {

                // reset channel numbers so that bottom antennas have ch 1-8
                channel_num = GetChannelNum8_LowAnt(string_i, antenna_i);

                if (antenna_i < 2) {
                    // only antenna 0, 1 which are bottom 2 antennas

                    // set channel_num as new value (antenna 0, 1 are only
                    // possible antennas for channel_num 1 - 8)

                    int trig_bin = 0;
                    while (trig_bin < trig_window_bin) {

                        double diode_noise_RMS = trigger->GetAntNoise_diodeRMS(channel_num-1, settings1);
                        // with threshold offset by chs
                        if (
                            trigger->Full_window[antenna_counter][this_bin + trig_bin] <
                            (   detector->GetThres(station_index, channel_num - 1, settings1)
                                * diode_noise_RMS
                                * detector->GetThresOffset(station_index, channel_num - 1, settings1)))
                        {
                            // if this channel passed the trigger!

                            station_r->strings[string_i].antennas[antenna_i].Trig_Pass = this_bin + trig_bin;
                            N_pass++;
                            if (station_d->strings[string_i].antennas[antenna_i].type == 0) {
                                // Vpol
                                N_pass_V++;
                            }
                            if (station_d->strings[string_i].antennas[antenna_i].type == 1) {
                                // Hpol
                                N_pass_H++;
                            }
                            if (last_trig_bin < this_bin + trig_bin)
                                last_trig_bin = this_bin + trig_bin;   // added for fixed V_mimic
                            trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
                            Passed_chs.push_back(antenna_counter);

                        }

                        trig_bin++;
                    }
                }
            }

            // other cases: use all possible chs for trigger analysis
            else {

                int trig_bin = 0;
                while (trig_bin < trig_window_bin) {

                    double diode_noise_RMS = trigger->GetAntNoise_diodeRMS(channel_num-1, settings1);
                    if( this_bin+offset+trig_bin >= settings1->DATA_BIN_SIZE )
                        break; //if trigger window hits wf end, cannot scan this channel further with this this_bin
                    if (
                        trigger->Full_window[antenna_counter][this_bin + trig_bin + offset] <
                        (   detector->GetThres(station_index, channel_num - 1, settings1)
                            * diode_noise_RMS
                            * detector->GetThresOffset(station_index, channel_num - 1, settings1))
                    ) {
                        // if this channel passed the trigger!

                        station_r->strings[string_i].antennas[antenna_i].Trig_Pass = this_bin + trig_bin + offset;
                        N_pass++;
                        if (station_d->strings[string_i].antennas[antenna_i].type == 0) { // Vpol
                            N_pass_V++;
                        }
                        if (station_d->strings[string_i].antennas[antenna_i].type == 1) { // Hpol
                            N_pass_H++;
                        }
                        if (last_trig_bin < this_bin + trig_bin)
                            last_trig_bin = this_bin + trig_bin + offset;   // added for fixed V_mimic
                        trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
                        Passed_chs.push_back(antenna_counter);

                    }

                    trig_bin++;
                }
            }

            antenna_counter++;   // if station not passed the trigger, just go to next channel

        }   // while antenna_counter < n_scanned_channels

        if (((trig_mode == 0) && (N_pass > settings1->N_TRIG - 1))  // trig_mode = 0 case!
            ||  // or
            ((trig_mode == 1) && ((N_pass_V > settings1->N_TRIG_V - 1) || (N_pass_H > settings1->N_TRIG_H - 1)))    // trig_mode = 1 case!
        ) {

            int check_ch = 0; // Indexer for the Passed_chs object
            station_r->Global_Pass = last_trig_bin; // Save the waveform bin when the station triggered and read out

            this_bin = max_total_bin; // also if we know this station is trigged, don't need to check rest of time window

            // Loop over all antennas on the station and save the data-like waveform, V_mimic
            for (int ch_loop = 0; ch_loop < n_scanned_channels; ch_loop++) {
                int string_i = detector->getStringfromArbAntID(station_index, ch_loop);
                int antenna_i = detector->getAntennafromArbAntID(station_index, ch_loop);

                int waveformLength = settings1->WAVEFORM_LENGTH;
                int waveformCenter = settings1->WAVEFORM_CENTER;

                station_r->strings[string_i].antennas[antenna_i].Likely_Sol[0] = -1;  // no likely init
                station_r->strings[string_i].antennas[antenna_i].Likely_Sol[1] = -1;  // no likely init

                // If this channel was one that passed, also determine which
                //   ray and interaction triggered the detector and print some info
                if (ch_loop == Passed_chs[check_ch] && check_ch < N_pass) {
                    // `check_ch < N_Pass` present to circumvent a bug in vector Passed_chs.clear()

                    // Determine which ray and interaciton triggered triggered
                    //   the station based on signal and trigger bins
                    station_r->strings[string_i].antennas[antenna_i].Find_Likely_Sol();
                    int likely_int = station_r->strings[string_i].antennas[antenna_i].Likely_Sol[0];
                    int likely_ray = station_r->strings[string_i].antennas[antenna_i].Likely_Sol[1];

                    if (settings1->TRIG_ONLY_LOW_CH_ON == 0) {

                        cout << endl
                             << "trigger passed at bin "
                             << station_r->strings[string_i].antennas[antenna_i].Trig_Pass
                             << "  passed ch : " << ch_loop
                             << " (" << station_d->strings[string_i].antennas[antenna_i].type
                             << "type) Direct dist btw posnu : "
                             << event->Nu_Interaction[0].posnu.Distance(station_d->strings[string_i].antennas[antenna_i])
                             << " noiseID : " << station_r->strings[string_i].antennas[antenna_i].noise_ID[0];
                        if (station_r->strings[string_i].antennas[antenna_i].Likely_Sol[0] != -1) {
                            cout << " ViewAngle : "
                                 << station_r->strings[string_i].antennas[antenna_i].view_ang[likely_int][likely_ray] *DEGRAD ;
                            cout << " LikelyTrigSignal : interaction "
                                 << station_r->strings[string_i].antennas[antenna_i].Likely_Sol[0];
                            cout << ", ray "
                                 << station_r->strings[string_i].antennas[antenna_i].Likely_Sol[1];
                        }

                    }
                    else if (settings1->TRIG_ONLY_LOW_CH_ON == 1) {

                        cout << endl
                             << "trigger passed at bin "
                             << station_r->strings[string_i].antennas[antenna_i].Trig_Pass
                             << "  passed ant: str[" << string_i << "].ant[" << antenna_i << "] ("
                             << station_d->strings[string_i].antennas[antenna_i].type
                             << "type) Direct dist btw posnu : "
                             << event->Nu_Interaction[0].posnu.Distance(station_d->strings[string_i].antennas[antenna_i])
                             << " noiseID : " << station_r->strings[string_i].antennas[antenna_i].noise_ID[0];
                        if (station_r->strings[string_i].antennas[antenna_i].Likely_Sol[0] != -1) {
                            cout << " ViewAngle : "
                                 << station_r->strings[string_i].antennas[antenna_i].view_ang[likely_int][likely_ray] *DEGRAD ;
                            cout << " LikelyTrigSignal : interaction "
                                 << station_r->strings[string_i].antennas[antenna_i].Likely_Sol[0];
                            cout << ", ray "
                                 << station_r->strings[string_i].antennas[antenna_i].Likely_Sol[1];
                        }

                    }

                    check_ch++; // Increase the Passed_chs index tracker

                    // now save the voltage waveform to V_mimic
                    for (int mimicbin = 0; mimicbin < waveformLength; mimicbin++) {

                        // new DAQ waveform writing mechanism test
                        if (V_mimic_mode == 0) {
                            // Global passed bin is the center of the window
                            const int thisBin = last_trig_bin + waveformCenter - waveformLength / 2 + mimicbin;
                            const int thisTimeBin = waveformCenter - waveformLength / 2 + mimicbin;

                            station_r->strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][thisBin]) *1.e3);   // save in mV
                            station_r->strings[string_i].antennas[antenna_i].time.push_back(thisBin);
                            station_r->strings[string_i].antennas[antenna_i].time_mimic.push_back(thisTimeBin *settings1->TIMESTEP *1.e9);    // save in ns
                        }
                        else if (V_mimic_mode == 1) {
                            // Global passed bin is the center of the window + delay to each chs from araGeom
                            const int thisBin = last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop]
                                                - detector->params.TestBed_BH_Mean_delay_bin) + waveformCenter - waveformLength / 2 + mimicbin;
                            const int thisTimeBin = -1*(detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin)
                                                    + waveformCenter - waveformLength / 2 + mimicbin;

                            station_r->strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][thisBin]) *1.e3);   // save in mV
                            station_r->strings[string_i].antennas[antenna_i].time.push_back(thisBin);
                            station_r->strings[string_i].antennas[antenna_i].time_mimic.push_back(thisTimeBin *settings1->TIMESTEP *1.e9);   // save in ns
                        }
                        else if (V_mimic_mode == 2) {
                            // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
                            const int thisBin = last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop]
                                                - detector->params.TestBed_BH_Mean_delay_bin
                                                + station_d->strings[string_i].antennas[antenna_i].manual_delay_bin)
                                                + waveformCenter - waveformLength / 2 + mimicbin;
                            const int thisTimeBin = -1*(detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin
                                                    + station_d->strings[string_i].antennas[antenna_i].manual_delay_bin)
                                                    + waveformCenter - waveformLength / 2 + mimicbin;

                            station_r->strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][thisBin]) *1.e3);    // save in mV
                            station_r->strings[string_i].antennas[antenna_i].time.push_back(thisBin);
                            station_r->strings[string_i].antennas[antenna_i].time_mimic.push_back(thisTimeBin *settings1->TIMESTEP *1.e9 + detector->params.TestBed_WFtime_offset_ns);    // save in ns
                        }
                        if (mimicbin == 0) {
                            ///< calculates time of center of each rays signal based on readout window time config
                            for (int interaction_idx=0; interaction_idx<station_r->strings[string_i].antennas[antenna_i].SignalBin.size(); interaction_idx++){
                                for (int m = 0; m < station_r->strings[string_i].antennas[antenna_i].SignalBin[interaction_idx].size(); m++) {
                                    double signal_center_offset = (double)(
                                        station_r->strings[string_i].antennas[antenna_i].SignalBin[interaction_idx][m] -
                                        station_r->strings[string_i].antennas[antenna_i].time[0]
                                    ) * settings1->TIMESTEP * 1.e9;
                                    double signal_center_time = signal_center_offset + station_r->strings[string_i].antennas[antenna_i].time_mimic[0];
                                    //! signal_center_offset: time offset between beginning of readout window and center of signal
                                    //! signal_center_time: time of center of signal based on readout window time config
                                    station_r->strings[string_i].antennas[antenna_i].SignalBinTime[interaction_idx].push_back(signal_center_time);
                                }
                            }
                        }
                    }

                    // set global_trig_bin values
                    if (V_mimic_mode == 0) {
                        // Global passed bin is the center of the window
                        station_r->strings[string_i].antennas[antenna_i].global_trig_bin = (waveformLength / 2 - waveformCenter);
                    }
                    else if (V_mimic_mode == 1) {
                        // Global passed bin is the center of the window + delay to each chs from araGeom
                        station_r->strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[ch_loop]
                                                                                             - detector->params.TestBed_BH_Mean_delay_bin) - waveformCenter + waveformLength / 2;
                    }
                    else if (V_mimic_mode == 2) {
                        // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
                        station_r->strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[ch_loop]
                                                                                             - detector->params.TestBed_BH_Mean_delay_bin
                                                                                             + station_d->strings[string_i].antennas[antenna_i].manual_delay_bin)
                                                                                             - waveformCenter + waveformLength / 2;
                    }
                }
                else {
                    station_r->strings[string_i].antennas[antenna_i].Trig_Pass = 0.;

                    // new DAQ waveform writing mechanism test
                    for (int mimicbin = 0; mimicbin < waveformLength; mimicbin++) {
                        if (V_mimic_mode == 0) {
                            // Global passed bin is the center of the window
                            const int thisBin = last_trig_bin + waveformCenter - waveformLength / 2 + mimicbin;
                            const int thisTimeBin = waveformCenter - waveformLength / 2 + mimicbin;

                            station_r->strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][thisBin]) *1.e3);   // save in mV
                            station_r->strings[string_i].antennas[antenna_i].time.push_back(thisBin);
                            station_r->strings[string_i].antennas[antenna_i].time_mimic.push_back(thisTimeBin *settings1->TIMESTEP *1.e9);    // save in ns
                        }
                        else if (V_mimic_mode == 1) {
                            // Global passed bin is the center of the window + delay to each chs from araGeom
                            const int thisBin = last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop]
                                                - detector->params.TestBed_BH_Mean_delay_bin)
                                                + waveformCenter - waveformLength / 2 + mimicbin;
                            const int thisTimeBin = -1*(detector->params.TestBed_Ch_delay_bin[ch_loop]
                                                    - detector->params.TestBed_BH_Mean_delay_bin)
                                                    + waveformCenter - waveformLength / 2 + mimicbin;

                            station_r->strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][thisBin]) *1.e3);   // save in mV
                            station_r->strings[string_i].antennas[antenna_i].time.push_back(thisBin);
                            station_r->strings[string_i].antennas[antenna_i].time_mimic.push_back(thisTimeBin *settings1->TIMESTEP *1.e9);   // save in ns
                        }
                        if (V_mimic_mode == 2) {
                            // Global passed bin is the center of the window + delay to each chs from araGeom + fitted delay
                            const int thisBin = last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop]
                                                - detector->params.TestBed_BH_Mean_delay_bin
                                                + station_d->strings[string_i].antennas[antenna_i].manual_delay_bin)
                                                + waveformCenter - waveformLength / 2 + mimicbin;
                            const int thisTimeBin = -1*(detector->params.TestBed_Ch_delay_bin[ch_loop]
                                                    - detector->params.TestBed_BH_Mean_delay_bin
                                                    + station_d->strings[string_i].antennas[antenna_i].manual_delay_bin)
                                                    + waveformCenter - waveformLength / 2 + mimicbin;

                            station_r->strings[string_i].antennas[antenna_i].V_mimic.push_back((trigger->Full_window_V[ch_loop][thisBin]) *1.e3);    // save in mV
                            station_r->strings[string_i].antennas[antenna_i].time.push_back(thisBin);
                            station_r->strings[string_i].antennas[antenna_i].time_mimic.push_back(thisTimeBin *settings1->TIMESTEP *1.e9 + detector->params.TestBed_WFtime_offset_ns);    // save in ns
                        }
                        if (mimicbin == 0) {
                            ///< calculates time of center of each rays signal based on readout window time config
                            for (int interaction_idx=0; interaction_idx<station_r->strings[string_i].antennas[antenna_i].SignalBin.size(); interaction_idx++){
                                for (int m = 0; m < station_r->strings[string_i].antennas[antenna_i].SignalBin[interaction_idx].size(); m++) {
                                    double signal_center_offset = (double)(
                                        station_r->strings[string_i].antennas[antenna_i].SignalBin[interaction_idx][m] -
                                        station_r->strings[string_i].antennas[antenna_i].time[0]
                                    ) * settings1->TIMESTEP * 1.e9;
                                    double signal_center_time = signal_center_offset + station_r->strings[string_i].antennas[antenna_i].time_mimic[0];
                                    //! signal_center_offset: time offset between beginning of readout window and center of signal
                                    //! signal_center_time: time of center of signal based on readout window time config
                                    station_r->strings[string_i].antennas[antenna_i].SignalBinTime[interaction_idx].push_back(signal_center_time);
                                }
                            }
                        }
                    }

                    // set global_trig_bin values
                    if (V_mimic_mode == 0) {
                        // Global passed bin is the center of the window
                        station_r->strings[string_i].antennas[antenna_i].global_trig_bin = waveformLength / 2 - waveformCenter;
                    }
                    else if (V_mimic_mode == 1) {
                        // Global passed bin is the center of the window + delay to each chs from araGeom
                        station_r->strings[string_i].antennas[antenna_i].global_trig_bin = (
                            (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin)
                            + waveformLength / 2 - waveformCenter);
                    }
                    else if (V_mimic_mode == 2) {
                        // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
                        station_r->strings[string_i].antennas[antenna_i].global_trig_bin = (
                            (detector->params.TestBed_Ch_delay_bin[ch_loop] - detector->params.TestBed_BH_Mean_delay_bin + station_d->strings[string_i].antennas[antenna_i].manual_delay_bin)
                            + waveformLength / 2 - waveformCenter);
                    }

                    // done V_mimic for non-triggered chs (done fixed V_mimic)
                }

            }

            // save trig search bin info
            station_r->total_trig_search_bin = station_r->Global_Pass + trig_window_bin - trig_search_init;

        }   // if global trig!
        else {
            this_bin++;   // also if station not passed the trigger, just go to next bin
        }
    }   // while this_bin

}

int Report::triggerCheckLoop(
    Settings *settings1, Detector *detector, Event *event, Trigger *trigger,
    int stationID, int trig_search_init, int max_total_bin, int trig_window_bin, int scan_mode
){

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

  }// for trig_j

  double *TDR_all_sorted_temp;
  double *TDR_Vpol_sorted_temp;
  double *TDR_Hpol_sorted_temp;

  if(scan_mode>1){

    TDR_all_sorted_temp=new double[numChan];
    TDR_Vpol_sorted_temp=new double[numChanVpol];
    TDR_Hpol_sorted_temp=new double[numChanHpol];

    // only need to initialize for scan_mode>1
    for(int trig_j=0;trig_j<numChan;trig_j++)
      TDR_all_sorted_temp[trig_j]=0;
    for(int trig_j=0;trig_j<numChanVpol;trig_j++)
      TDR_Vpol_sorted_temp[trig_j]=0;
    for(int trig_j=0;trig_j<numChanHpol;trig_j++)
      TDR_Hpol_sorted_temp[trig_j]=0;

  }// if scan_mode>1

  int global_pass_bit=0; // whether this event passes (in all windows)
  int check_TDR_configuration=0; // check if we need to reorder our TDR arrays
  int SCTR_cluster_bit[numChan];

  for(int trig_j=0;trig_j<numChan;trig_j++)
    SCTR_cluster_bit[trig_j]=0;

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
         (settings1->TRIG_ONLY_LOW_CH_ON==1 && settings1->DETECTOR!=3 && antenna_i<2 ) )
      { // channel filter: choose if to use lower/borehole channels or not

        int channel_num = detector->GetChannelfromStringAntenna ( i, string_i, antenna_i, settings1 );

        // assign Pthresh a value
        double diode_noise_RMS = trigger->GetAntNoise_diodeRMS(channel_num-1, settings1);
        Pthresh_value[trig_j]=trigger->Full_window[trig_j][trig_i]/(diode_noise_RMS * detector->GetThresOffset( i, channel_num-1,settings1) );

        // this is to count how many local trigger clusters there are
        if(Pthresh_value[trig_j]<powerthreshold){

          if(SCTR_cluster_bit[trig_j]==0)
            stations[i].strings[string_i].antennas[antenna_i].SingleChannelTriggers++;

          // records all the different Pthresh values that caused local trigger.
          if(settings1->TRIG_SCAN_MODE>2){

            if(SCTR_cluster_bit[trig_j]==0){// if first trigger in cluster

              stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.push_back(Pthresh_value[trig_j]);


            }
            else{// choose the highest trigger value (most negative) in cluster

              if(Pthresh_value[trig_j]<stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.back())
                stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.back()=Pthresh_value[trig_j];

            }

          }// trig scan mode > 2

          SCTR_cluster_bit[trig_j]=1;

        }// if local trigger
        else
          SCTR_cluster_bit[trig_j]=0;// if no local trigger, set zero to start a new cluster at next local trigger

        // and how many bins scanned
        stations[i].strings[string_i].antennas[antenna_i].TotalBinsScannedPerChannel++;





        // fill the buffers (if any changes occur mark check_TDR_configuration as non-zero)
        if(trig_i<trig_search_init+trig_window_bin)
          check_TDR_configuration+=buffer[trig_j]->fill(Pthresh_value[trig_j]);
        else
          check_TDR_configuration+=buffer[trig_j]->add(Pthresh_value[trig_j]);

        if(buffer[trig_j]->addToNPass>0){// if there is at least one value above threshold in the buffer, this is ++

          N_pass++;
          if(detector->stations[i].strings[string_i].antennas[antenna_i].type == 0)
            N_pass_V++;
          if(detector->stations[i].strings[string_i].antennas[antenna_i].type == 1)
            N_pass_H++;


         }// if addToNPass>0

      }// non-testbed case (i.e. use all channels)

    }// for trig_j (channel scan)

    // check if global trigger...

    if( (settings1->TRIG_MODE==0&&( N_pass >= settings1->N_TRIG )) ||
        (settings1->TRIG_MODE==1&&( N_pass_V >= settings1->N_TRIG_V ||
          N_pass_H >= settings1->N_TRIG_H ) )
      ){ // if there's a trigger !


      global_pass_bit=1;
      window_pass_bit=1;
      if(first_trigger==0){ // if this is the first trigger, mark this position and save event

        first_trigger=1;

        for(int trig_j=0;trig_j<numChan;trig_j++){

          int string_i = detector->getStringfromArbAntID( i, trig_j);
          int antenna_i = detector->getAntennafromArbAntID( i, trig_j);

          if(buffer[trig_j]->addToNPass>0)
            stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = trig_i-buffer[trig_j]->numBinsToOldestTrigger(); // mark the bin on which we triggered...
          else
            stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = 0.;

        }// for trig_j

        saveTriggeredEvent(settings1, detector, event, trigger, stationID, trig_search_init, max_total_bin, trig_window_bin, trig_i);

      }// first trigger

      if(scan_mode==1)
        return trig_i; //  if we aren't going to scan all the Pthresh values, just return
    }

    // if there's a trigger and anything changes in the buffers, restock the TDR arrays
    if( scan_mode>1 && check_TDR_configuration && window_pass_bit ){

      for(int trig_j=0;trig_j<numChan;trig_j++)
        TDR_all_sorted_temp[trig_j]=0;
      for(int trig_j=0;trig_j<numChanVpol;trig_j++)
        TDR_Vpol_sorted_temp[trig_j]=0;
      for(int trig_j=0;trig_j<numChanHpol;trig_j++)
        TDR_Hpol_sorted_temp[trig_j]=0;

      if(settings1->TRIG_MODE==0){ // for N out of all mode

        if(N_pass>=settings1->N_TRIG)
          for(int ii=0;ii<N_pass; ii++){// find the N_pass best channel's TDR and store them.

            double best_thresh=0;
            int best_chan=0;

            for(int trig_j=0;trig_j<numChan;trig_j++)
              if(buffer[trig_j]->temp_value<best_thresh){
                best_thresh=buffer[trig_j]->temp_value;
                best_chan=trig_j;
              }

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

            if(detector->stations[i].strings[string_i].antennas[antenna_i].type == 0 && buffer[trig_j]->temp_value<best_thresh){

              best_thresh=buffer[trig_j]->temp_value;
              best_chan=trig_j;

            }// if best
          }// for trig_j
          buffer[best_chan]->temp_value=0;

          TDR_Vpol_sorted_temp[ii]=best_thresh;

        }// for ii

        // for Hpol only
        if(N_pass_H>=settings1->N_TRIG_H)
          for(int ii=0;ii<N_pass_H; ii++){// find the N_pass best channel's TDR and store them.

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

      }// if trig_mode==1


      // check if temp TDR arrays improved, if so update TDR arrays:
      if(settings1->TRIG_MODE==0){

        if(N_pass>=settings1->N_TRIG) for(int ii=0;ii<N_pass;ii++) if(TDR_all_sorted_temp[ii]<stations[i].TDR_all_sorted[ii]) stations[i].TDR_all_sorted[ii]=TDR_all_sorted_temp[ii];

      }
      if(settings1->TRIG_MODE==1){

        if(N_pass_V>=settings1->N_TRIG_V) for(int ii=0;ii<N_pass_V;ii++) if(TDR_Vpol_sorted_temp[ii]<stations[i].TDR_Vpol_sorted[ii]) stations[i].TDR_Vpol_sorted[ii]=TDR_Vpol_sorted_temp[ii];
        if(N_pass_H>=settings1->N_TRIG_H) for(int ii=0;ii<N_pass_H;ii++) if(TDR_Hpol_sorted_temp[ii]<stations[i].TDR_Hpol_sorted[ii]) stations[i].TDR_Hpol_sorted[ii]=TDR_Hpol_sorted_temp[ii];

        // for this mode, can get TDR_all_sorted from these two arrays:
      }


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
        for(int p=0;p<80;p++)
          cout<<"*";
        cout<<"\n  ordering problem: "
            <<stations[i].TDR_all_sorted[0]
            <<" "
            <<stations[i].TDR_all_sorted[1]
            <<" "
            <<stations[i].TDR_all_sorted[2]
            <<"\n";
        for(int p=0;p<80;p++)
          cout<<"*";
        cout<<"\n";

      }// ordering problem


    }// trig mode 0


    if(settings1->TRIG_MODE==1){
      cout<<"\nPthresh best: ";
      cout<<"  Vpol: ";
      for(int ii=0;ii<stations[i].TDR_Vpol_sorted.size();ii++)
        cout<<" "<<stations[i].TDR_Vpol_sorted[ii];
      cout<<"  Hpol: ";
      for(int ii=0;ii<stations[i].TDR_Hpol_sorted.size();ii++)
        cout<<" "<<stations[i].TDR_Hpol_sorted[ii];
      cout<<"\n";


      // debug output:
      if(stations[i].TDR_Vpol_sorted[0]>stations[i].TDR_Vpol_sorted[1]
          || stations[i].TDR_Vpol_sorted[1]>stations[i].TDR_Vpol_sorted[2]){

        cout<<"\n";
        for(int p=0;p<80;p++) cout<<"*";
        cout<<"\n  ordering problem (final) Vpol: "
            <<stations[i].TDR_Vpol_sorted[0]
            <<" "<<stations[i].TDR_Vpol_sorted[1]
            <<" "<<stations[i].TDR_Vpol_sorted[2]
            <<"\n";
        for(int p=0;p<80;p++)
          cout<<"*";
        cout<<"\n";

      }// ordering problem

      // debug output:
      if(stations[i].TDR_Hpol_sorted[0]>stations[i].TDR_Hpol_sorted[1]
          || stations[i].TDR_Hpol_sorted[1]>stations[i].TDR_Hpol_sorted[2]){

        cout<<"\n";
        for(int p=0;p<80;p++)
          cout<<"*";
        cout<<"\n  ordering problem (final), Hpol: "
            <<stations[i].TDR_Hpol_sorted[0]
            <<" "<<stations[i].TDR_Hpol_sorted[1]
            <<" "<<stations[i].TDR_Hpol_sorted[2]
            <<"\n";
        for(int p=0;p<80;p++)
          cout<<"*";
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

int Report::saveTriggeredEvent(
    Settings *settings1, Detector *detector, Event *event, Trigger *trigger,
    int stationID, int trig_search_init, int max_total_bin, int trig_window_bin, int last_trig_bin
) {

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

      // Determine which ray and interaciton triggered triggered
      //   the station based on signal and trigger bins
      stations[i].strings[string_i].antennas[antenna_i].Find_Likely_Sol();
      int likely_int = stations[i].strings[string_i].antennas[antenna_i].Likely_Sol[0];
      int likely_ray = stations[i].strings[string_i].antennas[antenna_i].Likely_Sol[1];

        if ( settings1->TRIG_ONLY_LOW_CH_ON==0 ) {

          cout << endl
               << "trigger passed at bin "
               << stations[i].strings[string_i].antennas[antenna_i].Trig_Pass
               << "  passed ch : "<<trig_j
               << " ("<<detector->stations[i].strings[string_i].antennas[antenna_i].type
               << "type) Direct dist btw posnu : "
               << event->Nu_Interaction[0].posnu.Distance( detector->stations[i].strings[string_i].antennas[antenna_i] )
               << " noiseID : "
               << stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];
          if (stations[i].strings[string_i].antennas[antenna_i].Likely_Sol[0] != -1) {
              cout << " ViewAngle : "
                   << stations[i].strings[string_i].antennas[antenna_i].view_ang[likely_int][likely_ray] *DEGRAD ;
              cout << " LikelyTrigSignal : interaction "
                   << stations[i].strings[string_i].antennas[antenna_i].Likely_Sol[0];
              cout << ", ray "
                   << stations[i].strings[string_i].antennas[antenna_i].Likely_Sol[1];
          }
        }
        else if ( settings1->TRIG_ONLY_LOW_CH_ON==1 ) {

          cout << endl
               << "trigger passed at bin "
               << stations[i].strings[string_i].antennas[antenna_i].Trig_Pass
               <<"  passed ant: str[" << string_i << "].ant[" <<antenna_i<< "] ("
               << detector->stations[i].strings[string_i].antennas[antenna_i].type
               << "type) Direct dist btw posnu : "
               << event->Nu_Interaction[0].posnu.Distance( detector->stations[i].strings[string_i].antennas[antenna_i] )
               << " noiseID : "
               << stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];

          if (stations[i].strings[string_i].antennas[antenna_i].Likely_Sol[0] != -1) {
              cout << " ViewAngle : "
                   << stations[i].strings[string_i].antennas[antenna_i].view_ang[likely_int][likely_ray] *DEGRAD ;
              cout << " LikelyTrigSignal : interaction "
                   << stations[i].strings[string_i].antennas[antenna_i].Likely_Sol[0];
              cout << ", ray "
                   << stations[i].strings[string_i].antennas[antenna_i].Likely_Sol[1];
          }

        }

    }// if Trig_Pass

    // now save the voltage waveform to V_mimic
    for (int mimicbin=0; mimicbin < waveformLength/2; mimicbin++) {

        // new DAQ waveform writing mechanism test
        if (settings1->V_MIMIC_MODE == 0) { // Global passed bin is the center of the window

          const int thisBin = last_trig_bin - waveformLength/2 + waveformCenter + mimicbin;
          const int thisTimeBin = waveformLength/2 + waveformCenter + mimicbin;

          stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back( ( trigger->Full_window_V[trig_j][thisBin] )*1.e3 );// save in mV
          stations[i].strings[string_i].antennas[antenna_i].time.push_back(thisBin);
          stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back( thisTimeBin * settings1->TIMESTEP*1.e9  );// save in ns
        }
        else if (settings1->V_MIMIC_MODE == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom

          const int thisBin = last_trig_bin - (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin) - waveformLength/2 + waveformCenter + mimicbin;
          const int thisTimeBin = -(detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin) - waveformLength/2 + waveformCenter + mimicbin;

          stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back( ( trigger->Full_window_V[trig_j][thisBin] )*1.e3 );// save in mV
          stations[i].strings[string_i].antennas[antenna_i].time.push_back(thisBin);
          stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back( thisTimeBin * settings1->TIMESTEP*1.e9  );// save in ns

        }

        else if (settings1->V_MIMIC_MODE == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye

          const int thisBin = last_trig_bin - (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin
                              + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) - waveformLength/2 + waveformCenter + mimicbin;
          const int thisTimeBin = -1*(detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin
                                  + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin) - waveformLength/2 + waveformCenter + mimicbin;

          stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back( ( trigger->Full_window_V[trig_j][thisBin] )*1.e3 );// save in mV
          stations[i].strings[string_i].antennas[antenna_i].time.push_back(thisBin);
          stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back( thisTimeBin * settings1->TIMESTEP*1.e9 + detector->params.TestBed_WFtime_offset_ns );// save in ns

        }
        if (mimicbin == 0) {
            ///< calculates time of center of each rays signal based on readout window time config
            for (int interaction_idx=0; interaction_idx < stations[i].strings[string_i].antennas[antenna_i].SignalBin.size(); interaction_idx++){
                for (int m = 0; m < stations[i].strings[string_i].antennas[antenna_i].SignalBin[interaction_idx].size(); m++) {
                    double signal_center_offset = (double)(
                        stations[i].strings[string_i].antennas[antenna_i].SignalBin[interaction_idx][m] -
                        stations[i].strings[string_i].antennas[antenna_i].time[0]
                    ) * settings1->TIMESTEP * 1.e9;
                    double signal_center_time = signal_center_offset + stations[i].strings[string_i].antennas[antenna_i].time_mimic[0];
                    //! signal_center_offset: time offset between beginning of readout window and center of signal
                    //! signal_center_time: time of center of signal based on readout window time config
                    stations[i].strings[string_i].antennas[antenna_i].SignalBinTime[interaction_idx].push_back(signal_center_time);
                }
            }
         }

      }


       // set global_trig_bin values
       if (settings1->V_MIMIC_MODE == 0) { // Global passed bin is the center of the window
         stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = waveformLength/2 + waveformCenter ;
       }
       else if (settings1->V_MIMIC_MODE == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
         stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin)
                                                                             + waveformLength/2 + waveformCenter ;
       }
       else if (settings1->V_MIMIC_MODE == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
          stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin
                                                                              + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
                                                                              + waveformLength/2 + waveformCenter ;
       }

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
       for(int trig_i=0;trig_i<settings1->DATA_BIN_SIZE/2;trig_i++)
          gr[trig_j]->SetPoint(trig_i, (trig_search_init+trig_i), trigger->Full_window[trig_j][trig_i]);

        gr[trig_j]->SetNameTitle(Form("TDR_waveform%dC%02d",  settings1->OUTPUT_TDR_GRAPH, trig_j),
                                 Form("Tunnel diode response waveform %d, channel %02d, trig_pass= %d, P_{th}= %le; time bins; power after convolution with tunnel diode",
                                 settings1->OUTPUT_TDR_GRAPH, trig_j, stations[i].strings[string_i].antennas[antenna_i].Trig_Pass, thresh_value));
        gr[trig_j]->Write();
        delete gr[trig_j];
       }

       delete [] gr;


     }

  return 1;

}// saveTriggeredEvent

#ifdef ARA_UTIL_EXISTS
void Report::MakeUsefulEvent(Detector *detector, Settings *settings1, Trigger *trigger, int stationID, int stationIndex, UsefulIcrrStationEvent *theUsefulEvent) {
    if (stationID < detector->stations.size()){
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
            theUsefulEvent->fNumPointsRF[ch_loop] = UsefulEventBin;
        }
    }
}
#endif

#ifdef ARA_UTIL_EXISTS
void Report::MakeUsefulEvent(Detector *detector, Settings *settings1, Trigger *trigger, int stationID, int stationIndex, UsefulAtriStationEvent *theUsefulEvent) {

    int i = stationID;
	int stationID_AraRoot = settings1->DETECTOR_STATION_ARAROOT;
	cout << "StationID: " << stationID << endl;
	cout << "StationID_AraRoot: " << stationID_AraRoot << endl;
	theUsefulEvent->fNumChannels = 32;
	theUsefulEvent->stationId = stationID_AraRoot;

	int ch_limit;
	if (stationID == 0){
	  ch_limit = 14;
	} else {
	  ch_limit = 16;
	}

	int maxElecChans = 32;

	for (int ch_loop=0; ch_loop < ch_limit; ch_loop++) {
	  int elecChan = AraGeomTool::Instance()->getElecChanFromRFChan(ch_loop, stationID_AraRoot);
	  int string_i = 0;
	  int antenna_i = 0;
	  detector->GetSSAfromChannel(stationID, ch_loop, &antenna_i, &string_i, settings1);

	  int AraRootChannel = 0;
	  AraRootChannel = detector->GetChannelfromStringAntenna (stationID, string_i, antenna_i, settings1);

	  int UsefulEventBin;

	  UsefulEventBin = settings1->WAVEFORM_LENGTH;

	  vector < double > volts;
	  volts.resize(UsefulEventBin);
	  vector < double > times;
	  times.resize(UsefulEventBin);

	  for (int mimicbin=0; mimicbin<UsefulEventBin; mimicbin++) {
	    if (stations[stationIndex].Global_Pass > 0){
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

	  volts.clear();
	  times.clear();

    }
}
#endif


void Report::ClearUselessfromConnect(Detector *detector, Settings *settings1, Trigger *trigger){

    for (int i = 0; i< stations.size(); i++) {
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
    Event *event, Settings *settings1, Trigger *trigger, Detector *detector
){
    // For the provided antenna:
    // Calculate bin values for the signal and determine if how many rays fit in the readout window
    // Combine signals from rays and load noise waveforms
    // Pass the noise+signal waveform through the tunnel diode
    // Enforce voltage saturation

    // init : bin + maxt_diode_bin, fin : bin + NFOUR/2

    // Clear old signals
    vector < vector <int> > signal_bin; // the center of bin where signal should locate
    signal_bin.resize(event->Nu_Interaction.size());

    // Determine the bin (index) where the signal will arrive in output waveforms
    //   and flag which signals will fit into the same time window
    for (int interaction_idx=0; interaction_idx<event->Nu_Interaction.size(); interaction_idx++) {
        for (int m = 0; m < antenna->arrival_time[interaction_idx].size(); m++) { // loop over raysol numbers
            // Store the bin where the signal is located
            signal_bin[interaction_idx].push_back(
                (antenna->arrival_time[interaction_idx][m] - stations[station_number].min_arrival_time) / (settings1->TIMESTEP)
                + settings1->NFOUR *2 + trigger->maxt_diode_bin);
            antenna->SignalBin[interaction_idx].push_back(signal_bin[interaction_idx][m]);
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
    // Loop over ray solutions and get signals from each
    // Loop does not run if there are no ray solutions
    else {

      int this_signalbin = 0;
      int n_connected_rays = antenna->ray_sol_cnt;
      vector<double> V_signal;

      // Identify the indices of the last non-empty waveform (to apply final padding inside Combine_Waveforms)
      // Loop through all interaction indices
      int last_valid_interaction = -1;  // Will store the index of the last interaction with data
      int last_valid_m = -1;            // Will store the index of the last ray within that interaction

      for (int i = 0; i < antenna->V.size(); ++i) {
          // Loop over rays for this interaction
          for (int m = 0; m < signal_bin[i].size(); ++m) {
              // Check if the waveform is non-empty
              if (!antenna->V[i][m].empty()) {
                  // Update the last known valid (non-empty) waveform indices
                  last_valid_interaction = i;
                  last_valid_m = m;
              }
          }
      }

      for (int interaction_idx=0; interaction_idx<antenna->V.size(); interaction_idx++){
        for (int m = 0; m < signal_bin[interaction_idx].size(); ++m) {

            bool is_last = (interaction_idx == last_valid_interaction) && (m == last_valid_m); // reached the last waveform

            Combine_Waveforms(
                signal_bin[interaction_idx][m], this_signalbin,
                antenna->V[interaction_idx][m], V_signal,
                &this_signalbin, &V_signal, is_last);
        }
      }

      if(V_signal.size() >= settings1->DATA_BIN_SIZE){
        throw runtime_error("Full waveform trace is longer than DATA_BIN_SIZE.");
      }

      // If there are ray solutions, prepare signal wf array with 0s
      //   and save an array with the full noise-only waveform to the antenna

      // Initialize array_length = waveform length + (2*NFOUR + maxt_diode_bin)
      int waveform_length = (int)V_signal.size();
      int array_length = waveform_length + this_signalbin + 2*settings1->NFOUR + trigger->maxt_diode_bin;

      // Make vectors of 0s for array we're saving all ray signals to
      antenna->V_convolved.clear();
      antenna->V_noise.clear();
      antenna->V_convolved.resize(array_length, 0.0);
      antenna->V_noise.resize(array_length, 0.0);

      GetNoiseThenConvolve(
          antenna, V_signal,
          BINSIZE, this_signalbin, n_connected_rays,
          channel_index, station_number,
          settings1, trigger, detector);

    }

}

int Antenna_r::Get_Max_SignalBin(){
    int max_signal_bin = 0;
    for (int interaction_idx=0; interaction_idx<SignalBin.size(); interaction_idx++) {
        for (int ray_idx=0; ray_idx<SignalBin[interaction_idx].size(); ray_idx++) {
            if (SignalBin[interaction_idx][ray_idx] > max_signal_bin) {
                max_signal_bin = SignalBin[interaction_idx][ray_idx];
            }
        }
    }
    return max_signal_bin;
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

void Report::Combine_Waveforms(int signalbin_0, int signalbin_1,
    vector<double> V0, vector<double> V1,
    int* signalbin_combined, vector<double>* V_combined, bool pad_to_power_of_two
) {
  // adds waveforms V0 & V1 into combined vector V_combined
  // uses global signalbin indices to ensure vectors are properly combined
  // V_combined is not guaranteed to have a length which is a multiple of BINSIZE
  // signalbin_combined tracks global signalbin index of combined vector for use later

  vector<double>& V = *V_combined;
  int& signalbin = *signalbin_combined;

  // ensure voltage array is empty
  V.clear();

  // check if both input vectors are empty
  if(V0.empty() && V1.empty())
    throw runtime_error("Cannot combine two empty signal vectors!");

  // if one vector is empty, assign its signalbin to that of the other so nothing breaks
  if(V0.empty())
    signalbin_0 = signalbin_1;
  else if(V1.empty())
    signalbin_1 = signalbin_0;

  // resize to necessary length to combine vectors
  signalbin = min(signalbin_0, signalbin_1); // starting index in terms of global index
  const int maxsignalbin = max(signalbin_0+V0.size()-1, signalbin_1+V1.size()-1); // last index in terms of global index
  const int len = maxsignalbin - signalbin + 1; // length needed to combine both signals
  V.resize(len, 0);

  // add in the zeroth vector
  for(unsigned int bin = 0; bin < V0.size(); ++bin) {
    const int combined_bin = bin + (signalbin_0- signalbin);
    V[combined_bin] += V0[bin];
  }

  // add in the first vector
  for(unsigned int bin = 0; bin < V1.size(); ++bin) {
    const int combined_bin = bin + (signalbin_1- signalbin);
    V[combined_bin] += V1[bin];
  }

  // for the last added waveform if length is not a power of 2, zero pad
  if(pad_to_power_of_two && (V.size() & (V.size() - 1)) != 0) {
    const int newLen = int(pow(2, ceil(log2(V.size())))); // get next power of 2
    V.resize(newLen, 0.);
  }

  return;
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
    int min_wf_bin = 0;
    int max_wf_bin = 0;
    int offset = 0;
    vector <double> diode_response;
    if ( n_connected_rays > 1 ) { // multiple ray solutions in one window
        wf_length = V_signal.size(); // when using Select_Wave_Convlv_Exchange this is 2*BINSIZE
        offset = trigger->maxt_diode_bin;
        min_wf_bin = this_signalbin - BINSIZE/2 + offset;
        max_wf_bin = this_signalbin - BINSIZE/2 + wf_length;
        diode_response = detector->getDiodeModel(2*wf_length, settings1);
    }
    else if ( antenna->ray_sol_cnt == 0 ){ // No rays connected to this antenna
        this_signalbin = BINSIZE/2;
        wf_length = V_signal.size(); // when using Select_Wave_Convlv_Exchange this is BINSIZE
        offset = 0;
        min_wf_bin = 0;
        max_wf_bin = wf_length;
        diode_response = detector->getDiodeModel(2*wf_length, settings1);
    }
    else { // Only one ray signal in the window
        wf_length = V_signal.size(); // when using Select_Wave_Convlv_Exchange this is BINSIZE
        offset = trigger->maxt_diode_bin;
        min_wf_bin = this_signalbin - BINSIZE/2 + offset;
        max_wf_bin = this_signalbin - BINSIZE/2 + wf_length;
        diode_response = detector->getDiodeModel(2*wf_length, settings1);
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
        const int vBin = bin - (min_wf_bin - offset);
        trigger->Full_window[channel_index][bin] = V_total_forconvlv[vBin];
        trigger->Full_window_V[channel_index][bin] += V_signal[vBin];
        antenna->V_convolved[bin] += V_signal[vBin];
        antenna->V_noise[bin] += V_noise[vBin];

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

    int channel_num = detector->GetChannelfromStringAntenna ( StationIndex, string_num, ant_num, settings1 );

    for (int bin=0; bin<settings1->DATA_BIN_SIZE; bin++) {   // test for full window
        trigger->Full_window[ID][bin] = (
            trigger->Full_window[ID][bin] * detector->GetGainOffset( StationIndex, channel_num-1, settings1 )
            * detector->GetGainOffset( StationIndex, channel_num-1, settings1 )
        ); // offset in voltage factor so we need power (V^2 factor to diode response)
        trigger->Full_window_V[ID][bin] = (
            trigger->Full_window_V[ID][bin] * detector->GetGainOffset( StationIndex, channel_num-1, settings1 )
        ); // gain in voltage factor
    }

}



int Report::GetChNumFromArbChID( Detector *detector, int ID, int StationIndex, Settings *settings1 ) {

    int string_num = detector->getStringfromArbAntID( StationIndex, ID );
    int ant_num = detector->getAntennafromArbAntID( StationIndex, ID );
    int StationID = detector->stations[StationIndex].StationID;

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


void Report::GetParameters(
    Position &src, Position &trg, Vector &nnu, double &viewangle, double receive_angle,
    Vector &launch_vector, Vector &receive_vector, Vector &n_trg_slappy, Vector &n_trg_pokey
) {

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


void Report::ApplyAntFactors(
    double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector,
    int ant_type, double &pol_factor, double &vmmhz, double antenna_theta, double antenna_phi
) {  // vmmhz is input and output. output will have some antenna factors on it

    pol_factor = calculatePolFactor(Pol_vector, ant_type, antenna_theta, antenna_phi);

    // apply 3dB spliter, d nu to prepare FFT
    // now actually vmmhz is not V/m/MHz but V/m/Hz unit
    // sqrt(2) for 3dB spliter for TURF, SURF. 1/TIMESTEP moved to MakeArraysforFFT
    // 1/(1.E6) for V/MHz to V/Hz
    vmmhz = vmmhz/sqrt(2.)/(1.E6);

    // apply antenna effective height and 0.5 (to calculate power with heff), and polarization factor
    // not vmmhz is actually V/Hz unit
    vmmhz = vmmhz * 0.5 * heff * pol_factor;

    // now we have to use MakeArraysforFFT to change signal arrays for FFT
}


void Report::ApplyAntFactors_Tdomain(
    double AntPhase, double heff, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img,
    Settings *settings1, double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode, bool applyInverse
) {
    double heff_copol = heff;
    double heff_crosspol = 0.0;
    double phase_copol = AntPhase;
    double phase_crosspol = 0.0;

    ApplyAntFactors_Tdomain(
        phase_copol, phase_crosspol,
        heff_copol, heff_crosspol,
        Pol_vector, ant_type,
        pol_factor, vm_real, vm_img,
        settings1, antenna_theta, antenna_phi,
        freq, useInTransmitterMode, applyInverse);
}


void Report::ApplyAntFactors_Tdomain(double phase_copol, double phase_crosspol, double heff_copol, double heff_crosspol,
                                     Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img,
                                     Settings *settings1, double antenna_theta, double antenna_phi, double freq,
                                     bool useInTransmitterMode, bool applyInverse
) {
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
    //Do both co-pol and cross-pol polarization factors
    double pol_factor_copol = calculatePolFactor(Pol_vector, ant_type, antenna_theta, antenna_phi);
    double pol_factor_crosspol = calculatePolFactor(Pol_vector, 1 - ant_type, antenna_theta, antenna_phi);

    // Calculate amplitude amplifications from co-pol and cross-pol
    double v_amplification_copol = heff_copol * pol_factor_copol;
    double v_amplification_crosspol = 0.0;

    //turn on cross-pol for receiver
    if (settings1->CROSSPOL_RX){
        v_amplification_crosspol = heff_crosspol * pol_factor_crosspol;
    }

     // Add the contributions linearly for Rx
     double v_amplification = v_amplification_copol + v_amplification_crosspol;
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
        double v_amp_in  = sqrt(vm_real*vm_real + vm_img*vm_img);
        double v_amp = v_amp_in; //save the incomming voltage
        // Apply amplitude sign for inversion if necessary
        v_amp *= pow(v_amplification, amplitudeSign);
        //If in transmitter mode, we must apply additional frequency and impedance terms to the amplitude.
        if (useInTransmitterMode==true){
            phase_current += PI/2;
            double psi = settings1->CLOCK_ANGLE;
            double delta_psi = TMath::ATan(heff_crosspol / heff_copol); // Tx cross-pol tilt
            double theta = antenna_theta*PI/180;
            double phi = antenna_phi*PI/180;

            //turn on cross-pol Tx
            if (settings1->CROSSPOL_TX != 1){
                heff_crosspol = 0.0; //need to be turned off for the calculation of amplitude v_amp
                delta_psi = 0.0; // turn off cross-pol tx tilt
            }

            //Adjust polarization vector
            double newPol_vectorX = -cos(psi+delta_psi)*cos(theta)*cos(phi) + sin(psi+delta_psi)*sin(phi);
            double newPol_vectorY = -cos(psi+delta_psi)*cos(theta)*sin(phi) - sin(psi+delta_psi)*cos(phi);
            double newPol_vectorZ = cos(psi+delta_psi)*sin(theta);

            Pol_vector = Vector(newPol_vectorX, newPol_vectorY, newPol_vectorZ);

            //copol and cross-pol add quadratically in E-field
            v_amp *= pow(freq/CLIGHT*(Z0/Zr)/4/sqrt(2.)*(1/v_amplification)*(sqrt(heff_crosspol* heff_crosspol + heff_copol*heff_copol )), amplitudeSign);
        }
        // Calculate the combined real and imaginary terms from co-pol and cross-pol
        double vm_real_copol = v_amp_in * (v_amplification_copol * cos(phase_current + (phaseSign * phase_copol * RADDEG)));
        double vm_img_copol = v_amp_in * (v_amplification_copol * sin(phase_current + (phaseSign * phase_copol * RADDEG)));

        double vm_real_crosspol = v_amp_in * (v_amplification_crosspol * cos(phase_current + (phaseSign * phase_crosspol * RADDEG)));
        double vm_img_crosspol = v_amp_in * (v_amplification_crosspol * sin(phase_current + (phaseSign * phase_crosspol * RADDEG)));

        // Combine co-pol and cross-pol real and imaginary parts
        vm_real = vm_real_copol + vm_real_crosspol;
        vm_img = vm_img_copol + vm_img_crosspol;
    }
    else { // Amplitude-only mode

        vm_real *= pow(v_amplification, amplitudeSign);
        vm_img *= pow(v_amplification, amplitudeSign);
    }
}

void Report::ApplyAntFactors_Tdomain_FirstTwo(
    double heff, double heff_lastbin, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_bin0, double &vm_bin1,
    Settings *settings1, double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode, bool applyInverse
) {
    // Default to no cross-pol
    double heff_crosspol = 0.0;
    double heff_crosspol_lastbin = 0.0;

    ApplyAntFactors_Tdomain_FirstTwo(
        heff, heff_lastbin, heff_crosspol, heff_crosspol_lastbin,
        Pol_vector, ant_type, pol_factor, vm_bin0, vm_bin1,
        settings1, antenna_theta, antenna_phi, freq, useInTransmitterMode, applyInverse);
}


void Report::ApplyAntFactors_Tdomain_FirstTwo(double heff_copol, double heff_copol_lastbin, double heff_crosspol,
                                              double heff_crosspol_lastbin, Vector &Pol_vector, int ant_type,
                                              double &pol_factor, double &vm_bin0, double &vm_bin1, Settings *settings1,
                                              double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode,
                                              bool applyInverse
) {
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

    // Initialize cross-pol amplification factors
    double pol_factor_crosspol = 0.0;
    double v_amplification_copol_bin0 = heff_copol * pol_factor;
    double v_amplification_copol_bin1 = heff_copol_lastbin * pol_factor;
    double v_amplification_crosspol_bin0 = 0.0;
    double v_amplification_crosspol_bin1 = 0.0;

    // Handle cross-pol factors for receiver mode
    if (!useInTransmitterMode && settings1->CROSSPOL_RX) {
        pol_factor_crosspol = calculatePolFactor(Pol_vector, 1 - ant_type, antenna_theta, antenna_phi);
        v_amplification_crosspol_bin0 = heff_crosspol * pol_factor_crosspol;
        v_amplification_crosspol_bin1 = heff_crosspol_lastbin * pol_factor_crosspol;
    }

    // Combine co-pol and cross-pol amplifications
    double v_amplification_bin0 = v_amplification_copol_bin0 + v_amplification_crosspol_bin0;
    double v_amplification_bin1 = v_amplification_copol_bin1 + v_amplification_crosspol_bin1;

    // Apply amplification for receiver mode
    vm_bin0 *= pow(v_amplification_bin0, amplitudeSign);
    vm_bin1 *= pow(v_amplification_bin1, amplitudeSign);

    if (useInTransmitterMode) {
        // Handle cross-pol factors for transmitter mode
        if (settings1->CROSSPOL_TX) {
            pol_factor_crosspol = calculatePolFactor(Pol_vector, 1 - ant_type, antenna_theta, antenna_phi);
            v_amplification_crosspol_bin0 = heff_crosspol * pol_factor_crosspol;
            v_amplification_crosspol_bin1 = heff_crosspol_lastbin * pol_factor_crosspol;
        }

        // Co-pol and cross-pol add quadratically in E-field
        double tx_amplification_bin0 = freq / CLIGHT * (Z0 / Zr) / 4 / sqrt(2.0) * (1/v_amplification_bin0) *
                                       (sqrt(heff_crosspol * heff_crosspol + heff_copol * heff_copol));
        double tx_amplification_bin1 = freq / CLIGHT * (Z0 / Zr) / 4 / sqrt(2.0) * (1/v_amplification_bin1) *
                                       (sqrt(heff_crosspol_lastbin * heff_crosspol_lastbin + heff_copol_lastbin * heff_copol_lastbin));

        vm_bin0 *= pow(tx_amplification_bin0, amplitudeSign);
        vm_bin1 *= pow(tx_amplification_bin1, amplitudeSign);
    }

    // I don't understand why we had it twice in the original function -ASG 12/09/24
    //    else {
    //        vm_bin0 *= pow(heff * pol_factor, amplitudeSign);
    //        vm_bin1 *= pow(heff_lastbin * pol_factor, amplitudeSign);
    //    }

}

void Report::InvertAntFactors_Tdomain(
    double AntPhase, double heff, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img, Settings *settings1,
    double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode
) {
    ApplyAntFactors_Tdomain(
        AntPhase, heff, Pol_vector, ant_type,
        pol_factor, vm_real, vm_img,
        settings1, antenna_theta, antenna_phi,
        freq, useInTransmitterMode, true); // applyInverse = true
}

void Report::InvertAntFactors_Tdomain(double AntPhase_copol, double AntPhase_crosspol,
                                          double heff_copol, double heff_crosspol,
                                          Vector &Pol_vector, int ant_type,
                                          double &pol_factor, double &vm_real, double &vm_img,
                                          Settings *settings1, double antenna_theta,
                                          double antenna_phi, double freq,
                                          bool useInTransmitterMode) {
    /* Report::InvertAntFactors_Tdomain_new()
    Purpose: Inverts the antenna factors for both co-polarization and cross-polarization
             by simply calling ApplyAntFactors_Tdomain_new with the boolean applyInverse enabled.
    */
    ApplyAntFactors_Tdomain(AntPhase_copol, AntPhase_crosspol, heff_copol, heff_crosspol,
                                Pol_vector, ant_type, pol_factor, vm_real, vm_img,
                                settings1, antenna_theta, antenna_phi, freq, useInTransmitterMode, true);

}

void Report::InvertAntFactors_Tdomain_FirstTwo(
    double heff, double heff_lastbin, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_bin0, 
    double &vm_bin1, Settings *settings1, double antenna_theta, double antenna_phi, double freq, bool useInTransmitterMode
) {
    ApplyAntFactors_Tdomain_FirstTwo(
        heff, heff_lastbin,
        Pol_vector, ant_type, pol_factor, vm_bin0, vm_bin1,
        settings1, antenna_theta, antenna_phi, freq,
        useInTransmitterMode, true);
}

void Report::InvertAntFactors_Tdomain_FirstTwo(double heff_copol, double heff_copol_lastbin,
                                                   double heff_crosspol, double heff_crosspol_lastbin,
                                                   Vector &Pol_vector, int ant_type,
                                                   double &pol_factor, double &vm_bin0, double &vm_bin1,
                                                   Settings *settings1, double antenna_theta, double antenna_phi, double freq,
                                                   bool useInTransmitterMode) {
    /* Report::InvertAntFactors_Tdomain_FirstTwo_new()
    Purpose: Inverts the antenna factors for both co-polarization and cross-polarization
             by simply calling ApplyAntFactors_Tdomain_FirstTwo_new with the boolean applyInverse enabled.
    */

    ApplyAntFactors_Tdomain_FirstTwo(heff_copol, heff_copol_lastbin, heff_crosspol, heff_crosspol_lastbin,
                                         Pol_vector, ant_type, pol_factor, vm_bin0, vm_bin1,
                                         settings1, antenna_theta, antenna_phi, freq, useInTransmitterMode, true);

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
    // If using this function to invert the electronics response, we apply a minus sign the phase in order to undo the
    //   phase applied by the electronics.  We also divide the amplitude by the factor rather than multiply.
    //   To do this, we simply apply a -1 to the power of the factor applied to the amplitude so that it is divided out.
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

        // V amplitude, apply gain (unitless) to amplitude
        double v_amp  = sqrt(vm_real*vm_real + vm_img*vm_img) * pow(detector->GetElectGain_1D_OutZero( freq, gain_ch_no ),amplitudeSign);

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




void Report::ApplyElect_Tdomain_FirstTwo(
    double freq0, double freq1, Detector *detector, double &vm_bin0, double &vm_bin1,
    int gain_ch_no, Settings *settings1, bool applyInverse
) {  // read elect chain gain (unitless), phase (rad) and apply to V/m

    /*  Report::ApplyElect_Tdomain_FirstTwo()
    Purpose:
        Multiply (or divide) first two bins of voltage in Fourier space by the electronics gain and phase.

    For a more verbose explanation, see Report::ApplyElect_Tdomain.
    This simply applies the operation to the first and last bin, which much be done separately as a consequence of FFT's.
    Does not calculate phase.
    */
    double amplitudeSign = 1.;
    // If using this function to invert the electronics response, we apply a minus sign the phase in order to undo the
    //   phase applied by the electronics.  We also divide the amplitude by the factor rather than multiply.
    //   To do this, we simply apply a -1 to the power of the factor applied to the amplitude so that it is divided out.
    if(applyInverse==true){amplitudeSign*=-1;};

    vm_bin0 = vm_bin0 * pow(detector->GetElectGain_1D_OutZero( freq0 , gain_ch_no),amplitudeSign);
    vm_bin1 = vm_bin1 * pow(detector->GetElectGain_1D_OutZero( freq1 , gain_ch_no),amplitudeSign);

    //Apply power splitter/attenuator based on station.
    ApplySplitterFactor(vm_bin0, vm_bin1, detector, settings1, applyInverse);
}

void Report::ApplySplitterFactor(double &vm_real, double &vm_img, Detector* detector, Settings *settings1, bool applyInverse) {
    // Apply splitter/attenuation factor in the digitizer path based on station.

    // If we wish to invert the splitter factor, we must divide the factor out rather than multiply.  To accomplish this,
    //   we simply raise the factor to the -1 power in the multiplication step.
    double amplitudeSign=1;
    if (applyInverse==true){amplitudeSign*=-1;};

    const double factor = detector->GetSplitterFactor(settings1);

    vm_real *= pow(factor, amplitudeSign);
    vm_img *= pow(factor, amplitudeSign);

    return;
}

void Report::InvertElect_Tdomain(
    double freq, Detector *detector, double &vm_real, double &vm_img, int gain_ch_no, Settings *settings1
) {  // read elect chain gain (unitless), phase (rad) and apply to V/m

    /* Report::InvertElect_Tdomain()
    Purpose:  Inverts the electronics factors in a convenient way by simply calling ApplyElecFactors with the boolean applyInverse enabled.
    */
    ApplyElect_Tdomain(freq, detector, vm_real, vm_img, gain_ch_no, settings1, true);

}

void Report::InvertElect_Tdomain_FirstTwo(
    double freq0, double freq1, Detector *detector, double &vm_bin0, double &vm_bin1, int gain_ch_no, Settings *settings1
) {  // read elect chain gain (unitless), phase (rad) and apply to V/m

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

void Report::ApplyRFCM_OutZero(int ch, double freq, Detector *detector, double &vmmhz, double RFCM_OFFSET) {  // read RFCM gain in dB and apply unitless gain to vmmhz

    vmmhz = vmmhz * pow(10., ( detector->GetRFCMGain_OutZero(ch, freq) + RFCM_OFFSET )/20.);   // from dB to unitless gain for voltage

    return;
}

void Report::ApplyNoiseFig_databin(int ch, int bin_n, Detector *detector, double &vmmhz, Settings *settings1) {  // read noise figure and apply unitless gain to vmmhz

  double tempNoise = vmmhz*vmmhz;
  if(detector->GetNoiseFig_databin(ch,bin_n)>1.0){
      vmmhz = TMath::Sqrt( tempNoise*(detector->GetTransm_databin(ch, bin_n)) + tempNoise/(settings1->NOISE_TEMP)*220.0*(detector->GetNoiseFig_databin(ch,bin_n) - 1.0) ) ;
  }
  else{
      vmmhz = vmmhz;
  }

  return;
}

void Report::ApplyNoiseFig_OutZero(int ch, double freq, Detector *detector, double &vmmhz, Settings *settings1) {  // read noise figure and apply unitless gain to vmmhz

    double tempNoise = vmmhz*vmmhz;

    if(detector->GetNoiseFig_OutZero(ch, freq)>1.0){
        vmmhz = TMath::Sqrt( tempNoise*(detector->GetTransm_OutZero(ch, freq)) + tempNoise/(settings1->NOISE_TEMP)*220.0*(detector->GetNoiseFig_OutZero(ch, freq) - 1.0) ) ;
    }
    else{
      vmmhz = vmmhz;
    }

    return;
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

        Vfft_noise_after.clear();  // remove previous Vfft_noise values
        Vfft_noise_before.clear();  // remove previous Vfft_noise values

        double V_tmp; // copy original flat H_n [V] value
        double current_amplitude, current_phase;

        GetNoisePhase(settings1); // get random phase for noise

        // returned array vnoise currently have real value = imag value as no phase term applied

        for (int k=0; k<settings1->DATA_BIN_SIZE/2; k++) {

            V_tmp = v_noise / sqrt(2.) ; // copy original flat H_n [V] value, and apply 1/sqrt2 for SURF/TURF divide same as signal

            if (settings1->USE_TESTBED_RFCM_ON == 0) {
                const double dfreq = 1./(settings1->DATA_BIN_SIZE * settings1->TIMESTEP * (1.E6)); // from Hz to MHz;
                const double freq = k*dfreq;
                ApplyFilter_OutZero(freq, detector, V_tmp);
                ApplyPreamp_OutZero(freq, detector, V_tmp);
                ApplyFOAM_OutZero(freq, detector, V_tmp);
                if (settings1->APPLY_NOISE_FIGURE==1){
                    ApplyNoiseFig_OutZero(0, freq, detector, V_tmp, settings1);
                }
            }
            else if (settings1->USE_TESTBED_RFCM_ON == 1) {
                // apply RFCM gain
                // as this mode don't have different chs' noise waveform separately, just use ch0 RFCM gain...
                cerr<<"Trying not allowed mode : NOISE_CHANNEL_MODE=0 and USE_TESTBED_RFCM_ON=1 same time!"<<endl;
                break;
            }

            Vfft_noise_before.push_back( V_tmp );

            current_phase = noise_phase[k];

            // use real value array value, extra 1/1.177 to make total power same with "before random_rician".
            Tools::get_random_rician( 0., 0., sqrt(2./M_PI)/1.177 * V_tmp, current_amplitude, current_phase);

            // vnoise is currently noise spectrum (before fft, unit : V)
           vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
           vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);

            Vfft_noise_after.push_back( vnoise[2*k] );
            Vfft_noise_after.push_back( vnoise[2*k+1] );

            // inverse FFT normalization factor!
            vnoise[2 * k] *= 2./((double)settings1->DATA_BIN_SIZE);
            vnoise[2 * k + 1] *= 2./((double)settings1->DATA_BIN_SIZE);

        }

        // now vnoise is time domain waveform
        Tools::realft( vnoise, -1, settings1->DATA_BIN_SIZE);

    }
    else {  // currently there are no more options!
        cout<<"no noise option for NOISE = "<<settings1->NOISE<<endl;
    }


}



// generate DATA_BIN_SIZE sized noise waveform array in time domain
void Report::GetNoiseWaveforms_ch(Settings * settings1, Detector * detector, double v_noise, double * vnoise, int ch) {

    if (settings1 -> NOISE == 0) { // NOISE == 0 : flat thermal noise with Johnson-Nyquist noise

        Vfft_noise_after.clear(); // remove previous Vfft_noise values
        Vfft_noise_before.clear(); // remove previous Vfft_noise values

        Vfft_noise_after.resize(settings1->DATA_BIN_SIZE);
        Vfft_noise_before.resize(settings1->DATA_BIN_SIZE/2);

        double V_tmp; // copy original flat H_n [V] value
        double current_amplitude, current_phase;

        GetNoisePhase(settings1); // get random phase for noise

        for (int k = 0; k < settings1 -> DATA_BIN_SIZE / 2; k++) {

            V_tmp = v_noise / sqrt(2.); // copy original flat H_n [V] value, and apply 1/sqrt2 for SURF/TURF divide same as signal

            if (settings1 -> USE_TESTBED_RFCM_ON == 0) {
                const double dfreq = 1./(settings1->DATA_BIN_SIZE * settings1->TIMESTEP * (1.E6)); // from Hz to MHz;
                const double freq = k*dfreq;
                ApplyFilter_OutZero(freq, detector, V_tmp);
                ApplyPreamp_OutZero(freq, detector, V_tmp);
                ApplyFOAM_OutZero(freq, detector, V_tmp);
                if (settings1 -> APPLY_NOISE_FIGURE == 1) {
                    ApplyNoiseFig_OutZero(ch % 16, freq, detector, V_tmp, settings1);
                }
            }
            else if (settings1 -> USE_TESTBED_RFCM_ON == 1) {
                const double dfreq = 1./(settings1->DATA_BIN_SIZE * settings1->TIMESTEP * (1.E6)); // from Hz to MHz;
                const double freq = k*dfreq;
                // apply RFCM gain
                ApplyRFCM_OutZero(ch, freq, detector, V_tmp, settings1 -> RFCM_OFFSET);
            }

            Vfft_noise_before[k] = V_tmp;

            current_phase = noise_phase[k];

            // use real value array value, extra 1/1.177 to make total power same with "before random_rician".
            Tools::get_random_rician(0., 0., sqrt(2. / M_PI) / 1.177 * V_tmp, current_amplitude, current_phase);

            // vnoise is currently noise spectrum (before fft, unit : V)
            vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
            vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);

            Vfft_noise_after[2*k] = vnoise[2 * k];
            Vfft_noise_after[2*k + 1] = vnoise[2 * k + 1];

            // inverse FFT normalization factor!
            vnoise[2 * k] *= 2. / ((double) settings1 -> DATA_BIN_SIZE);
            vnoise[2 * k + 1] *= 2. / ((double) settings1 -> DATA_BIN_SIZE);

        }

        // now vnoise is time domain waveform
        Tools::realft(vnoise, -1, settings1 -> DATA_BIN_SIZE);

    } else if (settings1 -> NOISE == 1) { // NOISE == 1 : use Rayleigh dist. fits

        Vfft_noise_after.clear(); // remove previous Vfft_noise values
        Vfft_noise_before.clear(); // remove previous Vfft_noise values

        Vfft_noise_after.resize(settings1->DATA_BIN_SIZE);
        Vfft_noise_before.resize(settings1->DATA_BIN_SIZE/2);

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

                    Vfft_noise_before[k] = detector -> GetRayleighFit_databin(ch, k);

                    current_phase = noise_phase[k];

                    V_tmp = detector -> GetRayleighFit_databin(ch, k) * sqrt((double) settings1 -> DATA_BIN_SIZE / (double)(settings1 -> NFOUR / 2.));

                    // use real value array value, extra 1/1.177 to make total power same with "before random_rician".
                    Tools::get_random_rician(0., 0., V_tmp, current_amplitude, current_phase);

                } else {

                    V_tmp = v_noise / sqrt(2.); // copy original flat H_n [V] value, and apply 1/sqrt2 for SURF/TURF divide same as signal

                    if (settings1 -> USE_TESTBED_RFCM_ON == 0) {
                        const double dfreq = 1./(settings1->DATA_BIN_SIZE * settings1->TIMESTEP * (1.E6)); // from Hz to MHz;
                        const double freq = k*dfreq;
                        ApplyFilter_OutZero(freq, detector, V_tmp);
                        ApplyPreamp_OutZero(freq, detector, V_tmp);
                        ApplyFOAM_OutZero(freq, detector, V_tmp);

                    } else if (settings1 -> USE_TESTBED_RFCM_ON == 1) {
                        const double dfreq = 1./(settings1->DATA_BIN_SIZE * settings1->TIMESTEP * (1.E6)); // from Hz to MHz;
                        const double freq = k*dfreq;
                        // apply RFCM gain
                        ApplyRFCM_OutZero(ch, freq, detector, V_tmp, settings1 -> RFCM_OFFSET);
                    }

                    Vfft_noise_before[k] = V_tmp;

                    current_phase = noise_phase[k];

                    // use real value array value, extra 1/1.177 to make total power same with "before random_rician".
                    Tools::get_random_rician(0., 0., sqrt(2. / M_PI) / 1.177 * V_tmp, current_amplitude, current_phase);

                }

                // vnoise is currently noise spectrum (before fft, unit : V)
                vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
                vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);

                Vfft_noise_after[2*k] = vnoise[2 * k];
                Vfft_noise_after[2*k + 1] =  vnoise[2 * k + 1];

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

                Vfft_noise_before[k] = fits_for_this_station[ch][k];
                current_phase = noise_phase[k];
                V_tmp = fits_for_this_station[ch][k];

                // the right normalization factor is N * sqrt(deltaF)
                V_tmp *= double(settings1->DATA_BIN_SIZE);
                V_tmp *= sqrt(this_delta_f);

                // draw a random number from the distribution with this rayleigh fit parameter
                Tools::get_random_rician(0., 0., V_tmp, current_amplitude, current_phase);

                // set the real and imaginary components of the FFT
                vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
                vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);

                // stash those values
                Vfft_noise_after[2*k] = vnoise[2 * k];
                Vfft_noise_after[2*k+1] = vnoise[2 * k + 1];

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

                Vfft_noise_before[k] = fits_for_this_station[ch][k];
                current_phase = noise_phase[k];
                V_tmp = fits_for_this_station[ch][k];

                // the right normalization factor is N * sqrt(deltaF)
                V_tmp *= double(settings1->DATA_BIN_SIZE);
                V_tmp *= sqrt(this_delta_f);

                // draw a random number from the distribution with this rayleigh fit parameter
                Tools::get_random_rician(0., 0., V_tmp, current_amplitude, current_phase);

                // set the real and imaginary components of the FFT
                vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
                vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);

                // stash those values
                Vfft_noise_after[2*k] = vnoise[2 * k];
                Vfft_noise_after[2*k + 1] = vnoise[2 * k + 1];

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
    noise_phase.resize(settings1->DATA_BIN_SIZE/2);

    for (int i=0; i<settings1->DATA_BIN_SIZE/2; i++) {  // noise with DATA_BIN_SIZE bin array
        noise_phase[i] = 2*PI*gRandom->Rndm();  // get random phase for flat thermal noise
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
            int noise_pass_nogo = 1; // index for checking if any same noise_ID is used in different chs.
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
                    // last noise filled in Full_window
                }
                else
                {
                    // when it's not final noise waveform, full noise will fill in Full_window
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
                    // last noise filled in Full_window
                }
                else
                {
                    // when it's not final noise waveform, full noise will fill in Full_window
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
                        // last noise filled in Full_window
                    }
                    else
                    {
                        // when it's not final noise waveform, full noise will fill in Full_window
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
                        // last noise filled in Full_window
                    }
                    else
                    {
                        // when it's not final noise waveform, full noise will fill in Full_window
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

    Tools::Zero(vsignal_forfft,NFOUR/2);

    if (settings1->AVZ_NORM_FACTOR_MODE == 0) { // use previous normalization factors

        double previous_value_e_even=0.;
        double previous_value_e_odd=0.;
        int count_nonzero=0;
        int iprevious=0;
        int ifirstnonzero=-1;
        int ilastnonzero=2000;

        for (int i=0;i<detector->GetFreqBin();i++) {

            // freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
            // but there are only NFOUR/4 different values
            // it's the index among the NFOUR/4 that we're interested in
            int ifour=Tools::Getifreq(detector->GetFreq(i),detector->freq_forfft[0],detector->freq_forfft[NFOUR/2-1],NFOUR/4);

            if (ifour!=-1 && 2*ifour+1<NFOUR/2) {
                count_nonzero++;
                if (ifirstnonzero==-1)
                    ifirstnonzero=ifour;

                // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft
                vsignal_forfft[2*ifour]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP);
                vsignal_forfft[2*ifour+1]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP); // phase is 90 deg.
                // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting

                // how about we interpolate instead of doing a box average

                for (int j=iprevious+1;j<ifour;j++) {
                    vsignal_forfft[2*j]=previous_value_e_even+(vsignal_forfft[2*ifour]-previous_value_e_even)*(double)(j-iprevious)/(double)(ifour-iprevious);
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
        for (int j=0;j<NFOUR/4;j++) {
            vsignal_forfft[2*j]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
            vsignal_forfft[2*j+1]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
        }

        for (int ifour=0;ifour<NFOUR/4;ifour++) {

            vsignal_forfft[2*ifour]*=cos(settings1->PHASE*PI/180.);
            vsignal_forfft[2*ifour+1]*=sin(settings1->PHASE*PI/180.);

        }
    }
    else if (settings1->AVZ_NORM_FACTOR_MODE == 1) { // use new (fixed) normalization factors

        double dF = 1. / ((double)(NFOUR/2) * TIMESTEP); // in Hz

        double dF_org = detector->GetFreq(1) - detector->GetFreq(0); // in Hz

         double dF_factor = 1.;

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

    Tools::Zero(vsignal_forfft,NFOUR/2);

    if (settings1->AVZ_NORM_FACTOR_MODE == 0) { // use previous normalization factors


        double previous_value_e_even=0.;
        double previous_value_e_odd=0.;
        int count_nonzero=0;
        int iprevious=0;
        int ifirstnonzero=-1;
        int ilastnonzero=2000;

        for (int i=0;i<detector->GetFreqBin();i++) {

            // freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
            // but there are only NFOUR/4 different values
            // it's the index among the NFOUR/4 that we're interested in
            int ifour=Tools::Getifreq(detector->GetFreq(i),detector->freq_forfft[0],detector->freq_forfft[NFOUR/2-1],NFOUR/4);

            if (ifour!=-1 && 2*ifour+1<NFOUR/2) {
                count_nonzero++;
                if (ifirstnonzero==-1)
                    ifirstnonzero=ifour;

                // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft
                vsignal_forfft[2*ifour]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP);
                vsignal_forfft[2*ifour+1]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP); // phase is 90 deg.
                // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting

                // how about we interpolate instead of doing a box average

                for (int j=iprevious+1;j<ifour;j++) {
                    vsignal_forfft[2*j]=previous_value_e_even+(vsignal_forfft[2*ifour]-previous_value_e_even)*(double)(j-iprevious)/(double)(ifour-iprevious);
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
        for (int j=0;j<NFOUR/4;j++) {
            vsignal_forfft[2*j]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
            vsignal_forfft[2*j+1]*=1./sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
        }

        for (int ifour=0;ifour<NFOUR/4;ifour++) {
            vsignal_forfft[2*ifour]*=cos(settings1->PHASE*PI/180.);
            vsignal_forfft[2*ifour+1]*=sin(settings1->PHASE*PI/180.);
        }

    }
    else if (settings1->AVZ_NORM_FACTOR_MODE == 1) { // use new (fixed) normalization factors

        double dF = 1. / ((double)(NFOUR/2) * TIMESTEP); // in Hz

        double dF_org = detector->GetFreq(1) - detector->GetFreq(0); // in Hz

         double dF_factor = 1.;

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

    Tools::Zero(vsignal_forfft,NFOUR/2);

    double previous_value_e_even=0.;
    double previous_value_e_odd=0.;
    int count_nonzero=0;
    int iprevious=0;
    int ifirstnonzero=-1;
    int ilastnonzero=2000;

    for (int i=0;i<detector->GetFreqBin();i++) {

	// freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
	// but there are only NFOUR/4 different values
	// it's the index among the NFOUR/4 that we're interested in
	int ifour=Tools::Getifreq(detector->GetFreq(i),detector->freq_forfft[0],detector->freq_forfft[NFOUR/2-1],NFOUR/4);

	if (ifour!=-1 && 2*ifour+1<NFOUR/2) {
	    count_nonzero++;
	    if (ifirstnonzero==-1)
		ifirstnonzero=ifour;

        // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft
	    vsignal_forfft[2*ifour]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP);
	    vsignal_forfft[2*ifour+1]=vsignal_array[i]*2/((double)NFOUR/2)/(TIMESTEP); // phase is 90 deg.
	    // the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting

	    // how about we interpolate instead of doing a box average

	    for (int j=iprevious+1;j<ifour;j++) {
            vsignal_forfft[2*j]=previous_value_e_even+(vsignal_forfft[2*ifour]-previous_value_e_even)*(double)(j-iprevious)/(double)(ifour-iprevious);
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


double Report::FindPeak(vector< double > waveform, int n) {  // same with icemc, trigger-> AntTrigger::FindPeak
    double peak=0.;
    for (int i=0;i<n;i++) {
        if (fabs( (double) waveform[i])>peak)
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

    for (int i = 0; i< detector->stations.size(); i++) {

        for (int j=0; j< detector->stations[i].strings.size(); j++) {

            for (int k=0; k< detector->stations[i].strings[j].antennas.size(); k++) {

                if (stations[i].strings[j].antennas[k].ray_sol_cnt) {

                    bool nonzero_peak = false;

                    for (int n=0; n<stations[i].strings[j].antennas[k].V.size(); n++){
                        if (stations[i].strings[j].antennas[k].PeakV[n].size() > 0 &&
                            stations[i].strings[j].antennas[k].PeakV[n][0] != 0.      ) {
                            nonzero_peak = true;
                            break;
                        }
                    }

                    if (nonzero_peak){
                        stations[i].strings[j].antennas[k].Rank.push_back(current+1);  // set Rank as 1 (for now)
                    }
                    else{
                        stations[i].strings[j].antennas[k].Rank.push_back(0);  // rank 0, PeakV is 0, non-countable rank.
                    }

                }   // if ray_sol_cnt != 0
            }
        }
    }

    check = 1;
    pre_maxpeak = 1.E5; // unreasonably big value which real PeakV will never reach
    while (check!=0) {
        check=0;
        maxpeak = 0.;
        for (int i = 0; i< detector->stations.size(); i++) {

            for (int j=0; j< detector->stations[i].strings.size(); j++) {

                for (int k=0; k< detector->stations[i].strings[j].antennas.size(); k++) {

                    //lets start with 1st ray_sol only
                    if (stations[i].strings[j].antennas[k].ray_sol_cnt) {
                        if (stations[i].strings[j].antennas[k].Rank[0] != 0) {  // there is non-zero value! and ray_sol_cnt also non-zero!
                            if ( stations[i].strings[j].antennas[k].PeakV[0].size()!=0 && stations[i].strings[j].antennas[k].PeakV[0][0] < pre_maxpeak) {

                                if (maxpeak < stations[i].strings[j].antennas[k].PeakV[0][0] ) {
                                    maxpeak = stations[i].strings[j].antennas[k].PeakV[0][0];

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

}


int Report::GetChannelNum8_LowAnt(int string_num, int antenna_num) {

    int outputch = 16;

    // just give ch numbers 1-8 for antenna 0 - 1 while higher ch numbers for antenna 2 - 3
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

                for (int n=0; n<station.strings[s].antennas[a].V.size(); n++){
                    // Loop over interactions

                    for (int m=0; m<station.strings[s].antennas[a].V[n].size(); m++){
                        // Loop over ray solutions from this interaction

                        // Check if the max (abs) value is greater than
                        //   the signals from the other rays
                        double ray_max_signal = Tools::getMaxAbsoluteMagnitude(
                            station.strings[s].antennas[a].V[n][m]
                        );
                        if ( ray_max_signal > max_signal ){
                            max_signal = ray_max_signal;
                        }

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


bool Report::isTrigger(
    vector <double> waveform, int *brightest_event, double randNum,
    Antenna_r *antenna, int ch_ID, Settings *settings, Trigger *trigger) {

    double avgSnr = 0.;
    if(settings->TRIG_ANALYSIS_MODE == 2) { // Noise only triggers
        avgSnr = pa_force_trigger_snr;
    }
    else {
        // Estimate average SNR in topmost vpol
        if( waveform.size() > 0 ) {
            //   it as the noise WF to get_SNR() (since the RMS of an
            //   array with one element is the absolute value of that element)
            vector <double> tmp_noise_RMS;
            double ant_noise_voltage_RMS = trigger->GetAntNoise_voltageRMS(ch_ID, settings);
            tmp_noise_RMS.push_back( ant_noise_voltage_RMS );

            // Calculate SNR in this antenna using trace w/o noise since
            //   this SNR is used with a signal-only SNR efficiency curve
            //   to estimate trigger likelihood
            avgSnr = get_SNR(
                waveform,
                tmp_noise_RMS);
        }
        else {
            avgSnr = 0.0;
        }
    }

    // Scale SNR with respect to the antenna's viewing angle of the signal
    double viewangle = antenna->view_ang[brightest_event[0]][brightest_event[1]];
    viewangle = viewangle * 180.0/PI - 90.0;
    double snr_50 = interpolate(
        trigger->angle_PA, trigger->aSNR_PA, // x and y coordinates of curve to interpolate
        viewangle, // x value to interpolate y value for
        (*(&trigger->angle_PA+1) - trigger->angle_PA) - 1 // len(ang_data) - 1
    );
    avgSnr = avgSnr*2.0/snr_50;

    // Estimate the PA signal efficiency of for this SNR from curve of efficiency vs data
    double eff = interpolate(
        trigger->snr_PA, trigger->eff_PA, // x and y coordinates of curve to interpolate
        avgSnr, // x value to interpolate y value for
        (*(&trigger->snr_PA+1) - trigger->snr_PA) - 1 // len(snr_PA) - 1
    );

    if (avgSnr > 0.5){
        if (eff >= 1.0){
            cout<<" efficiency : "<<eff<<"  avg SNR : "<<avgSnr<<endl;
            return true;
        }
        randNum = gRandom->Rndm();
        if (randNum < eff){
            cout<<" efficiency : "<<eff<<"  avg SNR : "<<avgSnr<<"  random number : "<<randNum<<endl;
            return true;
        }
    }

    return false;

}

int Report::get_PA_trigger_bin(
    int ch_ID, Antenna_r *antenna, vector <double> waveform, double randNum,
    Settings *settings, Trigger *trigger
){
    // From the given waveform, find and return the starting index of the
    //   first 10.7 ns time window with an SNR that triggers the PA.
    // Return -1 if no 10.7 ns time window triggers the PA

    // Determine number of bins corresponding to a 10.7 ns integration window
    // 10.7 ns from page 12 of the PA design and instrumentation paper: https://arxiv.org/abs/1809.04573
    int trigger_window_bins = (10.7E-9) / settings->TIMESTEP;

    // Create a 10.7 ns window of waveform to analyze for a trigger-worthy signal
    //   and move this window along the full waveform, returning the bin
    //   where a strong enough signal exists, if it does.
    const int waveform_length = (int)waveform.size();
    if(waveform_length < trigger_window_bins) {
        throw runtime_error("Waveform length is shorter than 10.7 ns!");
    }

    for (int bin=0; bin < waveform_length-trigger_window_bins + 1; bin++) {

        // Get this 10.7ns snapshot of waveform and find the bin with greatest absolute value
        vector <double> waveform_window(waveform.begin()+bin, waveform.begin()+bin+trigger_window_bins);
        int bin_max_value =
            (TMath::MaxElement(trigger_window_bins, &waveform_window[0]) > -1* TMath::MinElement(trigger_window_bins, &waveform_window[0]))?
            TMath::LocMax(trigger_window_bins, &waveform_window[0])
            : TMath::LocMin(trigger_window_bins, &waveform_window[0]); // Bin (in waveform_window) with the greatest absolute value

        // Identify the ray solution that arrived closest to bin with the
        //   greatest absolute value in this 10.7ns window.
        int Likely_Sol[2] = {0, 0}; // Initialize to direct ray
        int mindBin = 1.e9; // Minimum difference between this bin and a ray solutions signal_bin.
                            // Initialized unphysically large.
        int dBin = 0; // Difference between this bin and some ray solution's signal bin.
        for (int n=0; n<antenna->SignalBin.size(); n++) { // loop over interactions
            for (int m=0; m<antenna->SignalBin[n].size(); m++){ // loop over raysol numbers
                if ( antenna->SignalExt[n][m] ) {
                    dBin = abs( antenna->SignalBin[n][m] - bin_max_value );
                    if ( dBin < mindBin ) {
                        Likely_Sol[0] = n;
                        Likely_Sol[1] = m;
                        mindBin = dBin;
                    }
                }
            }
        }

        // Scale SNR according to full phased array angular response
        double avgSnr = 0.;
        if(settings->TRIG_ANALYSIS_MODE == 2) { // Noise only triggers
            avgSnr = pa_force_trigger_snr;
        }
        else {
            // Estimate average SNR in topmost vpol
            if( waveform.size() > 0 ) {
                //   it as the noise WF to get_SNR() (since the RMS of an
                //   array with one element is the absolute value of that element)
                vector <double> tmp_noise_RMS;
                double ant_noise_voltage_RMS = trigger->GetAntNoise_voltageRMS(ch_ID, settings);
                tmp_noise_RMS.push_back( ant_noise_voltage_RMS );

                // Calculate SNR in this antenna using trace w/o noise since
                //   this SNR is used with a signal-only SNR efficiency curve
                //   to estimate trigger likelihood
                avgSnr = get_SNR(waveform_window, tmp_noise_RMS);
            }
            else {
                avgSnr = 0.0;
            }
        }

        // Scale SNR with respect to the antenna's viewing angle of the signal
        const double viewangle = antenna->theta_rec[Likely_Sol[0]][Likely_Sol[1]] * 180.0/PI - 90.0 ;
        double snr_50 = interpolate(
            trigger->angle_PA, trigger->aSNR_PA, // x and y coordinates of curve to interpolate
            viewangle, // x value to interpolate y value for
            (*(&trigger->angle_PA+1) - trigger->angle_PA) - 1 // len(ang_data) - 1
        );
        avgSnr = avgSnr*2.0/snr_50;

        // Estimate the PA signal efficiency of for this SNR from curve of efficiency vs data
        double eff = interpolate(
            trigger->snr_PA, trigger->eff_PA, // x and y coordinates of curve to interpolate
            avgSnr, // x value to interpolate y value for
            (*(&trigger->snr_PA+1) - trigger->snr_PA) - 1 // len(snr_PA) - 1
        );

        // If the PA triggers on this SNR, return the starting bin of this window
        if (avgSnr > 0.5){
            if (eff >= 1.0){
                cout<<" efficiency : "<<eff<<"  avg SNR : "<<avgSnr<<endl;
                return bin;
            }
            if (randNum < eff){
                cout<<" efficiency : "<<eff<<"  avg SNR : "<<avgSnr<<endl;
                return bin;
            }
        }

    }

    // If no windows triggered, return -1
    return -1;

}

void Report::checkPATrigger(
    int i, Detector *detector, Event *event, int evt, Trigger *trigger, Settings *settings1,
    int trig_search_init, int max_total_bin
){
    // Calculates max SNR in topmost PA Vpol
    //   and multplies it by viewing angle factors.
    // Then determines signal efficiency by interpolating oberved SNR from
    //   efficiency vs SNR data.
    // Triggers if signal efficiency is above certain treshold.
    // If triggered, saves relevant information.

    // Identify the antenna used to estimate PA trigger
    int trigger_antenna_string_i = 0;
    int trigger_antenna_i = 8;
    int trigger_antenna_number = 8;
    Antenna_r *trigger_antenna = &stations[i].strings[trigger_antenna_string_i].antennas[trigger_antenna_i];
    int trigger_ch_ID = GetChNumFromArbChID(detector, trigger_antenna_number, i, settings1) - 1;

    // If the antenna we use to trigger the PA doesn't have any ray solutions,
    //   do not perform the trigger check.
    if (trigger_antenna->ray_sol_cnt == 0){
        return;
    }

    // For phased array, waveform length is 680 ns, but
    // for PA trigger check only 20 ns around the signal bin.
    // This is to avoid getting the second ray
    // KAH and ARB are not sure where the 1200 number comes from
    // !!!HEY IF YOU ARE FIXING THIS!!! please update the corresponding check in Settings::CheckCompatibilitiesSettings
    int BINSIZE = 1200/(settings1->TIMESTEP*1.e9);  // Number of bins (aka datapoints) of data to save

    int waveformLength = settings1->WAVEFORM_LENGTH;
    int waveformCenter = settings1->WAVEFORM_CENTER;

    // Create the signal-only waveform that will be used for triggering
    vector<double>* V_convolved = &trigger_antenna->V_convolved;
    const int convolved_len = (int)V_convolved->size();
    vector<double> trigger_waveform(V_convolved->begin(), V_convolved->begin()+convolved_len-1);

    // Find and log the event and ray with the most signal
    int brightest_event[2];
    trigger_antenna->Get_Brightest_Interaction(&brightest_event);

    // Create an object to save the random number generated in isTrigger()
    //   for use in get_PA_trigger_bin().
    double random_number = 1.01;

    if( isTrigger(
        trigger_waveform, brightest_event, random_number,
        trigger_antenna, trigger_ch_ID, settings1, trigger
    )){
        cout<<endl<<"PA trigger ~~~  Event Number : "<<evt<<endl;

        int last_trig_bin = get_PA_trigger_bin(
            trigger_ch_ID, trigger_antenna, trigger_waveform, random_number,
            settings1, trigger
        );
        if (last_trig_bin == -1) {
            cerr<<"Waveform triggered but we cannot find a trigger bin."<<endl;
            last_trig_bin = trigger_antenna->SignalBin[brightest_event[0]][brightest_event[1]];
        }
        int my_ch_id = 0;
        stations[i].Global_Pass = last_trig_bin;
        trigger_antenna->Trig_Pass = last_trig_bin;
        for (size_t str = 0; str < detector->stations[i].strings.size(); str++) {
            for (size_t ant = 0; ant < detector->stations[i].strings[str].antennas.size(); ant++) {
                double peakvalue = 0;
                for (int bin=0; bin<waveformLength; bin++) {

                    int bin_value = last_trig_bin + waveformCenter - waveformLength/2 + bin;
                    stations[i].strings[str].antennas[ant].V_mimic.push_back(trigger->Full_window_V[my_ch_id][bin_value]);// save in V (KAH)
                    stations[i].strings[str].antennas[ant].time.push_back( bin_value );
                    stations[i].strings[str].antennas[ant].time_mimic.push_back( ( bin) * settings1->TIMESTEP*1.e9 );// save in ns
                    if (TMath::Abs(trigger->Full_window_V[ant][bin_value]) > peakvalue) {
                        peakvalue = TMath::Abs(trigger->Full_window_V[my_ch_id][bin_value]);
                    }

                }//end bin
                my_ch_id ++;
            }//end ant
        }//end detector

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
                            if(Pthresh_value[trig_j]<stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.back())
                                stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.back() = Pthresh_value[trig_j];
                        }
                    }// trig scan mode > 2

                    SCTR_cluster_bit[trig_j]=1;

                }// if local trigger
                else SCTR_cluster_bit[trig_j]=0;// if no local trigger, set zero to start a new cluster at next local trigger

                // and how many bins scanned
                stations[i].strings[string_i].antennas[antenna_i].TotalBinsScannedPerChannel++;

                // fill the buffers (if any changes occur mark check_TDR_configuration as non-zero)
                if(trig_i<trig_search_init+trig_window_bin)
                    check_TDR_configuration+=buffer[trig_j]->fill(Pthresh_value[trig_j]);
                else
                    check_TDR_configuration += buffer[trig_j]->add(Pthresh_value[trig_j]);

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

                        // mark the bin on which we triggered...
                        if(buffer[trig_j]->addToNPass>0)
                            stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = trig_i-buffer[trig_j]->numBinsToOldestTrigger();
                        else
                            stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = 0.;

                    }// for trig_j
                    bin_to_save_on = trig_i;
                }// first trigger

            } // end if N_Pass_Vpol > some value

            // if there's a trigger and anything changes in the buffers, restock the TDR arrays
            if(settings1->TRIG_SCAN_MODE>1&&check_TDR_configuration&&window_pass_bit){

                for(int trig_j=0;trig_j<numChan;trig_j++)
                    TDR_all_sorted_temp[trig_j]=0;
                for(int trig_j=0;trig_j<numChanVpol;trig_j++)
                    TDR_Vpol_sorted_temp[trig_j]=0;
                for(int trig_j=0;trig_j<numChanHpol;trig_j++)
                    TDR_Hpol_sorted_temp[trig_j]=0;

                // Changes TDR sorting and buffer[best_chan] for Vpol only:
                for(int ii=0;ii<numChanVpol; ii++){// find the best channel's TDR and store them.

                    double best_thresh=0;
                    int best_chan=0;

                    // Pick a new best_thresh and best_chan
                    for(int trig_j=0;trig_j<numChan;trig_j++){

                        int string_i = detector->getStringfromArbAntID( i, trig_j);
                        int antenna_i = detector->getAntennafromArbAntID( i, trig_j);
                        if(detector->stations[i].strings[string_i].antennas[antenna_i].type == 0
                            && buffer[trig_j]->temp_value<best_thresh)
                        {
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
                        if(detector->stations[i].strings[string_i].antennas[antenna_i].type==1
                            && buffer[trig_j]->temp_value<best_thresh)
                        {
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

        if ( bin_to_save_on == -1 )
            bin_to_save_on = trigger_antenna->SignalBin[brightest_event[0]][brightest_event[1]];

        // Do what saveTriggeredEvent() does
        for(int trig_j=0; trig_j<numChan;trig_j++){

            int string_i = detector->getStringfromArbAntID( i, trig_j);
            int antenna_i = detector->getAntennafromArbAntID( i, trig_j);

            // Determine which ray and interaciton triggered triggered
            //   the station based on signal and trigger bins
            stations[i].strings[string_i].antennas[antenna_i].Find_Likely_Sol(); // no likely init

            // set global_trig_bin values
            if (settings1->V_MIMIC_MODE == 0) { // Global passed bin is the center of the window
                stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = waveformLength/2 + waveformCenter ;
            }
            else if (settings1->V_MIMIC_MODE == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
                stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (
                (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin)
                + waveformLength/2 + waveformCenter );
            }
            else if (settings1->V_MIMIC_MODE == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
                stations[i].strings[string_i].antennas[antenna_i].global_trig_bin = (
                (detector->params.TestBed_Ch_delay_bin[trig_j] - detector->params.TestBed_BH_Mean_delay_bin
                        + detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
                + waveformLength/2 + waveformCenter );
            }

            stations[i].total_trig_search_bin = stations[i].Global_Pass + trig_window_bin - trig_search_init;

        } // end for trig_j in numchans

    } // end if isTrigger

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
