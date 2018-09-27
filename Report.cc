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

#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "TRandom3.h"
#include "TMath.h"

#include <cstdlib>

ClassImp (Report);
ClassImp (Antenna_r);
ClassImp (Surface_antenna_r);
ClassImp (String_r);
ClassImp (Station_r);

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
	connect_signals.clear(); // if ray_sol time delay is small enough to connect each other
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

	for (int i = 0; i < detector->params.number_of_stations; i++) {
		// vector stations
		stations.push_back(tmp_station);
		for (int j = 0; j < detector->stations[i].strings.size(); j++) {
			// vector strings
			stations[i].strings.push_back(tmp_string);
			for (int k = 0;
					k < detector->stations[i].strings[j].antennas.size(); k++) {
				// vector antennas
				stations[i].strings[j].antennas.push_back(tmp_antenna);
			}
		}
		for (int j = 0; j < detector->stations[i].surfaces.size(); j++) {
			// vector surface antennas
			stations[i].surfaces.push_back(tmp_surface);
		}

		if (settings1->TRIG_SCAN_MODE > 0) {    // scan Pthresh mode

			int numChan = 0;
			int numChanVpol = 0;
			int numChanHpol = 0;

			for (int j = 0; j < detector->stations[i].strings.size(); j++) {

				for (int k = 0;
						k < detector->stations[i].strings[j].antennas.size();
						k++) {

					int string_i = detector->getStringfromArbAntID(i, numChan);
					int antenna_i = detector->getAntennafromArbAntID(i,
							numChan);
					if (detector->stations[i].strings[string_i].antennas[antenna_i].type
							== 0)
						numChanVpol++;
					if (detector->stations[i].strings[string_i].antennas[antenna_i].type
							== 1)
						numChanHpol++;
					numChan++;
				}    // for k

			}    // for j

			stations[i].TDR_all.clear();
			for (int ch = 0; ch < numChan; ch++)
				stations[i].TDR_all.push_back(0);
			stations[i].TDR_all_sorted.clear();
			if (settings1->TRIG_MODE == 0)
				for (int ch = 0; ch < numChan; ch++)
					stations[i].TDR_all_sorted.push_back(0);

			stations[i].TDR_Vpol_sorted.clear();
			if (settings1->TRIG_MODE == 1)
				for (int ch = 0; ch < numChanVpol; ch++)
					stations[i].TDR_Vpol_sorted.push_back(0);

			stations[i].TDR_Hpol_sorted.clear();
			if (settings1->TRIG_MODE == 1)
				for (int ch = 0; ch < numChanHpol; ch++)
					stations[i].TDR_Hpol_sorted.push_back(0);

		}    // if TRIG_SCAN_MODE 

	}    // for i (number of stations)

}

void Antenna_r::clear() { // if any vector variable added in Antenna_r, need to be added here!

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
	//VHz_antfactor.clear();
	//VHz_filter.clear();
	Vfft.clear();
	Vfft_noise.clear();

	time.clear();
	time_mimic.clear();
	V_mimic.clear();
	Ax.clear();
	Ay.clear();
	Az.clear();

	V.clear();
	//V_noise.clear();
	//V_total.clear();
	//V_total_diode.clear();
	//V_total_timedelay.clear();

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
	SignalExt.clear();

	SCT_threshold_pass.clear();

}

void Antenna_r::clear_useless(Settings *settings1) { // to reduce the size of output AraOut.root, remove some information

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
		//V_noise.clear();
		//V_total.clear();
		//V_total_diode.clear();
		//V_total_timedelay.clear();

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

	} else if (settings1->DATA_SAVE_MODE == 2) {

		Heff.clear();
		//VHz_antfactor.clear();
		//VHz_filter.clear();
		Vfft.clear();
		Vfft_noise.clear();

		Ax.clear();
		Ay.clear();
		Az.clear();

		V.clear();
		//V_noise.clear();
		//V_total.clear();
		//V_total_diode.clear();
		//V_total_timedelay.clear();

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
		V_mimic.clear();

		// additional for before ant waveform
		//Vm_wo_antfactor.clear();
		Vm_zoom.clear();
		Vm_zoom_T.clear();

	}

}

void Report::clear_useless(Settings *settings1) { // to reduce the size of output AraOut.root, remove some information

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

//void Report::Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger) {
//void Report::Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger, UsefulIcrrStationEvent *theUsefulEvent) {
//void Report::Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger) {
//void Report::Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger, UsefulIcrrStationEvent *theUsefulEvent) {
//void Report::Connect_Interaction_Detector (Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger, UsefulIcrrStationEvent *theUsefulEvent, int evt) {
void Report::Connect_Interaction_Detector(Event *event, Detector *detector, RaySolver *raysolver, Signal *signal, IceModel *icemodel, Settings *settings1, Trigger *trigger, int evt) {

	int ray_sol_cnt;
	double viewangle;
	Position launch_vector; // direction of ray at the source
	Position receive_vector;    // direction of ray at the target antenna
	Vector n_trg_pokey;         // unit pokey vector at the target
	Vector n_trg_slappy;        // unit slappy vector at the target
	vector < vector<double> > ray_output;

	double vmmhz1m_tmp, vmmhz1m_sum, vmmhz1m_em; // currently not using vmmhz1m_em
	Position Pol_vector;                    // polarization vector at the source
	double mag; // magnification factor. it can vary in case of plane / spherical wave
	double fresnel;                                 // fresnel factor
	double Pol_factor;                              // polarization factor
	double tmp; // for non use return values

	double freq_tmp, heff, antenna_theta, antenna_phi; // values needed for apply antenna gain factor and prepare fft, trigger
	double volts_forfft[settings1->NFOUR / 2];       // array for fft
	//double T_forfft[settings1->NFOUR/2];       // array for fft
	double dT_forfft;
	//double V_total_forconvlv[settings1->DATA_BIN_SIZE];   //
	//
	double volts_forint[settings1->NFOUR / 2];       // array for interpolation
	double T_forint[settings1->NFOUR / 2];       // array for interpolation

	double dF_NFOUR = 1.
			/ ((double) (settings1->NFOUR / 2) * settings1->TIMESTEP); // in Hz

	int waveformLength = settings1->WAVEFORM_LENGTH;
	int waveformCenter = settings1->WAVEFORM_CENTER;

	//double dF_outbin;
	double dF_Nnew;

	double heff_lastbin;
	double freq_lastbin;

	int check_toomuch_Tdelay;   // return value from MixSignalNoise_Tdelay

	double min_arrival_time_tmp; // min arrival time between all antennas, raysolves
	double max_arrival_time_tmp; // max arrival time between all antennas, raysolves
	double max_PeakV_tmp;       // max PeakV of all antennas in the station

	int trig_window_bin = (int) (settings1->TRIG_WINDOW / settings1->TIMESTEP); // coincidence window bin for trigger

	RandomTshift = gRandom->Rndm();

	init_T = settings1->TIMESTEP * -1.e9
			* ((double) settings1->NFOUR / 4 + RandomTshift); // locate zero time at the middle and give random time shift

	for (int n = 0; n < settings1->NFOUR / 2; n++) {
		T_forint[n] = init_T + (double) n * settings1->TIMESTEP * 1.e9; // in ns
	}

	// decide whether debug mode or not
	int debugmode = 0;
	if (settings1->DEBUG_MODE_ON == 1 && evt < settings1->DEBUG_SKIP_EVT)
		debugmode = 1;
	//else if ( settings1->DEBUG_MODE_ON == 1 && evt >= settings1->DEBUG_SKIP_EVT ) cout<<"After DEBUG_SKIP_EVT"<<endl;
	else if (settings1->DEBUG_MODE_ON == 1 && evt >= settings1->DEBUG_SKIP_EVT)
		cout << evt << " " << endl;
	// skip most of computation intensive processes if debugmode == 1

	int N_pass; // number of trigger passed channels (antennas)
	int N_pass_V; // number of trigger passed channels (Vpol antennas)
	int N_pass_H; // number of trigger passed channels (Hpol antennas)

	for (int i = 0; i < detector->params.number_of_stations; i++) {

		min_arrival_time_tmp = 10.; // first min_arrival_time is unreasonably big value
		max_arrival_time_tmp = 0.; // first max_arrival_time is unreasonably small value
		max_PeakV_tmp = 0.;  // first max_PeakV_tmp is 0.

		stations[i].Total_ray_sol = 0; // initial Total_ray_sol value

		for (int j = 0; j < detector->stations[i].strings.size(); j++) {

			for (int k = 0; k < detector->stations[i].strings[j].antennas.size(); k++) {
				//    cout << i << " : " << j << " : " << k << endl;

				stations[i].strings[j].antennas[k].clear(); // clear data in antenna which stored in previous event

				// run ray solver, see if solution exist
				// if not, skip (set something like Sol_No = 0;
				// if solution exist, calculate view angle and calculate TaperVmMHz

				// added one more condition to run raysolver ( direct distance is less than 3km )
				//

				//    cout << event->Nu_Interaction[0].posnu.GetX() << " : " << event->Nu_Interaction[0].posnu.GetY() << " : " << event->Nu_Interaction[0].posnu.GetZ() << endl;
				//    cout << event->Nu_Interaction[0].pickposnu << " : " << event->Nu_Interaction[0].posnu.Distance( detector->stations[i].strings[j].antennas[k] ) << " : " << settings1->RAYSOL_RANGE << endl; 
				if (event->Nu_Interaction[0].pickposnu
						&& event->Nu_Interaction[0].posnu.Distance(
								detector->stations[i].strings[j].antennas[k])
								<= settings1->RAYSOL_RANGE) { // if posnu is selected inside the antarctic ic:"<<viewangle<<" th_em:"<<d_theta_em[l]<<" th_had:"<<d_theta_had[l]<<" emfrac:"<<emfrac<<" hadfrac:"<<hadfrac<<" vmmhz1m:"<<vmmhz1m[l]<<endl;e
						//if (event->Nu_Interaction[0].pickposnu && event->Nu_Interaction[0].posnu.Distance( detector->stations[i].strings[j].antennas[k] ) <= settings1->RAYSOL_RANGE && debugmode == 0 ) {    // if posnu is selected inside the antarctic ic:"<<viewangle<<" th_em:"<<d_theta_em[l]<<" th_had:"<<d_theta_had[l]<<" emfrac:"<<emfrac<<" hadfrac:"<<hadfrac<<" vmmhz1m:"<<vmmhz1m[l]<<endl;e
						//cout << i << " : " << j << " : " << k << endl;

					//raysolver->Solve_Ray(event->Nu_Interaction[0].posnu, detector->stations[i].strings[j].antennas[k], icemodel, ray_output, settings1);   // solve ray between source and antenna
					RayStep.clear(); // remove previous values
					raysolver->Solve_Ray(event->Nu_Interaction[0].posnu,
							detector->stations[i].strings[j].antennas[k],
							icemodel, ray_output, settings1, RayStep); // solve ray between source and antenna

					ray_sol_cnt = 0;

					if (raysolver->solution_toggle) { // if there are solution from raysolver
						//if (raysolver->solution_toggle && debugmode==0 ) {  // if there are solution from raysolver

						while (ray_sol_cnt < ray_output[0].size()) { // for number of soultions (could be 1 or 2)

							stations[i].strings[j].antennas[k].arrival_time.push_back(
									ray_output[4][ray_sol_cnt]);

							// get ice attenuation factor
							//
							double IceAttenFactor = 1.;
							if (settings1->USE_ARA_ICEATTENU == 1) { // use new ARA measured ice attenuation values

								double dx, dz, dl;
								for (int steps = 1;
										steps
												< (int) RayStep[ray_sol_cnt][0].size();
										steps++) {

									dx = RayStep[ray_sol_cnt][0][steps - 1]
											- RayStep[ray_sol_cnt][0][steps];
									dz = RayStep[ray_sol_cnt][1][steps - 1]
											- RayStep[ray_sol_cnt][1][steps];
									dl = sqrt((dx * dx) + (dz * dz));

									// use new ice model
									IceAttenFactor *=
											exp(
													-dl
															/ icemodel->GetARAIceAttenuLength(
																	-RayStep[ray_sol_cnt][1][steps]));
								}
								//cout<<"new iceattenfactor : "<<IceAttenFactor<<", old way : "<<exp(-ray_output[0][ray_sol_cnt]/icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0))<<endl;
							} else if (settings1->USE_ARA_ICEATTENU == 0) { // use old method

								IceAttenFactor =
										exp(
												-ray_output[0][ray_sol_cnt]
														/ icemodel->EffectiveAttenuationLength(
																settings1,
																event->Nu_Interaction[0].posnu,
																0));
							}

							if (debugmode == 0) {

								// set viewangle, launch_vector, receive vectors
								viewangle = ray_output[1][ray_sol_cnt];
								GetParameters(
										event->Nu_Interaction[0].posnu, // posnu
										detector->stations[i].strings[j].antennas[k], // trg antenna
										event->Nu_Interaction[0].nnu,     // nnu
										viewangle, // inputs launch_angle, returns viewangle
										ray_output[2][ray_sol_cnt], // receive_angle
										launch_vector, receive_vector,
										n_trg_slappy, n_trg_pokey);

								// check viewangle that if ray in near Cherenkov cone
								//
								if (viewangle * DEGRAD > 55.
										&& viewangle * DEGRAD < 57.) { // if viewangle is 56 deg +- 1 deg
										//cout<<"near cone! view angle : "<<viewangle*DEGRAD<<"  station["<<i<<"].string["<<j<<"].antenna["<<k<<"] with  ray_sol_cnt : "<<ray_sol_cnt<<endl;
								}

								// store information to report
								stations[i].strings[j].antennas[k].view_ang.push_back(
										viewangle);
								stations[i].strings[j].antennas[k].launch_ang.push_back(
										ray_output[1][ray_sol_cnt]);
								stations[i].strings[j].antennas[k].rec_ang.push_back(
										ray_output[2][ray_sol_cnt]);
								stations[i].strings[j].antennas[k].Dist.push_back(
										ray_output[0][ray_sol_cnt]);
								stations[i].strings[j].antennas[k].L_att.push_back(
										icemodel->EffectiveAttenuationLength(
												settings1,
												event->Nu_Interaction[0].posnu,
												0));
								//stations[i].strings[j].antennas[k].arrival_time.push_back(ray_output[4][ray_sol_cnt]);
								stations[i].strings[j].antennas[k].reflect_ang.push_back(
										ray_output[3][ray_sol_cnt]);

								stations[i].strings[j].antennas[k].vmmhz.resize(
										ray_sol_cnt + 1);

								stations[i].strings[j].antennas[k].Heff.resize(
										ray_sol_cnt + 1);

								stations[i].strings[j].antennas[k].Vm_zoom.resize(
										ray_sol_cnt + 1);
								stations[i].strings[j].antennas[k].Vm_zoom_T.resize(
										ray_sol_cnt + 1);

								//stations[i].strings[j].antennas[k].Vm_wo_antfactor.resize(ray_sol_cnt+1);

								//stations[i].strings[j].antennas[k].VHz_antfactor.resize(ray_sol_cnt+1);
								//stations[i].strings[j].antennas[k].VHz_filter.resize(ray_sol_cnt+1);
								stations[i].strings[j].antennas[k].Vfft.resize(
										ray_sol_cnt + 1);
								stations[i].strings[j].antennas[k].Vfft_noise.resize(
										ray_sol_cnt + 1);

								stations[i].strings[j].antennas[k].V.resize(
										ray_sol_cnt + 1);
								//stations[i].strings[j].antennas[k].V_noise.resize(ray_sol_cnt+1);
								//stations[i].strings[j].antennas[k].V_total.resize(ray_sol_cnt+1);
								//stations[i].strings[j].antennas[k].V_total_diode.resize(ray_sol_cnt+1);
								//stations[i].strings[j].antennas[k].V_total_timedelay.resize(ray_sol_cnt+1);
								//stations[i].strings[j].antennas[k].time.resize(ray_sol_cnt+1);

								stations[i].strings[j].antennas[k].SignalExt.resize(
										ray_sol_cnt + 1);

								// calculate the polarization vector at the source
								Pol_vector = GetPolarization(
										event->Nu_Interaction[0].nnu,
										launch_vector);

								icemodel->GetFresnel(
										ray_output[1][ray_sol_cnt], // launch_angle
										ray_output[2][ray_sol_cnt], // rec_angle
										ray_output[3][ray_sol_cnt], // reflect_angle
										event->Nu_Interaction[0].posnu,
										launch_vector, receive_vector,
										settings1, fresnel, mag, Pol_vector); // input src Pol and return Pol at trg

								if (ray_output[3][ray_sol_cnt] < PI / 2.) { // when not reflected at the surface, angle = 100
									stations[i].strings[j].antennas[k].reflection.push_back(
											1);  // say this is reflected ray
								} else {
									stations[i].strings[j].antennas[k].reflection.push_back(
											0); // say this is not reflected ray
								}

								stations[i].strings[j].antennas[k].Pol_vector.push_back(
										Pol_vector); // this Pol_vector is for the target antenna
								stations[i].strings[j].antennas[k].Mag.push_back(
										mag); // magnification factor
								stations[i].strings[j].antennas[k].Fresnel.push_back(
										fresnel); // Fresnel factor

								vmmhz1m_sum = 0;

								GetAngleAnt(receive_vector,
										detector->stations[i].strings[j].antennas[k],
										antenna_theta, antenna_phi); // get theta, phi for signal ray arrived at antenna
								//cout<<"antenna theta : "<<antenna_theta<<"  phi : "<<antenna_phi<<endl;  

								// old freq domain signal mode (AVZ model)
								if (settings1->SIMULATION_MODE == 0) {

									// initially give raysol has actual signal
									stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] =
											1;

									double vmmhz_filter[(int) (detector->GetFreqBin())];

									for (int l = 0; l < detector->GetFreqBin();
											l++) { // for detector freq bin numbers

										//cout<<"TaperVmMHz inputs VA:"<<viewangle<<" th_em:"<<event->Nu_Interaction[0].d_theta_em[l]<<" th_had:"<<event->Nu_Interaction[0].d_theta_had[l]<<" emfrac:"<<event->Nu_Interaction[0].emfrac<<" hadfrac:"<<event->Nu_Interaction[0].hadfrac<<" vmmhz1m:"<<event->Nu_Interaction[0].vmmhz1m[l]<<endl;
										//                                       switch (event->IsCalpulser){
										//                                           case 0:
										if (event->IsCalpulser > 0) {
											vmmhz1m_tmp =
													event->Nu_Interaction[0].vmmhz1m[l]
															* settings1->CALPUL_AMP;
											//vmmhz1m_tmp = event->Nu_Interaction[0].vmmhz1m[l];// calpulser -> let's use slight offset from cone
											//signal->TaperVmMHz( settings1->CALPUL_OFFCONE_ANGLE*RADDEG, event->Nu_Interaction[0].d_theta_em[l], event->Nu_Interaction[0].d_theta_had[l], event->Nu_Interaction[0].emfrac, event->Nu_Interaction[0].hadfrac, vmmhz1m_tmp, vmmhz1m_em);
										} else {
											vmmhz1m_tmp =
													event->Nu_Interaction[0].vmmhz1m[l];
											signal->TaperVmMHz(viewangle,
													event->Nu_Interaction[0].d_theta_em[l],
													event->Nu_Interaction[0].d_theta_had[l],
													event->Nu_Interaction[0].emfrac,
													event->Nu_Interaction[0].hadfrac,
													vmmhz1m_tmp, vmmhz1m_em);
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

										//                                       signal->TaperVmMHz( viewangle, event->Nu_Interaction[0].d_theta_em[l], event->Nu_Interaction[0].d_theta_had[l], event->Nu_Interaction[0].emfrac, event->Nu_Interaction[0].hadfrac, vmmhz1m_tmp, vmmhz1m_em);
										//cout<<"TaperVmMHz (1m at view angle) at "<<l<<"th bin : "<<vmmhz1m_tmp<<endl;

										// multiply all factors for traveling ice
										//vmmhz1m_tmp = vmmhz1m_tmp / ray_output[0][ray_sol_cnt] * exp(-ray_output[0][ray_sol_cnt]/icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0)) * mag * fresnel;  // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
										vmmhz1m_tmp = vmmhz1m_tmp
												/ ray_output[0][ray_sol_cnt]
												* IceAttenFactor * mag
												* fresnel; // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
										//cout<<"AttenLength : "<<icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0)<<endl;

										vmmhz1m_sum += vmmhz1m_tmp;

										stations[i].strings[j].antennas[k].vmmhz[ray_sol_cnt].push_back(
												vmmhz1m_tmp);

										freq_tmp = detector->GetFreq(l); // freq in Hz

										cout << "Check 1" << endl;
										/*
										 // Get ant gain with 2-D interpolation (may have bug?) 
										 //
										 heff = GaintoHeight(detector->stations[i].strings[j].antennas[k].GetG(detector, freq_tmp*1.E-6, // to MHz
										 antenna_theta, antenna_phi), 
										 freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]) );
										 */
										if (settings1->ALL_ANT_V_ON == 0) {
											if (settings1->ANTENNA_MODE == 0) {
												heff =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		antenna_theta,
																		antenna_phi,
																		detector->stations[i].strings[j].antennas[k].type),
																freq_tmp,
																icemodel->GetN(
																		detector->stations[i].strings[j].antennas[k]));
											}
											if (settings1->ANTENNA_MODE == 1) {
												heff =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		antenna_theta,
																		antenna_phi,
																		detector->stations[i].strings[j].antennas[k].type,
																		k),
																freq_tmp,
																icemodel->GetN(
																		detector->stations[i].strings[j].antennas[k]));
											}
										} else if (settings1->ALL_ANT_V_ON
												== 1) {
											if (settings1->ANTENNA_MODE == 0) {
												heff =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		antenna_theta,
																		antenna_phi,
																		0),
																freq_tmp,
																icemodel->GetN(
																		detector->stations[i].strings[j].antennas[k]));
											}
											if (settings1->ANTENNA_MODE == 1) {
												heff =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		antenna_theta,
																		antenna_phi,
																		0, k),
																freq_tmp,
																icemodel->GetN(
																		detector->stations[i].strings[j].antennas[k]));
											}
										}

										//cout<<"n_medium : "<<icemodel->GetN(detector->stations[i].strings[j].antennas[k])<<endl;
										//cout<<"gain : "<<detector->stations[i].strings[j].antennas[k].GetG(detector, freq_tmp*1.E-6, antenna_theta, antenna_phi)<<endl;
										//cout<<"heff : "<<heff<<endl;
										stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(
												heff);

										// apply pol factor, heff
										if (event->IsCalpulser == 1) {
											//cout<<"set signal pol as Hpol for Calpulser1 evts"<<endl;
											Pol_vector = n_trg_slappy;
										} else if (event->IsCalpulser == 2) {
											//cout<<"set signal pol as Vpol for Calpulser2 evts"<<endl;
											Pol_vector = n_trg_pokey;
										} else if (event->IsCalpulser == 3) {
											//cout<<"set signal pol as Hpol for Calpulser2 evts"<<endl;
											Pol_vector = n_trg_slappy;
										} else if (event->IsCalpulser == 4) {
											//cout<<"set signal pol as Vpol + Hpol for Calpulser2 evts"<<endl;
											Pol_vector = n_trg_slappy
													+ n_trg_pokey;
										}

										ApplyAntFactors(heff, n_trg_pokey,
												n_trg_slappy, Pol_vector,
												detector->stations[i].strings[j].antennas[k].type,
												Pol_factor, vmmhz1m_tmp);

										cout << "Check 2" << endl;
										//stations[i].strings[j].antennas[k].VHz_antfactor[ray_sol_cnt].push_back( vmmhz1m_tmp );

										// apply filter
										ApplyFilter(l, detector, vmmhz1m_tmp);

										// apply Preamp gain
										ApplyPreamp(l, detector, vmmhz1m_tmp);

										// apply FOAM gain
										ApplyFOAM(l, detector, vmmhz1m_tmp);

										//stations[i].strings[j].antennas[k].VHz_filter[ray_sol_cnt].push_back( vmmhz1m_tmp );
										vmmhz_filter[l] = vmmhz1m_tmp;

									}                        // end for freq bin

									stations[i].strings[j].antennas[k].Pol_factor.push_back(
											Pol_factor);

									//cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"].vmmhz1m["<<ray_sol_cnt<<"][0] : "<<stations[i].strings[j].antennas[k].vmmhz[ray_sol_cnt][0]<<endl;

									//MakeArraysforFFT(settings1, detector, i, stations[i].strings[j].antennas[k].VHz_filter[ray_sol_cnt], volts_forfft);
									MakeArraysforFFT(settings1, detector, i,
											vmmhz_filter, volts_forfft);

									// save freq domain array which is prepaired for realft
									for (int n = 0; n < settings1->NFOUR / 2;
											n++) {
										stations[i].strings[j].antennas[k].Vfft[ray_sol_cnt].push_back(
												volts_forfft[n]);
									}

									// now, after realft, volts_forfft is time domain signal at backend of antenna
									Tools::realft(volts_forfft, -1,
											settings1->NFOUR / 2);
									//Tools::realft(volts_forfft,1,settings1->NFOUR/2);

									//cout<<"Finished getting V signal part!!"<<endl;

									stations[i].strings[j].antennas[k].PeakV.push_back(
											FindPeak(volts_forfft,
													settings1->NFOUR / 2));

									// Vfft_noise_org is in fft freq bin!!
									// same unit with Vfft [V] but filter not applied

									Tools::NormalTimeOrdering(
											settings1->NFOUR / 2, volts_forfft);
									//cout<<"finished NormalTimeOrdering!!"<<endl;

									for (int n = 0; n < settings1->NFOUR / 2;
											n++) {

										if (settings1->TRIG_ANALYSIS_MODE
												!= 2) { // not pure noise mode (we need signal)
											stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(
													volts_forfft[n]);
										} else if (settings1->TRIG_ANALYSIS_MODE
												== 2) { // pure noise mode (set signal to 0)
											stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(
													0.);
										}

										//stations[i].strings[j].antennas[k].time[ray_sol_cnt].push_back( stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt] + (double)(n - settings1->NFOUR/4)* settings1->TIMESTEP );   // time at 0 s is when ray started at the posnu
									}
									//cout<<"finished push_back V, V_noise V_total, and time!!"<<endl;

								} // if SIMULATION_MODE = 0

								else if (settings1->SIMULATION_MODE == 1) {

									// if event is not calpulser
									if (event->IsCalpulser == 0) {

										if (settings1->EVENT_TYPE == 0) {
											// see if integrated shower profile LQ is non-zero
											//
											// and near the cone viewangle
											//
											if (event->Nu_Interaction[0].LQ > 0
													&& (fabs(
															viewangle
																	- signal->CHANGLE_ICE)
															<= settings1->OFFCONE_LIMIT
																	* RADDEG)) {

												// initially give raysol has actual signal
												stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] =
														1;

												// let's make NFOUR/2 bin of time domain pure signal part for now
												// later once we understand how to apply antenna phase, total electronics with phase, apply those

												//double atten_factor = 1. / ray_output[0][ray_sol_cnt] * exp(-ray_output[0][ray_sol_cnt]/icemodel->EffectiveAttenuationLength(settings1, event->Nu_Interaction[0].posnu, 0)) * mag * fresnel;  // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)
												double atten_factor =
														1.
																/ ray_output[0][ray_sol_cnt]
																* IceAttenFactor
																* mag * fresnel; // assume whichray = 0, now vmmhz1m_tmp has all factors except for the detector properties (antenna gain, etc)

												// signal before the antenna (get signal at 1m and apply atten factor)
												signal->GetVm_FarField_Tarray(
														event, settings1,
														viewangle, atten_factor,
														outbin, Tarray, Earray,
														stations[i].strings[j].antennas[k].skip_bins[ray_sol_cnt]);

												dT_forfft = Tarray[1]
														- Tarray[0]; // step in ns

												int Ntmp = settings1->TIMESTEP
														* 1.e9 / dT_forfft;
												stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] = 1;
												while (Ntmp > 1) {
													Ntmp = Ntmp / 2;
													stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] =
															stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																	* 2;
												}
												stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] =
														stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																* settings1->NFOUR
																/ 2;
												// now new NFOUR for zero padding

												// now we have to make NFOUR/2 number of bins with random init time
												// as a test, make first as it is and zero pad

												double V_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];
												double T_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];

												for (int n = 0; n < stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt];n++) {

													if (n < outbin) {
														stations[i].strings[j].antennas[k].Vm_zoom[ray_sol_cnt].push_back(
																Earray[n]);
														stations[i].strings[j].antennas[k].Vm_zoom_T[ray_sol_cnt].push_back(
																Tarray[n]);
													}

													// make Tarray, Earray located at the center of Nnew array
													//T_forfft[n] = Tarray[outbin/2] - (dT_forfft*(double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]/2)) + (double)n*dT_forfft;

													T_forfft[n] =
															Tarray[outbin / 2]
																	- (dT_forfft
																			* (double) (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																					/ 2
																					- n));

													if ((n
															>= stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																	/ 2
																	- outbin / 2)
															&& (n
																	< stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																			/ 2
																			+ outbin
																					/ 2))
													{
														V_forfft[n] =
																Earray[n
																		- (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																				/ 2
																				- outbin
																						/ 2)];
													} else{
														V_forfft[n] = 0.;
													}
													//stations[i].strings[j].antennas[k].Vm_wo_antfactor[ray_sol_cnt].push_back( V_forfft[n] );
												}

												// just get peak from the array
												stations[i].strings[j].antennas[k].PeakV.push_back(
														FindPeak(Earray,
																outbin));

												// this forward fft volts_forfft is now in unit of V at each freq we can just apply each bin's gain factor to each freq bins
												//
												// without any phase consideration,
												// apply same factor to both real, img parts
												//
												//
												//
												//
												//
												//
												//
												//
												// get spectrum with zero padded WF
												//Tools::realft(volts_forfft,1,settings1->NFOUR/2);
												Tools::realft(V_forfft, 1,
														stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

												//dF_outbin = 1./( (double)(outbin) * (Tarray[1]-Tarray[0])*1.e-9 ); // in Hz
												dF_Nnew =
														1.
																/ ((double) (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt])
																		* (dT_forfft)
																		* 1.e-9); // in Hz

												//freq_tmp = dF_Nnew*((double)stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]/2.+1.+0.5);// in Hz 0.5 to place the middle of the bin and avoid zero freq
												freq_tmp =
														dF_Nnew
																* ((double) stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																		/ 2.
																		+ 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

												freq_lastbin = freq_tmp;

												/*
												 // Get ant gain with 2-D interpolation (may have bug?) 
												 //
												 heff_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
												 antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
												 freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]) );
												 */
												if (settings1->ALL_ANT_V_ON
														== 0) {
													if (settings1->ANTENNA_MODE
															== 0) {
														heff_lastbin =
																GaintoHeight(
																		detector->GetGain_1D_OutZero(
																				freq_tmp
																						* 1.E-6, // to MHz
																				antenna_theta,
																				antenna_phi,
																				detector->stations[i].strings[j].antennas[k].type),
																		freq_tmp,
																		icemodel->GetN(
																				detector->stations[i].strings[j].antennas[k]));
													}
													if (settings1->ANTENNA_MODE
															== 1) {
														heff_lastbin =
																GaintoHeight(
																		detector->GetGain_1D_OutZero(
																				freq_tmp
																						* 1.E-6, // to MHz
																				antenna_theta,
																				antenna_phi,
																				detector->stations[i].strings[j].antennas[k].type,
																				k),
																		freq_tmp,
																		icemodel->GetN(
																				detector->stations[i].strings[j].antennas[k]));
													}
												} else if (settings1->ALL_ANT_V_ON
														== 1) {

													if (settings1->ANTENNA_MODE
															== 0) {
														heff_lastbin =
																GaintoHeight(
																		detector->GetGain_1D_OutZero(
																				freq_tmp
																						* 1.E-6, // to MHz
																				antenna_theta,
																				antenna_phi,
																				0),
																		freq_tmp,
																		icemodel->GetN(
																				detector->stations[i].strings[j].antennas[k]));
													}
													if (settings1->ANTENNA_MODE
															== 1) {
														heff_lastbin =
																GaintoHeight(
																		detector->GetGain_1D_OutZero(
																				freq_tmp
																						* 1.E-6, // to MHz
																				antenna_theta,
																				antenna_phi,
																				0,
																				k),
																		freq_tmp,
																		icemodel->GetN(
																				detector->stations[i].strings[j].antennas[k]));
													}

												}

												//
												//
												for (int n = 0;
														n
																< stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																		/ 2;
														n++) {

													freq_tmp =
															dF_Nnew
																	* ((double) n
																			+ 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

													/*
													 heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
													 antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
													 freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]) );
													 */
													if (settings1->ALL_ANT_V_ON
															== 0) {
														if (settings1->ANTENNA_MODE
																== 0) {
															heff =
																	GaintoHeight(
																			detector->GetGain_1D_OutZero(
																					freq_tmp
																							* 1.E-6, // to MHz
																					antenna_theta,
																					antenna_phi,
																					detector->stations[i].strings[j].antennas[k].type),
																			freq_tmp,
																			icemodel->GetN(
																					detector->stations[i].strings[j].antennas[k]));
														}

														if (settings1->ANTENNA_MODE
																== 1) {
															heff =
																	GaintoHeight(
																			detector->GetGain_1D_OutZero(
																					freq_tmp
																							* 1.E-6, // to MHz
																					antenna_theta,
																					antenna_phi,
																					detector->stations[i].strings[j].antennas[k].type,
																					k),
																			freq_tmp,
																			icemodel->GetN(
																					detector->stations[i].strings[j].antennas[k]));
														}
													} else if (settings1->ALL_ANT_V_ON
															== 1) {
														if (settings1->ANTENNA_MODE
																== 0) {
															heff =
																	GaintoHeight(
																			detector->GetGain_1D_OutZero(
																					freq_tmp
																							* 1.E-6, // to MHz
																					antenna_theta,
																					antenna_phi,
																					0),
																			freq_tmp,
																			icemodel->GetN(
																					detector->stations[i].strings[j].antennas[k]));
														}
														if (settings1->ANTENNA_MODE
																== 1) {
															heff =
																	GaintoHeight(
																			detector->GetGain_1D_OutZero(
																					freq_tmp
																							* 1.E-6, // to MHz
																					antenna_theta,
																					antenna_phi,
																					0,
																					k),
																			freq_tmp,
																			icemodel->GetN(
																					detector->stations[i].strings[j].antennas[k]));
														}
													}

													stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(
															heff);

													//
													// apply ant factors
													//
													if (n > 0) {

														if (settings1->ALL_ANT_V_ON
																== 0) {

															ApplyAntFactors_Tdomain(
																	detector->GetAntPhase_1D(
																			freq_tmp
																					* 1.e-6,
																			antenna_theta,
																			antenna_phi,
																			detector->stations[i].strings[j].antennas[k].type),
																	heff,
																	n_trg_pokey,
																	n_trg_slappy,
																	Pol_vector,
																	detector->stations[i].strings[j].antennas[k].type,
																	Pol_factor,
																	V_forfft[2
																			* n],
																	V_forfft[2
																			* n
																			+ 1],
																	settings1);
														} else if (settings1->ALL_ANT_V_ON
																== 1) {
															ApplyAntFactors_Tdomain(
																	detector->GetAntPhase_1D(
																			freq_tmp
																					* 1.e-6,
																			antenna_theta,
																			antenna_phi,
																			0),
																	heff,
																	n_trg_pokey,
																	n_trg_slappy,
																	Pol_vector,
																	detector->stations[i].strings[j].antennas[k].type,
																	Pol_factor,
																	V_forfft[2
																			* n],
																	V_forfft[2
																			* n
																			+ 1],
																	settings1);
														}

													} else {
														ApplyAntFactors_Tdomain_FirstTwo(
																heff,
																heff_lastbin,
																n_trg_pokey,
																n_trg_slappy,
																Pol_vector,
																detector->stations[i].strings[j].antennas[k].type,
																Pol_factor,
																V_forfft[2 * n],
																V_forfft[2 * n
																		+ 1]);
													}

													//
													// apply entire elect chain gain, phase
													//
													if (n > 0) {
														ApplyElect_Tdomain(
																freq_tmp
																		* 1.e-6,
																detector,
																V_forfft[2 * n],
																V_forfft[2 * n
																		+ 1],
																settings1);
													} else {
														ApplyElect_Tdomain_FirstTwo(
																freq_tmp
																		* 1.e-6,
																freq_lastbin
																		* 1.e-6,
																detector,
																V_forfft[2 * n],
																V_forfft[2 * n
																		+ 1]);
													}

												}            // end for freq bin

															 // now get time domain waveform back by inv fft
												Tools::realft(V_forfft, -1,
														stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

												// we need to do normal time ordering as we did zero padding(?)
												// If signal is located at the center, we don't need to do NormalTimeOrdering???
												//Tools::NormalTimeOrdering(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], V_forfft);

												// skip linear interpolation for now
												Tools::SimpleLinearInterpolation_OutZero(
														stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt],
														T_forfft, V_forfft,
														settings1->NFOUR / 2,
														T_forint, volts_forint);

												// check what we save as V[], volts_forint? or volts_forfft

												for (int n = 0;
														n < settings1->NFOUR / 2;
														n++) {
													if (settings1->TRIG_ANALYSIS_MODE
															!= 2) { // not pure noise mode (we need signal)
														//stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back( volts_forfft[n] );
														//stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back( V_forfft[n] * 2./(double)(settings1->NFOUR/2) ); // 2/N for inverse FFT normalization factor
														//stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back( volts_forint[n] * 2./(double)(settings1->NFOUR/2) ); // 2/N for inverse FFT normalization factor
														stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(
																volts_forint[n]
																		* 2.
																		/ (double) (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt])); // 2/N for inverse FFT normalization factor
														// cout<<"Push back this V value "<<stations[i].strings[j].antennas[k].V[ray_sol_cnt][n]<<endl;
													} else if (settings1->TRIG_ANALYSIS_MODE
															== 2) { // pure noise mode (set signal to 0)
														stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(
																0.);
													}

												}

											}
											// see if integrated shower profile LQ is non-zero
											//
											// and near the cone viewangle

											else { // no signal generating

												// initially give raysol has actual signal
												stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] =
														0;

												// if no signal, push_back 0 values (otherwise the value inside will remain as old value)
												for (int n = 0;
														n < settings1->NFOUR / 2;
														n++) {

													if (n < outbin) {
														stations[i].strings[j].antennas[k].Vm_zoom[ray_sol_cnt].push_back(
																0.);
														stations[i].strings[j].antennas[k].Vm_zoom_T[ray_sol_cnt].push_back(
																n);
													}
													//stations[i].strings[j].antennas[k].Vm_wo_antfactor[ray_sol_cnt].push_back( 0. );
													//stations[i].strings[j].antennas[k].VHz_antfactor[ray_sol_cnt].push_back( 0. );
													//stations[i].strings[j].antennas[k].VHz_filter[ray_sol_cnt].push_back( 0. );
													stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(
															0.);
												}

												stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] =
														settings1->NFOUR / 2;
												// just get peak from the array
												stations[i].strings[j].antennas[k].PeakV.push_back(
														0.);
											}

										} else if (settings1->EVENT_TYPE== 10) {

											// initially give raysol has actual signal
											stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] =
													1;

											int waveform_bin =
													(int) signal->ArbitraryWaveform_V.size();
											//          cout << waveform_bin << endl;

											//dT_forfft = Tarray[1] - Tarray[0]; // step in ns
											//dT_forfft = detector->CalPulserWF_ns[1] - detector->CalPulserWF_ns[0]; // step in ns
											dT_forfft =
													signal->ArbitraryWaveform_T[1]
															- signal->ArbitraryWaveform_T[0]; // step in ns

											//          cout << "dT_forfft: " << dT_forfft << endl;
											int Ntmp = settings1->TIMESTEP
													* 1.e9 / dT_forfft;
											//          cout << Ntmp << endl;

											stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] =
													1;
											while (Ntmp > 1) {
												Ntmp = Ntmp / 2;
												stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] =
														stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																* 2;
											}
											stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] =
													stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
															* settings1->NFOUR
															/ 2;
											// now new NFOUR for zero padding

											// now we have to make NFOUR/2 number of bins with random init time
											//
											// as a test, make first as it is and zero pad

											double V_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];
											double T_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];

											for (int n = 0;
													n
															< stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt];
													n++) {
												//cout << n << endl;
												if (n < waveform_bin) {
													stations[i].strings[j].antennas[k].Vm_zoom[ray_sol_cnt].push_back(
															signal->ArbitraryWaveform_V[n]);
													stations[i].strings[j].antennas[k].Vm_zoom_T[ray_sol_cnt].push_back(
															signal->ArbitraryWaveform_T[n]);
													//cout << signal->ArbitraryWaveform_T[n] << " : " << signal->ArbitraryWaveform_V[n] << endl;
												}

												// make Tarray, Earray located at the center of Nnew array
												//T_forfft[n] = Tarray[outbin/2] - (dT_forfft*(double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]/2)) + (double)n*dT_forfft;

												T_forfft[n] =
														signal->ArbitraryWaveform_T[waveform_bin
																/ 2]
																- (dT_forfft
																		* (double) (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																				/ 2
																				- n));

												if ((n
														>= stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																/ 2
																- waveform_bin
																		/ 2)
														&& (n
																< stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																		/ 2
																		+ waveform_bin
																				/ 2)) {
													V_forfft[n] =
															signal->ArbitraryWaveform_V[n
																	- (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																			/ 2
																			- waveform_bin
																					/ 2)];

												} else
													V_forfft[n] = 0.;

												//stations[i].strings[j].antennas[k].Vm_wo_antfactor[ray_sol_cnt].push_back( V_forfft[n] );

											}

											// just get peak from the array
											//stations[i].strings[j].antennas[k].PeakV.push_back( FindPeak(detector->CalPulserWF_V, CP_bin) );
											stations[i].strings[j].antennas[k].PeakV.push_back(
													-1.);        // just let -1.

											// get spectrum with zero padded WF
											//Tools::realft(volts_forfft,1,settings1->NFOUR/2);
											Tools::realft(V_forfft, 1,
													stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

											//dF_outbin = 1./( (double)(outbin) * (Tarray[1]-Tarray[0])*1.e-9 ); // in Hz
											dF_Nnew =
													1.
															/ ((double) (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt])
																	* (dT_forfft)
																	* 1.e-9); // in Hz

											//freq_tmp = dF_Nnew*((double)stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]/2.+1.+0.5);// in Hz 0.5 to place the middle of the bin and avoid zero freq
											freq_tmp =
													dF_Nnew
															* ((double) stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																	/ 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

											freq_lastbin = freq_tmp;

											/*
											 // heff last bin for transmitter ant
											 double heff_lastbin_trans;
											 double ant_theta_trans = ray_output[1][ray_sol_cnt] * DEGRAD; // from 0 to 180
											 //cout<<"ant theta trans : "<<ant_theta_trans<<"deg"<<endl;
											 if ( settings1->ALL_ANT_V_ON==0 ) {
											 heff_lastbin_trans = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
											 ant_theta_trans, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
											 freq_tmp, icemodel->GetN(event->Nu_Interaction[0].posnu) );
											 }
											 else if ( settings1->ALL_ANT_V_ON==1 ) {
											 heff_lastbin_trans = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
											 ant_theta_trans, antenna_phi, 0), 
											 freq_tmp, icemodel->GetN(event->Nu_Interaction[0].posnu) );
											 }
											 
											 
											 
											 
											 // heff last bin for receiver ant
											 /*
											 // Get ant gain with 2-D interpolation (may have bug?) 
											 //
											 heff_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
											 antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
											 freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]) );
											 */
											/*          
											 if ( settings1->ALL_ANT_V_ON==0 ) {
											 heff_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
											 antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
											 freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]) );
											 }
											 else if ( settings1->ALL_ANT_V_ON==1 ) {
											 heff_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
											 antenna_theta, antenna_phi, 0), 
											 freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]) );
											 }
											 */
											// Pol_vector = n_trg_slappy;
											Pol_vector = n_trg_pokey;

											//
											//
											for (int n = 0;
													n
															< stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																	/ 2; n++) {

												freq_tmp = dF_Nnew
														* ((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq
												/*
												 //
												 // apply ant factors (transmitter ant)
												 //
												 if ( settings1->ALL_ANT_V_ON==0 ) {
												 heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
												 ant_theta_trans, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
												 freq_tmp, icemodel->GetN(event->Nu_Interaction[0].posnu) );
												 }
												 else if ( settings1->ALL_ANT_V_ON==1 ) {
												 heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
												 ant_theta_trans, antenna_phi, 0), 
												 freq_tmp, icemodel->GetN(event->Nu_Interaction[0].posnu) );
												 }
												 //
												 if ( n > 0 ) {
												 
												 if ( settings1->ALL_ANT_V_ON==0 ) {
												 
												 ApplyAntFactors_Tdomain_Transmitter( detector->GetAntPhase_1D( freq_tmp*1.e-6, ant_theta_trans, antenna_phi, detector->stations[i].strings[j].antennas[k].type ),
												 heff, n_trg_pokey, n_trg_slappy, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2*n], V_forfft[2*n+1], settings1 );
												 }
												 else if ( settings1->ALL_ANT_V_ON==1 ) {
												 ApplyAntFactors_Tdomain_Transmitter( detector->GetAntPhase_1D( freq_tmp*1.e-6, ant_theta_trans, antenna_phi, 0 ),
												 heff, n_trg_pokey, n_trg_slappy, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2*n], V_forfft[2*n+1], settings1 );
												 }
												 
												 
												 
												 }
												 else {
												 ApplyAntFactors_Tdomain_FirstTwo( heff, heff_lastbin_trans, n_trg_pokey, n_trg_slappy, Pol_vector, detector->stations[i].strings[j].antennas[k].type, Pol_factor, V_forfft[2*n], V_forfft[2*n+1] );
												 }
												 
												 //
												 // apply ant factors (receiver ant)
												 //
												 /*
												 heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
												 antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
												 freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]) );
												 */

												if (settings1->ALL_ANT_V_ON
														== 0) {
													if (settings1->ANTENNA_MODE
															== 0) {
														heff =
																GaintoHeight(
																		detector->GetGain_1D_OutZero(
																				freq_tmp
																						* 1.E-6, // to MHz
																				antenna_theta,
																				antenna_phi,
																				detector->stations[i].strings[j].antennas[k].type),
																		freq_tmp,
																		icemodel->GetN(
																				detector->stations[i].strings[j].antennas[k]));
													}
													if (settings1->ANTENNA_MODE
															== 1) {
														heff =
																GaintoHeight(
																		detector->GetGain_1D_OutZero(
																				freq_tmp
																						* 1.E-6, // to MHz
																				antenna_theta,
																				antenna_phi,
																				detector->stations[i].strings[j].antennas[k].type,
																				k),
																		freq_tmp,
																		icemodel->GetN(
																				detector->stations[i].strings[j].antennas[k]));
													}
												}

												else if (settings1->ALL_ANT_V_ON
														== 1) {
													if (settings1->ANTENNA_MODE
															== 0) {
														heff =
																GaintoHeight(
																		detector->GetGain_1D_OutZero(
																				freq_tmp
																						* 1.E-6, // to MHz
																				antenna_theta,
																				antenna_phi,
																				0),
																		freq_tmp,
																		icemodel->GetN(
																				detector->stations[i].strings[j].antennas[k]));

													}
													if (settings1->ANTENNA_MODE
															== 1) {
														heff =
																GaintoHeight(
																		detector->GetGain_1D_OutZero(
																				freq_tmp
																						* 1.E-6, // to MHz
																				antenna_theta,
																				antenna_phi,
																				0,
																				k),
																		freq_tmp,
																		icemodel->GetN(
																				detector->stations[i].strings[j].antennas[k]));

													}
												}

												stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(
														heff);

												if (n > 0) {

													if (settings1->ALL_ANT_V_ON
															== 0) {

														ApplyAntFactors_Tdomain(
																detector->GetAntPhase_1D(
																		freq_tmp
																				* 1.e-6,
																		antenna_theta,
																		antenna_phi,
																		detector->stations[i].strings[j].antennas[k].type),
																heff,
																n_trg_pokey,
																n_trg_slappy,
																Pol_vector,
																detector->stations[i].strings[j].antennas[k].type,
																Pol_factor,
																V_forfft[2 * n],
																V_forfft[2 * n
																		+ 1],
																settings1);
													} else if (settings1->ALL_ANT_V_ON
															== 1) {
														ApplyAntFactors_Tdomain(
																detector->GetAntPhase_1D(
																		freq_tmp
																				* 1.e-6,
																		antenna_theta,
																		antenna_phi,
																		0),
																heff,
																n_trg_pokey,
																n_trg_slappy,
																Pol_vector,
																detector->stations[i].strings[j].antennas[k].type,
																Pol_factor,
																V_forfft[2 * n],
																V_forfft[2 * n
																		+ 1],
																settings1);
													}

												} else {
													ApplyAntFactors_Tdomain_FirstTwo(
															heff, heff_lastbin,
															n_trg_pokey,
															n_trg_slappy,
															Pol_vector,
															detector->stations[i].strings[j].antennas[k].type,
															Pol_factor,
															V_forfft[2 * n],
															V_forfft[2 * n + 1]);
												}

												//
												// apply entire elect chain gain, phase
												//
												if (n > 0) {
													ApplyElect_Tdomain(
															freq_tmp * 1.e-6,
															detector,
															V_forfft[2 * n],
															V_forfft[2 * n + 1],
															settings1);
												} else {
													ApplyElect_Tdomain_FirstTwo(
															freq_tmp * 1.e-6,
															freq_lastbin
																	* 1.e-6,
															detector,
															V_forfft[2 * n],
															V_forfft[2 * n + 1]);
												}

											}             // end for freq bin

														  // now get time domain waveform back by inv fft
											Tools::realft(V_forfft, -1,
													stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

											// we need to do normal time ordering as we did zero padding(?)
											// If signal is located at the center, we don't need to do NormalTimeOrdering???
											//Tools::NormalTimeOrdering(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], V_forfft);

											// skip linear interpolation for now
											Tools::SimpleLinearInterpolation_OutZero(
													stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt],
													T_forfft, V_forfft,
													settings1->NFOUR / 2,
													T_forint, volts_forint);

											// check what we save as V[], volts_forint? or volts_forfft

											for (int n = 0;
													n < settings1->NFOUR / 2;
													n++) {

												if (settings1->TRIG_ANALYSIS_MODE
														!= 2) { // not pure noise mode (we need signal)
													//stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back( volts_forfft[n] );
													//stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back( V_forfft[n] * 2./(double)(settings1->NFOUR/2) ); // 2/N for inverse FFT normalization factor
													//stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back( volts_forint[n] * 2./(double)(settings1->NFOUR/2) ); // 2/N for inverse FFT normalization factor
													stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(
															settings1->ARBITRARY_EVENT_ATTENUATION
																	* volts_forint[n]
																	* 2.
																	/ (double) (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt])); // 2/N for inverse FFT normalization factor
												} else if (settings1->TRIG_ANALYSIS_MODE
														== 2) { // pure noise mode (set signal to 0)
													stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(
															0.);
												}

											}

										} // Arbitrary Events

									} // if not calpulser event

									// if calpulser event
									else if (event->IsCalpulser > 0) {

										// initially give raysol has actual signal
										stations[i].strings[j].antennas[k].SignalExt[ray_sol_cnt] =
												1;

										int CP_bin =
												(int) detector->CalPulserWF_ns.size();

										//dT_forfft = Tarray[1] - Tarray[0]; // step in ns
										dT_forfft = detector->CalPulserWF_ns[1]
												- detector->CalPulserWF_ns[0]; // step in ns

										int Ntmp = settings1->TIMESTEP * 1.e9
												/ dT_forfft;
										stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] =
												1;
										while (Ntmp > 1) {
											Ntmp = Ntmp / 2;
											stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] =
													stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
															* 2;
										}
										stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt] =
												stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
														* settings1->NFOUR / 2;
										// now new NFOUR for zero padding

										// now we have to make NFOUR/2 number of bins with random init time
										//
										// as a test, make first as it is and zero pad

										double V_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];
										double T_forfft[stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]];

										for (int n = 0;
												n
														< stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt];
												n++) {

											if (n < CP_bin) {
												stations[i].strings[j].antennas[k].Vm_zoom[ray_sol_cnt].push_back(
														detector->CalPulserWF_V[n]);
												stations[i].strings[j].antennas[k].Vm_zoom_T[ray_sol_cnt].push_back(
														detector->CalPulserWF_ns[n]);

											}

											// make Tarray, Earray located at the center of Nnew array
											//T_forfft[n] = Tarray[outbin/2] - (dT_forfft*(double)(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]/2)) + (double)n*dT_forfft;

											T_forfft[n] =
													detector->CalPulserWF_ns[CP_bin
															/ 2]
															- (dT_forfft
																	* (double) (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																			/ 2
																			- n));

											if ((n
													>= stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
															/ 2 - CP_bin / 2)
													&& (n
															< stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																	/ 2
																	+ CP_bin / 2)) {
												V_forfft[n] =
														detector->CalPulserWF_V[n
																- (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																		/ 2
																		- CP_bin
																				/ 2)];
											} else
												V_forfft[n] = 0.;

											//stations[i].strings[j].antennas[k].Vm_wo_antfactor[ray_sol_cnt].push_back( V_forfft[n] );

										}

										// just get peak from the array
										//stations[i].strings[j].antennas[k].PeakV.push_back( FindPeak(detector->CalPulserWF_V, CP_bin) );
										stations[i].strings[j].antennas[k].PeakV.push_back(
												-1.);           // just let -1.

										// get spectrum with zero padded WF
										//Tools::realft(volts_forfft,1,settings1->NFOUR/2);
										Tools::realft(V_forfft, 1,
												stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

										//dF_outbin = 1./( (double)(outbin) * (Tarray[1]-Tarray[0])*1.e-9 ); // in Hz
										dF_Nnew =
												1.
														/ ((double) (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt])
																* (dT_forfft)
																* 1.e-9); // in Hz

										//freq_tmp = dF_Nnew*((double)stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]/2.+1.+0.5);// in Hz 0.5 to place the middle of the bin and avoid zero freq
										freq_tmp =
												dF_Nnew
														* ((double) stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																/ 2. + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

										freq_lastbin = freq_tmp;

										// heff last bin for transmitter ant
										double heff_lastbin_trans;
										double ant_theta_trans =
												ray_output[1][ray_sol_cnt]
														* DEGRAD; // from 0 to 180
										//cout<<"ant theta trans : "<<ant_theta_trans<<"deg"<<endl;
										if (settings1->ALL_ANT_V_ON == 0) {
											if (settings1->ANTENNA_MODE == 0) {
												heff_lastbin_trans =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		ant_theta_trans,
																		antenna_phi,
																		detector->stations[i].strings[j].antennas[k].type),
																freq_tmp,
																icemodel->GetN(
																		event->Nu_Interaction[0].posnu));
											}
											if (settings1->ANTENNA_MODE == 1) {
												heff_lastbin_trans =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		ant_theta_trans,
																		antenna_phi,
																		detector->stations[i].strings[j].antennas[k].type,
																		k),
																freq_tmp,
																icemodel->GetN(
																		event->Nu_Interaction[0].posnu));
											}
										} else if (settings1->ALL_ANT_V_ON
												== 1) {
											if (settings1->ANTENNA_MODE == 0) {
												heff_lastbin_trans =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		ant_theta_trans,
																		antenna_phi,
																		0),
																freq_tmp,
																icemodel->GetN(
																		event->Nu_Interaction[0].posnu));
											}
											if (settings1->ANTENNA_MODE == 1) {
												heff_lastbin_trans =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		ant_theta_trans,
																		antenna_phi,
																		0, k),
																freq_tmp,
																icemodel->GetN(
																		event->Nu_Interaction[0].posnu));
											}
										}

										// heff last bin for receiver ant
										/*
										 // Get ant gain with 2-D interpolation (may have bug?) 
										 //
										 heff_lastbin = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
										 antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
										 freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]) );
										 */
										if (settings1->ALL_ANT_V_ON == 0) {
											if (settings1->ANTENNA_MODE == 0) {
												heff_lastbin =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		antenna_theta,
																		antenna_phi,
																		detector->stations[i].strings[j].antennas[k].type),
																freq_tmp,
																icemodel->GetN(
																		detector->stations[i].strings[j].antennas[k]));
											}
											if (settings1->ANTENNA_MODE == 1) {
												heff_lastbin =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		antenna_theta,
																		antenna_phi,
																		detector->stations[i].strings[j].antennas[k].type,
																		k),
																freq_tmp,
																icemodel->GetN(
																		detector->stations[i].strings[j].antennas[k]));
											}
										} else if (settings1->ALL_ANT_V_ON
												== 1) {
											if (settings1->ANTENNA_MODE == 0) {
												heff_lastbin =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		antenna_theta,
																		antenna_phi,
																		0),
																freq_tmp,
																icemodel->GetN(
																		detector->stations[i].strings[j].antennas[k]));
											}
											if (settings1->ANTENNA_MODE == 1) {
												heff_lastbin =
														GaintoHeight(
																detector->GetGain_1D_OutZero(
																		freq_tmp
																				* 1.E-6, // to MHz
																		antenna_theta,
																		antenna_phi,
																		0, k),
																freq_tmp,
																icemodel->GetN(
																		detector->stations[i].strings[j].antennas[k]));
											}
										}

										// apply calpulser waveform
										// apply pol factor, heff
										if (event->IsCalpulser == 1) {
											//cout<<"set signal pol as Hpol for Calpulser1 evts"<<endl;
											Pol_vector = n_trg_slappy;
										} else if (event->IsCalpulser == 2) {
											//cout<<"set signal pol as Vpol for Calpulser2 evts"<<endl;
											Pol_vector = n_trg_pokey;
										} else if (event->IsCalpulser == 3) {
											//cout<<"set signal pol as Hpol for Calpulser2 evts"<<endl;
											Pol_vector = n_trg_slappy;
										} else if (event->IsCalpulser == 4) {
											//cout<<"set signal pol as Vpol + Hpol for Calpulser2 evts"<<endl;
											Pol_vector = n_trg_slappy
													+ n_trg_pokey;
										}

										//
										//
										for (int n = 0;
												n
														< stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]
																/ 2; n++) {

											freq_tmp = dF_Nnew
													* ((double) n + 0.5); // in Hz 0.5 to place the middle of the bin and avoid zero freq

											//
											// apply ant factors (transmitter ant)
											//
											if (settings1->ALL_ANT_V_ON == 0) {
												if (settings1->ANTENNA_MODE
														== 0) {
													heff =
															GaintoHeight(
																	detector->GetGain_1D_OutZero(
																			freq_tmp
																					* 1.E-6, // to MHz
																			ant_theta_trans,
																			antenna_phi,
																			detector->stations[i].strings[j].antennas[k].type),
																	freq_tmp,
																	icemodel->GetN(
																			event->Nu_Interaction[0].posnu));
												}
												if (settings1->ANTENNA_MODE
														== 1) {
													heff =
															GaintoHeight(
																	detector->GetGain_1D_OutZero(
																			freq_tmp
																					* 1.E-6, // to MHz
																			ant_theta_trans,
																			antenna_phi,
																			detector->stations[i].strings[j].antennas[k].type,
																			k),
																	freq_tmp,
																	icemodel->GetN(
																			event->Nu_Interaction[0].posnu));
												}
											} else if (settings1->ALL_ANT_V_ON
													== 1) {
												if (settings1->ANTENNA_MODE
														== 0) {
													heff =
															GaintoHeight(
																	detector->GetGain_1D_OutZero(
																			freq_tmp
																					* 1.E-6, // to MHz
																			ant_theta_trans,
																			antenna_phi,
																			0),
																	freq_tmp,
																	icemodel->GetN(
																			event->Nu_Interaction[0].posnu));
												}
												if (settings1->ANTENNA_MODE
														== 1) {
													heff =
															GaintoHeight(
																	detector->GetGain_1D_OutZero(
																			freq_tmp
																					* 1.E-6, // to MHz
																			ant_theta_trans,
																			antenna_phi,
																			0,
																			k),
																	freq_tmp,
																	icemodel->GetN(
																			event->Nu_Interaction[0].posnu));
												}
											}
											//
											if (n > 0) {

												if (settings1->ALL_ANT_V_ON
														== 0) {

													ApplyAntFactors_Tdomain_Transmitter(
															detector->GetAntPhase_1D(
																	freq_tmp
																			* 1.e-6,
																	ant_theta_trans,
																	antenna_phi,
																	detector->stations[i].strings[j].antennas[k].type),
															heff, n_trg_pokey,
															n_trg_slappy,
															Pol_vector,
															detector->stations[i].strings[j].antennas[k].type,
															Pol_factor,
															V_forfft[2 * n],
															V_forfft[2 * n + 1],
															settings1);
												} else if (settings1->ALL_ANT_V_ON
														== 1) {
													ApplyAntFactors_Tdomain_Transmitter(
															detector->GetAntPhase_1D(
																	freq_tmp
																			* 1.e-6,
																	ant_theta_trans,
																	antenna_phi,
																	0), heff,
															n_trg_pokey,
															n_trg_slappy,
															Pol_vector,
															detector->stations[i].strings[j].antennas[k].type,
															Pol_factor,
															V_forfft[2 * n],
															V_forfft[2 * n + 1],
															settings1);
												}

											} else {
												ApplyAntFactors_Tdomain_FirstTwo(
														heff,
														heff_lastbin_trans,
														n_trg_pokey,
														n_trg_slappy,
														Pol_vector,
														detector->stations[i].strings[j].antennas[k].type,
														Pol_factor,
														V_forfft[2 * n],
														V_forfft[2 * n + 1]);
											}

											//
											// apply ant factors (receiver ant)
											//
											/*
											 heff = GaintoHeight(detector->GetGain_1D_OutZero(freq_tmp*1.E-6, // to MHz
											 antenna_theta, antenna_phi, detector->stations[i].strings[j].antennas[k].type), 
											 freq_tmp, icemodel->GetN(detector->stations[i].strings[j].antennas[k]) );
											 */
											if (settings1->ALL_ANT_V_ON == 0) {
												if (settings1->ANTENNA_MODE
														== 0) {
													heff =
															GaintoHeight(
																	detector->GetGain_1D_OutZero(
																			freq_tmp
																					* 1.E-6, // to MHz
																			antenna_theta,
																			antenna_phi,
																			detector->stations[i].strings[j].antennas[k].type),
																	freq_tmp,
																	icemodel->GetN(
																			detector->stations[i].strings[j].antennas[k]));
												}
												if (settings1->ANTENNA_MODE
														== 1) {
													heff =
															GaintoHeight(
																	detector->GetGain_1D_OutZero(
																			freq_tmp
																					* 1.E-6, // to MHz
																			antenna_theta,
																			antenna_phi,
																			detector->stations[i].strings[j].antennas[k].type,
																			k),
																	freq_tmp,
																	icemodel->GetN(
																			detector->stations[i].strings[j].antennas[k]));
												}
											} else if (settings1->ALL_ANT_V_ON
													== 1) {
												if (settings1->ANTENNA_MODE
														== 0) {
													heff =
															GaintoHeight(
																	detector->GetGain_1D_OutZero(
																			freq_tmp
																					* 1.E-6, // to MHz
																			antenna_theta,
																			antenna_phi,
																			0),
																	freq_tmp,
																	icemodel->GetN(
																			detector->stations[i].strings[j].antennas[k]));
												}
												if (settings1->ANTENNA_MODE
														== 1) {
													heff =
															GaintoHeight(
																	detector->GetGain_1D_OutZero(
																			freq_tmp
																					* 1.E-6, // to MHz
																			antenna_theta,
																			antenna_phi,
																			0,
																			k),
																	freq_tmp,
																	icemodel->GetN(
																			detector->stations[i].strings[j].antennas[k]));
												}
											}

											stations[i].strings[j].antennas[k].Heff[ray_sol_cnt].push_back(
													heff);

											if (n > 0) {

												if (settings1->ALL_ANT_V_ON
														== 0) {

													ApplyAntFactors_Tdomain(
															detector->GetAntPhase_1D(
																	freq_tmp
																			* 1.e-6,
																	antenna_theta,
																	antenna_phi,
																	detector->stations[i].strings[j].antennas[k].type),
															heff, n_trg_pokey,
															n_trg_slappy,
															Pol_vector,
															detector->stations[i].strings[j].antennas[k].type,
															Pol_factor,
															V_forfft[2 * n],
															V_forfft[2 * n + 1],
															settings1);
												} else if (settings1->ALL_ANT_V_ON
														== 1) {
													ApplyAntFactors_Tdomain(
															detector->GetAntPhase_1D(
																	freq_tmp
																			* 1.e-6,
																	antenna_theta,
																	antenna_phi,
																	0), heff,
															n_trg_pokey,
															n_trg_slappy,
															Pol_vector,
															detector->stations[i].strings[j].antennas[k].type,
															Pol_factor,
															V_forfft[2 * n],
															V_forfft[2 * n + 1],
															settings1);
												}

											} else {
												ApplyAntFactors_Tdomain_FirstTwo(
														heff, heff_lastbin,
														n_trg_pokey,
														n_trg_slappy,
														Pol_vector,
														detector->stations[i].strings[j].antennas[k].type,
														Pol_factor,
														V_forfft[2 * n],
														V_forfft[2 * n + 1]);
											}

											//
											// apply entire elect chain gain, phase
											//
											if (n > 0) {
												ApplyElect_Tdomain(
														freq_tmp * 1.e-6,
														detector,
														V_forfft[2 * n],
														V_forfft[2 * n + 1],
														settings1);
											} else {
												ApplyElect_Tdomain_FirstTwo(
														freq_tmp * 1.e-6,
														freq_lastbin * 1.e-6,
														detector,
														V_forfft[2 * n],
														V_forfft[2 * n + 1]);
											}

										}             // end for freq bin

													  // now get time domain waveform back by inv fft
										Tools::realft(V_forfft, -1,
												stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt]);

										// we need to do normal time ordering as we did zero padding(?)
										// If signal is located at the center, we don't need to do NormalTimeOrdering???
										//Tools::NormalTimeOrdering(stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt], V_forfft);

										// skip linear interpolation for now
										Tools::SimpleLinearInterpolation_OutZero(
												stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt],
												T_forfft, V_forfft,
												settings1->NFOUR / 2, T_forint,
												volts_forint);

										// check what we save as V[], volts_forint? or volts_forfft

										for (int n = 0;
												n < settings1->NFOUR / 2; n++) {

											if (settings1->TRIG_ANALYSIS_MODE
													!= 2) { // not pure noise mode (we need signal)
												//stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back( volts_forfft[n] );
												//stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back( V_forfft[n] * 2./(double)(settings1->NFOUR/2) ); // 2/N for inverse FFT normalization factor
												//stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back( volts_forint[n] * 2./(double)(settings1->NFOUR/2) ); // 2/N for inverse FFT normalization factor
												stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(
														volts_forint[n] * 2.
																/ (double) (stations[i].strings[j].antennas[k].Nnew[ray_sol_cnt])); // 2/N for inverse FFT normalization factor
											} else if (settings1->TRIG_ANALYSIS_MODE
													== 2) { // pure noise mode (set signal to 0)
												stations[i].strings[j].antennas[k].V[ray_sol_cnt].push_back(
														0.);
											}

										}

									} // Calpulser events

								} // if SIMULATION_MODE = 1

								// check max_PeakV
								if (max_PeakV_tmp
										< stations[i].strings[j].antennas[k].PeakV[ray_sol_cnt]) {
									max_PeakV_tmp =
											stations[i].strings[j].antennas[k].PeakV[ray_sol_cnt];
								}

							} // if not debug mode

							// check min_arrival_time
							if (min_arrival_time_tmp
									> stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt]) {
								min_arrival_time_tmp =
										stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt];
							}

							// check max_arrival_time
							if (max_arrival_time_tmp
									< stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt]) {
								max_arrival_time_tmp =
										stations[i].strings[j].antennas[k].arrival_time[ray_sol_cnt];
							}

							ray_sol_cnt++;

						} // end while number of solutions

					} // end if solution exist

					else {

						//cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"].trg = "<<stations[i].strings[j].antennas[k].trg[ray_sol_cnt]<<"  No vmmhz1m data!"<<endl;
					}

				} // end if posnu selected

				else {
					//cout<<" No posnu!!!!!! No signals calculated at all!!"<<endl;
					ray_sol_cnt = 0;
				}

				stations[i].strings[j].antennas[k].ray_sol_cnt = ray_sol_cnt; // save number of RaySolver solutions

				stations[i].Total_ray_sol += ray_sol_cnt; // add ray_sol_cnt to Total_ray_sol

			} // for number_of_antennas_string

		} // for number_of_strings_station

		// set each station's min/max arrival time
		stations[i].min_arrival_time = min_arrival_time_tmp;
		stations[i].max_arrival_time = max_arrival_time_tmp;

		// set each station's max PeakV
		stations[i].max_PeakV = max_PeakV_tmp;

	} // for number_of_stations

	// clear ray_output info
	ray_output.clear();

	// do only if it's not in debugmode
	if (debugmode == 0) {
		// after all values are stored in Report, set ranking of signal between antennas
		SetRank(detector);
	}

	// now loop over all antennas again to make DATA_BIN_SIZE array for signal + noise. (with time delay)
	// with that signal + noise array, we'll do convolution with diode response.
	// with the convolution result, we'll do trigger check
	for (int i = 0; i < detector->params.number_of_stations; i++) {

		N_pass = 0;
		N_pass_V = 0;
		N_pass_H = 0;

		stations[i].Global_Pass = 0;
		int check_passed_global_trigger = 0; // this switch determines if station globally triggers (in all TRIG_SCAN_MODEs)

		// int nfour_station = settings1->NFOUR;
		// double timestep_station = settings1->TIMESTEP;
		// double trig_window_station = settings1->TRIG_WINDOW;
		// int trig_window_bin = int(trig_window_station/timestep_station);
		// cout << nfour_station << " : " << timestep_station << " : " << trig_window_station << " : " << trig_window_bin << endl;

		// int nfour_station = detector->stations[i].NFOUR;
		// double timestep_station = detector->stations[i].TIMESTEP;
		// double trig_window_station = detector->stations[i].TRIG_WINDOW;
		// int trig_window_bin = int(trig_window_station/timestep_station);

		if (stations[i].Total_ray_sol) { // if there is any ray_sol (don't check trigger if there is no ray_sol at all)
			//if (stations[i].Total_ray_sol && settings1->TRIG_ANALYSIS_MODE != 2) { // if there is any ray_sol (don't check trigger if there is no ray_sol at all) and TRIG_ANALYSIS_MODE is 0 (signal + noise mode)

			//cout<<"3"<<endl;

			// calculate total number of bins we need to do trigger check
			//max_total_bin = (stations[i].max_arrival_time - stations[i].min_arrival_time)/settings1->TIMESTEP + settings1->NFOUR/2 + trigger->maxt_diode_bin;
			//max_total_bin = (stations[i].max_arrival_time - stations[i].min_arrival_time)/settings1->TIMESTEP + settings1->NFOUR*2 + trigger->maxt_diode_bin; // make more time
			max_total_bin = (stations[i].max_arrival_time
					- stations[i].min_arrival_time) / settings1->TIMESTEP
					+ settings1->NFOUR * 3 + trigger->maxt_diode_bin; // make more time

			//stations[i].max_total_bin = max_total_bin;

			// test generating new noise waveform for only stations if there's any ray trace solutions
			if (settings1->NOISE_WAVEFORM_GENERATE_MODE == 0) { // noise waveforms will be generated for each evts

				// redefine DATA_BIN_SIZE
				int DATA_BIN_SIZE_tmp;
				  // long DATA_BIN_SIZE_tmp;
				//for (int DBS=10; DBS<15; DBS++) {
				for (int DBS = 10; DBS < 16; DBS++) {
				  // for (int DBS=10; DBS<17; DBS++) {
					DATA_BIN_SIZE_tmp = (int) pow(2., (double) DBS);
					//if (DATA_BIN_SIZE_tmp > max_total_bin) DBS = 15; // come out
					if (DATA_BIN_SIZE_tmp > max_total_bin)
						DBS = 16; // come out
					  // if (DATA_BIN_SIZE_tmp > max_total_bin) DBS = 17; // come out
					  // if (DATA_BIN_SIZE_tmp > max_total_bin) DATA_BIN_SIZE_tmp = (long)pow(2., 17.0);
								// if (DATA_BIN_SIZE_tmp > max_total_bin+settings1->NFOUR/2) DBS = 15; // come out
				}
				settings1->DATA_BIN_SIZE = DATA_BIN_SIZE_tmp;
				  // cout<<"new DATA_BIN_SIZE : "<<DATA_BIN_SIZE_tmp<<endl;
				  // cout<<"max_total_bin : "<<max_total_bin<<endl;

				// reset all filter values in Detector class
				detector->get_NewDiodeModel(settings1);
				detector->ReadFilter_New(settings1);
				detector->ReadPreamp_New(settings1);
				detector->ReadFOAM_New(settings1);

				detector->ReadElectChain_New(settings1);

				if (settings1->USE_TESTBED_RFCM_ON == 1) {
					detector->ReadRFCM_New(settings1);
				}
				if (settings1->NOISE == 1 && settings1->DETECTOR == 3) {
					detector->ReadRayleigh_New(settings1);
				}
				  // if ( settings1->NOISE==1 && settings1->DETECTOR==4) {
					  // detector->ReadRayleigh_Station(settings1);
				  // }

				// reset Trigger class noise temp values
				trigger->Reset_V_noise_freqbin(settings1, detector);

				// now call new noise waveforms with new DATA_BIN_SIZE
				trigger->GetNewNoiseWaveforms(settings1, detector, this);
				//      cout << "New noise waveforms gotten" << endl;
			}

			// do only if it's not in debugmode
			if (debugmode == 1)
				N_noise = 1;
			// now, check if DATA_BIN_SIZE is enough for total time delay between antennas
			else {

				N_noise = (int) (max_total_bin / settings1->DATA_BIN_SIZE) + 1;
			}
			//cout<<"N_noise : "<<N_noise<<endl;

			if (N_noise > 1)
				cout << "N_noise : " << N_noise << " max_total_bin : "
						<< max_total_bin << " might cause error!!" << endl;
			// mostly N_noise should be "1"

			// now, check the number of bins we need for portion of noise waveforms
			remain_bin = max_total_bin % settings1->DATA_BIN_SIZE;

			ch_ID = 0;

			for (int j = 0; j < detector->stations[i].strings.size(); j++) {
				//cout << j << endl;
				for (int k = 0; k < detector->stations[i].strings[j].antennas.size(); k++) {

					// select noise waveform from trigger class
					for (int l = 0; l < N_noise; l++) {

						// if we are sharing same noise waveform for all chs, make sure diff chs use diff noise waveforms
						if (settings1->NOISE_CHANNEL_MODE == 0) {

							// get random noise_ID and check if there are same noise_ID in different ch.
							// if there's same noise_ID, get noise_ID again until no noise_ID are same between chs
							noise_pass_nogo = 1;
							while (noise_pass_nogo) {
								noise_ID[l] = (int) (settings1->NOISE_EVENTS
										* gRandom->Rndm());
								noise_pass_nogo = 0;
								for (int j_sub = 0; j_sub < j + 1; j_sub++) {
									for (int k_sub = 0;
											k_sub
													< detector->stations[i].strings[j].antennas.size();
											k_sub++) {
										if (j_sub == j) { // if we are checking current string
											if (k_sub < k) { // check antennas before current antenna
												if (noise_ID[l]
														== stations[i].strings[j_sub].antennas[k_sub].noise_ID[0]) { // check only first one for now;;;
													noise_pass_nogo = 1;
												}
											}
										} else { // if we are checking previous string, check upto entire antennas
											if (noise_ID[l]
													== stations[i].strings[j_sub].antennas[k_sub].noise_ID[0]) { // check only first one for now;;;
												noise_pass_nogo = 1;
											}
										}
									}
								}
							} // while noise_pass_nogo
						} // if NOISE_CHANNEL_MODE = 0

						// if we are using diff noise waveform for diff chs, just pick any noise waveform
						else if (settings1->NOISE_CHANNEL_MODE == 1) {
							noise_ID[l] = (int) (settings1->NOISE_EVENTS
									* gRandom->Rndm());
							//      cout << "picking noise waveform" << endl;
							//cout << noise_ID[l]<< endl;
						}             // if NOISE_CHANNEL_MODE = 1

						// if we are using diff noise waveform for diff chs, just pick any noise waveform
						else if (settings1->NOISE_CHANNEL_MODE == 2) {
							noise_ID[l] = (int) (settings1->NOISE_EVENTS
									* gRandom->Rndm());
						}             // if NOISE_CHANNEL_MODE = 1

						// save noise ID
						stations[i].strings[j].antennas[k].noise_ID.push_back(
								noise_ID[l]);

						//                           cout<<"noise_ID for "<<l<<"th noisewaveform is : "<<noise_ID[l]<<"  N_noise : "<<N_noise<<" ray_sol_cnt : "<<stations[i].strings[j].antennas[k].ray_sol_cnt<<endl;

						// do only if it's not in debugmode
						if (debugmode == 0) {

							//          cout << l << " : " << N_noise-1 << " : " << ch_ID << endl;
							// if we are sharing same noise waveform for all chs, make sure diff chs use diff noise waveforms
							if (settings1->NOISE_CHANNEL_MODE == 0) {
								if (l == N_noise - 1) { // when it's final noise waveform
								//for (int bin=0; bin<remain_bin; bin++) {
									for (int bin = 0;
											bin < settings1->DATA_BIN_SIZE;
											bin++) {   // test for full window
										trigger->Full_window[ch_ID][bin] =
												(trigger->v_noise_timedomain_diode[noise_ID[l]][bin]);
										trigger->Full_window_V[ch_ID][bin] =
												(trigger->v_noise_timedomain[noise_ID[l]][bin]);
									}
									//cout<<"last noise filled in Full_window!"<<endl;
								} else {   // when it's not final noise waveform
									//cout<<"full noise will fill in Full_window!"<<endl;
									for (int bin = 0;
											bin < settings1->DATA_BIN_SIZE;
											bin++) {
										trigger->Full_window[ch_ID][bin] =
												(trigger->v_noise_timedomain_diode[noise_ID[l]][bin]);
										trigger->Full_window_V[ch_ID][bin] =
												(trigger->v_noise_timedomain[noise_ID[l]][bin]);
									}
								}
							}

							// if we are sharing same noise waveform for all chs, make sure diff chs use diff noise waveforms
							else if (settings1->NOISE_CHANNEL_MODE == 1) {
								if (l == N_noise - 1) { // when it's final noise waveform
								//for (int bin=0; bin<remain_bin; bin++) {
									for (int bin = 0;
											bin < settings1->DATA_BIN_SIZE;
											bin++) {   // test for full window
										//          cout << GetChNumFromArbChID(detector,ch_ID,i,settings1)-1 << " : " << noise_ID[l] << " : " << ch_ID << endl;
										trigger->Full_window[ch_ID][bin] =
												(trigger->v_noise_timedomain_diode_ch[GetChNumFromArbChID(
														detector, ch_ID, i,
														settings1) - 1][noise_ID[l]][bin]);
										trigger->Full_window_V[ch_ID][bin] =
												(trigger->v_noise_timedomain_ch[GetChNumFromArbChID(
														detector, ch_ID, i,
														settings1) - 1][noise_ID[l]][bin]);
									}
									//              cout << "getting full window: " << ch_ID << endl;
									//cout<<"last noise filled in Full_window!"<<endl;
								} else {   // when it's not final noise waveform
									//cout<<"full noise will fill in Full_window!"<<endl;
									for (int bin = 0;
											bin < settings1->DATA_BIN_SIZE;
											bin++) {

										trigger->Full_window[ch_ID][bin] =
												(trigger->v_noise_timedomain_diode_ch[GetChNumFromArbChID(
														detector, ch_ID, i,
														settings1) - 1][noise_ID[l]][bin]);
										trigger->Full_window_V[ch_ID][bin] =
												(trigger->v_noise_timedomain_ch[GetChNumFromArbChID(
														detector, ch_ID, i,
														settings1) - 1][noise_ID[l]][bin]);
									}
								}
							}

							// if we are sharing same noise waveform for first 8 chs and share same noisewaveforms for others, make sure diff chs use diff noise waveforms
							else if (settings1->NOISE_CHANNEL_MODE == 2) {

								if ((GetChNumFromArbChID(detector, ch_ID, i,
										settings1) - 1) < 8) {

									if (l == N_noise - 1) { // when it's final noise waveform
									//for (int bin=0; bin<remain_bin; bin++) {
										for (int bin = 0;
												bin < settings1->DATA_BIN_SIZE;
												bin++) { // test for full window
											trigger->Full_window[ch_ID][bin] =
													(trigger->v_noise_timedomain_diode_ch[GetChNumFromArbChID(
															detector, ch_ID, i,
															settings1) - 1][noise_ID[l]][bin]);
											trigger->Full_window_V[ch_ID][bin] =
													(trigger->v_noise_timedomain_ch[GetChNumFromArbChID(
															detector, ch_ID, i,
															settings1) - 1][noise_ID[l]][bin]);
										}
										//cout<<"last noise filled in Full_window!"<<endl;
									} else { // when it's not final noise waveform
										//cout<<"full noise will fill in Full_window!"<<endl;
										for (int bin = 0;
												bin < settings1->DATA_BIN_SIZE;
												bin++) {
											trigger->Full_window[ch_ID][bin] =
													(trigger->v_noise_timedomain_diode_ch[GetChNumFromArbChID(
															detector, ch_ID, i,
															settings1) - 1][noise_ID[l]][bin]);
											trigger->Full_window_V[ch_ID][bin] =
													(trigger->v_noise_timedomain_ch[GetChNumFromArbChID(
															detector, ch_ID, i,
															settings1) - 1][noise_ID[l]][bin]);
										}
									}

								}

								else {

									if (l == N_noise - 1) { // when it's final noise waveform
									//for (int bin=0; bin<remain_bin; bin++) {
										for (int bin = 0;
												bin < settings1->DATA_BIN_SIZE;
												bin++) { // test for full window
											trigger->Full_window[ch_ID][bin] =
													(trigger->v_noise_timedomain_diode_ch[8][noise_ID[l]][bin]);
											trigger->Full_window_V[ch_ID][bin] =
													(trigger->v_noise_timedomain_ch[8][noise_ID[l]][bin]);
										}
										//cout<<"last noise filled in Full_window!"<<endl;
									} else { // when it's not final noise waveform
										//cout<<"full noise will fill in Full_window!"<<endl;
										for (int bin = 0;
												bin < settings1->DATA_BIN_SIZE;
												bin++) {
											trigger->Full_window[ch_ID][bin] =
													(trigger->v_noise_timedomain_diode_ch[8][noise_ID[l]][bin]);
											trigger->Full_window_V[ch_ID][bin] =
													(trigger->v_noise_timedomain_ch[8][noise_ID[l]][bin]);
										}
									}

								}

							}

						} // if we are not in debug mode

					}

					// currently there is a initial spoiled bins (maxt_diode_bin) at the initial Full_window "AND" at the initial of connected noisewaveform (we can fix this by adding values but not accomplished yet)

					//cout<<"done filling noise diode arrays for trigger!"<<endl;

					// now we need to calculate bin values for the signal
					// with the bin values, grap noise voltage waveform (NFOUR/2) and do myconvlv
					// and replace by init : bin + maxt_diode_bin, fin : bin + NFOUR/2
					//
					// after we do this for all channels, do trigger check

					// calculate the bin values for the signal
					signal_bin.clear();
					signal_dbin.clear();
					connect_signals.clear();

					// do only if it's not in debugmode
					if (debugmode == 0) {

						for (int m = 0; m < stations[i].strings[j].antennas[k].ray_sol_cnt; m++) {   // loop over raysol numbers
							//signal_bin.push_back( (stations[i].strings[j].antennas[k].arrival_time[m] - stations[i].min_arrival_time)/(settings1->TIMESTEP) + settings1->NFOUR/4 );
							//signal_bin.push_back( (stations[i].strings[j].antennas[k].arrival_time[m] - stations[i].min_arrival_time)/(settings1->TIMESTEP) + settings1->NFOUR/2 + trigger->maxt_diode_bin );
							signal_bin.push_back(
									(stations[i].strings[j].antennas[k].arrival_time[m]
											- stations[i].min_arrival_time)
											/ (settings1->TIMESTEP)
											+ settings1->NFOUR * 2
											+ trigger->maxt_diode_bin);
							//signal_bin.push_back( (stations[i].strings[j].antennas[k].arrival_time[m] - stations[i].min_arrival_time)/(settings1->TIMESTEP) + settings1->NFOUR + trigger->maxt_diode_bin );

							// store signal located bin
							stations[i].strings[j].antennas[k].SignalBin.push_back(
									signal_bin[m]);

							if (m > 0) {
								signal_dbin.push_back(
										signal_bin[m] - signal_bin[m - 1]);
								//cout<<"  signal_dbin["<<m-1<<"] : "<<signal_dbin[m-1];
								if (signal_dbin[m - 1] < settings1->NFOUR / 2) { // if two ray_sol time delay is smaller than time window
									connect_signals.push_back(1);
									//cout<<"need two signal connection!!"<<endl;
								} else {
									connect_signals.push_back(0);
								}
							} else if (stations[i].strings[j].antennas[k].ray_sol_cnt
									== 1 && m == 0) { // if only one solution
								connect_signals.push_back(0);
							}
							//cout<<"\n";

						}
						//

						//cout<<"done calculating signal bins / connect or not"<<endl;

						// grap noise waveform (NFOUR/2 bins or NFOUR) for diode convlv
						for (int m = 0; m < stations[i].strings[j].antennas[k].ray_sol_cnt; m++) {   // loop over raysol numbers
							// when ray_sol_cnt == 0, this loop inside codes will not run

							if (m == 0) {    // if it's first sol

								if (connect_signals[m] == 1) {

									// do two convlv with double array m, m+1
									//

									Select_Wave_Convlv_Exchange(settings1,
											trigger, detector, signal_bin[m],
											signal_bin[m + 1],
											stations[i].strings[j].antennas[k].V[m],
											stations[i].strings[j].antennas[k].V[m
													+ 1], noise_ID, ch_ID, i);
								} else if (connect_signals[m] == 0) {
									// cout << noise_ID << " : " << ch_ID << " : " << i << " : " << endl;
									// do NFOUR/2 size array convlv (m)
									//
									Select_Wave_Convlv_Exchange(settings1,
											trigger, detector, signal_bin[m],
											stations[i].strings[j].antennas[k].V[m],
											noise_ID, ch_ID, i);

								}
							}

							else {   // if it's not the first sol

								if (m + 1 < stations[i].strings[j].antennas[k].ray_sol_cnt) { // if there is next raysol

									if (connect_signals[m] == 1) { // next raysol is connected

										if (connect_signals[m - 1] == 1) { // and previous raysol also connected

											// double size array with m-1, m, m+1 raysols all added
											//
											Select_Wave_Convlv_Exchange(
													settings1, trigger,
													detector, signal_bin[m - 1],
													signal_bin[m],
													signal_bin[m + 1],
													stations[i].strings[j].antennas[k].V[m
															- 1],
													stations[i].strings[j].antennas[k].V[m],
													stations[i].strings[j].antennas[k].V[m
															+ 1], noise_ID,
													ch_ID, i);
										}

										else if (connect_signals[m - 1] == 0) { // and previous raysol not connected

											// double size array with m , m+1 raysols
											//
											Select_Wave_Convlv_Exchange(
													settings1, trigger,
													detector, signal_bin[m],
													signal_bin[m + 1],
													stations[i].strings[j].antennas[k].V[m],
													stations[i].strings[j].antennas[k].V[m
															+ 1], noise_ID,
													ch_ID, i);

										}
									} else if (connect_signals[m] == 0) { // next raysol is not connected

										if (connect_signals[m - 1] == 1) { // and previous raysol is connected

											// skip the process as this should have done before
											//

										}

										else if (connect_signals[m - 1] == 0) { // and previous raysol not connected

											// single size array with only m raysol
											//
											Select_Wave_Convlv_Exchange(
													settings1, trigger,
													detector, signal_bin[m],
													stations[i].strings[j].antennas[k].V[m],
													noise_ID, ch_ID, i);

										}
									}
								}

								else { // there is no next raysol (this "m" is the last raysol)

									if (connect_signals[m - 1] == 1) { // and previous raysol is connected

										// skip the process as this should have done before
										//

									} else if (connect_signals[m - 1] == 0) { // and previous raysol is not connected

										// single size array with only m raysol
										//
										Select_Wave_Convlv_Exchange(settings1,
												trigger, detector,
												signal_bin[m],
												stations[i].strings[j].antennas[k].V[m],
												noise_ID, ch_ID, i);

									}
								}

							} // if not the first raysol (all other raysols)

						} // end loop over raysols

						//cout<<"done convlv for signal + noise"<<endl;
						//
						//
						// I think I have to apply gain difference factors in here...
						//
						if (settings1->USE_MANUAL_GAINOFFSET == 1)
							Apply_Gain_Offset(settings1, trigger, detector,
									ch_ID, i); // last i for stationID

					} // if it's not debugmode

					ch_ID++; // now to next channel

				} // for antennas

			} // for strings

			// do only if it's not in debugmode
			if (debugmode == 0) {

				//
				// before we move to next station, do trigger check here!!!
				//

				int trig_i, trig_j, trig_bin;
				int trig_search_init;

				// parts that are added for fixed non-trigger passed chs' V_mimic (fixed V_mimic)
				int last_trig_bin;   // stores last trigger passed bin number
				// mode select for non-trigger passed chs' V_mimic
				int V_mimic_mode = settings1->V_MIMIC_MODE;
				// 0 for orginal style (save the middle of trig_window)
				// 1 for saving waveform starting from last_trig_bin
				// 2 for saving waveform where last_trig_bin located on the middle of the waveform
				//
				//

				int trig_mode = settings1->TRIG_MODE;
				// global trigger mode
				// 0 for orginal N_TRIG out of 16 channels
				// 1 for new stations, N_TRIG_V out of Vpol channels or N_TRIG_H out of Hpol channels

				int check_ch;

				//trig_i = trigger->maxt_diode_bin;
				//trig_i = trigger->maxt_diode_bin + settings1->NFOUR; // give some time shift for mimicing force trig events
				trig_search_init = trigger->maxt_diode_bin + settings1->NFOUR; // give some time shift for mimicing force trig events
				trig_i = trig_search_init;

				// save trig search bin info default (if no global trig, this value saved)
				stations[i].total_trig_search_bin = max_total_bin
						- trig_search_init;

				if (settings1->TRIG_SCAN_MODE == 0) { // ******************** old mode left as-is ********************

					// avoid really long trig_window_bin case (change trig_window to check upto max_total_bin)
					if (max_total_bin - trig_window_bin <= trig_i)
						trig_window_bin = max_total_bin - trig_i - 1;

					//while (trig_i < settings1->DATA_BIN_SIZE - trig_window_bin ) {
					while (trig_i < max_total_bin - trig_window_bin) {

						N_pass = 0;
						N_pass_V = 0;
						N_pass_H = 0;
						last_trig_bin = 0;
						Passed_chs.clear();

						//for ( trig_j=0; trig_j<ch_ID; trig_j++) {    // loop over all channels
						trig_j = 0;
						while (trig_j < ch_ID) {

							int string_i = detector->getStringfromArbAntID(i,
									trig_j);
							int antenna_i = detector->getAntennafromArbAntID(i,
									trig_j);

							int channel_num =
									detector->GetChannelfromStringAntenna(i,
											string_i, antenna_i, settings1);
							// check if we want to use BH chs only for trigger analysis
							//if (settings1->TRIG_ONLY_BH_ON == 1) {
							if ((settings1->TRIG_ONLY_BH_ON == 1)
									&& (settings1->DETECTOR == 3)) { // trig by BH is only for TestBed case

								// check if this channel is BH ch (DAQchan)
								if (detector->stations[i].strings[string_i].antennas[antenna_i].DAQchan
										== 0) {

									trig_bin = 0;
									while (trig_bin < trig_window_bin) {

										//cout<<"trig_bin : "<<trig_bin<<endl;

										if (settings1->NOISE_CHANNEL_MODE
												== 0) {
											// with threshold offset by chs
											if (trigger->Full_window[trig_j][trig_i
													+ trig_bin]
													< (detector->GetThres(i,
															channel_num - 1,
															settings1)
															* trigger->rmsdiode
															* detector->GetThresOffset(
																	i,
																	channel_num
																			- 1,
																	settings1))) { // if this channel passed the trigger!
												//cout<<"trigger passed at bin "<<trig_i+trig_bin<<" ch : "<<trig_j<<endl;
												//stations[i].strings[(int)((trig_j)/4)].antennas[(int)((trig_j)%4)].Trig_Pass = trig_i+trig_bin;
												stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
														trig_i + trig_bin;
												N_pass++;
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 0) { // Vpol
													N_pass_V++;
												}
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 1) { // Hpol
													N_pass_H++;
												}
												if (last_trig_bin
														< trig_i + trig_bin)
													last_trig_bin = trig_i
															+ trig_bin; // added for fixed V_mimic
												trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
												Passed_chs.push_back(trig_j);
											}
										} else if (settings1->NOISE_CHANNEL_MODE
												== 1) {
											// with threshold offset by chs
											if (trigger->Full_window[trig_j][trig_i
													+ trig_bin]
													< (detector->GetThres(i,
															channel_num - 1,
															settings1)
															* trigger->rmsdiode_ch[channel_num
																	- 1]
															* detector->GetThresOffset(
																	i,
																	channel_num
																			- 1,
																	settings1))) { // if this channel passed the trigger!
												//cout<<"trigger passed at bin "<<trig_i+trig_bin<<" ch : "<<trig_j<<endl;
												//stations[i].strings[(int)((trig_j)/4)].antennas[(int)((trig_j)%4)].Trig_Pass = trig_i+trig_bin;
												stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
														trig_i + trig_bin;
												N_pass++;
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 0) { // Vpol
													N_pass_V++;
												}
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 1) { // Hpol
													N_pass_H++;
												}
												if (last_trig_bin
														< trig_i + trig_bin)
													last_trig_bin = trig_i
															+ trig_bin; // added for fixed V_mimic
												trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
												Passed_chs.push_back(trig_j);
											}
										} else if (settings1->NOISE_CHANNEL_MODE
												== 2) {
											// with threshold offset by chs
											// for TRIG_ONLY_BH_ON = 1 case, we are only using first 8 chs so don't worry about other chs
											if (trigger->Full_window[trig_j][trig_i
													+ trig_bin]
													< (detector->GetThres(i,
															channel_num - 1,
															settings1)
															* trigger->rmsdiode_ch[channel_num
																	- 1]
															* detector->GetThresOffset(
																	i,
																	channel_num
																			- 1,
																	settings1))) { // if this channel passed the trigger!
												//cout<<"trigger passed at bin "<<trig_i+trig_bin<<" ch : "<<trig_j<<endl;
												//stations[i].strings[(int)((trig_j)/4)].antennas[(int)((trig_j)%4)].Trig_Pass = trig_i+trig_bin;
												stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
														trig_i + trig_bin;
												N_pass++;
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 0) { // Vpol
													N_pass_V++;
												}
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 1) { // Hpol
													N_pass_H++;
												}
												if (last_trig_bin
														< trig_i + trig_bin)
													last_trig_bin = trig_i
															+ trig_bin; // added for fixed V_mimic
												trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
												Passed_chs.push_back(trig_j);
											}
										}

										trig_bin++;

									}
								}

							}

							// in this case just use first 8 chs' thres values
							else if ((settings1->TRIG_ONLY_LOW_CH_ON == 1)
									&& (settings1->DETECTOR != 3)) { // non-TestBed case, and only trig by lower 8 channels

								// reset channel numbers so that bottom antennas have ch 1-8
								channel_num = GetChannelNum8_LowAnt(string_i,
										antenna_i);

								if (antenna_i < 2) { // only antenna 0, 1 which are bottom 2 antennas

									// set channel_num as new value (antenna 0, 1 are only possible antennas for channel_num 1 - 8)

									trig_bin = 0;
									while (trig_bin < trig_window_bin) {

										if (settings1->NOISE_CHANNEL_MODE
												== 0) {
											// with threshold offset by chs
											if (trigger->Full_window[trig_j][trig_i
													+ trig_bin]
													< (detector->GetThres(i,
															channel_num - 1,
															settings1)
															* trigger->rmsdiode
															* detector->GetThresOffset(
																	i,
																	channel_num
																			- 1,
																	settings1))) { // if this channel passed the trigger!
												//cout<<"trigger passed at bin "<<trig_i+trig_bin<<" ch : "<<trig_j<<endl;
												//stations[i].strings[(int)((trig_j)/4)].antennas[(int)((trig_j)%4)].Trig_Pass = trig_i+trig_bin;
												stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
														trig_i + trig_bin;
												N_pass++;
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 0) { // Vpol
													N_pass_V++;
												}
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 1) { // Hpol
													N_pass_H++;
												}
												if (last_trig_bin
														< trig_i + trig_bin)
													last_trig_bin = trig_i
															+ trig_bin; // added for fixed V_mimic
												trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
												Passed_chs.push_back(trig_j);
											}
										} else if (settings1->NOISE_CHANNEL_MODE
												== 1) {
											// with threshold offset by chs
											if (trigger->Full_window[trig_j][trig_i
													+ trig_bin]
													< (detector->GetThres(i,
															channel_num - 1,
															settings1)
															* trigger->rmsdiode_ch[channel_num
																	- 1]
															* detector->GetThresOffset(
																	i,
																	channel_num
																			- 1,
																	settings1))) { // if this channel passed the trigger!
												stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
														trig_i + trig_bin;
												N_pass++;
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 0) { // Vpol
													N_pass_V++;
												}
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 1) { // Hpol
													N_pass_H++;
												}
												if (last_trig_bin
														< trig_i + trig_bin)
													last_trig_bin = trig_i
															+ trig_bin; // added for fixed V_mimic
												trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
												Passed_chs.push_back(trig_j);
											}
										} else if (settings1->NOISE_CHANNEL_MODE
												== 2) {
											// with threshold offset by chs
											// for TRIG_ONLY_BH_ON = 1 case, we are only using first 8 chs so don't worry about other chs
											if (channel_num - 1 < 8) {
												if (trigger->Full_window[trig_j][trig_i
														+ trig_bin]
														< (detector->GetThres(i,
																channel_num - 1,
																settings1)
																* trigger->rmsdiode_ch[channel_num
																		- 1]
																* detector->GetThresOffset(
																		i,
																		channel_num
																				- 1,
																		settings1))) { // if this channel passed the trigger!
													stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
															trig_i + trig_bin;
													N_pass++;
													if (detector->stations[i].strings[string_i].antennas[antenna_i].type
															== 0) { // Vpol
														N_pass_V++;
													}
													if (detector->stations[i].strings[string_i].antennas[antenna_i].type
															== 1) { // Hpol
														N_pass_H++;
													}
													if (last_trig_bin
															< trig_i + trig_bin)
														last_trig_bin = trig_i
																+ trig_bin; // added for fixed V_mimic
													trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
													Passed_chs.push_back(
															trig_j);
												}
											} else { // chs starting from 8 (counted from 0), uses same rmsdiode value
												if (trigger->Full_window[trig_j][trig_i
														+ trig_bin]
														< (detector->GetThres(i,
																channel_num - 1,
																settings1)
																* trigger->rmsdiode_ch[8]
																* detector->GetThresOffset(
																		i,
																		channel_num
																				- 1,
																		settings1))) { // if this channel passed the trigger!
													stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
															trig_i + trig_bin;
													N_pass++;
													if (detector->stations[i].strings[string_i].antennas[antenna_i].type
															== 0) { // Vpol
														N_pass_V++;
													}
													if (detector->stations[i].strings[string_i].antennas[antenna_i].type
															== 1) { // Hpol
														N_pass_H++;
													}
													if (last_trig_bin
															< trig_i + trig_bin)
														last_trig_bin = trig_i
																+ trig_bin; // added for fixed V_mimic
													trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
													Passed_chs.push_back(
															trig_j);
												}
											}
										}

										trig_bin++;

									}
								}

							}

							//else if (settings1->TRIG_ONLY_BH_ON == 0) {
							else { // other cases, use all possible chs for trigger analysis
								trig_bin = 0;
								while (trig_bin < trig_window_bin) {

									if (settings1->NOISE_CHANNEL_MODE == 0) {
										// with threshold offset by chs
										if (trigger->Full_window[trig_j][trig_i
												+ trig_bin]
												< (detector->GetThres(i,
														channel_num - 1,
														settings1)
														* trigger->rmsdiode
														* detector->GetThresOffset(
																i,
																channel_num - 1,
																settings1))) { // if this channel passed the trigger!
											//cout<<"trigger passed at bin "<<trig_i+trig_bin<<" ch : "<<trig_j<<endl;
											//stations[i].strings[(int)((trig_j)/4)].antennas[(int)((trig_j)%4)].Trig_Pass = trig_i+trig_bin;
											stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
													trig_i + trig_bin;
											N_pass++;
											if (detector->stations[i].strings[string_i].antennas[antenna_i].type
													== 0) { // Vpol
												N_pass_V++;
											}
											if (detector->stations[i].strings[string_i].antennas[antenna_i].type
													== 1) { // Hpol
												N_pass_H++;
											}
											if (last_trig_bin
													< trig_i + trig_bin)
												last_trig_bin = trig_i
														+ trig_bin; // added for fixed V_mimic
											trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
											Passed_chs.push_back(trig_j);
										}
									} else if (settings1->NOISE_CHANNEL_MODE
											== 1) {
										// with threshold offset by chs
										if (trigger->Full_window[trig_j][trig_i
												+ trig_bin]
												< (detector->GetThres(i,
														channel_num - 1,
														settings1)
														* trigger->rmsdiode_ch[channel_num
																- 1]
														* detector->GetThresOffset(
																i,
																channel_num - 1,
																settings1))) { // if this channel passed the trigger!
											stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
													trig_i + trig_bin;
											N_pass++;
											if (detector->stations[i].strings[string_i].antennas[antenna_i].type
													== 0) { // Vpol
												N_pass_V++;
											}
											if (detector->stations[i].strings[string_i].antennas[antenna_i].type
													== 1) { // Hpol
												N_pass_H++;
											}
											if (last_trig_bin
													< trig_i + trig_bin)
												last_trig_bin = trig_i
														+ trig_bin; // added for fixed V_mimic
											trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
											Passed_chs.push_back(trig_j);
										}
									} else if (settings1->NOISE_CHANNEL_MODE
											== 2) {
										// with threshold offset by chs
										// for TRIG_ONLY_BH_ON = 1 case, we are only using first 8 chs so don't worry about other chs
										if (channel_num - 1 < 8) {
											if (trigger->Full_window[trig_j][trig_i
													+ trig_bin]
													< (detector->GetThres(i,
															channel_num - 1,
															settings1)
															* trigger->rmsdiode_ch[channel_num
																	- 1]
															* detector->GetThresOffset(
																	i,
																	channel_num
																			- 1,
																	settings1))) { // if this channel passed the trigger!
												stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
														trig_i + trig_bin;
												N_pass++;
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 0) { // Vpol
													N_pass_V++;
												}
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 1) { // Hpol
													N_pass_H++;
												}
												if (last_trig_bin
														< trig_i + trig_bin)
													last_trig_bin = trig_i
															+ trig_bin; // added for fixed V_mimic
												trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
												Passed_chs.push_back(trig_j);
											}
										} else { // chs starting from 8 (counted from 0), uses same rmsdiode value
											if (trigger->Full_window[trig_j][trig_i
													+ trig_bin]
													< (detector->GetThres(i,
															channel_num - 1,
															settings1)
															* trigger->rmsdiode_ch[8]
															* detector->GetThresOffset(
																	i,
																	channel_num
																			- 1,
																	settings1))) { // if this channel passed the trigger!
												stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
														trig_i + trig_bin;
												N_pass++;
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 0) { // Vpol
													N_pass_V++;
												}
												if (detector->stations[i].strings[string_i].antennas[antenna_i].type
														== 1) { // Hpol
													N_pass_H++;
												}
												if (last_trig_bin
														< trig_i + trig_bin)
													last_trig_bin = trig_i
															+ trig_bin; // added for fixed V_mimic
												trig_bin = trig_window_bin; // if confirmed this channel passed the trigger, no need to do rest of bins
												Passed_chs.push_back(trig_j);
											}
										}
									}

									trig_bin++;

								}
							}

							// check all triggered channels not just 3
							//if ( N_pass > 2 ) trig_j += ch_ID;    // if the number of passed channels is 3 or more, no need to check other remaining channels as this station is trigged!
							//else trig_j++;   // if station not passed the trigger, just go to next channel
							trig_j++; // if station not passed the trigger, just go to next channel

						}    // while trig_j < ch_ID

						//if ( N_pass > settings1->N_TRIG-1 ) {  // now as global trigged!! = more or eq to N_TRIG triggered
						if (((trig_mode == 0)
								&& (N_pass > settings1->N_TRIG - 1)) // trig_mode = 0 case!
						||// or
								((trig_mode == 1)
										&& ((N_pass_V > settings1->N_TRIG_V - 1)
												|| (N_pass_H
														> settings1->N_TRIG_H
																- 1))) // trig_mode = 1 case!
								) {

							check_ch = 0;
							//stations[i].Global_Pass = trig_i;
							stations[i].Global_Pass = last_trig_bin; // where actually global trigger occured

							//trig_i = settings1->DATA_BIN_SIZE;    // also if we know this station is trigged, don't need to check rest of time window
							trig_i = max_total_bin; // also if we know this station is trigged, don't need to check rest of time window
							for (int ch_loop = 0; ch_loop < ch_ID; ch_loop++) {
								//cout << ch_loop << "/" << ch_ID << endl;
								int string_i = detector->getStringfromArbAntID(
										i, ch_loop);
								int antenna_i =
										detector->getAntennafromArbAntID(i,
												ch_loop);
								//      cout << "string:antenna: " << string_i << " : " << antenna_i << endl;
								stations[i].strings[string_i].antennas[antenna_i].Likely_Sol =
										-1; // no likely init
								if (ch_loop == Passed_chs[check_ch]
										&& check_ch < N_pass) { // added one more condition (check_ch<N_Pass) for bug in vector Passed_chs.clear()???

										// store which ray sol is triggered based on trig time
										//

									int mindBin = 1.e9; // init big value
									int dBin = 0;

									for (int m = 0;
											m
													< stations[i].strings[string_i].antennas[antenna_i].ray_sol_cnt;
											m++) {   // loop over raysol numbers

										if (stations[i].strings[string_i].antennas[antenna_i].SignalExt[m]) {

											dBin =
													abs(
															stations[i].strings[string_i].antennas[antenna_i].SignalBin[m]
																	- stations[i].strings[string_i].antennas[antenna_i].Trig_Pass);

											if (dBin < mindBin) {
												stations[i].strings[string_i].antennas[antenna_i].Likely_Sol =
														m; // store the ray sol number which is minimum difference between Trig_Pass bin
												mindBin = dBin;
											}

										}

									}

									//skip this passed ch as it already has bin info
									//cout<<"trigger passed at bin "<<stations[i].strings[(int)((ch_loop)/4)].antennas[(int)((ch_loop)%4)].Trig_Pass<<"  passed ch : "<<ch_loop<<" Direct dist btw posnu : "<<event->Nu_Interaction[0].posnu.Distance( detector->stations[i].strings[(int)((ch_loop)/4)].antennas[(int)((ch_loop)%4)] )<<endl;
									//cout<<endl<<"trigger passed at bin "<<stations[i].strings[string_i].antennas[antenna_i].Trig_Pass<<"  passed ch : "<<ch_loop<<" ("<<detector->stations[i].strings[string_i].antennas[antenna_i].type<<"type) Direct dist btw posnu : "<<event->Nu_Interaction[0].posnu.Distance( detector->stations[i].strings[string_i].antennas[antenna_i] )<<" noiseID : "<<stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];
									//cout<<endl<<"trigger passed at bin "<<stations[i].strings[string_i].antennas[antenna_i].Trig_Pass<<"  passed ch : "<<ch_loop<<" ("<<detector->stations[i].strings[string_i].antennas[antenna_i].type<<"type) Direct dist btw posnu : "<<event->Nu_Interaction[0].posnu.Distance( detector->stations[i].strings[string_i].antennas[antenna_i] )<<" noiseID : "<<stations[i].strings[string_i].antennas[antenna_i].noise_ID[0]<<" ViewAngle : "<<stations[i].strings[string_i].antennas[antenna_i].view_ang[0]*DEGRAD;

									if (settings1->TRIG_ONLY_LOW_CH_ON == 0) {
										cout << endl << "trigger passed at bin "
												<< stations[i].strings[string_i].antennas[antenna_i].Trig_Pass
												<< "  passed ch : " << ch_loop
												<< " ("
												<< detector->stations[i].strings[string_i].antennas[antenna_i].type
												<< "type) Direct dist btw posnu : "
												<< event->Nu_Interaction[0].posnu.Distance(
														detector->stations[i].strings[string_i].antennas[antenna_i])
												<< " noiseID : "
												<< stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];
										if (stations[i].strings[string_i].antennas[antenna_i].Likely_Sol
												!= -1) {
											cout << " ViewAngle : "
													<< stations[i].strings[string_i].antennas[antenna_i].view_ang[0]
															* DEGRAD
													<< " LikelyTrigSignal : "
													<< stations[i].strings[string_i].antennas[antenna_i].Likely_Sol;
										}
									}

									else if (settings1->TRIG_ONLY_LOW_CH_ON
											== 1) {
										cout << endl << "trigger passed at bin "
												<< stations[i].strings[string_i].antennas[antenna_i].Trig_Pass
												<< "  passed ant: str["
												<< string_i << "].ant["
												<< antenna_i << "] ("
												<< detector->stations[i].strings[string_i].antennas[antenna_i].type
												<< "type) Direct dist btw posnu : "
												<< event->Nu_Interaction[0].posnu.Distance(
														detector->stations[i].strings[string_i].antennas[antenna_i])
												<< " noiseID : "
												<< stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];
										if (stations[i].strings[string_i].antennas[antenna_i].Likely_Sol
												!= -1) {
											cout << " ViewAngle : "
													<< stations[i].strings[string_i].antennas[antenna_i].view_ang[0]
															* DEGRAD
													<< " LikelyTrigSignal : "
													<< stations[i].strings[string_i].antennas[antenna_i].Likely_Sol;
										}
									}

									//cout << "True station number: " << detector->InstalledStations[0].VHChannel[string_i][antenna_i] << endl;
									//cout << event->Nu_Interaction[0].posnu[0] << " : " <<  event->Nu_Interaction[0].posnu[1] << " : " << event->Nu_Interaction[0].posnu[2] << endl;
									check_ch++;

									// now save the voltage waveform to V_mimic
									//

									for (int mimicbin = 0;
											mimicbin < waveformLength;
											mimicbin++) {

										// new DAQ waveform writing mechanism test
										if (V_mimic_mode == 0) { // Global passed bin is the center of the window
											stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back(
													(trigger->Full_window_V[ch_loop][last_trig_bin
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin])
															* 1.e3); // save in mV
											stations[i].strings[string_i].antennas[antenna_i].time.push_back(
													last_trig_bin
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin);
											stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back(
													(waveformCenter
															- waveformLength / 2
															+ mimicbin)
															* settings1->TIMESTEP
															* 1.e9); // save in ns
										} else if (V_mimic_mode == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
											stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back(
													(trigger->Full_window_V[ch_loop][last_trig_bin
															- (detector->params.TestBed_Ch_delay_bin[ch_loop]
																	- detector->params.TestBed_BH_Mean_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin])
															* 1.e3); // save in mV
											stations[i].strings[string_i].antennas[antenna_i].time.push_back(
													last_trig_bin
															- (detector->params.TestBed_Ch_delay_bin[ch_loop]
																	- detector->params.TestBed_BH_Mean_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin);
											stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back(
													(-(detector->params.TestBed_Ch_delay_bin[ch_loop]
															- detector->params.TestBed_BH_Mean_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin)
															* settings1->TIMESTEP
															* 1.e9); // save in ns
										} else if (V_mimic_mode == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
											stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back(
													(trigger->Full_window_V[ch_loop][last_trig_bin
															- (detector->params.TestBed_Ch_delay_bin[ch_loop]
																	- detector->params.TestBed_BH_Mean_delay_bin
																	+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin])
															* 1.e3); // save in mV
											stations[i].strings[string_i].antennas[antenna_i].time.push_back(
													last_trig_bin
															- (detector->params.TestBed_Ch_delay_bin[ch_loop]
																	- detector->params.TestBed_BH_Mean_delay_bin
																	+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin);
											stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back(
													(-(detector->params.TestBed_Ch_delay_bin[ch_loop]
															- detector->params.TestBed_BH_Mean_delay_bin
															+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin)
															* settings1->TIMESTEP
															* 1.e9
															+ detector->params.TestBed_WFtime_offset_ns); // save in ns
										}
									}

									// set global_trig_bin values
									if (V_mimic_mode == 0) { // Global passed bin is the center of the window
										stations[i].strings[string_i].antennas[antenna_i].global_trig_bin =
												(waveformLength / 2
														- waveformCenter);
									} else if (V_mimic_mode == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
										stations[i].strings[string_i].antennas[antenna_i].global_trig_bin =
												(detector->params.TestBed_Ch_delay_bin[ch_loop]
														- detector->params.TestBed_BH_Mean_delay_bin)
														- waveformCenter
														+ waveformLength / 2;
									} else if (V_mimic_mode == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
										stations[i].strings[string_i].antennas[antenna_i].global_trig_bin =
												(detector->params.TestBed_Ch_delay_bin[ch_loop]
														- detector->params.TestBed_BH_Mean_delay_bin
														+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
														- waveformCenter
														+ waveformLength / 2;
									}

								} else {
									//stations[i].strings[(int)((ch_loop)/4)].antennas[(int)((ch_loop)%4)].Trig_Pass = stations[i].Global_Pass + settings1->NFOUR/4;   // so that global trig is 
									//stations[i].strings[(int)((ch_loop)/4)].antennas[(int)((ch_loop)%4)].Trig_Pass = 0.;
									stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
											0.;

									// new DAQ waveform writing mechanism test
									for (int mimicbin = 0;
											mimicbin < waveformLength;
											mimicbin++) {
										if (V_mimic_mode == 0) { // Global passed bin is the center of the window
											stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back(
													(trigger->Full_window_V[ch_loop][last_trig_bin
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin])
															* 1.e3); // save in mV
											stations[i].strings[string_i].antennas[antenna_i].time.push_back(
													last_trig_bin
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin);
											stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back(
													(waveformCenter
															- waveformLength / 2
															+ mimicbin)
															* settings1->TIMESTEP
															* 1.e9); // save in ns
										} else if (V_mimic_mode == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
											stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back(
													(trigger->Full_window_V[ch_loop][last_trig_bin
															- (detector->params.TestBed_Ch_delay_bin[ch_loop]
																	- detector->params.TestBed_BH_Mean_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin])
															* 1.e3); // save in mV
											stations[i].strings[string_i].antennas[antenna_i].time.push_back(
													last_trig_bin
															- (detector->params.TestBed_Ch_delay_bin[ch_loop]
																	- detector->params.TestBed_BH_Mean_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin);
											stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back(
													(-(detector->params.TestBed_Ch_delay_bin[ch_loop]
															- detector->params.TestBed_BH_Mean_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin)
															* settings1->TIMESTEP
															* 1.e9); // save in ns
										}
										if (V_mimic_mode == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted delay
											//stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back( trigger->Full_window_V[ch_loop][ last_trig_bin - (detector->params.TestBed_Ch_delay_bin[ch_loop]-detector->params.TestBed_BH_Mean_delay_bin) - settings1->NFOUR/4 + mimicbin ] );
											stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back(
													(trigger->Full_window_V[ch_loop][last_trig_bin
															- (detector->params.TestBed_Ch_delay_bin[ch_loop]
																	- detector->params.TestBed_BH_Mean_delay_bin
																	+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin])
															* 1.e3); // save in mV
											stations[i].strings[string_i].antennas[antenna_i].time.push_back(
													last_trig_bin
															- (detector->params.TestBed_Ch_delay_bin[ch_loop]
																	- detector->params.TestBed_BH_Mean_delay_bin
																	+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin);
											stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back(
													(-(detector->params.TestBed_Ch_delay_bin[ch_loop]
															- detector->params.TestBed_BH_Mean_delay_bin
															+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
															+ waveformCenter
															- waveformLength / 2
															+ mimicbin)
															* settings1->TIMESTEP
															* 1.e9
															+ detector->params.TestBed_WFtime_offset_ns); // save in ns
										}
									}

									// set global_trig_bin values
									if (V_mimic_mode == 0) { // Global passed bin is the center of the window
										stations[i].strings[string_i].antennas[antenna_i].global_trig_bin =
												waveformLength / 2
														- waveformCenter;
									} else if (V_mimic_mode == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
										stations[i].strings[string_i].antennas[antenna_i].global_trig_bin =
												(detector->params.TestBed_Ch_delay_bin[ch_loop]
														- detector->params.TestBed_BH_Mean_delay_bin)
														+ waveformLength / 2
														- waveformCenter;
									} else if (V_mimic_mode == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
										stations[i].strings[string_i].antennas[antenna_i].global_trig_bin =
												(detector->params.TestBed_Ch_delay_bin[ch_loop]
														- detector->params.TestBed_BH_Mean_delay_bin
														+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
														+ waveformLength / 2
														- waveformCenter;
									}

									// done V_mimic for non-triggered chs ( done fixed V_mimic )
								}

								double arrivtime =
										stations[i].strings[string_i].antennas[antenna_i].arrival_time[0];
								double X =
										detector->stations[i].strings[string_i].antennas[antenna_i].GetX();
								double Y =
										detector->stations[i].strings[string_i].antennas[antenna_i].GetY();
								double Z =
										detector->stations[i].strings[string_i].antennas[antenna_i].GetZ();
								//std::cout << "Arrival time:X:Y:Z " << arrivtime << " : " << X << " : " << Y << " : " << Z << std::endl;
								/*
								 int AraRootChannel = 0;
								 AraRootChannel = detector->GetChannelfromStringAntenna (i, string_i, antenna_i, settings1);

								 int UsefulEventBin;
								 if ( settings1->NFOUR/2 < EFFECTIVE_LAB3_SAMPLES*2) UsefulEventBin = settings1->NFOUR/2;
								 else UsefulEventBin = EFFECTIVE_LAB3_SAMPLES*2;

								 //for (int mimicbin=0; mimicbin<settings1->NFOUR/2; mimicbin++) {
								 for (int mimicbin=0; mimicbin<UsefulEventBin; mimicbin++) {
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
							stations[i].total_trig_search_bin =
									stations[i].Global_Pass + trig_window_bin
											- trig_search_init;

						} // if global trig!
						else {
							trig_i++; // also if station not passed the trigger, just go to next bin
							/*
							 for (int ch_loop=0; ch_loop<ch_ID; ch_loop++) {
							 int string_i = detector->getStringfromArbAntID( i, ch_loop);
							 int antenna_i = detector->getAntennafromArbAntID( i, ch_loop);
							 int AraRootChannel = 0;
							 AraRootChannel = detector->GetChannelfromStringAntenna (i, string_i, antenna_i, settings1);

							 int UsefulEventBin;
							 if ( settings1->NFOUR/2 < EFFECTIVE_LAB3_SAMPLES*2) UsefulEventBin = settings1->NFOUR/2;
							 else UsefulEventBin = EFFECTIVE_LAB3_SAMPLES*2;

							 //for (int mimicbin=0; mimicbin<settings1->NFOUR/2; mimicbin++) {
							 for (int mimicbin=0; mimicbin<UsefulEventBin; mimicbin++) {
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
					}    // while trig_i

				}    // if TRIG_SCAN_MODE==0

				else if (settings1->TRIG_SCAN_MODE > 0) {

					if (settings1->TRIG_SCAN_MODE == 9) { //hi-low trigger development

						//printf("RMS voltage: %.10e \n",trigger->rmsvoltage);

						/*
						before now, we have an iterator variable called trig_i which has been set to "trig_search_init"			
						this is where in the trace we should start scanning
						*/

						// cout<<"Max total bin is "<<max_total_bin<<endl;
						// cout<<"Trig window bin "<<trig_window_bin<<endl;

						//avoid really long trig_window_bin case (change trig_window to check up to max_total_bin)
						if(max_total_bin - trig_window_bin <= trig_i){
							trig_window_bin = max_total_bin - trig_i -1;
						}

						/*
						We want to scan over this TRIG_WINDOW worth of time
						and trig_window_bin is this TRIG_WINDOW in bin variables (so TRIG_WINDOW/TIMESTEP)
						so we will advance the trig_i variable, slowly scanning rightward on the trace
						we should stop when there isn't enough trace left for a full trig_window_bin
						ideally (if the previos programmer's were careful) the waveform will have been zero padded
						sufficiently much that this won't be a problem
						*/
						// vector <bool> channel_trig_this_coincidence;
						// channel_trig_this_coincidence.resize(ch_ID);						
						// for(int this_chan=0; this_chan<ch_ID;this_chan++) channel_trig_this_coincidence.push_back(false);

						int num_lpda_trigger = 0;
						int num_dipole_trigger = 0;
						int num_phased_trigger = 0;
						Passed_chs.clear();

						//these might be irrelevant in this triggering scheme; here for heritage
						// N_pass = 0;
						// N_pass_V = 0;
						// N_pass_H = 0;

						while(trig_i < max_total_bin - trig_window_bin){

							// printf("Coincidence Trig window start bin: %d \n", trig_i);

							int trig_k = trig_i;
							int hi_lo_window_bin = int(5.e-9/settings1->TIMESTEP); //this is the width of the hi-lo window
							//printf("Hi low bin size is %d \n",hi_lo_window_bin);

							// trig_j=0; //heritage code, not sure if this will be necessary long term

							//Now we go and loop over all of the samples in the 5ns hi-lo coincidence window of this 50 ns chunk
							// int num_hi_lo=0; //counter for how many channels in this multi-coincidence window have a hi+low
							// vector <int> passing_chans;
						
							//this is to keep track of if at any point in our 5ns bin sliding a channel goes hi
							vector <bool> channel_trig_this_hilo_bin;
							channel_trig_this_hilo_bin.resize(ch_ID);
							for(int this_chan=0; this_chan<ch_ID;this_chan++) channel_trig_this_hilo_bin.push_back(false);

							while(trig_k < trig_i + trig_window_bin){
								// printf("	Hi-Lo window start bin: %d \n",trig_k);

								//loop over channels and check each one for a hi+lo
								int check_chan = 0;
								while(check_chan < ch_ID){

									int string_i = detector->getStringfromArbAntID(i, check_chan);
									int antenna_i = detector->getAntennafromArbAntID(i, check_chan);
									int channel_num = detector->GetChannelfromStringAntenna(i, string_i, antenna_i, settings1);
									int ant_type = detector->stations[i].strings[string_i].antennas[antenna_i].type;
									double X = detector->stations[i].strings[string_i].antennas[antenna_i].GetX();
									double Y = detector->stations[i].strings[string_i].antennas[antenna_i].GetY();
									double Z = detector->stations[i].strings[string_i].antennas[antenna_i].GetZ();
									//printf("		String %d , Antenna %d, Channel Num %d \n",string_i, antenna_i, channel_num);

									bool has_hi=false;
									bool has_lo=false;
									double high_check_val = 100000.; //some crazy high value
									double lo_check_val = -100000.; //some crazy low value
									double noise_level = 0.;
									if(settings1->CUSTOM_ELECTRONICS==1 && settings1->TRIG_ANALYSIS_MODE==1)
									{
										//you should only be using CUSTOM_ELECTRONICS right now if the noise is OFF
										noise_level = 9.3e-6; //noise off
									}
									else if(settings1->CUSTOM_ELECTRONICS==0)
									{
										//no custom electronics, so use the RMS of the generated noise
										noise_level = trigger->rmsvoltage; 
									}

									double lpda_hi = 4.1 * noise_level;
									double lpda_lo = -4.1 * noise_level;
									double dipole_hi = 3. * noise_level;
									double dipole_lo = -3. * noise_level;
									double phased_array_hi = 2. * noise_level;
									double phased_array_lo = -2. * noise_level;
									
									
									//we agreed to make the phased array dipole stations[0].strings[0].antennas[0]
									bool should_i_check_this_chan_for_trigger=false;
									if(string_i==0 && antenna_i==0){
										high_check_val = phased_array_hi;
										lo_check_val = phased_array_lo;
										should_i_check_this_chan_for_trigger=true;
									}
									else if(ant_type==1){ //(lpdas are type==1 in this hacky trigger mode)
										high_check_val = lpda_hi;
										lo_check_val = lpda_lo;
										should_i_check_this_chan_for_trigger=true;
									}
									//it's a bicone and falls into the "shallow bicone" category
									else if(ant_type == 0 && Z < -8.){
										high_check_val = dipole_hi;
										lo_check_val = dipole_lo;
										should_i_check_this_chan_for_trigger=true;
									}
									if(should_i_check_this_chan_for_trigger==false) {check_chan++; continue;} //skip
									// printf("			Chan %d \n",check_chan);

									//Now, we loop over all the samples in this hi_lo_window_bin set
									trig_bin = 0;

									while(trig_bin < hi_lo_window_bin){
										//we're in the trig_i'th multi-channel coincidence spot, in the trig_k'th hi-lo coincidence sub-window, on sample trig_bin 
										//we can check it by looking at trig_k + trig_bin
										int bin_to_check = trig_k + trig_bin;
										// printf("				Trig sample %d, or real sample %d \n", trig_bin, bin_to_check);
										double val= trigger -> Full_window_V[check_chan][bin_to_check]; //get the voltage
										if(val >= high_check_val){
											has_hi=true; //it satisfies the high trigger
										}
										else if( val <= lo_check_val) {
											has_lo=true; //it satisfies the low trigger
										}
										trig_bin++; //advance which bin we're triggering on
									}
									/*
									if in this 5ns window, the channel achieved a hi+low trigger, then for this multi-channel coincidence window
									we get to say we have a channel with a hi+low trigger
									*/
									//see if any channels triggered in this particular 5ns section
									if(has_hi && has_lo){
										channel_trig_this_hilo_bin[check_chan]=true; //mark this channel as being hit
										// printf("		Hi-Lo Trigger!! Chan %d <-------------------\n");
										//this way, if once you've shifted the 5ns over and this channel triggers again, it just keeps getting set hi
										// passing_chans.push_back(check_chan);
									}
									check_chan++;
								}//loop over channels
								trig_k++; //advance the 5ns window
							}//single-channel hi-lo coincidence window

							/*now, for this specific multi-channel coincidence window, we need to see if enough channels
							have gone hi to warrant a global trigger
							so, we need to loop again over all channels
							*/
							int check_chan=0;
							while(check_chan<ch_ID){
								int string_i = detector->getStringfromArbAntID(i, check_chan);
								int antenna_i = detector->getAntennafromArbAntID(i, check_chan);
								int channel_num = detector->GetChannelfromStringAntenna(i, string_i, antenna_i, settings1);
								int ant_type = detector->stations[i].strings[string_i].antennas[antenna_i].type;
								double X = detector->stations[i].strings[string_i].antennas[antenna_i].GetX();
								double Y = detector->stations[i].strings[string_i].antennas[antenna_i].GetY();
								double Z = detector->stations[i].strings[string_i].antennas[antenna_i].GetZ();
								int what_type_am_i = 0;
								
								//we agreed to make the phased array dipole stations[0].strings[0].antennas[0]
								if(string_i==0 && antenna_i==0){
									what_type_am_i=50;
								}
								else if(ant_type==1){//lpda (lpdas are type==1 in this hacky trigger mode)
									what_type_am_i=30;
								}
								//it's a bicone and falls into the "shallow bicone" category
								else if(ant_type == 0 && Z < -8.){
									what_type_am_i=40;
								}

								if(what_type_am_i==30 && channel_trig_this_hilo_bin[check_chan]==true){
									num_lpda_trigger++;
									//stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = 1;
									Passed_chs.push_back(check_chan);
								}
								else if(what_type_am_i==40 && channel_trig_this_hilo_bin[check_chan]==true){
									num_dipole_trigger++;
									//stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = 1;
									Passed_chs.push_back(check_chan);
								}
								else if(what_type_am_i==50 && channel_trig_this_hilo_bin[check_chan]==true){
									num_phased_trigger++;
									//stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = 1;
									Passed_chs.push_back(check_chan);
								}
								check_chan++;							
							}

							//we allow for a global trigger if 2/4 LPDAs, or 4/4 dipoles, or 1 phased array
							if(num_lpda_trigger>1 || num_dipole_trigger==4 || num_phased_trigger>0){
								// if(num_lpda_trigger>1) printf("	It's the LPDA that triggerd! (%d chans) \n",num_lpda_trigger);
								// if(num_dipole_trigger==4) printf("	It's the Dipoles that triggerd! (%d chans)\n",num_dipole_trigger);
								// if(num_phased_trigger>0) printf("	It's the PA that triggered! (%d chans)\n",num_phased_trigger);
								stations[i].Global_Pass = 1; //if this is true in a 50ns window, we've achieved a global trigger!
								// printf("	Full station trigger!!  Coincidence window %d <-------------------------------------------------------\n",trig_i);
								break; //get out, we're done
							}
							// else{
							// 	//otherwise we didn't trigger, and we need to wipe this away for the next round
							// 	Passed_chs.clear();
							// }
							trig_i++; //advance the multi-channel coincidence window
						}//multi-channel coincdience  window loop
					} //mode 9, hi-low development
					else if (settings1->TRIG_SCAN_MODE == 10) {
						N_pass = 0;
						N_pass_V = 0;
						N_pass_H = 0;
						Passed_chs.clear();
						for (int ch_loop = 0; ch_loop < ch_ID; ch_loop++) {
							int string_i = detector->getStringfromArbAntID(i,
									ch_loop);
							int antenna_i = detector->getAntennafromArbAntID(i,
									ch_loop);
							//cout<<"On channel: "<<ch_loop<<", string: "<<string_i<<" antenna: "<<antenna_i<<endl;
							int num_this_ant = 0;
							int type =
									detector->stations[i].strings[string_i].antennas[antenna_i].type; //0 is VPol, 1 is Hpol
							if (type == 0) { //only trigger on Vpol
								for (int sols = 0;
										sols
												< stations[i].strings[string_i].antennas[antenna_i].ray_sol_cnt;
										sols++) {
									if (stations[i].strings[string_i].antennas[antenna_i].PeakV[0]
											> 40e-6) {
										//cout<<"   Channel "<<ch_loop<< "Peak V value is  "<<stations[i].strings[string_i].antennas[antenna_i].PeakV[sols]<<endl;
										num_this_ant++;
									}
								}
							}
							if (num_this_ant > 0) {
								N_pass_V++;
								Passed_chs.push_back(ch_loop);
								stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = 1; //set this channel as passing the trigger
							}
						}
						if (N_pass_V > 2) { //this will have it pass the Vpol trigger
							//cout<<"Pass vpol trigger!"<<endl;
							stations[i].Global_Pass = 1; //it deserves a global trigger
						}
					} //end trigger mode 10
					else {
						cout<<"In this else loop for some reason..>"<<endl;
						triggerCheckLoop(settings1, detector, event, trigger, i,
								trig_search_init, max_total_bin,
								trig_window_bin, settings1->TRIG_SCAN_MODE);
					}
				}
			} // if it's not debugmode

			// delete noise waveforms 
			if (settings1->NOISE_WAVEFORM_GENERATE_MODE == 0) { // noise waveforms will be generated for each evts
				// remove noise waveforms for next evt
				trigger->ClearNoiseWaveforms();
			}

		}           // if there is any ray_sol in the station
		/*
		 // now remove all information which are useless
		 for (int c_j=0; c_j< detector->stations[i].strings.size(); c_j++) {

		 for (int c_k=0; c_k< detector->stations[i].strings[c_j].antennas.size(); c_k++) {

		 stations[i].strings[c_j].antennas[c_k].clear_useless(settings1);  // clear data in antenna which stored in previous event
		 //stations[i].strings[c_j].antennas[c_k].clear();  // clear data in antenna which stored in previous event

		 }
		 }
		 // clear useless done
		 */
	} // for stations

	// also clear all vector info to reduce output root file size
	clear_useless(settings1); // to reduce the size of output AraOut.root, remove some information

}   // end Connect_Interaction_Detector

int Report::triggerCheckLoop(Settings *settings1, Detector *detector, Event *event, Trigger *trigger, int stationID, int trig_search_init, int max_total_bin, int trig_window_bin, int scan_mode) {

	class CircularBuffer {

	public:
		int i;
		int mode; // in mode 1 only check number of values above threshold, in >1 check what the best value is too
		int changelog; // if the best value changed, this is returned by add and fill
		int N; // size of buffer
		double *buffer;
		double pthresh; // general run's pthresh
		double best_value; // best value
		double temp_value; // copy of best value, updates each call (so we can zero it when sorting)
		double epsilon; // best value needs to be this close to current last value
		double last_value; // value leaving buffer
		int addToNPass; // number of values above pthresh inside buffer

		CircularBuffer(int size, double threshold, int scan_mode) :
				N(size), pthresh(threshold), mode(scan_mode) {
			i = 0;
			best_value = 0;
			temp_value = 0;
			last_value = 0;
			addToNPass = 0;
			epsilon = 1e-6;
			buffer = new double[N];
			for (int j = 0; j < N; j++)
				buffer[j] = 0;
		}
		~CircularBuffer() {
			delete[] buffer;
		}

		int add(double input_value) {
			changelog = 0;
			temp_value = best_value;

			if (buffer[i] < pthresh)
				addToNPass--; // if the value leaving the buffer is over threshold, we reduce the counter.
			if (mode > 1)
				last_value = buffer[i]; // in mode 1 we don't care about the values, just about the addToNPass

			if (input_value < pthresh) {
				addToNPass++; // if value entering buffer is over threshold we increase the counter.
				buffer[i] = input_value;
			} else
				buffer[i] = 0; // if under threshold just insert zero

			if (mode > 1 && buffer[i] < best_value) {
				best_value = buffer[i];
				temp_value = best_value;
				changelog = 1;
			} // improve best threshold value in the buffer
			if (mode > 1 && std::fabs(last_value - best_value) < epsilon) { // i.e. if last_value==best_value. this happens when the best value leaves buffer
				best_value = findBestValue(); // rescan whole buffer. this should be rather rare...
				temp_value = best_value;
				changelog = 1;
			}
			i++;
			if (i == N)
				i = 0;
			return changelog; // in mode>1: return value is zero if no changes, and non-zero if there are changes to best_value (for mode==1 it's always zero);
		} // add 

		int fill(double input_value) { // buffer is just filling, don't use last_value
			changelog = 0;
			temp_value = best_value;
			if (input_value < pthresh) {
				addToNPass++; // if value entering buffer is over threshold we increase the counter.
				buffer[i] = input_value;
				// cout<<"addToNPass++\n";
			} else
				buffer[i] = 0; // if under threshold add just zero
			if (mode > 1 && buffer[i] < best_value) {
				best_value = buffer[i];
				temp_value = best_value;
				changelog = 1;
			} // improve best threshold value in the buffer
			i++;
			if (i == N)
				i = 0;
			return changelog;
		} // fill

		double findBestValue() {
			double temp_best = 0;
			for (int ii = 0; ii < N; ii++)
				if (buffer[ii] < temp_best) {
					temp_best = buffer[ii];
				}
			return temp_best;
		} // find best value

		int numBinsToOldestTrigger() {
			if (addToNPass < 1) {
				cerr
						<< "ERROR: there are no triggers in this buffer right now!\n";
				return -1;
			}
			int j; // the number of bins after "i" where the first trigger is found 

			for (j = 1; j < N; j++) {
				// i-1 is b/c we did i++ in the last call to add/fill.
				if (buffer[(i - 1 + j) % N] < 0)
					break; // if there's any value here, its because it passed the threshold, so take it
			} // for j
			return N - j; // the backward count of how many bins between i and the earliest trigger... 
		} // numBinsToOldestTrigger

		int numBinsToLatestTrigger() {
			if (addToNPass < 1) {
				cerr
						<< "ERROR: there are no triggers in this buffer right now!\n";
				return -1;
			}
			int j; // the number of bins after "i" where the first trigger is found 
			for (j = 0; j < N; j++) {
				int bin = i - 1 - j;
				if (bin < 0)
					bin = N + (i - j);
				if (buffer[bin] < 0) { // if there's any value here, its because it passed the threshold, so take it
					break;
				}
			} // for j
			return j; // the count of how many bins between i and the latest trigger... 
		} // numBinsToLatestTrigger

	};
	//end circular buffer class

	int i = stationID;

	int numChan = stations[i].TDR_all.size();
	int numChanVpol = stations[i].TDR_Vpol_sorted.size();
	int numChanHpol = stations[i].TDR_Hpol_sorted.size();

	double powerthreshold = settings1->POWERTHRESHOLD;
	// this value is compared to POWERTHRESHOLD for local trigger.  

	int first_trigger = 0;

	double Pthresh_value[numChan];
	CircularBuffer **buffer = new CircularBuffer*[numChan];
	for (int trig_j = 0; trig_j < numChan; trig_j++) { // initialize Trig_Pass and buffers

		Pthresh_value[trig_j] = 0;
		buffer[trig_j] = new CircularBuffer(trig_window_bin, powerthreshold,
				scan_mode);

		int string_i = detector->getStringfromArbAntID(i, trig_j);
		int antenna_i = detector->getAntennafromArbAntID(i, trig_j);

		stations[i].strings[string_i].antennas[antenna_i].SingleChannelTriggers =
				0;
		stations[i].strings[string_i].antennas[antenna_i].TotalBinsScannedPerChannel =
				0;

	} // for trig_j

	double *TDR_all_sorted_temp;
	double *TDR_Vpol_sorted_temp;
	double *TDR_Hpol_sorted_temp;

	if (scan_mode > 1) {

		TDR_all_sorted_temp = new double[numChan];
		TDR_Vpol_sorted_temp = new double[numChanVpol];
		TDR_Hpol_sorted_temp = new double[numChanHpol];

		// only need to initialize for scan_mode>1
		for (int trig_j = 0; trig_j < numChan; trig_j++)
			TDR_all_sorted_temp[trig_j] = 0;
		for (int trig_j = 0; trig_j < numChanVpol; trig_j++)
			TDR_Vpol_sorted_temp[trig_j] = 0;
		for (int trig_j = 0; trig_j < numChanHpol; trig_j++)
			TDR_Hpol_sorted_temp[trig_j] = 0;

	} // if scan_mode>1

	int global_pass_bit = 0;
	int check_TDR_configuration = 0; // check if we need to reorder our TDR arrays
	int SCTR_cluster_bit[numChan];

	for (int trig_j = 0; trig_j < numChan; trig_j++)
		SCTR_cluster_bit[trig_j] = 0;

	for (int trig_i = trig_search_init; trig_i < max_total_bin; trig_i++) { // scan the different window positions

		// for trigger check:
		int N_pass = 0;
		int N_pass_V = 0;
		int N_pass_H = 0;

		check_TDR_configuration = 0;

		for (int trig_j = 0; trig_j < numChan; trig_j++) {

			int string_i = detector->getStringfromArbAntID(i, trig_j);
			int antenna_i = detector->getAntennafromArbAntID(i, trig_j);

			if (settings1->TRIG_ONLY_BH_ON == 0
					|| (settings1->TRIG_ONLY_BH_ON == 1
							&& settings1->DETECTOR == 3
							&& detector->stations[i].strings[string_i].antennas[antenna_i].DAQchan
									== 0)
					|| (settings1->TRIG_ONLY_LOW_CH_ON == 1
							&& settings1->DETECTOR != 3 && antenna_i < 2)) { // channel filter: choose if to use lower/borehole channels or not
				int channel_num = detector->GetChannelfromStringAntenna(i,
						string_i, antenna_i, settings1);

				// assign Pthresh a value 
				if (settings1->NOISE_CHANNEL_MODE == 0)
					Pthresh_value[trig_j] = trigger->Full_window[trig_j][trig_i]
							/ (trigger->rmsdiode
									* detector->GetThresOffset(i,
											channel_num - 1, settings1));
				if (settings1->NOISE_CHANNEL_MODE == 1)
					Pthresh_value[trig_j] = trigger->Full_window[trig_j][trig_i]
							/ (trigger->rmsdiode_ch[channel_num - 1]
									* detector->GetThresOffset(i,
											channel_num - 1, settings1));
				if (settings1->NOISE_CHANNEL_MODE == 2) {
					if (channel_num - 1 < 8)
						Pthresh_value[trig_j] =
								trigger->Full_window[trig_j][trig_i]
										/ (trigger->rmsdiode_ch[channel_num - 1]
												* detector->GetThresOffset(i,
														channel_num - 1,
														settings1));
					else
						Pthresh_value[trig_j] =
								trigger->Full_window[trig_j][trig_i]
										/ (trigger->rmsdiode_ch[8]
												* detector->GetThresOffset(i,
														channel_num - 1,
														settings1));
				}
				// this is to count how many local trigger clusters there are 
				if (Pthresh_value[trig_j] < powerthreshold) {

					if (SCTR_cluster_bit[trig_j] == 0)
						stations[i].strings[string_i].antennas[antenna_i].SingleChannelTriggers++;

					// records all the different Pthresh values that caused local trigger.
					if (settings1->TRIG_SCAN_MODE > 2) {

						if (SCTR_cluster_bit[trig_j] == 0) { // if first trigger in cluster

							stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.push_back(
									Pthresh_value[trig_j]);

						} else { // choose the highest trigger value (most negative) in cluster

							if (Pthresh_value[trig_j]
									< stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.back())
								stations[i].strings[string_i].antennas[antenna_i].SCT_threshold_pass.back() =
										Pthresh_value[trig_j];

						}

					} // trig scan mode > 2

					SCTR_cluster_bit[trig_j] = 1;

				} // if local trigger
				else
					SCTR_cluster_bit[trig_j] = 0; // if no local trigger, set zero to start a new cluster at next local trigger

				// and how many bins scanned
				stations[i].strings[string_i].antennas[antenna_i].TotalBinsScannedPerChannel++;

				// fill the buffers (if any changes occur mark check_TDR_configuration as non-zero)
				if (trig_i < trig_search_init + trig_window_bin)
					check_TDR_configuration += buffer[trig_j]->fill(
							Pthresh_value[trig_j]);
				else
					check_TDR_configuration += buffer[trig_j]->add(
							Pthresh_value[trig_j]);

				if (buffer[trig_j]->addToNPass > 0) { // if there is at least one value above threshold in the buffer, this is ++

					N_pass++;
					if (detector->stations[i].strings[string_i].antennas[antenna_i].type
							== 0)
						N_pass_V++;
					if (detector->stations[i].strings[string_i].antennas[antenna_i].type
							== 1)
						N_pass_H++;

   // if(last_trig_bin<trig_i) last_trig_bin=trig_i;// not sure we need this variable in this mode... 

				} // if addToNPass>0

			} // non-testbed case (i.e. use all channels)

		} // for trig_j (channel scan)

		// check if global trigger...

		if ((settings1->TRIG_MODE == 0 && (N_pass >= settings1->N_TRIG))
				|| (settings1->TRIG_MODE == 1
						&& (N_pass_V >= settings1->N_TRIG_V
								|| N_pass_H >= settings1->N_TRIG_H))
				//|| (settings1->TRIG_MODE==2&&( N_pass_0 >= settings1->N_TRIG_0 || N_pass_1 >= settings1->N_TRIG_1 ))
				) { // if there's a trigger !

			global_pass_bit = 1;
			if (first_trigger == 0) { // if this is the first trigger, mark this position and save event

				first_trigger = 1;

				for (int trig_j = 0; trig_j < numChan; trig_j++) {

					int string_i = detector->getStringfromArbAntID(i, trig_j);
					int antenna_i = detector->getAntennafromArbAntID(i, trig_j);
   // if(buffer[trig_j]->addToNPass>0) stations[i].strings[string_i].antennas[antenna_i].Trig_Pass = trig_i-buffer[trig_j]->numBinsToLatestTrigger(); // mark the bin on which we triggered...
					if (buffer[trig_j]->addToNPass > 0)
						stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
								trig_i
										- buffer[trig_j]->numBinsToOldestTrigger(); // mark the bin on which we triggered...
					else
						stations[i].strings[string_i].antennas[antenna_i].Trig_Pass =
								0.;

				} // for trig_j

				saveTriggeredEvent(settings1, detector, event, trigger,
						stationID, trig_search_init, max_total_bin,
						trig_window_bin, trig_i);
 // cout<<"\nPthresh value=";
 // if(scan_mode==1) for(int trig_j=0;trig_j<numChan; trig_j++) cout<<" "<<Pthresh_value[trig_j];
 // if(scan_mode>1)  for(int trig_j=0;trig_j<numChan; trig_j++) cout<<" "<<buffer[trig_j]->best_value;
 // cout<<"\n";

			} // first trigger

	  // if(scan_mode==1) return last_trig_bin; //  if we aren't going to scan all the Pthresh values, just return
			if (scan_mode == 1)
				return trig_i; //  if we aren't going to scan all the Pthresh values, just return
		} else
			global_pass_bit = 0; // no trigger

		if (scan_mode > 1 && check_TDR_configuration && global_pass_bit) { // if there's a trigger and anything changes in the buffers, restock the TDR arrays

			for (int trig_j = 0; trig_j < numChan; trig_j++)
				TDR_all_sorted_temp[trig_j] = 0;
			for (int trig_j = 0; trig_j < numChanVpol; trig_j++)
				TDR_Vpol_sorted_temp[trig_j] = 0;
			for (int trig_j = 0; trig_j < numChanHpol; trig_j++)
				TDR_Hpol_sorted_temp[trig_j] = 0;

			for (int trig_j = 0; trig_j < numChan; trig_j++) { // fill the TDR (unsorted) arrays if they improved... 

				if (buffer[trig_j]->best_value < stations[i].TDR_all[trig_j])
					stations[i].TDR_all[trig_j] = buffer[trig_j]->best_value;

			} // for trig_j

			if (settings1->TRIG_MODE == 0) { // for N out of 16 mode

				if (N_pass >= settings1->N_TRIG)
					for (int ii = 0; ii < N_pass; ii++) { // find the N_pass best channel's TDR and store them.

						double best_thresh = 0;
						int best_chan = 0;

						for (int trig_j = 0; trig_j < numChan; trig_j++)
							if (buffer[trig_j]->temp_value < best_thresh) {
								best_thresh = buffer[trig_j]->temp_value;
								best_chan = trig_j;
							}

						buffer[best_chan]->temp_value = 0;

						TDR_all_sorted_temp[ii] = best_thresh;

					} // for ii

					  // debug output:
				if (TDR_all_sorted_temp[0] > TDR_all_sorted_temp[1]
						|| TDR_all_sorted_temp[1] > TDR_all_sorted_temp[2]) {

					cout << "\n";
					for (int p = 0; p < 80; p++)
						cout << "*";
					cout << "\n  ordering problem: " << TDR_all_sorted_temp[0]
							<< " " << TDR_all_sorted_temp[1] << " "
							<< TDR_all_sorted_temp[2] << "\n";
					for (int p = 0; p < 80; p++)
						cout << "*";
					cout << "\n";

				}

			} // if trig_mode==0
			if (settings1->TRIG_MODE == 1) { // for N out of either polarization

				// for Vpol only:
				if (N_pass_V >= settings1->N_TRIG_V)
					for (int ii = 0; ii < N_pass_V; ii++) { // find the N_pass best channel's TDR and store them.

						double best_thresh = 0;
						int best_chan = 0;

						for (int trig_j = 0; trig_j < numChan; trig_j++) {

							int string_i = detector->getStringfromArbAntID(i,
									trig_j);
							int antenna_i = detector->getAntennafromArbAntID(i,
									trig_j);

							if (detector->stations[i].strings[string_i].antennas[antenna_i].type
									== 0
									&& buffer[trig_j]->temp_value
											< best_thresh) {

								best_thresh = buffer[trig_j]->temp_value;
								best_chan = trig_j;

							} // if best
						} // for trig_j
						buffer[best_chan]->temp_value = 0;

						TDR_Vpol_sorted_temp[ii] = best_thresh;

					} // for ii

			//        debug output:
   // if(TDR_Vpol_sorted_temp[0]>TDR_Vpol_sorted_temp[1]||TDR_Vpol_sorted_temp[1]>TDR_Vpol_sorted_temp[2]){
	
   //   cout<<"\n";
   //   for(int p=0;p<80;p++) cout<<"*";
   //   cout<<"\n  ordering problem, Vpol: "<<TDR_Vpol_sorted_temp[0]<<" "<<TDR_Vpol_sorted_temp[1]<<" "<<TDR_Vpol_sorted_temp[2]<<"\n";
   //   for(int p=0;p<80;p++) cout<<"*";
   //   cout<<"\n";
	 
   // }
			//      for Hpol only
				if (N_pass_H >= settings1->N_TRIG_H)
					for (int ii = 0; ii < N_pass_H; ii++) { // find the N_pass best channel's TDR and store them.

						double best_thresh = 0;
						int best_chan = 0;

						for (int trig_j = 0; trig_j < numChan; trig_j++) {

							int string_i = detector->getStringfromArbAntID(i,
									trig_j);
							int antenna_i = detector->getAntennafromArbAntID(i,
									trig_j);

							if (detector->stations[i].strings[string_i].antennas[antenna_i].type
									== 1
									&& buffer[trig_j]->temp_value
											< best_thresh) {

								best_thresh = buffer[trig_j]->temp_value;
								best_chan = trig_j;

							}  // if best
						}  // for trig_j

						buffer[best_chan]->temp_value = 0;

						TDR_Hpol_sorted_temp[ii] = best_thresh;

					}  // for ii

   //        // debug output:
   // if(TDR_Hpol_sorted_temp[0]>TDR_Hpol_sorted_temp[1]||TDR_Hpol_sorted_temp[1]>TDR_Hpol_sorted_temp[2]){
	
   //   cout<<"\n";
   //   for(int p=0;p<80;p++) cout<<"*";
   //   cout<<"\n  ordering problem, Hpol: "<<TDR_Hpol_sorted_temp[0]<<" "<<TDR_Hpol_sorted_temp[1]<<" "<<TDR_Hpol_sorted_temp[2]<<"\n";
   //   for(int p=0;p<80;p++) cout<<"*";
   //   cout<<"\n";
	 
   // }
			}  // if trig_mode==1

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
			 //    if(TDR_Vpol_sorted_temp[0]>TDR_Vpol_sorted_temp[1]||TDR_Vpol_sorted_temp[1]>TDR_Vpol_sorted_temp[2]){
			 //     
			 //      cout<<"\n";
			 //      for(int p=0;p<80;p++) cout<<"*";
			 //      cout<<"\n  ordering problem, Vpol: "<<TDR_Vpol_sorted_temp[0]<<" "<<TDR_Vpol_sorted_temp[1]<<" "<<TDR_Vpol_sorted_temp[2]<<"\n";
			 //      for(int p=0;p<80;p++) cout<<"*";
			 //      cout<<"\n";
			 //      
			 //    }
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
			if (settings1->TRIG_MODE == 0) {

				if (N_pass >= settings1->N_TRIG)
					for (int ii = 0; ii < N_pass; ii++)
						if (TDR_all_sorted_temp[ii]
								< stations[i].TDR_all_sorted[ii])
							stations[i].TDR_all_sorted[ii] =
									TDR_all_sorted_temp[ii];

			}
			if (settings1->TRIG_MODE == 1) {

				if (N_pass_V >= settings1->N_TRIG_V)
					for (int ii = 0; ii < N_pass_V; ii++)
						if (TDR_Vpol_sorted_temp[ii]
								< stations[i].TDR_Vpol_sorted[ii])
							stations[i].TDR_Vpol_sorted[ii] =
									TDR_Vpol_sorted_temp[ii];
				if (N_pass_H >= settings1->N_TRIG_H)
					for (int ii = 0; ii < N_pass_H; ii++)
						if (TDR_Hpol_sorted_temp[ii]
								< stations[i].TDR_Hpol_sorted[ii])
							stations[i].TDR_Hpol_sorted[ii] =
									TDR_Hpol_sorted_temp[ii];

				// for this mode, can get TDR_all_sorted from these two arrays:
			}

			/*
			 if(settings1->TRIG_MODE==2){
			 
			 if(N_pass_0>=settings1->N_TRIG_0) for(int ii=0;ii<N_pass_0;ii++) if(TDR_Vpol_sorted_temp[ii]<stations[i].TDR_Vpol_sorted[ii]) stations[i].TDR_Vpol_sorted[ii]=TDR_Vpol_sorted_temp[ii];
			 if(N_pass_1>=settings1->N_TRIG_1) for(int ii=0;ii<N_pass_1;ii++) if(TDR_Hpol_sorted_temp[ii]<stations[i].TDR_Hpol_sorted[ii]) stations[i].TDR_Hpol_sorted[ii]=TDR_Hpol_sorted_temp[ii];
			 
			 // for this mode, can get TDR_all_sorted from these two arrays:
			 }
			 */

		}  // if trigger and buffer changed

	}  // while trig_i

	if (scan_mode > 1 && stations[i].Global_Pass) {

		if (settings1->TRIG_MODE == 0) {

			cout << "\nPthresh best: ";
			for (int ii = 0; ii < 3; ii++)
				cout << " " << stations[i].TDR_all_sorted[ii];
			cout << "\n";

			// debug output:
			if (stations[i].TDR_all_sorted[0] > stations[i].TDR_all_sorted[1]
					|| stations[i].TDR_all_sorted[1]
							> stations[i].TDR_all_sorted[2]) {

				cout << "\n";
				for (int p = 0; p < 80; p++)
					cout << "*";
				cout << "\n  ordering problem: "
						<< stations[i].TDR_all_sorted[0] << " "
						<< stations[i].TDR_all_sorted[1] << " "
						<< stations[i].TDR_all_sorted[2] << "\n";
				for (int p = 0; p < 80; p++)
					cout << "*";
				cout << "\n";

			}  // ordering problem

		}  // trig mode 0

		if (settings1->TRIG_MODE == 1) {
			cout << "\nPthresh best: ";
			cout << "  Vpol: ";
			for (int ii = 0; ii < stations[i].TDR_Vpol_sorted.size(); ii++)
				cout << " " << stations[i].TDR_Vpol_sorted[ii];
			cout << "  Hpol: ";
			for (int ii = 0; ii < stations[i].TDR_Hpol_sorted.size(); ii++)
				cout << " " << stations[i].TDR_Hpol_sorted[ii];
			cout << "\n";

			// debug output:
			if (stations[i].TDR_Vpol_sorted[0] > stations[i].TDR_Vpol_sorted[1]
					|| stations[i].TDR_Vpol_sorted[1]
							> stations[i].TDR_Vpol_sorted[2]) {

				cout << "\n";
				for (int p = 0; p < 80; p++)
					cout << "*";
				cout << "\n  ordering problem (final) Vpol: "
						<< stations[i].TDR_Vpol_sorted[0] << " "
						<< stations[i].TDR_Vpol_sorted[1] << " "
						<< stations[i].TDR_Vpol_sorted[2] << "\n";
				for (int p = 0; p < 80; p++)
					cout << "*";
				cout << "\n";

			}  // ordering problem

			// debug output:
			if (stations[i].TDR_Hpol_sorted[0] > stations[i].TDR_Hpol_sorted[1]
					|| stations[i].TDR_Hpol_sorted[1]
							> stations[i].TDR_Hpol_sorted[2]) {

				cout << "\n";
				for (int p = 0; p < 80; p++)
					cout << "*";
				cout << "\n  ordering problem (final), Hpol: "
						<< stations[i].TDR_Hpol_sorted[0] << " "
						<< stations[i].TDR_Hpol_sorted[1] << " "
						<< stations[i].TDR_Hpol_sorted[2] << "\n";
				for (int p = 0; p < 80; p++)
					cout << "*";
				cout << "\n";

			}  // ordering problem

		}  // trig mode 1

	}  // if scan_mode>1 and global pass

	for (int trig_j = 0; trig_j < numChan; trig_j++)
		delete buffer[trig_j];
	delete[] buffer;

	if (scan_mode > 1) {

		delete[] TDR_all_sorted_temp;
		delete[] TDR_Vpol_sorted_temp;
		delete[] TDR_Hpol_sorted_temp;
	}

	return stations[i].Global_Pass;

}

int Report::saveTriggeredEvent(Settings *settings1, Detector *detector, Event *event, Trigger *trigger, int stationID, int trig_search_init, int max_total_bin, int trig_window_bin, int last_trig_bin) {

	int i = stationID;
	int numChan = stations[i].TDR_all.size();
	cout << "saving event" << endl;
	stations[i].Global_Pass = last_trig_bin;

	int waveformLength = settings1->WAVEFORM_LENGTH;
	int waveformCenter = settings1->WAVEFORM_CENTER;

	for (int trig_j = 0; trig_j < numChan; trig_j++) {

		int string_i = detector->getStringfromArbAntID(i, trig_j);
		int antenna_i = detector->getAntennafromArbAntID(i, trig_j);

		if (stations[i].strings[string_i].antennas[antenna_i].Trig_Pass) { // if this channel triggered

			stations[i].strings[string_i].antennas[antenna_i].Likely_Sol = -1; // no likely init

			int mindBin = 1.e9; // init big values
			int dBin = 0;

			for (int m = 0;
					m
							< stations[i].strings[string_i].antennas[antenna_i].ray_sol_cnt;
					m++) {   // loop over raysol numbers

				if (stations[i].strings[string_i].antennas[antenna_i].SignalExt[m]) {

					dBin =
							abs(
									stations[i].strings[string_i].antennas[antenna_i].SignalBin[m]
											- stations[i].strings[string_i].antennas[antenna_i].Trig_Pass);

					if (dBin < mindBin) {

						stations[i].strings[string_i].antennas[antenna_i].Likely_Sol =
								m; // store the ray sol number which is minimum difference between Trig_Pass bin
						mindBin = dBin;

					}
				}

			} // for m (ray sol numbers)

			if (settings1->TRIG_ONLY_LOW_CH_ON == 0) {

				cout << endl << "trigger passed at bin "
						<< stations[i].strings[string_i].antennas[antenna_i].Trig_Pass
						<< "  passed ch : " << trig_j << " ("
						<< detector->stations[i].strings[string_i].antennas[antenna_i].type
						<< "type) Direct dist btw posnu : "
						<< event->Nu_Interaction[0].posnu.Distance(
								detector->stations[i].strings[string_i].antennas[antenna_i])
						<< " noiseID : "
						<< stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];

				if (stations[i].strings[string_i].antennas[antenna_i].Likely_Sol
						!= -1) {

					cout << " ViewAngle : "
							<< stations[i].strings[string_i].antennas[antenna_i].view_ang[0]
									* DEGRAD << " LikelyTrigSignal : "
							<< stations[i].strings[string_i].antennas[antenna_i].Likely_Sol;

				}
			} else if (settings1->TRIG_ONLY_LOW_CH_ON == 1) {

				cout << endl << "trigger passed at bin "
						<< stations[i].strings[string_i].antennas[antenna_i].Trig_Pass
						<< "  passed ant: str[" << string_i << "].ant["
						<< antenna_i << "] ("
						<< detector->stations[i].strings[string_i].antennas[antenna_i].type
						<< "type) Direct dist btw posnu : "
						<< event->Nu_Interaction[0].posnu.Distance(
								detector->stations[i].strings[string_i].antennas[antenna_i])
						<< " noiseID : "
						<< stations[i].strings[string_i].antennas[antenna_i].noise_ID[0];

				if (stations[i].strings[string_i].antennas[antenna_i].Likely_Sol
						!= -1) {

					cout << " ViewAngle : "
							<< stations[i].strings[string_i].antennas[antenna_i].view_ang[0]
									* DEGRAD << " LikelyTrigSignal : "
							<< stations[i].strings[string_i].antennas[antenna_i].Likely_Sol;

				}

			}

		} // if Trig_Pass

	// now save the voltage waveform to V_mimic
		for (int mimicbin = 0; mimicbin < waveformLength / 2; mimicbin++) {

			// new DAQ waveform writing mechanism test
			if (settings1->V_MIMIC_MODE == 0) { // Global passed bin is the center of the window

				stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back(
						(trigger->Full_window_V[trig_j][last_trig_bin
								- waveformLength / 2 + waveformCenter + mimicbin])
								* 1.e3); // save in mV
				stations[i].strings[string_i].antennas[antenna_i].time.push_back(
						last_trig_bin - waveformLength / 2 + waveformCenter
								+ mimicbin);
				stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back(
						(waveformLength / 2 + waveformCenter + mimicbin)
								* settings1->TIMESTEP * 1.e9); // save in ns
			} else if (settings1->V_MIMIC_MODE == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
				stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back(
						(trigger->Full_window_V[trig_j][last_trig_bin
								- (detector->params.TestBed_Ch_delay_bin[trig_j]
										- detector->params.TestBed_BH_Mean_delay_bin)
								- waveformLength / 2 + waveformCenter + mimicbin])
								* 1.e3); // save in mV
				stations[i].strings[string_i].antennas[antenna_i].time.push_back(
						last_trig_bin
								- (detector->params.TestBed_Ch_delay_bin[trig_j]
										- detector->params.TestBed_BH_Mean_delay_bin)
								- waveformLength / 2 + waveformCenter
								+ mimicbin);
				stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back(
						(-(detector->params.TestBed_Ch_delay_bin[trig_j]
								- detector->params.TestBed_BH_Mean_delay_bin)
								- waveformLength / 2 + waveformCenter + mimicbin)
								* settings1->TIMESTEP * 1.e9); // save in ns

			}

			else if (settings1->V_MIMIC_MODE == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye

				stations[i].strings[string_i].antennas[antenna_i].V_mimic.push_back(
						(trigger->Full_window_V[trig_j][last_trig_bin
								- (detector->params.TestBed_Ch_delay_bin[trig_j]
										- detector->params.TestBed_BH_Mean_delay_bin
										+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
								- waveformLength / 2 + waveformCenter + mimicbin])
								* 1.e3); // save in mV
				stations[i].strings[string_i].antennas[antenna_i].time.push_back(
						last_trig_bin
								- (detector->params.TestBed_Ch_delay_bin[trig_j]
										- detector->params.TestBed_BH_Mean_delay_bin
										+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
								- waveformLength / 2 + waveformCenter
								+ mimicbin);
				stations[i].strings[string_i].antennas[antenna_i].time_mimic.push_back(
						(-(detector->params.TestBed_Ch_delay_bin[trig_j]
								- detector->params.TestBed_BH_Mean_delay_bin
								+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
								- waveformLength / 2 + waveformCenter + mimicbin)
								* settings1->TIMESTEP * 1.e9
								+ detector->params.TestBed_WFtime_offset_ns); // save in ns

			}

		}

		// set global_trig_bin values
		if (settings1->V_MIMIC_MODE == 0) { // Global passed bin is the center of the window
			stations[i].strings[string_i].antennas[antenna_i].global_trig_bin =
					waveformLength / 2 + waveformCenter;
		} else if (settings1->V_MIMIC_MODE == 1) { // Global passed bin is the center of the window + delay to each chs from araGeom
			stations[i].strings[string_i].antennas[antenna_i].global_trig_bin =
					(detector->params.TestBed_Ch_delay_bin[trig_j]
							- detector->params.TestBed_BH_Mean_delay_bin)
							+ waveformLength / 2 + waveformCenter;
		} else if (settings1->V_MIMIC_MODE == 2) { // Global passed bin is the center of the window + delay to each chs from araGeom + fitted by eye
			stations[i].strings[string_i].antennas[antenna_i].global_trig_bin =
					(detector->params.TestBed_Ch_delay_bin[trig_j]
							- detector->params.TestBed_BH_Mean_delay_bin
							+ detector->stations[i].strings[string_i].antennas[antenna_i].manual_delay_bin)
							+ waveformLength / 2 + waveformCenter;
		}

		double arrivtime =
				stations[i].strings[string_i].antennas[antenna_i].arrival_time[0];
		double X =
				detector->stations[i].strings[string_i].antennas[antenna_i].GetX();
		double Y =
				detector->stations[i].strings[string_i].antennas[antenna_i].GetY();
		double Z =
				detector->stations[i].strings[string_i].antennas[antenna_i].GetZ();

		stations[i].total_trig_search_bin = stations[i].Global_Pass
				+ trig_window_bin - trig_search_init;

	} // for trig_j

	if (settings1->OUTPUT_TDR_GRAPH > 0) {

		settings1->OUTPUT_TDR_GRAPH--;

		TGraph **gr = new TGraph*[numChan];

		for (int trig_j = 0; trig_j < numChan; trig_j++) {

			int string_i = detector->getStringfromArbAntID(i, trig_j);
			int antenna_i = detector->getAntennafromArbAntID(i, trig_j);
			int channel_num = detector->GetChannelfromStringAntenna(i, string_i,
					antenna_i, settings1);
			double thresh_value = 0;

			// assign Pthresh a value 
			if (settings1->NOISE_CHANNEL_MODE == 0)
				thresh_value = detector->GetThres(i, channel_num - 1, settings1)
						* trigger->rmsdiode
						* detector->GetThresOffset(i, channel_num - 1,
								settings1);
			if (settings1->NOISE_CHANNEL_MODE == 1)
				thresh_value = detector->GetThres(i, channel_num - 1, settings1)
						* trigger->rmsdiode_ch[channel_num - 1]
						* detector->GetThresOffset(i, channel_num - 1,
								settings1);
			if (settings1->NOISE_CHANNEL_MODE == 2) {

				if (channel_num - 1 < 8)
					thresh_value = detector->GetThres(i, channel_num - 1,
							settings1) * trigger->rmsdiode_ch[channel_num - 1]
							* detector->GetThresOffset(i, channel_num - 1,
									settings1);
				else
					thresh_value = detector->GetThres(i, channel_num - 1,
							settings1) * trigger->rmsdiode_ch[8]
							* detector->GetThresOffset(i, channel_num - 1,
									settings1);

			}

			gr[trig_j] = new TGraph();
			for (int trig_i = 0; trig_i < settings1->DATA_BIN_SIZE / 2;
					trig_i++)
				gr[trig_j]->SetPoint(trig_i, (trig_search_init + trig_i),
						trigger->Full_window[trig_j][trig_i]);
			gr[trig_j]->SetNameTitle(
					Form("TDR_waveform%dC%02d", settings1->OUTPUT_TDR_GRAPH,
							trig_j),
					Form(
							"Tunnel diode response waveform %d, channel %02d, trig_pass= %d, P_{th}= %le; time bins; power after convolution with tunnel diode",
							settings1->OUTPUT_TDR_GRAPH, trig_j,
							stations[i].strings[string_i].antennas[antenna_i].Trig_Pass,
							thresh_value));
			gr[trig_j]->Write();
			delete gr[trig_j];
		}

		delete[] gr;

	}

	return 1;

} // saveTriggeredEvent

#ifdef ARA_UTIL_EXISTS
void Report::MakeUsefulEvent(Detector *detector, Settings *settings1, Trigger *trigger, int stationID, int stationIndex, UsefulIcrrStationEvent *theUsefulEvent) {
	if (stationID < detector->params.number_of_stations) {
		int i = stationID;
		cout << stationID << endl;
		int ch_limit;
		if (stationID == 0) {
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
				if (stations[stationIndex].Global_Pass > 0) {
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
	cout << "StationID: " << stationID << endl;
	theUsefulEvent->fNumChannels = 32;
	theUsefulEvent->stationId = stationID;

	//  cout << endl << stationID << endl;

	int ch_limit;
	if (stationID == 0) {
		ch_limit = 14;
	} else {
		ch_limit = 16;
	}

	int maxElecChans = 32;

	for (int ch_loop=0; ch_loop < ch_limit; ch_loop++) {
		//    int elecChan = AraGeom->getElecChanFromRFChan(ch_loop, stationID);
		int elecChan = AraGeomTool::Instance()->getElecChanFromRFChan(ch_loop, stationID);
		int string_i = 0;
		int antenna_i = 0;
		detector->GetSSAfromChannel(stationID, ch_loop, &antenna_i, &string_i, settings1);

		//    cout << ch_loop << " : " << elecChan << " : " << string_i << " : " << antenna_i << endl;

		//    int string_i = detector->getStringfromArbAntID( stationIndex, ch_loop);
		//    int antenna_i = detector->getAntennafromArbAntID( stationIndex, ch_loop);
		int AraRootChannel = 0;
		AraRootChannel = detector->GetChannelfromStringAntenna (stationID, string_i, antenna_i, settings1);

		int UsefulEventBin;
		//    if ( settings1->NFOUR/2 < EFFECTIVE_LAB3_SAMPLES*2) UsefulEventBin = settings1->NFOUR/2;
		//    else UsefulEventBin = EFFECTIVE_LAB3_SAMPLES*2;

		UsefulEventBin = settings1->WAVEFORM_LENGTH;

		vector < double > volts;
		volts.resize(UsefulEventBin);
		vector < double > times;
		times.resize(UsefulEventBin);

		//    theUsefulEvent->fVolts[AraRootChannel-1].resize(UsefulEventBin);
		//    theUsefulEvent->fTimes[AraRootChannel-1].resize(UsefulEventBin);
		//    cout << string_i << " : " << antenna_i << endl;

		//for (int mimicbin=0; mimicbin<settings1->NFOUR/2; mimicbin++) {
		for (int mimicbin=0; mimicbin<UsefulEventBin; mimicbin++) {
			if (stations[stationIndex].Global_Pass > 0) {
				//        cout << "Test 1" << endl;
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

		//    cout << elecChan << " : " << theUsefulEvent->fTimes[elecChan][1] <<  " : " << stations[stationIndex].strings[string_i].antennas[antenna_i].time_mimic[1] << endl;
		//    cout << elecChan << " : " << theUsefulEvent->fVolts[elecChan][1] <<  " : " << stations[stationIndex].strings[string_i].antennas[antenna_i].V_mimic[1] << endl;

		volts.clear();
		times.clear();

	}
	//  }
}
#endif

void Report::ClearUselessfromConnect(Detector *detector, Settings *settings1, Trigger *trigger) {

	for (int i = 0; i < detector->params.number_of_stations; i++) {
		// now remove all information which are useless
		for (int c_j = 0; c_j < detector->stations[i].strings.size(); c_j++) {
			for (int c_k = 0;
					c_k < detector->stations[i].strings[c_j].antennas.size();
					c_k++) {
				stations[i].strings[c_j].antennas[c_k].clear_useless(settings1); // clear data in antenna which stored in previous event
				//stations[i].strings[c_j].antennas[c_k].clear();  // clear data in antenna which stored in previous event
			}
		}
		//clear_useless(settings1);
	}
}

// this one is for single signal
void Report::Select_Wave_Convlv_Exchange(Settings *settings1, Trigger *trigger, Detector *detector, int signalbin, vector<double> &V, int *noise_ID, int ID, int StationIndex) {

	int BINSIZE = settings1->NFOUR / 2;

   // int BINSIZE = detector->stations[StationIndex].NFOUR/2;
	int bin_value;
	//vector <double> V_total_forconvlv;   // total time domain waveform (noise + signal)

	V_total_forconvlv.clear();

	// first, fill the noise values
	for (int bin = 0; bin < BINSIZE; bin++) {   //BINSIZE should be NFOUR/2
		bin_value = signalbin - BINSIZE / 2 + bin;

		// save the noise + signal waveform

		if (settings1->NOISE_CHANNEL_MODE == 0) {
			V_total_forconvlv.push_back(
					trigger->v_noise_timedomain[noise_ID[(int) (bin_value
							/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
							% settings1->DATA_BIN_SIZE)] + V[bin]);
		} else if (settings1->NOISE_CHANNEL_MODE == 1) {
			V_total_forconvlv.push_back(
					trigger->v_noise_timedomain_ch[GetChNumFromArbChID(detector,
							ID, StationIndex, settings1) - 1][noise_ID[(int) (bin_value
							/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
							% settings1->DATA_BIN_SIZE)] + V[bin]);
		}

		else if (settings1->NOISE_CHANNEL_MODE == 2) {
			if ((GetChNumFromArbChID(detector, ID, StationIndex, settings1) - 1)
					< 8) {
				V_total_forconvlv.push_back(
						trigger->v_noise_timedomain_ch[GetChNumFromArbChID(
								detector, ID, StationIndex, settings1) - 1][noise_ID[(int) (bin_value
								/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
								% settings1->DATA_BIN_SIZE)] + V[bin]);
			} else {
				V_total_forconvlv.push_back(
						trigger->v_noise_timedomain_ch[8][noise_ID[(int) (bin_value
								/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
								% settings1->DATA_BIN_SIZE)] + V[bin]);
			}
		}

	}

	// do myconvlv and replace the diode response array
 //   trigger->myconvlv( V_total_forconvlv, detector->stations[StationIndex].NFOUR, detector->stations[StationIndex].TIMESTEP, detector->fdiode_real, V_total_forconvlv);
 //   trigger->myconvlv( V_total_forconvlv, detector->stations[StationIndex].NFOUR, detector->fdiode_real, V_total_forconvlv);

	trigger->myconvlv(V_total_forconvlv, BINSIZE, detector->fdiode_real,
			V_total_forconvlv);

	// do replace the part we get from noise + signal
	for (int bin = signalbin - BINSIZE / 2 + (trigger->maxt_diode_bin);
			bin < signalbin + BINSIZE / 2; bin++) {
		trigger->Full_window[ID][bin] = V_total_forconvlv[bin - signalbin
				+ BINSIZE / 2];
		trigger->Full_window_V[ID][bin] += V[bin - signalbin + BINSIZE / 2];

		//
		// electronics saturation effect
		if (trigger->Full_window_V[ID][bin] > settings1->V_SATURATION)
			trigger->Full_window_V[ID][bin] = settings1->V_SATURATION;
		else if (trigger->Full_window_V[ID][bin]
				< -1. * settings1->V_SATURATION)
			trigger->Full_window_V[ID][bin] = -1. * settings1->V_SATURATION;

	}

	//V_total_forconvlv.clear();

}

// this one is for two connected signals 
void Report::Select_Wave_Convlv_Exchange(Settings *settings1, Trigger *trigger, Detector *detector, int signalbin1, int signalbin2, vector<double> &V1, vector<double> &V2, int *noise_ID, int ID, int StationIndex) {

	int BINSIZE = settings1->NFOUR / 2;

   // int BINSIZE = detector->stations[StationIndex].NFOUR/2;
	int bin_value;
	int signal_dbin = signalbin2 - signalbin1;
	//vector <double> V_total_forconvlv;   // total time domain waveform (noise + signal)
	double V_tmp[BINSIZE * 2];
	for (int bin_tmp = 0; bin_tmp < BINSIZE * 2; bin_tmp++) {
		V_tmp[bin_tmp] = 0.;
	}

	V_total_forconvlv.clear();

	// first, fill the noise values
	for (int bin = 0; bin < BINSIZE * 2; bin++) {
		bin_value = signalbin1 - BINSIZE / 2 + bin;

		// save the noise waveform
		if (settings1->NOISE_CHANNEL_MODE == 0) {
			V_total_forconvlv.push_back(
					trigger->v_noise_timedomain[noise_ID[(int) (bin_value
							/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
							% settings1->DATA_BIN_SIZE)]);
		} else if (settings1->NOISE_CHANNEL_MODE == 1) {
			V_total_forconvlv.push_back(
					trigger->v_noise_timedomain_ch[GetChNumFromArbChID(detector,
							ID, StationIndex, settings1) - 1][noise_ID[(int) (bin_value
							/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
							% settings1->DATA_BIN_SIZE)]);
		} else if (settings1->NOISE_CHANNEL_MODE == 2) {
			if ((GetChNumFromArbChID(detector, ID, StationIndex, settings1) - 1)
					< 8) {
				V_total_forconvlv.push_back(
						trigger->v_noise_timedomain_ch[GetChNumFromArbChID(
								detector, ID, StationIndex, settings1) - 1][noise_ID[(int) (bin_value
								/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
								% settings1->DATA_BIN_SIZE)]);
			} else {
				V_total_forconvlv.push_back(
						trigger->v_noise_timedomain_ch[8][noise_ID[(int) (bin_value
								/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
								% settings1->DATA_BIN_SIZE)]);
			}
		}

		// exchange from pure noise to noise + signal
		if (bin < signal_dbin) {  // bins where only first signal is shown
			V_total_forconvlv[bin] = V_total_forconvlv[bin] + V1[bin];
			V_tmp[bin] = V1[bin];
		} else if (bin < BINSIZE) { // bins where first + second signal is shown
			V_total_forconvlv[bin] = V_total_forconvlv[bin] + V1[bin]
					+ V2[bin - signal_dbin];
			V_tmp[bin] = V1[bin] + V2[bin - signal_dbin];
		} else if (bin < BINSIZE + signal_dbin) { // bins where only second signal is shown
			V_total_forconvlv[bin] = V_total_forconvlv[bin]
					+ V2[bin - signal_dbin];
			V_tmp[bin] = V2[bin - signal_dbin];
		}

	}

	// do myconvlv and replace the diode response array
 //   trigger->myconvlv( V_total_forconvlv, detector->stations[StationIndex].NFOUR, detector->stations[StationIndex].TIMESTEP, detector->fdiode_real_double, V_total_forconvlv);
 //   trigger->myconvlv( V_total_forconvlv, detector->stations[StationIndex].NFOUR, detector->fdiode_real_double, V_total_forconvlv);

	trigger->myconvlv(V_total_forconvlv, BINSIZE * 2,
			detector->fdiode_real_double, V_total_forconvlv);

	// do replace the part we get from noise + signal
	for (int bin = signalbin1 - BINSIZE / 2 + (trigger->maxt_diode_bin);
			bin < signalbin1 + BINSIZE / 2 + BINSIZE; bin++) {
		trigger->Full_window[ID][bin] = V_total_forconvlv[bin - signalbin1
				+ BINSIZE / 2];
		trigger->Full_window_V[ID][bin] +=
				V_tmp[bin - signalbin1 + BINSIZE / 2];
		//
		// electronics saturation effect
		if (trigger->Full_window_V[ID][bin] > settings1->V_SATURATION)
			trigger->Full_window_V[ID][bin] = settings1->V_SATURATION;
		else if (trigger->Full_window_V[ID][bin]
				< -1. * settings1->V_SATURATION)
			trigger->Full_window_V[ID][bin] = -1. * settings1->V_SATURATION;
	}

	//V_total_forconvlv.clear();

}

// this one is for three connected signals 
void Report::Select_Wave_Convlv_Exchange(Settings *settings1, Trigger *trigger, Detector *detector, int signalbin0, int signalbin1, int signalbin2, vector<double> &V0, vector<double> &V1, vector<double> &V2, int *noise_ID, int ID, int StationIndex) {

	int BINSIZE = settings1->NFOUR / 2;

   // int BINSIZE = detector->stations[StationIndex].NFOUR/2;
	int bin_value;
	int signal_dbin = signalbin2 - signalbin1;
	int signal_dbin0 = signalbin1 - signalbin0;
	//vector <double> V_total_forconvlv;   // total time domain waveform (noise + signal)
	double V_tmp[BINSIZE * 2];
	for (int bin_tmp = 0; bin_tmp < BINSIZE * 2; bin_tmp++) {
		V_tmp[bin_tmp] = 0.;
	}

	V_total_forconvlv.clear();

	// first, fill the noise values
	for (int bin = 0; bin < BINSIZE * 2; bin++) {
		bin_value = signalbin1 - BINSIZE / 2 + bin;

		// save the noise waveform
		if (settings1->NOISE_CHANNEL_MODE == 0) {
			V_total_forconvlv.push_back(
					trigger->v_noise_timedomain[noise_ID[(int) (bin_value
							/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
							% settings1->DATA_BIN_SIZE)]);
		} else if (settings1->NOISE_CHANNEL_MODE == 1) {
			V_total_forconvlv.push_back(
					trigger->v_noise_timedomain_ch[GetChNumFromArbChID(detector,
							ID, StationIndex, settings1) - 1][noise_ID[(int) (bin_value
							/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
							% settings1->DATA_BIN_SIZE)]);
		} else if (settings1->NOISE_CHANNEL_MODE == 2) {
			if ((GetChNumFromArbChID(detector, ID, StationIndex, settings1) - 1)
					< 8) {
				V_total_forconvlv.push_back(
						trigger->v_noise_timedomain_ch[GetChNumFromArbChID(
								detector, ID, StationIndex, settings1) - 1][noise_ID[(int) (bin_value
								/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
								% settings1->DATA_BIN_SIZE)]);
			} else {
				V_total_forconvlv.push_back(
						trigger->v_noise_timedomain_ch[8][noise_ID[(int) (bin_value
								/ settings1->DATA_BIN_SIZE)]][(int) (bin_value
								% settings1->DATA_BIN_SIZE)]);
			}
		}

		// exchange from pure noise to noise + signal
		if (bin < signal_dbin) {  // bins where no second signal is shown
			if (signal_dbin0 + bin < BINSIZE) { // previous signal is also here!
				V_total_forconvlv[bin] = V_total_forconvlv[bin] + V1[bin]
						+ V0[signal_dbin0 + bin];
				V_tmp[bin] = V1[bin] + V0[signal_dbin0 + bin];
			} else {  // no previous signal, and next signal
				V_total_forconvlv[bin] = V_total_forconvlv[bin] + V1[bin];
				V_tmp[bin] = V1[bin];
			}
		} else if (bin < BINSIZE) { // bins where first + second signal is shown
			if (signal_dbin0 + bin < BINSIZE) { // previous signal is also here!
				V_total_forconvlv[bin] = V_total_forconvlv[bin] + V1[bin]
						+ V0[signal_dbin0 + bin] + V2[bin - signal_dbin];
				V_tmp[bin] = V1[bin] + V0[signal_dbin0 + bin]
						+ V2[bin - signal_dbin];
			} else {  // no previous signal, and next signal
				V_total_forconvlv[bin] = V_total_forconvlv[bin] + V1[bin]
						+ V2[bin - signal_dbin];
				V_tmp[bin] = V1[bin] + V2[bin - signal_dbin];
			}
		} else if (bin < BINSIZE + signal_dbin) { // bins where only second signal is shown
			V_total_forconvlv[bin] = V_total_forconvlv[bin]
					+ V2[bin - signal_dbin];
			V_tmp[bin] = V2[bin - signal_dbin];
		}

	}

	// do myconvlv and replace the diode response array
 //   trigger->myconvlv( V_total_forconvlv, detector->stations[StationIndex].NFOUR, detector->stations[StationIndex].TIMESTEP, detector->fdiode_real_double, V_total_forconvlv);
 //   trigger->myconvlv( V_total_forconvlv, detector->stations[StationIndex].NFOUR, detector->fdiode_real_double, V_total_forconvlv);

	trigger->myconvlv(V_total_forconvlv, BINSIZE * 2,
			detector->fdiode_real_double, V_total_forconvlv);

	// do replace the part we get from noise + signal
	for (int bin = signalbin1 - BINSIZE / 2 + (trigger->maxt_diode_bin);
			bin < signalbin1 + BINSIZE / 2 + BINSIZE; bin++) {
		trigger->Full_window[ID][bin] = V_total_forconvlv[bin - signalbin1
				+ BINSIZE / 2];
		trigger->Full_window_V[ID][bin] +=
				V_tmp[bin - signalbin1 + BINSIZE / 2];
		//
		// electronics saturation effect
		if (trigger->Full_window_V[ID][bin] > settings1->V_SATURATION)
			trigger->Full_window_V[ID][bin] = settings1->V_SATURATION;
		else if (trigger->Full_window_V[ID][bin]
				< -1. * settings1->V_SATURATION)
			trigger->Full_window_V[ID][bin] = -1. * settings1->V_SATURATION;
	}

	//V_total_forconvlv.clear();

}

void Report::Apply_Gain_Offset(Settings *settings1, Trigger *trigger, Detector *detector, int ID, int StationIndex) {

	int string_num = detector->getStringfromArbAntID(StationIndex, ID);
	int ant_num = detector->getAntennafromArbAntID(StationIndex, ID);

	//int channel_num = detector->GetChannelfromStringAntenna ( StationIndex, string_num, ant_num );
	int channel_num = detector->GetChannelfromStringAntenna(StationIndex,
			string_num, ant_num, settings1);
	//cout<<"station "<<StationIndex<<" ch"<<channel_num<<" applying gain offset "<<detector->GetGainOffset( StationIndex, channel_num, settings1 )<<" applying"<<endl;

	for (int bin = 0; bin < settings1->DATA_BIN_SIZE; bin++) { // test for full window

		trigger->Full_window[ID][bin] = (trigger->Full_window[ID][bin]
				* detector->GetGainOffset(StationIndex, channel_num - 1,
						settings1)
				* detector->GetGainOffset(StationIndex, channel_num - 1,
						settings1)); // offset in voltage factor so we need power (V^2 factor to diode response)
		trigger->Full_window_V[ID][bin] = (trigger->Full_window_V[ID][bin]
				* detector->GetGainOffset(StationIndex, channel_num - 1,
						settings1)); // gain in voltage factor
	}

}

int Report::GetChNumFromArbChID(Detector *detector, int ID, int StationIndex, Settings *settings1) {

	int string_num = detector->getStringfromArbAntID(StationIndex, ID);
	int ant_num = detector->getAntennafromArbAntID(StationIndex, ID);

	int StationID = detector->stations[StationIndex].StationID;
	//    cout << "Station ID: " << StationID << endl;
	//    cout << "string_num: " << string_num << endl;
	//    cout << "ant_num: " << ant_num << endl;
	//    cout << StationID << " : " << string_num << " : " << ant_num << endl;

	//int channel_num = detector->GetChannelfromStringAntenna ( StationIndex, string_num, ant_num );
	int channel_num = detector->GetChannelfromStringAntenna(StationID,
			string_num, ant_num, settings1);

	return channel_num;
}

Vector Report::GetPolarization(Vector &nnu, Vector &launch_vector) {
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
	if (nnu * launch_vector > 0 && n_pol * nnu > 0)
		cout << "error in GetPolarization.  \n";

	return n_pol;
} //GetPolarization

void Report::GetParameters(Position &src, Position &trg, Vector &nnu, double &viewangle, double receive_angle, Vector &launch_vector, Vector &receive_vector, Vector &n_trg_slappy, Vector &n_trg_pokey) {

	viewangle = PI / 2. - viewangle;  // viewangle was actually launch angle

	launch_vector = (trg.Cross(trg.Cross(src))).Rotate(viewangle,
			trg.Cross(src));
	launch_vector = launch_vector.Unit();
	viewangle = launch_vector.Angle(nnu);

	//cout<<"launch_vector angle between R1 (trg) : "<<launch_vector.Angle(trg)<<"\n";

	receive_vector = trg.Rotate(receive_angle, src.Cross(trg));
	receive_vector = receive_vector.Unit();

	n_trg_pokey = trg.Unit();
	n_trg_slappy = (trg.Cross(src)).Unit();

}

double Report::GaintoHeight(double gain, double freq, double n_medium) {

	// from gain=4*pi*A_eff/lambda^2
	// and h_eff=2*sqrt(A_eff*Z_rx/Z_air)
	// gain is unitless value

	return 2
			* sqrt(
					gain / 4 / PI * CLIGHT * CLIGHT
							/ (freq * freq * n_medium * n_medium) * Zr
							/ (Z0 / n_medium)); // n_medium parts are changed from icemc(I believe this is correct one; E. Hong)
}

void Report::ApplyAntFactors(double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &pol_factor, double &vmmhz) { // vmmhz is input and output. output will have some antenna factors on it

	//double pol_factor;

	if (ant_type == 0) {    // if v pol
		pol_factor = n_trg_pokey * Pol_vector;
	} else if (ant_type == 1) {   // if h pol
		pol_factor = n_trg_slappy * Pol_vector;
	}
	pol_factor = abs(pol_factor);

	// apply 3dB spliter, d nu to prepare FFT
	// now actually vmmhz is not V/m/MHz but V/m/Hz unit
	//vmmhz = vmmhz/sqrt(2.)/(settings1->TIMESTEP*1.E6); //sqrt(2) for 3dB spliter for TURF, SURF
	vmmhz = vmmhz / sqrt(2.) / (1.E6); //sqrt(2) for 3dB spliter for TURF, SURF. 1/TIMESTEP moved to MakeArraysforFFT
	// 1/(1.E6) for V/MHz to V/Hz

	// apply antenna effective height and 0.5 (to calculate power with heff), and polarization factor
	// not vmmhz is actually V/Hz unit
	vmmhz = vmmhz * 0.5 * heff * pol_factor;

	// now we have to use MakeArraysforFFT to change signal arrays for FFT
}

void Report::ApplyAntFactors_Tdomain(double AntPhase, double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img, Settings *settings1) { // vm is input and output. output will have some antenna factors on it

	//double pol_factor;

	if (ant_type == 0) {    // if v pol
		pol_factor = n_trg_pokey * Pol_vector;
	} else if (ant_type == 1) {   // if h pol
		pol_factor = n_trg_slappy * Pol_vector;
	}
	pol_factor = abs(pol_factor);

	if (settings1->PHASE_SKIP_MODE != 1) {

		double phase_current;

		if (vm_real != 0.) {

			phase_current = atan(vm_img / vm_real);

			// phase in +-PI range
			if (vm_real < 0.) {
				if (vm_img > 0.)
					phase_current += PI;
				else if (vm_img < 0.)
					phase_current -= PI;
			}
		} else {

			if (vm_img > 0.)
				phase_current = PI;
			else if (vm_img < 0.)
				phase_current = -PI;
			else
				phase_current = 0.;
		}

		// V amplitude
		double v_amp = sqrt(vm_real * vm_real + vm_img * vm_img) / sqrt(2.)
				* 0.5 * heff * pol_factor; // sqrt(2) for 3dB splitter for TURF, SURF, 0.5 to calculate power with heff

		// real, img terms with phase shift
		vm_real = v_amp * cos(phase_current + AntPhase * RADDEG);
		vm_img = v_amp * sin(phase_current + AntPhase * RADDEG);

		//vm_real = v_amp * cos( phase_current - AntPhase*RADDEG ); // subtract AntPhase for four1 function's equation definition (inverse in img values)
		//vm_img = v_amp * sin( phase_current - AntPhase*RADDEG );
	}

	else { // only amplitude

		vm_real = vm_real / sqrt(2.) * 0.5 * heff * pol_factor; // only amplitude

		vm_img = vm_img / sqrt(2.) * 0.5 * heff * pol_factor; // only amplitude
	}

}

void Report::ApplyAntFactors_Tdomain_Transmitter(double AntPhase, double heff, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_real, double &vm_img,Settings *settings1) { // vm is input and output. output will have some antenna factors on it

	//double pol_factor;

	if (ant_type == 0) {    // if v pol
		pol_factor = n_trg_pokey * Pol_vector;
	} else if (ant_type == 1) {   // if h pol
		pol_factor = n_trg_slappy * Pol_vector;
	}
	pol_factor = abs(pol_factor);

	if (settings1->PHASE_SKIP_MODE != 1) {

		double phase_current;

		if (vm_real != 0.) {

			phase_current = atan(vm_img / vm_real);

			// phase in +-PI range
			if (vm_real < 0.) {
				if (vm_img > 0.)
					phase_current += PI;
				else if (vm_img < 0.)
					phase_current -= PI;
			}
		} else {

			if (vm_img > 0.)
				phase_current = PI;
			else if (vm_img < 0.)
				phase_current = -PI;
			else
				phase_current = 0.;
		}

		// V amplitude
		double v_amp = sqrt(vm_real * vm_real + vm_img * vm_img) / sqrt(2.)
				* 0.5 * heff * pol_factor; // sqrt(2) for 3dB splitter for TURF, SURF, 0.5 to calculate power with heff

		// real, img terms with phase shift
		//vm_real = v_amp * cos( phase_current + AntPhase*RADDEG );
		//vm_img = v_amp * sin( phase_current + AntPhase*RADDEG );

		vm_real = v_amp * cos(phase_current - AntPhase * RADDEG); // subtract AntPhase for four1 function's equation definition (inverse in img values)
		vm_img = v_amp * sin(phase_current - AntPhase * RADDEG);
	}

	else { // only amplitude

		vm_real = vm_real / sqrt(2.) * 0.5 * heff * pol_factor; // only amplitude

		vm_img = vm_img / sqrt(2.) * 0.5 * heff * pol_factor; // only amplitude
	}

}

void Report::ApplyAntFactors_Tdomain_FirstTwo(double heff, double heff_lastbin, Vector &n_trg_pokey, Vector &n_trg_slappy, Vector &Pol_vector, int ant_type, double &pol_factor, double &vm_bin0, double &vm_bin1) { // vm is input and output. output will have some antenna factors on it

	//double pol_factor;

	if (ant_type == 0) {    // if v pol
		pol_factor = n_trg_pokey * Pol_vector;
	} else if (ant_type == 1) {   // if h pol
		pol_factor = n_trg_slappy * Pol_vector;
	}
	pol_factor = abs(pol_factor);

	vm_bin0 = vm_bin0 / sqrt(2.) * 0.5 * heff * pol_factor; // sqrt(2) for 3dB splitter for TURF, SURF, 0.5 to calculate power with heff
	vm_bin1 = vm_bin1 / sqrt(2.) * 0.5 * heff_lastbin * pol_factor; // sqrt(2) for 3dB splitter for TURF, SURF, 0.5 to calculate power with heff

}

void Report::ApplyFilter(int bin_n, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetFilterGain(bin_n)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyPreamp(int bin_n, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetPreampGain(bin_n)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyFOAM(int bin_n, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetFOAMGain(bin_n)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyFilter_NFOUR(int bin_n, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetFilterGain_NFOUR(bin_n)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyFilter_OutZero(double freq, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetFilterGain_1D_OutZero(freq)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyElect_Tdomain(double freq, Detector *detector, double &vm_real, double &vm_img, Settings *settings1) { // read elect chain gain (unitless), phase (rad) and apply to V/m

	if (settings1->PHASE_SKIP_MODE == 0) {

		double phase_current;

		if (vm_real != 0.) {

			phase_current = atan(vm_img / vm_real);

			// phase in +-PI range
			if (vm_real < 0.) {
				if (vm_img > 0.)
					phase_current += PI;
				else if (vm_img < 0.)
					phase_current -= PI;
			}
		} else {

			if (vm_img > 0.)
				phase_current = PI;
			else if (vm_img < 0.)
				phase_current = -PI;
			else
				phase_current = 0.;
		}

		// V amplitude
		double v_amp = sqrt(vm_real * vm_real + vm_img * vm_img)
				* detector->GetElectGain_1D_OutZero(freq); // apply gain (unitless) to amplitude

				// real, img terms with phase shift
				//vm_real = v_amp * cos( phase_current + detector->GetElectPhase_1D(freq) );
				//vm_img = v_amp * sin( phase_current + detector->GetElectPhase_1D(freq) );

		vm_real = v_amp * cos(phase_current - detector->GetElectPhase_1D(freq));
		vm_img = v_amp * sin(phase_current - detector->GetElectPhase_1D(freq));
	}

	else {

		vm_real = vm_real * detector->GetElectGain_1D_OutZero(freq); // only amplitude

		vm_img = vm_img * detector->GetElectGain_1D_OutZero(freq); // only amplitude
	}

}

void Report::ApplyElect_Tdomain_FirstTwo(double freq0, double freq1, Detector *detector, double &vm_bin0, double &vm_bin1) { // read elect chain gain (unitless), phase (rad) and apply to V/m

	vm_bin0 = vm_bin0 * detector->GetElectGain_1D_OutZero(freq0);
	vm_bin1 = vm_bin1 * detector->GetElectGain_1D_OutZero(freq1);

}

void Report::ApplyPreamp_NFOUR(int bin_n, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetPreampGain_NFOUR(bin_n)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyPreamp_OutZero(double freq, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetPreampGain_1D_OutZero(freq)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyFOAM_NFOUR(int bin_n, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetFOAMGain_NFOUR(bin_n)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyFOAM_OutZero(double freq, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetFOAMGain_1D_OutZero(freq)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyRFCM(int ch, int bin_n, Detector *detector, double &vmmhz, double RFCM_OFFSET) { // read RFCM gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz
			* pow(10., (detector->GetRFCMGain(ch, bin_n) + RFCM_OFFSET) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyFilter_databin(int bin_n, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetFilterGain_databin(bin_n)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyPreamp_databin(int bin_n, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetPreampGain_databin(bin_n)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyFOAM_databin(int bin_n, Detector *detector, double &vmmhz) { // read filter gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz * pow(10., (detector->GetFOAMGain_databin(bin_n)) / 20.); // from dB to unitless gain for voltage

}

void Report::ApplyRFCM_databin(int ch, int bin_n, Detector *detector, double &vmmhz, double RFCM_OFFSET) { // read RFCM gain in dB and apply unitless gain to vmmhz

	vmmhz = vmmhz
			* pow(10.,
					(detector->GetRFCMGain_databin(ch, bin_n) + RFCM_OFFSET)
							/ 20.);   // from dB to unitless gain for voltage

}

void Report::ApplyNoiseFig_databin(int ch, int bin_n, Detector *detector, double &vmmhz, Settings *settings1) { // read noise figure and apply unitless gain to vmmhz
	//cout<<"Entered ApplyNoiseFig_databin    ";
	//cout<<"vmmhz: "<<vmmhz;

	double tempNoise = vmmhz * vmmhz;
	//  cout<<"    1: "<<detector->GetNoiseFig_databin(ch,bin_n);
	//  cout<<"    2: "<<detector->GetTransm_databin(ch, bin_n);
	//  cout<<"    3: "<<settings1->NOISE_TEMP;
	//FILE *fp = fopen("noiseTemp_NOISE_CHANNEL_MODE_0.txt","a+");
	if (detector->GetNoiseFig_databin(ch, bin_n) > 1.0) {
		vmmhz = TMath::Sqrt(
				tempNoise * (detector->GetTransm_databin(ch, bin_n))
						+ tempNoise / (settings1->NOISE_TEMP) * 220.0
								* (detector->GetNoiseFig_databin(ch, bin_n)
										- 1.0));
		//cout<<"   vmmhz   "<<vmmhz;
		//cout<<"   bin: "<<bin_n<<" freq: "<<detector->GetFreq(bin_n) << " noise temp: " << 246.0*((detector->GetTransm_databin(ch, bin_n)) + 1.0/(settings1->NOISE_TEMP)*220.0*(detector->GetNoiseFig_databin(ch,bin_n) - 1.0) ) << endl;  
		//fprintf(fp, "%f   ",246.0*((detector->GetTransm_databin(ch, bin_n)) + 1.0/(settings1->NOISE_TEMP)*220.0*(detector->GetNoiseFig_databin(ch,bin_n) - 1.0) ));
	} else {
		vmmhz = vmmhz;
		//fprintf(fp, "%f   ",246.0);
	}

	//fclose(fp);
	//    if(ch==0 && bin_n==2050) cout << "Noise ch: "<< ch <<"   "<< detector->GetNoiseFig_databin(ch,bin_n) << "   " << detector->GetTransm_databin(ch, bin_n)
	//    << "   " << TMath::Sqrt( tempNoise ) << "   " << vmmhz << endl;
	//  if(ch==8 && bin_n==2050) cout << "Noise ch: "<< ch <<"   "<< detector->GetNoiseFig_databin(ch,bin_n) << "   " << detector->GetTransm_databin(ch, bin_n)
	//    << "   " << TMath::Sqrt( tempNoise ) << "   " << vmmhz << endl;
	//    if(ch==24 && bin_n==2000) cout << "Noise ch: "<< ch << "   " << detector->GetNoiseFig_databin(ch,bin_n)*220.0/246.0 << "   " << TMath::Sqrt(detector->GetTransm_databin(ch, bin_n) ) << endl;

}

void Report::GetAngleAnt(Vector &rec_vector, Position &antenna, double &ant_theta, double &ant_phi) { //ant_theta and ant_phi is in degree 

	// need to fix some parts.
	// currently phi is not correct.
	// I think we have to actually rotate (0,-1,0) to antenna location and get phi from rotated (0,-1,0)
	ant_theta = rec_vector.Angle(antenna) * DEGRAD;

	Vector unit0_10;
	Vector unit100;
	unit0_10.SetXYZ(0, -1, 0);
	unit100.SetXYZ(1, 0, 0);
	Vector proj;

	proj = antenna.Cross(rec_vector.Cross(antenna));

	if (proj * unit100 > 0.) { // rec_vector is pointing to positive x direction
		ant_phi = 360. - (proj.Angle(unit0_10)) * DEGRAD;
	} else {
		ant_phi = proj.Angle(unit0_10) * DEGRAD;
	}

}

// generate DATA_BIN_SIZE sized noise waveform array in time domain
void Report::GetNoiseWaveforms(Settings *settings1, Detector *detector, double v_noise, double *vnoise) {

	if (settings1->NOISE == 0) { // NOISE == 0 : flat thermal noise with Johnson-Nyquist noise
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
		for (int k = 0; k < settings1->DATA_BIN_SIZE / 2; k++) {

			V_tmp = v_noise / sqrt(2.); // copy original flat H_n [V] value, and apply 1/sqrt2 for SURF/TURF divide same as signal

			if (settings1->USE_TESTBED_RFCM_ON == 0) {
				ApplyFilter_databin(k, detector, V_tmp);
				ApplyPreamp_databin(k, detector, V_tmp);
				ApplyFOAM_databin(k, detector, V_tmp);
				if (settings1->APPLY_NOISE_FIGURE == 1) {
					ApplyNoiseFig_databin(0, k, detector, V_tmp, settings1);
				}
			} else if (settings1->USE_TESTBED_RFCM_ON == 1) {
				// apply RFCM gain
				// as this mode don't have different chs' noise waveform separately, just use ch0 RFCM gain...
				cerr
						<< "Trying not allowed mode : NOISE_CHANNEL_MODE=0 and USE_TESTBED_RFCM_ON=1 same time!"
						<< endl;
				break;
				//ApplyRFCM_databin(0, detector, V_tmp);
			}

			Vfft_noise_before.push_back(V_tmp);

			current_phase = noise_phase[k];

			//Tools::get_random_rician( 0., 0., sqrt(2./ M_PI) * V_tmp, current_amplitude, current_phase);    // use real value array value
			Tools::get_random_rician(0., 0., sqrt(2. / M_PI) / 1.177 * V_tmp,
					current_amplitude, current_phase); // use real value array value, extra 1/1.177 to make total power same with "before random_rician".

			// vnoise is currently noise spectrum (before fft, unit : V)
			//vnoise[2 * k] = sqrt(current_amplitude) * cos(noise_phase[k]);
			//vnoise[2 * k + 1] = sqrt(current_amplitude) * sin(noise_phase[k]);
			vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
			vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);

			//vnoise[2 * k] = (V_tmp) * cos(noise_phase[k]);
			//vnoise[2 * k + 1] = (V_tmp) * sin(noise_phase[k]);

			Vfft_noise_after.push_back(vnoise[2 * k]);
			Vfft_noise_after.push_back(vnoise[2 * k + 1]);

			// inverse FFT normalization factor!
			vnoise[2 * k] *= 2. / ((double) settings1->DATA_BIN_SIZE);
			vnoise[2 * k + 1] *= 2. / ((double) settings1->DATA_BIN_SIZE);

		}

		// now vnoise is time domain waveform
		Tools::realft(vnoise, -1, settings1->DATA_BIN_SIZE);

		// save timedomain noise to Report class
		/*
		 for (int k=0; k<settings1->DATA_BIN_SIZE; k++) {
		 V_noise_timedomain.push_back( vnoise[k] );
		 }
		 */

	} else {  // currently there are no more options!
		cout << "no noise option for NOISE = " << settings1->NOISE << endl;
	}

}

// generate DATA_BIN_SIZE sized noise waveform array in time domain
void Report::GetNoiseWaveforms_ch(Settings *settings1, Detector *detector, double v_noise, double *vnoise, int ch) {
	//  cout << "channel: " << ch << endl;
	if (settings1->NOISE == 0) { // NOISE == 0 : flat thermal noise with Johnson-Nyquist noise

		Vfft_noise_after.clear();  // remove previous Vfft_noise values
		Vfft_noise_before.clear();  // remove previous Vfft_noise values
		//V_noise_timedomain.clear(); // remove previous V_noise_timedomain values

		double V_tmp; // copy original flat H_n [V] value
		double current_amplitude, current_phase;

		GetNoisePhase(settings1); // get random phase for noise

		for (int k = 0; k < settings1->DATA_BIN_SIZE / 2; k++) {

			V_tmp = v_noise / sqrt(2.); // copy original flat H_n [V] value, and apply 1/sqrt2 for SURF/TURF divide same as signal

			if (settings1->USE_TESTBED_RFCM_ON == 0) {
				ApplyFilter_databin(k, detector, V_tmp);
				ApplyPreamp_databin(k, detector, V_tmp);
				ApplyFOAM_databin(k, detector, V_tmp);
				if (settings1->APPLY_NOISE_FIGURE == 1) {
					ApplyNoiseFig_databin(ch % 16, k, detector, V_tmp,
							settings1);
				}
			} else if (settings1->USE_TESTBED_RFCM_ON == 1) {
				// apply RFCM gain
				ApplyRFCM_databin(ch, k, detector, V_tmp,
						settings1->RFCM_OFFSET);
			}

			Vfft_noise_before.push_back(V_tmp);

			current_phase = noise_phase[k];

			//Tools::get_random_rician( 0., 0., sqrt(2./ M_PI) * V_tmp, current_amplitude, current_phase);    // use real value array value
			Tools::get_random_rician(0., 0., sqrt(2. / M_PI) / 1.177 * V_tmp,
					current_amplitude, current_phase); // use real value array value, extra 1/1.177 to make total power same with "before random_rician".

			// vnoise is currently noise spectrum (before fft, unit : V)
			//vnoise[2 * k] = sqrt(current_amplitude) * cos(noise_phase[k]);
			//vnoise[2 * k + 1] = sqrt(current_amplitude) * sin(noise_phase[k]);
			vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
			vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);

			Vfft_noise_after.push_back(vnoise[2 * k]);
			Vfft_noise_after.push_back(vnoise[2 * k + 1]);

			// inverse FFT normalization factor!
			vnoise[2 * k] *= 2. / ((double) settings1->DATA_BIN_SIZE);
			vnoise[2 * k + 1] *= 2. / ((double) settings1->DATA_BIN_SIZE);

		}

		//  cout << "After applying noise figure" << endl;

		// now vnoise is time domain waveform
		Tools::realft(vnoise, -1, settings1->DATA_BIN_SIZE);

		// save timedomain noise to Report class
		/*
		 for (int k=0; k<settings1->DATA_BIN_SIZE; k++) {
		 V_noise_timedomain.push_back( vnoise[k] );
		 }
		 */
		//  cout << "After fft" << endl;
	}

	else if (settings1->NOISE == 1) { // NOISE == 1 : use Rayleigh dist. fit from TestBed data

		Vfft_noise_after.clear();  // remove previous Vfft_noise values
		Vfft_noise_before.clear();  // remove previous Vfft_noise values
		//V_noise_timedomain.clear(); // remove previous V_noise_timedomain values

		double V_tmp; // copy original flat H_n [V] value
		double current_amplitude, current_phase;

		GetNoisePhase(settings1); // get random phase for noise

		for (int k = 0; k < settings1->DATA_BIN_SIZE / 2; k++) {

			if (ch < detector->RayleighFit_ch) {

				Vfft_noise_before.push_back(
						detector->GetRayleighFit_databin(ch, k));

				current_phase = noise_phase[k];

				V_tmp = detector->GetRayleighFit_databin(ch, k)
						* sqrt(
								(double) settings1->DATA_BIN_SIZE
										/ (double) (settings1->NFOUR / 2.));

				//Tools::get_random_rician( 0., 0., sqrt(2./ M_PI) * V_tmp, current_amplitude, current_phase);    // use real value array value
				//Tools::get_random_rician( 0., 0., sqrt(2./M_PI)/1.177 * V_tmp, current_amplitude, current_phase);    // use real value array value, extra 1/1.177 to make total power same with "before random_rician".
				Tools::get_random_rician(0., 0., V_tmp, current_amplitude,
						current_phase); // use real value array value, extra 1/1.177 to make total power same with "before random_rician".

			}

			else {

				V_tmp = v_noise / sqrt(2.); // copy original flat H_n [V] value, and apply 1/sqrt2 for SURF/TURF divide same as signal

				if (settings1->USE_TESTBED_RFCM_ON == 0) {
					ApplyFilter_databin(k, detector, V_tmp);
					ApplyPreamp_databin(k, detector, V_tmp);
					ApplyFOAM_databin(k, detector, V_tmp);

				} else if (settings1->USE_TESTBED_RFCM_ON == 1) {
					// apply RFCM gain
					ApplyRFCM_databin(ch, k, detector, V_tmp,
							settings1->RFCM_OFFSET);
				}

				Vfft_noise_before.push_back(V_tmp);

				current_phase = noise_phase[k];

				//Tools::get_random_rician( 0., 0., sqrt(2./ M_PI) * V_tmp, current_amplitude, current_phase);    // use real value array value
				Tools::get_random_rician(0., 0.,
						sqrt(2. / M_PI) / 1.177 * V_tmp, current_amplitude,
						current_phase); // use real value array value, extra 1/1.177 to make total power same with "before random_rician".

			}

			// vnoise is currently noise spectrum (before fft, unit : V)
			//vnoise[2 * k] = sqrt(current_amplitude) * cos(noise_phase[k]);
			//vnoise[2 * k + 1] = sqrt(current_amplitude) * sin(noise_phase[k]);
			vnoise[2 * k] = (current_amplitude) * cos(noise_phase[k]);
			vnoise[2 * k + 1] = (current_amplitude) * sin(noise_phase[k]);

			Vfft_noise_after.push_back(vnoise[2 * k]);
			Vfft_noise_after.push_back(vnoise[2 * k + 1]);

			// inverse FFT normalization factor!
			vnoise[2 * k] *= 2. / ((double) settings1->DATA_BIN_SIZE);
			vnoise[2 * k + 1] *= 2. / ((double) settings1->DATA_BIN_SIZE);

		}

		// now vnoise is time domain waveform
		Tools::realft(vnoise, -1, settings1->DATA_BIN_SIZE);

		// save timedomain noise to Report class
		/*
		 for (int k=0; k<settings1->DATA_BIN_SIZE; k++) {
		 V_noise_timedomain.push_back( vnoise[k] );
		 }
		 */
	}

	else {  // currently there are no more options!
		cout << "no noise option for NOISE = " << settings1->NOISE << endl;
	}

}

void Report::GetNoisePhase(Settings *settings1) {

	noise_phase.clear();    // remove all previous noise_phase vector values

	//for (int i=0; i<settings1->NFOUR/4; i++) {
	for (int i = 0; i < settings1->DATA_BIN_SIZE / 2; i++) { // noise with DATA_BIN_SIZE bin array
		noise_phase.push_back(2 * PI * gRandom->Rndm()); // get random phase for flat thermal noise
	}
}

void Report::MakeArraysforFFT(Settings *settings1, Detector *detector, int StationIndex, vector<double> &vsignal_array, double *vsignal_forfft) {

	// from icemc, anita class MakeArraysforFFT
	int NFOUR = settings1->NFOUR;
	double TIMESTEP = settings1->TIMESTEP;

   // int NFOUR = detector->stations[StationIndex].NFOUR;
   // int TIMESTEP = detector->stations[StationIndex].TIMESTEP;
	Tools::Zero(vsignal_forfft, NFOUR / 2);

	if (settings1->AVZ_NORM_FACTOR_MODE == 0) { // use previous normalization factors

		double previous_value_e_even = 0.;
		double previous_value_e_odd = 0.;
		int count_nonzero = 0;
		int iprevious = 0;
		int ifirstnonzero = -1;
		int ilastnonzero = 2000;
		//for (int i=0;i<NFREQ;i++) {
		for (int i = 0; i < detector->GetFreqBin(); i++) {

			// freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
			// but there are only NFOUR/4 different values
			// it's the index among the NFOUR/4 that we're interested in
			int ifour = Tools::Getifreq(detector->GetFreq(i),
					detector->freq_forfft[0],
					detector->freq_forfft[NFOUR / 2 - 1], NFOUR / 4);

			if (ifour != -1 && 2 * ifour + 1 < NFOUR / 2) {
				count_nonzero++;
				if (ifirstnonzero == -1)
					ifirstnonzero = ifour;

				vsignal_forfft[2 * ifour] = vsignal_array[i] * 2
						/ ((double) NFOUR / 2) / (TIMESTEP); // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft

				//      cout << "ifour, vsignal is " << ifour << " " << vsignal_e_forfft[2*ifour] << "\n";

				vsignal_forfft[2 * ifour + 1] = vsignal_array[i] * 2
						/ ((double) NFOUR / 2) / (TIMESTEP); // phase is 90 deg.
				// the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting

				// how about we interpolate instead of doing a box average

				for (int j = iprevious + 1; j < ifour; j++) {
					vsignal_forfft[2 * j] =
							previous_value_e_even
									+ (vsignal_forfft[2 * ifour]
											- previous_value_e_even)
											* (double) (j - iprevious)
											/ (double) (ifour - iprevious);
					//  cout << "j, vsignal is " << j << " " << vsignal_e_forfft[2*j] << "\n";

					vsignal_forfft[2 * j + 1] = previous_value_e_odd
							+ (vsignal_forfft[2 * ifour + 1]
									- previous_value_e_odd)
									* (double) (j - iprevious)
									/ (double) (ifour - iprevious);
				}

				ilastnonzero = ifour;
				iprevious = ifour;
				previous_value_e_even = vsignal_forfft[2 * ifour];
				previous_value_e_odd = vsignal_forfft[2 * ifour + 1];
			}

		} // end loop over nfreq

		// don't applying any extra factor for the change in array (change of bin size)
		// as change in the bin size doesn't matter for the total energy
		// total energy is just same as integral over frequency range and this frequency range will not change
		//
		//cout << "ifirstnonzero, ilastnonzero are " << ifirstnonzero << " " << ilastnonzero << "\n";
		//cout << "non zero count is " << count_nonzero << "\n";
		//cout << "ratio is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
		for (int j = 0; j < NFOUR / 4; j++) {
			vsignal_forfft[2 * j] *= 1.
					/ sqrt(
							(double) count_nonzero
									/ (double) (ilastnonzero - ifirstnonzero));
			vsignal_forfft[2 * j + 1] *= 1.
					/ sqrt(
							(double) count_nonzero
									/ (double) (ilastnonzero - ifirstnonzero));
		}

		//  Tools::InterpolateComplex(vsignal_e_forfft,NFOUR/4);
		//Tools::InterpolateComplex(vsignal_h_forfft,NFOUR/4);
		for (int ifour = 0; ifour < NFOUR / 4; ifour++) {

			vsignal_forfft[2 * ifour] *= cos(settings1->PHASE * PI / 180.);
			vsignal_forfft[2 * ifour + 1] *= sin(settings1->PHASE * PI / 180.);

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
	} else if (settings1->AVZ_NORM_FACTOR_MODE == 1) { // use new (fixed) normalization factors

		double dF = 1. / ((double) (NFOUR / 2) * TIMESTEP); // in Hz
		//cout<<"dF1 : "<<dF<<endl;

		double dF_org = detector->GetFreq(1) - detector->GetFreq(0); // in Hz

		//double dF_factor = sqrt( dF_org / dF );
		double dF_factor = 1.;
		//double dF_factor = sqrt( dF / dF_org );

		int Norg = detector->GetFreqBin();

		double FreqOrg[Norg + 1];
		double VmMHzOrg[Norg + 1];

		double FreqNFOUR[NFOUR / 4 + 1]; // one more bin
		double VmMHzNFOUR[NFOUR / 4 + 1]; // one more bin

		for (int i = 0; i < Norg + 1; i++) {

			if (i == 0) {

				VmMHzOrg[i] = 0.;
				FreqOrg[i] = 0.;
			} else {

				VmMHzOrg[i] = vsignal_array[i - 1];
				FreqOrg[i] = detector->GetFreq(i - 1);
			}

		}

		for (int ifour = 0; ifour < NFOUR / 4 + 1; ifour++) {

			FreqNFOUR[ifour] = dF * (double) ifour;
			VmMHzNFOUR[ifour] = 0.;

		}

		//Tools::SimpleLinearInterpolation_OutZero( Norg, FreqOrg, vsignal_array, NFOUR/4+1, FreqNFOUR, VmMHzNFOUR );
		Tools::SimpleLinearInterpolation_OutZero(Norg + 1, FreqOrg, VmMHzOrg,
				NFOUR / 4 + 1, FreqNFOUR, VmMHzNFOUR);

		for (int ifour = 1; ifour < NFOUR / 4; ifour++) {

			// same amplitude, 2/(NFOUR/2) for inverse FFT normalization factor
			vsignal_forfft[2 * ifour] = VmMHzNFOUR[ifour] * 2
					/ ((double) NFOUR / 2) * dF_factor / TIMESTEP; // change to V/Hz, apply norm factor, change to Hn
			vsignal_forfft[2 * ifour + 1] = VmMHzNFOUR[ifour] * 2
					/ ((double) NFOUR / 2) * dF_factor / TIMESTEP;

			// apply phase
			vsignal_forfft[2 * ifour] *= cos(settings1->PHASE * PI / 180.);
			vsignal_forfft[2 * ifour + 1] *= sin(settings1->PHASE * PI / 180.);
		}

		// first and last freq bin read values to 0, 1 bin
		vsignal_forfft[0] = VmMHzNFOUR[0] * 2 / ((double) NFOUR / 2) * dF_factor
				/ TIMESTEP;
		vsignal_forfft[1] = VmMHzNFOUR[NFOUR / 4] * 2 / ((double) NFOUR / 2)
				* dF_factor / TIMESTEP;
	}

}

void Report::MakeArraysforFFT(Settings *settings1, Detector *detector, int StationIndex, double *vsignal_array, double *vsignal_forfft) {

	// from icemc, anita class MakeArraysforFFT
	int NFOUR = settings1->NFOUR;
	double TIMESTEP = settings1->TIMESTEP;

   // int NFOUR = detector->stations[StationIndex].NFOUR;
   // int TIMESTEP = detector->stations[StationIndex].TIMESTEP;
	Tools::Zero(vsignal_forfft, NFOUR / 2);

	if (settings1->AVZ_NORM_FACTOR_MODE == 0) { // use previous normalization factors

		double previous_value_e_even = 0.;
		double previous_value_e_odd = 0.;
		int count_nonzero = 0;
		int iprevious = 0;
		int ifirstnonzero = -1;
		int ilastnonzero = 2000;
		//for (int i=0;i<NFREQ;i++) {
		for (int i = 0; i < detector->GetFreqBin(); i++) {

			// freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
			// but there are only NFOUR/4 different values
			// it's the index among the NFOUR/4 that we're interested in
			int ifour = Tools::Getifreq(detector->GetFreq(i),
					detector->freq_forfft[0],
					detector->freq_forfft[NFOUR / 2 - 1], NFOUR / 4);

			if (ifour != -1 && 2 * ifour + 1 < NFOUR / 2) {
				count_nonzero++;
				if (ifirstnonzero == -1)
					ifirstnonzero = ifour;

				vsignal_forfft[2 * ifour] = vsignal_array[i] * 2
						/ ((double) NFOUR / 2) / (TIMESTEP); // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft

				//      cout << "ifour, vsignal is " << ifour << " " << vsignal_e_forfft[2*ifour] << "\n";

				vsignal_forfft[2 * ifour + 1] = vsignal_array[i] * 2
						/ ((double) NFOUR / 2) / (TIMESTEP); // phase is 90 deg.
				// the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting

				// how about we interpolate instead of doing a box average

				for (int j = iprevious + 1; j < ifour; j++) {
					vsignal_forfft[2 * j] =
							previous_value_e_even
									+ (vsignal_forfft[2 * ifour]
											- previous_value_e_even)
											* (double) (j - iprevious)
											/ (double) (ifour - iprevious);
					//  cout << "j, vsignal is " << j << " " << vsignal_e_forfft[2*j] << "\n";

					vsignal_forfft[2 * j + 1] = previous_value_e_odd
							+ (vsignal_forfft[2 * ifour + 1]
									- previous_value_e_odd)
									* (double) (j - iprevious)
									/ (double) (ifour - iprevious);
				}

				ilastnonzero = ifour;
				iprevious = ifour;
				previous_value_e_even = vsignal_forfft[2 * ifour];
				previous_value_e_odd = vsignal_forfft[2 * ifour + 1];
			}

		} // end loop over nfreq

		// don't applying any extra factor for the change in array (change of bin size)
		// as change in the bin size doesn't matter for the total energy
		// total energy is just same as integral over frequency range and this frequency range will not change
		//
		//cout << "ifirstnonzero, ilastnonzero are " << ifirstnonzero << " " << ilastnonzero << "\n";
		//cout << "non zero count is " << count_nonzero << "\n";
		//cout << "ratio is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
		for (int j = 0; j < NFOUR / 4; j++) {
			vsignal_forfft[2 * j] *= 1.
					/ sqrt(
							(double) count_nonzero
									/ (double) (ilastnonzero - ifirstnonzero));
			vsignal_forfft[2 * j + 1] *= 1.
					/ sqrt(
							(double) count_nonzero
									/ (double) (ilastnonzero - ifirstnonzero));
		}

		//  Tools::InterpolateComplex(vsignal_e_forfft,NFOUR/4);
		//Tools::InterpolateComplex(vsignal_h_forfft,NFOUR/4);
		for (int ifour = 0; ifour < NFOUR / 4; ifour++) {

			vsignal_forfft[2 * ifour] *= cos(settings1->PHASE * PI / 180.);
			vsignal_forfft[2 * ifour + 1] *= sin(settings1->PHASE * PI / 180.);

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
	} else if (settings1->AVZ_NORM_FACTOR_MODE == 1) { // use new (fixed) normalization factors

		double dF = 1. / ((double) (NFOUR / 2) * TIMESTEP); // in Hz
		//cout<<"dF1 : "<<dF<<endl;

		double dF_org = detector->GetFreq(1) - detector->GetFreq(0); // in Hz

		//double dF_factor = sqrt( dF_org / dF );
		double dF_factor = 1.;
		//double dF_factor = sqrt( dF / dF_org );

		int Norg = detector->GetFreqBin();

		double FreqOrg[Norg + 1];
		double VmMHzOrg[Norg + 1];

		double FreqNFOUR[NFOUR / 4 + 1]; // one more bin
		double VmMHzNFOUR[NFOUR / 4 + 1]; // one more bin

		for (int i = 0; i < Norg + 1; i++) {

			if (i == 0) {

				VmMHzOrg[i] = 0.;
				FreqOrg[i] = 0.;
			} else {

				VmMHzOrg[i] = vsignal_array[i - 1];
				FreqOrg[i] = detector->GetFreq(i - 1);
			}

		}

		for (int ifour = 0; ifour < NFOUR / 4 + 1; ifour++) {

			FreqNFOUR[ifour] = dF * (double) ifour;
			VmMHzNFOUR[ifour] = 0.;

		}

		//Tools::SimpleLinearInterpolation_OutZero( Norg, FreqOrg, vsignal_array, NFOUR/4+1, FreqNFOUR, VmMHzNFOUR );
		Tools::SimpleLinearInterpolation_OutZero(Norg + 1, FreqOrg, VmMHzOrg,
				NFOUR / 4 + 1, FreqNFOUR, VmMHzNFOUR);

		for (int ifour = 1; ifour < NFOUR / 4; ifour++) {

			// same amplitude, 2/(NFOUR/2) for inverse FFT normalization factor
			vsignal_forfft[2 * ifour] = VmMHzNFOUR[ifour] * 2
					/ ((double) NFOUR / 2) * dF_factor / TIMESTEP; // change to V/Hz, apply norm factor, change to Hn
			vsignal_forfft[2 * ifour + 1] = VmMHzNFOUR[ifour] * 2
					/ ((double) NFOUR / 2) * dF_factor / TIMESTEP;

			// apply phase
			vsignal_forfft[2 * ifour] *= cos(settings1->PHASE * PI / 180.);
			vsignal_forfft[2 * ifour + 1] *= sin(settings1->PHASE * PI / 180.);
		}

		// first and last freq bin read values to 0, 1 bin
		vsignal_forfft[0] = VmMHzNFOUR[0] * 2 / ((double) NFOUR / 2) * dF_factor
				/ TIMESTEP;
		vsignal_forfft[1] = VmMHzNFOUR[NFOUR / 4] * 2 / ((double) NFOUR / 2)
				* dF_factor / TIMESTEP;
	}

}

void Report::MakeArraysforFFT_noise(Settings *settings1, Detector *detector, int StationIndex, vector<double> &vsignal_array, double *vsignal_forfft) {

	// from icemc, anita class MakeArraysforFFT
	int NFOUR = settings1->NFOUR;
	double TIMESTEP = settings1->TIMESTEP;

   // int NFOUR = detector->stations[StationIndex].NFOUR;
   // int TIMESTEP = detector->stations[StationIndex].TIMESTEP;

	Tools::Zero(vsignal_forfft, NFOUR / 2);

	double previous_value_e_even = 0.;
	double previous_value_e_odd = 0.;
	int count_nonzero = 0;
	int iprevious = 0;
	int ifirstnonzero = -1;
	int ilastnonzero = 2000;
	//for (int i=0;i<NFREQ;i++) {
	for (int i = 0; i < detector->GetFreqBin(); i++) {

		// freq_forfft has NFOUR/2 elements because it is meant to cover real and imaginary values
		// but there are only NFOUR/4 different values
		// it's the index among the NFOUR/4 that we're interested in
		int ifour = Tools::Getifreq(detector->GetFreq(i),
				detector->freq_forfft[0], detector->freq_forfft[NFOUR / 2 - 1],
				NFOUR / 4);

		if (ifour != -1 && 2 * ifour + 1 < NFOUR / 2) {
			count_nonzero++;
			if (ifirstnonzero == -1)
				ifirstnonzero = ifour;

			vsignal_forfft[2 * ifour] = vsignal_array[i] * 2
					/ ((double) NFOUR / 2) / (TIMESTEP); // inverse fft normalization factor (2/(N/2)), 1/dt for change integration fft form to discrete numerical fft

			//      cout << "ifour, vsignal is " << ifour << " " << vsignal_e_forfft[2*ifour] << "\n";

			vsignal_forfft[2 * ifour + 1] = vsignal_array[i] * 2
					/ ((double) NFOUR / 2) / (TIMESTEP); // phase is 90 deg.
			// the 2/(nfour/2) needs to be included since were using Tools::realft with the -1 setting

			// how about we interpolate instead of doing a box average

			for (int j = iprevious + 1; j < ifour; j++) {
				vsignal_forfft[2 * j] = previous_value_e_even
						+ (vsignal_forfft[2 * ifour] - previous_value_e_even)
								* (double) (j - iprevious)
								/ (double) (ifour - iprevious);
				//  cout << "j, vsignal is " << j << " " << vsignal_e_forfft[2*j] << "\n";

				vsignal_forfft[2 * j + 1] = previous_value_e_odd
						+ (vsignal_forfft[2 * ifour + 1] - previous_value_e_odd)
								* (double) (j - iprevious)
								/ (double) (ifour - iprevious);
			}

			ilastnonzero = ifour;
			iprevious = ifour;
			previous_value_e_even = vsignal_forfft[2 * ifour];
			previous_value_e_odd = vsignal_forfft[2 * ifour + 1];
		}

	} // end loop over nfreq

	// don't applying any extra factor for the change in array (change of bin size)
	// as change in the bin size doesn't matter for the total energy
	// total energy is just same as integral over frequency range and this frequency range will not change
	//
	//cout << "ifirstnonzero, ilastnonzero are " << ifirstnonzero << " " << ilastnonzero << "\n";
	//cout << "non zero count is " << count_nonzero << "\n";
	//cout << "ratio is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
	for (int j = 0; j < NFOUR / 4; j++) {
		vsignal_forfft[2 * j] *= 1.
				/ sqrt(
						(double) count_nonzero
								/ (double) (ilastnonzero - ifirstnonzero));
		vsignal_forfft[2 * j + 1] *= 1.
				/ sqrt(
						(double) count_nonzero
								/ (double) (ilastnonzero - ifirstnonzero));
	}

	// phase for noise will be applied in GetNoiseWaveforms

}

double Report::FindPeak(double *waveform, int n) { // same with icemc, trigger-> AntTrigger::FindPeak

	double peak = 0.;

	for (int i = 0; i < n; i++) {
		if (fabs(waveform[i]) > peak)
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
	for (int i = 0; i < detector->params.number_of_stations; i++) {

		for (int j = 0; j < detector->stations[i].strings.size(); j++) {

			for (int k = 0;
					k < detector->stations[i].strings[j].antennas.size(); k++) {

				if (stations[i].strings[j].antennas[k].ray_sol_cnt) {
					if (stations[i].strings[j].antennas[k].PeakV[0] == 0.) {
						stations[i].strings[j].antennas[k].Rank.push_back(0); // rank 0, PeakV is 0, non-countable rank.
					} else {
						stations[i].strings[j].antennas[k].Rank.push_back(
								current + 1); // elses, if PeakV is not 0, set Rank as 1 (for now)
					}
				}   // if ray_sol_cnt != 0
			}
		}
	}

	//cout<<"Start while loop for Ranking!!"<<endl;

	check = 1;
	pre_maxpeak = 1.E5; // unreasonably big value which real PeakV will never reach
	while (check != 0) {
		check = 0;
		maxpeak = 0.;
		for (int i = 0; i < detector->params.number_of_stations; i++) {

			for (int j = 0; j < detector->stations[i].strings.size(); j++) {

				for (int k = 0;
						k < detector->stations[i].strings[j].antennas.size();
						k++) {

					//for (int l=0; l< detector->stations[i].strings[j].antennas[k].ray_sol_cnt; l++) {
					//
					//lets start with 1st ray_sol only
					//

					if (stations[i].strings[j].antennas[k].ray_sol_cnt) {
						if (stations[i].strings[j].antennas[k].Rank[0] != 0) { // there is non-zero value! and ray_sol_cnt also non-zero!
							if (stations[i].strings[j].antennas[k].PeakV[0]
									< pre_maxpeak) {

								if (maxpeak
										< stations[i].strings[j].antennas[k].PeakV[0]) {
									maxpeak =
											stations[i].strings[j].antennas[k].PeakV[0];

									stations[i].strings[j].antennas[k].Rank[0] =
											current + 1;
									check++;
								} // is maxpeak < PeakV
								else {
									stations[i].strings[j].antennas[k].Rank[0] =
											current + 2;
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
	if (string_num == 0 && antenna_num == 0) {
		outputch = 1;
	} else if (string_num == 0 && antenna_num == 1) {
		outputch = 2;
	} else if (string_num == 1 && antenna_num == 0) {
		outputch = 3;
	} else if (string_num == 1 && antenna_num == 1) {
		outputch = 4;
	} else if (string_num == 2 && antenna_num == 0) {
		outputch = 5;
	} else if (string_num == 2 && antenna_num == 1) {
		outputch = 6;
	} else if (string_num == 3 && antenna_num == 0) {
		outputch = 7;
	} else if (string_num == 3 && antenna_num == 1) {
		outputch = 8;
	}

	// higher antennas
	else if (string_num == 0 && antenna_num == 2) {
		outputch = 9;
	} else if (string_num == 0 && antenna_num == 3) {
		outputch = 10;
	} else if (string_num == 1 && antenna_num == 2) {
		outputch = 11;
	} else if (string_num == 1 && antenna_num == 3) {
		outputch = 12;
	} else if (string_num == 2 && antenna_num == 2) {
		outputch = 13;
	} else if (string_num == 2 && antenna_num == 3) {
		outputch = 14;
	} else if (string_num == 3 && antenna_num == 2) {
		outputch = 15;
	} else if (string_num == 3 && antenna_num == 3) {
		outputch = 16;
	}

	return outputch;

}

TGraph *Report::getWaveform(Detector *detector, int ch, int station_i, int event_num, int run_num) {

	int string_i = detector->getStringfromArbAntID(station_i, ch);
	int antenna_i = detector->getAntennafromArbAntID(station_i, ch);

	TGraph *gr = new TGraph();

	int N =
			stations[station_i].strings[string_i].antennas[antenna_i].V_mimic.size();

	for (int i = 0; i < N; i++) {

		gr->SetPoint(i,
				stations[station_i].strings[string_i].antennas[antenna_i].time_mimic[i],
				stations[station_i].strings[string_i].antennas[antenna_i].V_mimic[i]);

	}    // for i

	gr->SetNameTitle(Form("WF%dS%02dC%02d", event_num, station_i, ch),
			Form(
					"Simulated waveform %d, station %d, channel %d;time (ns); voltage (mV)",
					event_num, station_i, ch));

	return gr;

}

vector<TGraph*> Report::getWaveformVector(Detector *detector, int station_i, int event_num, int run_num) {

	vector<TGraph*> waveforms;

	for (int ch = 0; ch < detector->stations[station_i].number_of_antennas;
			ch++) {

		int string_i = detector->getStringfromArbAntID(station_i, ch);
		int antenna_i = detector->getAntennafromArbAntID(station_i, ch);

		TGraph *gr = new TGraph();

		int N =
				stations[station_i].strings[string_i].antennas[antenna_i].V_mimic.size();

		for (int i = 0; i < N; i++) {

			gr->SetPoint(i,
					stations[station_i].strings[string_i].antennas[antenna_i].time_mimic[i],
					stations[station_i].strings[string_i].antennas[antenna_i].V_mimic[i]);

		}    // for i

		gr->SetNameTitle(Form("WF%dS%02dC%02d", event_num, station_i, ch),
				Form(
						"Simulated waveform %d, station %d, channel %d;time (ns); voltage (mV)",
						event_num, station_i, ch));

		waveforms.push_back(gr);

	}    // for ch

	return waveforms;

}

vector<TGraph*> Report::getWaveformVectorVpol(Detector *detector, int station_i, int event_num, int run_num) {

	vector<TGraph*> waveforms;

	for (int ch = 0; ch < detector->stations[station_i].number_of_antennas;
			ch++) {

		int string_i = detector->getStringfromArbAntID(station_i, ch);
		int antenna_i = detector->getAntennafromArbAntID(station_i, ch);

		if (detector->stations[station_i].strings[string_i].antennas[antenna_i].type
				!= 0)
			continue; // jump to next channel if this isn't Vpol

		TGraph *gr = new TGraph();

		int N =
				stations[station_i].strings[string_i].antennas[antenna_i].V_mimic.size();

		for (int i = 0; i < N; i++) {

			gr->SetPoint(i,
					stations[station_i].strings[string_i].antennas[antenna_i].time_mimic[i],
					stations[station_i].strings[string_i].antennas[antenna_i].V_mimic[i]);

		} // for i

		gr->SetNameTitle(Form("WF%dS%02dC%02d", event_num, station_i, ch),
				Form(
						"Simulated waveform %d, station %d, channel %d;time (ns); voltage (mV)",
						event_num, station_i, ch));

		waveforms.push_back(gr);

	} // for ch

	return waveforms;

}

vector<TGraph*> Report::getWaveformVectorHpol(Detector *detector, int station_i, int event_num, int run_num) {

	vector<TGraph*> waveforms;

	for (int ch = 0; ch < detector->stations[station_i].number_of_antennas;
			ch++) {

		int string_i = detector->getStringfromArbAntID(station_i, ch);
		int antenna_i = detector->getAntennafromArbAntID(station_i, ch);

		if (detector->stations[station_i].strings[string_i].antennas[antenna_i].type
				!= 1)
			continue; // jump to next channel if this isn't Hpol

		TGraph *gr = new TGraph();

		int N =
				stations[station_i].strings[string_i].antennas[antenna_i].V_mimic.size();

		for (int i = 0; i < N; i++) {

			gr->SetPoint(i,
					stations[station_i].strings[string_i].antennas[antenna_i].time_mimic[i],
					stations[station_i].strings[string_i].antennas[antenna_i].V_mimic[i]);

		} // for i

		gr->SetNameTitle(Form("WF%dS%02dC%02d", event_num, station_i, ch),
				Form(
						"Simulated waveform %d, station %d, channel %d;time (ns); voltage (mV)",
						event_num, station_i, ch));

		waveforms.push_back(gr);

	} // for ch

	return waveforms;

}

vector<double> Report::getHitTimesVector(Detector *detector, int station_i, int polarization) {

	vector<double> times;

	for (int ch = 0; ch < detector->stations[station_i].number_of_antennas;
			ch++) {

		int string_i = detector->getStringfromArbAntID(station_i, ch);
		int antenna_i = detector->getAntennafromArbAntID(station_i, ch);

		if (polarization >= 0
				&& detector->stations[station_i].strings[string_i].antennas[antenna_i].type
						!= polarization)
			continue; // jump to next channel if this isn't Hpol/Vpol

		if (stations[station_i].strings[string_i].antennas[antenna_i].arrival_time.size()) {

			times.push_back(
					1e9
							* stations[station_i].strings[string_i].antennas[antenna_i].arrival_time[0]); // get the direct beam arrival time

		}

		else
			times.push_back(-1000); // if there's no ray-trace solution just put something in.

	}

	return times;

}

vector<double> Report::getHitTimesVectorVpol(Detector *detector, int station_i) {
	return getHitTimesVector(detector, station_i, 0);
}

vector<double> Report::getHitTimesVectorHpol(Detector *detector, int station_i) {

	return getHitTimesVector(detector, station_i, 1);

}

