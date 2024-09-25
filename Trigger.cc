#include "Trigger.h"
#include "Detector.h"
#include "Settings.h"
#include "Report.h"
#include "Efficiencies.h"
#include "Tools.h"
#include "Constants.h"
#include <math.h>
#include <vector>

ClassImp(Trigger);

Trigger::Trigger() {

}


Trigger::~Trigger() {

    v_noise_timedomain.clear();
    v_noise_timedomain_diode.clear();

    v_noise_timedomain_ch.clear();
    v_noise_timedomain_diode_ch.clear();

    V_noise_freqbin_ch.clear();

     meandiode_ch.clear();
     rmsdiode_ch.clear();
     rmsvoltage_ch.clear();

}



void Trigger::ClearNoiseWaveforms() {

    v_noise_timedomain.clear();
    v_noise_timedomain_diode.clear();

    v_noise_timedomain_ch.clear();
    v_noise_timedomain_diode_ch.clear();

}


Trigger::Trigger(Detector *detector, Settings *settings1) {

    TIMESTEP = detector->TIMESTEP;
    maxt_diode = detector->maxt_diode;
    maxt_diode_bin = (int)(maxt_diode/TIMESTEP);
    NFOUR = detector->NFOUR;
    DATA_BIN_SIZE = settings1->DATA_BIN_SIZE;

    powerthreshold = settings1->POWERTHRESHOLD;

    iminbin = NFOUR/4 - detector->ibinshift + detector->idelaybeforepeak;
    imaxbin = NFOUR/4 - detector->ibinshift + detector->idelaybeforepeak + detector->iwindow;


    if (settings1->NOISE_CHANNEL_MODE == 0) {
    
        V_noise_freqbin = sqrt( (double)(settings1->DATA_BIN_SIZE) * 50. * KBOLTZ * settings1->NOISE_TEMP / (settings1->TIMESTEP * 2.) );
    }

    else if (settings1->NOISE_CHANNEL_MODE == 1) {
      cout << "Detector.params.numberofantennas: " << detector->params.number_of_antennas << endl;
        for (int ch=0; ch<detector->params.number_of_antennas; ch++) {
            //V_noise_freqbin_ch.push_back( sqrt( (double)(settings1->DATA_BIN_SIZE) * 50. * KBOLTZ * settings1->NOISE_TEMP / (settings1->TIMESTEP * 2.) ) );
            V_noise_freqbin_ch.push_back( sqrt( (double)(settings1->DATA_BIN_SIZE) * 50. * KBOLTZ * detector->GetTemp(0, ch, settings1) / (settings1->TIMESTEP * 2.) ) );
        }
    }

    // mode 2 : make separate noise waveforms for just first 8 chs
    else if (settings1->NOISE_CHANNEL_MODE == 2) {
    
        for (int ch=0; ch<8; ch++) {
            V_noise_freqbin_ch.push_back( sqrt( (double)(settings1->DATA_BIN_SIZE) * 50. * KBOLTZ * detector->GetTemp(0, ch, settings1) / (settings1->TIMESTEP * 2.) ) );
        }
        
        // total 9 array, last array for shared systemp
        V_noise_freqbin_ch.push_back( sqrt( (double)(settings1->DATA_BIN_SIZE) * 50. * KBOLTZ * settings1->NOISE_TEMP / (settings1->TIMESTEP * 2.) ) );

    }

    // Load in Phased Array Data
    if (settings1->TRIG_SCAN_MODE==5){

        cout << "Phased Array mode! Reading in data: " << endl;
    
        // Load Efficiency vs SNR curve: 
        // From arxiv.1809.04573, Fig 15, Dashed Pink 0.8Hz/beam curve
        ifstream infile;
        infile.open("data/nuphase_trig_effc.txt",ios::in);
        if(infile.fail()){ // checks to see if file opended
            cerr << "   error, file could not be opened" << endl;
        }
        int num = 0;
        while(!infile.eof()) { // reads file to end of *file*, not line
            infile >> snr_PA[num]  // read first column number
                   >> eff_PA[num]; // read second column number
            ++num;
        }
        cout<<"   number of triggering efficiency data points: "<<num<<endl;
        infile.close();

        //load Angle vs SNR curve:
        // KAH says this curve is created by converting Fig 13b from 
        //   arxiv.1809.04573 into digitized, linear units but ARB has not 
        //   been able to replicate this
        infile.clear();
        infile.open("data/nuphase_SNR_angle.txt",ios::in);
        if(infile.fail()) { // checks to see if file opended
            cerr << "   error, file could not be opened" << endl;
        }
        num = 0;
        while(!infile.eof()) { // reads file to end of *file*, not line
            infile >> angle_PA[num]   // read first column number
                    >>  aSNR_PA[num];  // read second column number
            ++num;
        }
        cout<<"   number of SNR vs angle data points: "<<num<<endl;
        infile.close();

    } // End if DETECTOR=5 (Phased Array sim)

}



void Trigger::Reset_V_noise_freqbin(Settings *settings1, Detector *detector) {


    V_noise_freqbin_ch.clear();


    if (settings1->NOISE_CHANNEL_MODE == 0) {
    
        V_noise_freqbin = sqrt( (double)(settings1->DATA_BIN_SIZE) * 50. * KBOLTZ * settings1->NOISE_TEMP / (settings1->TIMESTEP * 2.) );
    }

    else if (settings1->NOISE_CHANNEL_MODE == 1) {
    
        for (int ch=0; ch<detector->params.number_of_antennas; ch++) {
            //V_noise_freqbin_ch.push_back( sqrt( (double)(settings1->DATA_BIN_SIZE) * 50. * KBOLTZ * settings1->NOISE_TEMP / (settings1->TIMESTEP * 2.) ) );
            V_noise_freqbin_ch.push_back( sqrt( (double)(settings1->DATA_BIN_SIZE) * 50. * KBOLTZ * detector->GetTemp(0, ch, settings1) / (settings1->TIMESTEP * 2.) ) );
        }
    }

    // mode 2 : make separate noise waveforms for just first 8 chs
    else if (settings1->NOISE_CHANNEL_MODE == 2) {
    
        for (int ch=0; ch<8; ch++) {
            V_noise_freqbin_ch.push_back( sqrt( (double)(settings1->DATA_BIN_SIZE) * 50. * KBOLTZ * detector->GetTemp(0, ch, settings1) / (settings1->TIMESTEP * 2.) ) );
        }
        
        // total 9 array, last array for shared systemp
        V_noise_freqbin_ch.push_back( sqrt( (double)(settings1->DATA_BIN_SIZE) * 50. * KBOLTZ * settings1->NOISE_TEMP / (settings1->TIMESTEP * 2.) ) );

    }

}




void Trigger::SetMeanRmsDiode(Settings *settings1, Detector *detector, Report *report) {
//void Trigger::SetMeanRmsDiode(Detector *detector, vector <double> &v_noise_timedomain, vector <double> &v_noise_timedomain_diode) {

    // set meandiode
    // same with icemc -> anita -> Initialize
    cerr<<"Preparing Noise"<<endl;

    // if using default noise temp setting (same temp for all chs)
    if (settings1->NOISE_CHANNEL_MODE == 0) {
        
        int ngeneratedevents=settings1->NOISE_EVENTS;  // should this value read at Settings class
        double v_noise[settings1->DATA_BIN_SIZE];    // noise voltage time domain (with filter applied)
        
        // Prepare to save data
        meandiode = 0.;
        rmsdiode = 0.;
        rmsvoltage = 0.;
        
        // make the size of v_noise_timedomain_diode ngeneratedevents (this will be huge!)
        v_noise_timedomain.resize(ngeneratedevents);  
        v_noise_timedomain_diode.resize(ngeneratedevents);
        
        // Get all `ngeneratedevents` noise waveforms
        if ( ngeneratedevents > 10 ) {
            cerr<<"Generating noise waveforms"<<endl;
        }
        for (int i=0; i<ngeneratedevents; i++) {

            // Print updates to cerr
            if ( ngeneratedevents > 10 ){
                if ((i)%(ngeneratedevents/10) == 0) {
                    cerr<< (i/ngeneratedevents)*100 <<"% done"<<endl;
                }
            }

            // get v_noise array (noise voltage in time domain)
            report->GetNoiseWaveforms(settings1, detector, V_noise_freqbin, v_noise);

            // do normal time ordering (not sure if this is necessary)
            Tools::NormalTimeOrdering(settings1->DATA_BIN_SIZE, v_noise);

            // Convolve signal through the tunnel diode
            myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real_databin,v_noise_timedomain_diode[i]);

            // Save noise waveform and determine mean noise value in this channel
            for (int m=0; m<settings1->DATA_BIN_SIZE; m++) {
                
                // Add contribution from this bin to the mean diode value we're calculating
                if ( m>=(int)(maxt_diode/TIMESTEP) && m<settings1->DATA_BIN_SIZE ) {
                    meandiode += (
                        v_noise_timedomain_diode[i][m]
                        / ( (double)ngeneratedevents * ( (double)settings1->DATA_BIN_SIZE - maxt_diode/TIMESTEP ) ) );
                } // end loop over samples where diode function is fully contained 

                // save pure noise (not diode convlved) waveforms
                v_noise_timedomain[i].push_back(v_noise[m]);

                // save pure noise spectrum before applying Rayleigh distribution (only at first evt as they will be all same)
                if ( i == 0 && m<settings1->DATA_BIN_SIZE/2 ) {
                    Vfft_noise_before.push_back(report->Vfft_noise_before[m]);
                }

            }

        }   // get meandiode with `ngeneratedevents` noisewaveforms;
        if ( ngeneratedevents > 10 ) {
            cerr<< "100% done"<<endl;
        }
        
        // now as v_noise_timedomain_diode's waveforms are still stored, we can just calculate rms from them
        if ( ngeneratedevents > 10 ) {
            cerr<<"Calculating noise waveforms for RMS"<<endl;
        }
        for (int i=0; i<ngeneratedevents; i++) {

            // Print updates to cerr
            if ( ngeneratedevents > 10 ){
                if ((i)%(ngeneratedevents/10) == 0) {
                    cerr<< (i/ngeneratedevents)*100 <<"% done"<<endl;
                }
            }

            // we are going to reuse generated v_noise_timedomain_diode above.
            //report->GetNoiseWaveforms(settings1, detector, V_noise_freqbin, v_noise);
            //myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real,v_noise_timedomain_diode);

            // Add contribution from this bin to the mean RMS we're calculating
            for (int m=(int)(maxt_diode/TIMESTEP); m<settings1->DATA_BIN_SIZE; m++) {
                rmsdiode += (
                    (v_noise_timedomain_diode[i][m]-meandiode) * (v_noise_timedomain_diode[i][m]-meandiode)
                    /( (double)ngeneratedevents * ( (double)settings1->DATA_BIN_SIZE - maxt_diode/TIMESTEP ) ) );
                rmsvoltage += (
                    (v_noise_timedomain[i][m]) * (v_noise_timedomain[i][m])
                    /( (double)ngeneratedevents * ( (double)settings1->DATA_BIN_SIZE - maxt_diode/TIMESTEP ) ) );
            }

        }   // get rmsdiode with `ngeneratedevents` noisewaveforms
        cerr<< "100% done"<<endl;

        // Finish RMS calculation and print out
        rmsdiode=sqrt(rmsdiode);
        rmsvoltage=sqrt(rmsvoltage);
        cout << "From pure noise waveforms, diode responses" << "\n";
        cout << "mean, rms diode are " << meandiode << " " << rmsdiode << "\n";
        cout << "rms voltage is "<<rmsvoltage<<"\n";
        cout<<" DATA_BIN_SIZE : "<<settings1->DATA_BIN_SIZE<<"\n";

        // if we are doing pure signal trigger analysis, set all noise waveform values to 0
        if (settings1->TRIG_ANALYSIS_MODE == 1 ) {
            for (int i=0; i<ngeneratedevents; i++) {
                for (int m=0; m<settings1->DATA_BIN_SIZE; m++) {
                    v_noise_timedomain[i][m] = 0.;
                    v_noise_timedomain_diode[i][m] = 0.;
                }
            }
        }

    }// if NOISE_CHANNEL_MODE = 0

    // if using mode 1 noise temp setting (different temp for each chs)
    else if (settings1->NOISE_CHANNEL_MODE == 1) {

        int ngeneratedevents=settings1->NOISE_EVENTS;  // should this value read at Settings class
        double v_noise[settings1->DATA_BIN_SIZE];    // noise voltage time domain (with filter applied)

        int num_chs = detector->params.number_of_antennas;
        cout << "num chs: " << num_chs << endl;

        // Prepare to save data for all `num_chs` channels
        meandiode_ch.resize(num_chs);
        rmsdiode_ch.resize(num_chs);
        rmsvoltage_ch.resize(num_chs);
        v_noise_timedomain_ch.resize(num_chs);
        v_noise_timedomain_diode_ch.resize(num_chs);

        for (int i=0; i<num_chs; i++) {

            // make the size of v_noise_timedomain_diode[ch] ngeneratedevents (this will be huge!)
            v_noise_timedomain_ch[i].resize(ngeneratedevents);
            v_noise_timedomain_diode_ch[i].resize(ngeneratedevents);

            meandiode_ch[i] = 0.;
            rmsdiode_ch[i] = 0.;
            rmsvoltage_ch[i] = 0.;

        }

        // Get all `ngeneratedevents` noise waveforms for `num_chs` channels
        if ( ngeneratedevents > 10 ) {
            cerr<<"Generating noise waveforms"<<endl;
        }
        for (int ch=0; ch<num_chs; ch++) {

            for (int i=0; i<ngeneratedevents; i++) {

                // Print updates to cerr
                if ( ngeneratedevents > 10 ){
                    if ((i)%(ngeneratedevents/10) == 0) {
                        cerr<< (i/ngeneratedevents)*100 <<"% done for ch"<<ch<<endl;
                    }
                }

                // get v_noise array (noise voltage in time domain)
                report->GetNoiseWaveforms_ch(settings1, detector, V_noise_freqbin_ch[ch], v_noise, ch);

                // do normal time ordering (not sure if this is necessary)
                Tools::NormalTimeOrdering(settings1->DATA_BIN_SIZE, v_noise);

                // Convolve signal through the tunnel diode
                myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real_databin,v_noise_timedomain_diode_ch[ch][i]);

                // Save noise waveform and determine mean noise value in this channel
                for (int m=0; m<settings1->DATA_BIN_SIZE; m++) {

                    // Add contribution from this bin to the mean diode value we're calculating
                    if ( m>=(int)(maxt_diode/TIMESTEP) && m<settings1->DATA_BIN_SIZE ) {
                        meandiode_ch[ch] += (
                            v_noise_timedomain_diode_ch[ch][i][m]
                            /( (double)ngeneratedevents * ( (double)settings1->DATA_BIN_SIZE - maxt_diode/TIMESTEP ) ) );

                    } // end loop over samples where diode function is fully contained 

                    // save pure noise (not diode convlved) waveforms
                    v_noise_timedomain_ch[ch][i].push_back(v_noise[m]);

                    // save pure noise spectrum before applying Rayleigh distribution (only at first evt as they will be all same)
                    if ( i == 0 && m<settings1->DATA_BIN_SIZE/2 ) {
                        Vfft_noise_before.push_back(report->Vfft_noise_before[m]);
                    }

                }

            }   // get meandiode with `ngeneratedevents` noisewaveforms

        } // loop over chs
        if ( ngeneratedevents > 10 ) {
            cerr<< "100% done"<<endl;
        }

        // now as v_noise_timedomain_diode's waveforms are still stored, we can just calculate rms from them
        if ( ngeneratedevents > 10 ) {
            cerr<<"Calculating noise waveforms RMS"<<endl;
        }
        for (int ch=0; ch<num_chs; ch++) {

            for (int i=0; i<ngeneratedevents; i++) {

                // Print updates to cerr
                if ( ngeneratedevents > 10 ){
                    if ((i)%(ngeneratedevents/10) == 0) {
                        cerr<< (i/ngeneratedevents)*100 <<"% done for ch"<<ch<<endl;
                    }
                }

                // Add contribution from this bin to the mean RMS we're calculating
                for (int m=(int)(maxt_diode/TIMESTEP);m<settings1->DATA_BIN_SIZE;m++) {
                    rmsdiode_ch[ch] += (
                        (v_noise_timedomain_diode_ch[ch][i][m]-meandiode_ch[ch]) * (v_noise_timedomain_diode_ch[ch][i][m]-meandiode_ch[ch])
                        /( (double)ngeneratedevents * ( (double)settings1->DATA_BIN_SIZE - maxt_diode/TIMESTEP ) ) );
                    rmsvoltage_ch[ch]+=( 
                        (v_noise_timedomain_ch[ch][i][m]) * (v_noise_timedomain_ch[ch][i][m])
                        /( (double)ngeneratedevents * ( (double)settings1->DATA_BIN_SIZE - maxt_diode/TIMESTEP ) ) );
                }

            }   // get rmsdiode with `ngeneratedevents` noisewaveforms
        }
        if ( ngeneratedevents > 10 ) {
            cerr<< "100% done"<<endl;
        }

        // Finish RMS calculation and print out
        cout << "From pure noise waveforms, diode responses" << "\n";
        for (int ch=0; ch<num_chs; ch++) {
            rmsdiode_ch[ch]=sqrt(rmsdiode_ch[ch]);
            rmsvoltage_ch[ch]=sqrt(rmsvoltage_ch[ch]);
            cout << "For ch"<<ch<<" mean, rms diode are " << meandiode_ch[ch] << " " << rmsdiode_ch[ch] << " rms voltage is "<<rmsvoltage_ch[ch]<<"\n";
        }
        cout<<" DATA_BIN_SIZE : "<<settings1->DATA_BIN_SIZE<<"\n";

        // if we are doing pure signal trigger analysis, set all noise waveform values to 0
        if (settings1->TRIG_ANALYSIS_MODE == 1 ) {
            for (int ch=0; ch<num_chs; ch++) {
                for (int i=0; i<ngeneratedevents; i++) {
                    for (int m=0; m<settings1->DATA_BIN_SIZE; m++) {
                        v_noise_timedomain_ch[ch][i][m] = 0.;
                        v_noise_timedomain_diode_ch[ch][i][m] = 0.;
                    }
                }
            }
        }

    }// else if NOISE_CHANNEL_MODE = 1

    // if using mode 2 noise temp setting (different temp for first 8 chs)
    else if (settings1->NOISE_CHANNEL_MODE == 2) {

        int ngeneratedevents=settings1->NOISE_EVENTS;  // should this value read at Settings class
        double v_noise[settings1->DATA_BIN_SIZE];    // noise voltage time domain (with filter applied)

        //int num_chs = detector->params.number_of_antennas;
        int num_chs = 8+1;// 8 chs for separated systemp temp and one more for sharing temp

        // Prepare to save data for all `num_chs` channels
        meandiode_ch.resize(num_chs);
        rmsdiode_ch.resize(num_chs);
        rmsvoltage_ch.resize(num_chs);
        v_noise_timedomain_ch.resize(num_chs);
        v_noise_timedomain_diode_ch.resize(num_chs);

        for (int i=0; i<num_chs; i++) {

            // make the size of v_noise_timedomain_diode[ch] ngeneratedevents (this will be huge!)
            v_noise_timedomain_ch[i].resize(ngeneratedevents);
            v_noise_timedomain_diode_ch[i].resize(ngeneratedevents);

            meandiode_ch[i] = 0.;
            rmsdiode_ch[i] = 0.;
            rmsvoltage_ch[i] = 0.;

        }

        // Get all `ngeneratedevents` noise waveforms for `num_chs` channels
        if ( ngeneratedevents > 10 ) {
            cerr<<"Generating noise waveforms"<<endl;
        }
        for (int ch=0; ch<num_chs; ch++) {

            for (int i=0; i<ngeneratedevents; i++) {

                // Print updates to cerr
                if ( ngeneratedevents > 10 ){
                    if ((i)%(ngeneratedevents/10) == 0) {
                        cerr<< (i/ngeneratedevents)*100 <<"% done for ch"<<ch<<endl;
                    }
                }

                // get v_noise array (noise voltage in time domain)
                report->GetNoiseWaveforms_ch(settings1, detector, V_noise_freqbin_ch[ch], v_noise, ch);

                // do normal time ordering (not sure if this is necessary)
                Tools::NormalTimeOrdering(settings1->DATA_BIN_SIZE, v_noise);

                // Convolve signal through the tunnel diode
                myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real_databin,v_noise_timedomain_diode_ch[ch][i]);

                // Save noise waveform and determine mean noise value in this channel
                for (int m=0; m<settings1->DATA_BIN_SIZE; m++) {

                    // Add contribution from this bin to the mean diode value we're calculating
                    if ( m>=(int)(maxt_diode/TIMESTEP) && m<settings1->DATA_BIN_SIZE ) {
                        meandiode_ch[ch] += (
                            v_noise_timedomain_diode_ch[ch][i][m]
                            /( (double)ngeneratedevents * ( (double)settings1->DATA_BIN_SIZE - maxt_diode/TIMESTEP ) ) );

                    } // end loop over samples where diode function is fully contained 

                    // save pure noise (not diode convlved) waveforms
                    v_noise_timedomain_ch[ch][i].push_back(v_noise[m]);

                    // save pure noise spectrum before applying Rayleigh distribution (only at first evt as they will be all same)
                    if ( i == 0 && m<settings1->DATA_BIN_SIZE/2 ) {
                        Vfft_noise_before.push_back(report->Vfft_noise_before[m]);
                    }

                }

            }   // get meandiode with `ngeneratedevents` noisewaveforms

        } // loop over chs
        if ( ngeneratedevents > 10 ) {
            cerr<< "100% done"<<endl;
        }

        // now as v_noise_timedomain_diode's waveforms are still stored, we can just calculate rms from them
        if ( ngeneratedevents > 10 ) {
            cerr<<"Calculating noise waveforms RMS"<<endl;
        }
        for (int ch=0; ch<num_chs; ch++) {

            for (int i=0; i<ngeneratedevents; i++) {

                // Print updates to cerr
                if ( ngeneratedevents > 10 ){
                    if ((i)%(ngeneratedevents/10) == 0) {
                        cerr<< (i/ngeneratedevents)*100 <<"% done for ch"<<ch<<endl;
                    }
                }

                // Add contribution from this bin to the mean RMS we're calculating
                for (int m=(int)(maxt_diode/TIMESTEP);m<settings1->DATA_BIN_SIZE;m++) {
                    rmsdiode_ch[ch]+=(
                        (v_noise_timedomain_diode_ch[ch][i][m]-meandiode_ch[ch]) * (v_noise_timedomain_diode_ch[ch][i][m]-meandiode_ch[ch])
                        /( (double)ngeneratedevents * ( (double)settings1->DATA_BIN_SIZE - maxt_diode/TIMESTEP ) ) );
                    rmsvoltage_ch[ch]+=(
                        (v_noise_timedomain_ch[ch][i][m]) * (v_noise_timedomain_ch[ch][i][m])
                        /( (double)ngeneratedevents * ( (double)settings1->DATA_BIN_SIZE - maxt_diode/TIMESTEP ) ) );
                }
            }   // get rmsdiode with `ngeneratedevents` noisewaveforms
        }
        if ( ngeneratedevents > 10 ) {
            cerr<< "100% done"<<endl;
        }

        // Finish RMS calculation and print out
        cout << "From pure noise waveforms, diode responses" << "\n";
        for (int ch=0; ch<num_chs; ch++) {
            rmsdiode_ch[ch]=sqrt(rmsdiode_ch[ch]);
            rmsvoltage_ch[ch]=sqrt(rmsvoltage_ch[ch]);
            cout << "For ch"<<ch<<" mean, rms diode are " << meandiode_ch[ch] << " " << rmsdiode_ch[ch] << " rms voltage is "<<rmsvoltage_ch[ch]<<"\n";
        }
        cout<<" DATA_BIN_SIZE : "<<settings1->DATA_BIN_SIZE<<"\n";

        // if we are doing pure signal trigger analysis, set all noise waveform values to 0
        if (settings1->TRIG_ANALYSIS_MODE == 1 ) {
            for (int ch=0; ch<num_chs; ch++) {
                for (int i=0; i<ngeneratedevents; i++) {
                    for (int m=0; m<settings1->DATA_BIN_SIZE; m++) {
                        v_noise_timedomain_ch[ch][i][m] = 0.;
                        v_noise_timedomain_diode_ch[ch][i][m] = 0.;
                    }
                }
            }
        }

    }// else if NOISE_CHANNEL_MODE = 2

}

double Trigger::GetAntNoise_diodeMean(int ch_ID, Settings *settings1){
    // Return the mean value of noise in the diode for the requested antenna

    if (settings1->NOISE_CHANNEL_MODE == 0){
        // All antennas use the same noise WFs
        return meandiode;
    }
    else if (settings1->NOISE_CHANNEL_MODE == 1){
        // All antennas have their own noise WFs
        return meandiode_ch[ch_ID];
    }
    else if (settings1->NOISE_CHANNEL_MODE == 2){
        // Only deep channels (ch_ID<8) use their own noise WFs. The others share one.
        if (ch_ID < 8){
            return meandiode_ch[ch_ID];
        }
        else {
            return meandiode_ch[8];
        }
    }
    else{
        cerr<<"Unsupported NOISE_CHANNEL_MODE of "<<settings1->NOISE_CHANNEL_MODE;
        throw("Cannot return antenna's noise diode mean.");
    }

}

double Trigger::GetAntNoise_diodeRMS(int ch_ID, Settings *settings1){
    // Return the RMS value of noise in the diode for the requested antenna

    if (settings1->NOISE_CHANNEL_MODE == 0){
        // All antennas use the same noise WFs
        return rmsdiode;
    }
    else if (settings1->NOISE_CHANNEL_MODE == 1){
        // All antennas have their own noise WFs
        return rmsdiode_ch[ch_ID];
    }
    else if (settings1->NOISE_CHANNEL_MODE == 2){
        // Only deep channels (ch_ID<8) use their own noise WFs. The others share one.
        if (ch_ID < 8){
            return rmsdiode_ch[ch_ID];
        }
        else {
            return rmsdiode_ch[8];
        }
    }
    else{
        cerr<<"Unsupported NOISE_CHANNEL_MODE of "<<settings1->NOISE_CHANNEL_MODE;
        throw("Cannot return antenna's noise diode RMS.");
    }

}


double Trigger::GetAntNoise_voltageRMS(int ch_ID, Settings *settings1){
    // Return the RMS value of noise in volts for the requested antenna

    if (settings1->NOISE_CHANNEL_MODE == 0){
        // All antennas use the same noise WFs
        return rmsvoltage;
    }
    else if (settings1->NOISE_CHANNEL_MODE == 1){
        // All antennas have their own noise WFs
        return rmsvoltage_ch[ch_ID];
    }
    else if (settings1->NOISE_CHANNEL_MODE == 2){
        // Only deep channels (ch_ID<8) use their own noise WFs. The others share one.
        if (ch_ID < 8){
            return rmsvoltage_ch[ch_ID];
        }
        else {
            return rmsvoltage_ch[8];
        }
    }
    else{
        cerr<<"Unsupported NOISE_CHANNEL_MODE of "<<settings1->NOISE_CHANNEL_MODE;
        throw("Cannot return antenna's noise voltage RMS.");
    }
    
}


void Trigger::GetNewNoiseWaveforms(Settings *settings1, Detector *detector, Report *report)
{
    // get some new noise waveforms
    // if using default noise temp setting (same temp for all chs)
    if (settings1->NOISE_CHANNEL_MODE == 0)
    {

        int ngeneratedevents = settings1->NOISE_EVENTS;  // should this value read at Settings class
        double v_noise[settings1->DATA_BIN_SIZE];  // noise voltage time domain (with filter applied)

        v_noise_timedomain.resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)
        v_noise_timedomain_diode.resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)

        if (settings1->TRIG_ANALYSIS_MODE != 1)
        {

            for (int i = 0; i < ngeneratedevents; i++)
            {

                // get v_noise array (noise voltage in time domain)
                report->GetNoiseWaveforms(settings1, detector, V_noise_freqbin, v_noise);

                // do normal time ordering (not sure if this is necessary)
                Tools::NormalTimeOrdering(settings1->DATA_BIN_SIZE, v_noise);

                myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real_databin, v_noise_timedomain_diode[i]);

                for (int m = 0; m < settings1->DATA_BIN_SIZE; m++)
                {

                    v_noise_timedomain[i].push_back(v_noise[m]);  // save pure noise (not diode convlved) waveforms

                }
            }  // get meandiode with 1000 noisewaveforms
        }

        // if we are doing pure signal trigger analysis, set all noise waveform values to 0
        else if (settings1->TRIG_ANALYSIS_MODE == 1)
        {

            for (int i = 0; i < ngeneratedevents; i++)
            {

                for (int m = 0; m < settings1->DATA_BIN_SIZE; m++)
                {

                    v_noise_timedomain[i].push_back(0.);
                    v_noise_timedomain_diode[i].push_back(0.);
                }
            }
        }
    }  // if NOISE_CHANNEL_MODE = 0
    
    // if using mode 1 noise temp setting (different temp for each chs)
    else if (settings1->NOISE_CHANNEL_MODE == 1)
    {
        int ngeneratedevents = settings1->NOISE_EVENTS;  // should this value read at Settings class
        double v_noise[settings1->DATA_BIN_SIZE];  // noise voltage time domain (with filter applied)

        int num_chs = detector->params.number_of_antennas;

        v_noise_timedomain_ch.resize(num_chs);  // make the size of v_noise_timedomain_diode as number of chs
        v_noise_timedomain_diode_ch.resize(num_chs);  // make the size of v_noise_timedomain_diode as number of chs

        for (int i = 0; i < num_chs; i++)
        {
            v_noise_timedomain_ch[i].resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)
            v_noise_timedomain_diode_ch[i].resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)

        }

        if (settings1->TRIG_ANALYSIS_MODE != 1)
        {

            for (int ch = 0; ch < num_chs; ch++)
            {

                for (int i = 0; i < ngeneratedevents; i++)
                {

                    // get v_noise array (noise voltage in time domain)
                    report->GetNoiseWaveforms_ch(settings1, detector, V_noise_freqbin_ch[ch], v_noise, ch);

                    // cout << "After getting noise waveforms" << endl;

                    // do normal time ordering (not sure if this is necessary)
                    Tools::NormalTimeOrdering(settings1->DATA_BIN_SIZE, v_noise);

                    myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real_databin, v_noise_timedomain_diode_ch[ch][i]);

                    // cout << "After convolve" << endl;

                    for (int m = 0; m < settings1->DATA_BIN_SIZE; m++)
                    {

                        v_noise_timedomain_ch[ch][i].push_back(v_noise[m]);  // save pure noise (not diode convlved) waveforms

                    }

                    // cout << "after push back " << endl;
                }  // get meandiode with 1000 noisewaveforms

            }  // loop over chs

        }

        // if we are doing pure signal trigger analysis, set all noise waveform values to 0
        else if (settings1->TRIG_ANALYSIS_MODE == 1)
        {

            for (int ch = 0; ch < num_chs; ch++)
            {

                for (int i = 0; i < ngeneratedevents; i++)
                {

                    for (int m = 0; m < settings1->DATA_BIN_SIZE; m++)
                    {

                        v_noise_timedomain_ch[ch][i].push_back(0.);
                        v_noise_timedomain_diode_ch[ch][i].push_back(0.);
                    }
                }
            }
        }
    }  // else if NOISE_CHANNEL_MODE = 1

    // if using mode 2 noise temp setting (different temp for first 8 chs)
    else if (settings1->NOISE_CHANNEL_MODE == 2)
    {

        int ngeneratedevents = settings1->NOISE_EVENTS;  // should this value read at Settings class
        double v_noise[settings1->DATA_BIN_SIZE];  // noise voltage time domain (with filter applied)

        // int num_chs = detector->params.number_of_antennas;
        int num_chs = 8 + 1; // 8 chs for separated systemp temp and one more for sharing temp

        v_noise_timedomain_ch.resize(num_chs);  // make the size of v_noise_timedomain_diode as number of chs
        v_noise_timedomain_diode_ch.resize(num_chs);  // make the size of v_noise_timedomain_diode as number of chs

        for (int i = 0; i < num_chs; i++)
        {
            v_noise_timedomain_ch[i].resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)
            v_noise_timedomain_diode_ch[i].resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)

        }

        if (settings1->TRIG_ANALYSIS_MODE != 1)
        {

            for (int ch = 0; ch < num_chs; ch++)
            {

                for (int i = 0; i < ngeneratedevents; i++)
                {

                    // get v_noise array (noise voltage in time domain)
                    report->GetNoiseWaveforms_ch(settings1, detector, V_noise_freqbin_ch[ch], v_noise, ch);

                    // do normal time ordering (not sure if this is necessary)
                    Tools::NormalTimeOrdering(settings1->DATA_BIN_SIZE, v_noise);

                    myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real_databin, v_noise_timedomain_diode_ch[ch][i]);

                    for (int m = 0; m < settings1->DATA_BIN_SIZE; m++)
                    {

                        v_noise_timedomain_ch[ch][i].push_back(v_noise[m]); // save pure noise (not diode convlved) waveforms

                    }
                }  // get meandiode with 1000 noisewaveforms

            }  // loop over chs

        }

        // if we are doing pure signal trigger analysis, set all noise waveform values to 0
        else if (settings1->TRIG_ANALYSIS_MODE == 1)
        {

            for (int ch = 0; ch < num_chs; ch++)
            {

                for (int i = 0; i < ngeneratedevents; i++)
                {

                    for (int m = 0; m < settings1->DATA_BIN_SIZE; m++)
                    {

                        v_noise_timedomain_ch[ch][i].push_back(0.);
                        v_noise_timedomain_diode_ch[ch][i].push_back(0.);
                    }
                }
            }
        }
    }  // else if NOISE_CHANNEL_MODE = 2

}




    


// myconvlv in AraSim is actually different with myconvlv in icemc!!!

// moved mindiodeconvl, onediodeconvl, power_noise as I can't find their use.
//void Trigger::myconvlv(double *data,const int NFOUR,double *fdiode,double &mindiodeconvl,double &onediodeconvl,double *power_noise,double *diodeconv) {
void Trigger::myconvlv(vector <double> &data,const int DATA_BIN_SIZE,vector <double> &fdiode, vector <double> &diodeconv) {
    
    
    const int length=DATA_BIN_SIZE;
//    double data_copy[length];
    //double fdiode_real[length];

    // we are going to make double sized array for complete convolution
    double power_noise_copy[length*2];

/*    
    for (int i=0;i<NFOUR/2;i++) {
	data_copy[i]=data[i];
    }
    for (int i=NFOUR/2;i<length;i++) {
	data_copy[i]=0.;
    }
*/  


    // fill half of the array as power (actually energy) and another half (actually extanded part) with zero padding (Numerical Recipes 643 page)
    for (int i=0;i<length;i++) {
      power_noise_copy[i]=(data[i]*data[i])/Zr*TIMESTEP;
    }
    for (int i=length;i<length*2;i++) {
      power_noise_copy[i]=0.;
    }
    
    
    
    // do forward fft to get freq domain (energy of pure signal)
    Tools::realft(power_noise_copy,1,length*2);
    
    double ans_copy[length*2];
    
    
    
    // change the sign (from numerical recipes 648, 649 page)
    for (int j=1;j<length;j++) {
      ans_copy[2*j]=(power_noise_copy[2*j]*fdiode[2*j]-power_noise_copy[2*j+1]*fdiode[2*j+1])/((double)length);
      //ans_copy[2*j]=(power_noise_copy[2*j]*fdiode[2*j]+power_noise_copy[2*j+1]*fdiode[2*j+1])/((double)length/2);
      ans_copy[2*j+1]=(power_noise_copy[2*j+1]*fdiode[2*j]+power_noise_copy[2*j]*fdiode[2*j+1])/((double)length);
      //ans_copy[2*j+1]=(power_noise_copy[2*j+1]*fdiode[2*j]-power_noise_copy[2*j]*fdiode[2*j+1])/((double)length/2);
    }
    ans_copy[0]=power_noise_copy[0]*fdiode[0]/((double)length);
    ans_copy[1]=power_noise_copy[1]*fdiode[1]/((double)length);

    // 1/length is actually 2/(length * 2)
    //
    
    Tools::realft(ans_copy,-1,length*2);
    
    diodeconv.clear();  // remove previous values in diodeconv
    
    // only save first half of convolution result as last half is result from zero padding (and some spoiled bins)
    //for (int i=0;i<length;i++) {
    for (int i=0;i<length+maxt_diode_bin;i++) {
      //power_noise[i]=power_noise_copy[i];
      //diodeconv[i]=ans_copy[i];
      diodeconv.push_back( ans_copy[i] );
    }
    
    
    
    int iminsamp,imaxsamp; // find min and max samples such that
    // the diode response is fully overlapping with the noise waveform
    iminsamp=(int)(maxt_diode/TIMESTEP);
    // the noise waveform is NFOUR/2 long
    // then for a time maxt_diode/TIMESTEP at the end of that, the
    // diode response function is only partially overlappying with the
    // waveform in the convolution
    imaxsamp=NFOUR/2;
    
    //  cout << "iminsamp, imaxsamp are " << iminsamp << " " << imaxsamp << "\n";
    if (imaxsamp<iminsamp) {
      cout << "Noise waveform is not long enough for this diode response.\n";
      exit(1);
    }
    
    //int ibin=((iminsamp+imaxsamp)-(iminsamp+imaxsamp)%2)/2;
    //cout << "ibin is " << ibin << "\n";
    
    // return the 50th sample, right in the middle
   // onediodeconvl=diodeconv[ibin];
    
    
    //mindiodeconvl=0.;
    
    //for (int i=iminsamp;i<imaxsamp;i++) {
	
//	if (diodeconv[i]<mindiodeconvl)
//	    mindiodeconvl=diodeconv[i];
	
	
  //  }
    
}



// input data is not vector but double array
void Trigger::myconvlv(double *data,const int DATA_BIN_SIZE,vector <double> &fdiode, vector <double> &diodeconv) {
    
    // cout << "Trigger::myconvlv" << endl;
    const int length=DATA_BIN_SIZE;
//    double data_copy[length];
    //double fdiode_real[length];

    // we are going to make double sized array for complete convolution
    double power_noise_copy[length*2];

/*    
    for (int i=0;i<NFOUR/2;i++) {
	data_copy[i]=data[i];
    }
    for (int i=NFOUR/2;i<length;i++) {
	data_copy[i]=0.;
    }
*/  


    // fill half of the array as power (actually energy) and another half (actually extanded part) with zero padding (Numerical Recipes 643 page)
    for (int i=0;i<length;i++) {
	power_noise_copy[i]=(data[i]*data[i])/Zr*TIMESTEP;
    }
    for (int i=length;i<length*2;i++) {
	power_noise_copy[i]=0.;
    }
    
    
    
    // do forward fft to get freq domain (energy of pure signal)
    Tools::realft(power_noise_copy,1,length*2);
    
    double ans_copy[length*2];
    
    
    
    // change the sign (from numerical recipes 648, 649 page)
    for (int j=1;j<length;j++) {
	ans_copy[2*j]=(power_noise_copy[2*j]*fdiode[2*j]-power_noise_copy[2*j+1]*fdiode[2*j+1])/((double)length);
	//ans_copy[2*j]=(power_noise_copy[2*j]*fdiode[2*j]+power_noise_copy[2*j+1]*fdiode[2*j+1])/((double)length/2);
	ans_copy[2*j+1]=(power_noise_copy[2*j+1]*fdiode[2*j]+power_noise_copy[2*j]*fdiode[2*j+1])/((double)length);
	//ans_copy[2*j+1]=(power_noise_copy[2*j+1]*fdiode[2*j]-power_noise_copy[2*j]*fdiode[2*j+1])/((double)length/2);
    }
    ans_copy[0]=power_noise_copy[0]*fdiode[0]/((double)length);
    ans_copy[1]=power_noise_copy[1]*fdiode[1]/((double)length);

    // 1/length is actually 2/(length * 2)
    //
    
    Tools::realft(ans_copy,-1,length*2);
    
    diodeconv.clear();  // remove previous values in diodeconv
    
    // only save first half of convolution result as last half is result from zero padding (and some spoiled bins)
    //for (int i=0;i<length;i++) {
    for (int i=0;i<length+maxt_diode_bin;i++) {
	//power_noise[i]=power_noise_copy[i];
	//diodeconv[i]=ans_copy[i];
	diodeconv.push_back( ans_copy[i] );
    }
    
    
    
    int iminsamp,imaxsamp; // find min and max samples such that
    // the diode response is fully overlapping with the noise waveform
    iminsamp=(int)(maxt_diode/TIMESTEP);
    // the noise waveform is NFOUR/2 long
    // then for a time maxt_diode/TIMESTEP at the end of that, the
    // diode response function is only partially overlappying with the
    // waveform in the convolution
    imaxsamp=NFOUR/2;
    
    //  cout << "iminsamp, imaxsamp are " << iminsamp << " " << imaxsamp << "\n";
    if (imaxsamp<iminsamp) {
	cout << "Noise waveform is not long enough for this diode response.\n";
	exit(1);
    }
    
    //int ibin=((iminsamp+imaxsamp)-(iminsamp+imaxsamp)%2)/2;
    //cout << "ibin is " << ibin << "\n";
    
    // return the 50th sample, right in the middle
   // onediodeconvl=diodeconv[ibin];
    
    
    //mindiodeconvl=0.;
    
    //for (int i=iminsamp;i<imaxsamp;i++) {
	
//	if (diodeconv[i]<mindiodeconvl)
//	    mindiodeconvl=diodeconv[i];
	
	
  //  }
    
}








// test for myconvlv with NFOUR/2 version
void Trigger::myconvlv_half(vector <double> &data,const int NFOUR,vector <double> &fdiode, vector <double> &diodeconv) {
    
    
    const int length=NFOUR;
//    double data_copy[length];
    //double fdiode_real[length];
    double power_noise_copy[length];

/*    
    for (int i=0;i<NFOUR/2;i++) {
	data_copy[i]=data[i];
    }
    for (int i=NFOUR/2;i<length;i++) {
	data_copy[i]=0.;
    }
*/  


    for (int i=0;i<length;i++) {
	power_noise_copy[i]=(data[i]*data[i])/Zr*TIMESTEP;
    }
    
    
    
    Tools::realft(power_noise_copy,1,length);
    
    double ans_copy[length];
    
    
    
    // change the sign (from numerical recipes 648, 649 page)
    for (int j=1;j<length/2;j++) {
	ans_copy[2*j]=(power_noise_copy[2*j]*fdiode[2*j]-power_noise_copy[2*j+1]*fdiode[2*j+1])/((double)length/2);
	//ans_copy[2*j]=(power_noise_copy[2*j]*fdiode[2*j]+power_noise_copy[2*j+1]*fdiode[2*j+1])/((double)length/2);
	ans_copy[2*j+1]=(power_noise_copy[2*j+1]*fdiode[2*j]+power_noise_copy[2*j]*fdiode[2*j+1])/((double)length/2);
	//ans_copy[2*j+1]=(power_noise_copy[2*j+1]*fdiode[2*j]-power_noise_copy[2*j]*fdiode[2*j+1])/((double)length/2);
    }
    ans_copy[0]=power_noise_copy[0]*fdiode[0]/((double)length/2);
    ans_copy[1]=power_noise_copy[1]*fdiode[1]/((double)length/2);
    
    
    Tools::realft(ans_copy,-1,length);
    
    diodeconv.clear();  // remove previous values in diodeconv
    
    for (int i=0;i<length;i++) {
	//power_noise[i]=power_noise_copy[i];
	//diodeconv[i]=ans_copy[i];
	diodeconv.push_back( ans_copy[i] );
    }
    
    
    
    int iminsamp,imaxsamp; // find min and max samples such that
    // the diode response is fully overlapping with the noise waveform
    iminsamp=(int)(maxt_diode/TIMESTEP);
    // the noise waveform is NFOUR/2 long
    // then for a time maxt_diode/TIMESTEP at the end of that, the
    // diode response function is only partially overlappying with the
    // waveform in the convolution
    imaxsamp=NFOUR/2;
    
    //  cout << "iminsamp, imaxsamp are " << iminsamp << " " << imaxsamp << "\n";
    if (imaxsamp<iminsamp) {
	cout << "Noise waveform is not long enough for this diode response.\n";
	exit(1);
    }
    
    //int ibin=((iminsamp+imaxsamp)-(iminsamp+imaxsamp)%2)/2;
    //cout << "ibin is " << ibin << "\n";
    
    // return the 50th sample, right in the middle
   // onediodeconvl=diodeconv[ibin];
    
    
    //mindiodeconvl=0.;
    
    //for (int i=iminsamp;i<imaxsamp;i++) {
	
//	if (diodeconv[i]<mindiodeconvl)
//	    mindiodeconvl=diodeconv[i];
	
	
  //  }
    
}



void Trigger::myconvlv_half(double *data,const int NFOUR,vector <double> &fdiode, vector <double> &diodeconv) {
    
    
    const int length=NFOUR;
//    double data_copy[length];
    //double fdiode_real[length];
    double power_noise_copy[length];
  
/*    
    for (int i=0;i<NFOUR/2;i++) {
	data_copy[i]=data[i];
    }
    for (int i=NFOUR/2;i<length;i++) {
	data_copy[i]=0.;
    }
*/  


    for (int i=0;i<length;i++) {
	power_noise_copy[i]=(data[i]*data[i])/Zr*TIMESTEP;
    }
    
    
    
    Tools::realft(power_noise_copy,1,length);
    
    double ans_copy[length];
    
    
    
    for (int j=1;j<length/2;j++) {
	ans_copy[2*j]=(power_noise_copy[2*j]*fdiode[2*j]-power_noise_copy[2*j+1]*fdiode[2*j+1])/((double)length/2);
	//ans_copy[2*j]=(power_noise_copy[2*j]*fdiode[2*j]+power_noise_copy[2*j+1]*fdiode[2*j+1])/((double)length/2);
	ans_copy[2*j+1]=(power_noise_copy[2*j+1]*fdiode[2*j]+power_noise_copy[2*j]*fdiode[2*j+1])/((double)length/2);
	//ans_copy[2*j+1]=(power_noise_copy[2*j+1]*fdiode[2*j]-power_noise_copy[2*j]*fdiode[2*j+1])/((double)length/2);
    }
    ans_copy[0]=power_noise_copy[0]*fdiode[0]/((double)length/2);
    ans_copy[1]=power_noise_copy[1]*fdiode[1]/((double)length/2);
    
    
    Tools::realft(ans_copy,-1,length);
    
    diodeconv.clear();  // remove previous values in diodeconv
    
    for (int i=0;i<length;i++) {
	//power_noise[i]=power_noise_copy[i];
	//diodeconv[i]=ans_copy[i];
	diodeconv.push_back( ans_copy[i] );
    }
    
    
    
    int iminsamp,imaxsamp; // find min and max samples such that
    // the diode response is fully overlapping with the noise waveform
    iminsamp=(int)(maxt_diode/TIMESTEP);
    // the noise waveform is NFOUR/2 long
    // then for a time maxt_diode/TIMESTEP at the end of that, the
    // diode response function is only partially overlappying with the
    // waveform in the convolution
    imaxsamp=NFOUR/2;
    
    //  cout << "iminsamp, imaxsamp are " << iminsamp << " " << imaxsamp << "\n";
    if (imaxsamp<iminsamp) {
	cout << "Noise waveform is not long enough for this diode response.\n";
	exit(1);
    }
    
    //int ibin=((iminsamp+imaxsamp)-(iminsamp+imaxsamp)%2)/2;
    //cout << "ibin is " << ibin << "\n";
    
    // return the 50th sample, right in the middle
   // onediodeconvl=diodeconv[ibin];
    
    
    //mindiodeconvl=0.;
    
    //for (int i=iminsamp;i<imaxsamp;i++) {
	
//	if (diodeconv[i]<mindiodeconvl)
//	    mindiodeconvl=diodeconv[i];
	
	
  //  }
    
}





int Trigger::CheckChannelsPass( vector <double> &V_total_diode ) {

    int ch_pass_bin = 0;

    //for (int ibin=iminbin; ibin<imaxbin; ibin++) {
    for (int ibin=(int)(maxt_diode/TIMESTEP);ibin<DATA_BIN_SIZE;ibin++) {

        if ( V_total_diode[ibin] < powerthreshold * rmsdiode ) {
            ch_pass_bin = ibin;
            ibin = DATA_BIN_SIZE;   // stop loop
        }

    }

    return ch_pass_bin; // if there was any iwindow values pass the threshold, return = bin number where passed the trigger, if nothing, return = 0;

    
}





//Trigger::~Trigger() {

//}
  

