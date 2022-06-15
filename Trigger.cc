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
    //
    // same with icemc -> anita -> Initialize

    // if using default noise temp setting (same temp for all chs)
  if (settings1->NOISE_CHANNEL_MODE == 0) {
    
    int ngeneratedevents=settings1->NOISE_EVENTS;  // should this value read at Settings class
    double v_noise[settings1->DATA_BIN_SIZE];    // noise voltage time domain (with filter applied)
    
    meandiode = 0.;
    rmsdiode = 0.;
    
    rmsvoltage = 0.;
    
    v_noise_timedomain.resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)
    v_noise_timedomain_diode.resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)
    
    cerr<<"Generating noise waveforms"<<endl;
    int print_steps;
    if ( ngeneratedevents > 10 ) print_steps = ngeneratedevents/10;
    else print_steps = ngeneratedevents;
    

    for (int i=0; i<ngeneratedevents; i++) {

        if ((i)%print_steps == 0) cerr<< (i/print_steps) * 10 <<"% done"<<endl;

        // get v_noise array (noise voltage in time domain)
        report->GetNoiseWaveforms(settings1, detector, V_noise_freqbin, v_noise);

        // do normal time ordering (not sure if this is necessary)
        Tools::NormalTimeOrdering(settings1->DATA_BIN_SIZE, v_noise);


	myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real_databin,v_noise_timedomain_diode[i]);
	//myconvlv_new(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real,v_noise_timedomain_diode[i]);

	for (int m=0; m<settings1->DATA_BIN_SIZE; m++) {
	    if ( m>=(int)(maxt_diode/TIMESTEP) && m<settings1->DATA_BIN_SIZE ) {
	    //bwslice_meandiode[j]+=timedomain_output_e[j][m]/((double)ngeneratedevents*((double)NFOUR/2-maxt_diode/TIMESTEP));
	    meandiode+=v_noise_timedomain_diode[i][m]/((double)ngeneratedevents * ((double)settings1->DATA_BIN_SIZE-maxt_diode/TIMESTEP));
	    } // end loop over samples where diode function is fully contained 

            v_noise_timedomain[i].push_back(v_noise[m]);    // save pure noise (not diode convlved) waveforms

            // save pure noise spectrum before applying Rayleigh distribution (only at first evt as they will be all same)
            if ( i == 0 && m<settings1->DATA_BIN_SIZE/2 ) {
                Vfft_noise_before.push_back(report->Vfft_noise_before[m]);    // save pure noise (not diode convlved) waveforms
            }

        }


    }   // get meandiode with 1000 noisewaveforms
    cerr<< "100% done"<<endl;
    





    // now as v_noise_timedomain_diode's waveforms are still stored in it, we can just calculate rms from them
    //
    cerr<<"Calculating noise waveforms RMS"<<endl;
    for (int i=0; i<ngeneratedevents; i++) {

        if ((i)%print_steps == 0) cerr<< (i/print_steps) * 10 <<"% done"<<endl;

        // we are going to reuse generated v_noise_timedomain_diode above.
        //report->GetNoiseWaveforms(settings1, detector, V_noise_freqbin, v_noise);
	//myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real,v_noise_timedomain_diode);

        for (int m=(int)(maxt_diode/TIMESTEP);m<settings1->DATA_BIN_SIZE;m++) {
		
	    rmsdiode+=(v_noise_timedomain_diode[i][m]-meandiode)*(v_noise_timedomain_diode[i][m]-meandiode)/((double)ngeneratedevents * ((double)settings1->DATA_BIN_SIZE-maxt_diode/TIMESTEP));

	    rmsvoltage+=(v_noise_timedomain[i][m])*(v_noise_timedomain[i][m])/((double)ngeneratedevents * ((double)settings1->DATA_BIN_SIZE-maxt_diode/TIMESTEP));
		
	}
    }   // get rmsdiode with 1000 noisewaveforms
    cerr<< "100% done"<<endl;


    // now we can use stored v_noise_timedomain_diode noise waveforms anytime later.


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

    meandiode_ch.resize(num_chs);
    rmsdiode_ch.resize(num_chs);
    rmsvoltage_ch.resize(num_chs);


    v_noise_timedomain_ch.resize(num_chs);  // make the size of v_noise_timedomain_diode as number of chs
    v_noise_timedomain_diode_ch.resize(num_chs);  // make the size of v_noise_timedomain_diode as number of chs


    for (int i=0; i<num_chs; i++) {
        v_noise_timedomain_ch[i].resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)
        v_noise_timedomain_diode_ch[i].resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)

        meandiode_ch[i] = 0.;
        rmsdiode_ch[i] = 0.;
        rmsvoltage_ch[i] = 0.;

    }


    cerr<<"Generating noise waveforms"<<endl;
    //int print_steps = ngeneratedevents/10;
    int print_steps;
    if ( ngeneratedevents > 10 ) print_steps = ngeneratedevents/10;
    else print_steps = ngeneratedevents;


for (int ch=0; ch<num_chs; ch++) {

    for (int i=0; i<ngeneratedevents; i++) {

        if ((i)%print_steps == 0) cerr<< (i/print_steps) * 10 <<"% done for ch"<<ch<<endl;

        // get v_noise array (noise voltage in time domain)
        report->GetNoiseWaveforms_ch(settings1, detector, V_noise_freqbin_ch[ch], v_noise, ch);

        // do normal time ordering (not sure if this is necessary)
        Tools::NormalTimeOrdering(settings1->DATA_BIN_SIZE, v_noise);


	myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real_databin,v_noise_timedomain_diode_ch[ch][i]);
        //

	for (int m=0; m<settings1->DATA_BIN_SIZE; m++) {
	    if ( m>=(int)(maxt_diode/TIMESTEP) && m<settings1->DATA_BIN_SIZE ) {
	    //bwslice_meandiode[j]+=timedomain_output_e[j][m]/((double)ngeneratedevents*((double)NFOUR/2-maxt_diode/TIMESTEP));
	    meandiode_ch[ch] += v_noise_timedomain_diode_ch[ch][i][m]/((double)ngeneratedevents * ((double)settings1->DATA_BIN_SIZE-maxt_diode/TIMESTEP));

	    } // end loop over samples where diode function is fully contained 

            v_noise_timedomain_ch[ch][i].push_back(v_noise[m]);    // save pure noise (not diode convlved) waveforms

            // save pure noise spectrum before applying Rayleigh distribution (only at first evt as they will be all same)
            if ( i == 0 && m<settings1->DATA_BIN_SIZE/2 ) {
                Vfft_noise_before.push_back(report->Vfft_noise_before[m]);    // save pure noise (not diode convlved) waveforms
            }

        }


    }   // get meandiode with 1000 noisewaveforms

} // loop over chs

    cerr<< "100% done"<<endl;


    cerr<<"Calculating noise waveforms RMS"<<endl;

for (int ch=0; ch<num_chs; ch++) {

    for (int i=0; i<ngeneratedevents; i++) {

        if ((i)%print_steps == 0) cerr<< (i/print_steps) * 10 <<"% done for ch"<<ch<<endl;


        for (int m=(int)(maxt_diode/TIMESTEP);m<settings1->DATA_BIN_SIZE;m++) {
		
	    rmsdiode_ch[ch]+=(v_noise_timedomain_diode_ch[ch][i][m]-meandiode_ch[ch])*(v_noise_timedomain_diode_ch[ch][i][m]-meandiode_ch[ch])/((double)ngeneratedevents * ((double)settings1->DATA_BIN_SIZE-maxt_diode/TIMESTEP));

	    rmsvoltage_ch[ch]+=(v_noise_timedomain_ch[ch][i][m])*(v_noise_timedomain_ch[ch][i][m])/((double)ngeneratedevents * ((double)settings1->DATA_BIN_SIZE-maxt_diode/TIMESTEP));
		
	}
    }   // get rmsdiode with 1000 noisewaveforms
}
    cerr<< "100% done"<<endl;

    


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

    meandiode_ch.resize(num_chs);
    rmsdiode_ch.resize(num_chs);
    rmsvoltage_ch.resize(num_chs);



    v_noise_timedomain_ch.resize(num_chs);  // make the size of v_noise_timedomain_diode as number of chs
    v_noise_timedomain_diode_ch.resize(num_chs);  // make the size of v_noise_timedomain_diode as number of chs


    for (int i=0; i<num_chs; i++) {
        v_noise_timedomain_ch[i].resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)
        v_noise_timedomain_diode_ch[i].resize(ngeneratedevents);  // make the size of v_noise_timedomain_diode as ngeneratedevents (this will be huge!)

        meandiode_ch[i] = 0.;
        rmsdiode_ch[i] = 0.;
        rmsvoltage_ch[i] = 0.;

    }


    cerr<<"Generating noise waveforms"<<endl;
    //int print_steps = ngeneratedevents/10;
    int print_steps;
    if ( ngeneratedevents > 10 ) print_steps = ngeneratedevents/10;
    else print_steps = ngeneratedevents;


for (int ch=0; ch<num_chs; ch++) {

    for (int i=0; i<ngeneratedevents; i++) {

        if ((i)%print_steps == 0) cerr<< (i/print_steps) * 10 <<"% done for ch"<<ch<<endl;

        // get v_noise array (noise voltage in time domain)
        report->GetNoiseWaveforms_ch(settings1, detector, V_noise_freqbin_ch[ch], v_noise, ch);

        // do normal time ordering (not sure if this is necessary)
        Tools::NormalTimeOrdering(settings1->DATA_BIN_SIZE, v_noise);


	myconvlv(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real_databin,v_noise_timedomain_diode_ch[ch][i]);
        //

	for (int m=0; m<settings1->DATA_BIN_SIZE; m++) {
	    if ( m>=(int)(maxt_diode/TIMESTEP) && m<settings1->DATA_BIN_SIZE ) {
	    //bwslice_meandiode[j]+=timedomain_output_e[j][m]/((double)ngeneratedevents*((double)NFOUR/2-maxt_diode/TIMESTEP));
	    meandiode_ch[ch] += v_noise_timedomain_diode_ch[ch][i][m]/((double)ngeneratedevents * ((double)settings1->DATA_BIN_SIZE-maxt_diode/TIMESTEP));

	    } // end loop over samples where diode function is fully contained 

            v_noise_timedomain_ch[ch][i].push_back(v_noise[m]);    // save pure noise (not diode convlved) waveforms

            // save pure noise spectrum before applying Rayleigh distribution (only at first evt as they will be all same)
            if ( i == 0 && m<settings1->DATA_BIN_SIZE/2 ) {
                Vfft_noise_before.push_back(report->Vfft_noise_before[m]);    // save pure noise (not diode convlved) waveforms
            }

        }


    }   // get meandiode with 1000 noisewaveforms

} // loop over chs

    cerr<< "100% done"<<endl;


    cerr<<"Calculating noise waveforms RMS"<<endl;

for (int ch=0; ch<num_chs; ch++) {

    for (int i=0; i<ngeneratedevents; i++) {

        if ((i)%print_steps == 0) cerr<< (i/print_steps) * 10 <<"% done for ch"<<ch<<endl;


        for (int m=(int)(maxt_diode/TIMESTEP);m<settings1->DATA_BIN_SIZE;m++) {
		
	    rmsdiode_ch[ch]+=(v_noise_timedomain_diode_ch[ch][i][m]-meandiode_ch[ch])*(v_noise_timedomain_diode_ch[ch][i][m]-meandiode_ch[ch])/((double)ngeneratedevents * ((double)settings1->DATA_BIN_SIZE-maxt_diode/TIMESTEP));

	    rmsvoltage_ch[ch]+=(v_noise_timedomain_ch[ch][i][m])*(v_noise_timedomain_ch[ch][i][m])/((double)ngeneratedevents * ((double)settings1->DATA_BIN_SIZE-maxt_diode/TIMESTEP));
		
	}
    }   // get rmsdiode with 1000 noisewaveforms
}
    cerr<< "100% done"<<endl;

    


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
                //myconvlv_new(v_noise, settings1->DATA_BIN_SIZE, detector->fdiode_real,v_noise_timedomain_diode[i]);

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
  

