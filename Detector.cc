#include "Detector.h"
#include "Tools.h"
#include "Event.h"
#include "IceModel.h"
#include "Settings.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "Constants.h"
#include "TF1.h"

#ifdef ARA_UTIL_EXISTS
#include "AraGeomTool.h"
#endif

ClassImp(Detector);
ClassImp(Parameters);
ClassImp(InstalledStation);
ClassImp(Surface_antenna);
ClassImp(Antenna);
ClassImp(Antenna_string);
ClassImp(ARA_station);


Detector::Detector() {
    //Default constructor
}


Detector::Detector(Settings *settings1, IceModel *icesurface, string setupfile) {
    //Detector::Detector(int mode, IceModel *icesurface) {
    cout << "good" << endl;
    
    // set freq_forfft for later use
    //
    
    // set freq_forfft array
    // same with icemc anita class initialization function
    
    
    double freqstep=1./(double)(settings1->NFOUR/2)/(settings1->TIMESTEP);
    
    NFOUR = settings1->NFOUR;
    TIMESTEP = settings1->TIMESTEP;
    
    NoiseFig_numCh=16;

    //for (int i=0;i<HALFNFOUR/2;i++) {
    for (int i=0;i<settings1->NFOUR/4;i++) {
        //--------------------------------------------------
        // freq_forfft[2*i]=(double)i*freqstep;
        // freq_forfft[2*i+1]=(double)i*freqstep;
        //-------------------------------------------------- 
        freq_forfft.push_back( (double)i*freqstep );    // even numbers
        freq_forfft.push_back( (double)i*freqstep );    // odd numbers
        
    }
    for (int i=settings1->NFOUR/4;i<settings1->NFOUR/2;i++) {
        //--------------------------------------------------
        // freq_forfft[2*i]=(double)i*freqstep;
        // freq_forfft[2*i+1]=(double)i*freqstep;
        //-------------------------------------------------- 
        freq_forfft.push_back( (double)i*freqstep );    // even numbers
        freq_forfft.push_back( (double)i*freqstep );    // odd numbers
        
    }
    // end of settings freq_forfft
    
    
    
    
    //set mode ex) mode 0 = testbed,
    // mode 1 = ARA_1
    // mode 2 = ARA_2
    // mode 3 = ARA_37
    int mode = settings1->DETECTOR;
    Detector_mode = mode;
    
    int string_id = -1;
    //    int string_id = 0;
    int antenna_id = 0;
    
    
    //    Parameters params;
    //    vector <Antenna_string> strings;
    ARA_station temp_station;
    Antenna_string temp;
    Antenna_string temp_string;
    Antenna temp_antenna;
    Surface_antenna temp_surface;
    
    params.number_of_strings = 0;
    params.number_of_antennas = 0;
    
    //initialize few params values.
    params.freq_step = 60;
    params.ang_step = 2664;
    params.freq_width = 16.667;
    params.freq_init = 83.333;
    params.DeployedStations = 4;
    //end initialize
    
    //Parameters to use if using Arianna_WIPLD_hpol.dat
    
    if (settings1->ANTENNA_MODE == 2){
      params.freq_step = 238;//60
      params.ang_step = 2664;//2664;
      params.freq_width = 5; //16.667;
      params.freq_init = 15;//83.333;
      params.DeployedStations = 4;
    }
    
    //copy freq_width, freq_init in params to Detector freq_width, freq_init
    freq_step = params.freq_step;
    ang_step = params.ang_step;
    freq_width = params.freq_width;
    freq_init = params.freq_init;
    //end copy
    
    
    string testbed_file = "testbed_info.txt";
//--------------------------------------------------
//     string ARA_N_file = "ARA_N_info.txt";
//     string ARA37_file = "ARA37_info.txt";
//-------------------------------------------------- 
    string ARA_N_file = setupfile;
    string ARA37_file = setupfile;
    
    string line, label;

    // setup installed station information
    // setup actual installed staion information regardless of what DETECTOR mode is in use
    SetupInstalledStations();
    
    //    IceModel *icesurface = new IceModel;
    //cout<<"Ice surface at 0,0 : "<<icesurface->Geoid(0.)<<endl;
    
    
    ////////////////////////////////////////////////////////////////////////////////////    
    
    if (mode == 0) {
        cout<<"\n\tDector mode 0 : testbed"<<endl;
        ifstream testbed( testbed_file.c_str() );
        cout<<"We use "<<testbed_file.c_str()<<" as antenna info."<<endl;
        
        
        if ( testbed.is_open() ) {
            while (testbed.good() ) {
                getline (testbed, line);
                
                if (line[0] != "/"[0]) {
                    label = line.substr(0, line.find_first_of("=") );
                    
                    if (label == "number_of_strings") {
                        params.number_of_strings = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                        for (int i=0; i<(int) params.number_of_strings; i++) {
                            strings.push_back(temp);
                        }
                        cout<<"read numner_of_strings"<<endl;
                        //                        Parameters.number_of_strings = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                    }
                    else if (label == "antenna_string") {
                        string_id++;
                        antenna_id = 0;
                        cout<<"read antenna_string"<<endl;
                        //                        if (string_id + 1 >= params.number_of_strings) {
                        //                            cout<<"Error! Too many strings!"<<endl;
                        //                        }
                    }
                    else if (label == "x") {
                        //strings[string_id].x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        strings[string_id].SetX( atof( line.substr( line.find_first_of("=") + 1).c_str() ) );
                        cout<<"read x : "<<(double)strings[string_id].GetX()<<" string_id : "<<string_id<<endl;
                    }
                    else if (label == "y") {
                        //strings[string_id].y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        strings[string_id].SetY( atof( line.substr( line.find_first_of("=") + 1).c_str() ) );
                        cout<<"read y : "<<(double)strings[string_id].GetY()<<" string_id : "<<string_id<<endl;
                    }
                    else if (label == "z") {
                        strings[string_id].antennas.push_back(temp_antenna);
                        //strings[string_id].antennas[antenna_id].z = atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() );
                        strings[string_id].antennas[antenna_id].SetZ( atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() ) );
                        strings[string_id].antennas[antenna_id].type = atoi( line.substr( line.find_first_of(",") + 1).c_str() );
                        cout<<"read z : "<<(double)strings[string_id].antennas[antenna_id].GetZ()<<" string_id : "<<string_id<<" antenna_id : "<<antenna_id<<" type : "<<(int)strings[string_id].antennas[antenna_id].type<<endl;
                        antenna_id++;
                        params.number_of_antennas++;
                        //                        Parameters.number_of_antennas++;
                    }
                    
                    
                }
            }
            testbed.close();
        }
        
        else {
            cout<<"Unable to open antenna array file !"<<endl;
            //            return 1;
        }
        
        
        // testbed version of FlattoEarth_ARA 
        // strings and antennas on the strings use geoid surface!
        double Dist = 0.;   //for sqrt(x^2 + y^2)
        double R1 = icesurface->Surface(0.,0.); // from core of earth to surface at theta, phi = 0.
        double theta_tmp;
        double phi_tmp;
        
        // set same theta, phi to all antennas in same string
        for (int i=0; i<params.number_of_strings; i++) {
            
            Dist = sqrt( pow(strings[i].GetX(),2) + pow(strings[i].GetY(),2) );
            theta_tmp = Dist/R1;    // assume R1 is constant (which is not)
            phi_tmp = atan2(strings[i].GetY(),strings[i].GetX());
            
            if (phi_tmp<0.) phi_tmp += 2.*PI;
            
            // set theta, phi for strings.
            strings[i].SetThetaPhi(theta_tmp, phi_tmp);
            //set R for strings.
            strings[i].SetR( icesurface->Surface( strings[i].Lon(), strings[i].Lat()) );
            
            cout<<"R, Theta, Phi : "<<strings[i].R()<<" "<<strings[i].Theta()<<" "<<strings[i].Phi()<<endl;
            
            // set antennas r, theta, phi
            for (int j=0; j<antenna_id; j++) {
                strings[i].antennas[j].SetRThetaPhi( strings[i].R() + strings[i].antennas[j].GetZ() , strings[i].Theta(), strings[i].Phi() );
            }
        }
        
        
        
        
        
        
    }
    
    /////////////////////////////////////////////////////////////////////////////
    
    
    
    
    else if (mode == 1) {
//        cout<<"\n\tDector mode 1 : Specific number of stations (less than 7 stations) !"<<endl;
        ifstream ARA_N( ARA_N_file.c_str() );
//        cout<<"We use "<<ARA_N_file.c_str()<<" as antenna info."<<endl;
        
        
        // initialize info
        params.number_of_stations = 1;
        params.number_of_strings_station = 4;   // ARA-1 has 4 strings
        params.number_of_antennas_string = 4; // 4 antennas on each strings
        params.number_of_surfaces_station = 4;
        
	

        //double core_x = 0.; 
        //double core_y = 0.;
        params.core_x = 10000.; 
        params.core_y = 10000.;
        double R_string = 10.;  // all units are in meter
        double R_surface = 60.;
        double z_max = 200.;
        double z_btw = 10.;
        double z_btw_array[6]; // assume there will be less than 6 bore hole antennas at each string
	// these z_btw array will be used when settings->BH_ANT_SEP_DIST_ON=1 case
        for (int i=0; i<6; i++) {
            if (i==0) z_btw_array[i] = 0.;
            //else z_btw_array[i] = z_btw;
            else if (i==1) z_btw_array[i] = 2.;
            else if (i==2) z_btw_array[i] = 15.;
            else if (i==3) z_btw_array[i] = 2.;
            else z_btw_array[i] = z_btw;
        }
        double z_btw_total;
        params.stations_per_side = 4;       // total 37 stations
        params.station_spacing = 2000.;     // 2km spacing
        params.antenna_orientation = 0;     // all antenna facing x
        params.bore_hole_antenna_layout = settings1->BORE_HOLE_ANTENNA_LAYOUT;
        // finish initialization
        //
        
        
        
        
        
        // Read new parameters if there are...
        if ( ARA_N.is_open() ) {
            while (ARA_N.good() ) {
                getline (ARA_N, line);
                
                if (line[0] != "/"[0]) {
                    label = line.substr(0, line.find_first_of("=") );
                    
                    if (label == "core_x") {
                        params.core_x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_x"<<endl;
                    }
                    else if (label == "core_y") {
                        params.core_y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_y"<<endl;
                    }
                    else if (label == "R_string") {
                        R_string = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_string"<<endl;
                    }
                    else if (label == "R_surface") {
                        R_surface = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_surface"<<endl;
                    }
                    else if (label == "z_max") {
                        z_max = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_max"<<endl;
                    }
                    else if (label == "z_btw") {
                        z_btw = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw"<<endl;
                    }
                    else if (label == "z_btw01") {
                        z_btw_array[1] = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw bh ant0 and ant1"<<endl;
                    }
                    else if (label == "z_btw12") {
                        z_btw_array[2] = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw bh ant1 and ant2"<<endl;
                    }
                    else if (label == "z_btw23") {
                        z_btw_array[3] = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw bh ant2 and ant3"<<endl;
                    }
                    else if (label == "z_btw34") {
                        z_btw_array[4] = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw bh ant3 and ant4"<<endl;
                    }
                    else if (label == "z_btw45") {
                        z_btw_array[5] = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw bh ant4 and ant5"<<endl;
                    }
                    else if (label == "number_of_stations") {
                        params.number_of_stations = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read stations_per_side"<<endl;
                    }
                    else if (label == "station_spacing") {
                        params.station_spacing = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read station_spacting"<<endl;
                    }
                    else if (label == "antenna_orientation") {
                        params.antenna_orientation = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read antenna_orientation"<<endl;
                    }
                    else if (label == "number_of_strings_station") {
                        params.number_of_strings_station = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read number of strings"<<endl;
                    }
                    else if (label == "number_of_antennas_string") {
                        params.number_of_antennas_string = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read number of antennas per string"<<endl;
                    }
                }
            }
            ARA_N.close();
        }
        // finished reading new parameters
        
        
        
        // set number of antennas in a string
        if (params.bore_hole_antenna_layout == 0) { // VHVH layout
            params.number_of_antennas_string = 4;
        }
        else if (params.bore_hole_antenna_layout == 1) { // VHV layout
            params.number_of_antennas_string = 3;
        }
        else if (params.bore_hole_antenna_layout == 2) { // VHVV layout
            params.number_of_antennas_string = 4;
        }
        else if (params.bore_hole_antenna_layout == 3) { // VHHH layout
            params.number_of_antennas_string = 4;
        }
        else if (params.bore_hole_antenna_layout == 4) { // VHH layout
            params.number_of_antennas_string = 3;
        }
        else if (params.bore_hole_antenna_layout == 5) { // VVVV layout
            params.number_of_antennas_string = 4;
        }
        else if (params.bore_hole_antenna_layout == 6) { // VV layout
            params.number_of_antennas_string = 2;
        }
	else if (params.bore_hole_antenna_layout == 7) { // V layout
            params.number_of_antennas_string = 1;
        }

        
        
        
        
        //
        // caculate number of stations, strings, antennas 
        params.number_of_strings = params.number_of_stations * params.number_of_strings_station;
        params.number_of_antennas = params.number_of_strings * params.number_of_antennas_string;
        //
        //
        
        
        
        //
        // prepare vectors
        for (int i=0; i<params.number_of_stations; i++) {
            stations.push_back(temp_station);
            
            for (int j=0; j<params.number_of_surfaces_station; j++) {
                stations[i].surfaces.push_back(temp_surface);
            }
            
            for (int k=0; k<params.number_of_strings_station; k++) {
                stations[i].strings.push_back(temp_string);
                
                for (int l=0; l<params.number_of_antennas_string; l++) {
                    stations[i].strings[k].antennas.push_back(temp_antenna);
                }
                
            }
            
            
        }
        // end prepare vectors
        //
        
        
        
        
        
        //
        // for ARA-37 (or more than 1 station case), need code for setting position for all 37 stations here!
        //
        int station_count = 0;
        
        int side_step;

        double next_dir = PI*2./3;

        for (int istation = 0; istation < (int)params.number_of_stations; istation++) {
	  stations[istation].StationID=istation;
	
            if (station_count < (int)params.number_of_stations - 1) {


                if ( station_count < 6 ) { // first layer
                    stations[station_count].SetX( params.core_x + (double)params.station_spacing * cos( (PI/3.) * (double)station_count ) );
                    stations[station_count].SetY( params.core_y + (double)params.station_spacing * sin( (PI/3.) * (double)station_count ) );
                }
                else if ( station_count < 18 ) { // second layer

                    // if the first outter layer station
                    if ( station_count == 6 ) {
                        side_step = 2;
                        stations[station_count].SetX( params.core_x + (double)params.station_spacing * 2. );
                        stations[station_count].SetY( params.core_y );
                    }
                    // after first station
                    else { 
                        if ( side_step > 0 ) {
                            stations[station_count].SetX( stations[station_count-1].GetX() + (double)params.station_spacing * cos(next_dir) );
                            stations[station_count].SetY( stations[station_count-1].GetY() + (double)params.station_spacing * sin(next_dir) );
                            side_step--;
                        }
                        else {
                            side_step = 1;
                            next_dir+=PI/3.; // rotate
                            stations[station_count].SetX( stations[station_count-1].GetX() + (double)params.station_spacing * cos(next_dir) );
                            stations[station_count].SetY( stations[station_count-1].GetY() + (double)params.station_spacing * sin(next_dir) );
                        }
                    }
                }
                else if ( station_count < 36 ) { // third layer

                    // if the first outter layer station
                    if ( station_count == 6 ) {
                        side_step = 3;
                        stations[station_count].SetX( params.core_x + (double)params.station_spacing * 3. );
                        stations[station_count].SetY( params.core_y );
                    }
                    // after first station
                    else { 
                        if ( side_step > 0 ) {
                            stations[station_count].SetX( stations[station_count-1].GetX() + (double)params.station_spacing * cos(next_dir) );
                            stations[station_count].SetY( stations[station_count-1].GetY() + (double)params.station_spacing * sin(next_dir) );
                            side_step--;
                        }
                        else {
                            side_step = 2;
                            next_dir+=PI/3.; // rotate
                            stations[station_count].SetX( stations[station_count-1].GetX() + (double)params.station_spacing * cos(next_dir) );
                            stations[station_count].SetY( stations[station_count-1].GetY() + (double)params.station_spacing * sin(next_dir) );
                        }
                    }

                }

                station_count++;
            }
            else if (station_count < (int)params.number_of_stations) {
                //stations[station_count].x = core_x;
                //stations[station_count].y = core_y;
                stations[station_count].SetX( params.core_x );
                stations[station_count].SetY( params.core_y );
                station_count++;
            }
            else {
                cout<<"\n\tError, too many stations !"<<endl;
            }
        }
        // finished setting all stations' position
        

        
//        cout<<"total station_count : "<<station_count<<endl;
        if (station_count != (int)params.number_of_stations) cout<<"\n\tError, station number not match !"<<endl;       
        
        //
        // set antenna values from parameters
        // set station positions
        if (settings1->READGEOM == 0) { // use idealized geometry
            //SetupInstalledStations();

            for (int i=0; i<params.number_of_stations; i++) {
                
                //
                // set string postions based on station position
                stations[i].strings[0].SetX( stations[i].GetX() - (R_string * cos(PI/4.)) );
                stations[i].strings[0].SetY( stations[i].GetY() + (R_string * sin(PI/4.)) );
                
                stations[i].strings[1].SetX( stations[i].GetX() + (R_string * cos(PI/4.)) );
                stations[i].strings[1].SetY( stations[i].GetY() + (R_string * sin(PI/4.)) );
                
                stations[i].strings[2].SetX( stations[i].GetX() - (R_string * cos(PI/4.)) );
                stations[i].strings[2].SetY( stations[i].GetY() - (R_string * sin(PI/4.)) );
                
                stations[i].strings[3].SetX( stations[i].GetX() + (R_string * cos(PI/4.)) );
                stations[i].strings[3].SetY( stations[i].GetY() - (R_string * sin(PI/4.)) );
                
                
                //
                // set antenna postions in borehole
                // and set type (h or v pol antenna) and set orientation (facing x or y)
                if ( params.bore_hole_antenna_layout == 0 || params.bore_hole_antenna_layout == 1) {
                    
                    for (int j=0; j<params.number_of_strings_station; j++) {
                        for (int k=0; k<params.number_of_antennas_string; k++) {

                            if (settings1->BH_ANT_SEP_DIST_ON==0) 
                            stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw*k );

                            else if (settings1->BH_ANT_SEP_DIST_ON==1) {
                                z_btw_total = 0.;
                                for (int l=0; l<k+1; l++) {
                                    z_btw_total += z_btw_array[l];
                                }
                                stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw_total );
                            }
                            
                            if (k%2 == 0) {
                                stations[i].strings[j].antennas[k].type = 0;   // v-pol
                            }
                            else {
                                stations[i].strings[j].antennas[k].type = 1;   // h-pol
                            }
                            
                            if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                                stations[i].strings[j].antennas[k].orient = 0;
                            }
                            else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                                if ( j==0||j==3 ) {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                }
                                else {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                }
                                
                            }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                        }
                    }
                    
                } // end if bore hole antenna layout = 0 or 1 (where VHVH way but different numbers)
                
                
                else if ( params.bore_hole_antenna_layout == 2) {   // it's V-H-V-V
                    
                    for (int j=0; j<params.number_of_strings_station; j++) {
                        for (int k=0; k<params.number_of_antennas_string; k++) {

                            if (settings1->BH_ANT_SEP_DIST_ON==0) 
                            stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw*k );

                            else if (settings1->BH_ANT_SEP_DIST_ON==1) {
                                z_btw_total = 0.;
                                for (int l=0; l<k+1; l++) {
                                    z_btw_total += z_btw_array[l];
                                }
                                stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw_total );
                            }
                            
                            if (k == 1) {   // only the second antenna is H pol
                                stations[i].strings[j].antennas[k].type = 1;   // h-pol
                            }
                            else {  // other antennas are V pol
                                stations[i].strings[j].antennas[k].type = 0;   // v-pol
                            }
                            
                            if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                                stations[i].strings[j].antennas[k].orient = 0;
                            }
                            else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                                if ( j==0||j==3 ) {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                }
                                else {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                }
                                
                            }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                        }
                    }
                    
                } // end if bore hole antenna layout = 2 (where VHVV way but different numbers)
                
                
                else if ( params.bore_hole_antenna_layout == 3 || params.            
                         bore_hole_antenna_layout == 4 ) {   // it's V-H-H-H or V-H-H
                    
                    for (int j=0; j<params.number_of_strings_station; j++) {
                        for (int k=0; k<params.number_of_antennas_string; k++) {

                            if (settings1->BH_ANT_SEP_DIST_ON==0) 
                            stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw*k );

                            else if (settings1->BH_ANT_SEP_DIST_ON==1) {
                                z_btw_total = 0.;
                                for (int l=0; l<k+1; l++) {
                                    z_btw_total += z_btw_array[l];
                                }
                                stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw_total );
                            }

                            
                            if (k == 0) {   // only the first antenna is V pol
                                stations[i].strings[j].antennas[k].type = 0;   // v-pol
                            }
                            else {  // other antennas are H pol
                                stations[i].strings[j].antennas[k].type = 1;   // h-pol
                            }
                            
                            if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                                stations[i].strings[j].antennas[k].orient = 0;
                            }
                            else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                                if ( j==0||j==3 ) {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                }
                                else {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                }
                                
                            }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                        }
                    }
                    
                } // end if bore hole antenna layout = 3 (where VHHH way)

                else if ( params.bore_hole_antenna_layout == 5 || 
			  params.bore_hole_antenna_layout == 6 || 
			  params.bore_hole_antenna_layout == 7 ) {   // it's V-V-V-V or V-V or V
                    
                    for (int j=0; j<params.number_of_strings_station; j++) {
                        for (int k=0; k<params.number_of_antennas_string; k++) {
			  
			  if (settings1->BH_ANT_SEP_DIST_ON==0) 
                            stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw*k );
			  
			  else if (settings1->BH_ANT_SEP_DIST_ON==1) {
			    z_btw_total = 0.;
			    for (int l=0; l<k+1; l++) {
			      z_btw_total += z_btw_array[l];
			    }
			    stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw_total );
			  }

			  stations[i].strings[j].antennas[k].type = 0;   // all antennas v-pol			  
                            
			  if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
			    stations[i].strings[j].antennas[k].orient = 0;
			  }
			  else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
			    if ( j==0||j==3 ) {
			      if ( k==0||k==1 ) {
				stations[i].strings[j].antennas[k].orient = 0;
			      }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                }
                                else {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                }
                                
                            }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                        }
                    }
                    
                } // end if bore hole antenna layout = 5,6,7 (VVVV, VV, V)
                
                //
                // set surface antenna postions
                stations[i].surfaces[0].SetX( stations[i].GetX() + (R_surface * cos(PI/3.)) );
                stations[i].surfaces[0].SetY( stations[i].GetY() + (R_surface * sin(PI/3.)) );
                
                stations[i].surfaces[1].SetX( stations[i].GetX() + (R_surface * cos(-PI/3.)) );
                stations[i].surfaces[1].SetY( stations[i].GetY() + (R_surface * sin(-PI/3.)) );
                
                stations[i].surfaces[2].SetX( stations[i].GetX() + (R_surface * cos(PI)) );
                stations[i].surfaces[2].SetY( stations[i].GetY() );
                
                stations[i].surfaces[3].SetX( stations[i].GetX() );
                stations[i].surfaces[3].SetY( stations[i].GetY() );
                

                stations[i].number_of_antennas = params.number_of_strings_station * params.number_of_antennas_string;


            }// loop over stations i


            // for idealized geometry, number of antennas in a station is constant
            max_number_of_antennas_station = params.number_of_strings_station * params.number_of_antennas_string;




        } // if idealized geometry
#ifdef ARA_UTIL_EXISTS

        else { // non-idealized geometry

	  //SetupInstalledStations();        
            //for (int i=0; i<params.number_of_stations; i++) {
      
                //AraGeomTool *araGeom=AraGeomTool::Instance();
                AraGeomTool *araGeom = new AraGeomTool();
                cout<<"read AraGeomTool"<<endl;
                
                for (int i=0; i<params.number_of_stations; i++) {
                    for (int j = 0; j < params.number_of_strings_station; j++){
                        
                        double avgX, avgY;
                        
                        for (int k = 0; k < params.number_of_antennas_string; k++){
                            
                            //int chan = GetChannelfromStringAntenna (i+1,j,k);
                            int chan = GetChannelfromStringAntenna (i+1,j,k,settings1);
                            
                            stations[i].strings[j].antennas[k].SetX(stations[i].GetX()+araGeom->getStationInfo(i+1)->fAntInfo[chan-1].antLocation[0]);
                            stations[i].strings[j].antennas[k].SetY(stations[i].GetY()+araGeom->getStationInfo(i+1)->fAntInfo[chan-1].antLocation[1]);
                            //stations[i].strings[j].antennas[k].SetZ(araGeom->fStationInfo[i+1].fAntInfo[chan-1].antLocation[2]-double(settings1->DEPTH_CHANGE));
                            stations[i].strings[j].antennas[k].SetZ(araGeom->getStationInfo(i+1)->fAntInfo[chan-1].antLocation[2]);
                                                    cout <<
                             "DetectorStation:string:antenna:X:Y:Z:: " <<
                             i<< " : " <<
                             j<< " : " <<
                             k<< " : " <<
                             stations[i].strings[j].antennas[k].GetX() << " : " <<
                             stations[i].strings[j].antennas[k].GetY() << " : " <<
                             stations[i].strings[j].antennas[k].GetZ() << " : " <<
			     chan << " : " <<	
			     //araGeom->fStationInfo[i+1].fAntInfo[chan-1].antLocation[2]-double(settings1->DEPTH_CHANGE) << " : " <<
			     //double(settings1->DEPTH_CHANGE) << " : " <<
                             endl;
                             
                        }
                        
                        //int chanstring = GetChannelfromStringAntenna (i+1, j,2);
                        int chanstring = GetChannelfromStringAntenna (i+1, j,2,settings1);
                        
                        stations[i].strings[j].SetX(stations[i].GetX()+araGeom->getStationInfo(i+1)->fAntInfo[chanstring-1].antLocation[0]);
                        stations[i].strings[j].SetY(stations[i].GetY()+araGeom->getStationInfo(i+1)->fAntInfo[chanstring-1].antLocation[1]);
                        
                    }
                
                
                //
                // set antenna postions in borehole
                // and set type (h or v pol antenna) and set orientation (facing x or y)
                if ( params.bore_hole_antenna_layout == 0 || params.bore_hole_antenna_layout == 1) {
                    
                    for (int j=0; j<params.number_of_strings_station; j++) {
                        for (int k=0; k<params.number_of_antennas_string; k++) {

                            if (k%2 == 0) {
                                stations[i].strings[j].antennas[k].type = 0;   // v-pol
                            }
                            else {
                                stations[i].strings[j].antennas[k].type = 1;   // h-pol
                            }
                            
                            if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                                stations[i].strings[j].antennas[k].orient = 0;
                            }
                            else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                                if ( j==0||j==3 ) {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                }
                                else {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                }
                                
                            }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                        }
                    }
                    
                } // end if bore hole antenna layout = 0 or 1 (where VHVH way but different numbers)
                
                
                else if ( params.bore_hole_antenna_layout == 2) {   // it's V-H-V-V
                    
                    for (int j=0; j<params.number_of_strings_station; j++) {
                        for (int k=0; k<params.number_of_antennas_string; k++) {
                            
                            if (k == 1) {   // only the second antenna is H pol
                                stations[i].strings[j].antennas[k].type = 1;   // h-pol
                            }
                            else {  // other antennas are V pol
                                stations[i].strings[j].antennas[k].type = 0;   // v-pol
                            }
                            
                            if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                                stations[i].strings[j].antennas[k].orient = 0;
                            }
                            else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                                if ( j==0||j==3 ) {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                }
                                else {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                }
                                
                            }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                        }
                    }
                    
                } // end if bore hole antenna layout = 2 (where VHVV way but different numbers)
                
                
                else if ( params.bore_hole_antenna_layout == 3 || params.            
                         bore_hole_antenna_layout == 4 ) {   // it's V-H-H-H or V-H-H
                    
                    for (int j=0; j<params.number_of_strings_station; j++) {
                        for (int k=0; k<params.number_of_antennas_string; k++) {
                            
                            if (k == 0) {   // only the first antenna is V pol
                                stations[i].strings[j].antennas[k].type = 0;   // v-pol
                            }
                            else {  // other antennas are H pol
                                stations[i].strings[j].antennas[k].type = 1;   // h-pol
                            }
                            
                            if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                                stations[i].strings[j].antennas[k].orient = 0;
                            }
                            else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                                if ( j==0||j==3 ) {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                }
                                else {
                                    if ( k==0||k==1 ) {
                                        stations[i].strings[j].antennas[k].orient = 1;
                                    }
                                    else {
                                        stations[i].strings[j].antennas[k].orient = 0;
                                    }
                                }
                                
                            }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                        }
                    }
                    
                } // end if bore hole antenna layout = 3 (where VHHH way)
                
                
                
                
                //
                // set surface antenna postions
                stations[i].surfaces[0].SetX( stations[i].GetX() + (R_surface * cos(PI/3.)) );
                stations[i].surfaces[0].SetY( stations[i].GetY() + (R_surface * sin(PI/3.)) );
                
                stations[i].surfaces[1].SetX( stations[i].GetX() + (R_surface * cos(-PI/3.)) );
                stations[i].surfaces[1].SetY( stations[i].GetY() + (R_surface * sin(-PI/3.)) );
                
                stations[i].surfaces[2].SetX( stations[i].GetX() + (R_surface * cos(PI)) );
                stations[i].surfaces[2].SetY( stations[i].GetY() );
                
                stations[i].surfaces[3].SetX( stations[i].GetX() );
                stations[i].surfaces[3].SetY( stations[i].GetY() );

                
            } // end loop over stations i


            //}// end loop over stations i
            

            int antenna_count = 0;
            max_number_of_antennas_station = 0;
            // for non-idealized geometry, it's better to actually count number of stations
            for (int i=0; i<(int)(stations.size()); i++) {
            
                antenna_count = 0;
                for (int j=0; j<(int)(stations[i].strings.size()); j++) {
                    for (int k=0; k<(int)(stations[i].strings[j].antennas.size()); k++) {
                        antenna_count++;
                    }
                }
                stations[i].number_of_antennas = antenna_count;

                if (max_number_of_antennas_station < antenna_count) max_number_of_antennas_station = antenna_count;
            }

            
        }// if non-idealized geom
#endif
        
        

	if (settings1->ANTENNA_MODE == 0){
	  // test read V-pol gain file!!
	  ReadVgain("ARA_bicone6in_output.txt");
	  // test read H-pol gain file!!
	  ReadHgain("ARA_dipoletest1_output.txt");
	}
	else if (settings1->ANTENNA_MODE == 1) {
	  // test read V-pol gain file!!
	  ReadVgainSettings("ARA_bicone6in_output_updated2016.txt", settings1);
	  ReadVgainTopSettings("ARA_VPresult_topTrec.txt", settings1);
	  // test read H-pol gain file!!
	  ReadHgainSettings("ARA_dipoletest1_output_updated2016.txt", settings1);
	}
	else if (settings1->ANTENNA_MODE == 2){
	  // test read V-pol gain file!!
	  ReadVgain("Arianna_WIPLD_hpol.dat");
	  // test read H-pol gain file!!
	  ReadHgain("Arianna_WIPLD_hpol.dat");
	}
	
	//	if (settings1->NOISE == 2){
	  ReadNoiseFigure("./data/ARA02_noiseFig.txt", settings1);
	  //	}


        // read filter file!!
        ReadFilter("./data/filter.csv", settings1);
        // read preamp gain file!!
        ReadPreamp("./data/preamp.csv", settings1);
        // read FOAM gain file!!
        ReadFOAM("./data/FOAM.csv", settings1);
        // read gain offset for chs file!!
        ReadGainOffset_TestBed("./data/preamp_ch_gain_offset.csv", settings1);// only TestBed for now
        // read threshold offset for chs file!!
        ReadThresOffset_TestBed("./data/threshold_offset.csv", settings1);// only TestBed for now
	// read threshold values for chs file
	ReadThres_TestBed("./data/thresholds_TB.csv", settings1);// only TestBed for now
        // read system temperature for chs file!!
	cout << "check read testbed temp1" << endl;
        if (settings1->NOISE_CHANNEL_MODE != 0) {
	  
            ReadTemp_TestBed("./data/system_temperature.csv", settings1);// only TestBed for now

        }
        // read total elec. chain response file!!
        cout<<"start read elect chain"<<endl;
        if(settings1->CUSTOM_ELECTRONICS==0){
            //read the standard ARA electronics
            cout<<"     Reading standard ARA electronics response"<<endl;
             ReadElectChain("./data/ARA_Electronics_TotalGain_TwoFilters.txt", settings1);
            //ReadElectChain("./data/ARA_Electronics_TotalGainPhase.txt", settings1);
        }
        else if (settings1->CUSTOM_ELECTRONICS==1){
            //read a custom user defined electronics gain
            cout<<"     Reading custom electronics response"<<endl;
             ReadElectChain("./data/custom_electronics.txt", settings1);
        }
        cout<<"done read elect chain"<<endl;
        
        
    } // if mode == 1
    
    
    
    
    /////////////////////////////////////////////////////////////////////////////////    
    
    
    
    else if (mode == 2) {
        cout<<"\n\tDector mode 2 : Pentagon"<<endl;
        cout<<"\n\tBy default, ARA-37 is set"<<endl;
        ifstream ARA37( ARA37_file.c_str() );
        cout<<"We use "<<ARA37_file.c_str()<<" as antenna info."<<endl;

        //SetupInstalledStations();

        
        //
        // initialize info
        params.number_of_stations = 37;
        params.number_of_strings_station = 4;   // ARA-1 has 4 strings
        params.number_of_antennas_string = 4; // 4 antennas on each strings
        params.number_of_surfaces_station = 4;
        
        //double core_x = 0.;  // all units are in meter
        //double core_y = 0.;
        params.core_x = 10000.;  // all units are in meter
        params.core_y = 10000.;
        double R_string = 10.;
        double R_surface = 60.;
        double z_max = 200.;
        double z_btw = 10.;
        double z_btw_array[6]; // assume there will be less than 6 bore hole antennas at each string
	// these z_btw array will be used when settings->BH_ANT_SEP_DIST_ON=1 case
        for (int i=0; i<6; i++) {
            if (i==0) z_btw_array[i] = 0.;
            //else z_btw_array[i] = z_btw;
            else if (i==1) z_btw_array[i] = 2.;
            else if (i==2) z_btw_array[i] = 15.;
            else if (i==3) z_btw_array[i] = 2.;
            else z_btw_array[i] = z_btw;
        }
        double z_btw_total;
        params.stations_per_side = 4;       // total 37 stations
        params.station_spacing = 2000.;     // 2km spacing
        params.antenna_orientation = 0;     // all antenna facing x
        params.bore_hole_antenna_layout = settings1->BORE_HOLE_ANTENNA_LAYOUT;
        // finish initialization
        //
        
        
        
        
        
        
        // Read new parameters if there are...
        if ( ARA37.is_open() ) {
            while (ARA37.good() ) {
                getline (ARA37, line);
                
                if (line[0] != "/"[0]) {
                    label = line.substr(0, line.find_first_of("=") );
                    
                    if (label == "core_x") {
                        params.core_x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_x"<<endl;
                    }
                    else if (label == "core_y") {
                        params.core_y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_y"<<endl;
                    }
                    else if (label == "R_string") {
                        R_string = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_string"<<endl;
                    }
                    else if (label == "R_surface") {
                        R_surface = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_surface"<<endl;
                    }
                    else if (label == "z_max") {
                        z_max = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_max"<<endl;
                    }
                    else if (label == "z_btw") {
                        z_btw = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw"<<endl;
                    }
                    else if (label == "z_btw01") {
                        z_btw_array[1] = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw bh ant0 and ant1"<<endl;
                    }
                    else if (label == "z_btw12") {
                        z_btw_array[2] = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw bh ant1 and ant2"<<endl;
                    }
                    else if (label == "z_btw23") {
                        z_btw_array[3] = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw bh ant2 and ant3"<<endl;
                    }
                    else if (label == "z_btw34") {
                        z_btw_array[4] = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw bh ant3 and ant4"<<endl;
                    }
                    else if (label == "z_btw45") {
                        z_btw_array[5] = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw bh ant4 and ant5"<<endl;
                    }
                    else if (label == "stations_per_side") {
                        params.stations_per_side = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read stations_per_side"<<endl;
                    }
                    else if (label == "station_spacing") {
                        params.station_spacing = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read station_spacting"<<endl;
                    }
                    else if (label == "antenna_orientation") {
                        params.antenna_orientation = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read antenna_orientation"<<endl;
                    }
                }
            }
            ARA37.close();
        }
        // finished reading new parameters
        
        
        
        // set number of antennas in a string
        if (params.bore_hole_antenna_layout == 0) { // VHVH layout
            params.number_of_antennas_string = 4;
        }
        else if (params.bore_hole_antenna_layout == 1) { // VHV layout
            params.number_of_antennas_string = 3;
        }
        else if (params.bore_hole_antenna_layout == 2) { // VHVV layout
            params.number_of_antennas_string = 4;
        }
        else if (params.bore_hole_antenna_layout == 3) { // VHHH layout
            params.number_of_antennas_string = 4;
        }
        else if (params.bore_hole_antenna_layout == 4) { // VHH layout
            params.number_of_antennas_string = 3;
        }
        
        
        
        
        
        //
        // caculate number of stations, strings, antennas 
        params.number_of_stations = 1 + (3 * params.stations_per_side) * (params.stations_per_side - 1);
        
        params.number_of_strings = params.number_of_stations * params.number_of_strings_station;
        params.number_of_antennas = params.number_of_strings * params.number_of_antennas_string;
        // 
        
        
        
        //
        // prepare vectors
        for (int i=0; i<params.number_of_stations; i++) {
            stations.push_back(temp_station);
            
            for (int j=0; j<params.number_of_surfaces_station; j++) {
                stations[i].surfaces.push_back(temp_surface);
            }
            
            for (int k=0; k<params.number_of_strings_station; k++) {
                stations[i].strings.push_back(temp_string);
                
                for (int l=0; l<params.number_of_antennas_string; l++) {
                    stations[i].strings[k].antennas.push_back(temp_antenna);
                }
                
            }
            
            
        }
        // end perpare vectors
        //
        
        
        
        
        
        
        
        
        //
        // for ARA-37 (or more than 1 station case), need code for setting position for all 37 stations here!
        //
        //
        // here, this only works for pentagon shape!
        //
        double y_offset = (double)params.station_spacing * sqrt(3) / 2.;
        
        int station_count = 0;
        
        for (int irow = 0; irow < ((int)params.stations_per_side * 2)-1; irow++) {
            double current_y = y_offset * ( (double)params.stations_per_side - 1 - irow) + params.core_y;
            int stations_this_row = (2 * (int)params.stations_per_side - 1) - abs((int)params.stations_per_side - 1 - irow);
            
            for (int istation = 0; istation < stations_this_row; istation++) {
                if (station_count < (int)params.number_of_stations) {
                    stations[station_count].SetY( current_y );
                    stations[station_count].SetX( (double)params.station_spacing * ((double)istation - ((double)stations_this_row - 1.) / 2.) + params.core_x );
                    station_count++;
                }
                else {
                    cout<<"\n\tError, too many stations !"<<endl;
                }
            }
        }
        // finished setting all stations' position
        
        
        cout<<"total station_count : "<<station_count<<endl;
        if (station_count != (int)params.number_of_stations) cout<<"\n\tError, station number not match !"<<endl;        
        
        //
        // set antenna values from parameters
        // set station positions
        for (int i=0; i<params.number_of_stations; i++) {
            
            //
            // set string postions based on station position
            //            for (int j=0; j<params.number_of_strings_station; j++) {
            //            stations[i].string[0].x = stations[i].x - (R_string / 1.414);
            stations[i].strings[0].SetX( stations[i].GetX() - (R_string * cos(PI/4.)) );
            stations[i].strings[0].SetY( stations[i].GetY() + (R_string * sin(PI/4.)) );
            
            stations[i].strings[1].SetX( stations[i].GetX() + (R_string * cos(PI/4.)) );
            stations[i].strings[1].SetY( stations[i].GetY() + (R_string * sin(PI/4.)) );
            
            stations[i].strings[2].SetX( stations[i].GetX() - (R_string * cos(PI/4.)) );
            stations[i].strings[2].SetY( stations[i].GetY() - (R_string * sin(PI/4.)) );
            
            stations[i].strings[3].SetX( stations[i].GetX() + (R_string * cos(PI/4.)) );
            stations[i].strings[3].SetY( stations[i].GetY() - (R_string * sin(PI/4.)) );
            
            
            
            
            //
            // set antenna postions in borehole
            // and set type (h or v pol antenna) and set orientation (facing x or y)
            if ( params.bore_hole_antenna_layout == 0 || params.bore_hole_antenna_layout == 1) {
                for (int j=0; j<params.number_of_strings_station; j++) {
                    for (int k=0; k<params.number_of_antennas_string; k++) {

                        if (settings1->BH_ANT_SEP_DIST_ON==0) 
                        stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw*k );

                        else if (settings1->BH_ANT_SEP_DIST_ON==1) {
                            z_btw_total = 0.;
                            for (int l=0; l<k+1; l++) {
                                z_btw_total += z_btw_array[l];
                            }
                            stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw_total );
                        }
                        
                        if (k%2 == 0) {
                            stations[i].strings[j].antennas[k].type = 0;   // v-pol
                        }
                        else {
                            stations[i].strings[j].antennas[k].type = 1;   // h-pol
                        }
                        
                        if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                            stations[i].strings[j].antennas[k].orient = 0;
                        }
                        else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                            if ( j==0||j==3 ) {
                                if ( k==0||k==1 ) {
                                    stations[i].strings[j].antennas[k].orient = 0;
                                }
                                else {
                                    stations[i].strings[j].antennas[k].orient = 1;
                                }
                            }
                            else {
                                if ( k==0||k==1 ) {
                                    stations[i].strings[j].antennas[k].orient = 1;
                                }
                                else {
                                    stations[i].strings[j].antennas[k].orient = 0;
                                }
                            }
                            
                        }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                    }
                }
            } // end if bore hole antenna layout = 0 or 1 (where VHVH way but different numbers)
            
            
            else if ( params.bore_hole_antenna_layout == 2) {   // it's V-H-V-V
                for (int j=0; j<params.number_of_strings_station; j++) {
                    for (int k=0; k<params.number_of_antennas_string; k++) {

                        if (settings1->BH_ANT_SEP_DIST_ON==0) 
                        stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw*k );

                        else if (settings1->BH_ANT_SEP_DIST_ON==1) {
                            z_btw_total = 0.;
                            for (int l=0; l<k+1; l++) {
                                z_btw_total += z_btw_array[l];
                            }
                            stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw_total );
                        }
                        
                        if (k == 1) {   // only the second antenna is H pol
                            stations[i].strings[j].antennas[k].type = 1;   // h-pol
                        }
                        else {  // other antennas are V pol
                            stations[i].strings[j].antennas[k].type = 0;   // v-pol
                        }
                        
                        if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                            stations[i].strings[j].antennas[k].orient = 0;
                        }
                        else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                            if ( j==0||j==3 ) {
                                if ( k==0||k==1 ) {
                                    stations[i].strings[j].antennas[k].orient = 0;
                                }
                                else {
                                    stations[i].strings[j].antennas[k].orient = 1;
                                }
                            }
                            else {
                                if ( k==0||k==1 ) {
                                    stations[i].strings[j].antennas[k].orient = 1;
                                }
                                else {
                                    stations[i].strings[j].antennas[k].orient = 0;
                                }
                            }
                            
                        }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                    }
                }
            } // end if bore hole antenna layout = 0 or 1 (where VHVH way but different numbers)
            
            
            
            else if ( params.bore_hole_antenna_layout == 3 || params.bore_hole_antenna_layout == 4 ) {   // it's V-H-H-H or V-H-H
                for (int j=0; j<params.number_of_strings_station; j++) {
                    for (int k=0; k<params.number_of_antennas_string; k++) {

                        if (settings1->BH_ANT_SEP_DIST_ON==0) 
                        stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw*k );

                        else if (settings1->BH_ANT_SEP_DIST_ON==1) {
                            z_btw_total = 0.;
                            for (int l=0; l<k+1; l++) {
                                z_btw_total += z_btw_array[l];
                            }
                            stations[i].strings[j].antennas[k].SetZ( -z_max + z_btw_total );
                        }
                        
                        if (k == 0) {   // only the first antenna is V pol
                            stations[i].strings[j].antennas[k].type = 0;   // v-pol
                        }
                        else {  // other antennas are H pol
                            stations[i].strings[j].antennas[k].type = 1;   // h-pol
                        }
                        
                        if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                            stations[i].strings[j].antennas[k].orient = 0;
                        }
                        else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                            if ( j==0||j==3 ) {
                                if ( k==0||k==1 ) {
                                    stations[i].strings[j].antennas[k].orient = 0;
                                }
                                else {
                                    stations[i].strings[j].antennas[k].orient = 1;
                                }
                            }
                            else {
                                if ( k==0||k==1 ) {
                                    stations[i].strings[j].antennas[k].orient = 1;
                                }
                                else {
                                    stations[i].strings[j].antennas[k].orient = 0;
                                }
                            }
                            
                        }// end facing different. I know it only works with 4 strings, 4 antennas on each strings but couldn't find a better way than this. -Eugene
                    }
                }
            } // end if bore hole antenna layout = 3 (where VHHH way)
            
            
            
            
            
            
            
            //
            // set surface antenna postions
            stations[i].surfaces[0].SetX( stations[i].GetX() + (R_surface * cos(PI/3.)) );
            stations[i].surfaces[0].SetY( stations[i].GetY() + (R_surface * sin(PI/3.)) );
            
            stations[i].surfaces[1].SetX( stations[i].GetX() + (R_surface * cos(-PI/3.)) );
            stations[i].surfaces[1].SetY( stations[i].GetY() + (R_surface * sin(-PI/3.)) );
            
            stations[i].surfaces[2].SetX( stations[i].GetX() + (R_surface * cos(PI)) );
            //            stations[i].surfaces[2].y = stations[i].y + (R_surface * sin(PI));
            stations[i].surfaces[2].SetY( stations[i].GetY() );
            
            stations[i].surfaces[3].SetX( stations[i].GetX() );
            stations[i].surfaces[3].SetY( stations[i].GetY() );


            stations[i].number_of_antennas = params.number_of_strings_station * params.number_of_antennas_string;
            
        }// loop over stations i
        
        
        // for idealized geometry, number of antennas in a station is constant
        max_number_of_antennas_station = params.number_of_strings_station * params.number_of_antennas_string;
        

        

	if (settings1->ANTENNA_MODE == 0){
	  // test read V-pol gain file!!
	  ReadVgain("ARA_bicone6in_output.txt");
	  // test read H-pol gain file!!
	  ReadHgain("ARA_dipoletest1_output.txt");
	}
	else if (settings1->ANTENNA_MODE == 1) {
	  // test read V-pol gain file!!
	  ReadVgainSettings("ARA_bicone6in_output_updated2016.txt", settings1);
	  ReadVgainTopSettings("ARA_VPresult_topTrec.txt", settings1);
	  // test read H-pol gain file!!
	  ReadHgainSettings("ARA_dipoletest1_output_updated2016.txt", settings1);
	}
	else if (settings1->ANTENNA_MODE == 2){
	  // test read V-pol gain file!!
	  ReadVgain("Arianna_WIPLD_hpol.dat");
	  // test read H-pol gain file!!
	  ReadHgain("Arianna_WIPLD_hpol.dat");
	}

	//	if (settings1->NOISE==2){
	ReadNoiseFigure("./data/ARA02_noiseFig.txt", settings1);
	  //	}

        ReadFilter("./data/filter.csv", settings1);
        // read preamp gain file!!
        ReadPreamp("./data/preamp.csv", settings1);
        // read FOAM gain file!!
        ReadFOAM("./data/FOAM.csv", settings1);
        // read gain offset for chs file!!
        ReadGainOffset_TestBed("./data/preamp_ch_gain_offset.csv", settings1);// only TestBed for now
        // read threshold offset for chs file!!
        ReadThresOffset_TestBed("./data/threshold_offset.csv", settings1);// only TestBed for now
	// read threshold values for chs file
	ReadThres_TestBed("./data/thresholds_TB.csv", settings1);// only TestBed for now
        // read system temperature for chs file!!

	cout << "check read temp testbed 2" << endl;
       if (settings1->NOISE_CHANNEL_MODE != 0) {
            ReadTemp_TestBed("./data/system_temperature.csv", settings1);// only TestBed for now
        }
        // read total elec. chain response file!!
        cout<<"start read elect chain"<<endl;
        if(settings1->CUSTOM_ELECTRONICS==0){
            //read the standard ARA electronics
            cout<<"     Reading standard ARA electronics response"<<endl;
             ReadElectChain("./data/ARA_Electronics_TotalGain_TwoFilters.txt", settings1);
            //ReadElectChain("./data/ARA_Electronics_TotalGainPhase.txt", settings1);
        }
        else if (settings1->CUSTOM_ELECTRONICS==1){
            //read a custom user defined electronics gain
            cout<<"     Reading custom electronics response"<<endl;
             ReadElectChain("./data/custom_electronics.txt", settings1);
        }
        
        
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////    
    else if (mode == 3) {        //        cout<<"\n\tDector mode 3 : Testbed and eventual inclusion of a specific number of stations (less than 7 stations) !"<<endl;
        //        cout<<"We use "<<ARA_N_file.c_str()<<" as antenna info."<<endl;
        
        //SetupInstalledStations();        
        
        // initialize info
        params.number_of_stations = 1; //including Testbed
        params.number_of_strings_station = 4;   // ARA-1 has 4 strings
        params.number_of_antennas_string = 4; // 4 antennas on each strings
        params.number_of_surfaces_station = 4;
        params.number_of_channels = 20;
        
        //double core_x = 0.;
        //double core_y = 0.;
        params.core_x = 10000.;
        params.core_y = 10000.;
        double R_string = 10.;  // all units are in meter
        double R_surface = 60.;
        double z_max = 200.;
        double z_btw = 20.;
        params.stations_per_side = 4;       // total 37 stations
        params.station_spacing = 2000.;     // 2km spacing for borehole stations
        params.antenna_orientation = 0;     // all antenna facing x
        params.bore_hole_antenna_layout = settings1->BORE_HOLE_ANTENNA_LAYOUT;
        // finish initialization
        //
        

        
        // mode == 3 currently just use installed TestBed station geom information.
        // So don't need to read any more information
        
        // Read new parameters if there are...
        ifstream ARA_N( ARA_N_file.c_str() );
        if ( ARA_N.is_open() ) {
            while (ARA_N.good() ) {
                getline (ARA_N, line);
                
                if (line[0] != "/"[0]) {
                    label = line.substr(0, line.find_first_of("=") );
                    
                    if (label == "core_x") {
                        params.core_x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_x"<<endl;
                    }
                    else if (label == "core_y") {
                        params.core_y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_y"<<endl;
                    }
                    else if (label == "R_string") {
                        R_string = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_string"<<endl;
                    }
                    else if (label == "R_surface") {
                        R_surface = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_surface"<<endl;
                    }
                    else if (label == "z_max") {
                        z_max = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_max"<<endl;
                    }
                    else if (label == "z_btw") {
                        z_btw = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw"<<endl;
                    }
                    else if (label == "number_of_stations") {
                        params.number_of_stations = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read stations_per_side"<<endl;
                    }
                    else if (label == "station_spacing") {
                        params.station_spacing = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read station_spacting"<<endl;
                    }
                    else if (label == "antenna_orientation") {
                        params.antenna_orientation = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read antenna_orientation"<<endl;
                    }
                }
            }
            ARA_N.close();
        }
        // finished reading new parameters
        
        
        
        params.number_of_antennas_string = 4;

        
        // prepare vectors
        PrepareVectorsInstalled();
        // end prepare vectors
        
        //
        // for ARA-37 (or more than 1 station case), need code for setting position for all 37 stations here!
        //
        int station_count = 0;
        
        for (int istation = 0; istation < (int)params.number_of_stations; istation++) {
            if (station_count < (int)params.number_of_stations - 1) {
                //stations[station_count].x = core_x + (double)params.station_spacing * cos( (PI/3.) * (double)station_count );
                //stations[station_count].y = core_y + (double)params.station_spacing * sin( (PI/3.) * (double)station_count );
                stations[station_count].SetX( params.core_x + (double)params.station_spacing * cos( (PI/3.) * (double)station_count ) );
                stations[station_count].SetY( params.core_y + (double)params.station_spacing * sin( (PI/3.) * (double)station_count ) );
                station_count++;
            }
            else if (station_count < (int)params.number_of_stations) {
                //stations[station_count].x = core_x;
                //stations[station_count].y = core_y;
                stations[station_count].SetX( params.core_x );
                stations[station_count].SetY( params.core_y );
                station_count++;
            }
            else {
                cout<<"\n\tError, too many stations !"<<endl;
            }
        }
        // finished setting all stations' position
        
        
        //        cout<<"total station_count : "<<station_count<<endl;
        if (station_count != (int)params.number_of_stations) cout<<"\n\tError, station number not match !"<<endl;
        
        //
        // set antenna values from parameters
        // set station positions
        //cout << "READGEOM:" << settings1->READGEOM << endl;
        
#ifdef ARA_UTIL_EXISTS
        UseAntennaInfo(0, settings1);
#endif
//            UseAntennaInfo(1, settings1);
        for (int i = 0; i < (int)params.number_of_stations; i++){
            stations[i].StationID = i;
            if (settings1->USE_INSTALLED_TRIGGER_SETTINGS == 0){
                stations[i].NFOUR = 1024;
                stations[i].TIMESTEP = 1./2.6*1.E-9;
                stations[i].TRIG_WINDOW = 2.5E-7;
                stations[i].DATA_BIN_SIZE = settings1->DATA_BIN_SIZE;
            }
            else if (settings1->USE_INSTALLED_TRIGGER_SETTINGS == 1){
                if (stations[i].StationID == 0){
                    stations[i].NFOUR = 1024;
                    stations[i].TIMESTEP = 1./2.6*1.E-9;
                    stations[i].TRIG_WINDOW = 2.5E-7;
                    stations[i].DATA_BIN_SIZE = settings1->DATA_BIN_SIZE;
                }
               if (stations[i].StationID == 1){
                    stations[i].NFOUR = 1024;
                    stations[i].TIMESTEP = 1./2.6*1.E-9;
                    stations[i].TRIG_WINDOW = 2.5E-7;
                    stations[i].DATA_BIN_SIZE = settings1->DATA_BIN_SIZE;
                }
            }
        }
        
        params.number_of_antennas = 0;

            cout<<"DETECTOR=3 TB station geom info"<<endl;
        
            for (int j = 0; j < stations[0].strings.size(); j++){
                for (int k = 0; k < stations[0].strings[j].antennas.size(); k++){
                    
                     cout <<
                     "DetectorStation2:string:antenna:X:Y:Z:chno :: " <<
                     j<< " : " <<
                     k<< " : " <<
                     stations[0].strings[j].antennas[k].GetX() << " : " <<
                     stations[0].strings[j].antennas[k].GetY() << " : " <<
                     stations[0].strings[j].antennas[k].GetZ() << " : \t" <<
                     //GetChannelfromStringAntenna ( 0, j, k)<<
                     GetChannelfromStringAntenna ( 0, j, k, settings1)<<
                     endl;

                     params.number_of_antennas++;
                }
            }
        



            cout<<"after FlattoEarth, station0 location"<<endl;
    for (int j = 0; j < stations[0].strings.size(); j++){
        for (int k = 0; k < stations[0].strings[j].antennas.size(); k++){

             cout <<
             "Detector:station:string:antenna:X:Y:Z:R:Theta:Phi:: " <<
             "0" << " : " <<
             j<< " : " <<
             k<< " : " <<
             stations[0].strings[j].antennas[k].GetX() << " : " <<
             stations[0].strings[j].antennas[k].GetY() << " : " <<
             stations[0].strings[j].antennas[k].GetZ() << " : " <<
             stations[0].strings[j].antennas[k].R() << " : " <<
             stations[0].strings[j].antennas[k].Theta() << " : " <<
             stations[0].strings[j].antennas[k].Phi() << " : " <<
             icesurface->Surface(stations[0].strings[j].antennas[k].Lon(), stations[0].strings[j].antennas[k].Lat()) << " : " <<
//             icesurface->Surface(stations[0].strings[j].antennas[k].Lat(), stations[0].strings[j].antennas[k].Lon()) << " : " <<
             endl;
                 
                     
        }
    }



            int antenna_count = 0;
            max_number_of_antennas_station = 0;
            // for non-idealized geometry, it's better to actually count number of stations
            for (int i=0; i<(int)(stations.size()); i++) {
            
                antenna_count = 0;
                for (int j=0; j<(int)(stations[i].strings.size()); j++) {
                    for (int k=0; k<(int)(stations[i].strings[j].antennas.size()); k++) {
                        antenna_count++;
                    }
                }
                stations[i].number_of_antennas = antenna_count;

                if (max_number_of_antennas_station < antenna_count) max_number_of_antennas_station = antenna_count;
            }


            
            



	    if (settings1->ANTENNA_MODE == 0){
	      // test read V-pol gain file!!
	      ReadVgain("ARA_bicone6in_output.txt");
	      // test read H-pol gain file!!
	      ReadHgain("ARA_dipoletest1_output.txt");
	    }
	    else if (settings1->ANTENNA_MODE == 1) {
	      // test read V-pol gain file!!
	      ReadVgainSettings("ARA_bicone6in_output_updated2016.txt", settings1);
	      ReadVgainTopSettings("ARA_VPresult_topTrec.txt", settings1);
	      // test read H-pol gain file!!
	      ReadHgainSettings("ARA_dipoletest1_output_updated2016.txt", settings1);
	    }
	    else if (settings1->ANTENNA_MODE == 2) {
	      // test read V-pol gain file!!
	      ReadVgain("Arianna_WIPLD_hpol.dat");
	      // test read H-pol gain file!!
	      ReadHgain("Arianna_WIPLD_hpol.dat");
	    }
	    

	    ReadNoiseFigure("./data/ARA02_noiseFig.txt", settings1);
	      
            // read filter file!!
            ReadFilter("./data/filter.csv", settings1);
            // read preamp gain file!!
            ReadPreamp("./data/preamp.csv", settings1);
            // read FOAM gain file!!
            ReadFOAM("./data/FOAM.csv", settings1);



            if ( settings1->NOISE==1) {
                // read Rayleigh fit for freq range, bh channels
                //ReadRayleighFit_TestBed("data/RayleighFit_TB.csv", settings1, settings1->RAYLEIGH_CH); // read and save RFCM gain
                ReadRayleighFit_TestBed("data/RayleighFit_TB.csv", settings1); // read and save RFCM gain
            }

            if ( settings1->USE_TESTBED_RFCM_ON==1) {
                // read RFCM gain file!! (measured value in ICL)
                ReadRFCM_TestBed("data/TestBed_RFCM/R1C1.csv", settings1); // read and save RFCM gain for ch1
                ReadRFCM_TestBed("data/TestBed_RFCM/R1C2.csv", settings1); // read and save RFCM gain for ch2
                ReadRFCM_TestBed("data/TestBed_RFCM/R1C3.csv", settings1); // read and save RFCM gain for ch3
                ReadRFCM_TestBed("data/TestBed_RFCM/R1C4.csv", settings1); // read and save RFCM gain for ch4
                ReadRFCM_TestBed("data/TestBed_RFCM/R2C5.csv", settings1); // read and save RFCM gain for ch5
                ReadRFCM_TestBed("data/TestBed_RFCM/R2C6.csv", settings1); // read and save RFCM gain for ch6
                ReadRFCM_TestBed("data/TestBed_RFCM/R2C7.csv", settings1); // read and save RFCM gain for ch7
                ReadRFCM_TestBed("data/TestBed_RFCM/R2C8.csv", settings1); // read and save RFCM gain for ch8
                ReadRFCM_TestBed("data/TestBed_RFCM/R3C9.csv", settings1); // read and save RFCM gain for ch9
                ReadRFCM_TestBed("data/TestBed_RFCM/R3C10.csv", settings1); // read and save RFCM gain for ch10
                ReadRFCM_TestBed("data/TestBed_RFCM/R3C11.csv", settings1); // read and save RFCM gain for ch11
                ReadRFCM_TestBed("data/TestBed_RFCM/R3C12.csv", settings1); // read and save RFCM gain for ch12
                ReadRFCM_TestBed("data/TestBed_RFCM/R4C13.csv", settings1); // read and save RFCM gain for ch13
                ReadRFCM_TestBed("data/TestBed_RFCM/R4C14.csv", settings1); // read and save RFCM gain for ch14
                ReadRFCM_TestBed("data/TestBed_RFCM/R4C15.csv", settings1); // read and save RFCM gain for ch15
                ReadRFCM_TestBed("data/TestBed_RFCM/R4C16.csv", settings1); // read and save RFCM gain for ch16
            }

            // read gain offset for chs file!!
            ReadGainOffset_TestBed("./data/preamp_ch_gain_offset.csv", settings1);// only TestBed for now
            // read threshold offset for chs file!!
            ReadThresOffset_TestBed("./data/threshold_offset.csv", settings1);// only TestBed for now
	    // read threshold values for chs file
	    ReadThres_TestBed("./data/thresholds_TB.csv", settings1);// only TestBed for now
            // read system temperature for chs file!!
	cout << "check read temp testbed 3" << endl;
            if (settings1->NOISE_CHANNEL_MODE!=0) {
                ReadTemp_TestBed("./data/system_temperature.csv", settings1);// only TestBed for now
            }
	    
            // read total elec. chain response file!!
	    cout<<"start read elect chain"<<endl;
        if(settings1->CUSTOM_ELECTRONICS==0){
            //read the standard ARA electronics
            cout<<"     Reading standard ARA electronics response"<<endl;
             ReadElectChain("./data/ARA_Electronics_TotalGain_TwoFilters.txt", settings1);
            //ReadElectChain("./data/ARA_Electronics_TotalGainPhase.txt", settings1);
        }
        else if (settings1->CUSTOM_ELECTRONICS==1){
            //read a custom user defined electronics gain
            cout<<"     Reading custom electronics response"<<endl;
             ReadElectChain("./data/custom_electronics.txt", settings1);
        }
	    cout<<"done read elect chain"<<endl;
	    

	    
	    // if calpulser case
            if (settings1->CALPULSER_ON > 0) {
	      // read TestBed Calpulser waveform measured (before pulser)
	      ReadCalPulserWF("./data/CalPulserWF.txt", settings1);
            }
	    
	    

    }// if mode == 3

    else if (mode == 4) {        
      //        cout<<"\n\tDector mode 4 : Single installed station determined by DETECTOR_STATION !"<<endl;
      //        cout<<"We use "<<ARA_N_file.c_str()<<" as antenna info."<<endl;
        
        //SetupInstalledStations();        
        
        // initialize info
        params.number_of_stations = 1; //including Testbed
        params.number_of_strings_station = 4;   // ARA-1 has 4 strings
        params.number_of_antennas_string = 4; // 4 antennas on each strings
        params.number_of_surfaces_station = 4;
        params.number_of_channels = 20;
        
        //double core_x = 0.;
        //double core_y = 0.;
        params.core_x = 10000.;
        params.core_y = 10000.;
        double R_string = 10.;  // all units are in meter
        double R_surface = 60.;
        double z_max = 200.;
        double z_btw = 20.;
        params.stations_per_side = 4;       // total 37 stations
        params.station_spacing = 2000.;     // 2km spacing for borehole stations
        params.antenna_orientation = 0;     // all antenna facing x
        params.bore_hole_antenna_layout = settings1->BORE_HOLE_ANTENNA_LAYOUT;
        // finish initialization
        //

        // mode == 4 currently just use installed TestBed station geom information.
        // So don't need to read any more information
        
        // Read new parameters if there are...
        ifstream ARA_N( ARA_N_file.c_str() );
        if ( ARA_N.is_open() ) {
            while (ARA_N.good() ) {
                getline (ARA_N, line);
                
                if (line[0] != "/"[0]) {
                    label = line.substr(0, line.find_first_of("=") );
                    
                    if (label == "core_x") {
                        params.core_x = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_x"<<endl;
                    }
                    else if (label == "core_y") {
                        params.core_y = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read core_y"<<endl;
                    }
                    else if (label == "R_string") {
                        R_string = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_string"<<endl;
                    }
                    else if (label == "R_surface") {
                        R_surface = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read R_surface"<<endl;
                    }
                    else if (label == "z_max") {
                        z_max = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_max"<<endl;
                    }
                    else if (label == "z_btw") {
                        z_btw = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read z_btw"<<endl;
                    }
                    else if (label == "number_of_stations") {
                        params.number_of_stations = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read stations_per_side"<<endl;
                    }
                    else if (label == "station_spacing") {
                        params.station_spacing = atof( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read station_spacting"<<endl;
                    }
                    else if (label == "antenna_orientation") {
                        params.antenna_orientation = atoi( line.substr( line.find_first_of("=") + 1).c_str() );
                        cout<<"read antenna_orientation"<<endl;
                    }
                }
            }
            ARA_N.close();
        }
        // finished reading new parameters
        
        params.number_of_antennas_string = 4;
        
        // prepare vectors
        PrepareVectorsInstalled(settings1->DETECTOR_STATION);
        // end prepare vectors
        
        //
        // for ARA-37 (or more than 1 station case), need code for setting position for all 37 stations here!
        //
        int station_count = 0;


	stations[0].SetX( params.core_x );
	stations[0].SetY( params.core_y );

        
        //        cout<<"total station_count : "<<station_count<<endl;
        if (station_count != (int)params.number_of_stations) cout<<"\n\tError, station number not match !"<<endl;
        
        //
        // set antenna values from parameters
        // set station positions
        //cout << "READGEOM:" << settings1->READGEOM << endl;
        
#ifdef ARA_UTIL_EXISTS
        ImportStationInfo(settings1, 0, settings1->DETECTOR_STATION);
#endif
//            UseAntennaInfo(1, settings1);
	int stationID = settings1->DETECTOR_STATION;
        for (int i = 0; i < (int)params.number_of_stations; i++){
	  stations[i].StationID = settings1->DETECTOR_STATION;
	  if (settings1->USE_INSTALLED_TRIGGER_SETTINGS == 0){
	    stations[i].NFOUR = 1024;
	    stations[i].TIMESTEP = 1./2.6*1.E-9;
	    stations[i].TRIG_WINDOW = 2.5E-7;
	    stations[i].DATA_BIN_SIZE = settings1->DATA_BIN_SIZE;
	  }
	  else if (settings1->USE_INSTALLED_TRIGGER_SETTINGS == 1){
	    if (stations[i].StationID == 0){
	      stations[i].NFOUR = 1024;
	      stations[i].TIMESTEP = 1./2.6*1.E-9;
	      stations[i].TRIG_WINDOW = 2.5E-7;
	      stations[i].DATA_BIN_SIZE = settings1->DATA_BIN_SIZE;
	    }
	    if (stations[i].StationID == 1){
	      stations[i].NFOUR = 1024;
	      stations[i].TIMESTEP = 1./2.6*1.E-9;
	      stations[i].TRIG_WINDOW = 2.5E-7;
	      stations[i].DATA_BIN_SIZE = settings1->DATA_BIN_SIZE;
	    }
	  }
        }
        
        params.number_of_antennas = 0;
	
	cout<<"DETECTOR=4 imported station geom info"<<endl;
        
	for (int j = 0; j < stations[0].strings.size(); j++){
	  for (int k = 0; k < stations[0].strings[j].antennas.size(); k++){
	    
	    cout <<
	      "DetectorStation2:string:antenna:X:Y:Z:chno :: " <<
	      j<< " : " <<
	      k<< " : " <<
	      stations[0].strings[j].antennas[k].GetX() << " : " <<
	      stations[0].strings[j].antennas[k].GetY() << " : " <<
	      stations[0].strings[j].antennas[k].GetZ() << " : \t" <<
	      //GetChannelfromStringAntenna ( 0, j, k)<<
	      GetChannelfromStringAntenna ( stationID, j, k, settings1)<<
	      endl;
	    
	    params.number_of_antennas++;
	  }
	}
        



	cout<<"after FlattoEarth, station0 location"<<endl;
	for (int j = 0; j < stations[0].strings.size(); j++){
	  for (int k = 0; k < stations[0].strings[j].antennas.size(); k++){
	    
	    cout <<
             "Detector:station:string:antenna:X:Y:Z:R:Theta:Phi:: " <<
             "0" << " : " <<
             j<< " : " <<
             k<< " : " <<
             stations[0].strings[j].antennas[k].GetX() << " : " <<
             stations[0].strings[j].antennas[k].GetY() << " : " <<
             stations[0].strings[j].antennas[k].GetZ() << " : " <<
             stations[0].strings[j].antennas[k].R() << " : " <<
             stations[0].strings[j].antennas[k].Theta() << " : " <<
             stations[0].strings[j].antennas[k].Phi() << " : " <<
             icesurface->Surface(stations[0].strings[j].antennas[k].Lon(), stations[0].strings[j].antennas[k].Lat()) << " : " <<
//             icesurface->Surface(stations[0].strings[j].antennas[k].Lat(), stations[0].strings[j].antennas[k].Lon()) << " : " <<
             endl;
                 
                     
        }
	}



            int antenna_count = 0;
            max_number_of_antennas_station = 0;
            // for non-idealized geometry, it's better to actually count number of stations
            for (int i=0; i<(int)(stations.size()); i++) {
            
                antenna_count = 0;
                for (int j=0; j<(int)(stations[i].strings.size()); j++) {
                    for (int k=0; k<(int)(stations[i].strings[j].antennas.size()); k++) {
                        antenna_count++;
                    }
                }
                stations[i].number_of_antennas = antenna_count;

                if (max_number_of_antennas_station < antenna_count) max_number_of_antennas_station = antenna_count;
            }


            
            


	    if (settings1->ANTENNA_MODE == 0){
	      // test read V-pol gain file!!
	      ReadVgain("ARA_bicone6in_output.txt");
	      // test read H-pol gain file!!
	      ReadHgain("ARA_dipoletest1_output.txt");
	    }
	    else if (settings1->ANTENNA_MODE == 1) {
	      // test read V-pol gain file!!
	      ReadVgainSettings("ARA_bicone6in_output_updated2016.txt", settings1);
	      ReadVgainTopSettings("ARA_VPresult_topTrec.txt", settings1);
	      // test read H-pol gain file!!
	      ReadHgainSettings("ARA_dipoletest1_output_updated2016.txt", settings1);
	    }
	    else if (settings1->ANTENNA_MODE == 2){
	      // test read V-pol gain file!!
	      ReadVgain("Arianna_WIPLD_hpol.dat");
	      // test read H-pol gain file!!
	      ReadHgain("Arianna_WIPLD_hpol.dat");
	    }

	    //	    if (settings1->NOISE == 2){
	      //Read the noise figures
	      ReadNoiseFigure("./data/ARA02_noiseFig.txt", settings1);
	      //	    }

            // read filter file!!
            ReadFilter("./data/filter.csv", settings1);
            // read preamp gain file!!
            ReadPreamp("./data/preamp.csv", settings1);
            // read FOAM gain file!!
            ReadFOAM("./data/FOAM.csv", settings1);

	    if (settings1->NOISE_CHANNEL_MODE!=0) {
	      ReadTemp_TestBed("./data/system_temperature.csv", settings1);// only TestBed for now
	    }
	    


	    if (settings1->DETECTOR_STATION == 0){
	      if ( settings1->NOISE==1) {
                // read Rayleigh fit for freq range, bh channels
                //ReadRayleighFit_TestBed("data/RayleighFit_TB.csv", settings1, settings1->RAYLEIGH_CH); // read and save RFCM gain
                ReadRayleighFit_TestBed("data/RayleighFit_TB.csv", settings1); // read and save RFCM gain
	      }
	      
	      if ( settings1->USE_TESTBED_RFCM_ON==1) {
                // read RFCM gain file!! (measured value in ICL)
                ReadRFCM_TestBed("data/TestBed_RFCM/R1C1.csv", settings1); // read and save RFCM gain for ch1
                ReadRFCM_TestBed("data/TestBed_RFCM/R1C2.csv", settings1); // read and save RFCM gain for ch2
                ReadRFCM_TestBed("data/TestBed_RFCM/R1C3.csv", settings1); // read and save RFCM gain for ch3
                ReadRFCM_TestBed("data/TestBed_RFCM/R1C4.csv", settings1); // read and save RFCM gain for ch4
                ReadRFCM_TestBed("data/TestBed_RFCM/R2C5.csv", settings1); // read and save RFCM gain for ch5
                ReadRFCM_TestBed("data/TestBed_RFCM/R2C6.csv", settings1); // read and save RFCM gain for ch6
                ReadRFCM_TestBed("data/TestBed_RFCM/R2C7.csv", settings1); // read and save RFCM gain for ch7
                ReadRFCM_TestBed("data/TestBed_RFCM/R2C8.csv", settings1); // read and save RFCM gain for ch8
                ReadRFCM_TestBed("data/TestBed_RFCM/R3C9.csv", settings1); // read and save RFCM gain for ch9
                ReadRFCM_TestBed("data/TestBed_RFCM/R3C10.csv", settings1); // read and save RFCM gain for ch10
                ReadRFCM_TestBed("data/TestBed_RFCM/R3C11.csv", settings1); // read and save RFCM gain for ch11
                ReadRFCM_TestBed("data/TestBed_RFCM/R3C12.csv", settings1); // read and save RFCM gain for ch12
                ReadRFCM_TestBed("data/TestBed_RFCM/R4C13.csv", settings1); // read and save RFCM gain for ch13
                ReadRFCM_TestBed("data/TestBed_RFCM/R4C14.csv", settings1); // read and save RFCM gain for ch14
                ReadRFCM_TestBed("data/TestBed_RFCM/R4C15.csv", settings1); // read and save RFCM gain for ch15
                ReadRFCM_TestBed("data/TestBed_RFCM/R4C16.csv", settings1); // read and save RFCM gain for ch16
	      }

	      // read gain offset for chs file!!
	      ReadGainOffset_TestBed("./data/preamp_ch_gain_offset.csv", settings1);// only TestBed for now
	      // read threshold offset for chs file!!
	      ReadThresOffset_TestBed("./data/threshold_offset.csv", settings1);// only TestBed for now
	      // read threshold values for chs file
	      ReadThres_TestBed("./data/thresholds_TB.csv", settings1);// only TestBed for now
	      // read system temperature for chs file!!
	      cout << "check read temp testbed 4" << endl;
	      if (settings1->NOISE_CHANNEL_MODE!=0) {
		ReadTemp_TestBed("./data/system_temperature.csv", settings1);// only TestBed for now
	      }
	    }
	    /*
        if (settings1->DETECTOR_STATION == 2){
            if ( settings1->NOISE==1) {
                // read Rayleigh fit for freq range, bh channels
                //ReadRayleighFit_TestBed("data/RayleighFit_TB.csv", settings1, settings1->RAYLEIGH_CH); // read and save RFCM gain
                ReadRayleighFit_TestBed("data/RayleighFit_TB.csv", settings1); // read and save RFCM gain
            }
            
            /*
            // read gain offset for chs file!!
            ReadGainOffset_TestBed("./data/preamp_ch_gain_offset.csv", settings1);// only TestBed for now
            // read threshold offset for chs file!!
            ReadThresOffset_TestBed("./data/threshold_offset.csv", settings1);// only TestBed for now
            // read threshold values for chs file
            ReadThres_TestBed("./data/thresholds_TB.csv", settings1);// only TestBed for now
            // read system temperature for chs file!!
            
            if (settings1->NOISE_CHANNEL_MODE!=0) {
                ("./data/system_temperature.csv", settings1);// only TestBed for now
            }
             */
	    //        }
	    /*
        if (settings1->DETECTOR_STATION == 3){
            if ( settings1->NOISE==1) {
                // read Rayleigh fit for freq range, bh channels
                //ReadRayleighFit_TestBed("data/RayleighFit_TB.csv", settings1, settings1->RAYLEIGH_CH); // read and save RFCM gain
                ReadRayleighFit_TestBed("data/RayleighFit_TB.csv", settings1); // read and save RFCM gain
            }
            /*
            // read gain offset for chs file!!
            ReadGainOffset_TestBed("./data/preamp_ch_gain_offset.csv", settings1);// only TestBed for now
            // read threshold offset for chs file!!
            ReadThresOffset_TestBed("./data/threshold_offset.csv", settings1);// only TestBed for now
            // read threshold values for chs file
            ReadThres_TestBed("./data/thresholds_TB.csv", settings1);// only TestBed for now
            // read system temperature for chs file!!
            if (settings1->NOISE_CHANNEL_MODE!=0) {
                ("./data/system_temperature.csv", settings1);// only TestBed for now
            }
             */
	    //        }
        
        
            // read total elec. chain response file!!
	    cout<<"start read elect chain"<<endl;
        if(settings1->CUSTOM_ELECTRONICS==0){
            //read the standard ARA electronics
            cout<<"     Reading standard ARA electronics response"<<endl;
             ReadElectChain("./data/ARA_Electronics_TotalGain_TwoFilters.txt", settings1);
            //ReadElectChain("./data/ARA_Electronics_TotalGainPhase.txt", settings1);
        }
        else if (settings1->CUSTOM_ELECTRONICS==1){
            //read a custom user defined electronics gain
            cout<<"     Reading custom electronics response"<<endl;
             ReadElectChain("./data/custom_electronics.txt", settings1);
        }
	    cout<<"done read elect chain"<<endl;
    
	    // if calpulser case
            if (settings1->CALPULSER_ON > 0) {
	      // read TestBed Calpulser waveform measured (before pulser)
	      ReadCalPulserWF("./data/CalPulserWF.txt", settings1);
            }
	    


    }// if mode == 4


    /////////////////////////////////////////////////////////////////////////////////    



            
        

    // add additional depth if it's on
    AddAdditional_Depth(settings1);

    
            
    // change coordinate from flat surface to curved Earth surface
    //FlattoEarth_ARA(icesurface);
    FlattoEarth_ARA_sharesurface(icesurface);   // this one will share the lowest surface at each station.

    
    
    //cout<<"done settings detectors, gain, filters"<<endl;
    
    getDiodeModel(settings1);    // set diode_real and fdiode_real values.
    
    
    //    return 0;
    
    //cout<<"test2"<<endl;
}




inline void Detector::ReadVgain(string filename) {
    ifstream NecOut( filename.c_str() );
    
    string line;
    if ( NecOut.is_open() ) {
        while (NecOut.good() ) {

            for (int i=0; i<freq_step; i++) {
                getline (NecOut, line);
                if ( line.substr(0, line.find_first_of(":")) == "freq ") {
		  Freq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
		  cout<<"freq["<<i<<"] = "<<Freq[i]<<" MHz"<<endl;
		  getline (NecOut, line); //read SWR
		  
		  getline (NecOut, line); //read names
		  
		  for (int j=0; j<ang_step; j++) {
		    getline (NecOut, line); //read data line
		    //Vgain[i][j] = atof( line.substr( 18 ).c_str() );  // read gain (not dB)
		    Vgain[i][j] = atof( line.substr( 20, 33 ).c_str() );  // read gain (not dB)
		    Vphase[i][j] = atof( line.substr( 34 ).c_str() );  // read gain (not dB)
		    //	cout << atof( line.substr( 20, 33 ).c_str() ) << endl;
		    
		    //	cout<<"VGain : "<<Vgain[i][j]<<", VPhase : "<<Vphase[i][j]<<endl;
                        
                    }// end ang_step
                    
                }// end check freq label
                
            }// end freq_step
            
        }// end while NecOut.good
        NecOut.close();
    }// end if file open

}// end ReadVgain

inline void Detector::ReadVgainSettings(string filename, Settings *settings1) {
    ifstream NecOut( filename.c_str() );
    const int N = freq_step;
    double Transm[N];
    
    string line;
    if ( NecOut.is_open() ) {
        while (NecOut.good() ) {

            for (int i=0; i<freq_step; i++) {
                getline (NecOut, line);
                if ( line.substr(0, line.find_first_of(":")) == "freq ") {
                    Freq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
                    //                    cout<<"freq["<<i<<"] = "<<Freq[i]<<" MHz"<<endl;
                    getline (NecOut, line); //read SWR
           		     Transm[i] = atof(line.substr(5,11).c_str()); //What says "SWR" in "ARA_bicone6in_output_updated2016.txt" is actually trasnmission coefficient - MYL 01/23/20
                    getline (NecOut, line); //read names

                    for (int j=0; j<ang_step; j++) {
                        getline (NecOut, line); //read data line
                        //Vgain[i][j] = atof( line.substr( 18 ).c_str() );  // read gain (not dB)
                        Vgain[i][j] = /*Transm[i] **/ atof( line.substr( 20, 33 ).c_str() );  // read gain (not dB) //Transm[i] commented for in the antenna gain files "realized gain" is provided and so the transmission coefficient does not need to be multiplied one more time - MLY 01/23/20
                        Vphase[i][j] = atof( line.substr( 34 ).c_str() );  // read gain (not dB)
                                                
                        //cout<<"VGain : "<<Vgain[i][j]<<", VPhase : "<<Vphase[i][j]<<endl;
                        
                    }// end ang_step
                    
                }// end check freq label
                
            }// end freq_step
            
        }// end while NecOut.good
        NecOut.close();
    }// end if file open
    double xfreq[N];
    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double trans_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    // now below are values that shared in all channels
    for (int i=0;i<freq_step;i++) { // copy values
      xfreq[i] = Freq[i];
      
    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
      xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    
    
    Tools::SimpleLinearInterpolation( freq_step-1, xfreq, Transm, settings1->DATA_BIN_SIZE/2, xfreq_databin, trans_databin );
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
      transV_databin.push_back(trans_databin[i]); // from Hz to MHz
    }
    
}// end ReadVgain


inline void Detector::ReadVgainTopSettings(string filename, Settings *settings1) {
    ifstream NecOut( filename.c_str() );
    const int N = freq_step;

    double Transm[N]; 
    string line;
    if ( NecOut.is_open() ) {
        while (NecOut.good() ) {

            for (int i=0; i<freq_step; i++) {
                getline (NecOut, line);
                if ( line.substr(0, line.find_first_of(":")) == "freq ") {
                    Freq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
                    //                    cout<<"freq["<<i<<"] = "<<Freq[i]<<" MHz"<<endl;
                    getline (NecOut, line); //read SWR
		              Transm[i] = atof( line.substr(5, 11).c_str() ); //What say "SWR" in "ARA_VPresult_topTrec.txt" is actually transmission coefficient - MYL 01/23/20
		              //cerr << "Vpol Transm: " << Transm[i] << endl;

                    getline (NecOut, line); //read names

                    for (int j=0; j<ang_step; j++) {
                        getline (NecOut, line); //read data line
                        VgainTop[i][j] = /*Transm[i]**/atof( line.substr( 20, 33 ).c_str() ) ; //Transm[i] commented for in the antenna gain files "realized gain" is provided and so the transmission coefficient does not need to be multiplied one more time - MLY 01/23/20
                        VphaseTop[i][j] = atof( line.substr( 34 ).c_str() );    // + 180.0/TMath::Pi()*TMath::ATan(-Freq[i]/500.0);  // read gain (not dB)
                                                
                        //cout<<"VGain : "<<Vgain[i][j]<<", VPhase : "<<Vphase[i][j]<<endl;
                        
                    }// end ang_step
                    
                }// end check freq label
                
            }// end freq_step
            
        }// end while NecOut.good
        NecOut.close();
    }// end if file open
    double xfreq[N];
    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double trans_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    // now below are values that shared in all channels
    for (int i=0;i<freq_step;i++) { // copy values
        xfreq[i] = Freq[i];

    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    

        Tools::SimpleLinearInterpolation( freq_step-1, xfreq, Transm, settings1->DATA_BIN_SIZE/2, xfreq_databin, trans_databin );
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        transVTop_databin.push_back(trans_databin[i]); // from Hz to MHz
    }



}// end ReadVgainTop


inline void Detector::ReadHgain(string filename) {
  ifstream NecOut( filename.c_str() );
    
    string line;
    
    if ( NecOut.is_open() ) {
      //  cout << "OK"  <<endl;
        while (NecOut.good() ) {
            
            for (int i=0; i<freq_step; i++) {
                getline (NecOut, line);
                if ( line.substr(0, line.find_first_of(":")) == "freq ") {
                    Freq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
		    //	    cout<<"freq["<<i<<"] = "<<Freq[i]<<" MHz"<<endl;
		    getline (NecOut, line); //read SWR
                    getline (NecOut, line); //read names
                    
                    for (int j=0; j<ang_step; j++) {
                        getline (NecOut, line); //read data line
                        //Hgain[i][j] = atof( line.substr( 20 ).c_str() );  // read gain (not dB)
                        //Hgain[i][j] = atof( line.substr( 18, 25 ).c_str() );  // read gain (not dB)
			Hgain[i][j] = atof( line.substr( 20, 33 ).c_str() );  // read gain (not dB)
			Hphase[i][j] = atof( line.substr( 34 ).c_str() );  // read gain (not dB)

			//      cout<<"HGain : "<<Hgain[i][j]<<", HPhase : "<<Hphase[i][j]<<endl;
                        
                    }// end ang_step
                    
                }// end check freq label
                
            }// end freq_step
            
        }// end while NecOut.good
        NecOut.close();
    }// end if file open
    
}// end ReadHgain

inline void Detector::ReadHgainSettings(string filename, Settings *settings1) {

    ifstream NecOut( filename.c_str() );    
    string line;
    
    const int N = freq_step;
    //string line;
    double Transm[N]; 
    if ( NecOut.is_open() ) {
        while (NecOut.good() ) {
            
            for (int i=0; i<freq_step; i++) {
                getline (NecOut, line);
                if ( line.substr(0, line.find_first_of(":")) == "freq ") {
                    Freq[i] = atof( line.substr(6, line.find_first_of("M")).c_str() );
                    //                    cout<<"freq["<<i<<"] = "<<Freq[i]<<" MHz"<<endl;
                    getline (NecOut, line); //read SWR
            	     Transm[i] = atof( line.substr(5, 11).c_str() ); //What says "SWR" in "ARA_dipoletest1_output_updated2016.txt" is actually transmission coefficient - MYL 01/23/20 
                    getline (NecOut, line); //read names
                    
                    for (int j=0; j<ang_step; j++) {
                        getline (NecOut, line); //read data line
                        //Hgain[i][j] = atof( line.substr( 20 ).c_str() );  // read gain (not dB)
                        //Hgain[i][j] = atof( line.substr( 18, 25 ).c_str() );  // read gain (not dB)
                        Hgain[i][j] = /*Transm[i]**/atof( line.substr( 20, 33 ).c_str() );  // read gain (not dB) //Transm[i] commented for in the antenna gain files "realized gain" is provided and so the transmission coefficient does not need to be multiplied one more time - MLY 01/23/20                      
                        Hphase[i][j] = atof( line.substr( 34 ).c_str() );  // read gain (not dB)

                        //cout<<"HGain : "<<Hgain[i][j]<<", HPhase : "<<Hphase[i][j]<<endl;
                        
                    }// end ang_step
                    
                }// end check freq label
                
            }// end freq_step
            
        }// end while NecOut.good
         NecOut.close();
    }// end if file open
    double xfreq[N];
    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double trans_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    // now below are values that shared in all channels
    for (int i=0;i<freq_step;i++) { // copy values
        xfreq[i] = Freq[i];

    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    

        Tools::SimpleLinearInterpolation( freq_step-1, xfreq, Transm, settings1->DATA_BIN_SIZE/2, xfreq_databin, trans_databin );
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        transH_databin.push_back(trans_databin[i]); // from Hz to MHz
    }
 
}// end ReadHgain


double Detector::GetGain(double freq, double theta, double phi, int ant_m, int ant_o) { // using Interpolation on multidimentions!
    //double GetGain(double freq, double theta, double phi, int ant_m, int ant_o) { // using Interpolation on multidimentions!
    
    //Parameters params;
    
    // change antenna facing orientation
    if (ant_o == 0) {
        // no change...
    }
    else if (ant_o == 1) {
        if (phi - 90. >= 0.) {
            phi = phi - 90.;
        }
        else {
            phi = 360. + phi - 90.;
        }
    }
    else if (ant_o == 2) {
        if (phi - 180. >= 0.) {
            phi = phi - 180.;
        }
        else {
            phi = 360. + phi - 180.;
        }
    }
    else if (ant_o == 3) {
        if (phi - 270. >= 0.) {
            phi = phi - 270.;
        }
        else {
            phi = 360. + phi - 270.;
        }
    }
    else {
        cout<<"Wrong option selected for antenna orientation "<<ant_o<<" !!"<<endl;
        cout<<"ant_o will be replaced from "<<ant_o<<" to 0"<<endl;
    }
    // end changing antenna orientation
    
    
    int i = (int)(theta/5.);
    int j = (int)(phi/5.);
    
    double thetai = 5.*( (int)(theta/5.) );
    double thetai1 = 5.*( (int)(theta/5.) + 1.);
    double phij = 5.*( (int)(phi/5.) );
    double phij1 = 5.*( (int)(phi/5.) + 1.);
    
    double t = (theta - thetai)/(thetai1 - thetai);
    double u = (phi - phij)/(phij1 - phij);
    
    // in case when freq is out of nec2 freq range. use nearest min/max freq bin value. 
    if ( freq < freq_init ) {
        //cout<<"Frequency value is smaller than frequency range with Gain."<<endl;
        //cout<<"Frequency value "<<freq<<" will be replaced to minimum frequency value "<<freq_init<<endl;
        freq = freq_init;
    }
    else if ( freq > (freq_init + freq_width*((double)freq_step-1.) ) ) {
        //cout<<"Frequency value is bigger than frequency range with Gain."<<endl;
        //cout<<"Frequency value "<<freq<<" will be replaced to maximum frequency value "<< freq_init + freq_width*((double)freq_step-1.) - 0.01 <<endl;
        freq = freq_init + freq_width*((double)freq_step-1.) - 0.01;
    }
    
    
    //    int fx1 = (int)( (freq + (freq_width/2.) - freq_init)/freq_width );
    int fx1 = (int)( (freq - freq_init)/freq_width );
    int fx2 = fx1 + 1;
    //    cout<<"fx1 : "<<fx1<<endl;
    //    cout<<"fx2 : "<<fx2<<endl;
    
    double Gij, Gi1j, Gij1, Gi1j1, Gout1, Gout2, Gout;
    
    if (ant_m == 0) {   // for V pol antenna!!
        Gij = Vgain[fx1][(int)(37*j+i)];
        Gi1j = Vgain[fx1][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx1][(int)(i)];
            Gi1j1 = Vgain[fx1][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx1][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx1][(int)(37*(j+1)+i+1)];
        }
        
        Gout1 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest smaller freq bin
        
        Gij = Vgain[fx2][(int)(37*j+i)];
        Gi1j = Vgain[fx2][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx2][(int)(i)];
            Gi1j1 = Vgain[fx2][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx2][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx2][(int)(37*(j+1)+i+1)];
        }
        
        Gout2 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest higher freq bin
        
    }
    
    else if (ant_m == 1) {   // for H pol antenna!!
        Gij = Hgain[fx1][(int)(37*j+i)];
        Gi1j = Hgain[fx1][(int)(37*j+i+1)];
        if ( j == 71 ) {
            Gij1 = Hgain[fx1][(int)(i)];
            Gi1j1 = Hgain[fx1][(int)(i+1)];
        }
        else {
            Gij1 = Hgain[fx1][(int)(37*(j+1)+i)];
            Gi1j1 = Hgain[fx1][(int)(37*(j+1)+i+1)];
        }
        
        Gout1 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest smaller freq bin
        
        Gij = Vgain[fx2][(int)(37*j+i)];
        Gi1j = Vgain[fx2][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx2][(int)(i)];
            Gi1j1 = Vgain[fx2][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx2][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx2][(int)(37*(j+1)+i+1)];
        }
        
        Gout2 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest higher freq bin
    }
    
    else {
        cout<<"There is no antenna type : "<<ant_m<<" !!"<<endl;
        cout<<"Will return Gain = 0 !!"<<endl;
        Gout1 = 0.;
        Gout2 = 0.;
    }
    
    Gout = ((Gout2 - Gout1)/freq_width) * ( freq - (freq_init + fx1*freq_width) ) + Gout1; // get linear interpolation between two nearest freq bin.
    
    //  cout<<Gout<<endl;
    
    return Gout;
    
    // ant_o face x = 0, y = 1, -x = 2, -y = 3
    
}


double Detector::GetGain(double freq, double theta, double phi, int ant_m) {
    //double GetGain(double freq, double theta, double phi, int ant_m) {
    
    //Parameters params;
    
    int i = (int)(theta/5.);
    int j = (int)(phi/5.);
    
    double thetai = 5.*( (int)(theta/5.) );
    double thetai1 = 5.*( (int)(theta/5.) + 1.);
    double phij = 5.*( (int)(phi/5.) );
    double phij1 = 5.*( (int)(phi/5.) + 1.);
    
    double t = (theta - thetai)/(thetai1 - thetai);
    double u = (phi - phij)/(phij1 - phij);
    
    
    // in case when freq is out of nec2 freq range. use nearest min/max freq bin value. 
    if ( freq < freq_init ) {
        //cout<<"Frequency value is smaller than frequency range with Gain."<<endl;
        //cout<<"Frequency value "<<freq<<" will be replaced to minimum frequency value "<<freq_init<<endl;
        freq = freq_init;
    }
    else if ( freq > (freq_init + freq_width*((double)freq_step - 1.) ) ) {
        //cout<<"Frequency value is bigger than frequency range with Gain."<<endl;
        //cout<<"Frequency value "<<freq<<" will be replaced to maximum frequency value "<< freq_init + freq_width*((double)freq_step-1.) - 0.01 <<endl;
        freq = freq_init + freq_width*((double)freq_step-1.) - 0.01;
    }
    
    
    //    int fx1 = (int)( (freq + (freq_width/2.) - freq_init)/freq_width );
    int fx1 = (int)( (freq - freq_init)/freq_width );
    int fx2 = fx1 + 1;
    
    double Gij, Gi1j, Gij1, Gi1j1, Gout1, Gout2, Gout;
    
    if (ant_m == 0) {   // for V pol antenna!!
        Gij = Vgain[fx1][(int)(37*j+i)];
        Gi1j = Vgain[fx1][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx1][(int)(i)];
            Gi1j1 = Vgain[fx1][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx1][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx1][(int)(37*(j+1)+i+1)];
        }
        
        Gout1 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest smaller freq bin
        
        Gij = Vgain[fx2][(int)(37*j+i)];
        Gi1j = Vgain[fx2][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vgain[fx2][(int)(i)];
            Gi1j1 = Vgain[fx2][(int)(i+1)];
        }
        else {
            Gij1 = Vgain[fx2][(int)(37*(j+1)+i)];
            Gi1j1 = Vgain[fx2][(int)(37*(j+1)+i+1)];
        }
        
        Gout2 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest higher freq bin
        
    }
    
    else if (ant_m == 1) {   // for H pol antenna!!
        Gij = Hgain[fx1][(int)(37*j+i)];
        Gi1j = Hgain[fx1][(int)(37*j+i+1)];
        if ( j == 71 ) {
            Gij1 = Hgain[fx1][(int)(i)];
            Gi1j1 = Hgain[fx1][(int)(i+1)];
        }
        else {
            Gij1 = Hgain[fx1][(int)(37*(j+1)+i)];
            Gi1j1 = Hgain[fx1][(int)(37*(j+1)+i+1)];
        }
        
        Gout1 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest smaller freq bin
        
        Gij = Hgain[fx2][(int)(37*j+i)];
        Gi1j = Hgain[fx2][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Hgain[fx2][(int)(i)];
            Gi1j1 = Hgain[fx2][(int)(i+1)];
        }
        else {
            Gij1 = Hgain[fx2][(int)(37*(j+1)+i)];
            Gi1j1 = Hgain[fx2][(int)(37*(j+1)+i+1)];
        }
        
        Gout2 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //Gain at nearest higher freq bin
    }
    
    else {
        cout<<"There is no antenna type : "<<ant_m<<" !!"<<endl;
        cout<<"Will return Gain = 0 !!"<<endl;
        Gout1 = 0.;
        Gout2 = 0.;
    }
    
    Gout = ((Gout2 - Gout1)/freq_width) * ( freq - (freq_init + fx1*freq_width) ) + Gout1; // get linear interpolation between two nearest freq bin.
    
    cout<<Gout<<endl;
    return Gout;
    
    
}



double Detector::GetAntPhase( double freq, double theta, double phi, int ant_m ) {



    int i = (int)(theta/5.);
    int j = (int)(phi/5.);
    
    double thetai = 5.*( (int)(theta/5.) );
    double thetai1 = 5.*( (int)(theta/5.) + 1.);
    double phij = 5.*( (int)(phi/5.) );
    double phij1 = 5.*( (int)(phi/5.) + 1.);
    
    double t = (theta - thetai)/(thetai1 - thetai);
    double u = (phi - phij)/(phij1 - phij);
    
    
    // in case when freq is out of nec2 freq range. use nearest min/max freq bin value. 
    if ( freq < freq_init ) {
        //cout<<"Frequency value is smaller than frequency range with phase."<<endl;
        //cout<<"Frequency value "<<freq<<" will be replaced to minimum frequency value "<<freq_init<<endl;
        freq = freq_init;
    }
    else if ( freq > (freq_init + freq_width*((double)freq_step - 1.) ) ) {
        //cout<<"Frequency value is bigger than frequency range with phase."<<endl;
        //cout<<"Frequency value "<<freq<<" will be replaced to maximum frequency value "<< freq_init + freq_width*((double)freq_step-1.) - 0.01 <<endl;
        freq = freq_init + freq_width*((double)freq_step-1.) - 0.01;
    }
    
    
    //    int fx1 = (int)( (freq + (freq_width/2.) - freq_init)/freq_width );
    int fx1 = (int)( (freq - freq_init)/freq_width );
    int fx2 = fx1 + 1;
    
    double Gij, Gi1j, Gij1, Gi1j1, Gout1, Gout2, Gout;
    
    if (ant_m == 0) {   // for V pol antenna!!

        Gij = Vphase[fx1][(int)(37*j+i)];
        Gi1j = Vphase[fx1][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vphase[fx1][(int)(i)];
            Gi1j1 = Vphase[fx1][(int)(i+1)];
        }
        else {
            Gij1 = Vphase[fx1][(int)(37*(j+1)+i)];
            Gi1j1 = Vphase[fx1][(int)(37*(j+1)+i+1)];
        }
        
        Gout1 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //phase at nearest smaller freq bin
        
        Gij = Vphase[fx2][(int)(37*j+i)];
        Gi1j = Vphase[fx2][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Vphase[fx2][(int)(i)];
            Gi1j1 = Vphase[fx2][(int)(i+1)];
        }
        else {
            Gij1 = Vphase[fx2][(int)(37*(j+1)+i)];
            Gi1j1 = Vphase[fx2][(int)(37*(j+1)+i+1)];
        }
        
        Gout2 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //phase at nearest higher freq bin
        
    }
    
    else if (ant_m == 1) {   // for H pol antenna!!
        Gij = Hphase[fx1][(int)(37*j+i)];
        Gi1j = Hphase[fx1][(int)(37*j+i+1)];
        if ( j == 71 ) {
            Gij1 = Hphase[fx1][(int)(i)];
            Gi1j1 = Hphase[fx1][(int)(i+1)];
        }
        else {
            Gij1 = Hphase[fx1][(int)(37*(j+1)+i)];
            Gi1j1 = Hphase[fx1][(int)(37*(j+1)+i+1)];
        }
        
        Gout1 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //phase at nearest smaller freq bin
        
        Gij = Hphase[fx2][(int)(37*j+i)];
        Gi1j = Hphase[fx2][(int)(37*j+i+1)];
        if ( j == 71 ) {    // doing this as maximum phi is 355 deg
            Gij1 = Hphase[fx2][(int)(i)];
            Gi1j1 = Hphase[fx2][(int)(i+1)];
        }
        else {
            Gij1 = Hphase[fx2][(int)(37*(j+1)+i)];
            Gi1j1 = Hphase[fx2][(int)(37*(j+1)+i+1)];
        }
        
        Gout2 = (1.-t)*(1.-u)*Gij + t*(1.-u)*Gi1j + t*u*Gi1j1 + (1.-t)*u*Gij1;  //phase at nearest higher freq bin
    }
    
    else {
        //cout<<"There is no antenna type : "<<ant_m<<" !!"<<endl;
        cout<<"Will return phase = 0 !!"<<endl;
        Gout1 = 0.;
        Gout2 = 0.;
    }
    
    Gout = ((Gout2 - Gout1)/freq_width) * ( freq - (freq_init + fx1*freq_width) ) + Gout1; // get linear interpolation between two nearest freq bin.
    
    
    return Gout;




}




double Detector::GetGain_1D( double freq, double theta, double phi, int ant_m ) {
    // find nearest theta, phi bin
    //
    //int i = (int)(theta/5.);
    //int j = (int)(phi/5.);

    // check if angles range actually theta 0-180, phi 0-360
    int i = (int)( (theta+2.5)/5. );
    int j = (int)( (phi+2.5)/5. );

    if ( j == 72 ) j = 0;

    int angle_bin = 37*j+i;

    // now just do linear interpolation at that angle
    //

    double slope_1, slope_2; // slope of init, final part

    double Gout;

    int bin = (int)( (freq - freq_init) / freq_width )+1;

    // Vpol
    if ( ant_m == 0 ) {

        slope_1 = (Vgain[1][angle_bin] - Vgain[0][angle_bin]) / (Freq[1] - Freq[0]);
        slope_2 = (Vgain[freq_step-1][angle_bin] - Vgain[freq_step-2][angle_bin]) / (Freq[freq_step-1] - Freq[freq_step-2]);


        // if freq is lower than freq_init
        if ( freq < freq_init ) {

            Gout = slope_1 * (freq - Freq[0]) + Vgain[0][angle_bin];
        }
        // if freq is higher than last freq
        else if ( freq > Freq[freq_step-1] ) {

            Gout = slope_2 * (freq - Freq[freq_step-1]) + Vgain[freq_step-1][angle_bin];
        }

        else {

            Gout = Vgain[bin-1][angle_bin] + (freq-Freq[bin-1])*(Vgain[bin][angle_bin]-Vgain[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
        } // not outside the Freq[] range
    
    } // Vpol case

    // Hpol
    else if ( ant_m == 1 ) {

        slope_1 = (Hgain[1][angle_bin] - Hgain[0][angle_bin]) / (Freq[1] - Freq[0]);
        slope_2 = (Hgain[freq_step-1][angle_bin] - Hgain[freq_step-2][angle_bin]) / (Freq[freq_step-1] - Freq[freq_step-2]);


        // if freq is lower than freq_init
        if ( freq < freq_init ) {

            Gout = slope_1 * (freq - Freq[0]) + Hgain[0][angle_bin];
        }
        // if freq is higher than last freq
        else if ( freq > Freq[freq_step-1] ) {

            Gout = slope_2 * (freq - Freq[freq_step-1]) + Hgain[freq_step-1][angle_bin];
        }

        else {

            Gout = Hgain[bin-1][angle_bin] + (freq-Freq[bin-1])*(Hgain[bin][angle_bin]-Hgain[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
        } // not outside the Freq[] range
    
    } // Hpol case


    return Gout;

}





// set outside value as 0
double Detector::GetGain_1D_OutZero( double freq, double theta, double phi, int ant_m ) {


    // find nearest theta, phi bin
    //
    //int i = (int)(theta/5.);
    //int j = (int)(phi/5.);

    // check if angles range actually theta 0-180, phi 0-360
    int i = (int)( (theta+2.5)/5. );
    int j = (int)( (phi+2.5)/5. );

    if ( j == 72 ) j = 0;

    int angle_bin = 37*j+i;

    // now just do linear interpolation at that angle
    //

    double slope_1; // slope of init part

    double Gout;

    int bin = (int)( (freq - freq_init) / freq_width )+1;

    // Vpol
    if ( ant_m == 0 ) {

        slope_1 = (Vgain[1][angle_bin] - Vgain[0][angle_bin]) / (Freq[1] - Freq[0]);


        // if freq is lower than freq_init
        if ( freq < freq_init ) {

            Gout = slope_1 * (freq - Freq[0]) + Vgain[0][angle_bin];
        }
        // if freq is higher than last freq
        else if ( freq > Freq[freq_step-1] ) {

            //Gout = slope_2 * (freq - Freq[freq_step-1]) + Vgain[freq_step-1][angle_bin];
            Gout = 0.;
        }

        else {

            Gout = Vgain[bin-1][angle_bin] + (freq-Freq[bin-1])*(Vgain[bin][angle_bin]-Vgain[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
        } // not outside the Freq[] range
    
    } // Vpol case

    // Hpol
    else if ( ant_m == 1 ) {

        slope_1 = (Hgain[1][angle_bin] - Hgain[0][angle_bin]) / (Freq[1] - Freq[0]);


        // if freq is lower than freq_init
        if ( freq < freq_init ) {

            Gout = slope_1 * (freq - Freq[0]) + Hgain[0][angle_bin];
        }
        // if freq is higher than last freq
        else if ( freq > Freq[freq_step-1] ) {

            //Gout = slope_2 * (freq - Freq[freq_step-1]) + Hgain[freq_step-1][angle_bin];
            Gout = 0.;
        }

        else {

            Gout = Hgain[bin-1][angle_bin] + (freq-Freq[bin-1])*(Hgain[bin][angle_bin]-Hgain[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
        } // not outside the Freq[] range
    
    } // Hpol case


    return Gout;

}


double Detector::GetGain_1D_OutZero( double freq, double theta, double phi, int ant_m, int ant_number) {


    // find nearest theta, phi bin
    //
    //int i = (int)(theta/5.);
    //int j = (int)(phi/5.);

    // check if angles range actually theta 0-180, phi 0-360
    int i = (int)( (theta+2.5)/5. );
    int j = (int)( (phi+2.5)/5. );

    if ( j == 72 ) j = 0;

    int angle_bin = 37*j+i;

    // now just do linear interpolation at that angle
    //

    double slope_1; // slope of init part

    double Gout;

    int bin = (int)( (freq - freq_init) / freq_width )+1;

    // Vpol
    if ( ant_m == 0 ) {

      if(ant_number==0){//bottom Vpol
        slope_1 = (Vgain[1][angle_bin] - Vgain[0][angle_bin]) / (Freq[1] - Freq[0]);


        // if freq is lower than freq_init
        if ( freq < freq_init ) {

            Gout = slope_1 * (freq - Freq[0]) + Vgain[0][angle_bin];
        }
        // if freq is higher than last freq
        else if ( freq > Freq[freq_step-1] ) {

            //Gout = slope_2 * (freq - Freq[freq_step-1]) + Vgain[freq_step-1][angle_bin];
            Gout = 0.;
        }

        else {

            Gout = Vgain[bin-1][angle_bin] + (freq-Freq[bin-1])*(Vgain[bin][angle_bin]-Vgain[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
        } // not outside the Freq[] range
      }//bottom Vpol
      else if(ant_number==2){//Top Vpol
//	cerr << "Does it happen: yes it happens! " << ant_number << endl;
        slope_1 = (VgainTop[1][angle_bin] - VgainTop[0][angle_bin]) / (Freq[1] - Freq[0]);


        // if freq is lower than freq_init
        if ( freq < freq_init ) {

            Gout = slope_1 * (freq - Freq[0]) + VgainTop[0][angle_bin];
        }
        // if freq is higher than last freq
        else if ( freq > Freq[freq_step-1] ) {

            //Gout = slope_2 * (freq - Freq[freq_step-1]) + Vgain[freq_step-1][angle_bin];
            Gout = 0.;
        }

        else {

            Gout = VgainTop[bin-1][angle_bin] + (freq-Freq[bin-1])*(VgainTop[bin][angle_bin]-VgainTop[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
	}
      }//Top Vpol
    
    } // Vpol case

    // Hpol
    else if ( ant_m == 1 ) {

        slope_1 = (Hgain[1][angle_bin] - Hgain[0][angle_bin]) / (Freq[1] - Freq[0]);


        // if freq is lower than freq_init
        if ( freq < freq_init ) {

            Gout = slope_1 * (freq - Freq[0]) + Hgain[0][angle_bin];
        }
        // if freq is higher than last freq
        else if ( freq > Freq[freq_step-1] ) {

            //Gout = slope_2 * (freq - Freq[freq_step-1]) + Hgain[freq_step-1][angle_bin];
            Gout = 0.;
        }

        else {

            Gout = Hgain[bin-1][angle_bin] + (freq-Freq[bin-1])*(Hgain[bin][angle_bin]-Hgain[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
        } // not outside the Freq[] range
    
    } // Hpol case


    if ( Gout < 0. ) // gain can not go below 0
        Gout = 0.;

    return Gout;

}






double Detector::GetAntPhase_1D( double freq, double theta, double phi, int ant_m ) {


    // find nearest theta, phi bin
    //
    //int i = (int)(theta/5.);
    //int j = (int)(phi/5.);

    // check if angles range actually theta 0-180, phi 0-360
    int i = (int)( (theta+2.5)/5. );
    int j = (int)( (phi+2.5)/5. );

    if ( j == 72 ) j = 0;

    int angle_bin = 37*j+i;

    // now just do linear interpolation at that angle
    //

    double slope_1, slope_2; // slope of init, final part
    double slope_t1, slope_t2; // slope of pre, after the freq bin

    double phase;

    int bin = (int)( (freq - freq_init) / freq_width )+1;

    // Vpol
    if ( ant_m == 0 ) {

        slope_1 = (Vphase[1][angle_bin] - Vphase[0][angle_bin]) / (Freq[1] - Freq[0]);
        slope_2 = (Vphase[freq_step-1][angle_bin] - Vphase[freq_step-2][angle_bin]) / (Freq[freq_step-1] - Freq[freq_step-2]);


        // if freq is lower than freq_init
        if ( freq < freq_init ) {

            phase = slope_1 * (freq - Freq[0]) + Vphase[0][angle_bin];

            if ( phase > 180. ) {
                while ( phase > 180. ) {
                    phase = phase - 360.;
                }
            }
            else if ( phase < -180. ) {
                while ( phase < -180. ) {
                    phase = phase + 360.;
                }
            }
        }
        // if freq is higher than last freq
        else if ( freq > Freq[freq_step-1] ) {

            phase = slope_2 * (freq - Freq[freq_step-1]) + Vphase[freq_step-1][angle_bin];

            if ( phase > 180. ) {
                while ( phase > 180. ) {
                    phase = phase - 360.;
                }
            }
            else if ( phase < -180. ) {
                while ( phase < -180. ) {
                    phase = phase + 360.;
                }
            }
        }

        else {

            // not at the first two bins
            if ( bin<freq_step-1 && bin>1 ) {

                slope_t1 = (Vphase[bin-1][angle_bin] - Vphase[bin-2][angle_bin]) / (Freq[bin-1] - Freq[bin-2]);
                slope_t2 = (Vphase[bin+1][angle_bin] - Vphase[bin][angle_bin]) / (Freq[bin+1] - Freq[bin]);

                // down going case
                if ( slope_t1 * slope_t2 > 0. && Vphase[bin][angle_bin] - Vphase[bin-1][angle_bin] > 180. ) {

                    phase = Vphase[bin-1][angle_bin] + (freq-Freq[bin-1])*(Vphase[bin][angle_bin]-360.-Vphase[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
                }

                // up going case
                else if ( slope_t1 * slope_t2 > 0. && Vphase[bin][angle_bin] - Vphase[bin-1][angle_bin] < -180. ) {
                    phase = Vphase[bin-1][angle_bin] + (freq-Freq[bin-1])*(Vphase[bin][angle_bin]+360.-Vphase[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
                }

                // neither case
                else {
                    phase = Vphase[bin-1][angle_bin] + (freq-Freq[bin-1])*(Vphase[bin][angle_bin]-Vphase[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
                }

                // if outside the range, put inside
                if ( phase > 180. ) {
                    while ( phase > 180. ) {
                        phase = phase - 360.;
                    }
                }
                else if ( phase < -180. ) {
                    while ( phase < -180. ) {
                        phase = phase + 360.;
                    }
                }

            }// not first two bins

            else {
                phase = Vphase[bin-1][angle_bin] + (freq-Freq[bin-1])*(Vphase[bin][angle_bin]-Vphase[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
            }

            // if outside the range, put inside
            if ( phase > 180. ) {
                while ( phase > 180. ) {
                    phase = phase - 360.;
                }
            }
            else if ( phase < -180. ) {
                while ( phase < -180. ) {
                    phase = phase + 360.;
                }
            }

        } // not outside the Freq[] range
    
    } // Vpol case

    // Hpol
    else if ( ant_m == 1 ) {

        slope_1 = (Hphase[1][angle_bin] - Hphase[0][angle_bin]) / (Freq[1] - Freq[0]);
        slope_2 = (Hphase[freq_step-1][angle_bin] - Hphase[freq_step-2][angle_bin]) / (Freq[freq_step-1] - Freq[freq_step-2]);


        // if freq is lower than freq_init
        if ( freq < freq_init ) {

            phase = slope_1 * (freq - Freq[0]) + Hphase[0][angle_bin];

            if ( phase > 180. ) {
                while ( phase > 180. ) {
                    phase = phase - 360.;
                }
            }
            else if ( phase < -180. ) {
                while ( phase < -180. ) {
                    phase = phase + 360.;
                }
            }
        }
        // if freq is higher than last freq
        else if ( freq > Freq[freq_step-1] ) {

            phase = slope_2 * (freq - Freq[freq_step-1]) + Hphase[freq_step-1][angle_bin];

            if ( phase > 180. ) {
                while ( phase > 180. ) {
                    phase = phase - 360.;
                }
            }
            else if ( phase < -180. ) {
                while ( phase < -180. ) {
                    phase = phase + 360.;
                }
            }
        }

        else {

            // not at the first two bins
            if ( bin<freq_step-1 && bin>1 ) {

                slope_t1 = (Hphase[bin-1][angle_bin] - Hphase[bin-2][angle_bin]) / (Freq[bin-1] - Freq[bin-2]);
                slope_t2 = (Hphase[bin+1][angle_bin] - Hphase[bin][angle_bin]) / (Freq[bin+1] - Freq[bin]);

                // down going case
                if ( slope_t1 * slope_t2 > 0. && Hphase[bin][angle_bin] - Hphase[bin-1][angle_bin] > 180. ) {

                    phase = Hphase[bin-1][angle_bin] + (freq-Freq[bin-1])*(Hphase[bin][angle_bin]-360.-Hphase[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
                }

                // up going case
                else if ( slope_t1 * slope_t2 > 0. && Hphase[bin][angle_bin] - Hphase[bin-1][angle_bin] < -180. ) {
                    phase = Hphase[bin-1][angle_bin] + (freq-Freq[bin-1])*(Hphase[bin][angle_bin]+360.-Hphase[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
                }

                // neither case
                else {
                    phase = Hphase[bin-1][angle_bin] + (freq-Freq[bin-1])*(Hphase[bin][angle_bin]-Hphase[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
                }

                // if outside the range, put inside
                if ( phase > 180. ) {
                    while ( phase > 180. ) {
                        phase = phase - 360.;
                    }
                }
                else if ( phase < -180. ) {
                    while ( phase < -180. ) {
                        phase = phase + 360.;
                    }
                }

            }// not first two bins

            else {
                phase = Hphase[bin-1][angle_bin] + (freq-Freq[bin-1])*(Hphase[bin][angle_bin]-Hphase[bin-1][angle_bin])/(Freq[bin]-Freq[bin-1]);
            }

            // if outside the range, put inside
            if ( phase > 180. ) {
                while ( phase > 180. ) {
                    phase = phase - 360.;
                }
            }
            else if ( phase < -180. ) {
                while ( phase < -180. ) {
                    phase = phase + 360.;
                }
            }

        } // not outside the Freq[] range
    
    } // Hpol case


    return phase;

}



// set outside value as 0
double Detector::GetFilterGain_1D_OutZero( double freq ) {


    double slope_1; // slope of init part

    double Gout;

    int bin = (int)( (freq - freq_init) / freq_width )+1;


    slope_1 = (FilterGain[1] - FilterGain[0]) / (Freq[1] - Freq[0]);


    // if freq is lower than freq_init
    if ( freq < freq_init ) {

        Gout = slope_1 * (freq - Freq[0]) + FilterGain[0];
    }
    // if freq is higher than last freq
    else if ( freq > Freq[freq_step-1] ) {

        //Gout = slope_2 * (freq - Freq[freq_step-1]) + FilterGain[freq_step-1];
        Gout = 0.;
    }

    else {

        Gout = FilterGain[bin-1] + (freq-Freq[bin-1])*(FilterGain[bin]-FilterGain[bin-1])/(Freq[bin]-Freq[bin-1]);
    } // not outside the Freq[] range
    

    
    return Gout;

}



// set outside value as 0
double Detector::GetPreampGain_1D_OutZero( double freq ) {


    double slope_1; // slope of init part

    double Gout;

    int bin = (int)( (freq - freq_init) / freq_width )+1;


    slope_1 = (PreampGain[1] - PreampGain[0]) / (Freq[1] - Freq[0]);


    // if freq is lower than freq_init
    if ( freq < freq_init ) {

        Gout = slope_1 * (freq - Freq[0]) + PreampGain[0];
    }
    // if freq is higher than last freq
    else if ( freq > Freq[freq_step-1] ) {

        //Gout = slope_2 * (freq - Freq[freq_step-1]) + PreampGain[freq_step-1];
        Gout = 0.;
    }

    else {

        Gout = PreampGain[bin-1] + (freq-Freq[bin-1])*(PreampGain[bin]-PreampGain[bin-1])/(Freq[bin]-Freq[bin-1]);
    } // not outside the Freq[] range
    


    return Gout;

}




// set outside value as 0
double Detector::GetFOAMGain_1D_OutZero( double freq ) {


    double slope_1; // slope of init part

    double Gout;

    int bin = (int)( (freq - freq_init) / freq_width )+1;


    slope_1 = (FOAMGain[1] - FOAMGain[0]) / (Freq[1] - Freq[0]);


    // if freq is lower than freq_init
    if ( freq < freq_init ) {

        Gout = slope_1 * (freq - Freq[0]) + FOAMGain[0];
    }
    // if freq is higher than last freq
    else if ( freq > Freq[freq_step-1] ) {

        //Gout = slope_2 * (freq - Freq[freq_step-1]) + FOAMGain[freq_step-1];
        Gout = 0.;
    }

    else {

        Gout = FOAMGain[bin-1] + (freq-Freq[bin-1])*(FOAMGain[bin]-FOAMGain[bin-1])/(Freq[bin]-Freq[bin-1]);
    } // not outside the Freq[] range
    


    return Gout;

}




// set outside value as 0
double Detector::GetElectGain_1D_OutZero( double freq ) {


    double slope_1; // slope of init part

    double Gout;

    int bin = (int)( (freq - freq_init) / freq_width )+1;


    slope_1 = (ElectGain[1] - ElectGain[0]) / (Freq[1] - Freq[0]);


    // if freq is lower than freq_init
    if ( freq < freq_init ) {

        //Gout = slope_1 * (freq - Freq[0]) + ElectGain[0];
        Gout = 0.;
    }
    // if freq is higher than last freq
    else if ( freq > Freq[freq_step-1] ) {

        //Gout = slope_2 * (freq - Freq[freq_step-1]) + FOAMGain[freq_step-1];
        Gout = 0.;
    }

    else {

        Gout = ElectGain[bin-1] + (freq-Freq[bin-1])*(ElectGain[bin]-ElectGain[bin-1])/(Freq[bin]-Freq[bin-1]);
    } // not outside the Freq[] range
    


    return Gout;

}



// set outside value as 0
double Detector::GetElectPhase_1D( double freq ) {


    double slope_1, slope_2; // slope of init, final part
    double slope_t1, slope_t2; // slope of pre, after the freq bin

    double phase;

    int bin = (int)( (freq - freq_init) / freq_width )+1;

    // ElectPhase are in rad (not deg)

    slope_1 = (ElectPhase[1] - ElectPhase[0]) / (Freq[1] - Freq[0]);
    slope_2 = (ElectPhase[freq_step-1] - ElectPhase[freq_step-2]) / (Freq[freq_step-1] - Freq[freq_step-2]);


    // if freq is lower than freq_init
    if ( freq < freq_init ) {

        phase = slope_1 * (freq - Freq[0]) + ElectPhase[0];

        if ( phase > PI ) {
            while ( phase > PI ) {
                phase = phase - 2*PI;
            }
        }
        else if ( phase < -PI ) {
            while ( phase < -PI ) {
                phase = phase + 2*PI;
            }
        }
    }
    // if freq is higher than last freq
    else if ( freq > Freq[freq_step-1] ) {

        phase = slope_2 * (freq - Freq[freq_step-1]) + ElectPhase[freq_step-1];

        if ( phase > PI ) {
            while ( phase > PI ) {
                phase = phase - 2*PI;
            }
        }
        else if ( phase < -PI ) {
            while ( phase < -PI ) {
                phase = phase + 2*PI;
            }
        }
    }

    else {

        // not at the first two bins
        if ( bin<freq_step-1 && bin>1 ) {

            slope_t1 = (ElectPhase[bin-1] - ElectPhase[bin-2]) / (Freq[bin-1] - Freq[bin-2]);
            slope_t2 = (ElectPhase[bin+1] - ElectPhase[bin]) / (Freq[bin+1] - Freq[bin]);

            // down going case
            if ( slope_t1 * slope_t2 > 0. && ElectPhase[bin] - ElectPhase[bin-1] > PI ) {

                phase = ElectPhase[bin-1] + (freq-Freq[bin-1])*(ElectPhase[bin]-2*PI-ElectPhase[bin-1])/(Freq[bin]-Freq[bin-1]);
            }

            // up going case
            else if ( slope_t1 * slope_t2 > 0. && ElectPhase[bin] - ElectPhase[bin-1] < -PI ) {
                phase = ElectPhase[bin-1] + (freq-Freq[bin-1])*(ElectPhase[bin]+2*PI-ElectPhase[bin-1])/(Freq[bin]-Freq[bin-1]);
            }

            // neither case
            else {
                phase = ElectPhase[bin-1] + (freq-Freq[bin-1])*(ElectPhase[bin]-ElectPhase[bin-1])/(Freq[bin]-Freq[bin-1]);
            }

            // if outside the range, put inside
            if ( phase > PI ) {
                while ( phase > PI ) {
                    phase = phase - 2*PI;
                }
            }
            else if ( phase < -PI ) {
                while ( phase < -PI ) {
                    phase = phase + 2*PI;
                }
            }

        }// not first two bins

        else {
            phase = ElectPhase[bin-1] + (freq-Freq[bin-1])*(ElectPhase[bin]-ElectPhase[bin-1])/(Freq[bin]-Freq[bin-1]);
        }

        // if outside the range, put inside
        if ( phase > PI ) {
            while ( phase > PI ) {
                phase = phase - 2*PI;
            }
        }
        else if ( phase < -PI ) {
            while ( phase < -PI ) {
                phase = phase + 2*PI;
            }
        }

    } // not outside the Freq[] range
    


    return phase;


}





double Antenna::GetG(Detector *D, double freq, double theta, double phi) {
    
    return D->GetGain(freq, theta, phi, type, orient);
}




double Surface_antenna::GetG(Detector *D, double freq, double theta, double phi) {
    
    return D->GetGain(freq, theta, phi, type, orient);
}


inline void Detector::FlattoEarth_ARA(IceModel *icesurface) {
    
    double Dist = 0.;   //for sqrt(x^2 + y^2)
    double R1 = icesurface->Surface(0.,0.); // from core of earth to surface at theta, phi = 0.
    //--------------------------------------------------
    //     double R1 = icesurface->Geoid(0.); // from core of earth to surface at theta, phi = 0.
    //-------------------------------------------------- 
    double theta_tmp;
    double phi_tmp;
    
    // stations
    // stations, strings, and borehole antennas use geoid surface !!
    for (int i=0; i<params.number_of_stations; i++) {
        
        Dist = sqrt( pow(stations[i].GetX(),2) + pow(stations[i].GetY(),2) );
        theta_tmp = Dist/R1;
        phi_tmp = atan2(stations[i].GetY(),stations[i].GetX());
        
        if (phi_tmp<0.) phi_tmp += 2.*PI;
        
        // set theta, phi for stations.
        stations[i].SetThetaPhi(theta_tmp, phi_tmp);
        //set R for stations.
        stations[i].SetR( icesurface->Surface( stations[i].Lon(), stations[i].Lat()) );
        
        
        // strings
        for (int j=0; j<params.number_of_strings_station; j++) {
            Dist = sqrt( pow(stations[i].strings[j].GetX(),2) + pow(stations[i].strings[j].GetY(),2) );
            theta_tmp = Dist/R1;
            phi_tmp = atan2(stations[i].strings[j].GetY(),stations[i].strings[j].GetX());
            
            if (phi_tmp<0.) phi_tmp += 2.*PI;
            
            stations[i].strings[j].SetThetaPhi(theta_tmp, phi_tmp);
            // string Vector points the position where string meets the ice surface!
            stations[i].strings[j].SetR( icesurface->Surface( stations[i].strings[j].Lon(), stations[i].strings[j].Lat()) );
            
            
            
            // borehole antennas
            for (int k=0; k<params.number_of_antennas_string; k++) {
                stations[i].strings[j].antennas[k].SetRThetaPhi( stations[i].strings[j].R() + stations[i].strings[j].antennas[k].GetZ() , stations[i].strings[j].Theta(), stations[i].strings[j].Phi() );
            }
            
            
        }
        
        // surface antennas
        // surface antennas are on actual ice surface (not geoid surface)
        for (int l=0; l<params.number_of_surfaces_station; l++) {
            Dist = sqrt( pow(stations[i].surfaces[l].GetX(),2) + pow(stations[i].surfaces[l].GetY(),2) );
            theta_tmp = Dist/R1;
            phi_tmp = atan2(stations[i].surfaces[l].GetY(),stations[i].surfaces[l].GetX());
            
            if (phi_tmp<0.) phi_tmp += 2.*PI;
            
            stations[i].surfaces[l].SetThetaPhi(theta_tmp, phi_tmp);
            stations[i].surfaces[l].SetR( icesurface->Surface( stations[i].surfaces[l].Lon(), stations[i].surfaces[l].Lat()) );
        }
        
        
    }
    
    
}


inline void Detector::FlattoEarth_ARA_sharesurface(IceModel *icesurface) {    // each station share the lowest surface
    
    double Dist = 0.;   //for sqrt(x^2 + y^2)
    double R1 = icesurface->Surface(0.,0.); // from core of earth to surface at theta, phi = 0.
    //--------------------------------------------------
    //     double R1 = icesurface->Geoid(0.); // from core of eart to surface at theta, phi = 0.
    //-------------------------------------------------- 
    double theta_tmp;
    double phi_tmp;
    
    double lowest_surface;  // lowest surface of the string among the station
    
    // stations
    // stations, strings, and borehole antennas use geoid surface !!
    for (int i=0; i<int(stations.size()); i++) {
        
        Dist = sqrt( pow(stations[i].GetX(),2) + pow(stations[i].GetY(),2) );
        theta_tmp = Dist/R1;
        phi_tmp = atan2(stations[i].GetY(),stations[i].GetX());
        
        if (phi_tmp<0.) phi_tmp += 2.*PI;
        
        // set theta, phi for stations.
        stations[i].SetThetaPhi(theta_tmp, phi_tmp);
        //set R for stations.
        stations[i].SetR( icesurface->Surface( stations[i].Lon(), stations[i].Lat()) );
        
        lowest_surface = 1.E7;  // much bigger than the surface (approx radius of earth 6.E6)
        
        // strings
        for (int j=0; j<int(stations[i].strings.size()); j++) {
            Dist = sqrt( pow(stations[i].strings[j].GetX(),2) + pow(stations[i].strings[j].GetY(),2) );
            theta_tmp = Dist/R1;
            phi_tmp = atan2(stations[i].strings[j].GetY(),stations[i].strings[j].GetX());
            
            if (phi_tmp<0.) phi_tmp += 2.*PI;
            
            stations[i].strings[j].SetThetaPhi(theta_tmp, phi_tmp);
            // string Vector points the position where string meets the ice surface!
            stations[i].strings[j].SetR( icesurface->Surface( stations[i].strings[j].Lon(), stations[i].strings[j].Lat()) );
            
            
            // find the lowest surface among strings in a station
            if ( lowest_surface > stations[i].strings[j].R() ) {
                lowest_surface = stations[i].strings[j].R();
            }
            
            
            
        }
        
        
        // string loop again for borehole antennas
        for (int j=0; j<int(stations[i].strings.size()); j++) {
            // borehole antennas
            for (int k=0; k<int(stations[i].strings[j].antennas.size()); k++) {
                stations[i].strings[j].antennas[k].SetRThetaPhi( lowest_surface + stations[i].strings[j].antennas[k].GetZ() , stations[i].strings[j].Theta(), stations[i].strings[j].Phi() );
            }
        } // end string loop for borehole antennas
        
        
        
        // surface antennas
        // surface antennas are on actual ice surface (not geoid surface)
        for (int l=0; l<int(stations[i].surfaces.size()); l++) {
            Dist = sqrt( pow(stations[i].surfaces[l].GetX(),2) + pow(stations[i].surfaces[l].GetY(),2) );
            theta_tmp = Dist/R1;
            phi_tmp = atan2(stations[i].surfaces[l].GetY(),stations[i].surfaces[l].GetX());
            
            if (phi_tmp<0.) phi_tmp += 2.*PI;
            
            stations[i].surfaces[l].SetThetaPhi(theta_tmp, phi_tmp);
            stations[i].surfaces[l].SetR( icesurface->Surface( stations[i].surfaces[l].Lon(), stations[i].surfaces[l].Lat()) );
        }
        
        
    } // end loop over stations
    
    for (int i = 0; i < int(stations.size()); i++){
        for (int j = 0; j < int(stations[i].strings.size()); j++){
            for (int k = 0; k < int(stations[i].strings[j].antennas.size()); k++){
        cout << "Sharesurface: " <<
//        GetChannelfromStringAntenna( i+1 , j , k) << " : " <<
        i << " : " <<
        j << " : " <<
        k << " : " <<
        stations[i].strings[j].antennas[k].GetX() << " : " <<
        stations[i].strings[j].antennas[k].GetY() << " : " <<
        stations[i].strings[j].antennas[k].GetZ() << " : " <<
        endl;
            }
        }
    }
    
    
    
}





inline void Detector::AddAdditional_Depth(Settings *settings1) {    // each station share the lowest surface
    

    if (settings1->ADDITIONAL_DEPTH_ON == 1) {

        // stations
        for (int i=0; i<int(stations.size()); i++) {

            // strings
            for (int j=0; j<int(stations[i].strings.size()); j++) {

                // antennas
                for (int k = 0; k < int(stations[i].strings[j].antennas.size()); k++){

                    //stations[i].strings[j].antennas[k].SetR( stations[i].strings[j].antennas[k].GetR() - settings1->ADDITIONAL_DEPTH );
                    stations[i].strings[j].antennas[k].SetZ( stations[i].strings[j].antennas[k].GetZ() - settings1->ADDITIONAL_DEPTH );
                }
            }
        }

    }


}














inline void Detector::ReadFilter(string filename, Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain
    
    ifstream Filter( filename.c_str() );
    
    string line;
    
    int N=-1;
    
    vector <double> xfreq_tmp;
    vector <double> ygain_tmp;
    
    if ( Filter.is_open() ) {
        while (Filter.good() ) {
            
            getline (Filter, line);
            //xfreq[N] = atof( line.substr(0, line.find_first_of(",")).c_str() );
            xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) );
            //xfreq.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) * 1.E6 );  // from MHz to Hz
            
            //xfreq[N] = xfreq[N] * 1.E6; // from MHz to Hz
            
            //ygain[N] = atof( line.substr(line.find_first_of(",") + 1).c_str() );
            ygain_tmp.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) );
            //ygain.push_back( pow(pow(10,atof( line.substr(line.find_first_of(",") + 1).c_str() ) /10.),0.5) );  // from dB to unitless gain for voltage
            
            //ygain[N] = pow(pow(10,yFilter[i]/10.0),0.5);    // from dB to field strength unitless gain
            
            N++;
            
        }
        Filter.close();
    }
    
    else cout<<"Filter file can not opened!!"<<endl;
    
    double xfreq[N], ygain[N];  // need array for Tools::SimpleLinearInterpolation
    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double ygain_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    for (int i=0;i<N;i++) { // copy values
        xfreq[i] = xfreq_tmp[i];
        ygain[i] = ygain_tmp[i];
    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    
    
    // Tools::SimpleLinearInterpolation will return Filter array (in dB)
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, freq_step, Freq, FilterGain );
    
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->DATA_BIN_SIZE/2, xfreq_databin, ygain_databin );
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        FilterGain_databin.push_back( ygain_databin[i] );
    }
    
    
    // for NFOUR/2 t domain array
    double xfreq_NFOUR[settings1->NFOUR/4+1];   // array for FFT freq bin
    double ygain_NFOUR[settings1->NFOUR/4+1];   // array for gain in FFT bin
    
    df_fft = 1./ ( (double)(settings1->NFOUR/2) * settings1->TIMESTEP );

    for (int i=0;i<settings1->NFOUR/4+1;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_NFOUR[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->NFOUR/4+1, xfreq_NFOUR, ygain_NFOUR );
    
    for (int i=0;i<settings1->NFOUR/4+1;i++) {
        FilterGain_NFOUR.push_back( ygain_NFOUR[i] );
    }
    
}


void Detector::ReadFilter_New(Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain

    // We can use FilterGain array as a original array

    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double ygain_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );

    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    Tools::SimpleLinearInterpolation( freq_step, Freq, FilterGain, settings1->DATA_BIN_SIZE/2, xfreq_databin, ygain_databin );
        
    FilterGain_databin.clear();
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        FilterGain_databin.push_back( ygain_databin[i] );
    }

}


inline void Detector::ReadPreamp(string filename, Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain
    
    ifstream Preampgain( filename.c_str() );
    
    string line;
    
    int N=-1;
    
    vector <double> xfreq_tmp;
    vector <double> ygain_tmp;
    
    if ( Preampgain.is_open() ) {
        while (Preampgain.good() ) {
            
            getline (Preampgain, line);
            //xfreq[N] = atof( line.substr(0, line.find_first_of(",")).c_str() );
            xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) );
            //xfreq.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) * 1.E6 );  // from MHz to Hz
            
            //xfreq[N] = xfreq[N] * 1.E6; // from MHz to Hz
            
            //ygain[N] = atof( line.substr(line.find_first_of(",") + 1).c_str() );
            ygain_tmp.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) );
            //ygain.push_back( pow(pow(10,atof( line.substr(line.find_first_of(",") + 1).c_str() ) /10.),0.5) );  // from dB to unitless gain for voltage
            
            //ygain[N] = pow(pow(10,yFilter[i]/10.0),0.5);    // from dB to field strength unitless gain
            
            N++;
            
        }
        Preampgain.close();
    }
    
    else cout<<"Preamgain file can not opened!!"<<endl;
    
    double xfreq[N], ygain[N];  // need array for Tools::SimpleLinearInterpolation
    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double ygain_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    for (int i=0;i<N;i++) { // copy values
        xfreq[i] = xfreq_tmp[i];
        ygain[i] = ygain_tmp[i];
    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    
    
    // Tools::SimpleLinearInterpolation will return Preampgain array (in dB)
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, freq_step, Freq, PreampGain );
    
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->DATA_BIN_SIZE/2, xfreq_databin, ygain_databin );
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        PreampGain_databin.push_back( ygain_databin[i] );
    }
    
    
    // for NFOUR/2 t domain array
    double xfreq_NFOUR[settings1->NFOUR/4+1];   // array for FFT freq bin
    double ygain_NFOUR[settings1->NFOUR/4+1];   // array for gain in FFT bin
    
    df_fft = 1./ ( (double)(settings1->NFOUR/2) * settings1->TIMESTEP );

    for (int i=0;i<settings1->NFOUR/4+1;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_NFOUR[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->NFOUR/4+1, xfreq_NFOUR, ygain_NFOUR );
    
    for (int i=0;i<settings1->NFOUR/4+1;i++) {
        PreampGain_NFOUR.push_back( ygain_NFOUR[i] );
    }
    
    
}



void Detector::ReadPreamp_New(Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain

    // We can use FilterGain array as a original array

    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double ygain_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );

    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    Tools::SimpleLinearInterpolation( freq_step, Freq, PreampGain, settings1->DATA_BIN_SIZE/2, xfreq_databin, ygain_databin );
        
    PreampGain_databin.clear();
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        PreampGain_databin.push_back( ygain_databin[i] );
    }

}



inline void Detector::ReadFOAM(string filename, Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain
    
    ifstream FOAMgain( filename.c_str() );
    
    string line;
    
    int N=-1;
    
    vector <double> xfreq_tmp;
    vector <double> ygain_tmp;
    
    if ( FOAMgain.is_open() ) {
        while (FOAMgain.good() ) {
            
            getline (FOAMgain, line);
            //xfreq[N] = atof( line.substr(0, line.find_first_of(",")).c_str() );
            xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) );
            //xfreq.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) * 1.E6 );  // from MHz to Hz
            
            //xfreq[N] = xfreq[N] * 1.E6; // from MHz to Hz
            
            //ygain[N] = atof( line.substr(line.find_first_of(",") + 1).c_str() );
            ygain_tmp.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) );
            //ygain.push_back( pow(pow(10,atof( line.substr(line.find_first_of(",") + 1).c_str() ) /10.),0.5) );  // from dB to unitless gain for voltage
            
            //ygain[N] = pow(pow(10,yFilter[i]/10.0),0.5);    // from dB to field strength unitless gain
            
            N++;
            
        }
        FOAMgain.close();
    }
    
    else cout<<"Preamgain file can not opened!!"<<endl;
    
    double xfreq[N], ygain[N];  // need array for Tools::SimpleLinearInterpolation
    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double ygain_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    for (int i=0;i<N;i++) { // copy values
        xfreq[i] = xfreq_tmp[i];
        ygain[i] = ygain_tmp[i];
    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    
    
    // Tools::SimpleLinearInterpolation will return FOAMgain array (in dB)
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, freq_step, Freq, FOAMGain );
    
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->DATA_BIN_SIZE/2, xfreq_databin, ygain_databin );
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        FOAMGain_databin.push_back( ygain_databin[i] );
    }
    

    // for NFOUR/2 t domain array
    double xfreq_NFOUR[settings1->NFOUR/4+1];   // array for FFT freq bin
    double ygain_NFOUR[settings1->NFOUR/4+1];   // array for gain in FFT bin
    
    df_fft = 1./ ( (double)(settings1->NFOUR/2) * settings1->TIMESTEP );

    for (int i=0;i<settings1->NFOUR/4+1;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_NFOUR[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->NFOUR/4+1, xfreq_NFOUR, ygain_NFOUR );
    
    for (int i=0;i<settings1->NFOUR/4+1;i++) {
        FOAMGain_NFOUR.push_back( ygain_NFOUR[i] );
    }
}

void Detector::ReadNoiseFigure(string filename, Settings *settings1)
{
cout<<"In ReadNoiseFigure"<<endl;    
    ifstream nfFile( filename.c_str() );
    
    string line;
    
    int N=0;
    
    vector< vector<double> > all_chNF;
//    all_chNF.resize(17);

    vector <double> xfreq_tmp;
    vector <double> yNF_tmp;
    
    if ( nfFile.is_open() ) {
        while (nfFile.good() ) {
            
            getline (nfFile, line);

	    istringstream iss(line);

	    while( iss )
	    {
	        string sub;
	        iss >> sub;
		yNF_tmp.push_back( atof(sub.c_str()) );
	    }
	    cout << "freq: " << yNF_tmp[0];
	    all_chNF.push_back( yNF_tmp );
	    yNF_tmp.clear();
	    cout << "   ch1: " << all_chNF[N][1] << endl;
            
            N++;
            
        }
        nfFile.close();
    }
    else cout<<"Noise Figure file can not be opened!!"<<endl;
     
    double xfreq[N];  // need array for Tools::SimpleLinearInterpolation
    double NoiseFig[N];

    int ch_no = 16;// all_chNF[1].size()-2;
    cerr << "The number of channels: " << ch_no << endl;

    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double NoiseFig_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
   cout<<"DATA_BIN_SIZE: "<<settings1->DATA_BIN_SIZE<<" TIMESTEP: "<<settings1->TIMESTEP<<" df_fft: "<<df_fft<<endl; 
    // now below are values that shared in all channels
    for (int i=0;i<N;i++) { // copy values
        xfreq[i] = all_chNF[i][0];

        /*
        for (int ch=0; ch<ch_no; ch++) {
            ygain[i] = ygain_tmp[i];
        }
        */
    }

//    FILE *fp=fopen("xfreq_databin_file_NOISE_CHANNEL_MODE_1.txt","a+");

cout<<"number of f bins: "<<settings1->DATA_BIN_SIZE/2<<endl;


    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
      
//fprintf(fp,"%f   ",xfreq_databin[i]);

    }

//fprintf(fp,"\n");
//fclose(fp);

    // set vector array size to number of chs
    NoiseFig_databin_ch.resize(ch_no);
    //    cout<<"ch_no: "<<ch_no<<endl;        

    // now loop over channels and do interpolation
    for (int ch=0; ch<ch_no; ch++) {

        // copy fit values
      //        cout << "Test reading " << ch;
        for (int i=0;i<N;i++) {
            NoiseFig[i] = (all_chNF[i][ch+1]);
//	    if(i>30 && i<60) cout << xfreq[i] << "   " << NoiseFig[i] << "\t";
        }
	//	cout << endl;

        // Tools::SimpleLinearInterpolation will return NoiseFig array (in dB)
        Tools::SimpleLinearInterpolation( N-1, xfreq, NoiseFig, freq_step, Freq, NoiseFig_ch[ch] );
	//	cout << "2nd Test reading 2: ";
//        for (int i=0;i<30;i++){
//		cout << Freq[i] << "   " <<  NoiseFig_ch[ch][i] << "\t";
//	} 
//	cout << endl;
	Tools::SimpleLinearInterpolation( N-1, xfreq, NoiseFig, settings1->DATA_BIN_SIZE/2, xfreq_databin, NoiseFig_databin );
	//	cout << "Third test reading next ";
//        for (int i=0;i<30;i++){
//		cout << xfreq_databin[i+3000] << "   " << NoiseFig_databin[i+3000] << "\t";
//	} 
//	cout << endl;

    
        for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
/*if(i%100==0)*/ //cout<<"NoiseFig_databin: "<<NoiseFig_databin[i]<<endl;
            NoiseFig_databin_ch[ch].push_back( NoiseFig_databin[i] );
        }
    }
}


void Detector::ReadFOAM_New(Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain

    // We can use FilterGain array as a original array

    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double ygain_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );

    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    Tools::SimpleLinearInterpolation( freq_step, Freq, FOAMGain, settings1->DATA_BIN_SIZE/2, xfreq_databin, ygain_databin );
        
    FOAMGain_databin.clear();
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        FOAMGain_databin.push_back( ygain_databin[i] );
    }

}





inline void Detector::ReadCalPulserWF(string filename, Settings *settings1 ) {    // will store calpulser waveform array
    
    ifstream CalPulWF( filename.c_str() );
    
    string line;
    
    int N=-1;
    
    //vector <double> CalPulserWF_ns;
    //vector <double> CalPulserWF_V;
    CalPulserWF_ns.clear();
    CalPulserWF_V.clear();
    
    int firstread = 1;

    if ( CalPulWF.is_open() ) {
        while (CalPulWF.good() ) {

            // skip first two lines
            if ( firstread == 1 ) {
                getline (CalPulWF, line);
                getline (CalPulWF, line);
                firstread++;
            }
            
            getline (CalPulWF, line);
            //CalPulserWF_ns.push_back( atof( line.substr(0, 4).c_str() ) );
            //CalPulserWF_V.push_back( atof( line.substr(5).c_str() ) * settings1->CALPUL_AMP );
            CalPulserWF_ns.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) );
            CalPulserWF_V.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) * settings1->CALPUL_AMP );
            
            N++;


        }
        CalPulWF.close();
    }
    
    else cout<<"CalPulserWF file can not opened!!"<<endl;
    
    // remove last element
    CalPulserWF_ns.pop_back();
    CalPulserWF_V.pop_back();

    /*
    for (int i=0; i<(int)CalPulserWF_ns.size(); i++) {
            cout<<"CalPulser["<<i<<"] "<<CalPulserWF_ns[i]<<"ns, "<<CalPulserWF_V[i]<<"V"<<endl;
    }
    */

    cout<<"done reading CalPulserWF file"<<endl;
    
    
}














inline void Detector::ReadElectChain(string filename, Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain
    
    ifstream Elect( filename.c_str() );
    
    string line;
    string line2;
    string line3;
    
    int N=-1;
    
    vector <double> xfreq_tmp;
    vector <double> ygain_tmp;
    vector <double> phase_tmp;

    int skipline = 3;
    int first_time = 1;

    
    if ( Elect.is_open() ) {
        while (Elect.good() ) {
            
            if ( first_time == 1 ) {
                for (int sl=0; sl<skipline; sl++) {
                    getline (Elect, line);
                }
                first_time = 0;
            }
                    

            getline (Elect, line);

            N++;

            xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) );
            cout<<"freq : "<<xfreq_tmp[N]<<"\t";

            line2 = line.substr( line.find_first_of(",")+1);

            ygain_tmp.push_back( atof( line2.substr(0, line2.find_first_of(",")).c_str() ) );
            cout<<"gain : "<<ygain_tmp[N]<<"\t";

            line3 = line2.substr( line2.find_first_of(",")+1);

            phase_tmp.push_back( atof( line3.substr(0).c_str() ) );
            cout<<"phase : "<<phase_tmp[N]<<" N : "<<N<<endl;

            /*
            xfreq_tmp.push_back( atof( line.substr(0, 10).c_str() ) );
            cout<<"freq : "<<xfreq_tmp[N]<<"\t";
            
            ygain_tmp.push_back( atof( line.substr(11, 19).c_str() ) );
            cout<<"gain : "<<ygain_tmp[N]<<"\t";

            phase_tmp.push_back( atof( line.substr(31).c_str() ) );
            cout<<"phase : "<<phase_tmp[N]<<" N : "<<N<<endl;
            */
            
        }
        Elect.close();
    }
    
    else cout<<"Elect file can not opened!!"<<endl;

    cout<<"N : "<<N<<endl;
    
    double xfreq[N], ygain[N], phase[N];  // need array for Tools::SimpleLinearInterpolation
    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double ygain_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double phase_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    for (int i=0;i<N;i++) { // copy values
        xfreq[i] = xfreq_tmp[i];
        ygain[i] = ygain_tmp[i];
        phase[i] = phase_tmp[i];
    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    
    
    // Tools::SimpleLinearInterpolation will return Filter array (in dB)
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, freq_step, Freq, ElectGain );
    
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->DATA_BIN_SIZE/2, xfreq_databin, ygain_databin );
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        ElectGain_databin.push_back( ygain_databin[i] );
    }

    Tools::SimpleLinearInterpolation( N, xfreq, phase, freq_step, Freq, ElectPhase );
    
    Tools::SimpleLinearInterpolation( N, xfreq, phase, settings1->DATA_BIN_SIZE/2, xfreq_databin, phase_databin );
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        ElectPhase_databin.push_back( phase_databin[i] );
    }

    
    
    // for NFOUR/2 t domain array
    double xfreq_NFOUR[settings1->NFOUR/4+1];   // array for FFT freq bin
    double ygain_NFOUR[settings1->NFOUR/4+1];   // array for gain in FFT bin
    double phase_NFOUR[settings1->NFOUR/4+1];   // array for gain in FFT bin
    
    df_fft = 1./ ( (double)(settings1->NFOUR/2) * settings1->TIMESTEP );

    for (int i=0;i<settings1->NFOUR/4+1;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_NFOUR[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->NFOUR/4+1, xfreq_NFOUR, ygain_NFOUR );
    
    for (int i=0;i<settings1->NFOUR/4+1;i++) {
        ElectGain_NFOUR.push_back( ygain_NFOUR[i] );
    }

    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->NFOUR/4+1, xfreq_NFOUR, phase_NFOUR );
    
    for (int i=0;i<settings1->NFOUR/4+1;i++) {
        ElectPhase_NFOUR.push_back( phase_NFOUR[i] );
    }

    
}



void Detector::ReadElectChain_New(Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain

    // We can use FilterGain array as a original array

    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double ygain_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double phase_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin

    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );

    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    Tools::SimpleLinearInterpolation( freq_step, Freq, ElectGain, settings1->DATA_BIN_SIZE/2, xfreq_databin, ygain_databin );
        
    ElectGain_databin.clear();
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        ElectGain_databin.push_back( ygain_databin[i] );
    }

    Tools::SimpleLinearInterpolation( freq_step, Freq, ElectPhase, settings1->DATA_BIN_SIZE/2, xfreq_databin, phase_databin );
        
    ElectPhase_databin.clear();
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        ElectPhase_databin.push_back( phase_databin[i] );
    }


}






inline void Detector::ReadGainOffset_TestBed(string filename, Settings *settings1) {    // will return gain offset (unit in voltage) for different chs

    if (settings1->USE_MANUAL_GAINOFFSET == 1) {

        if (settings1->USE_CH_GAINOFFSET == 1) {

            ifstream GainOffset( filename.c_str() );
            
            string line;
            
            int N=0;
            
            if ( GainOffset.is_open() ) {
                while (GainOffset.good() ) {
                    
                    getline (GainOffset, line);
                    
                    GainOffset_TB_ch.push_back( atof( line.c_str() ) );
                    cout<<"GainOffset ch"<<N<<" : "<<GainOffset_TB_ch[N]<<endl;
                    
                    N++;
                    
                }
                GainOffset.close();
            }
        }

        else {// not using gain offset file, but just use same gainoffset value to all channels
            for (int N=0; N<16; N++) {
                GainOffset_TB_ch.push_back( settings1->MANUAL_GAINOFFSET_VALUE );
            }
        }

    }
    // use manual gain offset from setup file
    else if (settings1->USE_MANUAL_GAINOFFSET == 0) {
        for (int i=0; i<params.number_of_antennas; i++) {
            GainOffset_TB_ch.push_back( 1. );
        }
    }


}


inline void Detector::ReadThres_TestBed( string filename, Settings *settings1){
  if (settings1->TRIG_THRES_MODE != 1){
    for (int N=0; N<16; N++) {
      Thres_TB_ch.push_back(settings1->POWERTHRESHOLD);
      cout<<"Thres ch"<<N<<" : "<<Thres_TB_ch[N]<<endl;
    }
  }
  else if (settings1->TRIG_THRES_MODE == 1){
    ifstream Thres( filename.c_str() );
        
    string line;
        
    int N=0;
        
    if ( Thres.is_open() ) {
      while (Thres.good() ) {
	
	getline (Thres, line);
        
	Thres_TB_ch.push_back( atof( line.c_str() ) );
	cout<<"Thres ch" << N << " : " << Thres_TB_ch[N] << endl;
        
	N++;
        
      }
      Thres.close();
    }
  }
}


inline void Detector::ReadThresOffset_TestBed(string filename, Settings *settings1) {    // will return gain offset (unit in voltage) for different chs

    // no threshold offset (just use 1)
    if ( settings1->TRIG_THRES_MODE == 0 ) {
        for (int N=0; N<16; N++) {
            ThresOffset_TB_ch.push_back(1.);
            cout<<"ThresOffset ch"<<N<<" : "<<ThresOffset_TB_ch[N]<<endl;
        }
    }


    else if ( settings1->TRIG_THRES_MODE == 1 ) {
        for (int N=0; N<16; N++) {
            ThresOffset_TB_ch.push_back(1.);
            cout<<"ThresOffset ch"<<N<<" : "<<ThresOffset_TB_ch[N]<<endl;
        }
    }
    else if ( settings1->TRIG_THRES_MODE == 2 ) {
        for (int N=0; N<16; N++) {
            ThresOffset_TB_ch.push_back( pow(GainOffset_TB_ch[N],2) );
            cout<<"ThresOffset ch"<<N<<" : "<<ThresOffset_TB_ch[N]<<endl;
        }
    }

}


inline void Detector::ReadTemp_TestBed(string filename, Settings *settings1) {    // will return gain offset (unit in voltage) for different chs


    ifstream Temp( filename.c_str() );
    
    string line;
    
    int N=0;
    
    if ( Temp.is_open() ) {
        while (Temp.good() ) {
            
            getline (Temp, line);
            
            Temp_TB_ch.push_back( atof( line.c_str() ) );
            cout<<"System temp ch"<<N<<" : "<<Temp_TB_ch[N]<<endl;
            
            N++;
            
        }
        Temp.close();
    }

}



inline void Detector::ReadRFCM_TestBed(string filename, Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain
    
    ifstream RFCM( filename.c_str() );
    
    string line;
    
    int N=-1;
    
    vector <double> xfreq_tmp;
    vector <double> ygain_tmp;
    
    if ( RFCM.is_open() ) {
        while (RFCM.good() ) {
            
            getline (RFCM, line);
            xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() )*1.e-6 ); // from Hz to MHz
            
            ygain_tmp.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) );
            
            N++;
        }
        RFCM.close();
    }
    
    else cout<<"RFCM file can not opened!!"<<endl;
    
    double xfreq[N], ygain[N];  // need array for Tools::SimpleLinearInterpolation
    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double ygain_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    for (int i=0;i<N;i++) { // copy values
        xfreq[i] = xfreq_tmp[i];
        ygain[i] = ygain_tmp[i];
    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    

    // check if there's pre assigned chs
    int ch_no = RFCM_TB_databin_ch.size();
    
    // Tools::SimpleLinearInterpolation will return RFCM array (in dB)
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, freq_step, Freq, RFCM_TB_ch[ch_no] );
    
    Tools::SimpleLinearInterpolation( N, xfreq, ygain, settings1->DATA_BIN_SIZE/2, xfreq_databin, ygain_databin );


    // set vector array size to number of chs
    RFCM_TB_databin_ch.resize(ch_no+1);
    
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
        RFCM_TB_databin_ch[ch_no].push_back( ygain_databin[i] );
    }
    
    
    
    
}



void Detector::ReadRFCM_New(Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain

    // We can use FilterGain array as a original array

    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double ygain_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );

    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    int RFCM_ch = RFCM_TB_databin_ch.size();

    for (int ch=0; ch<RFCM_ch; ch++) {
        Tools::SimpleLinearInterpolation( freq_step, Freq, RFCM_TB_ch[ch], settings1->DATA_BIN_SIZE/2, xfreq_databin, ygain_databin );
            
        RFCM_TB_databin_ch[ch].clear();
        
        for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
            RFCM_TB_databin_ch[ch].push_back( ygain_databin[i] );
        }
    }

}


            
inline void Detector::ReadRayleighFit_TestBed(string filename, Settings *settings1, int ch_no) {    // will read Rayleigh fit result from the file

    ifstream Rayleigh_file( filename.c_str() );
    
    string line;
    
    //int N=-1;
    int init = 1;
    int ch_loop = 0;
    
    vector <double> xfreq_tmp;
    vector <vector <double> > fit_tmp; // 2d array for ch
    fit_tmp.resize(ch_no);
    int ch_tmp;
    double freq_tmp_tmp;
    
    if ( Rayleigh_file.is_open() ) {
        while (Rayleigh_file.good() ) {
            
            if (init == 1) { // ok, skip first line
                getline (Rayleigh_file, line);
                init++;
            }
            else { // from second line, read
                
                //getline (Rayleigh_file, line);
                getline (Rayleigh_file, line, ',');

                //xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) ); // freq in MHz
                //line_no1 = line.find_first_of(",");
                //freq_tmp_tmp = atof( line.substr(0, line.find_first_of(",")).c_str() ); // freq in MHz
                freq_tmp_tmp = atof( line.c_str() ); // freq in MHz

                getline (Rayleigh_file, line, ',');

                //chan_tmp.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) ); // channel number
                //ch_tmp = atof( line.substr(line.find_first_of(",") + 1).c_str() ); // channel number (skip)
                ch_tmp = atof( line.c_str() ); // channel number (skip)

                getline (Rayleigh_file, line, ',');
                
                //fit_tmp[ch_tmp].push_back( atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() ) ); // fit result
                fit_tmp[ch_tmp].push_back( atof( line.c_str() ) ); // fit result

                if (ch_tmp == 0) xfreq_tmp.push_back( freq_tmp_tmp );

                getline (Rayleigh_file, line, '\n');
                
                //N++;
            }

        }
        Rayleigh_file.close();
    }
    
    else cout<<"Rayleigh file can not opened!!"<<endl;

    RayleighFit_ch = ch_no;

    //int N = (int)xfreq_tmp.size();
    int N = (int)xfreq_tmp.size() - 1;
    
    double xfreq[N];  // need array for Tools::SimpleLinearInterpolation
    double Rayleigh[N];

    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double Rayleigh_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    // now below are values that shared in all channels
    for (int i=0;i<N;i++) { // copy values
        xfreq[i] = xfreq_tmp[i];

        /*
        for (int ch=0; ch<ch_no; ch++) {
            ygain[i] = ygain_tmp[i];
        }
        */
    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    // set vector array size to number of chs
    Rayleigh_TB_databin_ch.resize(ch_no);
    

    // now loop over channels and do interpolation
    for (int ch=0; ch<ch_no; ch++) {

        // copy fit values
        for (int i=0;i<N;i++) {
            Rayleigh[i] = fit_tmp[ch][i];
        }


        // Tools::SimpleLinearInterpolation will return Rayleigh array (in dB)
        Tools::SimpleLinearInterpolation( N, xfreq, Rayleigh, freq_step, Freq, Rayleigh_TB_ch[ch] );
        
        Tools::SimpleLinearInterpolation( N, xfreq, Rayleigh, settings1->DATA_BIN_SIZE/2, xfreq_databin, Rayleigh_databin );

    
        for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
            Rayleigh_TB_databin_ch[ch].push_back( Rayleigh_databin[i] );
        }
    
    }
    
    
    
}




inline void Detector::ReadRayleighFit_TestBed(string filename, Settings *settings1) {    // will read Rayleigh fit result from the file

    ifstream Rayleigh_file( filename.c_str() );
    
    string line;
    
    //int N=-1;
    int init = 1;
    int ch_loop = 0;
    
    vector <double> xfreq_tmp;
    vector <vector <double> > fit_tmp; // 2d array for ch

    fit_tmp.resize(16); // start with max number of chs
    int ch_no=0; // this is actual number of chs from file (will be obtained)
    int total_line = 0;

    int ch_tmp;
    double fit_tmp_tmp;
    double freq_tmp_tmp;

    //cout<<"Reading RayleighFit file!"<<endl;
    
    if ( Rayleigh_file.is_open() ) {
        while (Rayleigh_file.good() ) {
            
            if (init == 1) { // ok, skip first line
                getline (Rayleigh_file, line);
                init++;
            }
            else { // from second line, read


                //getline (Rayleigh_file, line);
                getline (Rayleigh_file, line, ',');

                //xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) ); // freq in MHz
                //line_no1 = line.find_first_of(",");
                //freq_tmp_tmp = atof( line.substr(0, line.find_first_of(",")).c_str() ); // freq in MHz
                freq_tmp_tmp = atof( line.c_str() ); // freq in MHz

                getline (Rayleigh_file, line, ',');

                //chan_tmp.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) ); // channel number
                //ch_tmp = atof( line.substr(line.find_first_of(",") + 1).c_str() ); // channel number (skip)
                ch_tmp = atof( line.c_str() ); // channel number (skip)

                getline (Rayleigh_file, line, ',');
                
                //fit_tmp[ch_tmp].push_back( atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() ) ); // fit result
                fit_tmp_tmp = atof( line.c_str() ); // fit result
                fit_tmp[ch_tmp].push_back( fit_tmp_tmp ); // fit result

                if (ch_tmp == 0) xfreq_tmp.push_back( freq_tmp_tmp );

                getline (Rayleigh_file, line, '\n');


                
                /*
                getline (Rayleigh_file, line);
                //xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) ); // freq in MHz
                freq_tmp_tmp = atof( line.substr(0, line.find_first_of(",")).c_str() ); // freq in MHz

                //chan_tmp.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) ); // channel number
                ch_tmp = atof( line.substr(line.find_first_of(",") + 1).c_str() ); // channel number (skip)
                
                //fit_tmp[ch_tmp].push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) ); // fit result
                fit_tmp_tmp = atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() ); // fit result
                fit_tmp[ch_tmp].push_back( fit_tmp_tmp ); // fit result

                if (ch_tmp == 0) xfreq_tmp.push_back( freq_tmp_tmp );
                */
                
                total_line++;

                //cout<<freq_tmp_tmp<<"\t"<<ch_tmp<<"\t"<<fit_tmp_tmp<<endl;
            }

        }
        Rayleigh_file.close();
    }
    
    else cout<<"Rayleigh file can not opened!!"<<endl;

    //int N = (int)xfreq_tmp.size();
    int N = (int)xfreq_tmp.size() - 1;
    total_line = total_line - 1;
    //cout<<"freq bin : "<<N<<endl;
    //cout<<"Total data lines : "<<total_line<<endl;
    ch_no = total_line / N;

    //cout<<"number of ch from RayleighFit file : "<<ch_no<<endl;

    fit_tmp.resize(ch_no); // now resize (no data part will be removed)

    RayleighFit_ch = ch_no;

    
    double xfreq[N];  // need array for Tools::SimpleLinearInterpolation
    double Rayleigh[N];

    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double Rayleigh_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    // now below are values that shared in all channels
    for (int i=0;i<N;i++) { // copy values
        xfreq[i] = xfreq_tmp[i];

        /*
        for (int ch=0; ch<ch_no; ch++) {
            ygain[i] = ygain_tmp[i];
        }
        */
    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    // set vector array size to number of chs
    Rayleigh_TB_databin_ch.resize(ch_no);
    

    // now loop over channels and do interpolation
    for (int ch=0; ch<ch_no; ch++) {

        // copy fit values
        for (int i=0;i<N;i++) {
            Rayleigh[i] = fit_tmp[ch][i];
        }


        // Tools::SimpleLinearInterpolation will return Rayleigh array (in dB)
        Tools::SimpleLinearInterpolation( N, xfreq, Rayleigh, freq_step, Freq, Rayleigh_TB_ch[ch] );
        
        Tools::SimpleLinearInterpolation( N, xfreq, Rayleigh, settings1->DATA_BIN_SIZE/2, xfreq_databin, Rayleigh_databin );

    
        for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
            Rayleigh_TB_databin_ch[ch].push_back( Rayleigh_databin[i] );
        }
    
    }
    
    
    
}


inline void Detector::ReadRayleighFit(string filename, Settings *settings1) {    // will read Rayleigh fit result from the file
    
    ifstream Rayleigh_file( filename.c_str() );
    
    string line;
    
    //int N=-1;
    int init = 1;
    int ch_loop = 0;
    
    vector <double> xfreq_tmp;
    vector <vector <double> > fit_tmp; // 2d array for ch
    
    fit_tmp.resize(16); // start with max number of chs
    int ch_no=0; // this is actual number of chs from file (will be obtained)
    int total_line = 0;
    
    int ch_tmp;
    double fit_tmp_tmp;
    double freq_tmp_tmp;
    
    //cout<<"Reading RayleighFit file!"<<endl;
    
    if ( Rayleigh_file.is_open() ) {
        while (Rayleigh_file.good() ) {
            
            if (init == 1) { // ok, skip first line
                getline (Rayleigh_file, line);
                init++;
            }
            else { // from second line, read
                
                
                //getline (Rayleigh_file, line);
                getline (Rayleigh_file, line, ',');
                
                //xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) ); // freq in MHz
                //line_no1 = line.find_first_of(",");
                //freq_tmp_tmp = atof( line.substr(0, line.find_first_of(",")).c_str() ); // freq in MHz
                freq_tmp_tmp = atof( line.c_str() ); // freq in MHz
                
                getline (Rayleigh_file, line, ',');
                
                //chan_tmp.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) ); // channel number
                //ch_tmp = atof( line.substr(line.find_first_of(",") + 1).c_str() ); // channel number (skip)
                ch_tmp = atof( line.c_str() ); // channel number (skip)
                
                getline (Rayleigh_file, line, ',');
                
                //fit_tmp[ch_tmp].push_back( atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() ) ); // fit result
                fit_tmp_tmp = atof( line.c_str() ); // fit result
                fit_tmp[ch_tmp].push_back( fit_tmp_tmp ); // fit result
                
                if (ch_tmp == 0) xfreq_tmp.push_back( freq_tmp_tmp );
                
                getline (Rayleigh_file, line, '\n');
                
                
                
                /*
                 getline (Rayleigh_file, line);
                 //xfreq_tmp.push_back( atof( line.substr(0, line.find_first_of(",")).c_str() ) ); // freq in MHz
                 freq_tmp_tmp = atof( line.substr(0, line.find_first_of(",")).c_str() ); // freq in MHz
                 
                 //chan_tmp.push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) ); // channel number
                 ch_tmp = atof( line.substr(line.find_first_of(",") + 1).c_str() ); // channel number (skip)
                 
                 //fit_tmp[ch_tmp].push_back( atof( line.substr(line.find_first_of(",") + 1).c_str() ) ); // fit result
                 fit_tmp_tmp = atof( line.substr( line.find_first_of("=") + 1, line.find_first_of(",") ).c_str() ); // fit result
                 fit_tmp[ch_tmp].push_back( fit_tmp_tmp ); // fit result
                 
                 if (ch_tmp == 0) xfreq_tmp.push_back( freq_tmp_tmp );
                 */
                
                total_line++;
                
                //cout<<freq_tmp_tmp<<"\t"<<ch_tmp<<"\t"<<fit_tmp_tmp<<endl;
            }
            
        }
        Rayleigh_file.close();
    }
    
    else cout<<"Rayleigh file can not opened!!"<<endl;
    
    //int N = (int)xfreq_tmp.size();
    int N = (int)xfreq_tmp.size() - 1;
    total_line = total_line - 1;
    //cout<<"freq bin : "<<N<<endl;
    //cout<<"Total data lines : "<<total_line<<endl;
    ch_no = total_line / N;
    
    //cout<<"number of ch from RayleighFit file : "<<ch_no<<endl;
    
    fit_tmp.resize(ch_no); // now resize (no data part will be removed)
    
    RayleighFit_ch = ch_no;
    
    
    double xfreq[N];  // need array for Tools::SimpleLinearInterpolation
    double Rayleigh[N];
    
    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double Rayleigh_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );
    
    // now below are values that shared in all channels
    for (int i=0;i<N;i++) { // copy values
        xfreq[i] = xfreq_tmp[i];
        
        /*
         for (int ch=0; ch<ch_no; ch++) {
         ygain[i] = ygain_tmp[i];
         }
         */
    }
    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }
    // set vector array size to number of chs
    Rayleigh_TB_databin_ch.resize(ch_no);
    
    
    // now loop over channels and do interpolation
    for (int ch=0; ch<ch_no; ch++) {
        
        // copy fit values
        for (int i=0;i<N;i++) {
            Rayleigh[i] = fit_tmp[ch][i];
        }
        
        
        // Tools::SimpleLinearInterpolation will return Rayleigh array (in dB)
        Tools::SimpleLinearInterpolation( N, xfreq, Rayleigh, freq_step, Freq, Rayleigh_TB_ch[ch] );
        
        Tools::SimpleLinearInterpolation( N, xfreq, Rayleigh, settings1->DATA_BIN_SIZE/2, xfreq_databin, Rayleigh_databin );
        
        
        for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
            Rayleigh_databin_ch[ch].push_back( Rayleigh_databin[i] );
        }
    }
}



void Detector::ReadNoiseFig_New(Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain

    // We can use FilterGain array as a original array
cout<<"In ReadNoiseFig_New"<<endl;

    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double NoiseFig_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );

    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    int NoiseFig_numCh = NoiseFig_databin_ch.size();
cout<<"NoiseFig_numCh: "<<NoiseFig_numCh<<endl;

    for (int ch=0; ch<NoiseFig_numCh; ch++) {

        Tools::SimpleLinearInterpolation( freq_step, Freq, NoiseFig_ch[ch], settings1->DATA_BIN_SIZE/2, xfreq_databin, NoiseFig_databin );
            
        NoiseFig_databin_ch[ch].clear();
        
        for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
            if(i%100==0) cout<<"NoiseFig_databin: "<<NoiseFig_databin[i]<<endl;
            NoiseFig_databin_ch[ch].push_back( NoiseFig_databin[i] );
        }
    }

}




void Detector::ReadRayleigh_New(Settings *settings1) {    // will return gain (dB) with same freq bin with antenna gain

    // We can use FilterGain array as a original array

    double xfreq_databin[settings1->DATA_BIN_SIZE/2];   // array for FFT freq bin
    double Rayleigh_databin[settings1->DATA_BIN_SIZE/2];   // array for gain in FFT bin
    double df_fft;
    
    df_fft = 1./ ( (double)(settings1->DATA_BIN_SIZE) * settings1->TIMESTEP );

    for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {    // this one is for DATA_BIN_SIZE
        xfreq_databin[i] = (double)i * df_fft / (1.E6); // from Hz to MHz
    }

    int Rayleigh_ch = Rayleigh_TB_databin_ch.size();

    for (int ch=0; ch<Rayleigh_ch; ch++) {

        Tools::SimpleLinearInterpolation( freq_step, Freq, Rayleigh_TB_ch[ch], settings1->DATA_BIN_SIZE/2, xfreq_databin, Rayleigh_databin );
            
        Rayleigh_TB_databin_ch[ch].clear();
        
        for (int i=0;i<settings1->DATA_BIN_SIZE/2;i++) {
            Rayleigh_TB_databin_ch[ch].push_back( Rayleigh_databin[i] );
        }
    }

}






double Detector::GetGainOffset( int StationID, int ch, Settings *settings1 ) { // returns voltage factor for specific channel gain off set

    if ( (StationID == 0)&&(settings1->DETECTOR==3) ) { // if TestBed, we have offset values
        return GainOffset_TB_ch[ch];
    }
    else if ( settings1->TRIG_ONLY_LOW_CH_ON==1 ) { // if only triggered by bottom chs with TB info
        return GainOffset_TB_ch[ch];
    }
    else {
        return 1.;// other stations, just return 1.
    }
}


double Detector::GetThresOffset( int StationID, int ch, Settings *settings1 ) { // returns voltage factor for specific channel gain off set

    if ( (StationID == 0)&&(settings1->DETECTOR==3) ) { // if TestBed, we have offset values
        return ThresOffset_TB_ch[ch];
    }
    else if ( settings1->TRIG_ONLY_LOW_CH_ON==1 ) { // if only triggered by bottom chs with TB info
        return ThresOffset_TB_ch[ch];
    }
    else {
        return 1.;// other stations, just return 1.
    }
}

double Detector::GetThres( int StationID, int ch, Settings *settings1 ){
  if ( (StationID == 0) && (settings1->DETECTOR==3) ) {
    //    cout << "TBChannel:Threshold: " << ch << " : " << Thres_TB_ch[ch] << endl;
    return Thres_TB_ch[ch];
  }

  else if ( settings1->TRIG_ONLY_LOW_CH_ON==1 ) { // if only triggered by bottom chs with TB info
    //    cout << "TBChannel:Threshold: " << ch << " : " << Thres_TB_ch[ch] << endl;
    return Thres_TB_ch[ch];
  }

  else {
    return settings1->POWERTHRESHOLD;
  }
}


        
double Detector::GetTemp( int StationID, int ch, Settings *settings1 ) {  // returns system temp for specific channel

    if ( (StationID == 0)&&(settings1->DETECTOR==3) ) { // if TestBed, we have offset values
        return Temp_TB_ch[ch];
    }
    else if ( settings1->TRIG_ONLY_LOW_CH_ON==1 ) { // if only triggered by bottom chs with TB info
        return Temp_TB_ch[ch];
    }
    else {
        return settings1->NOISE_TEMP;// other stations, just return setup NOISE_TEMP
    }
}



/*
void Detector::ShiftLocationFlat(int stationNum){
    for (int i = 0; i < int(stations[stationNum].strings.size()); i++){
        for (int j = 0; j < int(stations[stationNum].strings[i].antennas.size()); j++){
            stations.[stationNum].strings[i].antennas[j].SetX(stations.[stationNum].strings[i].antennas[j].GetX() + 0);
            stations.[stationNum].strings[i].antennas[j].SetY(stations.[stationNum].strings[i].antennas[j].GetY() + 0);
            stations.[stationNum].strings[i].antennas[j].SetZ(stations.[stationNum].strings[i].antennas[j].GetZ() + 0);
        }
    }
    for (int i = 0; i < int(stations[stationNum].surface.size()); i++){
            stations.[stationNum].surface[i].SetX(stations.[stationNum].surface[i].GetX() + 0);
            stations.[stationNum].surface[i].SetY(stations.[stationNum].surface[i].GetY() + 0);
            stations.[stationNum].surface[i].SetZ(stations.[stationNum].surface[i].GetZ() + 0);
    }
    
}
*/

void Detector::getDiodeModel(Settings *settings1) {
    
    
    //  this is our homegrown diode response function which is a downgoing gaussian followed by an upward step function
    TF1 *fdown1=new TF1("fl_down1","[3]+[0]*exp(-1.*(x-[1])*(x-[1])/(2*[2]*[2]))",-300.E-9,300.E-9);
    fdown1->SetParameter(0,-0.8);
    //  fdown1->SetParameter(1,15.E-9);
    fdown1->SetParameter(1,15.E-9);
    fdown1->SetParameter(2,2.3E-9);
    //fdown1->SetParameter(2,0.5E-9);
    fdown1->SetParameter(3,0.);
    
    TF1 *fdown2=new TF1("fl_down2","[3]+[0]*exp(-1.*(x-[1])*(x-[1])/(2*[2]*[2]))",-300.E-9,300.E-9);
    fdown2->SetParameter(0,-0.2);
    //  fdown2->SetParameter(1,15.E-9);
    fdown2->SetParameter(1,15.E-9);
    fdown2->SetParameter(2,4.0E-9);
    //fdown2->SetParameter(2,0.5E-9);
    fdown2->SetParameter(3,0.);
    
    /*
     // commented for 5 different banding as in ARA, we only need full band
     maxt_diode=70.E-9;
     idelaybeforepeak[0]=(int)(5.E-9/TIMESTEP);
     iwindow[0]=(int)(20.E-9/TIMESTEP);
     idelaybeforepeak[1]=(int)(5.E-9/TIMESTEP);
     iwindow[1]=(int)(20.E-9/TIMESTEP);
     idelaybeforepeak[2]=(int)(5.E-9/TIMESTEP);
     iwindow[2]=(int)(20.E-9/TIMESTEP);
     idelaybeforepeak[3]=(int)(5.E-9/TIMESTEP);
     iwindow[3]=(int)(20.E-9/TIMESTEP);
     idelaybeforepeak[4]=(int)(13.E-9/TIMESTEP);
     iwindow[4]=(int)(4.E-9/TIMESTEP);
     */
    
    //maxt_diode=70.E-9;
    //idelaybeforepeak=(int)(13.E-9/TIMESTEP);
    //iwindow=(int)(4.E-9/TIMESTEP);
    
    maxt_diode= settings1->MAXT_DIODE;
    maxt_diode_bin = (int)( maxt_diode / TIMESTEP );
    idelaybeforepeak= settings1->IDELAYBEFOREPEAK_DIODE;
    iwindow= settings1->IWINDOW_DIODE;
    ibinshift = NFOUR/4 - (int)( maxt_diode / TIMESTEP );
    
    //fdown1->Copy(fdiode);
    
    TF1 *f_up=new TF1("f_up","[0]*([3]*(x-[1]))^2*exp(-(x-[1])/[2])",-200.E-9,100.E-9);
    
    f_up->SetParameter(2,7.0E-9);
    f_up->SetParameter(0,1.);
    f_up->SetParameter(1,18.E-9);
    f_up->SetParameter(3,1.E9);
    
    
    double sum=0.;
	
    f_up->SetParameter(0,-1.*sqrt(2.*PI)*(fdown1->GetParameter(0)*fdown1->GetParameter(2)+fdown2->GetParameter(0)*fdown2->GetParameter(2))/(2.*pow(f_up->GetParameter(2),3.)*1.E18));
	
    for (int i=0;i<NFOUR/2;i++) {
        
        diode_real.push_back(0.);   // first puchback 0. value  (this is actually not standard way though works fine)
	    
        //if (time[i]>0. && time[i]<maxt_diode) {
        if (i<(int)(maxt_diode/TIMESTEP)) { // think this is same with above commented if
            
            diode_real[i]=fdown1->Eval((double)i*TIMESTEP)+fdown2->Eval((double)i*TIMESTEP);
            if (i>(int)(f_up->GetParameter(1)/TIMESTEP))
                diode_real[i]+=f_up->Eval((double)i*TIMESTEP);
            
            sum+=diode_real[i];
        }
        /*
         // as we set default as 0 above, we dont need to set 0 with extra step
         else {
         diode_real[i]=0.;  
         } 
         */
    }
    
    //cout<<"done settings diode_real arrays"<<endl;
    
    
    // diode_real is the time domain response of the diode
    //
    // now get f domain response with realft
    
    double diode_real_fft[settings1->DATA_BIN_SIZE*2];  // double sized array for myconvlv
    //double diode_real_fft[settings1->DATA_BIN_SIZE + 512];  // DATA_BIN_SIZE + 512 bin (zero padding) for myconvlv
    double diode_real_fft_half[NFOUR];    // double sized array for NFOUR/2
    double diode_real_fft_double[NFOUR*2];    // test with NFOUR*2 array
    
    
    //for (int i=0; i<settings1->DATA_BIN_SIZE + 512; i++) {  // 512 bin added for zero padding
    for (int i=0; i<settings1->DATA_BIN_SIZE*2; i++) {  // 512 bin added for zero padding
        if ( i<(int)(maxt_diode/TIMESTEP) ) {
            diode_real_fft[i] = diode_real[i];
        }
        else {
            diode_real_fft[i] = 0.;
        }
        
    }
    
    
    for (int i=0; i<NFOUR; i++) {
        if ( i<(int)(maxt_diode/TIMESTEP) ) {
            diode_real_fft_half[i] = diode_real[i];
        }
        else {
            diode_real_fft_half[i] = 0.;
        }
    }
    
    
    // test for double size array
    for (int i=0; i<NFOUR*2; i++) {
        if ( i<(int)(maxt_diode/TIMESTEP) ) {
            diode_real_fft_double[i] = diode_real[i];
        }
        else {
            diode_real_fft_double[i] = 0.;
        }
    }
    
    
    //cout<<"start realft diode_real_fft"<<endl;
    
    // forward FFT
    //Tools::realft(diode_real_fft,1,settings1->DATA_BIN_SIZE+512);
    Tools::realft(diode_real_fft,1,settings1->DATA_BIN_SIZE*2);
    
    // forward FFT for half size array
    Tools::realft(diode_real_fft_half,1,NFOUR);
    
    // forward FFT for double size array
    Tools::realft(diode_real_fft_double,1,NFOUR*2);
    
    
    //cout<<"done realft diode_real_fft"<<endl;
    
    
    fdiode_real_databin.clear();
    fdiode_real.clear();
    fdiode_real_double.clear();
    
    // save f domain diode response in fdiode_real
    //for (int i=0; i<settings1->DATA_BIN_SIZE+512; i++) {
    for (int i=0; i<settings1->DATA_BIN_SIZE*2; i++) {
        fdiode_real_databin.push_back( diode_real_fft[i] );
    }
    
    
    // save f domain diode response in fdiode_real_half
    //for (int i=0; i<NFOUR/2; i++) {
    for (int i=0; i<NFOUR; i++) {
        fdiode_real.push_back( diode_real_fft_half[i] );
    }
    
    // save f domain diode response in fdiode_real_double
    //for (int i=0; i<NFOUR; i++) {
    for (int i=0; i<NFOUR*2; i++) {
        fdiode_real_double.push_back( diode_real_fft_double[i] );
    }
    
    
    
}



// this is a test version for getting new noise waveforms for each event
// for a best performance, we can just set a new reasonable DATA_BIN_SIZE and make new values for those
void Detector::get_NewDiodeModel(Settings *settings1) {

    double diode_real_fft[settings1->DATA_BIN_SIZE*2];  // double sized array for myconvlv

    for (int i=0; i<settings1->DATA_BIN_SIZE*2; i++) {  // 512 bin added for zero padding
        if ( i<(int)(maxt_diode/TIMESTEP) ) {
            diode_real_fft[i] = diode_real[i];
        }
        else {
            diode_real_fft[i] = 0.;
        }
        
    }

    // forward FFT
    Tools::realft(diode_real_fft,1,settings1->DATA_BIN_SIZE*2);

    // clear previous data as we need new diode response array for new DATA_BIN_SIZE
    fdiode_real_databin.clear();

    // save f domain diode response in fdiode_real
    //for (int i=0; i<settings1->DATA_BIN_SIZE+512; i++) {
    for (int i=0; i<settings1->DATA_BIN_SIZE*2; i++) {
        fdiode_real_databin.push_back( diode_real_fft[i] );
    }

}



void Detector::PrepareVectorsInstalled(){
    ARA_station temp_station;
    Antenna_string temp;
    Antenna_string temp_string;
    Antenna temp_antenna;
    Surface_antenna temp_surface;
    
    // prepare vectors
    for (int i=0; i<params.number_of_stations; i++) {
        stations.push_back(temp_station);
        
        for (int j = 0; j < InstalledStations[i].nSurfaces; j++) {
            stations[i].surfaces.push_back(temp_surface);
        }
        
        for (int k = 0; k < InstalledStations[i].nStrings; k++) {
            
            stations[i].strings.push_back(temp_string);
            
            for (int l = 0; l < InstalledStations[i].VHChannel[k].size(); l++){
                stations[i].strings[k].antennas.push_back(temp_antenna);
            }
        }
    }
    
    
}

void Detector::PrepareVectorsInstalled(int importedStation) {

  ARA_station temp_station;
  Antenna_string temp;
  Antenna_string temp_string;
  Antenna temp_antenna;
  Surface_antenna temp_surface;
  
  // prepare vectors
  for (int i=0; i<params.number_of_stations; i++) {
    stations.push_back(temp_station);
    
    for (int j = 0; j < InstalledStations[importedStation].nSurfaces; j++) {
      stations[i].surfaces.push_back(temp_surface);
    }
    
    for (int k = 0; k < InstalledStations[importedStation].nStrings; k++) {
      
      stations[i].strings.push_back(temp_string);
      
      for (int l = 0; l < InstalledStations[importedStation].VHChannel[k].size(); l++){
	stations[i].strings[k].antennas.push_back(temp_antenna);
      }
    }
  }
}


void Detector::SetupInstalledStations(){
 
    int number_of_installed_stations = 4;
    
    InstalledStations.resize(number_of_installed_stations);
    
    std::vector < int > Antennas;
    
    if (InstalledStations.size() > 0){ // Testbed
        
        Antennas.push_back(4);Antennas.push_back(1);
        InstalledStations[0].VHChannel.push_back(Antennas); // Make string 0
        Antennas.clear();
        Antennas.push_back(2);Antennas.push_back(7);
        InstalledStations[0].VHChannel.push_back(Antennas); // Make string 1
        Antennas.clear();
        Antennas.push_back(6);Antennas.push_back(3);
        InstalledStations[0].VHChannel.push_back(Antennas); // Make string 2
        Antennas.clear();
        Antennas.push_back(5);Antennas.push_back(8);
        InstalledStations[0].VHChannel.push_back(Antennas); // Make string 3
        Antennas.clear();
        Antennas.push_back(12);Antennas.push_back(9);
        InstalledStations[0].VHChannel.push_back(Antennas); // Make string 4
        Antennas.clear();
        Antennas.push_back(14);Antennas.push_back(13);
        InstalledStations[0].VHChannel.push_back(Antennas); // Make string 5
        Antennas.clear();
        Antennas.push_back(10);
        InstalledStations[0].VHChannel.push_back(Antennas); // Make string 6
        Antennas.clear();
        Antennas.push_back(11);
        InstalledStations[0].VHChannel.push_back(Antennas); // Make string 7
        Antennas.clear();
        
        InstalledStations[0].nStrings = InstalledStations[0].VHChannel.size();

        InstalledStations[0].surfaceChannels.push_back(15);
        InstalledStations[0].surfaceChannels.push_back(16); 
        
        InstalledStations[0].nSurfaces = InstalledStations[0].surfaceChannels.size();
        
        InstalledStations[0].nChannels = 16;
        InstalledStations[1].nChannelsVH = 14;

    }
    
    if ( InstalledStations.size() > 1 ){ // Station 1 
      //        Antennas.push_back(5);Antennas.push_back(9);Antennas.push_back(1);Antennas.push_back(17);
      Antennas.push_back(5);Antennas.push_back(9);Antennas.push_back(1);Antennas.push_back(13);
      InstalledStations[1].VHChannel.push_back(Antennas); // Make string 0
      Antennas.clear();
      //        Antennas.push_back(6);Antennas.push_back(10);Antennas.push_back(2);Antennas.push_back(18);
      Antennas.push_back(6);Antennas.push_back(10);Antennas.push_back(2);Antennas.push_back(14);
      InstalledStations[1].VHChannel.push_back(Antennas); // Make string 1
      Antennas.clear();
      //        Antennas.push_back(7);Antennas.push_back(11);Antennas.push_back(3);Antennas.push_back(19);
      Antennas.push_back(7);Antennas.push_back(11);Antennas.push_back(3);Antennas.push_back(15);
      InstalledStations[1].VHChannel.push_back(Antennas); // Make string 2
      Antennas.clear();
      //        Antennas.push_back(8);Antennas.push_back(12);Antennas.push_back(4);Antennas.push_back(20);
      Antennas.push_back(4);Antennas.push_back(8);Antennas.push_back(0);Antennas.push_back(12);
      InstalledStations[1].VHChannel.push_back(Antennas); // Make string 3
      Antennas.clear();
      
      InstalledStations[1].nStrings = InstalledStations[1].VHChannel.size();
      
      /*
        InstalledStations[1].surfaceChannels.push_back(13);
        InstalledStations[1].surfaceChannels.push_back(14);
        InstalledStations[1].surfaceChannels.push_back(15);
        InstalledStations[1].surfaceChannels.push_back(16);
      */
      
      InstalledStations[1].surfaceChannels.push_back(17);
      InstalledStations[1].surfaceChannels.push_back(18);
      InstalledStations[1].surfaceChannels.push_back(19);
      InstalledStations[1].surfaceChannels.push_back(20);
      
      InstalledStations[1].nSurfaces = InstalledStations[1].surfaceChannels.size();
        
      InstalledStations[1].nChannels = 20;
      InstalledStations[1].nChannelsVH = 16;
      
    }
    
    if ( InstalledStations.size() > 2 ){ // Station 2
	
	//        Antennas.push_back(5);Antennas.push_back(9);Antennas.push_back(1);Antennas.push_back(17);
        Antennas.push_back(5);Antennas.push_back(13);Antennas.push_back(1);Antennas.push_back(9);
        InstalledStations[2].VHChannel.push_back(Antennas); // Make string 0
        Antennas.clear();
	//        Antennas.push_back(6);Antennas.push_back(10);Antennas.push_back(2);Antennas.push_back(18);
        Antennas.push_back(6);Antennas.push_back(14);Antennas.push_back(2);Antennas.push_back(10);
        InstalledStations[2].VHChannel.push_back(Antennas); // Make string 1
        Antennas.clear();
	//        Antennas.push_back(7);Antennas.push_back(11);Antennas.push_back(3);Antennas.push_back(19);
        Antennas.push_back(7);Antennas.push_back(15);Antennas.push_back(3);Antennas.push_back(11);
        InstalledStations[2].VHChannel.push_back(Antennas); // Make string 2
        Antennas.clear();
	//        Antennas.push_back(8);Antennas.push_back(12);Antennas.push_back(4);Antennas.push_back(20);
        Antennas.push_back(4);Antennas.push_back(12);Antennas.push_back(0);Antennas.push_back(8);
        InstalledStations[2].VHChannel.push_back(Antennas); // Make string 3
        Antennas.clear();
        
        InstalledStations[2].nStrings = InstalledStations[2].VHChannel.size();
        
	/*
	  InstalledStations[1].surfaceChannels.push_back(13);
	  InstalledStations[1].surfaceChannels.push_back(14);
	  InstalledStations[1].surfaceChannels.push_back(15);
	  InstalledStations[1].surfaceChannels.push_back(16);
	*/
	
        InstalledStations[2].surfaceChannels.push_back(16);
        InstalledStations[2].surfaceChannels.push_back(17);
        InstalledStations[2].surfaceChannels.push_back(18);
        InstalledStations[2].surfaceChannels.push_back(19);
	
        InstalledStations[2].nSurfaces = InstalledStations[2].surfaceChannels.size();
        
        InstalledStations[2].nChannels = 20;
        InstalledStations[2].nChannelsVH = 16;
      }


    if ( InstalledStations.size() > 3 ){ // Station 3
	//        Antennas.push_back(5);Antennas.push_back(9);Antennas.push_back(1);Antennas.push_back(17);
        Antennas.push_back(5);Antennas.push_back(13);Antennas.push_back(1);Antennas.push_back(9);
        InstalledStations[3].VHChannel.push_back(Antennas); // Make string 0
        Antennas.clear();
	//        Antennas.push_back(6);Antennas.push_back(10);Antennas.push_back(2);Antennas.push_back(18);
        Antennas.push_back(6);Antennas.push_back(14);Antennas.push_back(2);Antennas.push_back(10);
        InstalledStations[3].VHChannel.push_back(Antennas); // Make string 1
        Antennas.clear();
	//        Antennas.push_back(7);Antennas.push_back(11);Antennas.push_back(3);Antennas.push_back(19);
        Antennas.push_back(7);Antennas.push_back(15);Antennas.push_back(3);Antennas.push_back(11);
        InstalledStations[3].VHChannel.push_back(Antennas); // Make string 2
        Antennas.clear();
	//        Antennas.push_back(8);Antennas.push_back(12);Antennas.push_back(4);Antennas.push_back(20);
        Antennas.push_back(4);Antennas.push_back(12);Antennas.push_back(0);Antennas.push_back(8);
        InstalledStations[3].VHChannel.push_back(Antennas); // Make string 3
        Antennas.clear();
        
        InstalledStations[3].nStrings = InstalledStations[3].VHChannel.size();
        
	/*
	  InstalledStations[1].surfaceChannels.push_back(13);
	  InstalledStations[1].surfaceChannels.push_back(14);
	  InstalledStations[1].surfaceChannels.push_back(15);
	  InstalledStations[1].surfaceChannels.push_back(16);
	*/
	
        InstalledStations[3].surfaceChannels.push_back(16);
        InstalledStations[3].surfaceChannels.push_back(17);
        InstalledStations[3].surfaceChannels.push_back(18);
        InstalledStations[3].surfaceChannels.push_back(19);
	
        InstalledStations[3].nSurfaces = InstalledStations[3].surfaceChannels.size();
        
        InstalledStations[3].nChannels = 20;
        InstalledStations[3].nChannelsVH = 16;
      }
}
    
    
    


/*
void Detector::SetChannelStringAntennaMap()
{
    std::vector < int > Antennas;
    std::vector < vector < int > > Strings;
    
    // Station 0 = Testbed
    Antennas.push_back(4);Antennas.push_back(1);
    Strings.push_back(Antennas); // Make string 0
    Antennas.clear();
    Antennas.push_back(2);Antennas.push_back(7);
    Strings.push_back(Antennas); // Make string 1
    Antennas.clear();
    Antennas.push_back(6);Antennas.push_back(3);
    Strings.push_back(Antennas); // Make string 2
    Antennas.clear();
    Antennas.push_back(5);Antennas.push_back(8);
    Strings.push_back(Antennas); // Make string 3
    Antennas.clear();
    Antennas.push_back(12);Antennas.push_back(9);
    Strings.push_back(Antennas); // Make string 4
    Antennas.clear();
    Antennas.push_back(14);Antennas.push_back(13);
    Strings.push_back(Antennas); // Make string 5
    Antennas.clear();
    Antennas.push_back(10);
    Strings.push_back(Antennas); // Make string 6
    Antennas.clear();
    Antennas.push_back(11);
    Strings.push_back(Antennas); // Make string 7
    Antennas.clear();
    Antennas.push_back(15);    Antennas.push_back(16);
    Strings.push_back(Antennas); // Make string 8 = surface antennas
    Antennas.clear();
    ChannelfromStringAntenna.push_back(Strings);
    Strings.clear();

    //Station 1
    Antennas.push_back(5);Antennas.push_back(9);Antennas.push_back(1);Antennas.push_back(17);
    Strings.push_back(Antennas); // Make string 0
    Antennas.clear();
    Antennas.push_back(6);Antennas.push_back(10);Antennas.push_back(2);Antennas.push_back(18);
    Strings.push_back(Antennas); // Make string 1
    Antennas.clear();
    Antennas.push_back(7);Antennas.push_back(11);Antennas.push_back(3);Antennas.push_back(19);
    Strings.push_back(Antennas); // Make string 2
    Antennas.clear();
    Antennas.push_back(8);Antennas.push_back(12);Antennas.push_back(4);Antennas.push_back(20);
    Strings.push_back(Antennas); // Make string 3
    Antennas.clear();
    Antennas.push_back(13);Antennas.push_back(14);Antennas.push_back(15);Antennas.push_back(16);
    Strings.push_back(Antennas); // Make string 4 = surface antennas
    Antennas.clear();
    ChannelfromStringAntenna.push_back(Strings);
    Strings.clear();
    
    cout << "Stations in Channel - Antenna/string map:" << int(ChannelfromStringAntenna.size()) << endl;
    cout << "Strings in station 0 for Channel - Antenna/string map:" << " : " << int(ChannelfromStringAntenna[0].size()) <<endl;
    cout << "Antennas in station 0 for Channel - Antenna/string map:" << " : " << int(ChannelfromStringAntenna[0][0].size()) <<endl;
    cout << "Strings in station 1 for Channel - Antenna/string map:" << " : " << int(ChannelfromStringAntenna[1].size()) <<endl;
    cout << "Antennas in station 1 for Channel - Antenna/string map:" << " : " << int(ChannelfromStringAntenna[1][0].size()) <<endl;

    
}
*/

int Detector::GetChannelfromStringAntenna ( int stationNum, int stringnum, int antennanum){
    int ChannelNum;
    if (stationNum < int(InstalledStations.size())){
        if (stringnum < int(InstalledStations[stationNum].VHChannel.size())){
            if (antennanum < int(InstalledStations[stationNum].VHChannel[stringnum].size())){
                ChannelNum = InstalledStations[stationNum].VHChannel[stringnum][antennanum];
                return ChannelNum;
            }
            else {
                cerr << "Invalid request for station channel map: antenna number" << endl;
                return -1;
            }
        }
        else {
            cerr << "Invalid request for station channel map: string number" << endl;
            return -1;
        }
    }
    else {
        cerr << "Invalid request for station channel map: station number" << endl;
        cout << stationNum << " : " <<  int(InstalledStations.size()) << endl;
        return -1;
    }
}


// more general used function
int Detector::GetChannelfromStringAntenna ( int stationNum, int stringnum, int antennanum, Settings *settings1) {

    int ChannelNum;

    // for the cases when actual installed TestBed stations geom info is in use
    if ( settings1->DETECTOR==3 ) {
        if (stationNum < int(InstalledStations.size())){
            if (stringnum < int(InstalledStations[stationNum].VHChannel.size())){
                if (antennanum < int(InstalledStations[stationNum].VHChannel[stringnum].size())){
                    ChannelNum = InstalledStations[stationNum].VHChannel[stringnum][antennanum];
                    return ChannelNum;
                }
                else {
                    cerr << "Invalid request for station channel map: antenna number - 3" << endl;
                    return -1;
                }
            }
            else {
                cerr << "Invalid request for station channel map: string number" << endl;
                return -1;
            }
        }
        else {
            cerr << "Invalid request for station channel map: station number" << endl;
            cout << stationNum << " : " <<  int(InstalledStations.size()) << endl;
            return -1;
        }
    }

    else if(settings1->DETECTOR==4){
      int stationId=settings1->DETECTOR_STATION;
      //      cout << settings1->DETECTOR << endl;
      //      int stationId=stationNum;
      if (stationId < int(InstalledStations.size())){
	if (stringnum < int(InstalledStations[stationId].VHChannel.size())){
	  if (antennanum < int(InstalledStations[stationId].VHChannel[stringnum].size())){
	    ChannelNum = InstalledStations[stationId].VHChannel[stringnum][antennanum];
	    return ChannelNum+1;
	  }
	  else {
	    cerr << "Invalid request for station channel map: antenna number - 4" << endl;
	    //cerr << stationId << " : " << stringnum << " : " << antennanum << endl;
	    return -1;
	  }
	}
	else {
	  cerr << "Invalid request for station channel map: string number" << endl;
	  return -1;
	}
      }
      else {
	cerr << "Invalid request for station channel map: station number" << endl;
	cout << stationNum << " : " <<  int(InstalledStations.size()) << endl;
	return -1;
      }


    }
    // if only ideal stations are in use and also installed ARA1a (use ARA1a ch mapping for now)
    else {
      if (stringnum < int(InstalledStations[1].VHChannel.size())){
	if (antennanum < int(InstalledStations[1].VHChannel[stringnum].size())){
	  ChannelNum = InstalledStations[1].VHChannel[stringnum][antennanum];
	  return ChannelNum;
	  //return ChannelNum+1; // Lu 06/24/2017
	}
	else {
	  cerr << "Invalid request for station channel map: antenna number - test" << endl;
	  return -1;
	}
      }
    }
    // we can add new case when mixed installed and ideal stations case later



}



void Detector::GetSSAfromChannel ( int stationID, int channelNum, int * antennaNum, int * stringNum) {
    *stringNum = -1;
    *antennaNum = -1;
        for (int i = 0; i < int(InstalledStations[stationID].VHChannel.size()); i++){
            for (int j = 0; j < int(InstalledStations[stationID].VHChannel[i].size()); j++){
                if (channelNum == InstalledStations[stationID].VHChannel[i][j]){
                    *stringNum = i;
                    *antennaNum = j;
                }
            }
        }
    
    if (*stringNum == -1){
        cerr << "No string/antenna matches the channel number" << endl;
    }
    
    return;
}


void Detector::GetSSAfromChannel ( int stationID, int channelNum, int * antennaNum, int * stringNum, Settings *settings1) {
    *stringNum = -1;
    *antennaNum = -1;

    // for the cases when actual installed TestBed stations geom info is in use
    if ( settings1->DETECTOR==3 ) {
      for (int i = 0; i < int(InstalledStations[stationID].VHChannel.size()); i++){
	for (int j = 0; j < int(InstalledStations[stationID].VHChannel[i].size()); j++){
	  if (channelNum == InstalledStations[stationID].VHChannel[i][j]){
	    *stringNum = i;
	    *antennaNum = j;
	  }
	}
      }
      
      if (*stringNum == -1){
	cerr << "No string/antenna matches the channel number" << endl;
      }
    }
    else if (settings1->DETECTOR==4){

      for (int i = 0; i < int(InstalledStations[stationID].VHChannel.size()); i++){
	for (int j = 0; j < int(InstalledStations[stationID].VHChannel[i].size()); j++){
	  if (channelNum == InstalledStations[stationID].VHChannel[i][j]){
	    *stringNum = i;
	    *antennaNum = j;
	  }
	}
      }
      
      if (*stringNum == -1){
	cerr << "No string/antenna matches the channel number" << endl;
      }

    }
    // if only ideal stations are in use and also installed ARA1a (use ARA1a ch mapping for now)
    else {

        for (int i = 0; i < int(InstalledStations[1].VHChannel.size()); i++){
            for (int j = 0; j < int(InstalledStations[1].VHChannel[i].size()); j++){
                if (channelNum == InstalledStations[1].VHChannel[i][j]){
                    *stringNum = i;
                    *antennaNum = j;
                }
            }
        }
    
        if (*stringNum == -1){
            cerr << "No string/antenna matches the channel number" << endl;
        }
    }
    
    return;
}


#ifdef ARA_UTIL_EXISTS

void Detector::UseAntennaInfo(int stationNum, Settings *settings1){
    
    //AraGeomTool *araGeom=AraGeomTool::Instance();
    AraGeomTool *araGeom = new AraGeomTool();

    if (stationNum == 0) params.TestBed_BH_Mean_delay = 0.;
    //cout<<"No of chs in station "<<stationNum<<" : "<<InstalledStations[stationNum].nChannels+1<<endl;
    
    for ( int chan = 1; chan < InstalledStations[stationNum].nChannels+1; chan++){
        
        double avgX, avgY;
        
        int antennaNum, stringNum;
        //GetSSAfromChannel(stationNum, chan, &antennaNum, &stringNum);
        GetSSAfromChannel(stationNum, chan, &antennaNum, &stringNum, settings1);
	
        if (araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].polType != AraAntPol::kSurface){

            stations[stationNum].strings[stringNum].antennas[antennaNum].SetX(stations[stationNum].GetX()+araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].antLocation[0]);
            stations[stationNum].strings[stringNum].antennas[antennaNum].SetY(stations[stationNum].GetY()+araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].antLocation[1]);
            //stations[stationNum].strings[stringNum].antennas[antennaNum].SetZ(araGeom->fStationInfo[stationNum].fAntInfo[chan-1].antLocation[2]-double(settings1->DEPTH_CHANGE));
            stations[stationNum].strings[stringNum].antennas[antennaNum].SetZ(araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].antLocation[2]);
            
            /*
             cout <<
             "DetectorStation:string:antenna:X:Y:Z:: " <<
             i<< " : " <<
             j<< " : " <<
             k<< " : " <<
             stations[i].strings[j].antennas[k].GetX() << " : " <<
             stations[i].strings[j].antennas[k].GetY() << " : " <<
             stations[i].strings[j].antennas[k].GetZ() << " : " <<
             endl;
             */
            
            stations[stationNum].strings[stringNum].antennas[antennaNum].type = int(araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].polType);  //set polarization to match the deployed information
            
            stations[stationNum].strings[stringNum].SetX(stations[stationNum].GetX()+araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].antLocation[0]);
            stations[stationNum].strings[stringNum].SetY(stations[stationNum].GetY()+araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].antLocation[1]);
            
            
            
            if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                stations[stationNum].strings[stringNum].antennas[antennaNum].orient = 0;
            }
            else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                if ( stringNum==0||stringNum==3 ) {
                    if ( antennaNum==0||antennaNum==1 ) {
                        stations[stationNum].strings[stringNum].antennas[antennaNum].orient = 0;
                    }
                    else {
                        stations[stationNum].strings[stringNum].antennas[antennaNum].orient = 1;
                    }
                }
                else {
                    if ( antennaNum==0||antennaNum==1 ) {
                        stations[stationNum].strings[stringNum].antennas[antennaNum].orient = 1;
                    }
                    else {
                        stations[stationNum].strings[stringNum].antennas[antennaNum].orient = 0;
                    }
                }
                
            } //end orientation selection


            // put DAQ channel type 
            if (araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].daqChanType == AraDaqChanType::kDisconeChan) { // BH chs
                stations[stationNum].strings[stringNum].antennas[antennaNum].DAQchan = 0;
            }
            else if (araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].daqChanType == AraDaqChanType::kBatwingChan) { // not BH chs
                stations[stationNum].strings[stringNum].antennas[antennaNum].DAQchan = 1;
            }

            
            //cout << "Borehole ch: " << chan << " station: " << stationNum << " string: " << stringNum << " ant: " << antennaNum << " X: " << stations[stationNum].strings[stringNum].antennas[antennaNum].GetX() << " Y: " << stations[stationNum].strings[stringNum].antennas[antennaNum].GetY() << " Z: " << stations[stationNum].strings[stringNum].antennas[antennaNum].GetZ() << " Type: " << stations[stationNum].strings[stringNum].antennas[antennaNum].type << endl;
            cout << "Borehole ch: " << chan << " station: " << stationNum << " string: " << stringNum << " ant: " << antennaNum << " Type: " << stations[stationNum].strings[stringNum].antennas[antennaNum].type << endl;



            if (stationNum == 0) {

                //cout<<"TestBed ch"<<chan-1<<" delay : "<<araGeom->fStationInfo[stationNum].fAntInfo[chan-1].debugTotalCableDelay<<endl;
                params.TestBed_Ch_delay[chan-1] = araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].debugTotalCableDelay;
                params.TestBed_Ch_delay_bin[chan-1] = params.TestBed_Ch_delay[chan-1]/(settings1->TIMESTEP * 1.e9); // change TIMESTEP s to ns
                //cout<<"TestBed ch"<<chan-1<<" delay bin : "<<params.TestBed_Ch_delay_bin[chan-1]<<endl;
                if (chan<9) params.TestBed_BH_Mean_delay += params.TestBed_Ch_delay[chan-1];

                // give manual delay time for BH chs (for TestBed)
                if (chan == 2) stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay = 50.; // additional delay in ns
                else if (chan == 3) stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay = -10.; // additional delay in ns
                else if (chan == 4) stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay = 20.; // additional delay in ns
                else if (chan == 5) stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay = 30.; // additional delay in ns
                else if (chan == 6) stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay = 20.; // additional delay in ns
                else if (chan == 7) stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay = -10.; // additional delay in ns
                else if (chan == 8) stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay = 10.; // additional delay in ns
                else stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay = 0.;
            }
            else stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay = 0.;// no manual delay for other stations (not known yet)

            // set manual delay bin
            stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay_bin = stations[stationNum].strings[stringNum].antennas[antennaNum].manual_delay / (settings1->TIMESTEP * 1.e9);





        }// end polarization (antenna type) selection
        else {
            
            int antPolNum = araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].antPolNum;
            // set surface antenna postions
            
            stations[stationNum].surfaces[antPolNum].SetX( stations[stationNum].GetX()+araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].antLocation[0]);
            stations[stationNum].surfaces[antPolNum].SetY( stations[stationNum].GetY()+araGeom->getStationInfo(stationNum)->fAntInfo[chan-1].antLocation[1]);
            
            //                cout << "Surface: " << chan << " : " << stationNum << " : " << stringNum << " : " << antennaNum << " : " << stations[stationNum].surfaces[antPolNum].GetX() << " : " << stations[stationNum].surfaces[antPolNum].GetY() << " : " << stations[stationNum].surfaces[antPolNum].GetZ() << " : " << stations[stationNum].surfaces[antPolNum].type << endl;
            
            
            
        }
        
        

        
        
    } // end channel loop
    
    if (stationNum == 0) {
        params.TestBed_BH_Mean_delay /= 8.;
        params.TestBed_BH_Mean_delay_bin = params.TestBed_BH_Mean_delay/(settings1->TIMESTEP * 1.e9); // change TIMESTEP s to ns
        //cout<<"TestBed Mean BH chs delay : "<<params.TestBed_BH_Mean_delay<<endl;
        //cout<<"TestBed Mean BH chs delay bin : "<<params.TestBed_BH_Mean_delay_bin<<endl;

        params.TestBed_WFtime_offset_ns = -20.;
    }
    
}
#endif

#ifdef ARA_UTIL_EXISTS

void Detector::ImportStationInfo(Settings *settings1, int StationIndex, int StationID){
    
    //AraGeomTool *araGeom=AraGeomTool::Instance();
    AraGeomTool *araGeom = new AraGeomTool();

    //    int stationNum = settings1->DETECTOR_STATION;

    if (StationID == 0) params.TestBed_BH_Mean_delay = 0.;
    //cout<<"No of chs in station "<<stationNum<<" : "<<InstalledStations[stationNum].nChannels+1<<endl;
    
    for ( int chan = 0; chan < InstalledStations[StationID].nChannels; chan++){
        
      int antId;
      if (StationID == 0){ antId = chan+1;} 
      else { antId = chan; }

        double avgX, avgY;
        
        int antennaNum, stringNum;
        //GetSSAfromChannel(stationNum, chan, &antennaNum, &stringNum);
        GetSSAfromChannel(StationID, chan, &antennaNum, &stringNum, settings1);

        if (araGeom->getStationInfo(StationID)->fAntInfo[antId].polType != AraAntPol::kSurface){

            stations[StationIndex].strings[stringNum].antennas[antennaNum].SetX(stations[StationIndex].GetX()+araGeom->getStationInfo(StationID)->fAntInfo[antId].antLocation[0]);
            stations[StationIndex].strings[stringNum].antennas[antennaNum].SetY(stations[StationIndex].GetY()+araGeom->getStationInfo(StationID)->fAntInfo[antId].antLocation[1]);
            //stations[stationNum].strings[stringNum].antennas[antennaNum].SetZ(araGeom->fStationInfo[stationNum].fAntInfo[chan-1].antLocation[2]-double(settings1->DEPTH_CHANGE));
            stations[StationIndex].strings[stringNum].antennas[antennaNum].SetZ(araGeom->getStationInfo(StationID)->fAntInfo[antId].antLocation[2]);
            
            /*
             cout <<
             "DetectorStation:string:antenna:X:Y:Z:: " <<
             i<< " : " <<
             j<< " : " <<
             k<< " : " <<
             stations[i].strings[j].antennas[k].GetX() << " : " <<
             stations[i].strings[j].antennas[k].GetY() << " : " <<
             stations[i].strings[j].antennas[k].GetZ() << " : " <<
             endl;
             */
            
            stations[StationIndex].strings[stringNum].antennas[antennaNum].type = int(araGeom->getStationInfo(StationID)->fAntInfo[antId].polType);  //set polarization to match the deployed information
            
            stations[StationIndex].strings[stringNum].SetX(stations[StationIndex].GetX()+araGeom->getStationInfo(StationID)->fAntInfo[antId].antLocation[0]);
            stations[StationIndex].strings[stringNum].SetY(stations[StationIndex].GetY()+araGeom->getStationInfo(StationID)->fAntInfo[antId].antLocation[1]);
            
            
            
            if ( params.antenna_orientation == 0 ) {    // all borehole antennas facing same x
                stations[StationIndex].strings[stringNum].antennas[antennaNum].orient = 0;
            }
            else if ( params.antenna_orientation == 1 ) {   // borehole antennas one next facing different way
                if ( stringNum==0||stringNum==3 ) {
                    if ( antennaNum==0||antennaNum==1 ) {
                        stations[StationIndex].strings[stringNum].antennas[antennaNum].orient = 0;
                    }
                    else {
                        stations[StationIndex].strings[stringNum].antennas[antennaNum].orient = 1;
                    }
                }
                else {
                    if ( antennaNum==0||antennaNum==1 ) {
                        stations[StationIndex].strings[stringNum].antennas[antennaNum].orient = 1;
                    }
                    else {
                        stations[StationIndex].strings[stringNum].antennas[antennaNum].orient = 0;
                    }
                }
                
            } //end orientation selection


            // put DAQ channel type 
            if (araGeom->getStationInfo(StationID)->fAntInfo[antId].daqChanType == AraDaqChanType::kDisconeChan) { // BH chs
                stations[StationIndex].strings[stringNum].antennas[antennaNum].DAQchan = 0;
            }
            else if (araGeom->getStationInfo(StationID)->fAntInfo[antId].daqChanType == AraDaqChanType::kBatwingChan) { // not BH chs
                stations[StationIndex].strings[stringNum].antennas[antennaNum].DAQchan = 1;
            }

            
            //cout << "Borehole ch: " << chan << " station: " << stationNum << " string: " << stringNum << " ant: " << antennaNum << " X: " << stations[stationNum].strings[stringNum].antennas[antennaNum].GetX() << " Y: " << stations[stationNum].strings[stringNum].antennas[antennaNum].GetY() << " Z: " << stations[stationNum].strings[stringNum].antennas[antennaNum].GetZ() << " Type: " << stations[stationNum].strings[stringNum].antennas[antennaNum].type << endl;
            cout << "Borehole ch: " << chan << " inserted station: " << StationIndex << "station: " << StationID << " string: " << stringNum << " ant: " << antennaNum << " Type: " << stations[StationIndex].strings[stringNum].antennas[antennaNum].type << endl;



            if (StationID == 0) {

                //cout<<"TestBed ch"<<chan-1<<" delay : "<<araGeom->fStationInfo[stationNum].fAntInfo[chan-1].debugTotalCableDelay<<endl;
                params.TestBed_Ch_delay[chan] = araGeom->getStationInfo(StationID)->fAntInfo[antId].debugTotalCableDelay;
                params.TestBed_Ch_delay_bin[chan] = params.TestBed_Ch_delay[chan]/(settings1->TIMESTEP * 1.e9); // change TIMESTEP s to ns
                //cout<<"TestBed ch"<<chan-1<<" delay bin : "<<params.TestBed_Ch_delay_bin[chan-1]<<endl;
                if (chan<8) params.TestBed_BH_Mean_delay += params.TestBed_Ch_delay[chan];

                // give manual delay time for BH chs (for TestBed)
                if (chan == 1) stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay = 50.; // additional delay in ns
                else if (chan == 2) stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay = -10.; // additional delay in ns
                else if (chan == 3) stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay = 20.; // additional delay in ns
                else if (chan == 4) stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay = 30.; // additional delay in ns
                else if (chan == 5) stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay = 20.; // additional delay in ns
                else if (chan == 6) stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay = -10.; // additional delay in ns
                else if (chan == 7) stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay = 10.; // additional delay in ns
                else stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay = 0.;
            }
            else stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay = 0.;// no manual delay for other stations (not known yet)

            // set manual delay bin
            stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay_bin = stations[StationIndex].strings[stringNum].antennas[antennaNum].manual_delay / (settings1->TIMESTEP * 1.e9);

        }// end polarization (antenna type) selection
        else {
            
            int antPolNum = araGeom->getStationInfo(StationID)->fAntInfo[antId].antPolNum;
            // set surface antenna postions
            
            stations[StationIndex].surfaces[antPolNum].SetX( stations[StationIndex].GetX()+araGeom->getStationInfo(StationID)->fAntInfo[antId].antLocation[0]);
            stations[StationIndex].surfaces[antPolNum].SetY( stations[StationIndex].GetY()+araGeom->getStationInfo(StationID)->fAntInfo[antId].antLocation[1]);
            
            //                cout << "Surface: " << chan << " : " << stationNum << " : " << stringNum << " : " << antennaNum << " : " << stations[stationNum].surfaces[antPolNum].GetX() << " : " << stations[stationNum].surfaces[antPolNum].GetY() << " : " << stations[stationNum].surfaces[antPolNum].GetZ() << " : " << stations[stationNum].surfaces[antPolNum].type << endl;
            
        }
        
    } // end channel loop
    
    if (StationID == 0) {
        params.TestBed_BH_Mean_delay /= 8.;
        params.TestBed_BH_Mean_delay_bin = params.TestBed_BH_Mean_delay/(settings1->TIMESTEP * 1.e9); // change TIMESTEP s to ns
        //cout<<"TestBed Mean BH chs delay : "<<params.TestBed_BH_Mean_delay<<endl;
        //cout<<"TestBed Mean BH chs delay bin : "<<params.TestBed_BH_Mean_delay_bin<<endl;

        params.TestBed_WFtime_offset_ns = -20.;
    }
    
}
#endif




void Detector::SetupIdealStations(){
        
    IdealStations.resize(2);
    
    std::vector < int > Antennas;
    
    if (IdealStations.size() > 1){ // Station 1
        int stationID = 1;
        IdealStations[stationID].nChannels = 20;
        IdealStations[stationID].nChannelsVH = 16;
        
        Antennas.push_back(5);Antennas.push_back(9);Antennas.push_back(1);Antennas.push_back(17);
        IdealStations[stationID].VHChannel.push_back(Antennas); // Make string 0
        Antennas.clear();
        Antennas.push_back(6);Antennas.push_back(10);Antennas.push_back(2);Antennas.push_back(18);
        IdealStations[stationID].VHChannel.push_back(Antennas); // Make string 1
        Antennas.clear();
        Antennas.push_back(7);Antennas.push_back(11);Antennas.push_back(3);Antennas.push_back(19);
        IdealStations[stationID].VHChannel.push_back(Antennas); // Make string 2
        Antennas.clear();
        Antennas.push_back(8);Antennas.push_back(12);Antennas.push_back(4);Antennas.push_back(20);
        IdealStations[stationID].VHChannel.push_back(Antennas); // Make string 3
        Antennas.clear();
        
        IdealStations[stationID].nStrings = IdealStations[stationID].VHChannel.size();

        IdealStations[stationID].surfaceChannels.push_back(13);
        IdealStations[stationID].surfaceChannels.push_back(14);
        IdealStations[stationID].surfaceChannels.push_back(15);
        IdealStations[stationID].surfaceChannels.push_back(16);
        
        IdealStations[stationID].nSurfaces = IdealStations[stationID].surfaceChannels.size();
        
        Antennas.push_back(0);Antennas.push_back(1);Antennas.push_back(2);Antennas.push_back(3);
        IdealStations[stationID].VHID.push_back(Antennas); // Make string 0
        Antennas.clear();
        Antennas.push_back(4);Antennas.push_back(5);Antennas.push_back(6);Antennas.push_back(7);
        IdealStations[stationID].VHID.push_back(Antennas); // Make string 1
        Antennas.clear();
        Antennas.push_back(8);Antennas.push_back(9);Antennas.push_back(10);Antennas.push_back(11);
        IdealStations[stationID].VHID.push_back(Antennas); // Make string 2
        Antennas.clear();
        Antennas.push_back(12);Antennas.push_back(13);Antennas.push_back(14);Antennas.push_back(15);
        IdealStations[stationID].VHID.push_back(Antennas); // Make string 3
        Antennas.clear();
                
        IdealStations[stationID].surfaceID.push_back(0);
        IdealStations[stationID].surfaceID.push_back(1);
        IdealStations[stationID].surfaceID.push_back(2);
        IdealStations[stationID].surfaceID.push_back(3);
        
        for (int BHAntID = 0; BHAntID < IdealStations[stationID].nChannelsVH; BHAntID++){
             for (int i = 0; i < IdealStations[stationID].VHID.size(); i++){
                for (int j = 0; j < IdealStations[stationID].VHID[i].size(); j++){
                     if (IdealStations[stationID].VHID[i][j] == BHAntID){
                         IdealStations[stationID].IDString.push_back(i);
                         IdealStations[stationID].IDAntenna.push_back(j);
                     }
                 }
             }
        }
        for (int surfAntID = 0; surfAntID < IdealStations[stationID].nSurfaces; surfAntID++){
            for (int i = 0; i < IdealStations[stationID].surfaceID.size(); i++){
                if (IdealStations[stationID].surfaceID[i] == surfAntID){
                    IdealStations[stationID].IDSurface.push_back(i);
                }
            }
        }
    }    
}

int Detector::getStringfromArbAntID( int stationID, int ant_ID){
    for (int i = 0; i < stations[stationID].strings.size(); i++){
        if (ant_ID < stations[stationID].strings[i].antennas.size()) {
            return i;
        } else {
            ant_ID = ant_ID - stations[stationID].strings[i].antennas.size();
        }
    }
}

int Detector::getAntennafromArbAntID( int stationID, int ant_ID){
    for (int i = 0; i < stations[stationID].strings.size(); i++){
        if (ant_ID < stations[stationID].strings[i].antennas.size()) {
            return ant_ID;
        } else {
            ant_ID = ant_ID - stations[stationID].strings[i].antennas.size();
        }
    }
}

Detector::~Detector() {
    cout<<"Destruct class Detector"<<endl;
}



