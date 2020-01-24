////////////////////////////////////////////////////////////////////////////////////////////////
//class Detector:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef DETECTOR_H
#define DETECTOR_H


//#include "TObject.h"
#include <vector>
#include "Trigger.h"
#include "Position.h"
#include "Vector.h"
//#include "IceModel.h"


using namespace std;

//#include "Trigger.h"
//class Event;
//class Efficiencies;
//class Position;

class Detector;
class IceModel;
class Settings;
class TF1;

//struct Parameters {
//class Parameters : public TObject {
class Parameters {

    private:

    public:
    int stations_per_side;           //number of stations on one side of pentagon.
    double station_spacing;          // space between stations.

    int antenna_orientation;        //antenna orientation setting.

    int number_of_stations;         //total stations
    int number_of_strings_station;  // strings per station
    int number_of_antennas_string;  //antennas per string
    int number_of_surfaces_station;    //surface antennas per station
    int number_of_channels; // number of channels for each regular (non-TestBed) station
    
    int number_of_strings_station_TB; //number of strings in the TestBed
    int number_of_antennas_string_TB; //number of antennas in each string in the Testbed
    int number_of_surfaces_station_TB; //number of surface stations for the TestBed
    int number_of_channels_TB; // number of channels for the TestBed
    
    int num_of_channels[2];
    
    int bore_hole_antenna_layout;   // bore hole antenna layout, 0 : VHVH, 1 : VHV, 2 : VHVV

    int number_of_strings;
    int number_of_antennas;

    double core_x;
    double core_y;
    
    int DeployedStations;
    
        
//    static const int freq_step = 60;
//    static const int ang_step = 2664;
//    static const double freq_width = 16.667;  // this value could be changed when Nec2 condition changes!
//    static const double freq_init = 83.333;  // this value could be changed when Nec2 condition changes!
//    static const int freq_width = 16.667;  // this value could be changed when Nec2 condition changes!
//    static const int freq_init = 83.333;  // this value could be changed when Nec2 condition changes!
//    static const float freq_width = 16.667;  // this value could be changed when Nec2 condition changes!
//    static const float freq_init = 83.333;  // this value could be changed when Nec2 condition changes!

//    int freq_step = 60;
//    int ang_step = 2664;
//    double freq_width = 16.667;  // this value could be changed when Nec2 condition changes!
//    double freq_init = 83.333;  // this value could be changed when Nec2 condition changes!

    int freq_step;  // this value will be obtained through settings class and copy to Detector->freq_step
    int ang_step;   // also copy to Detector->ang_step
    double freq_width;  // this value could be changed when Nec2 condition changes! copy to Detector->freq_width
    double freq_init;  // this value could be changed when Nec2 condition changes! copy to Detector->freq_init


    double TestBed_Ch_delay[16];
    int TestBed_Ch_delay_bin[16]; // in bin
    double TestBed_BH_Mean_delay;
    int TestBed_BH_Mean_delay_bin; // in bin

    double TestBed_WFtime_offset_ns; // waveform time offset to match TestBed data waveform


    ClassDef(Parameters,1);


};





//struct Surface_antenna : public Position, public TObject {
//class Surface_antenna : public Position, public TObject {
class Surface_antenna : public Position {

    public:

//    double x, y;
    int type;   // need to be defined
    int orient; //0 : facing x, 1 : y, 2 : -x, 3 : -y.

    double GetG(Detector *D, double freq, double theta, double phi);    // read gain value from Detector class

    ClassDef(Surface_antenna,1);

};


    
//struct Antenna : public Position, public TObject {
//class Antenna : public Position, public TObject {
class Antenna : public Position {
    public:

//    double z;
    int type;  // type 0 : v-pol (bicone), type 1 : h-pol (bowtie for testbed, QSC for ARA)
    int orient; // 0 : facing x, 1 : y, 2 : -x, 3 : -y.

    int DAQchan;    // DAQ channel type (from AraGeomTools). 0 : discone (BH chs), 1 : BAT (shallow, surf; not first 8 chs)

    int manual_delay;   // to fit the waveform to actual TestBed waveform, added manual delay time

    int manual_delay_bin;   // to fit the waveform to actual TestBed waveform, added manual delay bin
    
    double GetG(Detector *D, double freq, double theta, double phi);    // read gain value from Detector class, return 2-D interpolated value


    //ClassDef(Antenna,1);
    ClassDef(Antenna,3);
    
};



//struct Antenna_string : public Position, public TObject {
//class Antenna_string : public Position, public TObject {
class Antenna_string : public Position {
//    double x, y;
//    int number_of_antennas;
    public:
    vector <Antenna> antennas;

    ClassDef(Antenna_string,1);
};

//struct ARA_station : public Position, public TObject {
//class ARA_station : public Position, public TObject {
class ARA_station : public Position {
//    double x, y;
    public:
    vector <Antenna_string> strings;
    vector <Surface_antenna> surfaces;
    int StationID;
    double TRIG_WINDOW; // in ns, the size of the trigger window used for the antennas
    int NFOUR; // 2 X nbins for readout waveform - for fourier tranform
    double TIMESTEP; // trigger and readout timestep
    double DATA_BIN_SIZE;

    int number_of_antennas; // total number of antennas for each stations
    
    ClassDef(ARA_station,1);
};



class InstalledStation {

    public:
    int nSurfaces;
    int nStrings;
    vector < int > surfaceChannels;
    vector < vector < int > > VHChannel;
    int nChannels;
    int nChannelsVH;
    vector < vector < int > > VHID;
    vector < int > surfaceID;

    ClassDef(InstalledStation,1);

};



//class Detector : public TObject {
      
class IdealStation{
      
    public:
      
        int nSurfaces;
        int nStrings;
        vector < int > surfaceChannels;
        vector < vector < int > > VHChannel;
        int nChannels;
        int nChannelsVH;
        vector < vector < int > > VHID;
        vector < int > surfaceID;
        vector < int > IDSurface;
        vector < int > IDAntenna;
        vector < int > IDString;
        
	ClassDef(IdealStation, 1);
	
};
  
class Detector {
    private:
        static const int freq_step_max = 60;
        static const int ang_step_max = 2664;
        void ReadVgain(string filename);
        void ReadVgainSettings(string filename, Settings *settings1);
	void ReadVgainTopSettings(string filename, Settings *settings1);
        void ReadHgain(string filename);
	void ReadHgainSettings(string filename, Settings *settings1);
        double Vgain[freq_step_max][ang_step_max];
        double Vphase[freq_step_max][ang_step_max];
        double VgainTop[freq_step_max][ang_step_max];
        double VphaseTop[freq_step_max][ang_step_max];
        double Hgain[freq_step_max][ang_step_max];
        double Hphase[freq_step_max][ang_step_max];
        double Freq[freq_step_max];


	

        void ReadFilter(string filename, Settings *settings1);
        double FilterGain[freq_step_max];   // Filter gain (dB) for Detector freq bin array
        vector <double> FilterGain_databin;   // Filter gain (dB) for DATA_BIN_SIZE bin array
        vector <double> FilterGain_NFOUR;   // Filter gain (dB) for NFOUR bin array

        void ReadPreamp(string filename, Settings *settings1);
        double PreampGain[freq_step_max];   // Filter gain (dB) for Detector freq bin array
        vector <double> PreampGain_databin;   // Filter gain (dB) for DATA_BIN_SIZE bin array
        vector <double> PreampGain_NFOUR;   // Filter gain (dB) for NFOUR bin array


        void ReadFOAM(string filename, Settings *settings1);
        double FOAMGain[freq_step_max];   // Filter gain (dB) for Detector freq bin array
        vector <double> FOAMGain_databin;   // Filter gain (dB) for DATA_BIN_SIZE bin array
        vector <double> FOAMGain_NFOUR;   // Filter gain (dB) for NFOUR bin array



        void ReadCalPulserWF(string filename, Settings *settings1 );  // will store calpulser waveform array



        void ReadElectChain(string filename, Settings *settings1);
        double ElectGain[freq_step_max];   // Elect chain gain (unitless) for Detector freq bin array
        double ElectPhase[freq_step_max];   // Elect chain phase (rad) for Detector freq bin array
        vector <double> ElectGain_databin;   // Elect gain (unitless) for DATA_BIN_SIZE bin array
        vector <double> ElectPhase_databin;   // Elect phase (rad) for DATA_BIN_SIZE bin array
        vector <double> ElectGain_NFOUR;   // Elect gain (unitless) for NFOUR bin array
        vector <double> ElectPhase_NFOUR;   // Elect phase (rad) for NFOUR bin array



        void ReadGainOffset_TestBed(string filename, Settings *settings1);
        vector <double> GainOffset_TB_ch;   // constant gain offset for the TestBed chs 

        void ReadThresOffset_TestBed(string filename, Settings *settings1);
        vector <double> ThresOffset_TB_ch;   // constant gain offset for the TestBed chs 

        void ReadThres_TestBed(string filename, Settings *settings1);
        vector <double> Thres_TB_ch;   // Threshold values for the TestBed chs 

        void ReadTemp_TestBed(string filename, Settings *settings1);
        //vector <double> Temp_TB_ch;   // constant gain offset for the TestBed chs 

        void ReadRFCM_TestBed(string filename, Settings *settings1);
        double RFCM_TB_ch[16][freq_step_max];   // Filter gain (dB) for Detector freq bin array
        vector < vector <double> > RFCM_TB_databin_ch;   // RFCM gain measured value for the TestBed (for each ch)


        void ReadRayleighFit_TestBed(string filename, Settings *settings1, int ch_no); // will read Rayleigh fit result from the file
        void ReadRayleighFit_TestBed(string filename, Settings *settings1); // will read Rayleigh fit result from the file
        double Rayleigh_TB_ch[16][freq_step_max];   // Filter gain (dB) for Detector freq bin array
        vector < vector <double> > Rayleigh_TB_databin_ch;   // RFCM gain measured value for the TestBed (for each ch)

        void ReadRayleighFit(string filename, Settings *settings1); // will read Rayleigh fit result from the file
        vector < vector <double> > Rayleigh_databin_ch;   // RFCM gain measured value for the TestBed (for each ch)

	void ReadNoiseFigure(string filename, Settings *settings1); 
	double NoiseFig_ch[16][freq_step_max];
	vector < vector < double > > NoiseFig_databin_ch;


	vector <double> transV_databin;
	vector <double> transVTop_databin;
	vector <double> transH_databin;


        void FlattoEarth_ARA(IceModel *icesurface);
        void FlattoEarth_ARA_sharesurface(IceModel *icesurface);  // each station share the lowest surface


        void AddAdditional_Depth(Settings *settings1); // each station share the lowest surface



        int freq_step;
        int ang_step;
        double freq_width;
        double freq_init;
        int Detector_mode;

    public:
        Parameters params;
        Detector ();    //default constructor
        Detector (Settings *settings1, IceModel *icesurface, string setupfile);
        //Detector (int mode, IceModel *icesurface);
        vector <ARA_station> stations;
        vector <Antenna_string> strings;

	int NoiseFig_numCh;

        vector <double> freq_forfft;

        double GetGain(double freq, double theta, double phi, int ant_m, int ant_o);    //read antenna gain at certain angle, certain type, and certain orientation
        double GetGain(double freq, double theta, double phi, int ant_m);   //read antenna gain at certain angle, certain type. (orientation : default)

        double GetGain_1D(double freq, double theta, double phi, int ant_m);   //read antenna gain at certain angle, certain type. (orientation : default) and use 1-D interpolation to get gain

        double GetGain_1D_OutZero(double freq, double theta, double phi, int ant_m);   //read antenna gain at certain angle, certain type. (orientation : default) and use 1-D interpolation to get gain, if freq bigger than freq range, return 0 gain

        double GetGain_1D_OutZero(double freq, double theta, double phi, int ant_m, int ant_number);   //read antenna gain at certain angle, certain type. (orientation : default) and use 1-D interpolation to get gain, if freq bigger than freq range, return 0 gain

	


        double GetAntPhase(double freq, double theta, double phi, int ant_m); // return antenna phase with 2-D interpolation

        double GetAntPhase_1D(double freq, double theta, double phi, int ant_m); // return antenna phase with 1-D interpolation


        double GetFilterGain(int bin) { return FilterGain[bin]; }   // same bin with Vgain, Hgain
        double GetFilterGain_databin(int bin) { return FilterGain_databin[bin]; }   // bin for FFT
        double GetFilterGain_NFOUR(int bin) { return FilterGain_NFOUR[bin]; }   // bin for FFT

        double GetFilterGain_1D_OutZero(double freq); // interpolated output, with outside band returns zero

        double GetPreampGain(int bin) { return PreampGain[bin]; }   // same bin with Vgain, Hgain
        double GetPreampGain_databin(int bin) { return PreampGain_databin[bin]; }   // bin for FFT
        double GetPreampGain_NFOUR(int bin) { return PreampGain_NFOUR[bin]; }   // bin for FFT
        double GetPreampGain_1D_OutZero(double freq);

        double GetFOAMGain(int bin) { return FOAMGain[bin]; }   // same bin with Vgain, Hgain
        double GetFOAMGain_databin(int bin) { return FOAMGain_databin[bin]; }   // bin for FFT
        double GetFOAMGain_NFOUR(int bin) { return FOAMGain_NFOUR[bin]; }   // bin for FFT
        double GetFOAMGain_1D_OutZero(double freq);

        double GetElectGain(int bin) { return ElectGain[bin]; }   // same bin with Vgain, Hgain
        double GetElectGain_databin(int bin) { return ElectGain_databin[bin]; }   // bin for FFT
        double GetElectGain_NFOUR(int bin) { return ElectGain_NFOUR[bin]; }   // bin for FFT
        double GetElectGain_1D_OutZero(double freq);
        double GetElectPhase_1D(double freq);



        double GetRFCMGain(int ch, int bin) { return RFCM_TB_ch[ch][bin]; }   // same bin with Vgain, Hgain
        double GetRFCMGain_databin(int ch, int bin) { return RFCM_TB_databin_ch[ch][bin]; }   // bin for FFT


        double GetRayleighFit(int ch, int bin) { return Rayleigh_TB_ch[ch][bin]; }   // same bin with Vgain, Hgain
        double GetRayleighFit_databin(int ch, int bin) { return Rayleigh_TB_databin_ch[ch][bin]; }   // bin for FFT


        double GetNoiseFig_databin(int ch, int bin) { return NoiseFig_databin_ch[ch%16][bin]; }   // bin for FFT
        double GetTransm_databin(int ch, int bin) {	if(ch%16<4){return transVTop_databin[bin];}
							else if(ch%16<8){return transV_databin[bin];} 
							else{return transH_databin[bin];} }   // bin for FFT
        void ReadNoiseFig_New(Settings *settings1); // get noise Figure array with new DATA_BIN_SIZE

        
        double GetGainOffset( int StationID, int ch, Settings *settings1 );  // returns voltage factor for specific channel gain off set

        double GetThresOffset( int StationID, int ch, Settings *settings1 );  // returns voltage factor for specific channel gain off set

        double GetThres( int StationID, int ch, Settings *settings1 );  // returns voltage factor threshold for specific channel 

        double GetTemp( int StationID, int ch, Settings *settings1 );  // returns voltage factor for specific channel gain off set
        vector <double> Temp_TB_ch;   // constant gain offset for the TestBed chs 


        double Getfreq_init() {return freq_init;}

        int Get_mode() {return Detector_mode;}

        int GetFreqBin() {return freq_step;}
        double GetFreq(int bin) {return Freq[bin]*1.e6;} //from MHz to Hz

        vector <double> diode_real; // NFOUR/2 array of t domain tunnel diode response. same with icemc -> anita -> diode_real  but only full bandwidth array 4
        vector <double> fdiode_real_databin;    // NFOUR array of f domain tunnel diode response (FFT of diode_real). also same with icemc -> anita -> fdiode_real  but only full bandwidth array 4
        vector <double> fdiode_real;    // NFOUR/2 array of f domain tunnel diode response (FFT of diode_real). also same with icemc -> anita -> fdiode_real  but only full bandwidth array 4
        vector <double> fdiode_real_double;    // NFOUR array of f domain tunnel diode response (FFT of diode_real). also same with icemc -> anita -> fdiode_real  but only full bandwidth array 4
        
        double TIMESTEP;    // will copy TIMESTEP from Settings
        int NFOUR;          // also will copy NFOUR from Settings

        //TF1 fdiode;   // removed. so not exactly same as icemc, but I believe it doesn't matter
        double maxt_diode;
        int maxt_diode_bin; // above maxt_diode in bin
        int idelaybeforepeak;
        int iwindow;
        int ibinshift;

        int RayleighFit_ch; // number of chs from RayleighFit

        int max_number_of_antennas_station; // maximum number of antennas in a station

        void getDiodeModel(Settings *settings1);   // similar with icemc -> anita -> getDiodeModel().  set diode_real and fdiode_real values.


        // this is a test version for getting new noise waveforms for each event
        // for a best performance, we can just set a new reasonable DATA_BIN_SIZE and make new values for those
        void get_NewDiodeModel(Settings *settings1);

        void ReadFilter_New(Settings *settings1);    // get filter vector array with new DATA_BIN_SIZE 

        void ReadPreamp_New(Settings *settings1);    // get filter vector array with new DATA_BIN_SIZE 

        void ReadFOAM_New(Settings *settings1);    // get filter vector array with new DATA_BIN_SIZE 

        void ReadElectChain_New(Settings *settings1);    // get elect chain vector array with new DATA_BIN_SIZE 

        void ReadRFCM_New(Settings *settings1);    // get filter vector array with new DATA_BIN_SIZE 

        void ReadRayleigh_New(Settings *settings1); // get Rayleigh fit array with new DATA_BIN_SIZE


    
//    vector < vector < vector < int > > > ChannelfromStringAntenna;
//    void SetChannelStringAntennaMap();
    int GetChannelfromStringAntenna (int stationNum, int stringnum, int antennanum );
    void GetSSAfromChannel ( int stationNum, int channelNum, int * antennaNum, int * stringNum );
    
#ifdef ARA_UTIL_EXISTS
    void UseAntennaInfo (int stationNum, Settings *settings1);
    void ImportStationInfo (Settings *settings1, int StationIndex, int StationID);
#endif

// more general used function
    void GetSSAfromChannel ( int stationID, int channelNum, int * antennaNum, int * stringNum, Settings *settings1);
    int GetChannelfromStringAntenna ( int stationID, int stringnum, int antennanum, Settings *settings1);
   
    
    /*
    struct InstalledStation {
        int nSurfaces;
        int nStrings;
        vector < int > surfaceChannels;
        vector < vector < int > > VHChannel;
        int nChannels;
        int nChannelsVH;
        vector < vector < int > > VHID;
        vector < int > surfaceID;

    };
    */
    
    vector < InstalledStation > InstalledStations;
    
//     class IdealStation{
//       
//     public:
//       
//         int nSurfaces;
//         int nStrings;
//         vector < int > surfaceChannels;
//         vector < vector < int > > VHChannel;
//         int nChannels;
//         int nChannelsVH;
//         vector < vector < int > > VHID;
//         vector < int > surfaceID;
//         vector < int > IDSurface;
//         vector < int > IDAntenna;
//         vector < int > IDString;
//         
//     };
    
    vector < IdealStation > IdealStations;
    

    void SetupInstalledStations();
    void PrepareVectorsInstalled();
    void PrepareVectorsInstalled(int importedStation);

    void SetupIdealStations();

    int getAntennafromArbAntID( int stationID, int ant_ID);
    int getStringfromArbAntID( int stationID, int ant_ID);
    
    
    vector <double> CalPulserWF_ns;
    vector <double> CalPulserWF_V;


        ~Detector();    //destructor

        ClassDef(Detector,1);
        
    
    
    
};



#endif //DETECTOR_H
