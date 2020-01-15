#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "Settings.h"
#include "Detector.h"

bool AraUtilExists = false;

#ifdef ARA_UTIL_EXISTS
#include "AraRootVersion.h"
#endif

ClassImp(Settings);

using namespace std;

Settings::Settings() {
    Initialize();
//    ReadFile();

}

Settings::~Settings() {
    //default destructor
}


void Settings::Initialize() {


// below : values from icemc Settings class
  NDISCONES_PASS=3;

  DEBUG=false;                   // debugging option
outputdir="outputs"; // directory where outputs go
 FREQ_LOW_SEAVEYS=200.E6;
 FREQ_HIGH_SEAVEYS=1200.E6;
 BW_SEAVEYS=FREQ_HIGH_SEAVEYS-FREQ_LOW_SEAVEYS;
 SIGMAPARAM=1;  // Connolly et al. 2011 default cross section parametrization
 SIGMA_FACTOR=1.;   // default sigma factor : 1
 YPARAM=1;  // 1: Connolly et al. 2011 default y parametrization, 2: Set ELAST_Y yourself
 ELAST_Y = 0.0;
 UNBIASED_SELECTION=1.; // (0) pick neutrino interaction in the ice and neutrino from any direction or (1) choose neutrino interaction point in the horizon on the balloon in the ice and neutrino direction on the cerenkov cone

 SIGMA_SELECT=0; // when in SIGMAPARAM=1 case, 0 : (default) use mean value, 1 : use upper bound, 2 : use lower bound


// end of values from icemc

 ARASIM_VERSION_MAJOR = ARASIM_MAJOR;
 ARASIM_VERSION_MINOR = ARASIM_MINOR;
 ARASIM_VERSION_SUBMINOR = ARASIM_SUBMINOR;
 ARASIM_VERSION = (double)ARASIM_VERSION_MAJOR + (double)ARASIM_VERSION_MINOR * 0.001 + (double)ARASIM_VERSION_SUBMINOR * 0.000001;
 
 ARAROOT_VERSION = 0.;

 ARAUTIL_EXISTS = false;
#ifdef ARA_UTIL_EXISTS
 ARAUTIL_EXISTS = true;
 ARAROOT_VERSION = (double)ARA_ROOT_MAJOR + (double)ARA_ROOT_MINOR * 0.01;
#endif
 

  NNU=100;

  // NEED TO FIGURE OUT A GOOD WAY TO READ THIS IN AND STORE THEM.
  // INPUT FILE AGAIN?  SOMETHING ELSE?
  //These were moved here from IceModel under the new compilation scheme
  ICE_MODEL=0; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP.
  NOFZ=1; // 1=depth dependent index of refraction,0=off
  CONSTANTCRUST=0; // set crust density and thickness to constant values.
  CONSTANTICETHICKNESS=0; // set ice thickness to constant value
  FIXEDELEVATION=0; // fix the elevation to the thickness of ice.
  MOOREBAY=0; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data
  USE_ARA_ICEATTENU=1; // use ARA measured ice attenuation value
  
  EXPONENT=19.; // 10^19 eV neutrinos only

  DETECTOR=1;   //ARA layout with small number of stations

  INTERACTION_MODE=1;   //PickNear mode (0: Aeff mode using sphere surface around station, 1: Veff mode using cylinder volume around station)

  POSNU_RADIUS=3000;    //radius for PickNear method

  WHICHPARAMETERIZATION=0;  //

  SIMULATION_MODE=1;    // default freq domain simulation

  EVENT_TYPE=0;         // default neutrino only events

  WAVE_TYPE=0;          // default wave type : plane wave (inside the ice)

  LPM=1;                //default : enable LPM effect

  SECONDARIES=1;        //default : enable secondary interactions

  TAUDECAY=1;           //default : let taudecay as secondary interactions

  TIMESTEP=(0.625)*1.E-9;  // default, in sec (old default: 0.5E-9, new default 0.625E-9

  PHASE=90.;            // default : 90 deg phase (it means all imaginary values)

  NFOUR=1024;           // default : 1024, same as in icemc
    
  NOISE=0;              // degault : 0, flat thermal noise, 1 : for TestBed (DETECTOR=3), use Rayleigh distribution fitted for borehole channels

  ATMOSPHERE=1;         // default : 1, include atmosphere

  TRIG_SCAN_MODE=0;	// default 0 (old mode) 1: new mode (faster) 2: scan all Pthresh values 3: scan also all N out of 8 
  
  POWERTHRESHOLD=-6.06; // old default : -6.15, new default: -6.06

  MAXT_DIODE=70.E-9;    // default : 70 ns

  IDELAYBEFOREPEAK_DIODE=(int)(13.E-9 / TIMESTEP);    // default : 13.e-9/TIMESTEP = 33

  IWINDOW_DIODE=(int)(4.E-9 / TIMESTEP);           // default : 4.e-9 / TIMESTEP = 10

  DATA_BIN_SIZE=16384;   // default : 16384

  NOISE_TEMP=325.;      // default : 325 K

  TRIG_ANALYSIS_MODE=0;    // default : 0, signal + noise

  TRIG_TIMEOUT=1.E-6;       // default : 1us

  TRIG_WINDOW=1.7E-7;       // old default : 110 ns, new default: 170 ns

  NOISE_EVENTS=16;        // default : 16 events

  DATA_SAVE_MODE=0;         // default : 0 (full mode)

  N_TRIG=3;                 // default : 3 (3 out of all channels in a station)

  RANDOM_MODE=1;            // default : 1 (seed is unique in time/space)

  SEED=1; // default: 1, only applies if RANDOM_MODE=0, provides base seed value and run number taken from arguments is added to this value in order to submit multiple repeatable runs instead of only one single long repeatable run
    
  BORE_HOLE_ANTENNA_LAYOUT=0;   // default : 0 (VHVH)

  DATA_LIKE_OUTPUT=1; //default : 0 (doesn't write out data-like events)
    
  RAYSOL_RANGE=5000; // default : 5000 m

  PICK_POSNU_DEPTH=0;     //default : 0 pick posnu depth from 0 to ice depth

  MAX_POSNU_DEPTH=200.;     // default : 200m depth max

  NNU_THIS_THETA=0;         // default 0: nnu angle pure random, 1: set a specific theta

  NNU_THETA=0.785;          // default : nnu theta : 45 deg

  NNU_D_THETA=0.0873;       // default : nnu d_theta : 5 deg

  NNU_THIS_PHI=0;//default 0: random phi, 1: a specific phi

  NNU_PHI=0.785;// default : nnu phi : 45 deg

  NNU_D_PHI=0.0873;// default : nnu_d_phi : 5 deg

    
    CALPULSER_ON=0; // default : calpulsers off
    
    TESTBED_ON=0; // default : 0 stations[0] is ARA1 not Testbed
    
    READGEOM=0; // default : 0 : use idealized geometry and do not read in from sqlite database
    
    V_MIMIC_MODE = 0; // default : 0 - write out all chs where global triggered bin is center of the window
                        // 1 - same as above 0 mode but apply TestBed ch delay - average BH ch delay
                        // 2 - same as above 0 mode but apply TestBed ch delay - average BH ch delay + additional delay to match with actual TestBed data waveforms
    
    USE_INSTALLED_TRIGGER_SETTINGS = 0; // default : 0 - use idealized settings for the trigger
    
    NUM_INSTALLED_STATIONS = 4;

    CALPUL_OFFCONE_ANGLE = 35.;

    CALPUL_AMP = 100.;

    TRIG_ONLY_BH_ON = 0;    // default trigger will occur with all chs (1 will do trigger analysis with BH chs only)

    TRIG_THRES_MODE = 0;    // default trigger threshold (0) will use 1 as offset (so no offset), (1) will use data/threshold_offset.csv as threshold off set factor

    NOISE_CHANNEL_MODE = 0;    //default noise temp setting (just same temp for all chs), 1 : all chs have different systemp, 2 : only first 8 chs have different systemp

    USE_TESTBED_RFCM_ON = 0;    // use RFCM measurement for testbed or not

    RFCM_OFFSET = 80.;  // if above USE_TESTBED_RFCM_ON = 1, we need RFCM attenuator factor cancel

    CONST_MEANDIODE = -6.5e-15; // just from one run

    CONST_RMSDIODE = 1.346e-13; // also from one run

    USE_MANUAL_GAINOFFSET = 0; //if use gain offset file to read values or just use constant gain offset from setup file (default 0 : use file)
            
    MANUAL_GAINOFFSET_VALUE = 1.; // gain offset value

    NOISE_WAVEFORM_GENERATE_MODE = 0; // mode 0 (default) will generate noise waveforms newly for each events. other values will use first generated noise waveforms for later events (huge mem usage)

    USE_CH_GAINOFFSET = 0; // if use gain offset for different channels. (default 0 : not using gain offset). mode 1 is only availbale for installed TestBed so far.

    // removed GETCHORD_MODE. This parameter is merged into INTERACTION_MODE
    //GETCHORD_MODE = 0; // which Getchord function to use. default 0 : old Getchord function (not correct tau weight, weight don't have ice inside interaction probability in it). 1 : new Getchord from icemc. This has new tau weight calculation and ice interaction probability applied to weight factor.

    taumodes = 0; // no tau created in the rock

    BH_ANT_SEP_DIST_ON = 1; // 0 : use constant borehole antenna distance. default 1 : use separate antenna distance. use z_btw01, z_btw12, ... in ARA_N_info.txt or ARA37_info.txt

    TRIG_MODE = 1; // default 1 : if any antennas got passed N_TRIG or more, global trig. 1 : either Vpol or Hpol antennas got passed N_TRIG_V or N_TRIG_H respectively, global trig.

    N_TRIG_V=3;                 // default : 3 (3 out of Vpolchannels in a station)

    N_TRIG_H=3;                 // default : 3 (3 out of Hpol channels in a station)

    FILL_TREE_MODE = 0; // default 0 : fill tree for all events, 1 : fill tree only usable posnu events, 2 : fill tree only trigger passed events

    ONLY_PASSED_EVENTS = 0;
    NNU_PASSED = 0;



    SHOWER_MODE = 2; // EM (0) or HAD (1) shower in t-domain signal. or either one which is bigger (3) or both EM and HAD (2)  default : 2, both EM and HAD showers

    SHOWER_STEP = 0.001; // step size in generating shower profile. default 0.001 m

    SHOWER_PARAM_MODEL = 0; // choose shower profile parameters (by Jaime fit = 0, or Carl's fit = 1). default = 0

    OFFCONE_LIMIT = 10.; // offcone angle (deg) limit to calculate time domain signal. Increasing this value will result in drametically increase computation time

    ALL_ANT_V_ON = 1; // use Vpol antenna gain for both Vpol and Hpol = 1, use Hpol gain for Hpol model = 0

    PHASE_SKIP_MODE = 0; // skip applying phase in t-domain mode (SIMULATION_MODE = 1). default 0 : don't skip (apply all phase), 1 : only upto Askaryan radiation, 2 : only upto antenna



    DEBUG_MODE_ON = 0; // 0 : off (do as usual), 1 : on (skip most of intensive computational process which don't have random generations)

    DEBUG_SKIP_EVT = 0; // when DEBUG_MODE_ON = 1, skip upto this number and then do as DEBUG_MODE_ON = 0


    V_SATURATION = 1.; // saturated voltage +-V_SATURATION

    ADDITIONAL_DEPTH_ON = 0; // whether add more depth to each antenas

    ADDITIONAL_DEPTH = 100.; // default additional depth value



    TRIG_ONLY_LOW_CH_ON = 0;    // default trigger will occur with all chs (1 will do trigger analysis with lower 8 chs; bottom 4 Vpol & bottom 4 Hpols)


    ACCUM_TRIG_SEARCH_BINS_STATION0 = 0.; // not actually setting value but gives us how much trigger searched bins there were in the run for station0

    NU_NUBAR_SELECT_MODE = 3; // default : 3 = random nu_nubar based on arXiv:1108.3163, section 3, 0 = just nu, 1 = just nubar 


    SELECT_FLAVOR = 0; // default : 0 = randomly 1:1:1 ratio, 1 : el. 2 : mu, 3 : tau
    SELECT_CURRENT = 2; // default: 2:random, 0:nc, 1:cc

    OUTPUT_TDR_GRAPH = 0;// saves a few example graphs of the tunnel diode response for a triggered event


    AVZ_NORM_FACTOR_MODE = 1; // default : 1 : don't apply sqrt(2) (actually applied but cancel that) as realft assume Hn as double-sided spectrum (invFFT normalization factor 2/N) and also remove dF binning factor in MakeArraysforFFT function, 0 : use normalization factors like in old version

    number_of_stations = 1;

    RAY_TRACE_ICE_MODEL_PARAMS=0; // Default: South Pole values fitted from RICE data

    WAVEFORM_LENGTH = 64/2*20; // Default: 64 digitization samples per block / 2 samples per waveform value * 20 blocks (value used for 2013-2016)
    
    WAVEFORM_CENTER = 0; // Default: 0, no offset in waveform centering

    POSNU_R = 1000.;
    POSNU_THETA=-3.1415926535/4.;
    POSNU_PHI=0.;

    ARBITRARY_EVENT_ATTENUATION = 1.0;
    PICK_ABOVE_HEIGHT = 3000;

    EVENT_GENERATION_MODE = 0;//default: 0: not event mode, 1: event mode
    //    EVENT_NUM = 10;//read in event number in EVENT_GENERATION_MODE=1, no more than 100 events
    ANTENNA_MODE=0; //default: 0 - old antenna model information
    APPLY_NOISE_FIGURE=0; // default: 0 - don't use new noise figure information

    CUSTOM_ELECTRONICS=0; //default: 0 -- don't use custom electronics, load regular "ARA_Electronics_TotalGain_TwoFilter.tst"


    /*
//arrays for saving read in event features in EVENT_GENERATION_MODE=1
    EVID[100] = {0};
    NUFLAVORINT[100] = {0};
    NUBAR[100] = {0};
    PNU[100] = {0};
    CURRENTINT[100] = {0};
    IND_POSNU_R[100] = {0};
    IND_POSNU_THETA[100] = {0};
    IND_POSNU_PHI[100] = {0};
    IND_NNU_THETA[100] = {0};
    IND_NNU_PHI[100] = {0};
    */
}

void Settings::ReadFile(string setupfile) {

  ifstream setFile (setupfile.c_str());
  
  string line, label;

  if ( setFile.is_open() ) {
      while (setFile.good() ) {
          getline (setFile,line);

          if (line[0] != "/"[0]) {
              label = line.substr(0, line.find_first_of("="));
              
              if (label == "NNU") {
                  NNU = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "ICE_MODEL") {
                  ICE_MODEL = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NOFZ") {
                  NOFZ = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "CONSTANTCRUST") {
                  CONSTANTCRUST = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "CONSTANTICETHICKNESS") {
                  CONSTANTICETHICKNESS = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "FIXEDELEVATION") {
                  FIXEDELEVATION = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "MOOREBAY") {
                  MOOREBAY = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "EXPONENT") {
                  EXPONENT = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "DETECTOR") {
                  DETECTOR = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "DETECTOR_STATION") {
                  DETECTOR_STATION = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "INTERACTION_MODE") {
                  INTERACTION_MODE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "POSNU_RADIUS") {
                  POSNU_RADIUS = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "WHICHPARAMETERIZATION") {
                  WHICHPARAMETERIZATION = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "SIMULATION_MODE") {
                  SIMULATION_MODE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "EVENT_TYPE") {
                  EVENT_TYPE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "WAVE_TYPE") {
                  WAVE_TYPE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "LPM") {
                  LPM = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "SECONDARIES") {
                  SECONDARIES = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TAUDECAY") {
                  TAUDECAY = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TIMESTEP") {
                  TIMESTEP = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "PHASE") {
                  PHASE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NFOUR") {
                  NFOUR = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NOISE") {
                  NOISE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "ATMOSPHERE") {
                  ATMOSPHERE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if(label == "TRIG_SCAN_MODE"){
                  TRIG_SCAN_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
	      }
              else if (label == "POWERTHRESHOLD") {
                  POWERTHRESHOLD = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "DATA_BIN_SIZE") {
                  DATA_BIN_SIZE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NOISE_TEMP") {
                  NOISE_TEMP = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TRIG_ANALYSIS_MODE") {
                  TRIG_ANALYSIS_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TRIG_TIMEOUT") {
                  TRIG_TIMEOUT = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TRIG_WINDOW") {
                  TRIG_WINDOW = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NOISE_EVENTS") {
                  NOISE_EVENTS = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "DATA_SAVE_MODE") {
                  DATA_SAVE_MODE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "N_TRIG") {
                  N_TRIG = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "RANDOM_MODE") {
                  RANDOM_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "SEED") {
                  SEED = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "BORE_HOLE_ANTENNA_LAYOUT") {
                  BORE_HOLE_ANTENNA_LAYOUT = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "RAYSOL_RANGE") {
                  RAYSOL_RANGE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "CALPULSER_ON") {
                  CALPULSER_ON = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TESTBED_ON") {
                  TESTBED_ON = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "READGEOM") {
                  cout << "Read in READGEOM" << endl;
                  READGEOM = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "PICK_POSNU_DEPTH") {
                  PICK_POSNU_DEPTH = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "MAX_POSNU_DEPTH") {
                  MAX_POSNU_DEPTH = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NNU_THIS_THETA") {
                  NNU_THIS_THETA = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NNU_THETA") {
                  NNU_THETA = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NNU_D_THETA") {
                  NNU_D_THETA = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NNU_THIS_PHI") {
                  NNU_THIS_PHI = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NNU_PHI") {
                  NNU_PHI = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NNU_D_PHI") {
                  NNU_D_PHI = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }

              else if (label == "DATA_LIKE_OUTPUT") {
                  DATA_LIKE_OUTPUT = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "V_MIMIC_MODE") {
                  V_MIMIC_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "USE_INSTALLED_TRIGGER_SETTINGS") {
                  USE_INSTALLED_TRIGGER_SETTINGS = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NUM_INSTALLED_STATIONS") {
                  NUM_INSTALLED_STATIONS = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }              
              else if (label == "CALPUL_OFFCONE_ANGLE") {
                  CALPUL_OFFCONE_ANGLE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "CALPUL_AMP") {
                  CALPUL_AMP = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TRIG_ONLY_BH_ON") {
                  TRIG_ONLY_BH_ON = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }              
              else if (label == "TRIG_THRES_MODE") {
                  TRIG_THRES_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }              
              else if (label == "NOISE_CHANNEL_MODE") {
                  NOISE_CHANNEL_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "CONST_MEANDIODE") {
                  CONST_MEANDIODE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "CONST_RMSDIODE") {
                  CONST_RMSDIODE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "USE_TESTBED_RFCM_ON") {
                  USE_TESTBED_RFCM_ON = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "RFCM_OFFSET") {
                  RFCM_OFFSET = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "USE_MANUAL_GAINOFFSET") {
                  USE_MANUAL_GAINOFFSET = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "MANUAL_GAINOFFSET_VALUE") {
                  MANUAL_GAINOFFSET_VALUE = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NOISE_WAVEFORM_GENERATE_MODE") {
                  NOISE_WAVEFORM_GENERATE_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "USE_CH_GAINOFFSET") {
                  USE_CH_GAINOFFSET = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              //else if (label == "GETCHORD_MODE") {
                  //GETCHORD_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              //}
              else if (label == "taumodes") {
                  taumodes = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "BH_ANT_SEP_DIST_ON") {
                  BH_ANT_SEP_DIST_ON = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TRIG_MODE") {
                  TRIG_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "N_TRIG_V") {
                  N_TRIG_V = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "N_TRIG_H") {
                  N_TRIG_H = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "FILL_TREE_MODE") {
                  FILL_TREE_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "ONLY_PASSED_EVENTS") {
                  ONLY_PASSED_EVENTS = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NNU_PASSED") {
                  NNU_PASSED = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "PICKNEARUNBIASED_R") {
                  PICKNEARUNBIASED_R = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "SHOWER_MODE") {
                  SHOWER_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "SHOWER_STEP") {
                  SHOWER_STEP = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "SHOWER_PARAM_MODEL") {
                  SHOWER_PARAM_MODEL = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "ALL_ANT_V_ON") {
                  ALL_ANT_V_ON = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "PHASE_SKIP_MODE") {
                  PHASE_SKIP_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "DEBUG_MODE_ON") {
                  DEBUG_MODE_ON = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "DEBUG_SKIP_EVT") {
                  DEBUG_SKIP_EVT = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "V_SATURATION") {
                  V_SATURATION = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "OFFCONE_LIMIT") {
                  OFFCONE_LIMIT = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "ADDITIONAL_DEPTH_ON") {
                  ADDITIONAL_DEPTH_ON = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "ADDITIONAL_DEPTH") {
                  ADDITIONAL_DEPTH = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "TRIG_ONLY_LOW_CH_ON") {
                  TRIG_ONLY_LOW_CH_ON = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }              
              else if (label == "USE_ARA_ICEATTENU") {
                  USE_ARA_ICEATTENU = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "SIGMA_SELECT") {
                  SIGMA_SELECT = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }              
              else if (label == "SIGMAPARAM") {
                  SIGMAPARAM = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }              
              else if (label == "SIGMA_FACTOR") {
                  SIGMA_FACTOR = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "YPARAM") {
                  YPARAM = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }              
              else if (label == "ELAST_Y") {
                 ELAST_Y = atof( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "NU_NUBAR_SELECT_MODE") {
                  NU_NUBAR_SELECT_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }              
              else if (label == "SELECT_FLAVOR") {
                  SELECT_FLAVOR = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "SELECT_CURRENT") {
                  SELECT_CURRENT = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }
              else if (label == "OUTPUT_TDR_GRAPH") {
                  OUTPUT_TDR_GRAPH = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }       
              else if (label == "AVZ_NORM_FACTOR_MODE") {
                  AVZ_NORM_FACTOR_MODE = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
              }              
	      else if (label == "number_of_stations") {
		number_of_stations = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
	      }
	      else if (label == "RAY_TRACE_ICE_MODEL_PARAMS") {
		RAY_TRACE_ICE_MODEL_PARAMS = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
	      }
	      else if (label == "WAVEFORM_LENGTH") {
		WAVEFORM_LENGTH = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
	      }
	      else if (label == "WAVEFORM_CENTER") {
		WAVEFORM_CENTER = atoi( line.substr(line.find_first_of("=") + 1).c_str() );
	      }
	      else if (label == "POSNU_R") {
		POSNU_R = atof( line.substr(line.find_first_of("=") + 1).c_str() );
	      }
	      else if (label == "POSNU_THETA") {
		POSNU_THETA = atof( line.substr(line.find_first_of("=") + 1).c_str() );
	      }
	      else if (label == "POSNU_PHI") {
		POSNU_PHI = atof( line.substr(line.find_first_of("=") + 1).c_str() );
	      }
	      else if (label == "ARBITRARY_EVENT_ATTENUATION") {
		ARBITRARY_EVENT_ATTENUATION = atof( line.substr(line.find_first_of("=") + 1).c_str() );
	      }
	      else if (label == "PICK_ABOVE_HEIGHT") {
		PICK_ABOVE_HEIGHT = atof( line.substr(line.find_first_of("=") + 1).c_str() );
	      }
              else if (label == "EVENT_GENERATION_MODE"){
                  EVENT_GENERATION_MODE = atoi(line.substr(line.find_first_of("=") + 1).c_str());
              }
	      //              else if (label == "EVENT_NUM"){
	      //                  EVENT_NUM = atoi(line.substr(line.find_first_of("=") + 1).c_str());
	      //              }
              else if (label == "ANTENNA_MODE"){
                  ANTENNA_MODE = atoi(line.substr(line.find_first_of("=") + 1).c_str());
              }
              else if (label == "APPLY_NOISE_FIGURE"){
                  APPLY_NOISE_FIGURE = atoi(line.substr(line.find_first_of("=") + 1).c_str());
              }
              else if (label == "CUSTOM_ELECTRONICS"){
              	   CUSTOM_ELECTRONICS = atoi(line.substr(line.find_first_of("=") + 1).c_str());
              }




          }
      }
      setFile.close();
  }
  else cout<<"Unable to open "<<setupfile<<" file!"<<endl;
  return;
}

void Settings::ReadEvtFile(string evtfile){
    ifstream evtFile(evtfile.c_str());

    std::string line;
    int l = 0;

    if ( evtFile.is_open() ) {
        while (evtFile.good() ) {
            getline(evtFile, line);
            if (line[0] != "/"[0]) {
                std::stringstream iss(line);
                int a, b, c, e;
                double d, f, g, h, i, j, k;
                if (!(iss >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k))
                    break;
/*                EVID[i] = atoi(a.c_str());
                NUFLAVORINT[i] = atoi(b.c_str());
                NUBAR[i] = atoi(c.c_str());
                PNU[i] = atof(d.c_str());
                CURRENTINT[i] = atoi(e.c_str());
                X[i] = atof(f.c_str());
                Y[i] = atof(g.c_str());
                Z[i] = atof(h.c_str());
                THETA[i] = atof(i.c_str());
                PHI[i] = atof(j.c_str());
*/
/*
                EVID[l] = a;
                NUFLAVORINT[l] = b;
                NUBAR[l] = c;
                PNU[l] = d;
                CURRENTINT[l] = e;
                IND_POSNU_R[l] = f;
                IND_POSNU_THETA[l] = g;
                IND_POSNU_PHI[l] = h;
                IND_NNU_THETA[l] = i;
                IND_NNU_PHI[l] = j;
*/
                EVID.push_back(a);
                NUFLAVORINT.push_back(b);
                NUBAR.push_back(c);
                PNU.push_back(d);
                CURRENTINT.push_back(e);
                IND_POSNU_R.push_back(f);
                IND_POSNU_THETA.push_back(g);
                IND_POSNU_PHI.push_back(h);
                IND_NNU_THETA.push_back(i);
                IND_NNU_PHI.push_back(j);
                ELAST.push_back(k);


                l++;
            }
        }
        evtFile.close();
	if (NNU == 0 || NNU > EVID.size()){
	  //	  EVENT_NUM = EVID.size();
	  NNU = EVID.size();
	}
	
    }
    else
        cout << "Unable to open " << evtfile << " file!" << endl;
    return;
}

int Settings::CheckCompatibilitiesDetector(Detector *detector) {

    int num_err = 0;

    // if there's something not going to work, count thoes settings

    if (DETECTOR==1 && READGEOM==1 && detector->params.number_of_stations>1) { // currently only ARA1a one station is possible
        cerr<<"DETECTOR=1, READGEOM=1 is currently only availble with number_of_stations=1 in ARA_N_info.txt file!"<<endl;
        num_err++;
    }

    if (DETECTOR==2 && READGEOM==1) { // currently READGEOM (using actual installed stations info) is not available in DETECTOR=2
        cerr<<"DETECTOR=2 and READGEOM=1 is currently not availble! Only ideal stations are available in DETECTOR=2!"<<endl;
        num_err++;
    }

    // check reasonable number of noise waveforms
    if (NOISE_WAVEFORM_GENERATE_MODE == 0) { // if generating new noise waveforms for every events
        if (NOISE_CHANNEL_MODE == 0) {// share all noise waveforms same with other channels
            //if (NOISE_EVENTS < detector->params.number_of_antennas) { // this is too low number of events!
            if (NOISE_EVENTS < detector->max_number_of_antennas_station) { // this is too low number of events!
                cerr<<"NOISE_EVENTS too less! At least use "<<detector->max_number_of_antennas_station<<"!"<<endl;
                num_err++;
            }
        }
        else if (NOISE_CHANNEL_MODE == 1) {// each chs will have separate noise waveforms
            if (NOISE_EVENTS > 1) { // this case 1 waveform is enough for each channels
                cerr<<"NOISE_EVENTS too many! With NOISE_WAVEFORM_GENERATE_MODE==0 and NOISE_CHANNEL_MODE==1, just use NOISE_EVENTS=1"<<endl;
                num_err++;
            }
        }
    }
    if (NOISE_WAVEFORM_GENERATE_MODE == 1) { // if generating noise waveforms in the begining and keep use them
        if (NOISE_CHANNEL_MODE == 0) {// share all noise waveforms same with other channels
            //if (NOISE_EVENTS < detector->params.number_of_antennas) { // this is too low number of events!
            if (NOISE_EVENTS < detector->max_number_of_antennas_station) { // this is too low number of events!
                cerr<<"NOISE_EVENTS too less! At least use "<<detector->max_number_of_antennas_station<<"!"<<endl;
                num_err++;
            }
        }
    }

    // check if there's enough system temperature values prepared for NOISE_CHANNEL_MODE=1
    if (NOISE_CHANNEL_MODE==1) {// use different system temperature values for different chs
      if (DETECTOR==3 && (detector->params.number_of_antennas > (int)(detector->Temp_TB_ch.size())) ) {
	  cout << detector->params.number_of_antennas << " : " <<(int)(detector->Temp_TB_ch.size()) << endl;
            cerr<<"System temperature values are not enough for all channels! Check number of channels you are using and numbers in data/system_temperature.csv"<<endl;
            num_err++;
        }
    }

    return num_err;
}

int Settings::CheckCompatibilitiesSettings() {

    int num_err = 0;


    /*
    // if INTERACTION_MODE is 0 (sphere area and obtain Aeff), make sure using GETCHORD_MODE=1
    if (INTERACTION_MODE==0) { // picknear_sphere mode
        if (GETCHORD_MODE==0) { // but use old getchord mode (not working!)
            cerr<<"In INTERACTION_MODE=0, you have to use GETCHORD_MODE=1"<<endl; 
            num_err++;
        }
    }
    else if (INTERACTION_MODE==1) { // picknear_cylinder mode
        if (GETCHORD_MODE==1) { // but use new getchord mode (not working!)
            cerr<<"In INTERACTION_MODE=1, you have to use GETCHORD_MODE=0"<<endl; 
            num_err++;
        }
    }
    */

    // if BH_ANT_SEP_DIST_ON=1, we can't use READGEOM=1 (actual installed geom)
    if (BH_ANT_SEP_DIST_ON==1 && READGEOM==1) {
        cerr<<"BH_ANT_SEP_DIST_ON=1 is only available in ideal station geom (READGEOM=0)!"<<endl; 
        num_err++;
    }

    // TRIG_MODE=1 (Vpol, Hpol separated) will not work with testbed station mode
    if (TRIG_MODE==1 && DETECTOR==3) {
        cerr<<"TRIG_MODE=1 is not available in TestBed mode (DETECTOR=3)!"<<endl; 
        num_err++;
    }



    // check modes which will only work for actual installed TestBed case
    //
    if (TRIG_ONLY_BH_ON==1 && DETECTOR!=3) {
        cerr<<"TRIG_ONLY_BH_ON=1 only works with DETECTOR=3!"<<endl;
        num_err++;
    }

    /*
    //if (NOISE_CHANNEL_MODE==1 && DETECTOR!=3) {
    if (NOISE_CHANNEL_MODE==1 && DETECTOR!=3 && TRIG_ONLY_LOW_CH_ON!=1) {
        //cerr<<"NOISE_CHANNEL_MODE=1 only works with DETECTOR=3!"<<endl;
        cerr<<"NOISE_CHANNEL_MODE=1 only works with DETECTOR=3 or TRIG_ONLY_LOW_CH_ON=1"<<endl;
        num_err++;
    }
    */

    //if (NOISE_CHANNEL_MODE==2 && DETECTOR!=3) {
    if (NOISE_CHANNEL_MODE==2 && DETECTOR!=3 && TRIG_ONLY_LOW_CH_ON!=1) {
        //cerr<<"NOISE_CHANNEL_MODE=2 only works with DETECTOR=3!"<<endl;
        cerr<<"NOISE_CHANNEL_MODE=2 only works with DETECTOR=3 or TRIG_ONLY_LOW_CH_ON=1"<<endl;
        num_err++;
    }

    if (DETECTOR==3 && READGEOM==0) {
        cerr<<"DETECTOR=3 will always need READGEOM=1"<<endl;
        num_err++;
    }

    if (USE_TESTBED_RFCM_ON==1 && DETECTOR!=3) {
        cerr<<"USE_TESTBED_RFCM_ON=1 only works with DETECTOR=3!"<<endl;
        num_err++;
    }

    if (USE_CH_GAINOFFSET==1 && DETECTOR!=3) {
        cerr<<"USE_CH_GAINOFFSET=1 only works with DETECTOR=3!"<<endl;
        num_err++;
    }

    if (TRIG_THRES_MODE==1 && DETECTOR!=3 && TRIG_ONLY_LOW_CH_ON!=1) {
        //cerr<<"TRIG_THRES_MODE=1 and 2 only works with DETECTOR=3!"<<endl;
        cerr<<"TRIG_THRES_MODE=1 only works with DETECTOR=3 or TRIG_ONLY_LOW_CH_ON=1"<<endl;
        num_err++;
    }

    if (TRIG_THRES_MODE==2 && DETECTOR!=3 && TRIG_ONLY_LOW_CH_ON!=1) {
        //cerr<<"TRIG_THRES_MODE=1 and 2 only works with DETECTOR=3!"<<endl;
        cerr<<"TRIG_THRES_MODE=2 only works with DETECTOR=3 or TRIG_ONLY_LOW_CH_ON=1"<<endl;
        num_err++;
    }

    if (USE_MANUAL_GAINOFFSET==1 && DETECTOR!=3) {
        cerr<<"USE_MANUAL_GAINOFFSET=1 only works with DETECTOR=3!"<<endl;
        num_err++;
    }

    if (CALPULSER_ON!=0 && DETECTOR!=3) {
        cerr<<"CALPULSER_ON=1 and above only works with DETECTOR=3!"<<endl;
        num_err++;
    }

    if ((V_MIMIC_MODE==1||V_MIMIC_MODE==2) && DETECTOR!=3) {
        cerr<<"V_MIMIC_MODE=1 and 2 only works with DETECTOR=3!"<<endl;
        num_err++;
    }
    
    if (USE_MANUAL_GAINOFFSET==1 && USE_CH_GAINOFFSET==1) {
        cerr<<"Can not use USE_MANUAL_GAINOFFSET=1 and USE_CH_GAINOFFSET=1 same time!"<<endl;
        num_err++;
    }

    //if (NOISE==1 && DETECTOR!=3) {
    if (NOISE==1 && DETECTOR!=3 && TRIG_ONLY_LOW_CH_ON!=1) {
        cerr<<"NOISE=1 only works with DETECTOR=3 or TRIG_ONLY_LOW_CH_ON=1"<<endl;
        num_err++;
    }

    if (NOISE==1 && USE_TESTBED_RFCM_ON==1) {
        cerr<<"NOISE=1 only works with USE_TESTBED_RFCM_ON=0!"<<endl;
        num_err++;
    }

    if (NOISE==1 && NOISE_CHANNEL_MODE==0) {
        cerr<<"NOISE=1 don't work with NOISE_CHANNEL_MODE=0!"<<endl;
        num_err++;
    }


    // This is for only ideal stations
    if (TRIG_ONLY_LOW_CH_ON==1 && DETECTOR==3) {
        cerr<<"TRIG_ONLY_LOW_CH_ON=1 doesn't work with DETECTOR=3!"<<endl;
        num_err++;
    }

    if (DATA_LIKE_OUTPUT != 0 && (DETECTOR==0 || DETECTOR==1 || DETECTOR==2)) {
        cerr<<"DATA_LIKE_OUTPUT=1,2 doesn't work with DETECTOR=0,1,2"<<endl;
        cerr<<"DATA_LIKE_OUTPUT controls data-like output into UsefulAtriStationEvent format; without a real station selected (using DETECTOR==3,4), the mapping to the data-like output will not function correctly"<<endl;
        num_err++;
    }

    if (DATA_LIKE_OUTPUT != 0 && (DETECTOR_STATION>3)) {
        cerr<<"DATA_LIKE_OUTPUT=1,2 doesn't work with DETECTOR_STATION>3"<<endl;
        cerr<<"DATA_LIKE_OUTPUT controls data-like output into UsefulAtriStationEvent format; without a real station selected (using DETECTOR==3,4), the mapping to the data-like output will not function correctly"<<endl;
        num_err++;
    }
    /*
    if(NOISE == 2 && (DETECTOR==0 || DETECTOR==1 || DETECTOR==2)){
      cerr << "CALIBRATION_MODE=1 doesn't work with DETECTOR=0,1,2" << endl;
      cerr << "CALIBRATION_MODE controls the response from installed stations 2,3 and thus has no real relevance for DETECTOR=0,1,2" << endl;
      num_err++;
    }
    */

    // This is for installed stations
    if (DETECTOR == 4 ) {
      if (ARAUTIL_EXISTS == false){
	cerr << "DETECTOR=4 only works with an installation of AraRoot" << endl;
	num_err++;
      } else {
	
	cerr << "DETECTOR is set to 4" << endl; 
	cerr << "Setting READGEOM to 1" << endl;
	READGEOM=1;

	cerr << "Setting number_of_stations to 1" << endl;
	number_of_stations = 1;

	if (DETECTOR_STATION <0 || DETECTOR_STATION >= NUM_INSTALLED_STATIONS){
	  cerr << "DETECTOR_STATION is not set to a valid station number" << endl;
	  num_err++;
	}
		
      }
    }

   //Check that DETECTOR_STATION=0 is only used with DETECTOR=3
   if (DETECTOR_STATION==0 && DETECTOR!=3){
      cerr << " DETECTOR_STATION=0 doesn't work with DETECTOR!=3. If you want to work with TestBed, use DETECTOR=3 & DETECTOR_STATION=0" << endl;
      num_err++;
   }


    return num_err;

}
