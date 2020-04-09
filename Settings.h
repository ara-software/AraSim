#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <string>
#include "TObject.h"
#include "AraSimVersion.h"

class Detector;

using namespace std;
//using std::string;

class Settings 
{
    protected:

    public:

        Settings();
        ~Settings();
        void Initialize();
        void ReadFile(string setupfile);
        void ReadEvtFile(string evtfile);

        int CheckCompatibilitiesSettings();// check if settings are not compatible to each other
        int CheckCompatibilitiesDetector(Detector *detector);// check if settings are not compatible to each other. checking against initialized Detector *detector object

	int ARASIM_VERSION_MAJOR;
        int ARASIM_VERSION_MINOR;
        int ARASIM_VERSION_SUBMINOR;
        double ARASIM_VERSION;

        double ARAROOT_VERSION;
	bool ARAUTIL_EXISTS;
        int NNU;

  // NEED TO FIGURE OUT A GOOD WAY TO READ THIS IN AND STORE THEM.
  // INPUT FILE AGAIN?  SOMETHING ELSE?
  //These were moved here from IceModel under the new compilation scheme
        int ICE_MODEL; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP.
        int NOFZ; // 1=depth dependent index of refraction,0=off
        int CONSTANTCRUST; // set crust density and thickness to constant values.
        int CONSTANTICETHICKNESS; // set ice thickness to constant value
        int FIXEDELEVATION; // fix the elevation to the thickness of ice.
        int MOOREBAY; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data
        int USE_ARA_ICEATTENU; // 0 : use old ice attenuation factor with one depth info, 1 : (default) use ARA measured ice attenuation factor with depth from ray steps

        double EXPONENT; // 10^19 eV neutrinos only

        int DETECTOR;   // choose detector layout

	int DETECTOR_STATION; // for DETECTOR=4, indicates the single station to be simulated
	                      // 0 = testbed, 1 = A1, 2 = A2, 3 = A3

	int number_of_stations; // the number of stations to be used in the simulation

        int INTERACTION_MODE;   // method to choose interaction point posnu. 0 : PickUnbiased, 1 : PickNear, 2 : PickExact, 3 : PickAboveIce

        double POSNU_RADIUS;    //PickNear radius in meter

        int WHICHPARAMETERIZATION;  //

        int SIMULATION_MODE;    // 0 : old freq domain mode, 1: new time domain mode

        int EVENT_TYPE;         // 0 : neutrino only events,  1 : blackhole evnet? ... etc

        int WAVE_TYPE;          // 0 : plane wave,  1 : spherical wave

        int LPM;                // 0 : no LPM effect, 1 : with LPM effect(default)

        int SECONDARIES;        // 0 : no secondary interaction, 1 : with secondary interactions

        int TAUDECAY;           // 0 : don't account tau decay as secondaries, 1 : let tau decay as secondaries

        double TIMESTEP;        // time step after fft. t domain bin width

        double PHASE;           // phase factor for fft. default 90 deg (in deg !!)

        int NFOUR;              // number of total bins for FFT. has to be power of 2 values

        int NOISE;              // noise condition settings degault 0 ( : thermal flat noise), 1 : Rayleigh dist. fit for installed TestBed geom, 2: Noise figure values for station 2 (??) from Thomas Meures 2015/2016

        int ATMOSPHERE;         // include atmosphere 1, no 0

        int TRIG_SCAN_MODE;
        
        double POWERTHRESHOLD;  // power threshold value. default -4.41 (same with icemc powerthreshold)

        double MAXT_DIODE;      // diode model max time, default : 70.e-9s

        int IDELAYBEFOREPEAK_DIODE; // diode model bins before the peak, default : (int) 13.e-9 / TIMESTEP = 33

        int IWINDOW_DIODE;            // diode model interesting bins after the peak, default : (int) 4.e-9 / TIMESTEP = 10

        int DATA_BIN_SIZE;          // bin size which mimic the data (time delay between antennas), default : 16384

        double NOISE_TEMP;          // noise temperature (default : 325 K = 230K (ice) + 95K (receiver), from Peter's spreadsheet)

        int TRIG_ANALYSIS_MODE; // trigger mode 0 : signal + noise, 1 : only pure signal, 2 : only pure noise, default : 0

        double TRIG_TIMEOUT;    // time out after the trigger (we have to wait this amount of time to trig next event), default : 1us = 1.E-6
        double TRIG_WINDOW;     // coincidence time window for trigger default : 250 ns

        int NOISE_EVENTS;       // number of noise events which will be stored in Trigger class for later use. This will also used to calculate mean, rms noise (with diode convlv). default : 1000

        int DATA_SAVE_MODE;     // 0 : save all information which are generated during the processes. 1 : light mode, remove most of data before the final global trigger (most of data except geometric info, final data V_mimic will be removed). 2 : light mode 2, except for physics information (energy, geometric info), all waveform information are removed (include noise waveforms in trigger class)

        int N_TRIG;         // number of coincidence channels triggered in TRIG_WINDOW default : 3

        int RANDOM_MODE;    // 0 : same random number generate, 1 : TRandom3(0) used. purely randomly generated (seed is guaranteed to be unique in space and time)
    
        int SEED; // default: 1, only applies if RANDOM_MODE=0, provides base seed value and run number taken from arguments is added to this value in order to submit multiple repeatable runs instead of only one single long repeatable run

        int BORE_HOLE_ANTENNA_LAYOUT;   // 0 = (V-H-V-H), 1 = (V-H-V), 2 = (V-H-V-V), 3 = (V-H-H-H), 4 = (V-H-H) default : 0
    
        int DATA_LIKE_OUTPUT; // Formerly WRITE_ALL_EVENTS, the mode numbering has changed as well
	                      // 0 - don't write any information to the data-like output tree
                              // 1 - only write globally triggered events,
                              // 2 - Write all event events including events that are not globally triggered
        // Note: NNU is the number of neutrinos that have been thrown in total, not just globally triggered events
        // When writing all events, the waveform stored in UsefulAraSimEvent->VoltsRF[] is just the untriggered, noiseless waveform of the initial signal, and it has not propagated to the antenna yet.

        double RAYSOL_RANGE;    // direct distance limit to do raysolver. If distance between posnu and antenna is bigger than RAYSOL_RANGE, AraSim will not try RaySol for that event. Default : 3000 m

        int CALPULSER_ON; // 0: no calpulsers in event list, 1: only throws calpulser 1 events, 2: only throws calpulser 2 events, 3: throws both calpulser 1 and 2 events alternating between them, 4: throws calpulser 1 and 2 events integrated with the simulated data (not yet implemented)
    
        int TESTBED_ON;
    
        int READGEOM;
    
        int V_MIMIC_MODE; // default : 0 - write out all chs where global triggered bin is center of the window
                        // 1 - same as above 0 mode but apply TestBed ch delay - average BH ch delay
                        // 2 - same as above 0 mode but apply TestBed ch delay - average BH ch delay + additional delay to match with actual TestBed data waveforms
    
        int USE_INSTALLED_TRIGGER_SETTINGS; // default : 0 - use idealized settings for the trigger
        //other options:  1 - use trigger settings for installed stations, i.e. trigger window, etc.
    
        int NUM_INSTALLED_STATIONS; // the number of stations including the testbed that are in reality installed in the ice and have position and electronics information
    
        int PICK_POSNU_DEPTH;  // whether use MAX_POSNU_DEPTH or not. 0 : pick posnu depth full ice depth, 1 : pick posnu depth only MAX_POSNU_DEPTH

        double MAX_POSNU_DEPTH;  // maximum posnu depth when above PICK_POSNU_DEPTH=1

        int NNU_THIS_THETA;     // if nnu theta angle will be selected randomly from [0, PI] (default=0) or set nnu theta to near some angle (=1). Theta = 0 corresponds to an upgoing event and Theta = pi corresponds to a downgoing event.

        double NNU_THETA;       // nnu theta when NNU_THIS_THETA=1. In radians.

        double NNU_D_THETA;     // nnu theta variation from NNU_THETA, when NNU_THIS_THETA=1 case.  Use radians.

        int NNU_THIS_PHI; // default 0: nnu angle pure random, 1: set a specific phi. Phi is in the range [0,2 Pi], where 0 is aligned with the iceflow and it's goes counterclockwise.

        double NNU_PHI;// default : nnu phi: 45 deg. Use radians.

        double NNU_D_PHI;// default : nnu_d_phi : 5 deg  Use radians.

        double CALPUL_OFFCONE_ANGLE;    // for calpulser events, what's the offcone angle value?

        double CALPUL_AMP;    // for calpulser events, how strong the calpulser waveforms?

        int TRIG_ONLY_BH_ON;    // if trigger will occur with all chs (0, default) or only borehole chs (1)

        int TRIG_THRES_MODE;    // if trigger threshold will use no specific offset value (0, default). or use data/threshold_offset.csv file values as a threshold offset for each chs

        int NOISE_CHANNEL_MODE;    // if using same noise temp for all chs (0, default), using diff temp for all chs (1), using diff temp for first 8 chs and share same temp for other chs (2), using same noise temp but allowing to use noise

        double CONST_MEANDIODE;
    
        double CONST_RMSDIODE;  // in case NOISE_CHANNEL_MODE = 1, just using this CONST_RMSDIODE value for threshold


        int USE_TESTBED_RFCM_ON;    // use RFCM measurement for testbed or not (default 0)
    
        double RFCM_OFFSET;  // if above USE_TESTBED_RFCM_ON = 1, we need RFCM attenuator factor cancel (default 80)


        int USE_MANUAL_GAINOFFSET; //whether use additional gain offset after all ant, elect chain. default 0 : don't apply any additional gain, 1 :  use gain offset
    
        int USE_CH_GAINOFFSET; // only if USE_MANUAL_GAINOFFSET = 1, default 0 : apply constant gain to all channels which is MANUAL_GAINOFFSET_VALUE, 1 : apply ch gain offset by using data/preamp_gain_offset.csv file (only installed TestBed mode available)


        double MANUAL_GAINOFFSET_VALUE; // gain offset value only if USE_MANUAL_GAINOFFSET = 1, and USE_CH_GAINOFFSET = 0


        int NOISE_WAVEFORM_GENERATE_MODE; // default 0 : generate new noise waveforms for each evts, if you set other values, noise waveforms will be generated once and use them for all evts



        int GETCHORD_MODE; // which Getchord function to use. default 0 : old Getchord function (not correct tau weight, weight don't have ice inside interaction probability in it). 1 : new Getchord from icemc. This has new tau weight calculation and ice interaction probability applied to weight factor.

        int taumodes; // taumodes 1 : tau created in the rock

        int BH_ANT_SEP_DIST_ON; // if we are going to use separate bore hole antenna distance or not. By default it's 0 (don't use separate dist)

        int TRIG_MODE; // default 0 : if any antennas got passed N_TRIG or more, global trig. 1 : either Vpol or Hpol antennas got passed N_TRIG_V or N_TRIG_H respectively, global trig.

        int N_TRIG_V; // default : 3 (3 out of Vpolchannels in a station)
    
        int N_TRIG_H; // default : 3 (3 out of Hpol channels in a station)

        int FILL_TREE_MODE; // default 0 : fill tree for all events, 1 : fill tree only usable posnu events, 2 : fill tree only trigger passed events

	int ONLY_PASSED_EVENTS; // 0: NNU represents the number of neutrinos total thrown in the simulation, 1: NNU represents the number of globally triggered neutrinos desired in the output file in the end
	int NNU_PASSED; // Number of neutrinos allowed to pass the trigger - for ONLY_PASSE_EVENTS, loop until this is reached, otherwise just a storage value

        double PICKNEARUNBIASED_R; // radius of the sphere surrounding the detector for INTERACTION_MODE=3, current default value is 5000 m


        int SHOWER_MODE; // in time domain caes (SIMULATION_MODE=1), if we just want EM shower portion, SHOWER_MODE = 0, HAD shower only = 1, if we want summed EM and HAD = 2, either one which is dominant = 3 (default = 2)

        double SHOWER_STEP; // step size (meter) in generating shower profile, default = 0.001

        int SHOWER_PARAM_MODEL; // shower profile parameter, 0 = fit function from Jaime, 1 = fit from Carl (OSU), default 0

        double OFFCONE_LIMIT; // offcone cut angel (deg) to reduce computation time 


        int ALL_ANT_V_ON; // default 1 : use both Vpol, Hpol ant gain from Vpol ant model but for Hpol ant case response is in Hpol, not Vpol, 0 : use Vpol ant gain from Vpol model, Hpol ant gain from Hpol model (where Hpol ant model is not reliable).

    
        int PHASE_SKIP_MODE; // skip applying phase in t-domain mode (SIMULATION_MODE = 1). default 0 : don't skip (apply all phase), 1 : only upto Askaryan radiation, 2 : only upto antenna


        int DEBUG_MODE_ON; // default 0 : do as usual, 1 : debug mode on -> skip noise generation, signal generating, trig analysis portion upto DEBUG_SKIP_EVT

        int DEBUG_SKIP_EVT; // default 0, number of events that skip noise, signal generating, and trigger analysis part


        double V_SATURATION; // default 1., output voltage saturation due to amplifier limit. by default it's +-1. V


    
        int ADDITIONAL_DEPTH_ON; // whether add more depth to each antenas, default : 0 (not on)
    
        double ADDITIONAL_DEPTH; // additional depth value, default 100, only applied when ADDITIONAL_DEPTH_ON=1


        
        int TRIG_ONLY_LOW_CH_ON;    // if trigger will occur with all chs (0, default) or only lower 8 chs (1)

    
        double ACCUM_TRIG_SEARCH_BINS_STATION0; // not actually setting value but gives us how much trigger searched bins there were in the run


	int OUTPUT_TDR_GRAPH;// saves a few example graphs of the tunnel diode response for a triggered event
	
	
    
        int NU_NUBAR_SELECT_MODE; // default : 3 = random nu_nubar based on arXiv:1108.3163, section 3, 0 = just nu, 1 = just nubar 


        int SELECT_FLAVOR; // default : 0 = random 1:1:1, 1: e, 2: mu, 3: tau

        int SELECT_CURRENT; //default:2:random, 0:nc, 1:cc


        int AVZ_NORM_FACTOR_MODE; // default : 1 : don't apply sqrt(2) (actually applied but cancel that) as realft assume Hn as double-sided spectrum (invFFT normalization factor 2/N) and also remove dF binning factor in MakeArraysforFFT function, 0 : use normalization factors like in old version

	int RAY_TRACE_ICE_MODEL_PARAMS; // which parameter set is used for the exponential ice model (defined in RayTrace_IceModel.cc) 
	//0 : default, South Pole model (fitted from RICE data
	//1 : South Pole model fitted from RICE #2
	//2 : South Pole (Eisen (2003))
	//3 : South Pole (Gow)
	//10 : Moore's Bay Model 1
	//11 : Moore's Bay Model 2
	//20 : Byrd (Ebimuna (1983))
	//30 : Mizuho (Ebimuna (1983))

	int WAVEFORM_LENGTH; // the number of samples in the waveform length for V_mimic and UsefulAtriStationEvent, default: 64/2*20 = 640

	int WAVEFORM_CENTER; // the relative location of the center of the write-out window with respect to the last triggered bin (which is laced at the center of the window by default), this effectively provides a global delay in the write-out window across all channels: positive values shift the write-out window to later times in the waveform, negative values shift the window to earlier times, default: 0

	double POSNU_R; // default: 1000; meters from station center
	double POSNU_THETA; // default: -PI/4; elevation angle from station center coordinates
	double POSNU_PHI; // default: 0; azimuth angle from station center coordinates
	// Only using INTERACTION_MODE = 2, these values set the position of the interaction point; can be used in EVENT_TYPE = 0 and 10

	double ARBITRARY_EVENT_ATTENUATION; // attenuation of the waveform intensity for arbitrary event waveforms
	double PICK_ABOVE_HEIGHT;

        int EVENT_GENERATION_MODE;//default: 0: not event mode, 1: event mode
	//        int EVENT_NUM;//number of events to read in, using EVENT_GENERATION_MODE=1, default initial 0: read in the number of events in the input file, resets to total number of events read in

	int ANTENNA_MODE; // 0: old default antenna models bicone/rotated dipole
	                   // 1: using different antenna response for the top Vpol antennas, otherwise same as old default

	int APPLY_NOISE_FIGURE; // 0: do not apply new noise figure from Thomas Meures 2016
	                        // 1: apply new noise figure to data

	int CUSTOM_ELECTRONICS; //0 (default): use the regular "ARA_Electronics_TotalGain_TwoFilter.txt" file
							//1 : load a custom electronics file, stored as "custom_electronics.txt" in the `data` directory


//arrays for saving read in event features in EVENT_GENERATION_MODE=1
        vector<int> EVID;
        vector<int> NUFLAVORINT;
        vector<int> NUBAR;
        vector<double> PNU;
        vector<int> CURRENTINT;
        vector<double> IND_POSNU_R;
        vector<double> IND_POSNU_THETA;
        vector<double> IND_POSNU_PHI;
        vector<double> IND_NNU_THETA;
        vector<double> IND_NNU_PHI;
        vector<double> ELAST;



    // below : values from icemc
    
    
    int UNBIASED_SELECTION;
    int WHICH; // which payload to use 0=Anita-lite,1=Ross,2=Smex,3=make your own
    int CYLINDRICALSYMMETRY; // is it cylindrically symmetric =1 if which=1,2, =0 if which=0
    // if which=3 then 0 or 1
    double SIGMA_FACTOR; // factor to multiply cross section by for error analysis
    int SIGMAPARAM; // 0=Reno, 1=Connolly et al. 2011 for cross section parametrization
    int YPARAM; // 0=Reno, 1=Connolly et al. 2011 for cross section parametrization, 2 = a specific elast_y
    double ELAST_Y;

    int SIGMA_SELECT; // when in SIGMAPARAM=1 case, 0 : (default) use mean value, 1 : use upper bound, 2 : use lower bound

    int SIGNAL_FLUCT;  // 1=add noise fluctuation to signal or 0=do not
    int TRIGGERSCHEME;  // frequency domain voltage, frequency domain energy, time domain diode integration
    int ZEROSIGNAL;  // zero the signal to see how many of our hits are noise hits
    int REMOVEPOLARIZATION; //Disable polarizations
    
    int EVENTSMAP;//whether draw the events distribution map
    
    int WHICHRAYS;  // how many rays to look at (1) direct only (2) direct and down-going.

// trigger
int LCPRCP; // 1 for circular polarization trigger, 0 for V and H
int JUSTVPOL; // 0 for both polarizations, 1 for just V polarization
// doesn't allow for both LCPRCP=1 and JUSTVPOL=1
//int FIFTHBAND; // 1 to include 0.2-1.2 GHz as a frequency band if JUSTVPOL==1
//int NFOLD=3;  // how many channels must pass the trigger - in old mechanism - only used for anita-lite
int NFOLD;  // how many channels must pass the trigger - in old mechanism - only used for anita-lite


//int CHMASKING=1; // whether or not to include channel masking
//int PHIMASKING=1; // whether or not to include phi masking  
int CHMASKING; // whether or not to include channel masking
int PHIMASKING; // whether or not to include phi masking  

//int NLAYERS=0;
//int NANTENNAS=0;

int NLAYERS;
int NANTENNAS;

/* int ONLYFINAL=1; // only write to final histogram */
/* int HIST_MAX_ENTRIES=10000; //maximum number of events to put in histograms */
/* int HIST=1;          //write to histograms  */

int ONLYFINAL; // only write to final histogram
int HIST_MAX_ENTRIES; //maximum number of events to put in histograms
int HIST;          //write to histograms 
double BW; // BANDWIDTH
//int DISCONES=1; // whether or not to use discones
int DISCONES; // whether or not to use discones

//double NDISCONES_PASS=3; // number of discones needed to pass
double NDISCONES_PASS; // number of discones needed to pass

int BORESIGHTS; // whether to loop over boresights
int SLAC; // whether or not we are simulating the slac run
double SLACSLOPE; // slope of the ice
double SLACICELENGTH;  // length of the block of ice
double SLAC_HORIZDIST; // horizontal distance from interaction to center of payload at slac beam test
double SLAC_DEPTH; // vertical depth of interaction at slac beam test
double SLAC_HORIZ_DEPTH; // horizontal depth of interaction at slac

int ROUGHNESS; // include effects of surface roughness
int FIRN; // whether or not to include the firn

//int SLOPEY=1; // 1=slopeyness on, 0=slopeyness off
//double SLOPEYSIZE=0.012; // This determines size of the slopeyness (0.10=5.4, 0.20=7.4 deg mean)

int SLOPEY; // 1=slopeyness on, 0=slopeyness off
double SLOPEYSIZE; // This determines size of the slopeyness (0.10=5.4, 0.20=7.4 deg mean)

 bool DEBUG;
string outputdir; // directory where outputs go

//double THERMALNOISE_FACTOR=1.0; // factor to multiply thermal noise for error analysis
double THERMALNOISE_FACTOR; // factor to multiply thermal noise for error analysis

//double FREQ_LOW_SEAVEYS=200.E6; // min frequency for seaveys
//const double FREQ_HIGH_SEAVEYS=1200.E6; // max frequency for seaveys

double FREQ_LOW_SEAVEYS; // min frequency for seaveys
double FREQ_HIGH_SEAVEYS; // max frequency for seaveys
 double BW_SEAVEYS;
 //int FORSECKEL=1; // Make array of strength of signal across frequencies for different viewing angles.
int FORSECKEL; // Make array of strength of signal across frequencies for different viewing angles.

double ROUGHSIZE; // roughness size


/* int ICE_MODEL=0; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP. */
/* int NOFZ=1; // 1=depth dependent index of refraction,0=off */
/* int CONSTANTCRUST=0; // set crust density and thickness to constant values. */
/* int CONSTANTICETHICKNESS=0; // set ice thickness to constant value */
/* int FIXEDELEVATION=0; // fix the elevation to the thickness of ice. */
/* int MOOREBAY=0; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data */
/* int USEPOSITIONWEIGHTS=1;// whether or not to restrict the neutrino position so it is within the horizon of the balloon */
/* int WRITE_FILE=0; //Select whether or not to write a new input file for CreateHorizons */

int USEPOSITIONWEIGHTS;// whether or not to restrict the neutrino position so it is within the horizon of the balloon
int WRITE_FILE; //Select whether or not to write a new input file for CreateHorizons

int MINRAY;
int MAXRAY;

int horizontal_banana_points;
  int vertical_banana_points;




 // end of values from icemc


  ClassDef(Settings,1);


};
#endif

