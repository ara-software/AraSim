/*! \file araIcrrDefines.h
    \brief The standard ARA include file.
    
    It will end up including all those definitions that are needed 
    all over the shop. Hopefully there won't be too many of these floating 
    around.

    
    June 2010  ped@m2-consult.com , 
    Nov 2010 rjn@hep.ucl.ac.uk
    Copyright: M2 Consulting, 2010, sod off is it Copyright M2 Consulting I wrote most of this for ANITA
*/


#ifndef ARA_DEFINES_H
#define ARA_DEFINES_H

// Hardware stuff
#define RUBIDIUM_FREQUENCY 280000000 //280MHz
#define ICRR1_CLOCK_PERIOD 25 //ns
#define ICRR1_CLOCK_CHANNEL 8 // counting from 0
#define CHANNELS_PER_LAB3 9
#define MAX_NUMBER_SAMPLES_LAB3 260
#define LAB3_PER_ICRR 3
#define NUM_DIGITIZED_ICRR_CHANNELS (LAB3_PER_ICRR*CHANNELS_PER_LAB3)
#define EFFECTIVE_LAB3_SAMPLES 256
#define ADC_MAX 4096


//New stuff rjn added
#define ANTS_PER_ICRR 16
#define RFCHANS_PER_ICRR ANTS_PER_ICRR
#define MAX_RFCHANS_PER_ICRR 20

#define TRANS_PER_ICRR 6
#define TOTAL_ANTS_PER_ICRR ANTS_PER_ICRR+TRANS_PER_ICRR

//New stuff added by jpd
#define RFCHANS_TESTBED 16
#define RFCHANS_STATION1 20

//Legacy defines that will eventually disappear
#define CHANNELS_PER_CHIP CHANNELS_PER_LAB3 
#define MAX_NUMBER_SAMPLES MAX_NUMBER_SAMPLES_LAB3 
#define ACTIVE_CHIPS LAB3_PER_ICRR 
#define EFFECTIVE_SAMPLES EFFECTIVE_LAB3_SAMPLES 
#define ANTS_PER_STATION ANTS_PER_ICRR 
#define RFCHANS_PER_STATION RFCHANS_PER_ICRR 
#define TRANS_PER_STATION TRANS_PER_ICRR 
#define TOTAL_ANTS TOTAL_ANTS_PER_ICRR 
#define NUM_DIGITIZED_CHANNELS NUM_DIGITIZED_ICRR_CHANNELS 

//To account for multiple ICRR type stations jd
#define ICRR_NO_STATIONS 2


#endif // ARA_DEFINES_H
