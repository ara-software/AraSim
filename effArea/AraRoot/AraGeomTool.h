//////////////////////////////////////////////////////////////////////////////
/////  AraGeomTool.h       ARA Geometry tool                             /////
/////                                                                    /////
/////  Description:                                                      /////
/////     The Ara class working out what is where                        /////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#ifndef ARAGEOMTOOL_H
#define ARAGEOMTOOL_H

//Includes
#include <TObject.h>
#include <TMath.h>
#include "araIcrrStructures.h"
#include "araIcrrDefines.h"
#include "AraAntennaInfo.h"
#include "AraStationInfo.h"

//sqlite includes
#ifndef __CINT__
#include <sqlite3.h>
#endif

//!  AraGeomTool -- The Ara Geometry and numbering tool
/*!
  The Ara geometry and numbering tool
  \ingroup rootclasses
*/
class AraGeomTool
{
 public:
   AraGeomTool(); ///< Default constructor
   ~AraGeomTool(); ///< Destructor

   //   AraAntennaInfo *getAntByRfChan(int chan);
   //   AraAntennaInfo *getAntByPolAndAnt(AraAntPol::AraAntPol_t antPol, int antNum);
   int getChanIndex(AraLabChip::AraLabChip_t chip, int chan) {return chip*CHANNELS_PER_LAB3 +chan;}

   //something wrong with this guy
   AraLabChip::AraLabChip_t getLabChipForChan(int chan, int stationId) {return fStationInfo[stationId].fAntInfo[chan].labChip;}

   int getNumLabChansForChan(int chan, int stationId) { return fStationInfo[stationId].fAntInfo[chan].numLabChans;}
   int getFirstLabChanForChan(int chan, int stationId) { return fStationInfo[stationId].fAntInfo[chan].labChans[0];}
   int getSecondLabChanForChan(int chan, int stationId) { return fStationInfo[stationId].fAntInfo[chan].labChans[1];}


   //  int getFirstLabChanIndexForChan(int chan) { return getChanIndex(getLabChipForChan(chan),getFirstLabChanForChan(chan));}


   int getFirstLabChanIndexForChan(int chan, int stationId) { return getChanIndex(getLabChipForChan(chan, stationId),getFirstLabChanForChan(chan, stationId));}


   int getSecondLabChanIndexForChan(int chan, int stationId) { return getChanIndex(getLabChipForChan(chan, stationId),getSecondLabChanForChan(chan, stationId));}

   //jpd helperfunction for diplexed channels
   int isDiplexed(int chan, int stationId) {return fStationInfo[stationId].fAntInfo[chan].isDiplexed;}

   Double_t getLowPassFilter(int chan, int stationId) { return fStationInfo[stationId].fAntInfo[chan].lowPassFilterMhz; }

   Double_t getHighPassFilter(int chan, int stationId) { return fStationInfo[stationId].fAntInfo[chan].highPassFilterMhz; }



   int getRFChanByPolAndAnt(AraAntPol::AraAntPol_t antPol, int antNum, int stationId);


   //jpd this is a hack to try and get AraCanvasMaker.cxx to work 
   int getRFChanByPolAndAnt(AraAntPol::AraAntPol_t antPol, int antNum);
   
   

   Double_t calcDeltaTInfinity(Double_t ant1[3], Double_t ant2[3],Double_t phiWave, Double_t thetaWave);
   Double_t calcDeltaTR(Double_t ant1[3], Double_t ant2[3], Double_t phiWave, Double_t thetaWave,Double_t R);

   Double_t calcDeltaTInfinity(Int_t chan1, Int_t chan2,Double_t phiWave, Double_t thetaWave, int stationId);
   Double_t calcDeltaTR(Int_t chan1, Int_t chan2, Double_t phiWave, Double_t thetaWave,Double_t R, int stationId);

   //Instance generator
   static AraGeomTool*  Instance();

   AraStationInfo fStationInfo[ICRR_NO_STATIONS]; //station info contains the antenna info and station information
   int fAntLookupTable[ICRR_NO_STATIONS][3][8]; //At some point should lose the magic numbers
   
   //Some variables to do with ice properties
   static Double_t nTopOfIce;

  
 protected:
   static AraGeomTool *fgInstance;  
   // protect against multiple instances

 private:
   //jpd this will be the implementation that will load from the sql DB
   void readChannelMapDb(Int_t stationId);

};


#endif //ARAGEOMTOOL_H
