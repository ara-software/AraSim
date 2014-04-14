//////////////////////////////////////////////////////////////////////////////
/////  AraAntennaInfo.h       ARA Antenna Information                    /////
/////                                                                    /////
/////  Description:                                                      /////
/////     A simple class for storing information about an Ara Antenna    /////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#ifndef ARAANTENNAINFO_H
#define ARAANTENNAINFO_H

//Includes
#include <TObject.h>
#include "araIcrrStructures.h"
#include "araIcrrDefines.h"
//#include "RawIcrrStationEvent.h"
//#include "AraEventCalibrator.h"

//!  AraAntennaInfo -- ARA Antenna Information
/*!
  A simple class for storing information about an Ara Antenna
  \ingroup rootclasses
*/

namespace AraAntType {
  typedef enum EAraAntType {
    kBicone = 1,
    kBowtieSlot = 2,
    kDiscone = 3,
    kBatwing = 4,
    kFatDipole =5,
    kQuadSlot = 6
  } AraAntType_t;
  const char *antTypeAsString(AraAntType::AraAntType_t antType); ///< Returns the antenna type as a string

}

namespace AraAntPol {
  typedef enum EAraAntPol {
    kVertical = 0,
    kHorizontal = 1,
    kSurface = 2
  } AraAntPol_t;
  const char *antPolAsString(AraAntPol::AraAntPol_t antPol); ///<Returns the antenna polarisation as a string
}

namespace AraDaqChanType {
  typedef enum EAraDaqChanType {
    kDisconeChan =1,
    kBatwingChan =2
  } AraDaqChanType_t;
}

namespace AraLabChip {
  typedef enum EAraLabChip {
    kA = 0,
    kB = 1,
    kC = 2
  } AraLabChip_t;
  const char *labChipAsString(AraLabChip::AraLabChip_t labChip);
}

namespace AraAntDir {
  typedef enum EAraAntDir {
    kReceiver = 1,
    kTransmitter = 2
  } AraAntDir_t;
}
    
namespace AraSurfaceOrientation {
  typedef enum EAraSurfaceOrientation {
    kNorthSouth =1,
    kEastWest =2
  } AraSurfaceOrientation_t;
}

class AraAntennaInfo: public TObject
{
 public:
   AraAntennaInfo(); ///< Default constructor
   ~AraAntennaInfo(); ///< Destructor

   void printAntennaInfo();
   const char *getDaqBoxChan();
   Int_t chanNum;
   AraDaqChanType::AraDaqChanType_t daqChanType;
   Int_t daqChanNum;
   Double_t highPassFilterMhz;
   Double_t lowPassFilterMhz;
   Int_t daqTrigChan;
   Int_t numLabChans;
   AraLabChip::AraLabChip_t labChip;
   Int_t labChans[2]; ///<These will count from 0
   Int_t isDiplexed; ///jpd
   Int_t diplexedChans[2]; ///jpd


   Int_t preAmpNum;
   Double_t avgNoiseFigure;
   Int_t rcvrNum;
   char designator[3];
   Int_t antPolNum;
   AraAntType::AraAntType_t antType;
   AraAntPol::AraAntPol_t polType;
   char locationName[4];
   Double_t antLocation[3]; ///< x,y,z in m
   Double_t cableDelay; ///< In ns
   AraAntDir::AraAntDir_t antDir;
   AraSurfaceOrientation::AraSurfaceOrientation_t antOrient; ///<Only for surface antennas
   
   Double_t debugHolePosition[3]; ////< x,y,z in m
   Double_t debugPreAmpDz; ///< in m
   Double_t debugHolePositionZft; ///< in ft
   Double_t debugHolePositionZm; ///< in m
   Double_t debugTrueAsBuiltPositon[3]; ///< x,y,z in m
   Double_t debugCableDelay2; //in ns
   Double_t debugFeedPointDelay; //in ns
   Double_t debugTotalCableDelay; //in ns
     
       

  ClassDef(AraAntennaInfo,1);
};


#endif //ARAANTENNAINFO_H
