//////////////////////////////////////////////////////////////////////////////
/////  AraAntennaInfo.h       ARA Antenna Information                    /////
/////                                                                    /////
/////  Description:                                                      /////
/////     A simple class for storing information about an Ara Antenna    /////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#include "AraAntennaInfo.h"
#include <iostream>
#include <fstream>
#include <cstring>
ClassImp(AraAntennaInfo);

const char *AraAntType::antTypeAsString(AraAntType::AraAntType_t antType)
{
  switch (antType) {
  case kBicone: return "Bicone";
  case kBowtieSlot: return "Bowtie-Slotted Cylinder";
  case kDiscone: return "Discone";
  case kBatwing: return "Batwing";
  case kFatDipole: return "Fat Dipole";
  case kQuadSlot: return "Quad-Slot Cylinder";
  default: return "Unknown";
  } 
}


const char *AraAntPol::antPolAsString(AraAntPol::AraAntPol_t antPol)
{
  switch (antPol) {
  case kVertical: return "Vertical";
  case kHorizontal: return "Horizontal";
  case kSurface: return "Surface";
  default: return "Unknown";
  } 
}


const char *AraLabChip::labChipAsString(AraLabChip::AraLabChip_t labChip)
{
  switch (labChip) {
  case kA: return "A";
  case kB: return "B";
  case kC: return "C";
  default: return "Unknown";
  } 
}



AraAntennaInfo::AraAntennaInfo() 
{
   //Default Constructor
}

AraAntennaInfo::~AraAntennaInfo() {
   //Default Destructor
}

void AraAntennaInfo::printAntennaInfo()
{
  std::cout << "*************************************************************\n";
  std::cout << "Antenna Info for Channel " << chanNum << "\n";
  std::cout << designator << " at " << locationName << "\n";
  std::cout << "DAQ Chan : " << getDaqBoxChan() << "\n";
  if(numLabChans==2)
    std::cout << "Lab Chans : " << AraLabChip::labChipAsString(labChip) << labChans[0]+1 << "," << AraLabChip::labChipAsString(labChip) << labChans[1]+1 << "\n";
  else
    std::cout << "Lab Chan : " << AraLabChip::labChipAsString(labChip) << labChans[0]+1 << "\n";
    
  std::cout << "Filters " << highPassFilterMhz << "-" << lowPassFilterMhz << " MHz\n";
  std::cout << AraAntType::antTypeAsString(antType) << "\t" << AraAntPol::antPolAsString(polType) << " polarisation\n";
  std::cout << antLocation[0] << "," << antLocation[1] << "," << antLocation[2] << " m\n";
  std::cout << "Delay " << cableDelay << " ns\n";
  std::cout << "Average Noise Figure " << avgNoiseFigure << " K\n";
  std::cout << "*************************************************************\n";

}


const char *AraAntennaInfo::getDaqBoxChan()
{
  static char boxChan[5];
  if(daqChanType==AraDaqChanType::kDisconeChan)
    sprintf(boxChan,"DIS%d",daqChanNum);
  else if(daqChanType==AraDaqChanType::kBatwingChan)
    sprintf(boxChan,"Bat%d",daqChanNum);
  return boxChan;
}

