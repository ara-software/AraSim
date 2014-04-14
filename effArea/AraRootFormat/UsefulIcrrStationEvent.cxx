//////////////////////////////////////////////////////////////////////////////
/////  UsefulIcrrStationEvent.cxx        ARA header reading class                  /////
/////                                                                    /////
/////  Description:                                                      /////
/////     A simple class that reads in useful ARA headers and produces     ///// 
/////   calibrated time and voltage stuff                                /////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#include "UsefulIcrrStationEvent.h"
//#include "FFTtools.h"
//#include "AraGeomTool.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
#include <cstring>
ClassImp(UsefulIcrrStationEvent);

UsefulIcrrStationEvent::UsefulIcrrStationEvent() 
{
   //Default Constructor
//  fCalibrator=0;
}


UsefulIcrrStationEvent::~UsefulIcrrStationEvent() {
   //Default Destructor
}


//UsefulIcrrStationEvent::UsefulIcrrStationEvent(RawIcrrStationEvent *rawEvent, AraCalType::AraCalType_t calType)
//  :RawIcrrStationEvent(*rawEvent)
//{
//  fCalibrator=AraEventCalibrator::Instance();
//  fCalibrator->calibrateEvent(this,calType);  
//}
/*
TGraph *UsefulIcrrStationEvent::getGraphFromElecChan(int chan)
{
  if(chan<0 || chan>=NUM_DIGITIZED_ICRR_CHANNELS)
    return NULL;
  return new TGraph(fNumPoints[chan],fTimes[chan],fVolts[chan]);
}

TGraph *UsefulIcrrStationEvent::getGraphFromRFChan(int chan)
{
  if(chan<0 || chan>=RFCHANS_PER_ICRR)
    return NULL;
  return new TGraph(fNumPointsRF[chan],fTimesRF[chan],fVoltsRF[chan]);
}
/*
TGraph *UsefulIcrrStationEvent::getFFTForRFChan(int chan)
{

   //   static AraGeomTool *fGeomTool = AraGeomTool::Instance();
   TGraph *gr = getGraphFromRFChan(chan);
   if(!gr) return NULL;
   Double_t newX[512],newY[512];
   Double_t intSample=1;
   Int_t maxSamps=256;
   // if(fGeomTool->getNumLabChansForChan(chan)==2) {
   intSample=0.5;
   maxSamps=512;
   //   }
   

   TGraph *grInt = FFTtools::getInterpolatedGraph(gr,intSample);
   

   Int_t numSamps  = grInt->GetN();
   Double_t *xVals = grInt->GetX();
   Double_t *yVals = grInt->GetY();
   for(int i=0;i<maxSamps;i++) {
      if(i<numSamps) {
	 newX[i]=xVals[i];
	 newY[i]=yVals[i];
      }
      else {
	 newX[i]=newX[i-1]+intSample;  
	 newY[i]=0;
      }      
  }
   TGraph *grNew = new TGraph(maxSamps,newX,newY);
   TGraph *grFFT = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(grNew);
   delete gr;
   delete grNew;
   delete grInt;
   return grFFT;
}

TH1D *UsefulIcrrStationEvent::getFFTHistForRFChan(int chan)
{

   Int_t numBins=256;
   Double_t minX=0.0;
   Double_t maxX=1000.0;
   minX = minX - ( (maxX-minX)/numBins/2.0 ); // adjust histogram edges so that the bin centers
   maxX = maxX + ( (maxX-minX)/numBins/2.0 ); // of the histograms are aligned with graph [add bdf]
   numBins++;
   char histName[180];
//   sprintf(histName,"%s_ffthist",this->GetName());
//   TH1D *histFFT = new TH1D(histName,histName,numBins,minX,maxX);
//   if(fillFFTHistoForRFChan(chan,histFFT)==0)
//      return histFFT;
   return NULL;

}
      
int UsefulIcrrStationEvent::fillFFTHistoForRFChan(int chan, TH1D *histFFT) 
{

   TGraph *grFFT =getFFTForRFChan(chan);
   if(!grFFT) return -1;
   Double_t *xVals=grFFT->GetX();
   Double_t *yVals=grFFT->GetY();
   Int_t numPoints=grFFT->GetN();
   for(int i=0;i<numPoints;i++) {    
      histFFT->Fill(xVals[i],yVals[i]);
   }
   delete grFFT;
   return 0;

}
 */
/*
int UsefulIcrrStationEvent::guessRCO(int chanIndex)
{
  ///< Looks at clock channel to try and guess which RCO phase we are in.
  int chip=chanIndex/CHANNELS_PER_LAB3;
  if(chip<0 || chip>=LAB3_PER_ICRR) return -1;
  static int firstTime=1;
  static unsigned int lastEventNumber=0; //random number
  static unsigned int lastUnixTimeUs=0; //random number
  static int fRcoGuess[LAB3_PER_ICRR]={0};
  fCalibrator=AraEventCalibrator::Instance();
  if(lastEventNumber!=this->head.eventNumber || firstTime || lastUnixTimeUs!=this->head.unixTimeUs) {
    fCalibrator->fillRCOGuessArray(this,fRcoGuess);
    lastEventNumber=this->head.eventNumber;
    lastUnixTimeUs=this->head.unixTimeUs;
    firstTime=0;
  }
  return fRcoGuess[chip];    
}

TGraph *UsefulIcrrStationEvent::getFFTForClock(int clock_number)
{

   if ( (clock_number<0) || ( clock_number>2) ) {
      return NULL;
   }

   int channel_number = clock_number*9+8;
   //   static AraGeomTool *fGeomTool = AraGeomTool::Instance();
   TGraph *gr = getGraphFromElecChan(channel_number);
   if(!gr) return NULL;
   Double_t newX[512],newY[512];
   Double_t intSample=1.0;
   Int_t maxSamps=256;
   intSample = 0.5;                   // interpolation spacing in ns                
   maxSamps  = 512;                   // number of samples in interpolated waveform 

   TGraph *grInt = FFTtools::getInterpolatedGraph(gr,intSample);

   Int_t numSamps  = grInt->GetN();
   Double_t *xVals = grInt->GetX();
   Double_t *yVals = grInt->GetY();
// printf("UsefulIcrrStationEvent::getFFTForClock() - intSample = %f - maxSamps = %d - numSamps = %d [%d]\n",intSample,maxSamps,numSamps,grInt->GetN());
   for(int i=0;i<maxSamps;i++) {
//    if ( i<numSamps ) {
//       printf("UsefulIcrrStationEvent::getFFTForClock() - i = %d - (x,y) = %.2f,%.2f",i,xVals[i],yVals[i]);
//    } else {
//       printf("UsefulIcrrStationEvent::getFFTForClock() - i = %d - i>numSamps ",i);
//    }
      if( i<numSamps ) {
//       printf(" - use point\n");
	 newX[i]=xVals[i];
	 newY[i]=yVals[i];
      }
      else {
	 newX[i]=newX[i-1]+intSample; /// change from +(1.);
         newY[i]=0;
//       printf(" - zero fill - (x,y) = %.2f,%.2f - previous_x = %.2f\n",newX[i],newY[i],newX[i-1]);
      }
   }
   TGraph *grNew = new TGraph(maxSamps,newX,newY);
   TGraph *grFFT = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(grNew);
   delete gr;
   delete grNew;
   delete grInt;
   return grFFT;

} // end of function UsefulIcrrStationEvent::getFFTForClock()

TH1D *UsefulIcrrStationEvent::getFFTHistForClock(int clock_number)
{

   Int_t numBins=256;
   Double_t minX=0.0;
   Double_t maxX=1000.0;   
   minX = minX - ( (maxX-minX)/numBins/2.0 ); // adjust histogram edges so that the bin centers
   maxX = maxX + ( (maxX-minX)/numBins/2.0 ); // of the histograms are aligned with graph [add bdf]
   numBins++;
   char histName[180];
   sprintf(histName,"%s_clock_%2.2dffthist",this->GetName(),clock_number);
   TH1D *histFFT = new TH1D(histName,histName,numBins,minX,maxX);
   if(fillFFTHistoForClock(clock_number,histFFT)==0)
      return histFFT;
   return NULL;
  
}
      
int UsefulIcrrStationEvent::fillFFTHistoForClock(int clock_number, TH1D *histFFT) 
{
   TGraph *grFFT = getFFTForClock(clock_number);
   if(!grFFT) return -1;
   Double_t *xVals=grFFT->GetX();
   Double_t *yVals=grFFT->GetY();
   Int_t numPoints=grFFT->GetN();
   for(int i=0;i<numPoints;i++) {    
      histFFT->Fill(xVals[i],yVals[i]);
   }
   delete grFFT;
   return 0;
}

bool UsefulIcrrStationEvent::isCalPulserEvent( )
{

   bool retcode = false;

   IcrrTriggerMonitor *trig = &(this->trig);
   unsigned short msw_clock_counter = trig->rovdd[0];
   unsigned short lsw_clock_counter = trig->rovdd[1];

// printf("UsefulIcrrStationEvent::isCalPulserEvent() - DEBUG - %d %d",msw_clock_counter,lsw_clock_counter);
   if ( ( msw_clock_counter == 0 ) &&                                                                // Rb peak is at msw=0 && lsw=5801+/-5 
        ( TMath::Abs(lsw_clock_counter-5801) < 5 ) ) {
//    printf(" <- is cal pulser\n");
      retcode = true;
   } else {
//    printf(" <- is *not* cal pulser\n");
      retcode = false;
   }

   return retcode;

} // end of UsefulIcrrStationEvent::isCalPulserEvent() member function
*/
