////////////////////////////////////////////////////////////////////////////////////////////////
//class Event:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EVENT_H
#define EVENT_H

#include <vector>
#include "Vector.h"
#include <string>
#include "Primaries.h"
#include "Report.h"
//#include "Neutrino.h"


using namespace std;

class Settings;
class Spectra;
class Primaries;
class Interaction;
//class Report;

class Event {

  public:

      int Event_type;   // 0 : neutrino event,  1 : blackhole event,  2 : monopole event,... etc

      int inu_thrown; // event number. in case we save triggered events only, this event number could be useful

      int inu_passed; // event number. in case we save triggered events only, this event number could be useful

    
      double pnu;   // energy of neutrino
      //Vector nnu;   // direction of neutrino
      string nuflavor;  // flavor of neutrino
      int nuflavorint;  // 1 : nue,  2 : numu,  3 : nutau
      int n_interactions;    // number of interaction inside the ice
      int IsCalpulser;

      int nu_nubar; // 0 : nu, 1 : nu_bar

      
    void Initialize ();

      void Choose_Evt_Type(Settings *settings1);    // choose the event type depend on the settings->EVENT_TYPE value

      std::vector <Interaction> Nu_Interaction;

      //vector <Report> test_report;  

      Event (); // default constructor
      //Event (Settings *settings1, Spectra *spectra1, Primaries *primary1, IceModel *icemodel, Detector *detector, Signal *signal, Secondaries *sec1 );
      Event (Settings *settings1, Spectra *spectra1, Primaries *primary1, IceModel *icemodel, Detector *detector, Signal *signal, Secondaries *sec1, int event_num );

      ~Event(); // destructor


      void delete_all(); // delete all vectors







//--------------------------------------------------
//   Event(int, Spectra*);
//   //  Event(double);
//   ~Event();
//   int iev; // index for this event
//   int part_type; // incident particle type (neutrino, magnetic monopole)
//   vector<double> n_incident;  // unit vector pointing in direction of incident particle trajectory (like a neutrino)
//   double e_incident;  // energy of particle incident at the interaction point
// 
//   // something that describes the radiation that emerges from the shower
//   // this may be a root file or text file that is generated from jaime's program
//   vector<double> eField1m; // E field - 3d vector as a function of frequency at 1 m from the interaction
//   
//   vector<double> getEField1m();
// 
//  
//  Event(double) : pos_sh(),int_type(0),iev(0),part_type(0),n_incident(3), e_incident(1.E18), eField1m(3) {}
//-------------------------------------------------- 

      ClassDef(Event,1);

};

#endif //EVENT_H








/* class Event { */
/* public: */

/* }; */
 
/* class NeutrinoEvent : public Event { */
/* public: */

/* }; */
 
/* class BlackHoleEvent : public Event { */
/* public: */

/* }; */
 
/* class MonopoleEvent : public Event { */
/* public: */

/* }; */
//class NeutrinoEvent;

/* template<class T, class U> */
/*   class NeutrinoEvent : public Event { */
/*  public: */
/*   NeutrinoEvent(T this_ev,U spectra1): Event(this_ev,spectra1) { */

/*   } */
/* }; */








/*
 * Create all available events and print their prices
 */
/* void event_information( EventFactory::EventType eventtype,int this_ev,Spectra spectra1 ) */
/* { */
/*   Event* event = EventFactory::createEvent(eventtype,this_ev,spectra1); */
/* 	//        std::cout << "Price of " << eventtype << " is " << event->getPrice() << std::endl; */
/*         delete event; */
/* } */
 






