#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#pragma link C++ class Detector+;
#pragma link C++ class Parameters+;
#pragma link C++ class InstalledStation+;
#pragma link C++ class IdealStation+;
#pragma link C++ class Trigger+;

#pragma link C++ class ARA_station+;
#pragma link C++ class Antenna_string+;
#pragma link C++ class Antenna+;
#pragma link C++ class Surface_antenna+;

#pragma link C++ class Event+;
#pragma link C++ class Report+;
#pragma link C++ class Antenna_r+;
#pragma link C++ class Surface_antenna_r+;
#pragma link C++ class String_r+;
#pragma link C++ class Station_r+;

#pragma link C++ class Position+;
#pragma link C++ class Vector+;
#pragma link C++ class Settings+;
#pragma link C++ class Spectra+;
#pragma link C++ class IceModel+;
#pragma link C++ class EarthModel+;
#pragma link C++ class Y+;
#pragma link C++ class Primaries+;
#pragma link C++ class Interaction+;

#ifdef ARA_UTIL_EXISTS
// #pragma link C++ class AtriEventHkData+;
// #pragma link C++ class AtriSensorHkData+;
// #pragma link C++ class RawAraStationEvent+;
// #pragma link C++ class RawIcrrStationHeader+;
// #pragma link C++ class RawIcrrStationEvent+;
// #pragma link C++ class RawAtriSimpleStationEvent+;
// #pragma link C++ class AraRawIcrrRFChannel+;
// #pragma link C++ class IcrrHkData+;
// #pragma link C++ class FullIcrrHkEvent+;
// #pragma link C++ class IcrrTriggerMonitor+;
// #pragma link C++ class UsefulAraStationEvent+;
// #pragma link C++ class UsefulIcrrStationEvent+;
// #pragma link C++ class UsefulAtriStationEvent+;
// #pragma link C++ class AraAntennaInfo+;
// #pragma link C++ class AraStationInfo+;
// #pragma link C++ class AraGeomTool+;
// #pragma link C++ class AraEventCalibrator+;
// #pragma link C++ class RawAtriStationBlock+;
// #pragma link C++ class RawAtriStationEvent+;
// #pragma link C++ class RawAraGenericHeader+;
// 
// #pragma link C++ typedef AraDataStructureType_t;


#pragma link C++ namespace     AraCalType;
#pragma link C++ enum          AraCalType::EAraCalType;
#pragma link C++ nestedtypedef AraCalType::AraCalType_t;

#pragma link C++ namespace     AraAntType;
#pragma link C++ enum          AraAntType::EAraAntType;
#pragma link C++ nestedtypedef AraAntType::AraAntType_t;


#pragma link C++ namespace     AraAntPol;
#pragma link C++ enum          AraAntPol::EAraAntPol;
#pragma link C++ nestedtypedef AraAntPol::AraAntPol_t;


#pragma link C++ namespace     AraDaqChanType;
#pragma link C++ enum          AraDaqChanType::EAraDaqChanType;
#pragma link C++ nestedtypedef AraDaqChanType::AraDaqChanType_t;


#pragma link C++ namespace     AraLabChip;
#pragma link C++ enum          AraLabChip::EAraLabChip;
#pragma link C++ nestedtypedef AraLabChip::AraLabChip_t;


#pragma link C++ namespace     AraAntDir;
#pragma link C++ enum          AraAntDir::EAraAntDir;
#pragma link C++ nestedtypedef AraAntDir::AraAntDir_t;

#endif

#endif

