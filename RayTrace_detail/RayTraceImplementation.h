// ------------  IMPORTANT!  ------------ //
//  Do not include this header directly!  //
// (It is for internal use by RayTrace.h) //
// -------------------------------------- //

#include <cmath>
#include <stdexcept>
#include <boost/type_traits.hpp>

namespace RayTrace{
	
	minimalRayPosition abs(const minimalRayPosition& p);
	rayPosition abs(const rayPosition& p);
	fullRayPosition abs(const fullRayPosition& p);
	
	///\param theta The angle of incidence (between the incoming ray and the surface normal)
	///\param polarization The angle of polarization, zero for completely perpedicular to the plane of reflection
	///\param n1 The index of refraction of the region the ray is leaving
	///\param n2 The index of refraction of the region the ray is entering
	double fresnelReflect(double theta, double& polarization, double n1, double n2);
	double fresnelTransmit(double theta, double& polarization, double n1, double n2);
	
	template <typename positionType>
	void correctAmplitudeReflect(positionType& pos, double polarization, double n1, double n2){}
	
	template <>
	void correctAmplitudeReflect<fullRayPosition>(fullRayPosition& pos, double polarization, double n1, double n2);
	
	template <>
	void correctAmplitudeReflect<positionRecordingWrapper<fullRayPosition> >(positionRecordingWrapper<fullRayPosition>& pos, double polarization, double n1, double n2);
	
	template <typename positionType>
	void correctAmplitudeTransmit(positionType& pos, double polarization, double n1, double n2){}
	
	template <>
	void correctAmplitudeTransmit<fullRayPosition>(fullRayPosition& pos, double polarization, double n1, double n2);
	
	template <>
	void correctAmplitudeTransmit<positionRecordingWrapper<fullRayPosition> >(positionRecordingWrapper<fullRayPosition>& pos, double polarization, double n1, double n2);
	
	//Callback dispatch machinery
	//object version
	template<typename positionType, typename callbackType>
	void callCallback(callbackType cb, const positionType& p, RKStepType s){
		cb(p,s);
	}
	
	//implementation for function pointers
	template<typename positionType, typename callbackType>
	void callCallback_impl(callbackType* cb, const positionType& p, RKStepType s, const boost::true_type&){
		if(cb!=NULL)
			cb(p,s);
	}
	
	//implementation for object pointers
	template<typename positionType, typename callbackType>
	void callCallback_impl(callbackType* cb, const positionType& p, RKStepType s, const boost::false_type&){
		if(cb!=NULL)
			(*cb)(p,s);
	}
	
	//pointer version
	template<typename positionType, typename callbackType>
	void callCallback(callbackType* cb, const positionType& p, RKStepType s){
		callCallback_impl(cb,p,s,typename boost::is_function<callbackType>());
	}
	//End of callback dispatch machinery


#undef DO_TRACE_WITH_CALLBACK	
#include "DoTraceImplementation.h"
#define DO_TRACE_WITH_CALLBACK 1
#include "DoTraceImplementation.h"	
#undef DO_TRACE_WITH_CALLBACK
	
} //namespace RayTrace
