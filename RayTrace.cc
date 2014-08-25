#include "RayTrace.h"
#include <cmath>
#include <iterator>
#include <algorithm>
#include <limits>
#include <queue>

#include <boost/scoped_ptr.hpp>

namespace RayTrace{
	
	const double TraceRecord::noReflection(100.0);
	
	//========== minimalRayPosition ==========//
	
	minimalRayPosition::minimalRayPosition():
	x(0.0),z(0.0),theta(0.0){}
	
	minimalRayPosition::minimalRayPosition(double x_, double z_, double theta_):
	x(x_),z(z_),theta(theta_){}
	
	minimalRayPosition& minimalRayPosition::operator +=(const minimalRayPosition& p){
		x+=p.x;
		z+=p.z;
		theta+=p.theta;
		return(*this);
	}
	
	minimalRayPosition& minimalRayPosition::operator -=(const minimalRayPosition& p){
		x-=p.x;
		z-=p.z;
		theta-=p.theta;
		return(*this);
	}
	
	minimalRayPosition& minimalRayPosition::operator *=(double m){
		x*=m;
		z*=m;
		theta*=m;
		return(*this);
	}
	
	const minimalRayPosition minimalRayPosition::operator +(const minimalRayPosition& p) const{
		minimalRayPosition res=*this;
		return(res+=p);
	}
	
	const minimalRayPosition minimalRayPosition::operator -(const minimalRayPosition& p) const{
		minimalRayPosition res=*this;
		return(res-=p);
	}
	
	const minimalRayPosition minimalRayPosition::operator *(double m) const{
		minimalRayPosition res=*this;
		return(res*=m);
	}
	
	const minimalRayPosition operator *(double m, const minimalRayPosition& p){
		return(p*m);
	}
	
	minimalRayPosition abs(const minimalRayPosition& p){
		return(minimalRayPosition(std::abs(p.x), std::abs(p.z), std::abs(p.theta)));
	}
	
	void minimalRayPosition::makeTiny(double s){
		x=s;
		z=s;
		theta=s;
	}
	
	void minimalRayPosition::giveData(TraceRecord& trace) const{
		trace.pathTime=-1.0;
		trace.receiptAngle=theta;
		trace.attenuation=-1.0;
	}
	
	//========== rayPosition ==========//
	
	rayPosition::rayPosition():
	x(0.0),z(0.0),theta(0.0),time(0.0){}
	
	rayPosition::rayPosition(double x_, double z_, double theta_, double time_):
	x(x_),z(z_),theta(theta_),time(time_){}
	
	rayPosition& rayPosition::operator +=(const rayPosition& p){
		x+=p.x;
		z+=p.z;
		theta+=p.theta;
		time+=p.time;
		return(*this);
	}
	
	rayPosition& rayPosition::operator -=(const rayPosition& p){
		x-=p.x;
		z-=p.z;
		theta-=p.theta;
		time-=p.time;
		return(*this);
	}
	
	rayPosition& rayPosition::operator *=(double m){
		x*=m;
		z*=m;
		theta*=m;
		time*=m;
		return(*this);
	}
	
	const rayPosition rayPosition::operator +(const rayPosition& p) const{
		rayPosition res=*this;
		return(res+=p);
	}
	
	const rayPosition rayPosition::operator -(const rayPosition& p) const{
		rayPosition res=*this;
		return(res-=p);
	}
	
	const rayPosition rayPosition::operator *(double m) const{
		rayPosition res=*this;
		return(res*=m);
	}
	
	const rayPosition operator *(double m, const rayPosition& p){
		return(p*m);
	}
	
	rayPosition abs(const rayPosition& p){
		return(rayPosition(std::abs(p.x), std::abs(p.z), std::abs(p.theta), std::abs(p.time)));
	}
	
	void rayPosition::makeTiny(double s){
		x=s;
		z=s;
		theta=s;
		time=s;
	}
	
	void rayPosition::giveData(TraceRecord& trace) const{
		trace.pathTime=time;
		trace.receiptAngle=theta;
		trace.attenuation=-1.0;
	}

	//========== fullRayPosition ==========//
	
	fullRayPosition::fullRayPosition():
	x(0.0),z(0.0),theta(0.0),time(0.0),attenuation(1.0){}
	
	fullRayPosition::fullRayPosition(double x_, double z_, double theta_, double time_, double attenuation_):
	x(x_),z(z_),theta(theta_),time(time_),attenuation(attenuation_){}
	
	fullRayPosition& fullRayPosition::operator +=(const fullRayPosition& p){
		x+=p.x;
		z+=p.z;
		theta+=p.theta;
		time+=p.time;
		attenuation+=p.attenuation;
		return(*this);
	}
	
	fullRayPosition& fullRayPosition::operator -=(const fullRayPosition& p){
		x-=p.x;
		z-=p.z;
		theta-=p.theta;
		time-=p.time;
		attenuation-=p.attenuation;
		return(*this);
	}
	
	fullRayPosition& fullRayPosition::operator *=(double m){
		x*=m;
		z*=m;
		theta*=m;
		time*=m;
		attenuation*=m;
		return(*this);
	}
	
	const fullRayPosition fullRayPosition::operator +(const fullRayPosition& p) const{
		fullRayPosition res=*this;
		return(res+=p);
	}
	
	const fullRayPosition fullRayPosition::operator -(const fullRayPosition& p) const{
		fullRayPosition res=*this;
		return(res-=p);
	}
	
	const fullRayPosition fullRayPosition::operator *(double m) const{
		fullRayPosition res=*this;
		return(res*=m);
	}
	
	const fullRayPosition operator *(double m, const fullRayPosition& p){
		return(p*m);
	}
	
	fullRayPosition abs(const fullRayPosition& p){
		return(fullRayPosition(std::abs(p.x), std::abs(p.z), std::abs(p.theta), std::abs(p.time), std::abs(p.attenuation)));
	}
	
	void fullRayPosition::makeTiny(double s){
		x=s;
		z=s;
		theta=s;
		time=s;
		attenuation=s;
	}
	
	void fullRayPosition::giveData(TraceRecord& trace) const{
		trace.pathTime=time;
		trace.receiptAngle=theta;
		trace.attenuation=attenuation;
	}

	//========== TraceFinder ==========//

	//const double TraceFinder::maximum_ice_depth = -2850.0;
	const double TraceFinder::maximum_ice_depth = -3500.0; // changed for testing max ice depth
	
	///\brief Computes the derivatives of the position coordinates with respect to path length
	///
	///This specialization computes the minimum usable set of coordinates: the radial and vertical 
	///position and the angle. 
	///
	///\param pos The position at which the derivatives are to be computed
	///\param der The object into which to record the calculated derivatives
	///\param frequency The frequency of the signal being propogated, in GHz
	template<>
	void TraceFinder::computeRayDerivatives<minimalRayPosition>(const minimalRayPosition& pos, minimalRayPosition& der, double frequency) const{
		//put dn/dz into der.theta
		double n;
		rModel->indexOfRefractionWithDerivative(pos.z,n,der.theta);
		der.x = sin(pos.theta);
		der.z = cos(pos.theta);
		der.theta = -(der.x * der.theta)/n;
	}
	
	///\brief Computes the derivatives of the position coordinates with respect to path length
	///
	///This specialization computes an intermediate set of coordinates: the radial and 
	///vertical position, the angle, and the time of flight. 
	///
	///\param pos The position at which the derivatives are to be computed
	///\param der The object into which to record the calculated derivatives
	///\param frequency The frequency of the signal being propogated, in GHz
	template<>
	void TraceFinder::computeRayDerivatives<rayPosition>(const rayPosition& pos, rayPosition& der, double frequency) const{
		//put n int der.time, put dn/dz into der.theta
		rModel->indexOfRefractionWithDerivative(pos.z,der.time,der.theta);
		der.x = sin(pos.theta);
		der.z = cos(pos.theta);
		der.theta = -(der.x * der.theta)/der.time;
		der.time /= speedOfLight;
	}
	
	///\brief Computes the derivatives of the position coordinates with respect to path length
	///
	///This specialization computes a complete set of coordinates: the radial and vertical 
	///position, the angle, the time of flight, and the amplitude correction. 
	///
	///\param pos The position at which the derivatives are to be computed
	///\param der The object into which to record the calculated derivatives
	///\param frequency The frequency of the signal being propogated, in GHz
	template<>
	void TraceFinder::computeRayDerivatives<fullRayPosition>(const fullRayPosition& pos, fullRayPosition& der, double frequency) const{
		//put n int der.time, put dn/dz into der.theta
		rModel->indexOfRefractionWithDerivative(pos.z,der.time,der.theta);
		der.x = sin(pos.theta);
		der.z = cos(pos.theta);
		der.theta = -(der.x * der.theta)/der.time;
		der.time /= speedOfLight;
		der.attenuation=-pos.attenuation/aModel->attenuationLength(pos.z,frequency);
		//std::cout << "  eval with z=" << pos.z << " att=" << pos.attenuation << " yields " << der.attenuation << std::endl;
	}

	template<typename positionType>
	void TraceFinder::rkStep(const positionType& pos, const typename positionType::derivativeType& der, const double h, positionType& newPos, positionType& errors, const double frequency) const{
		//const double a[5] = {0.2, 0.3, 0.6, 1.0, 0.875};
		const double b[15] = {0.2, 0.075, 0.3, -11.0/54.0, 1631.0/55296.0, 
			0.225, -0.9, 2.5, 175.0/512.0, 
			1.2, -70.0/27.0, 575.0/13824.0, 
			35.0/27.0, 44275.0/110592.0, 
			253.0/4096.0};
		const double c[6] = {37.0/378.0, 0, 250.0/621.0, 125.0/594.0, 0, 512.0/1771.0};
		const double cp[6] = {c[0]-2825.0/27648.0, 0, c[2]-18575.0/48384.0, c[3]-13525.0/55296.0, c[4]-277.0/14336.0, c[5]-0.25};
		
		positionType temp;
		typename positionType::derivativeType ak2, ak3, ak4, ak5, ak6; //derivatives
		//first step
		temp = pos + h*b[0]*der;
		//second step
		computeRayDerivatives(temp,ak2,frequency);
		temp = pos + h*(b[1]*der + b[5]*ak2);
		//third step
		computeRayDerivatives(temp,ak3,frequency);
		temp = pos + h*(b[2]*der + b[6]*ak2 + b[9]*ak3);
		//fourth step
		computeRayDerivatives(temp,ak4,frequency);
		temp = pos + h*(b[3]*der + b[7]*ak2 + b[10]*ak3 + b[12]*ak4);
		//fifth step
		computeRayDerivatives(temp,ak5,frequency);
		temp = pos + h*(b[4]*der + b[8]*ak2 + b[11]*ak3 + b[13]*ak4 + b[14]*ak5);
		//sixth step
		computeRayDerivatives(temp,ak6,frequency);
		newPos = pos + h*(c[0]*der + c[2]*ak3 + c[3]*ak4 + c[5]*ak6);
		//estimate error
		errors = h*(cp[0]*der + cp[2]*ak3 + cp[3]*ak4 + cp[4]*ak5 + cp[5]*ak6);
	}
	
	template<typename positionType>
	void TraceFinder::rkStep(const positionRecordingWrapper<positionType>& pos, const typename positionRecordingWrapper<positionType>::derivativeType& der, const double h, positionRecordingWrapper<positionType>& newPos, positionRecordingWrapper<positionType>& errors, const double frequency) const{
		//const double a[5] = {0.2, 0.3, 0.6, 1.0, 0.875};
		const double b[15] = {0.2, 0.075, 0.3, -11.0/54.0, 1631.0/55296.0, 
			0.225, -0.9, 2.5, 175.0/512.0, 
			1.2, -70.0/27.0, 575.0/13824.0, 
			35.0/27.0, 44275.0/110592.0, 
			253.0/4096.0};
		const double c[6] = {37.0/378.0, 0, 250.0/621.0, 125.0/594.0, 0, 512.0/1771.0};
		const double cp[6] = {c[0]-2825.0/27648.0, 0, c[2]-18575.0/48384.0, c[3]-13525.0/55296.0, c[4]-277.0/14336.0, c[5]-0.25};
		
		positionRecordingWrapper<positionType> temp;
		typename positionRecordingWrapper<positionType>::derivativeType ak2, ak3, ak4, ak5, ak6; //derivatives
		//first step
		temp = pos + h*b[0]*der;
		temp.recordStep(1);
		//second step
		computeRayDerivatives(temp,ak2,frequency);
		temp = pos + h*(b[1]*der + b[5]*ak2);
		temp.recordStep(2);
		//third step
		computeRayDerivatives(temp,ak3,frequency);
		temp = pos + h*(b[2]*der + b[6]*ak2 + b[9]*ak3);
		temp.recordStep(3);
		//fourth step
		computeRayDerivatives(temp,ak4,frequency);
		temp = pos + h*(b[3]*der + b[7]*ak2 + b[10]*ak3 + b[12]*ak4);
		temp.recordStep(4);
		//fifth step
		computeRayDerivatives(temp,ak5,frequency);
		temp = pos + h*(b[4]*der + b[8]*ak2 + b[11]*ak3 + b[13]*ak4 + b[14]*ak5);
		temp.recordStep(5);
		//sixth step
		computeRayDerivatives(temp,ak6,frequency);
		newPos = pos + h*(c[0]*der + c[2]*ak3 + c[3]*ak4 + c[5]*ak6);
		newPos.takeStepData(temp);
		//estimate error
		errors = h*(cp[0]*der + cp[2]*ak3 + cp[3]*ak4 + cp[4]*ak5 + cp[5]*ak6);
	}
	
	void TraceFinder::replayRkStep(const rkStepRecord& step, const double atten, const double& attenDer, double& newAtten, const double frequency) const{
		const double b[15] = {0.2, 0.075, 0.3, -11.0/54.0, 1631.0/55296.0, 
			0.225, -0.9, 2.5, 175.0/512.0, 
			1.2, -70.0/27.0, 575.0/13824.0, 
			35.0/27.0, 44275.0/110592.0, 
			253.0/4096.0};
		const double c[6] = {37.0/378.0, 0, 250.0/621.0, 125.0/594.0, 0, 512.0/1771.0};
		
		double temp;
		double ak2, ak3, ak4, ak5, ak6; //derivatives
		
		//first step
		temp = atten + step.length*b[0]*attenDer;
		//second step
		ak2=-temp/aModel->attenuationLength(step.z[1],frequency);
		//std::cout << "  eval with z=" << step.z[1] << " att=" << temp << " yields " << ak2 << std::endl;
		temp = atten + step.length*(b[1]*attenDer + b[5]*ak2);
		//third step
		ak3=-temp/aModel->attenuationLength(step.z[2],frequency);
		//std::cout << "  eval with z=" << step.z[2] << " att=" << temp << " yields " << ak3 << std::endl;
		temp = atten + step.length*(b[2]*attenDer + b[6]*ak2 + b[9]*ak3);
		//fourth step
		ak4=-temp/aModel->attenuationLength(step.z[3],frequency);
		//std::cout << "  eval with z=" << step.z[3] << " att=" << temp << " yields " << ak4 << std::endl;
		temp = atten + step.length*(b[3]*attenDer + b[7]*ak2 + b[10]*ak3 + b[12]*ak4);
		//fifth step
		ak5=-temp/aModel->attenuationLength(step.z[4],frequency);
		//std::cout << "  eval with z=" << step.z[4] << " att=" << temp << " yields " << ak5 << std::endl;
		temp = atten + step.length*(b[4]*attenDer + b[8]*ak2 + b[11]*ak3 + b[13]*ak4 + b[14]*ak5);
		//sixth step
		ak6=-temp/aModel->attenuationLength(step.z[5],frequency);
		//std::cout << "  eval with z=" << step.z[5] << " att=" << temp << " yields " << ak6 << std::endl;
		newAtten = atten + step.length*(c[0]*attenDer + c[2]*ak3 + c[3]*ak4 + c[5]*ak6);
	}
	
	///\brief Computes the maximum of a set of coordinate errors, scaled by given factors
	///
	///This specialization computes the minimum usable set of coordinates: the radial and 
	///vertical position and the angle. 
	///
	///\param errors The coordinate errors
	///\param scale The scaling factors for the errors
	///\return The maximum error divided by its associated scale
	template<>
	double TraceFinder::maxError<minimalRayPosition>(const minimalRayPosition& errors, const minimalRayPosition& scale) const{
		double result=std::abs(errors.x/scale.x);
		result=std::max(result,std::abs(errors.z/scale.z));
		result=std::max(result,std::abs(errors.theta/scale.theta));
		return(result);
	}
	
	///\brief Computes the maximum of a set of coordinate errors, scaled by given factors
	///
	///This specialization computes an intermediate set of coordinates: the radial and 
	///vertical position, the angle, and the time of flight. 
	///
	///\param errors The coordinate errors
	///\param scale The scaling factors for the errors
	///\return The maximum error divided by its associated scale
	template<>
	double TraceFinder::maxError<rayPosition>(const rayPosition& errors, const rayPosition& scale) const{
		double result=std::abs(errors.x/scale.x);
		result=std::max(result,std::abs(errors.z/scale.z));
		result=std::max(result,std::abs(errors.theta/scale.theta));
		result=std::max(result,std::abs(errors.time/scale.time));
		return(result);
	}
	
	///\brief Computes the maximum of a set of coordinate errors, scaled by given factors
	///
	///This specialization computes a complete set of coordinates: the radial and vertical 
	///position, the angle, the time of flight, and the amplitude correction. 
	///
	///\param errors The coordinate errors
	///\param scale The scaling factors for the errors
	///\return The maximum error divided by its associated scale
	template<>
	double TraceFinder::maxError<fullRayPosition>(const fullRayPosition& errors, const fullRayPosition& scale) const{
		double result=std::abs(errors.x/scale.x);
		result=std::max(result,std::abs(errors.z/scale.z));
		result=std::max(result,std::abs(errors.theta/scale.theta));
		result=std::max(result,std::abs(errors.time/scale.time));
		result=std::max(result,std::abs(errors.attenuation/scale.attenuation));
		return(result);
	}
	
	template<typename positionType>
	void confirmFinalStepPosition(positionType& pos, const positionType& temp, const double& length){
		pos=temp;
	}
	
	template<typename positionType>
	void confirmFinalStepPosition(positionRecordingWrapper<positionType>& pos, const positionRecordingWrapper<positionType>& temp, const double& length){
		double startZ=pos.z;
		pos=temp;
		pos.setFirstStep(startZ);
		pos.setStepLength(length);
		//pos.dumpData();
	}

	//requires htry>0 !!!
	template<typename positionType>
	void TraceFinder::rkStepControl(double& length, double frequency, positionType& pos, typename positionType::derivativeType& der, const positionType& scale, const double htry, const double eps, double& hdid, double& hnext) const{
		const double SAFETY=0.9;
		//const double growthPower=-0.2;
		//const double shrinkPower=-0.25;
		const double growLimit=1.89e-4; //=(5.0/SAFETY)**(1/growthPower)
		const double shrinkLimit=6561.0; //=(.1/SAFETY)**(1/shrinkPower)
		
		//std::cout << "\tTrying RK(CK) step of size " << htry << std::endl;
		
		double h=htry, errMax;
		positionType errors, temp;
                //std::cerr<<"begin while loop"<<std::endl;
		while(true){
			//take step
			rkStep(pos,der,h,temp,errors,frequency);
			errMax=maxError(errors,scale)/eps;
			if(errMax<=1.0) {
                                //std::cout<<"errMax<=1, break"<<std::endl;
				break; //error is tolerable
                        }
			//need to decrease step size
			if(errMax < shrinkLimit) {
                                //std::cout<<"errMax ("<<errMax<<")<shrinkLimit ("<<shrinkLimit<<"), h = "<<h;
				h=SAFETY*h/sqrt(sqrt(errMax)); // 1/sqrt(sqrt(errMax)) == errMax**shrinkPower
                                //std::cout<<" to "<<h<<std::endl;
                        }
			else {
                                //std::cout<<"errMax ("<<errMax<<")>=shrinkLimit ("<<shrinkLimit<<"), h = "<<h;
				h*=0.1;
                                //std::cout<<" to "<<h<<std::endl;
                        }

                        // test errMax !!
                        if( errMax != errMax ) { // error... nan errMax
                            //std::cout<<"errMax weird!"<<std::endl;
                            break;
                        }

			//std::cout << "\tError too big, reducing step size to " << h << std::endl;
			if((length+h)==length)
				throw std::runtime_error("TraceFinder::rkStepControl: stepsize underflow");

		}
		//increase step size
		if(errMax > growLimit){
			//Here, we use a taylor expansion (about 1) for errMax**growthPower. This 
			//expansion always underestimates the true function (in the domain [0,1])
			//so it will give conservative advice for the next step's size
			errMax-=1.0;
			hnext = SAFETY*h* (1.+errMax*(-.2+errMax*(.12+errMax*(-.088+errMax*.0704))));
		}
		else
			hnext = 5.0*h;
		//std::cout << "\tError acceptable, plan next step size to be " << hnext << std::endl;
		length+=(hdid=h);
		confirmFinalStepPosition(pos,temp,hdid);//pos=temp;
	}
	
	double fresnelReflect(double theta, double& polarization, double n1, double n2){
		double s=sin(theta);
		if(n1*s > n2){ //total internal reflection!
			//std::cout << "Total internal refraction" << std::endl;
			return(1.0);
		}
		double c=cos(theta);
		double d=sqrt(n2*n2 - n1*n1*s*s);
		double rPerpedicular=(n1*c-d)/(n1*c+d);
		double rParallel=(n2*n2*c-n1*d)/(n2*n2*c+n1*d);
		double cp=cos(polarization);
		double sp=sin(polarization);
		double R=sqrt(rPerpedicular*rPerpedicular*cp*cp + rParallel*rParallel*sp*sp);
		polarization=atan((rParallel*sp)/(rPerpedicular*cp));
		return(R);
	}
	
	double fresnelTransmit(double theta, double& polarization, double n1, double n2){
		double s=sin(theta);
		if(n1*s > n2){ //total internal reflection!
			//std::cout << "Total internal refraction" << std::endl;
			return(0.0);
		}
		double c=cos(theta);
		double d=sqrt(n2*n2 - n1*n1*s*s);
		double tPerpedicular=(2.*n1*c)/(n1*c+d);
		double tParallel=(2*n1*n2*c)/(n2*n2*c+n1*d);
		double cp=cos(polarization);
		double sp=sin(polarization);
		double T=sqrt(tPerpedicular*tPerpedicular*cp*cp + tParallel*tParallel*sp*sp);
		polarization=atan((tParallel*sp)/(tPerpedicular*cp));
		return(T);
	}
	
	template <>
	void correctAmplitudeReflect<fullRayPosition>(fullRayPosition& pos, double polarization, double n1, double n2){
		pos.attenuation*=fresnelReflect(pos.theta, polarization, n1, n2);
	}
	
	template <>
	void correctAmplitudeTransmit<fullRayPosition>(fullRayPosition& pos, double polarization, double n1, double n2){
		pos.attenuation*=fresnelTransmit(pos.theta, polarization, n1, n2);
	}
	
	template <>
	void correctAmplitudeReflect<positionRecordingWrapper<fullRayPosition> >(positionRecordingWrapper<fullRayPosition>& pos, double polarization, double n1, double n2){
		pos.attenuation*=fresnelReflect(pos.theta, polarization, n1, n2);
	}
	
	template <>
	void correctAmplitudeTransmit<positionRecordingWrapper<fullRayPosition> >(positionRecordingWrapper<fullRayPosition>& pos, double polarization, double n1, double n2){
		pos.attenuation*=fresnelTransmit(pos.theta, polarization, n1, n2);
	}
	
	double TraceFinder::recalculateAmplitude(const traceReplayRecord& trace, double frequency, double polarization){
		double attenuation=1.0, attenuationDerivative,nextAttenuation;
		//std::cout << " att: " << attenuation << std::endl;
		for(std::vector<stepRecord>::const_iterator step=trace.steps.begin(), end=trace.steps.end(); step!=end; ++step){
			switch(step->stepType){
				case RK_FIRST_STEP:
					//do nothing; should never occur
					break;
				case RK_AIR_STEP:
					//TODO: implement this!
					attenuation*=fresnelTransmit(step->angle, polarization, 1.0, rModel->indexOfRefraction(0.0));
					break;
				case RK_STEP:
					attenuationDerivative=-attenuation/aModel->attenuationLength(step->rkData.z[0],frequency);
					//std::cout << "  eval with z=" << step->rkData.z[0] << " att=" << attenuation << " yields " << attenuationDerivative << std::endl;
					replayRkStep(step->rkData, attenuation, attenuationDerivative, nextAttenuation, frequency);
					attenuation=nextAttenuation;
					break;
				case RK_REFLECT_STEP:
					if(step->angle<(pi/2.)){ //ray was reflected down from the ice surface
						//use the angle _after_ reflection
						attenuation*=fresnelReflect(step->angle, polarization, rModel->indexOfRefraction(maximum_ice_depth),1.0);
					}
					else{ //ray was reflected up from bedrock
						//use the angle _before_ reflection
						attenuation*=fresnelReflect(pi-step->angle, polarization, rModel->indexOfRefraction(0.0),1.0);
					}
					break;
			}
			//std::cout << " att: " << attenuation << std::endl;
		}
		return(attenuation);
	}

	std::pair<bool, double> TraceFinder::traceMax(double emit_depth, const rayTargetRecord& target, int &sol_error, double left, double right ) const{

            // test sol_error
            //int sol_error;

		//std::cout << "traceMax within domain [" << left << ',' << right << ']' << std::endl;
		//unsigned int steps=1;
		std::pair<bool,double> result;
		double a=left; //the left boundary of the search interval
		double b=right; //the right boundary of the search interval
		const double w=(3.0-sqrt(5.0))/2.0;
		double c=a+w*(b-a); //the best estimate point
		//frequency and polarization are unimportant
		//double cy=doTrace<minimalRayPosition>(emit_depth,c,target,SurfaceReflection,0.0,0.0).miss;
		double cy=doTrace<minimalRayPosition>(emit_depth,c,target,SurfaceReflection,0.0,0.0, sol_error).miss;
		//std::cout << "f(" << c << ")=" << cy << std::endl;
		double d,dy;
		const double tol=1e-3;
		while(cy<0.0 && fabs(b-a)>tol*c){
			//std::cout << "\tsearch domain is now [" << a << ',' << b << ']' << std::endl;
			if((c-a)>=(b-c)){ //the left subinterval is larger, so search within it
				d=c-w*(c-a);
				dy=doTrace<minimalRayPosition>(emit_depth,d,target,SurfaceReflection,0.0,0.0, sol_error ).miss; //evaluate the miss distance
				//std::cout << "\t " << d << ' ' << dy << std::endl;
				if(dy>cy){
					b=c;
					c=d;
					cy=dy;
				}
				else
					a=d;
			}
			else{ //otherwise seach within the right subinterval
				d=c+w*(b-c);
				dy=doTrace<minimalRayPosition>(emit_depth,d,target,SurfaceReflection,0.0,0.0, sol_error ).miss; //evaluate the miss distance
				//std::cout << "\t " << d << ' ' << dy << std::endl;
				if(dy>cy){
					a=c;
					c=d;
					cy=dy;
				}
				else
					b=d;
			}
			//steps++;
		}
		//std::cout << "traceMax stopped at c=" << c << " f(c)=" << cy << std::endl;
		result.first = cy>0.0; //whether the best point indicates the existance of roots
		result.second=c;
		//std::cout << "traceMax took " << steps << " steps" << std::endl;
		return(result);
	}

	std::pair<bool, double> TraceFinder::traceMin(double emit_depth, const rayTargetRecord& target, int &sol_error, double left, double right ) const{

            // test sol_error
            //int sol_error;

		//std::cout << "traceMin within domain [" << left << ',' << right << ']' << std::endl;
		//unsigned int steps=1;
		std::pair<bool,double> result;
		double a=left; //the left boundary of the search interval
		double b=right; //the right boundary of the search interval
		const double w=(3.0-sqrt(5.0))/2.0;
		double c=a+w*(b-a); //the best estimate point
		//frequency and polarization are unimportant
		double cy=doTrace<minimalRayPosition>(emit_depth,c,target,BedrockReflection,0.0,0.0, sol_error ).miss;
		//std::cout << "f(" << c << ")=" << cy << std::endl;
		double d,dy;
		const double tol=1e-3;
		while(cy>0.0 && fabs(b-a)>tol*c){
			//std::cout << "\tsearch domain is now [" << a << ',' << b << ']' << std::endl;
			if((c-a)>=(b-c)){ //the left subinterval is larger, so search within it
				//std::cout << "\t left sub-interval\n";
				d=c-w*(c-a);
				dy=doTrace<minimalRayPosition>(emit_depth,d,target,BedrockReflection,0.0,0.0, sol_error ).miss; //evaluate the miss distance
				//std::cout << "\t " << d << ' ' << dy << std::endl;
				if(dy<cy){
					b=c;
					c=d;
					cy=dy;
				}
				else
					a=d;
			}
			else{ //otherwise seach within the right subinterval
				//std::cout << "\t right sub-interval\n";
				d=c+w*(b-c);
				dy=doTrace<minimalRayPosition>(emit_depth,d,target,BedrockReflection,0.0,0.0, sol_error ).miss; //evaluate the miss distance
				//std::cout << "\t " << d << ' ' << dy << std::endl;
				if(dy<cy){
					a=c;
					c=d;
					cy=dy;
				}
				else
					b=d;
			}
			//steps++;
		}
		//std::cout << "traceMin stopped at c=" << c << " f(c)=" << cy << std::endl;
		result.first = cy<0.0; //whether the best point indicates the existance of roots
		result.second=c;
		//std::cout << "traceMin took " << steps << " steps" << std::endl;
		return(result);
	}

	std::pair<double,double> TraceFinder::traceRootImpl(double emit_depth, const rayTargetRecord& target, bool rising, unsigned short allowedReflections, double requiredAccuracy, 
									  double a, TraceRecord& aTrace, double c, TraceRecord& cTrace, double angle, int &sol_error ) const{

            // test sol_error
            //int sol_error;

		const double miss_eps=requiredAccuracy/100.0;
		const double angle_eps=1e-10;
		double lastMiss;
		TraceRecord trace;
		
		double e=0.0,ep=0.0,p,q,r,s,t;
		
		//initialize these so that the miss change check won't fail accidentally
		trace.miss=1.0e3*miss_eps;
		lastMiss=-1.0e3*miss_eps;
		
		while((c-a)>angle_eps && fabs(trace.miss-lastMiss)>miss_eps){
                    // changed
			//std::cout << "\tAngular range=" << (c-a) << ", miss change=" << fabs(trace.miss-lastMiss) << std::endl;
			lastMiss=trace.miss;
			//std::cout << " Angular range is now [" << a << ',' << c << "]\n";
			//std::cout << " Trying angle=" << angle << std::endl;
			trace=doTrace<minimalRayPosition>(emit_depth,angle,target,allowedReflections,0.0,0.0, sol_error );
			//steps++;
			//std::cout << "  miss distance was " << trace.miss << '\n';
			if(std::abs(trace.miss) < requiredAccuracy){
				//std::cout << "traceRoot took " << steps << " steps" << std::endl;
				//doTrace<fullRayPosition>(emit_depth,angle,target,allowedReflections,frequency,polarization)
				return(std::make_pair(angle,trace.miss));
			}
			
			//calculate the quadratic interpolation
			r=trace.miss/cTrace.miss;
			s=trace.miss/aTrace.miss;
			t=aTrace.miss/cTrace.miss;
			p=s*(t*(r-t)*(c-angle)+(r-1)*(angle-a));
			q=(r-1)*(s-1)*(t-1);
			
			//collapse the interval
			if((trace.miss>0.0) != rising){
				a=angle;
				aTrace=trace;
			}
			else{
				c=angle;
				cTrace=trace;
			}
			
			//accept the interpolation only if it falls within the current boundaries
			if((angle+(p/q))>=a && (angle+(p/q))<=c){
				//would like to interpolate, but only do so if recent convergence has been suitably rapid
				if(std::abs(p/q) >= 0.5*std::abs(e)){ //it has not; bisect instead
					//std::cout << " will bisect" << std::endl;
					e=ep;
					ep=0.5*(a+c)-angle;
					angle=0.5*(a+c);
				}
				else{ //it has; use the interpolation
					//std::cout << " will interpolate" << std::endl;
					angle+=p/q;
					e=ep;
					ep=p/q;
				}
			}
			else{ //otherwise, bisect
				//std::cout << " will bisect" << std::endl;
				e=ep;
				ep=0.5*(a+c)-angle;
				angle=0.5*(a+c);
			}
		}
		//std::cout << "traceRoot took " << steps << " steps and bailed out" << std::endl;
		return(std::make_pair(angle,trace.miss));
	}
	
	//rough implementation of Brent's method
	//This implementation neglects some of the checks made by Brent to ensure 
	//rapidity of convergence, and does not bother to use secant interpolation. 
	//In practice it seems not to matter, as it is still capable of the same 
	//best-case behavior as the full Brent's method. 
	//This implementation is more conveniently controllable for this purpose 
	//than NR's zbrent, with regard to the stopping conditions
	std::pair<double,double> TraceFinder::traceRoot(double emit_depth, const rayTargetRecord& target, double minAngle, double maxAngle, bool rising, unsigned short allowedReflections, double requiredAccuracy, int &sol_error ) const{

            // test sol_error
            //int sol_error;


		double a=minAngle,c=maxAngle;
		TraceRecord aTrace, cTrace;
		//frequency and polarization are irrelevant on all traces except the final one
		aTrace = doTrace<minimalRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0, sol_error );
		cTrace = doTrace<minimalRayPosition>(emit_depth,c,target,allowedReflections,0.0,0.0, sol_error );
		double angle=0.5*(minAngle+maxAngle);
		//std::cout << "Attempting to find root in range [" << minAngle << ',' << maxAngle << ']' << std::endl;
		//std::cout << " end point miss values are " << aTrace.miss << ',' << cTrace.miss << std::endl;
		if((aTrace.miss>0) == (cTrace.miss>0)) //didn't bracket a root. Return a nonsense angle and a huge error
			return(std::make_pair(-1.,1e6));
		return(traceRootImpl(emit_depth, target, rising, allowedReflections, requiredAccuracy, a, aTrace, c, cTrace, angle, sol_error ));
	}
	
        /*
	std::pair<double,double> TraceFinder::refineRoot(double emit_depth, const rayTargetRecord& target, const TraceRecord& seed, bool rising, unsigned short allowedReflections, double requiredAccuracy, int &sol_error ) const{

            // test sol_error
            //int sol_error;


		//bracket the root
		//std::cout << "Attempting to refine root near theta=" << seed.launchAngle << std::endl;
		double a, c;
		TraceRecord aTrace, cTrace;
		double testDisp=.001; //TODO: this is a problem since there could be a second root arbitrarily close to the one we seek
		//If we start with testDisp larger than the distance betwee these roots we can bracket them both and incorrectly conclude that we've bracketed nothing
		const unsigned int maxTests=(unsigned int)std::ceil(0.5*(sqrt(1.0+(8*pi/testDisp))-1.0));
		if((seed.miss>0.0)==rising)
			testDisp*=-1;
		unsigned int i;
		a=seed.launchAngle;
		for(i=0; i<maxTests; i++){
			a+=testDisp;
			testDisp*=2.0;
			if(a<0.0)
				a=0.0;
			else if(a>pi)
				a=pi;
			//std::cout << "\ttesting theta=" << a << std::endl;
			aTrace=doTrace<minimalRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0, sol_error );
                        // changed
			std::cout << "\t (miss=" << aTrace.miss << ")" << std::endl;
			if(std::abs(aTrace.miss) < requiredAccuracy) {
                                // changed
                                std::cout << "miss<requiredAccuracy" << std::endl;
				return(std::make_pair(a,aTrace.miss));
                        }
				//return(doTrace<fullRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0));
			if((aTrace.miss*seed.miss)<0.0) {
                                // changed
                                std::cout << "trace.miss * seed.miss < 0, break" << std::endl;
				break;
                        }
		}
		if(i==maxTests){
                        // changed
                        std::cout << "refineRoot, reached maxTests" << std::endl;

			//std::cerr << "Last attempt at bracketing was: [" << seed.launchAngle << "->" << seed.miss << ',' << a << "->" << aTrace.miss << ']' << std::endl;
			//std::cerr << "\trising edge = " << std::boolalpha << rising << std::endl;
			throw std::runtime_error("TraceFinder::refineRoot: exceeded maximum allowed number of steps for bracketing.");
		}
		if((seed.miss>0.0)==rising){
			c=seed.launchAngle;
			cTrace=seed;
		}
		else{
			c=a;
			cTrace=aTrace;
			a=seed.launchAngle;
			aTrace=seed;
		}
		
		//std::cout << "Endpoint miss distances are " << aTrace.miss << " and " << cTrace.miss << std::endl;
		double angle=0.5*(a+c);
		return(traceRootImpl(emit_depth, target, rising, allowedReflections, requiredAccuracy, a, aTrace, c, cTrace, angle, sol_error ));
	}
        */



        /*
        // changed
        // test obtaining min value case
	std::pair<double,double> TraceFinder::refineRoot(double emit_depth, const rayTargetRecord& target, const TraceRecord& seed, bool rising, unsigned short allowedReflections, double requiredAccuracy, int &sol_error ) const{

            // test sol_error
            //int sol_error;

            double min_angle, min_angle_plus, min_angle_minus;
            double min_miss = 1e10; // start with big value


		//bracket the root
		//std::cout << "Attempting to refine root near theta=" << seed.launchAngle << std::endl;
		double a, c;
		TraceRecord aTrace, cTrace;
		double testDisp=.001; //TODO: this is a problem since there could be a second root arbitrarily close to the one we seek
		//If we start with testDisp larger than the distance betwee these roots we can bracket them both and incorrectly conclude that we've bracketed nothing
		const unsigned int maxTests=(unsigned int)std::ceil(0.5*(sqrt(1.0+(8*pi/testDisp))-1.0));
		if((seed.miss>0.0)==rising)
			testDisp*=-1;
		unsigned int i;
		a=seed.launchAngle;
		for(i=0; i<maxTests; i++){
			a+=testDisp;
			testDisp*=2.0;
			if(a<0.0)
				a=0.0;
			else if(a>pi)
				a=pi;
			//std::cout << "\ttesting theta=" << a << std::endl;
			aTrace=doTrace<minimalRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0, sol_error );

                        if ( min_miss > std::abs(aTrace.miss) ) {
                            min_miss = std::abs(aTrace.miss);
                            min_angle = a;
		
                            if((seed.miss>0.0)==rising) {
                                min_angle_plus = a - testDisp;
                                min_angle_minus = a + testDisp;
                            }
                            else {
                                min_angle_plus = a + testDisp;
                                min_angle_minus = a - testDisp;
                            }
                        }

                        // changed
			std::cout << "\t (miss=" << aTrace.miss << ") at angle: "<<a<< std::endl;
			if(std::abs(aTrace.miss) < requiredAccuracy) {
                                // changed
                                std::cout << "miss<requiredAccuracy" << std::endl;
				return(std::make_pair(a,aTrace.miss));
                        }
				//return(doTrace<fullRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0));
			if((aTrace.miss*seed.miss)<0.0) {
                                // changed
                                std::cout << "trace.miss * seed.miss < 0, break" << std::endl;
				break;
                        }
		}
		if(i==maxTests){
                        // changed
                        std::cout << "refineRoot, reached maxTests" << std::endl;

			//std::cerr << "Last attempt at bracketing was: [" << seed.launchAngle << "->" << seed.miss << ',' << a << "->" << aTrace.miss << ']' << std::endl;
			//std::cerr << "\trising edge = " << std::boolalpha << rising << std::endl;
			throw std::runtime_error("TraceFinder::refineRoot: exceeded maximum allowed number of steps for bracketing.");
		}
		if((seed.miss>0.0)==rising){
			c=seed.launchAngle;
			cTrace=seed;
		}
		else{
			c=a;
			cTrace=aTrace;
			a=seed.launchAngle;
			aTrace=seed;
		}
		
                // changed
                std::cout<<"min miss : "<<min_miss<<" at angle : "<<min_angle<<", test w/ min angle : "<<min_angle_minus<<" max angle : "<<min_angle_plus<<std::endl;

		//std::cout << "Endpoint miss distances are " << aTrace.miss << " and " << cTrace.miss << std::endl;
		double angle=0.5*(a+c);
                // changed
		//return(traceRootImpl(emit_depth, target, rising, allowedReflections, requiredAccuracy, a, aTrace, c, cTrace, angle, sol_error ));
		return(evenmore_refineRoot(emit_depth, target, seed, min_angle_minus, min_angle_plus, rising, allowedReflections, requiredAccuracy, sol_error ));
	}
        */

        // changed
        // test obtaining min value case more carefully
	std::pair<double,double> TraceFinder::evenmore_refineRoot(double emit_depth, const rayTargetRecord& target, const TraceRecord& seed, double angle_min, double angle_max, bool rising, unsigned short allowedReflections, double requiredAccuracy, int &sol_error ) const{

            // test sol_error
            //int sol_error;

            double min_angle;
            double min_miss = 1e10; // start with big value


		//bracket the root
		//std::cout << "Attempting to refine root near theta=" << seed.launchAngle << std::endl;
		double a, c;
		TraceRecord aTrace, cTrace;
		//double testDisp=.001; //TODO: this is a problem since there could be a second root arbitrarily close to the one we seek
		double testDisp=.00001; //TODO: this is a problem since there could be a second root arbitrarily close to the one we seek
		//If we start with testDisp larger than the distance betwee these roots we can bracket them both and incorrectly conclude that we've bracketed nothing
		//const unsigned int maxTests=(unsigned int)std::ceil(0.5*(sqrt(1.0+(8*pi/testDisp))-1.0));
		const unsigned int maxTests=(unsigned int)((angle_max - angle_min)/testDisp);

		//if((seed.miss>0.0)==rising)
			//testDisp*=-1;
		unsigned int i;
		//a=seed.launchAngle;
		a=angle_min;
		for(i=0; i<maxTests; i++){
			a+=testDisp;
			testDisp*=2.0;
			if(a<0.0)
				a=0.0;
			else if(a>pi)
				a=pi;
			//std::cout << "\ttesting theta=" << a << std::endl;
			aTrace=doTrace<minimalRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0, sol_error );

                        if ( min_miss > std::abs(aTrace.miss) ) {
                            min_miss = std::abs(aTrace.miss);
                            min_angle = a;
                        }

                        // changed
			//std::cout << "\t second (miss=" << aTrace.miss << ") at angle: "<<a<< std::endl;
			if(std::abs(aTrace.miss) < requiredAccuracy) {
                                // changed
                                //std::cout << "miss<requiredAccuracy" << std::endl;
				return(std::make_pair(a,aTrace.miss));
                        }
				//return(doTrace<fullRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0));
			if((aTrace.miss*seed.miss)<0.0) {
                                // changed
                                //std::cout << "second trace.miss * seed.miss < 0, break" << std::endl;
				break;
                        }
		}





		if(i==maxTests){
                        // changed
                        //std::cout << "refineRoot, reached maxTests" << std::endl;

			//std::cerr << "Last attempt at bracketing was: [" << seed.launchAngle << "->" << seed.miss << ',' << a << "->" << aTrace.miss << ']' << std::endl;
			//std::cerr << "\trising edge = " << std::boolalpha << rising << std::endl;
			throw std::runtime_error("TraceFinder::refineRoot: exceeded maximum allowed number of steps for bracketing.");
		}
		if((seed.miss>0.0)==rising){
			c=seed.launchAngle;
			cTrace=seed;
		}
		else{
			c=a;
			cTrace=aTrace;
			a=seed.launchAngle;
			aTrace=seed;
		}
		
                // changed
                //std::cout<<"min miss : "<<min_miss<<" at angle : "<<min_angle<<std::endl;

		//std::cout << "Endpoint miss distances are " << aTrace.miss << " and " << cTrace.miss << std::endl;
		double angle=0.5*(a+c);
                // changed
		return(traceRootImpl(emit_depth, target, rising, allowedReflections, requiredAccuracy, a, aTrace, c, cTrace, angle, sol_error ));
	}


        // changed
        // test obtaining min value case
        // loosen required accuracy and see if they pass
	std::pair<double,double> TraceFinder::refineRoot(double emit_depth, const rayTargetRecord& target, const TraceRecord& seed, bool rising, unsigned short allowedReflections, double requiredAccuracy, int &sol_error ) const{

            // test sol_error
            //int sol_error;

            double min_angle, min_angle_plus, min_angle_minus;
            double min_miss = 1e10; // start with big value


		//bracket the root
		//std::cout << "Attempting to refine root near theta=" << seed.launchAngle << std::endl;
		double a, c;
		TraceRecord aTrace, cTrace, minmissTrace;
		double testDisp=.001; //TODO: this is a problem since there could be a second root arbitrarily close to the one we seek
		//If we start with testDisp larger than the distance betwee these roots we can bracket them both and incorrectly conclude that we've bracketed nothing
		const unsigned int maxTests=(unsigned int)std::ceil(0.5*(sqrt(1.0+(8*pi/testDisp))-1.0));
		if((seed.miss>0.0)==rising)
			testDisp*=-1;
		unsigned int i;
		a=seed.launchAngle;
		for(i=0; i<maxTests; i++){
			a+=testDisp;
			testDisp*=2.0;
			if(a<0.0)
				a=0.0;
			else if(a>pi)
				a=pi;
			//std::cout << "\ttesting theta=" << a << std::endl;
			aTrace=doTrace<minimalRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0, sol_error );

                        if ( min_miss > std::abs(aTrace.miss) ) {
                            min_miss = std::abs(aTrace.miss);
                            min_angle = a;
                            minmissTrace = aTrace;
		
                            if((seed.miss>0.0)==rising) {
                                min_angle_plus = a - testDisp;
                                min_angle_minus = a + testDisp;
                            }
                            else {
                                min_angle_plus = a + testDisp;
                                min_angle_minus = a - testDisp;
                            }
                        }

                        // changed
			//std::cout << "\t (miss=" << aTrace.miss << ") at angle: "<<a<< std::endl;
			if(std::abs(aTrace.miss) < requiredAccuracy) {
                                // changed
                                //std::cout << "miss<requiredAccuracy" << std::endl;
				return(std::make_pair(a,aTrace.miss));
                        }
				//return(doTrace<fullRayPosition>(emit_depth,a,target,allowedReflections,0.0,0.0));
			if((aTrace.miss*seed.miss)<0.0) {
                                // changed
                                //std::cout << "trace.miss * seed.miss < 0, break" << std::endl;
				break;
                        }
		}
		if(i==maxTests){
                        // changed
                        //std::cout << "refineRoot, reached maxTests" << std::endl;

			//std::cerr << "Last attempt at bracketing was: [" << seed.launchAngle << "->" << seed.miss << ',' << a << "->" << aTrace.miss << ']' << std::endl;
			//std::cerr << "\trising edge = " << std::boolalpha << rising << std::endl;
			throw std::runtime_error("TraceFinder::refineRoot: exceeded maximum allowed number of steps for bracketing.");
		}
		if((seed.miss>0.0)==rising){
			c=seed.launchAngle;
			cTrace=seed;
		}
		else{
			c=a;
			cTrace=aTrace;
			a=seed.launchAngle;
			aTrace=seed;
		}
		
                // changed
                //std::cout<<"min miss : "<<min_miss<<" at angle : "<<min_angle<<", test w/ min angle : "<<min_angle_minus<<" max angle : "<<min_angle_plus<<std::endl;


                /*
                // changed
                // test loosen required accuray 
                double loosen_cut = 5.; // const value for now
                if ( min_miss < loosen_cut ) {

                    std::cout<<"min miss : "<<min_miss<<" at angle : "<<min_angle<<", passed loosen cut "<<loosen_cut<<std::endl;
                    return(std::make_pair(min_angle,min_miss));
                }
                */



		//std::cout << "Endpoint miss distances are " << aTrace.miss << " and " << cTrace.miss << std::endl;
		double angle=0.5*(a+c);
                // changed
		//return(traceRootImpl(emit_depth, target, rising, allowedReflections, requiredAccuracy, a, aTrace, c, cTrace, angle, sol_error ));
		return(evenmore_refineRoot(emit_depth, target, seed, min_angle_minus, min_angle_plus, rising, allowedReflections, requiredAccuracy, sol_error ));
	}













	
	bool shorterPath(const TraceRecord& a, const TraceRecord& b){
		return(a.pathLen < b.pathLen);
	}
	
	///A quick-and-dirty function to simultaneously sort two sequences according to the values in the first sequence. 
	///That is, the first sequence is sorted in the normal way, and the second sequence undergoes the same permutation. 
	///This is only expected to be used on very short sequences, so it uses the inefficient but simple insertion sort algorithm. 
	///
	///\param begin1 An iterator to the beginning of the first sequence
	///\param end1 An iterator to the end of the first sequence
	///\param begin2 An iterator to the beginning of the second sequence 
	////		(which must be the same length as the first sequence)
	///\param comp The comparison function object used to determine the 
	///		proper ordering of the elements in the first sequence
	template<typename RandomAccessIterator1, typename RandomAccessIterator2, class Compare>
	void dual_insertion_sort(RandomAccessIterator1 begin1, RandomAccessIterator1 end1, RandomAccessIterator2 begin2, Compare comp){
		RandomAccessIterator2 frontier2=begin2+1;
		for(RandomAccessIterator1 frontier1=begin1+1; frontier1!=end1; ++frontier1,++frontier2){
			RandomAccessIterator1 it1=frontier1;
			RandomAccessIterator2 it2=frontier2;
			while(it1!=begin1){
				RandomAccessIterator1 nit=it1-1;
				if(comp(*it1,*nit)){
					std::iter_swap(it1,nit);
					std::iter_swap(it2,it2-1);
					it1=nit;
					--it2;
					continue;
				}
				break;
			}
		}
	}






















	//std::vector<TraceRecord> TraceFinder::findPaths(Vector sourcePos, Vector targetPos, double frequency, double polarization, unsigned short allowedReflections, double requiredAccuracy, std::vector<traceReplayRecord>* replayBuffer) const{
	std::vector<TraceRecord> TraceFinder::findPaths(Vector sourcePos, Vector targetPos, double frequency, double polarization, int &sol_cnt, int &sol_error, unsigned short allowedReflections, double requiredAccuracy, std::vector<traceReplayRecord>* replayBuffer) const{
	//std::vector<TraceRecord> TraceFinder::findPaths(Vector sourcePos, Vector targetPos, double frequency, double polarization, int &sol_cnt, int &sol_error, int mode, unsigned short allowedReflections, double requiredAccuracy, std::vector<traceReplayRecord>* replayBuffer) const{

                sol_cnt = 0;
                sol_error = 0;
                

                // changed
                //std::cout<<"begin findPaths"<<std::endl;
                /*
                std::cout<<"findPaths mode : "<<mode<<std::endl;
                if (mode == 1) {    // mode 1 for src, trg exchanged
                    Vector Pos_tmp;
                    Pos_tmp = sourcePos;
                    sourcePos = targetPos;
                    targetPos = Pos_tmp;
                }
                */




		const double initialSafety=1.e-2;
		const double minimumSafety=1.e-5;
		
		std::vector<TraceRecord> results;
		boost::scoped_ptr<pathRecorder<fullRayPosition> > recorder;
		if(replayBuffer!=NULL)
			recorder.reset(new pathRecorder<fullRayPosition>());
		//std::cout << "Finding paths from (" << sourcePos.GetX() << ',' << sourcePos.GetY() << ',' << sourcePos.GetZ() << ") to (" << targetPos.GetX() << ',' << targetPos.GetY() << ',' << targetPos.GetZ() << ')' << std::endl;
		if(sourcePos.GetZ()<maximum_ice_depth || targetPos.GetZ()<maximum_ice_depth)
                        {
			return(results);
                        }
		bool fullyContained=(sourcePos.GetZ()<=0.0 && targetPos.GetZ()<=0.0); //whether the ray is entirely inside the ice
		double dist = sqrt((targetPos.GetX()-sourcePos.GetX())*(targetPos.GetX()-sourcePos.GetX())+(targetPos.GetY()-sourcePos.GetY())*(targetPos.GetY()-sourcePos.GetY()));
		rayTargetRecord target(targetPos.GetZ(),dist);
		bool dualDirect=false; //whether there are knwon to be two direct solution rays
		
		//special case for pesky near-vertical rays
		if(dist<=requiredAccuracy){
			
                    // changed
                    //std::cout << "Computing direct ray (vertical)" << std::endl;

			if(replayBuffer==NULL)
                                {
				results.push_back(doVerticalTrace<fullRayPosition>(sourcePos.GetZ(), (sourcePos.GetZ()<targetPos.GetZ()?0.0:pi), target, NoReflection, frequency, polarization));
                                sol_cnt++;
                                //std::cout<<"\n\tSOL_CNT ADDED! 1\n"<<std::endl;
                                }
			else{
				results.push_back(doVerticalTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(), (sourcePos.GetZ()<targetPos.GetZ()?0.0:pi), target, NoReflection, frequency, polarization, recorder.get()));
				replayBuffer->push_back(recorder->getData());
				recorder->clearData();
                                sol_cnt++;
                                //std::cout<<"\n\tSOL_CNT ADDED! 2\n"<<std::endl;
			}
			if((allowedReflections & SurfaceReflection) && fullyContained){
				//std::cout << "Computing surface-reflected ray" << std::endl;
				if(replayBuffer==NULL)
                                        {
					results.push_back(doVerticalTrace<fullRayPosition>(sourcePos.GetZ(), 0.0, target, SurfaceReflection, frequency, polarization));
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 3\n"<<std::endl;
                                        }
				else{
					results.push_back(doVerticalTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(), 0.0, target, SurfaceReflection, frequency, polarization, recorder.get()));
					replayBuffer->push_back(recorder->getData());
					recorder->clearData();
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 4\n"<<std::endl;
				}
			}
			if(allowedReflections & BedrockReflection){
				//std::cout << "Computing bedrock-reflected ray" << std::endl;
				if(replayBuffer==NULL)
                                        {
					results.push_back(doVerticalTrace<fullRayPosition>(sourcePos.GetZ(), pi, target, BedrockReflection, frequency, polarization));
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 5\n"<<std::endl;
                                        }
				else{
					results.push_back(doVerticalTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(), pi, target, BedrockReflection, frequency, polarization, recorder.get()));
					replayBuffer->push_back(recorder->getData());
					recorder->clearData();
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 6\n"<<std::endl;
				}
			}
			std::sort(results.begin(),results.end(),&shorterPath);
			return(results);
		} // if special case for pesky near-vertical rays
		

                // if either src or trg not in the ice
		if(!fullyContained){
			
                    // changed
                    //std::cout << "Ray not fully contained" << std::endl;

			if(replayBuffer==NULL)
                                {
				results.push_back(findUncontainedFast(sourcePos, targetPos, frequency, polarization, requiredAccuracy, sol_error ));
                                sol_cnt++;
                                //std::cout<<"\n\tSOL_CNT ADDED! 7\n"<<std::endl;
                                }
			else{
				results.push_back(findUncontainedFast(sourcePos, targetPos, frequency, polarization, requiredAccuracy, sol_error, recorder.get()));
				replayBuffer->push_back(recorder->getData());
                                sol_cnt++;
                                //std::cout<<"\n\tSOL_CNT ADDED! 8\n"<<std::endl;
			}
			return(results);
		} // if either src or trg not in the ice


                // changed 
                // both cout
		//std::cout << "Looking for estimate" << std::endl;
		indexOfRefractionModel::RayEstimate est=rModel->estimateRayAngle(sourcePos.GetZ(), targetPos.GetZ(), dist); //TODO: put this back!!!
		//std::cout << "Estimate result : " << est.status << " refractionmodel solution : "<< indexOfRefractionModel::SOLUTION << std::endl;


		if(est.status==indexOfRefractionModel::SOLUTION){ //got a solution estimate
			
                    // changed
                    //std::cout << "Got fast solution, est.angle: " << est.angle << std::endl;

			//std::cout << "\t(theta=" << est.angle << ')' << std::endl;
			if(replayBuffer==NULL)
                                {
                                //std::cout << "replayBuffer==NULL" <<std::endl;
				results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(), est.angle, target, NoReflection, frequency, polarization, sol_error )); //this is supposed to be a direct ray, so no reflections should be needed
                                        
                                // changed always (if)?
                                //std::cout << "done 1st doTrace (if)\n" << std::endl;
                                //std::cout << "results pathLen: "<<results.front().pathLen<<" miss: "<<results.front().miss<<"\n" << std::endl;
                                sol_cnt++;
                                //std::cout<<"\n\tSOL_CNT ADDED! 9\n"<<std::endl;
                                }
			else{
                                //std::cout << "replayBuffer!=NULL" <<std::endl;
				results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(), est.angle, target, NoReflection, frequency, polarization, sol_error, recorder.get())); //this is supposed to be a direct ray, so no reflections should be needed

                                // changed
                                //std::cout << "done 1st doTrace (else)\n" << std::endl;
				replayBuffer->push_back(recorder->getData());
				recorder->clearData();
                                sol_cnt++;
                                //std::cout<<"\n\tSOL_CNT ADDED! 10\n"<<std::endl;
			}
			if(std::abs(results.front().miss) > requiredAccuracy){
                            // changed
                            //std::cout << "Fast Solution not close enough (" << results.front().miss << "), pathLen: "<<results.front().pathLen<<"; refining" << std::endl;

				est.angle=refineRoot(sourcePos.GetZ(), target, results.front(), false, NoReflection, requiredAccuracy, sol_error ).first;
                            // changed
                            //std::cout << "done refining" << std::endl;

				if(replayBuffer==NULL)
                                        {
                                        //std::cout << "Begin 2nd doTrace" << std::endl;
					results.front()=doTrace<fullRayPosition>(sourcePos.GetZ(), est.angle, target, NoReflection, frequency, polarization, sol_error );
                                
                                        // changed
                                        //std::cout << "done 2nd doTrace (if)\n" << std::endl;
                                        //std::cout << "results pathLen: "<<results.front().pathLen<<" miss: "<<results.front().miss<<"\n" << std::endl;

                                        // changed
                                        // if still don't satisfy accuracy, remove
                                        if ( std::abs(results.front().miss) > requiredAccuracy ) {
                                            //std::cout << "missed too much, pop back!" << std::endl;
                                            results.pop_back();
                                            sol_cnt--;
                                        }
                                        else if ( std::isnan(results.front().miss) ) {
                                            //std::cout << "nan info! pop back!" << std::endl;
                                            results.pop_back();
                                            sol_cnt--;
                                        }
                                        //else
                                            //sol_cnt++; // sol_cnt already done

                                        //std::cout<<"\n\tSOL_CNT ADDED! 11\n"<<std::endl;
                                        }
				else{
                                        //std::cout << "Begin 2nd doTrace" << std::endl;
					results.front()=doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(), est.angle, target, NoReflection, frequency, polarization, sol_error, recorder.get());
                                        
                                        // changed
                                        //std::cout << "done 2nd doTrace (else)\n" << std::endl;
                                        //sol_cnt++;

                                        // changed
                                        // if still don't satisfy accuracy, remove
                                        if ( std::abs(results.front().miss) > requiredAccuracy ) {
                                            //std::cout << "missed too much, pop back!" << std::endl;
                                            results.pop_back();
                                            recorder->clearData();
                                            replayBuffer->pop_back();
                                            sol_cnt--;
                                        }
                                        else if ( std::isnan(results.front().miss) ) {
                                            //std::cout << "nan info! pop back!" << std::endl;
                                            results.pop_back();
                                            recorder->clearData();
                                            replayBuffer->pop_back();
                                            sol_cnt--;
                                        }
                                        else {
                                            replayBuffer->front()=recorder->getData();
                                            recorder->clearData();
                                            //sol_cnt++; // sol_cnt already done
                                        }
                                        //std::cout<<"\n\tSOL_CNT ADDED!  12\n"<<std::endl;
				}
			}
			//if(allowedReflections & SurfaceReflection){
			if((allowedReflections & SurfaceReflection) && (sol_cnt>0) ){ // try reflect solution if there was 1st solution existing
				
                            // changed
                            //std::cout << "Looking for surface reflected solution, prev est.angle: " << est.angle << std::endl;

				//double angle=traceRoot(sourcePos.GetZ(),target,0.0,est.angle-1e-1,true,SurfaceReflection,requiredAccuracy);
				
				//Here, we want to find another root, potentially near the direct solution, but distinct from it. 
				//We expect this root to be in the damain [0,est.angle], but we don't want the root finding algorithm
				//to get distracted by the root that we know is at est.angle, so we subtract a small angle from the upper 
				//bound. However, the roots may be arbitrarily close together, so any guess for this safety value may be too large
				//and we may end up exlusing the part of the domain which contains the root we want. So, if the root finding
				//doesn't manage to find a suitable root we decrease the safety value. 
				double safety=initialSafety;
				while(est.angle-safety < 0)
					safety/=10.;
				std::pair<double,double> root(-1.,10.*requiredAccuracy); //first element is the launch angle, second is miss distance
				do{
					root=traceRoot(sourcePos.GetZ(),target,
								   0.0, //minimum angle
								   est.angle-safety, //maximum angle
								   true,SurfaceReflection,requiredAccuracy, sol_error );
					safety/=10.;
				}while(safety>=minimumSafety && std::abs(root.second)>requiredAccuracy);
				//TODO: check here that root.first is within [0.0,est.angle-safety]
				if(replayBuffer==NULL)
                                        {
                                        // I guess here occurs one kind of error (log_ray_solver), (not occur when NO surface reflection)
                                        // error which shows
                                        // errMax>shrinkLimit, h = inf
                                        //
                                        // (infinite loop shows above)
                                        //
                                        //std::cout << "Begin 3rd doTrace (replayBuffer==NULL)" << std::endl;
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),root.first,target,allowedReflections,frequency,polarization, sol_error ));
                                        // changed
                                        //std::cout << "Done 3rd doTrace (if)\n" << std::endl;
                                        //std::cout << "results pathLen: "<<results[1].pathLen<<" miss: "<<results[1].miss<<"\n" << std::endl;
                                        //sol_cnt++;

                                        // changed
                                        // if still don't satisfy accuracy, remove
                                        if ( std::abs(results[1].miss) > requiredAccuracy ) {
                                            //std::cout << "missed too much, pop back!" << std::endl;
                                            results.pop_back();
                                        }
                                        else if ( std::isnan(results[1].miss) ) {
                                            //std::cout << "nan info! pop back!" << std::endl;
                                            results.pop_back();
                                        }
                                        else
                                            sol_cnt++;

                                        //std::cout<<"\n\tSOL_CNT ADDED! 13\n"<<std::endl;
                                        }
				else{
                                        //std::cout << "Begin 3rd doTrace (replayBuffer!=NULL)" << std::endl;
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),root.first,target,allowedReflections,frequency,polarization,sol_error, recorder.get()));
                                        // changed
                                        //std::cout << "Done 3rd doTrace (else)\n" << std::endl;
					//replayBuffer->push_back(recorder->getData());
					//recorder->clearData();
                                        //sol_cnt++;

                                        // changed
                                        // if still don't satisfy accuracy, remove
                                        if ( std::abs(results[1].miss) > requiredAccuracy ) {
                                            //std::cout << "missed too much, pop back!" << std::endl;
                                            results.pop_back();
                                            recorder->clearData();
                                            replayBuffer->pop_back();
                                        }
                                        else if ( std::isnan(results[1].miss) ) {
                                            //std::cout << "nan info! pop back!" << std::endl;
                                            results.pop_back();
                                            recorder->clearData();
                                            replayBuffer->pop_back();
                                        }
                                        else {
                                            replayBuffer->front()=recorder->getData();
                                            recorder->clearData();
                                            sol_cnt++;
                                        }

                                        //std::cout<<"\n\tSOL_CNT ADDED! 14\n"<<std::endl;
				}
			} // if allowedreflections & SurfaceReflection
                        
                        // changed
                        // bedrock reflect part not modified yet!
			if(allowedReflections & BedrockReflection){
				
                            // changed
                            //std::cout << "Looking for bedrock reflected solution" << std::endl;

				//double angle=traceRoot(sourcePos.GetZ(),target,est.angle+1e-4,pi,true,BedrockReflection,requiredAccuracy).first;
				
				//this is the same logic as for the surface reflected case, but the seacrh is in the domain [est.angle,pi]
				double safety=initialSafety;
				while(est.angle+safety > pi)
					safety/=10;
				std::pair<double,double> root(-1.,10.*requiredAccuracy); //first element is the launch angle, second is miss distance
				do{
					root=traceRoot(sourcePos.GetZ(),target,
								   est.angle+safety, //minimum angle
								   pi, //maximum angle
								   true,BedrockReflection,requiredAccuracy, sol_error );
					safety/=10.;
				}while(safety>=minimumSafety && std::abs(root.second)>requiredAccuracy);
				if(replayBuffer==NULL)
                                        {
                                        //std::cout << "Begin 4th doTrace (replayBuffer==NULL)" << std::endl;
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),root.first,target,allowedReflections,frequency,polarization, sol_error ));
                                        //std::cout << "Done 4th doTrace" << std::endl;
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 15\n"<<std::endl;
                                        }
				else{
                                        //std::cout << "Begin 4th doTrace (replayBuffer!=NULL)" << std::endl;
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),root.first,target,allowedReflections,frequency,polarization,sol_error, recorder.get()));
                                        //std::cout << "Done 4th doTrace" << std::endl;
					replayBuffer->push_back(recorder->getData());
					recorder->clearData();
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 16\n"<<std::endl;
				}
			}
		} //got a solution estimate

		else{ //did not get a solution estimate
			
                    
                    // changed
                    //std::cout << "Didn't get fast solution" << std::endl;

			double a=0.0,b=pi; //the endpoints of a range bracketing a root
			std::pair<bool, double> minRes, maxRes, altMinRes;
			switch(est.status){
				case indexOfRefractionModel::SOLUTION:
					//impossible, this case was handled above
					break;
				case indexOfRefractionModel::UPPER_LIMIT:
					//std::cout << "UPPER_LIMIT, Got only upper limit of " << est.angle << std::endl;
					b=est.angle;//std::min(pi,est.angle+.1); //add a small amount in case the solution is right at the boundary given by the estimate
					minRes = traceMin(sourcePos.GetZ(),target, sol_error, b,pi );
					maxRes = traceMax(sourcePos.GetZ(),target, sol_error, 0.0,b );
					break;
				case indexOfRefractionModel::LOWER_LIMIT:
					//std::cout << "LOWER_LIMIT, Got only lower limit of " << est.angle << std::endl;
					a=est.angle;//std::max(0.0,est.angle-.1);
					minRes = traceMin(sourcePos.GetZ(),target, sol_error, a,pi );
					maxRes = traceMax(sourcePos.GetZ(),target, sol_error, 0.0,a );
					break;
				case indexOfRefractionModel::NO_SOLUTION:
					//std::cout << "NO_SOLUTION, Got no solution" << std::endl;
					//Note that getting indexOfRefractionModel::NO_SOLUTION implies no _direct_ solution exists
					//if we were instructed to find reflected rays, we still need to do the work to look for them
					if((allowedReflections & SurfaceReflection) || (allowedReflections & BedrockReflection)){
						maxRes = traceMax(sourcePos.GetZ(),target, sol_error, 0.0,pi );
						minRes = traceMin(sourcePos.GetZ(),target, sol_error, 0.0,pi );
					}
					break;
				case indexOfRefractionModel::UNKNOWN:
                                        // I guess here occurs one kind of error (log_ray_solver2), (occur whether or not surface reflection)
                                        // error which shows
                                        // terminate called throwing an exceptionAbort trap: 6
                                        //
					//std::cout << "UNKNOWN, Got no information" << std::endl;
					maxRes = traceMax(sourcePos.GetZ(),target, sol_error, 0.0,pi );
					//std::cout << "maxRes : " << maxRes.first << std::endl;
					minRes = traceMin(sourcePos.GetZ(),target, sol_error, 0.0,pi );
					//std::cout << "minRes : " << minRes.first << std::endl;
					break;
			}
                        //std::cout << "done est.status switch" << std::endl;

			if(minRes.first && maxRes.first){
                            //std::cout<<" minRes & maxRes true" << std::endl;
				if(maxRes.second>minRes.second){
					//watch out; there must be two direct solutions
					dualDirect=true;
					altMinRes=minRes;
					minRes = traceMin(sourcePos.GetZ(),target, sol_error, maxRes.second,pi ); //look for the other minimum which should exist
				}
				else{ //we need to waste some time checking for another minimum
					altMinRes = traceMin(sourcePos.GetZ(),target, sol_error, 0.0,maxRes.second );
					dualDirect=altMinRes.first; //if a minimum is found, which is negative, then there is second direct solution
				}
			}
			
			if(minRes.first && (allowedReflections & BedrockReflection)){
				//std::cout << "Looking for bedrock reflected solution" << std::endl;
				double angle=traceRoot(sourcePos.GetZ(),target,minRes.second,pi,true,BedrockReflection,requiredAccuracy, sol_error ).first;
				if(replayBuffer==NULL)
                                        {
                                        //std::cout << "Begin 1st doTrace (replayBuffer==NULL)" << std::endl;
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization, sol_error ));
                                        //std::cout << "Done 1st doTrace" << std::endl;
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 17\n"<<std::endl;
                                        }
				else{
                                        //std::cout << "Begin 1st doTrace (replayBuffer!=NULL)" << std::endl;
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization,sol_error, recorder.get()));
                                        //std::cout << "Done 1st doTrace" << std::endl;
					replayBuffer->push_back(recorder->getData());
					recorder->clearData();
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 18\n"<<std::endl;
				}
			}
			if(maxRes.first && (allowedReflections & SurfaceReflection) && fullyContained && !dualDirect){ //can't reflect frpom the surface if passing through the ice or if two direct solutions exist
				//std::cout << "Looking for surface reflected solution" << std::endl;
				double angle=traceRoot(sourcePos.GetZ(),target,0.0,maxRes.second,true,SurfaceReflection,requiredAccuracy, sol_error ).first;
				if(replayBuffer==NULL)
                                        {
                                        //std::cout << "Begin 2nd doTrace (replayBuffer==NULL)" << std::endl;
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization, sol_error ));
                                        //std::cout << "Done 2nd doTrace" << std::endl;
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 19\n"<<std::endl;
                                        }
				else{
                                        //std::cout << "Begin 2nd doTrace (replayBuffer!=NULL)" << std::endl;
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization,sol_error, recorder.get()));
                                        //std::cout << "Done 2nd doTrace" << std::endl;
					replayBuffer->push_back(recorder->getData());
					recorder->clearData();
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 20\n"<<std::endl;
				}
			}
			if((minRes.first || maxRes.first) && est.status!=indexOfRefractionModel::NO_SOLUTION){
				if(minRes.first && minRes.second<b)
					b=minRes.second;
				if(maxRes.first && maxRes.second>a)
					a=maxRes.second;
				//std::cout << "Looking for direct solution" << std::endl;
				double angle=traceRoot(sourcePos.GetZ(),target,a,b,false,NoReflection,requiredAccuracy, sol_error ).first;
				//std::cout << "Angle for direct ray is " << angle << std::endl;
				if(replayBuffer==NULL)
                                        {
                                        //std::cout << "Begin 3rd doTrace (replayBuffer==NULL)" << std::endl;
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization, sol_error ));
                                        //std::cout << "Done 3rd doTrace" << std::endl;
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 21\n"<<std::endl;
                                        }
				else{
                                        //std::cout << "Begin 3rd doTrace (replayBuffer!=NULL)" << std::endl;
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),angle,target,allowedReflections,frequency,polarization,sol_error, recorder.get()));
                                        //std::cout << "Done 3rd doTrace" << std::endl;
					replayBuffer->push_back(recorder->getData());
					recorder->clearData();
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 22\n"<<std::endl;
				}
				//TODO: should this next part be kept?
				if(std::abs(results.back().miss) > 10.0*requiredAccuracy){ //if it wasn't really a solution, throw it away
					//std::cout << "Supposed solution missed by " << results.back().miss << " meters, rejecting" << std::endl;
					results.pop_back();
					if(replayBuffer!=NULL)
						replayBuffer->pop_back();
				}
			}
			if(dualDirect){
				//std::cout << "Looking for secondary direct ray" << std::endl;
				double angle=traceRoot(sourcePos.GetZ(),target,altMinRes.second,maxRes.second,true,NoReflection,requiredAccuracy, sol_error ).first;
				//std::cout << "Angle for second direct ray is " << angle << std::endl;
				if(replayBuffer==NULL)
                                        {
                                        //std::cout << "Begin 4th doTrace (replayBuffer==NULL)" << std::endl;
					results.push_back(doTrace<fullRayPosition>(sourcePos.GetZ(),angle,target,NoReflection,frequency,polarization, sol_error ));
                                        //std::cout << "Done 4th doTrace" << std::endl;
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 23\n"<<std::endl;
                                        }
				else{
                                        //std::cout << "Begin 4th doTrace (replayBuffer!=NULL)" << std::endl;
					results.push_back(doTrace<positionRecordingWrapper<fullRayPosition> >(sourcePos.GetZ(),angle,target,NoReflection,frequency,polarization, sol_error, recorder.get()));
                                        //std::cout << "Done 4th doTrace" << std::endl;
					replayBuffer->push_back(recorder->getData());
					recorder->clearData();
                                        sol_cnt++;
                                        //std::cout<<"\n\tSOL_CNT ADDED! 24\n"<<std::endl;
				}
			}
		} //did not get a solution estimate
		//std::cout << "Have " << results.size() << " solution" << (results.size()!=1?"s":"") << std::endl;
                //

		if(results.size() > 1){
			if(replayBuffer==NULL)
                                {
                                //std::cout << "replayBuffer==NULL" <<std::endl;
				std::sort(results.begin(),results.end(),&shorterPath);
                                }
			else
                                {
                                //std::cout << "replayBuffer!=NULL" <<std::endl;
				dual_insertion_sort(results.begin(),results.end(),replayBuffer->begin(),&shorterPath);
                                }
		}
                else {  // no solution case
                    sol_cnt = 0;
                    //std::cout<<"\n\tSOL_CNT 0 as no sol\n"<<std::endl;
                }
		return(results);
	}



















	
	TraceRecord TraceFinder::findUncontainedFast(Vector sourcePos, Vector targetPos, double frequency, double polarization, double requiredAccuracy, int &sol_error, pathRecorder<fullRayPosition>* recorder ) const{

            // test sol_error
            //int sol_error;


		TraceRecord trace;
		//TODO: if both points are above surface, give obvious answers
		if(sourcePos.GetZ()<maximum_ice_depth || targetPos.GetZ()<maximum_ice_depth)
			return(trace);
		//the index of refraction of the ice at the surface
		const double surfaceN = rModel->indexOfRefraction(0.0);
		double wholeDist = sqrt((targetPos.GetX()-sourcePos.GetX())*(targetPos.GetX()-sourcePos.GetX())+(targetPos.GetY()-sourcePos.GetY())*(targetPos.GetY()-sourcePos.GetY()));
		double upper=std::max(sourcePos.GetZ(),targetPos.GetZ());
		double lower=std::min(sourcePos.GetZ(),targetPos.GetZ());
		rayTargetRecord target(lower,wholeDist);
		double minDist=0.0, maxDist=wholeDist;
		double dist=(minDist+maxDist)/2.; //the distance which will be covered by the ray while in the ice
		double thetaA=0.0,thetaB=0.0;
		while((maxDist-minDist)>requiredAccuracy/10.){
			//std::cout << "Trying dist=" << dist << std::endl;
			//thetaA is the angle of the straight line ray above the ice
			thetaA=pi-atan((wholeDist-dist)/upper);
			//thetaB is the angle the straight line ray would have below the ice
			thetaB=pi-asin(sin(pi-thetaA)/surfaceN);
			//std::cout << "Angles are " << thetaA << ' ' << thetaB << std::endl;
			
			indexOfRefractionModel::RayEstimate est=rModel->estimateRayAngle(0.0, lower, dist);
			if(est.status == indexOfRefractionModel::SOLUTION){
				//std::cout << "Got estimate at angle " << est.angle << std::endl;
				if(est.angle<thetaB)
					maxDist=dist;
				else if(est.angle>thetaB)
					minDist=dist;
				else //we got lucky
					break;
			}
			else{
				target.distance = dist;
				trace = doTrace<minimalRayPosition>(0.0,thetaB,target,NoReflection,0.0,0.0, sol_error );
				//std::cout << "Trace missed by " << trace.miss << std::endl;
				if(std::abs(trace.miss) <= requiredAccuracy)
					break;
				if(trace.miss<0.0)
					maxDist=dist;
				else if(trace.miss>0.0)
					minDist=dist;
			}
			
			dist=(minDist+maxDist)/2.;
		}
		
		//compute the path length the above ice segment
		double extraLength=sqrt(upper*upper+(wholeDist-dist)*(wholeDist-dist));
		//compute attenuation due to passing through the ice surface
		double transmissionAttenuation;
		if(sourcePos.GetZ() > targetPos.GetZ())
			transmissionAttenuation=fresnelTransmit(pi-thetaA, polarization, 1.0, surfaceN);
		else
			transmissionAttenuation=fresnelTransmit(pi-thetaB, polarization, surfaceN, 1.0);
		//do the final version trace
		target.distance = dist;	
		if(recorder==NULL)
			trace = doTrace<fullRayPosition>(0.0,thetaB,target,NoReflection,frequency,polarization, sol_error );
		else{
			callCallback(recorder,positionRecordingWrapper<fullRayPosition>(fullRayPosition(0.0,upper,((sourcePos.GetZ() > targetPos.GetZ())?thetaA:thetaB)-pi,extraLength/speedOfLight,transmissionAttenuation)),RK_AIR_STEP);
			trace = doTrace<positionRecordingWrapper<fullRayPosition> >(0.0,thetaB,target,NoReflection,frequency,polarization, sol_error, recorder);
		}
		//include surface crossing attenuation
		trace.attenuation*=transmissionAttenuation;
		//add on path length the above ice segment
		trace.pathLen+=extraLength;
		trace.pathTime+=extraLength/speedOfLight;
		//fix up the angles
		if(sourcePos.GetZ() > targetPos.GetZ()){ //upper is the source depth
			trace.launchAngle=thetaA;
			//receipt angle is as calculated by doTrace
			trace.reflectionAngle=-thetaA;
		}
		else{ //lower is the source depth
			trace.launchAngle=pi-trace.receiptAngle;
			trace.receiptAngle=pi-thetaA;
			trace.reflectionAngle=thetaB-pi;
		}
		return(trace);
	}
	
	double TraceFinder::signalStrength(const TraceRecord& ray, const Vector& src, const Vector& trg, unsigned short allowedReflections, int &sol_error ) const{

            // test sol_error
            //int sol_error;


		const double changeThresh=1e-3; //10^(-n) should give about n digits of accuracy
		const double startDelta=.002;
		//std::cout << "Considering focusing from " << src << " to " << trg << " at angle " << ray.launchAngle << std::endl;
		double dist = sqrt((trg.GetX()-src.GetX())*(trg.GetX()-src.GetX())+(trg.GetY()-src.GetY())*(trg.GetY()-src.GetY()));
		if(ray.launchAngle<.01 || ray.launchAngle>(pi-.01)){
			//std::cout << "ray too close to vertical" << std::endl;
			return(-1.0);
		}
		
		//std::cout << "Original ray miss distance: " << ray.miss << std::endl;
		unsigned int steps=1;
		double scaledChange=1e10;
		double lastChange=0.0;
		unsigned int scaleIncrease=0;
		const double multInc=2.0;
		double ratio,accel;
		std::vector<double> terms;
		std::queue<double> past;
		for(double deltaTheta=startDelta; scaledChange>changeThresh; deltaTheta/=multInc){
			//std::cout << "\tdeltaTheta = " << deltaTheta << std::endl;
			//std::cout << "\tlaunch angle will be " << ray.launchAngle-deltaTheta << std::endl;
			TraceRecord testRay = doTrace<minimalRayPosition>(src.GetZ(), ray.launchAngle+deltaTheta, rayTargetRecord(trg.GetZ(),dist), allowedReflections, 0.0, 0.0, sol_error );
			//std::cout << "\t\ttest ray miss distance: " << testRay.miss << std::endl;
			ratio=std::abs(deltaTheta/(ray.miss-testRay.miss));
			//std::cout << "\t\tdeltaTheta/deltaZ =  " << ratio << std::endl;
			terms.push_back(ratio);
			if(!past.empty()){
				double old=past.back();
				double mult=multInc;
				accel=ratio;
				past.push(ratio);
				for(unsigned int i=0; i<steps-1; i++){
					mult*=multInc;
					accel=(mult*accel-past.front())/(mult-1.0);
					past.pop();
					past.push(accel);
				}
				//std::cout << "\t\textrapolated value = " << accel << std::endl;
				scaledChange=std::abs(old-accel)/accel;
				//std::cout << "\t\tscaled change = " << scaledChange << std::endl;
				if(scaledChange>lastChange){
					//std::cout << "\t\tWARNING: SCALE INCREASE!" << std::endl;
					scaleIncrease++;
					if(scaleIncrease>=3){
						//std::cout << "\t\tWARNING: scale increased over " << scaleIncrease <<  " consecutive steps" << std::endl;
						//scaleIncrease>=n means terms.size()>=n+1
						//accel=(8**(terms.end()-2)-6**(terms.end()-1)+terms.back())/3.0; //n=2 case
						accel=(64.**(terms.end()-3)-56.**(terms.end()-2)+14.**(terms.end()-1)-terms.back())/21.0; //n=3 case
						break;
					}
					//break;
				}
				else
					scaleIncrease=0;
				lastChange=scaledChange;
			}
			else{
				past.push(ratio);
				accel=ratio;
			}
			steps++;
		}
		
		double ir=rModel->indexOfRefraction(trg.GetZ())/rModel->indexOfRefraction(src.GetZ());
		return(sqrt(ir*accel/dist));
	}
	
	//explicitly instantiate
	template void TraceFinder::rkStepControl<minimalRayPosition>(double&, double, minimalRayPosition&, minimalRayPosition&, const minimalRayPosition&, const double, const double, double&, double&) const;

} //namespace RayTrace
