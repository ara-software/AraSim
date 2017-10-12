#ifndef RAYTRACE_H
#define RAYTRACE_H

#include <ostream>
#include <vector>
#ifndef __CINT__
#include <boost/shared_ptr.hpp>
#endif
#include "Vector.h"

namespace RayTrace{
	
	///Pi (unitless)
	const double pi = 3.1415926535897932384626433832795028841971;
	
	///The speed of light in vacuum (meters per second).
	const double speedOfLight = 299792458;

	///\brief A structure to describe the location of a target for a ray-trace
	struct rayTargetRecord{
		///The vertical coordinate of the target, in meters
		double depth;
		
		///The horizontal distance form the source to the target, in meters
		double distance;
		
		///Constructs a rayTargetRecord
		///\param dep The depth of the target below the ice surface, in meters
		///\param dist The horizontal distance from the source to the target, in meters
		rayTargetRecord(double dep, double dist):depth(dep),distance(dist){}
	};
	
	///A constant used to indicate that no reflections should be allowed when performing a ray-trace
	const unsigned short NoReflection = 0;
	///A constant used to indicate that a ray-trace should be allowed to reflect off of the ice/air interface
	const unsigned short SurfaceReflection = 1;
	///A constant used to indicate that a ray-trace should be allowed to reflect off of the ice/bedrock interface
	const unsigned short BedrockReflection = 2;
	///A constant used to indicate that a ray-trace should be allowed to reflect off of either of the ice boundaries
	const unsigned short AllReflections = SurfaceReflection | BedrockReflection;

	///\brief A record of the results of a ray-trace. 
	///
	///Certain fields may not be computed by some functions which return these objects. 
	struct TraceRecord{
		///The length traversed along the path, in meters
		double pathLen;
		
		///\brief The time for light to traverse the path, in seconds
		///
		///This value may not be computed in some calculations, in this case it 
		///will be left with the default value of -1.0. 
		double pathTime;
		
		///The angle at which the path left the source, in radians
		double launchAngle;
		
		///The angle at which the path reached the target
		double receiptAngle;
		
		///\brief The angle, if any, at which the ray was reflected (or refracted), in radians
		///
		///In the case of no reflection, this value will be equal to the constant noReflection (100.0). 
		///For actual reflections, this value will be the angle of the ray when it struck 
		///the surface. For refraction, this value will be the negative of the angle of the ray when 
		///it struck the surface
		double reflectionAngle;
		
		///The distance by which the ray missed the target vertically, in meters
		double miss;
		
		///\brief The factor by which the amplitude of a signal traversing the path is attenuated
		///
		///This includes the effects of the attenuation model and fresnel effects if any reflections occur, 
		///but does _not_ include geometric spreading of the signal (1/r effect). 
		///This value may not be computed in some calculations, in this case it 
		///will be left with the default value of -1.0. 
		double attenuation;
		
		///\brief The angle of the ray's polarization
		///
		///This anlge is defined so that a value of zero corresponds to polarization in the direction 
		///perpedicular to the plane of the ray's propagation
		double polarization;



                // Eugene added below parameters for debug
                int sol_no;
                int sol_error;

		
		///\brief Constructs a TraceRecord with default values
		///
		///All values are set to 0.0 except pathTime, reflectionAngle, and amplitude. 
		//TraceRecord():pathLen(0.0),pathTime(-1.0),launchAngle(0.0),receiptAngle(0.0),reflectionAngle(noReflection),miss(0.0),attenuation(-1.0),polarization(0.0){}
		TraceRecord():pathLen(0.0),pathTime(-1.0),launchAngle(0.0),receiptAngle(0.0),reflectionAngle(noReflection),miss(0.0),attenuation(-1.0),polarization(0.0),sol_no(0),sol_error(0){}
		
		///\brief Constructs a TraceRecord with the given parameters
		///
		///\param len The length traversed along the path, in meters
		///\param time The time for light to traverse the path, in seconds
		///\param launch The angle at which the path left the source, in radians
		///\param receipt The angle at which the path reached the target
		///\param refl The angle, if any, at which the ray was reflected (or refracted), in radians
		///\param mis The distance by which the ray missed the target vertically, in meters
		///\param att The factor by which the amplitude of a signal traversing the path is decreased
		///\param pol The polarization angle of the ray, in radians
		TraceRecord(double len, double time, double launch, double receipt, double refl, double mis, double att, double pol):
		pathLen(len),pathTime(time),launchAngle(launch),receiptAngle(receipt),reflectionAngle(refl),miss(mis),attenuation(att),polarization(pol){}
		
		static const double noReflection;



	};

	///\brief A class providing a description of the index of refraction of the ice. 
	///
	///Subclasses must implement the three main functions to provide the index of refraction 
	///as a function of depth, and the derivative with respect to depth. They may also override 
	///estimateRayAngle if they have some knowledge that allows them to help the ray-tracer 
	///guess the correct launch angle for a path joining two points. 
	class indexOfRefractionModel{
	public:
		virtual ~indexOfRefractionModel(){};
		
		///Computes the index of refraction as a function of depth. 
		///
		///\param z The vertical coordinate of the point at which the index of refraction is to be calculated. 
		///\return The index of refraction at the given position
		virtual double indexOfRefraction(double z) const=0;
		
		///Computes the derivative of index of refraction with respect to depth  as a function of depth. 
		///
		///\param z The vertical coordinate of the point at which the index of refraction derivative is to be calculated. 
		///\return The index of refraction derivative at the given position
		virtual double indexOfRefractionDerivative(double z) const=0;
		
		///Computes the index of refraction as a function of depth and its derivative with respect to 
		///depth at the same time. 
		///
		///This is expected to be more efficient, since sub-calculations are often shared. 
		///\param z The vertical coordinate of the point at which the index of refraction is to be calculated. 
		///\param n The variable into which to place the calculated index of refraction. 
		///\param dndz The variable into which to place the calculated derivative. 
		virtual void indexOfRefractionWithDerivative(double z, double& n, double& dndz) const=0;
		
		///Used to describe the status of an estimate resulting from estimateRayAngle. 
		enum EstimateStatus{
			///The result is an exact solution
			SOLUTION,
			///There is no solution
			NO_SOLUTION,
			///There is no solution above the result value
			UPPER_LIMIT,
			///There is no solution below the result value
			LOWER_LIMIT,
			///Nothing can be determined about the existance of a solution
			UNKNOWN
		};
		struct RayEstimate{
			double angle;
			EstimateStatus status;
			RayEstimate():angle(0.0),status(UNKNOWN){}
			RayEstimate(EstimateStatus stat, double a):angle(a),status(stat){}
		};
		
		///Uses knowledge of the behavior of rays in the particular ice model to make estimates of the 
		///launch angles of rays to connect source and target positions
		///
		///Subclasses may choose not to implement this function, leaving in place the default 
		///implementation, which provides no advice. 
		///
		///\param sourceDepth The vertical coordinate of the source position, in meters. 
		///\param receiverDepth The vertical coordinate of the target position, in meters. 
		///\param distance The horizontal distance between the source and target positions, in meters. 
		///\return A pair containing an EstimateStatus value and an angle whose interpretation depends on the status. 
		virtual RayEstimate estimateRayAngle(double sourceDepth, double receiverDepth, double distance) const{
			return(RayEstimate(UNKNOWN,0.0));
		}
	};
	
	///\brief A class providing a description of the attenuation length for radio signals in the ice. 
	class attenuationModel{
	public:
		virtual ~attenuationModel(){};
		
		///Computes the attenuation length as a function of depth. 
		///
		///\param z The vertical coordinate of the point at which the attenuation length is to be calculated. 
		///\param frequency The signal frequency, in GHz. 
		///\return The attenuation length at the given depth, in meters. 
		virtual double attenuationLength(double z, double frequency) const=0;
	};

	///\brief A structure for tracking intermediate coordiantes during a ray-trace. 
	///
	///This version keeps track of only vertical and radial position coordiates, and the angle of the ray. 
	struct minimalRayPosition{
		///The radial coordinate of the ray position, in meters
		double x;
		
		///The vertical coordinate of the ray position, in meters
		double z;
		
		///The angle of the ray, in radians (from the upward vertical)
		double theta;
		
		///Default constructs a minimalRayPosition with all coordinates initialized to zero. 
		minimalRayPosition();
		
		///Constructs a minimalRayPosition with the given values. 
		///
		///\param x_ The radial coordinate of the ray position, in meters
		///\param z_ The vertical coordinate of the ray position, in meters
		///\param theta_ The angle of the ray, in radians (from the upward vertical)
		minimalRayPosition(double x_, double z_, double theta_);
		
		///Overloaded increment operator
		///\param p The position whose coordinates to add to this posttion's coordinates
		minimalRayPosition& operator +=(const minimalRayPosition& p);
		
		///Overloaded decrement operator
		///\param p The position whose coordinates to subtract this position's coordinates
		minimalRayPosition& operator -=(const minimalRayPosition& p);
		
		///Overloaded scaling operator
		///\param m The factor by which to scale this position's coordinates
		minimalRayPosition& operator *=(double m);
		
		///Overloaded addition operator
		///\param p The position whose coordinates to add to this position's coordinates
		const minimalRayPosition operator +(const minimalRayPosition& p) const;
		
		///Overloaded subtraction operator
		///\param p The position whose coordinates to subtract from this position's coordinates
		const minimalRayPosition operator -(const minimalRayPosition& p) const;
		
		///Overloaded multiplication operator
		///\param m The factor by which to scale this position's coordinates
		const minimalRayPosition operator *(double m) const;
		
		///\brief Replaces all of the coordinate values with the given value. 
		///\param s The value to use for all coordinates. 
		void makeTiny(double s);
		
		///\brief Copies information from the position object to a TraceRecord object
		///\param trace The trace record object into which to place data from this position. 
		void giveData(TraceRecord& trace) const;
		
		///\brief Adds flight time to the position. 
		///This function does nothing, as this position object does not record time. 
		///\param t The time, in seconds, to add. 
		void addTime(double t){}
		
		typedef minimalRayPosition derivativeType;
	};
	///Overloaded mutiplication operator for scaling a minimalRayPosition. 
	///\param m The factor by which to scale the coordinates of the position
	///\param p The position to scale. 
	const minimalRayPosition operator *(double m, const minimalRayPosition& p);
	
	///\brief A structure for tracking intermediate coordiantes during a ray-trace. 
	///
	///This version keeps track of vertical and radial position coordiates, the angle of the ray, 
	///and the time to traverse the path. 
	struct rayPosition{
		///The radial coordinate of the ray position, in meters
		double x;
		
		///The vertical coordinate of the ray position, in meters
		double z;
		
		///The angle of the ray, in radians (from the upward vertical)
		double theta;
		
		///The time to traverse the path, in seconds
		double time;
		
		///Default constructs a rayPosition with all coordinates initialized to zero. 
		rayPosition();
		
		///Constructs a rayPosition with the given values. 
		///
		///\param x_ The radial coordinate of the ray position, in meters
		///\param z_ The vertical coordinate of the ray position, in meters
		///\param theta_ The angle of the ray, in radians (from the upward vertical)
		///\param time_ The time to traverse the path, in seconds
		rayPosition(double x_, double z_, double theta_, double time_);
		
		///Overloaded increment operator
		///\param p The position whose coordinates to add to this position's coordinates
		rayPosition& operator +=(const rayPosition& p);
		
		///Overloaded decrement operator
		///\param p The position whose coordinates to subtract this position's coordinates
		rayPosition& operator -=(const rayPosition& p);
		
		///Overloaded scaling operator
		///\param m The factor by which to scale this position's coordinates
		rayPosition& operator *=(double m);
		
		///Overloaded addition operator
		///\param p The position whose coordinates to add to this position's coordinates
		const rayPosition operator +(const rayPosition& p) const;
		
		///Overloaded subtraction operator
		///\param p The position whose coordinates to subtract from this position's coordinates
		const rayPosition operator -(const rayPosition& p) const;
		
		///Overloaded multiplication operator
		///\param m The factor by which to scale this position's coordinates
		const rayPosition operator *(double m) const;
		
		///\brief Replaces all of the coordinate values with the given value. 
		///\param s The value to use for all coordinates. 
		void makeTiny(double s);
		
		///\brief Copies information from the position object to a TraceRecord object
		///\param trace The trace record object into which to place data from this position. 
		void giveData(TraceRecord& trace) const;
		
		///\brief Adds flight time to the position. 
		///\param t The time, in seconds, to add. 
		void addTime(double t){ time+=t; }
		
		typedef rayPosition derivativeType;
	};
	///Overloaded mutiplication operator for scaling a rayPosition. 
	///\param m The factor by which to scale the coordinates of the position
	///\param p The position to scale. 
	const rayPosition operator *(double m, const rayPosition& p);
	
	///\brief A structure for tracking intermediate coordiantes during a ray-trace. 
	///
	///This version keeps track of vertical and radial position coordiates, the angle of the ray, 
	///the time to traverse the path, and the factor by which signal aplitude is decreases along the path. 
	struct fullRayPosition{
		///The radial coordinate of the ray position, in meters
		double x;
		
		///The vertical coordinate of the ray position, in meters
		double z;
		
		///The angle of the ray, in radians (from the upward vertical)
		double theta;
		
		///The time to traverse the path, in seconds
		double time;
		
		///The factor by which signal amplitude is decreased along the path
		double attenuation;
		
		///Default constructs a fullRayPosition with all coordinates initialized to zero, except the amplitude, 
		///which is initiallized to 1.0. 
		fullRayPosition();
		
		///Constructs a fullRayPosition with the given values. 
		///
		///\param x_ The radial coordinate of the ray position, in meters
		///\param z_ The vertical coordinate of the ray position, in meters
		///\param theta_ The angle of the ray, in radians (from the upward vertical)
		///\param time_ The time to traverse the path, in seconds
		///\param attenuation_ The factor by which signal amplitude is decreased along the path
		fullRayPosition(double x_, double z_, double theta_, double time_, double attenuation_);
		
		///Overloaded increment operator
		///\param p The position whose coordinates to add to this position's coordinates
		fullRayPosition& operator +=(const fullRayPosition& p);
		
		///Overloaded decrement operator
		///\param p The position whose coordinates to subtract this position's coordinates
		fullRayPosition& operator -=(const fullRayPosition& p);
		
		///Overloaded scaling operator
		///\param m The factor by which to scale this position's coordinates
		fullRayPosition& operator *=(double m);
		
		///Overloaded addition operator
		///\param p The position whose coordinates to add to this position's coordinates
		const fullRayPosition operator +(const fullRayPosition& p) const;
		
		///Overloaded subtraction operator
		///\param p The position whose coordinates to subtract from this position's coordinates
		const fullRayPosition operator -(const fullRayPosition& p) const;
		
		///Overloaded multiplication operator
		///\param m The factor by which to scale this position's coordinates
		const fullRayPosition operator *(double m) const;
		
		///\brief Replaces all of the coordinate values with the given value. 
		///\param s The value to use for all coordinates. 
		void makeTiny(double s);
		
		///\brief Copies information from the position object to a TraceRecord object
		///\param trace The trace record object into which to place data from this position. 
		void giveData(TraceRecord& trace) const;
		
		///\brief Adds flight time to the position. 
		///\param t The time, in seconds, to add. 
		void addTime(double t){ time+=t; }
		
		typedef fullRayPosition derivativeType;
	};
	///Overloaded mutiplication operator for scaling a fullRayPosition. 
	///\param m The factor by which to scale the coordinates of the position
	///\param p The position to scale. 
	const fullRayPosition operator *(double m, const fullRayPosition& p);
	
	//Stuff related to saving path information
	
	///For internal use by TraceFinder
	struct rkStepRecord{
		double length;
		double z[6];
	};
	
	///For internal use by TraceFinder
	template<typename positionType>
	struct positionRecordingWrapper : public positionType{
		rkStepRecord rec;
		
		positionRecordingWrapper<positionType>(){
			rec.length=0.0;
			std::fill(&rec.z[0],&rec.z[6],0.0);
		}
		positionRecordingWrapper<positionType>(const positionType& pos):positionType(pos){
			rec.length=0.0;
			std::fill(&rec.z[0],&rec.z[6],0.0);
		}
		
		void recordStep(unsigned int step){
			rec.z[step]=positionType::z;
		}
		void setFirstStep(double z){
			rec.z[0]=z;
		}
		void setStepLength(double len){
			rec.length=len;
		}
		void takeStepData(const positionRecordingWrapper<positionType>& p){
			rec=p.rec;
		}
		
		positionType& operator=(const positionType& p){
			return(positionType::operator=(p));
		}
		
		typedef positionType derivativeType;
	};
	
	enum RKStepType{RK_FIRST_STEP, RK_AIR_STEP, RK_STEP, RK_REFLECT_STEP};
	
	struct stepRecord{
		RKStepType stepType;
		union{
			rkStepRecord rkData;
			double angle;
		};
		explicit stepRecord(const rkStepRecord& rkStep):stepType(RK_STEP){
			rkData=rkStep;
		}
		stepRecord(RKStepType type, double a):stepType(type),angle(a){}
	};
	
	class traceReplayRecord{
	private:
		std::vector<stepRecord> steps;
	public:
		template<typename positionType> friend class pathRecorder;
		friend class TraceFinder;
	};
	
	template<typename positionType>
	class pathRecorder{
	private:
		traceReplayRecord record;
	public:
		pathRecorder(){}
		void operator()(const positionRecordingWrapper<positionType>& p, RKStepType stepType){
			switch(stepType){
				case RK_FIRST_STEP:
					//do nothing
					break;
				case RK_AIR_STEP:
					//negate theta since we have other means of knowing that this is refraction:
					record.steps.push_back(stepRecord(RK_AIR_STEP,-p.theta));
					break;
				case RK_STEP:
					record.steps.push_back(stepRecord(p.rec));
					break;
				case RK_REFLECT_STEP:
					record.steps.push_back(stepRecord(RK_REFLECT_STEP,p.theta));
					break;
			}
		}
		const traceReplayRecord& getData(){
			return(record);
		}
		void clearData(){
			record.steps.clear();
		}
	};

	///\brief A class which uses models of ice parameters to calculate ray paths through the ice. 
	///
	///
	class TraceFinder{
	protected:
        #ifndef __CINT__
		///The model of the ice index of refraction used for path calculations
		boost::shared_ptr<const indexOfRefractionModel> rModel;
		
		///The model of the attenuation of radio signals used for path calculations
		boost::shared_ptr<const attenuationModel> aModel;
        #endif
		
		///\brief Computes the derivatives of the position coordinates with respect to path length
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. 
		///
		///\param pos The position at which the derivatives are to be computed
		///\param der The object into which to record the calculated derivatives
		///\param frequency The frequency of the signal being propogated, in GHz
		template<typename positionType>
		void computeRayDerivatives(const positionType& pos, typename positionType::derivativeType& der, double frequency) const;
		
		///\brief Computes the derivatives of the position coordinates with respect to path length
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. (This overload is necessary because function 
		///templates canot be partially specialized.)
		///
		///\param pos The position at which the derivatives are to be computed
		///\param der The object into which to record the calculated derivatives
		///\param frequency The frequency of the signal being propogated, in GHz
		template<typename positionType>
		void computeRayDerivatives(const positionRecordingWrapper<positionType>& pos, typename positionRecordingWrapper<positionType>::derivativeType& der, double frequency) const{
			computeRayDerivatives(static_cast<const positionType&>(pos),der,frequency);
		}
		
		///\brief Computes the maximum of a set of coordinate errors, scaled by given factors
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. 
		///
		///\param errors The coordinate errors
		///\param scale The scaling factors for the errors
		///\return The maximum error divided by its associated scale
		template<typename positionType>
		double maxError(const positionType& errors, const positionType& scale) const;
		
		///\brief Computes the maximum of a set of coordinate errors, scaled by given factors
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. (This overload is necessary because function 
		///templates canot be partially specialized.)
		///
		///\param errors The coordinate errors
		///\param scale The scaling factors for the errors
		///\return The maximum error divided by its associated scale
		template<typename positionType>
		double maxError(const positionRecordingWrapper<positionType>& errors, const positionRecordingWrapper<positionType>& scale) const{
			return(maxError(static_cast<const positionType&>(errors),static_cast<const positionType&>(scale)));
		}
		
		///\brief Takes one Runge-Kutta step from the given current position, reporting estimated errors on the result position. 
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. 
		///
		///\param length The path length traversed so far, which will be updated with the length of the step
		///\param pos The current position
		///\param der The derivatives of the position coordinates at the current position
		///\param h The length of the step to take
		///\param newPos The resulting position
		///\param errors The estimated errors on the step
		///\param frequency The frequency of the signal being propagated, in GHz
		template<typename positionType>
		void rkStep(const positionType& pos, const typename positionType::derivativeType& der, const double h, positionType& newPos, positionType& errors, const double frequency) const;
		
		///\brief Takes one Runge-Kutta step from the given current position, reporting estimated errors on the result position. 
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. This overload also records information so that the
		///step can be reproduced later
		///
		///\param length The path length traversed so far, which will be updated with the length of the step
		///\param pos The current position
		///\param der The derivatives of the position coordinates at the current position
		///\param h The length of the step to take
		///\param newPos The resulting position
		///\param errors The estimated errors on the step
		///\param frequency The frequency of the signal being propagated, in GHz
		template<typename positionType>
		void rkStep(const positionRecordingWrapper<positionType>& pos, const typename positionRecordingWrapper<positionType>::derivativeType& der, const double h, positionRecordingWrapper<positionType>& newPos, positionRecordingWrapper<positionType>& errors, const double frequency) const;
		
		void replayRkStep(const rkStepRecord& step, const double atten, const double& attenDer, double& newAtten, const double frequency) const;
		
		///\brief Takes one Runge-Kutta step choosing step size to control estimated errors. 
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. 
		///
		///\param length The path length traversed so far, which will be updated with the length of the step
		///\param frequency The frequency of the signal being propagated, in GHz
		///\param pos The current position
		///\param der The derivatives of the position coordinates at the current position
		///\param scale The maximum sizes allowable for the errors on each of the coordinates
		///\param htry The length of the step to attempt
		///\param eps Overall scale for estimated errors
		///\param hdid The variable into which to record the length of the step which was actually taken
		///\param hnext The varaible into thich to record the recommended length for the next step
		template<typename positionType>
		void rkStepControl(double& length, double frequency, positionType& pos, typename positionType::derivativeType& der, const positionType& scale, const double htry, const double eps, double& hdid, double& hnext) const;
		
		///\brief Core implementation shared by traceRoot and refineRoot
		///
		///\param emit_depth The depth of the source in the ice, in meters
		///\param target The description of the position of the target
		///\param minAngle The minimum launch angle to be considered, in radians
		///\param maxAngle The maximum launch angle to be considered, in radians
		///\param rising Whether the root being sought is on a rising slope
		///\param allowedReflections Which types of reflections are allowed
		///\param a The left bound dof the search interval
		///\param aTrace The trace computed at angle a
		///\param c The right bound dof the search interval
		///\param cTrace The trace computed at angle c
		///\param requiredAccuracy The maximum distance by which the result may miss the target
		///\return The initial angle of the best solution ray path
		//std::pair<double,double> traceRootImpl(double emit_depth, const rayTargetRecord& target, bool rising, unsigned short allowedReflections, double requiredAccuracy, 
		std::pair<double,double> traceRootImpl(double emit_depth, const rayTargetRecord& target, bool rising, unsigned short allowedReflections, double requiredAccuracy, 
							 double a, TraceRecord& aTrace, double c, TraceRecord& cTrace, double angle, int &sol_error ) const;
		
		///\brief Attempts to determine a path which misses the given target vertically by the smallest possible amount. 
		///
		///This function seeks roots of the vertical miss as a function of launch angle using a simplified version of Brent's method. 
		///
		///\param emit_depth The vertical coordinate of the source, in meters
		///\param target The description of the position of the target
		///\param minAngle The minimum launch angle to be considered, in radians
		///\param maxAngle The maximum launch angle to be considered, in radians
		///\param rising Whether the root being sought is on a rising slope
		///\param allowedReflections Which types of reflections are allowed
		///\param requiredAccuracy The maximum distance by which the result may miss the target
		///\return The initial angle of the best solution ray path
		//std::pair<double,double> traceRoot(double emit_depth, const rayTargetRecord& target, double minAngle, double maxAngle, bool rising, unsigned short allowedReflections, double requiredAccuracy) const;
		std::pair<double,double> traceRoot(double emit_depth, const rayTargetRecord& target, double minAngle, double maxAngle, bool rising, unsigned short allowedReflections, double requiredAccuracy, int &sol_error ) const;
		
		///\brief Attempts to improve an existing estimate of a path which misses the given target vertically by the smallest possible amount. 
		///
		///Given a path solution which is close to the target, but not as close as desired, this function tries to find a better one. 
		///It steps outward from the seed solution to try to bracket the root (solution with zero miss distance), and then uses a 
		///simplified version of Brent's method to seek the root. 
		///
		///\param emit_depth The vertical coordinate of the source, in meters
		///\param target The description of the position of the target
		///\param seed The existing trace which is insufficiently close to the target
		///\param rising Whether the root being sought is on a rising slope
		///\param allowedReflections Which types of reflections are allowed
		///\param requiredAccuracy The maximum distance by which the result may miss the target
		///\return The initial angle of the best solution ray path
		std::pair<double,double> refineRoot(double emit_depth, const rayTargetRecord& target, const TraceRecord& seed, bool rising, unsigned short allowedReflections, double requiredAccuracy, int &sol_error ) const;

		std::pair<double,double> evenmore_refineRoot(double emit_depth, const rayTargetRecord& target, const TraceRecord& seed, double angle_min, double angle_max, bool rising, unsigned short allowedReflections, double requiredAccuracy, int &sol_error ) const;

		
		///\brief Attempts to locate the ray path with passes the farthest above the given target
		///
		///\param emit_depth The vertical coordinate of the source, in meters
		///\param target The description of the position of the target
		///\param left The left bound of the angular search space
		///\param right The right bound of the angular search space
		//std::pair<bool, double> traceMax(double emit_depth, const rayTargetRecord& target, double left=0.0, double right=pi) const;
		std::pair<bool, double> traceMax(double emit_depth, const rayTargetRecord& target, int &sol_error, double left=0.0, double right=pi ) const;
		
		///\brief Attempts to locate the ray path with passes the farthest below the given target
		///
		///\param emit_depth The vertical coordinate of the source, in meters
		///\param target The description of the position of the target
		///\param left The left bound of the angular search space
		///\param right The right bound of the angular search space
		//std::pair<bool, double> traceMin(double emit_depth, const rayTargetRecord& target, double left=0.0, double right=pi) const;
		std::pair<bool, double> traceMin(double emit_depth, const rayTargetRecord& target, int &sol_error, double left=0.0, double right=pi ) const;
		
		///\brief Solves for rays which pass through the ice surface using fast estimates if possible
		///
		///\param sourcePos The position of the signal source (one endpoint of the path)
		///\param targetPos The position of the signal receiver (the otehr endpoint of the path)
		///\param frequency The frequency of the signal, in GHz
		///\param polarization The direction of polarization of the signal, where zero corresponds 
		///		to polarization perpeduicular to the plane of propagation
		///\param requiredAccuracy The maximum distance by which the result may miss the target
		//TraceRecord findUncontainedFast(Vector sourcePos, Vector targetPos, double frequency, double polarization, double requiredAccuracy, pathRecorder<fullRayPosition>* recorder=NULL) const;
		TraceRecord findUncontainedFast(Vector sourcePos, Vector targetPos, double frequency, double polarization, double requiredAccuracy, int &sol_error,  pathRecorder<fullRayPosition>* recorder=NULL) const;
		
	public:
		///The thickness of the ice sheet, in meters
		static const double maximum_ice_depth;
		
		///\brief Constructs a TraceFinder with given ice models
		///
		///\param rm The index of refraction model to be used
		///\param am The attenuation model to be used

#ifndef __CINT__

		TraceFinder(boost::shared_ptr<const indexOfRefractionModel> rm,boost::shared_ptr<const attenuationModel> am):rModel(rm),aModel(am){}
#endif		
		
		///\brief Attempts to determine all paths which pass between the given source and target.
		///
		///Paths are returned in order of increasing length. 
		///
		///\param sourcePos The position of the signal source (one endpoint of the path)
		///\param targetPos The position of the signal receiver (the otehr endpoint of the path)
		///\param frequency The frequency of the signal, in GHz
		///\param polarization The polarization angle of the signal, in radians
		///\param allowedReflections Which types of reflections are allowed
		///\param requiredAccuracy The maximum distance by which the result may miss the target



		//std::vector<TraceRecord> findPaths(Vector sourcePos, Vector targetPos, double frequency, double polarization, unsigned short allowedReflections=NoReflection, double requiredAccuracy=0.1, std::vector<traceReplayRecord>* replayBuffer=NULL) const;
		//std::vector<TraceRecord> findPaths(Vector sourcePos, Vector targetPos, double frequency, double polarization, int &sol_cnt, unsigned short allowedReflections=NoReflection, double requiredAccuracy=0.1, std::vector<traceReplayRecord>* replayBuffer=NULL) const;
		std::vector<TraceRecord> findPaths(Vector sourcePos, Vector targetPos, double frequency, double polarization, int &sol_cnt, int &sol_error, unsigned short allowedReflections=NoReflection, double requiredAccuracy=0.1, std::vector<traceReplayRecord>* replayBuffer=NULL) const;
		//std::vector<TraceRecord> findPaths(Vector sourcePos, Vector targetPos, double frequency, double polarization, int &sol_cnt, int &sol_error, int mode, unsigned short allowedReflections=NoReflection, double requiredAccuracy=0.1, std::vector<traceReplayRecord>* replayBuffer=NULL) const;
		



		///\brief Traces one path and reports its properties. 
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. Using minimalRayPosition will result in a minimum 
		///of calculations, but the results will include only valid path length, launch angle, receipt 
		///angle, reflection angle, and miss distance values, while using fullRayPosition will give 
		///all possible results. 
		///
		///At present, this function only handles rays which are entirely inside the ice, or propagate into the ice from above. 
		///
		///This function does not work properly for very nearly vertical paths, at the moment doVerticalTrace must be used instead. 
		///
		///\param depth The vertical coordinate of the starting point of the path
		///\param theta The starting angle of the path
		///\param target The description of the position of the target
		///\param allowedReflections Which types of reflections are allowed
		///\param frequency The frequency of the signal, in GHz
		///\param polarization The polarization angle of the signal, in radians
		template<typename positionType>
		//TraceRecord doTrace(double depth, double theta, const rayTargetRecord& target, unsigned short allowedReflections, double frequency, double polarization) const;
		TraceRecord doTrace(double depth, double theta, const rayTargetRecord& target, unsigned short allowedReflections, double frequency, double polarization, int &sol_error ) const;
		
		///\brief Traces one path and reports its properties. 
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. Using minimalRayPosition will result in a minimum 
		///of calculations, but the results will include only valid path length, launch angle, receipt 
		///angle, reflection angle, and miss distance values, while using fullRayPosition will give 
		///all possible results. 
		///
		///At present, this function only handles rays which are entirely inside the ice, or propagate into the ice from above. 
		///
		///This function does not work properly for very nearly vertical paths, at the moment doVerticalTrace must be used instead. 
		///
		///This overload is templated to accept a callback, which can be a function pointer, function object, 
		///or function object pointer. After each step in the trace the callback will be invoked with the current
		///position (whose type will be positionType) and an RKStepType value descibing what kind of step has 
		///just been taken. Note that the callback is passed by value, so it will be copied. This means that 
		///caution should be used if the callback is a stateful function object (as the accumulated state in 
		///the copy will not be available in the calling context). 
		///
		///\param depth The vertical coordinate of the starting point of the path
		///\param theta The starting angle of the path
		///\param target The description of the position of the target
		///\param allowedReflections Which types of reflections are allowed
		///\param frequency The frequency of the signal, in GHz
		///\param polarization The polarization angle of the signal, in radians
		///\param callback A function pointer, function object or function object 
		///		pointer which will be caled after each step in the trace with the 
		///		current trace position and the type of the step
		template<typename positionType, typename callbackType>
		//TraceRecord doTrace(double depth, double theta, const rayTargetRecord& target, unsigned short allowedReflections, double frequency, double polarization, callbackType callback) const;
		TraceRecord doTrace(double depth, double theta, const rayTargetRecord& target, unsigned short allowedReflections, double frequency, double polarization, int &sol_error, callbackType callback) const;
		
		double recalculateAmplitude(const traceReplayRecord& trace, double frequency, double polarization);
		
		///\brief Traces one vertical path and reports its properties. 
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. Using minimalRayPosition will result in a minimum 
		///of calculations, but the results will include only valid path length, launch angle, receipt 
		///angle, reflection angle, and miss distance values, while using fullRayPosition will give 
		///all possible results. 
		///
		///This function is intended only to handle paths which are very close to vertical. 
		///
		///\param depth The depth of the starting point of the path
		///\param theta The starting angle of the path
		///\param target The description of the position of the target
		///\param allowedReflections Which types of reflections are allowed
		///\param frequency The frequency of the signal, in GHz
		///\param polarization The polarization angle of the signal, in radians
		template<typename positionType>
		TraceRecord doVerticalTrace(double depth, double theta, const rayTargetRecord& target, unsigned short allowedReflections, double frequency, double polarization) const;
		
		///\brief Traces one vertical path and reports its properties. 
		///
		///This function is templated in the positionType it uses in order to facilitate selecting 
		///exactly which coordinates are computed. Using minimalRayPosition will result in a minimum 
		///of calculations, but the results will include only valid path length, launch angle, receipt 
		///angle, reflection angle, and miss distance values, while using fullRayPosition will give 
		///all possible results. 
		///
		///This function is intended only to handle paths which are very close to vertical. 
		///
		///This overload is templated to accept a callback, which can be a function pointer, function object, 
		///or function object pointer. After each step in the trace the callback will be invoked with the current
		///position (whose type will be positionType) and an RKStepType value descibing what kind of step has 
		///just been taken. Note that the callback is passed by value, so it will be copied. This means that 
		///caution should be used if the callback is a stateful function object (as the accumulated state in 
		///the copy will not be available in the calling context). 
		///
		///\param depth The depth of the starting point of the path
		///\param theta The starting angle of the path
		///\param target The description of the position of the target
		///\param allowedReflections Which types of reflections are allowed
		///\param frequency The frequency of the signal, in GHz
		///\param polarization The polarization angle of the signal, in radians
		///\param callback A function pointer, function object or function object 
		///		pointer which will be caled after each step in the trace with the 
		///		current trace position and the type of the step
		template<typename positionType, typename callbackType>
		TraceRecord doVerticalTrace(double depth, double theta, const rayTargetRecord& target, unsigned short allowedReflections, double frequency, double polarization, callbackType callback) const;
		
		///Gets the index of refraction model used by this TraceFinder
		///\return The index of refraction model

#ifndef __CINT__
		boost::shared_ptr<const indexOfRefractionModel> getRefractionModel() const{
			return(rModel);
		}
		
		///Gets the attenuation model used by this TraceFinder
		///\return The attenutation model
		boost::shared_ptr<const attenuationModel> getAttenuationModel() const{
			return(aModel);
		}
#endif		
		
		///\brief Computes the relative signal strength for a signal traveling between two points. 
		///
		///The result of this function describes the change in signal strength along the given path 
		///due to the spreading (or focusing) of rays. In a uniform medium without attenuation, this 
		///is the usual 1/r dependence. 
		///
		///To compute an actual detected signal strength, the result of 
		///this function should be multiplied by the source intensity, the source angular dependence 
		///(evaluated at the initial angle of the ray), the receiver angular dependence (evaluated 
		///at the final angle of the ray), and the attenuation along the ray path. 
		///
		///\param ray An already computed ray which passes from src to trg
		///\param src The starting position of the path
		///\param trg The ending position of the path
		///\param allowedReflections Which types of reflections are allowed
		///\return The unitless relative electric field strength arriving at the target
		double signalStrength(const TraceRecord& ray, const Vector& src, const Vector& trg, unsigned short allowedReflections, int &sol_error ) const;
	};

} //namespace RayTrace

#include "RayTrace_detail/RayTraceImplementation.h"

#endif
