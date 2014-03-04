#ifndef RAYTRACE_ICEMODELS_H
#define RAYTRACE_ICEMODELS_H

#include "RayTrace.h"

///\brief Implements an exponential parameterization for the ice index of refraction. 
class exponentialRefractiveIndex : public RayTrace::indexOfRefractionModel{
protected:
	const double A,B,C;
public:
	///Constructs the index of refraction model with the given parameters.
	///
	///\param n_surface The value of the index of refraction at the surface of the ice. 
	///\param n_deep The value the index of refraction approaches as the depth goes to infinity. 
	///\param transition The scale of the exponential transition, in inverse meters. 
	exponentialRefractiveIndex(double n_surface, double n_deep, double transition);
	
	virtual double indexOfRefraction(double z) const;
	virtual double indexOfRefractionDerivative(double z) const;
	virtual void indexOfRefractionWithDerivative(double z, double& n, double& dndz) const;
	virtual RayTrace::indexOfRefractionModel::RayEstimate estimateRayAngle(double sourceDepth, double receiverDepth, double distance) const;
};

class inverseExponentialRefractiveIndex : public RayTrace::indexOfRefractionModel{
protected:
	const double A,B,C;
public:
	///Constructs the index of refraction model with the given parameters.
	///
	///\param n_surface The value of the index of refraction at the surface of the ice. 
	///\param n_deep The value the index of refraction approaches as the depth goes to infinity. 
	///\param transition The scale of the exponential transition, in inverse meters.
	inverseExponentialRefractiveIndex(double n_surface, double n_deep, double transition);
	
	virtual double indexOfRefraction(double z) const;
	virtual double indexOfRefractionDerivative(double z) const;
	virtual void indexOfRefractionWithDerivative(double z, double& n, double& dndz) const;
    virtual RayTrace::indexOfRefractionModel::RayEstimate estimateRayAngle(double sourceDepth, double receiverDepth, double distance) const;
};

class simpleExponentialRefractiveIndex : public RayTrace::indexOfRefractionModel{
protected:
	const double A,B;
public:
	simpleExponentialRefractiveIndex(double a, double b);
	
	virtual double indexOfRefraction(double z) const;
	virtual double indexOfRefractionDerivative(double z) const;
	virtual void indexOfRefractionWithDerivative(double z, double& n, double& dndz) const;
};

class quadraticRefractiveIndex : public RayTrace::indexOfRefractionModel{
protected:
	const double A,B,C;
	const double maxPoint, maxVal;
public:
	quadraticRefractiveIndex(double a, double b, double c);
	
	virtual double indexOfRefraction(double z) const;
	virtual double indexOfRefractionDerivative(double z) const;
	virtual void indexOfRefractionWithDerivative(double z, double& n, double& dndz) const;
};

double todorDensity(double z);
double approxTodorDensity(double z);
double approxTodorDensityDerivative(double z);

class todorLinearRefractiveIndex : public RayTrace::indexOfRefractionModel{
protected:
	const double nref0,rho0;
public:
	todorLinearRefractiveIndex(double a, double b);
	
	virtual double indexOfRefraction(double z) const;
	virtual double indexOfRefractionDerivative(double z) const;
	virtual void indexOfRefractionWithDerivative(double z, double& n, double& dndz) const;
};

class todorChiRefractiveIndex : public RayTrace::indexOfRefractionModel{
protected:
	const double nref0,rho0;
public:
	todorChiRefractiveIndex(double a, double b);
	
	virtual double indexOfRefraction(double z) const;
	virtual double indexOfRefractionDerivative(double z) const;
	virtual void indexOfRefractionWithDerivative(double z, double& n, double& dndz) const;
};

class todorLLRefractiveIndex : public RayTrace::indexOfRefractionModel{
protected:
	const double nref0,rho0;
public:
	todorLLRefractiveIndex(double a, double b);
	
	virtual double indexOfRefraction(double z) const;
	virtual double indexOfRefractionDerivative(double z) const;
	virtual void indexOfRefractionWithDerivative(double z, double& n, double& dndz) const;
};

///\brief Implements an attenuation model with no attenuation. 
class negligibleAttenuationModel : public RayTrace::attenuationModel{
public:
	virtual double attenuationLength(double z, double frequency) const;
};

class basicAttenuationModel : public RayTrace::attenuationModel{
public:
	///Computes the ice temperature at a given depth. 
	///
	///\param z The depth at which to calculate the temperature. 
	///\return The ice temperature, in degrees Celsius. 
	double temperature(double z) const;
	virtual double attenuationLength(double z, double frequency) const;
};

class constantRefractiveIndex : public RayTrace::indexOfRefractionModel{
protected:
    const double n_fixed;
public:
    constantRefractiveIndex(double n);
    virtual double indexOfRefraction(double z) const;
    virtual double indexOfRefractionDerivative(double z) const;
	virtual void indexOfRefractionWithDerivative(double z, double& n, double& dndz) const;
};

#endif
