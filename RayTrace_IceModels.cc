#include "RayTrace_IceModels.h"

#include <limits>
#include <stdexcept>

exponentialRefractiveIndex::exponentialRefractiveIndex(double n_surface, double n_deep, double transition):
A(n_deep),B(n_surface-n_deep),C(transition){}

double exponentialRefractiveIndex::indexOfRefraction(double z) const{
	if(z>0.0)
		return(1.0);
	return(A+B*exp(C*z));
}

double exponentialRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z>0.0)
		return(0.0);
	return(B*C*exp(C*z));
}

void exponentialRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z>0.0){
		n=1.0;
		dndz=0.0;
	}
	else{
		n=B*exp(C*z);
		dndz=C*n;
		n+=A;
	}
}

RayTrace::indexOfRefractionModel::RayEstimate exponentialRefractiveIndex::estimateRayAngle(double sourceDepth, double receiverDepth, double distance) const {
	if(B==0.0){
		//in this degenerate case, n(z)==A for all z
		//which causes problems when we divide by (n-A) and (n0-A) below
		//however, this case is really simple, so we can handle it directly
		double theta=atan(distance/(sourceDepth-receiverDepth));
		if(theta<0.0)
			theta+=RayTrace::pi;
		return(RayEstimate(SOLUTION,theta)); // return status: SOLUTION (which is 0), angle: theta
	}
            
	double n0=A+B*exp(C*sourceDepth);
        //std::cout<<"source n0 : "<<n0<<std::endl; // changed
	double n=A+B*exp(C*receiverDepth);
        //std::cout<<"receiver n : "<<n<<std::endl; // changed
	bool swap=(n<n0);
	if(swap) {
            std::swap(n,n0);
            //std::cout<<"swap!"<<std::endl; // changed (so that src is always near surface?)
        }
	double s1=1e-10;
	double s2=1.0;
	
	double sDiff=C*distance;
	
	double a,b,c,s;
	
	double f;
	a=s1*s1*n0*n0;
	b=sqrt(A*A-a);
	c=A/b;
	f=log((((sqrt(n*n-a)+b)/(n-A))+c)/(((sqrt(n0*n0-a)+b)/(n0-A))+c))-((b*sDiff)/(s1*n0));
	
	double fmid;
	a=s2*s2*n0*n0;
	b=sqrt(A*A-a);
	c=A/b;
	fmid=log((((sqrt(n*n-a)+b)/(n-A))+c)/(((sqrt(n0*n0-a)+b)/(n0-A))+c))-((b*sDiff)/(s2*n0));
	
	if(f*fmid>=0.0 || std::isnan(f) || std::isnan(fmid)) {
            //std::cout<<"f*fmin>=0 | f nan | fmid nan"<<std::endl; // changed
            return(RayEstimate());
        }
	
	double ds;
	double rtb=(f<0.0?(ds=s2-s1,s1):(ds=s1-s2,s2));
        // changed
        //std::cout<<"f:"<<f<<std::endl;
        //std::cout<<"s1:"<<s1<<", s2:"<<s2<<std::endl;
        //std::cout<<"s2-s1:"<<s2-s1<<", s1-s2:"<<s1-s2<<std::endl;
        //std::cout<<"ds:"<<ds<<", rtb:"<<rtb<<std::endl;
	const unsigned int maxIter=40;
	unsigned int i;
	for(i=0; i<maxIter; i++){
		ds*=0.5;
		s=rtb+ds;
		a=s*s*n0*n0;
		b=sqrt(A*A-a);
		c=A/b;
		fmid=log((((sqrt(n*n-a)+b)/(n-A))+c)/(((sqrt(n0*n0-a)+b)/(n0-A))+c))-((b*sDiff)/(s*n0));
		if(fmid<=0.0)
			rtb=s;
		if(std::abs(fmid)<1.0e-4)
			break;
	}
	if(i==maxIter) {
                //std::cout<<"i=maxIter"<<std::endl; // changed
                return(RayEstimate());
        }
	if(swap){
		if(sourceDepth>receiverDepth) {
                        //std::cout<<"swap & src_depth("<<sourceDepth<<") > rec_depth("<<receiverDepth<<")"<<std::endl; // changed
			return(RayEstimate(SOLUTION,RayTrace::pi-asin((n0/n)*s)));
                }
                //std::cout<<"swap & src_depth("<<sourceDepth<<") <= rec_depth("<<receiverDepth<<")"<<std::endl; // changed
		return(RayEstimate(SOLUTION,asin((n0/n)*s)));
	}
	if(sourceDepth<receiverDepth) {
                //std::cout<<"src_depth("<<sourceDepth<<") < rec_depth("<<receiverDepth<<")"<<std::endl; // changed
		return(RayEstimate(SOLUTION,asin(s)));
        }
        //std::cout<<"just return from RayEstimate"<<std::endl; // changed
	return(RayEstimate(SOLUTION,RayTrace::pi-asin(s)));
}


inverseExponentialRefractiveIndex::inverseExponentialRefractiveIndex(double n_surface, double n_deep, double transition):
A(2*n_surface-n_deep),B(2.*(n_deep-n_surface)),C(transition){}

double inverseExponentialRefractiveIndex::indexOfRefraction(double z) const{
	if(z>0.0)
		return(1.0);
	return(A+B/(1.+exp(C*z)));
}

double inverseExponentialRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z>0.0)
		return(0.0);
	return(-B*C*exp(C*z)/((1.+exp(C*z))*(1.+exp(C*z))));
}

void inverseExponentialRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z>0.0){
		n=1.0;
		dndz=0.0;
	}
	else{
		double e=exp(C*z);
		n=A+B/(1.+e);
		dndz=-B*C*e/((1.+e)*(1.+e));
	}
}

RayTrace::indexOfRefractionModel::RayEstimate inverseExponentialRefractiveIndex::estimateRayAngle(double sourceDepth, double receiverDepth, double distance) const{
	if(B==0.0){
		//in this degenerate case, n(z)==A for all z
		//which causes problems when we divide by (n-A) and (n0-A) below
		//however, this case is really simple, so we can handle it directly
		double theta=atan(distance/(sourceDepth-receiverDepth));
		if(theta<0.0)
			theta+=RayTrace::pi;
		return(RayEstimate(SOLUTION,theta));
	}
	
    double n0=A+B/(1.+exp(C*sourceDepth));
    double n=A+B/(1.+exp(C*receiverDepth));
//	double n0=A+B*exp(C*sourceDepth);
//	double n=A+B*exp(C*receiverDepth);
	bool swap=(n<n0);
	if(swap)
		std::swap(n,n0);
	double s1=1e-10;
	double s2=1.0;
	
	double sDiff=C*distance;
	
	double a,b,c,s;
	
	double f;
	a=s1*s1*n0*n0;
	b=sqrt(A*A-a);
	c=A/b;
	f=log((((sqrt(n*n-a)+b)/(n-A))+c)/(((sqrt(n0*n0-a)+b)/(n0-A))+c))-((b*sDiff)/(s1*n0));
	
	double fmid;
	a=s2*s2*n0*n0;
	b=sqrt(A*A-a);
	c=A/b;
	fmid=log((((sqrt(n*n-a)+b)/(n-A))+c)/(((sqrt(n0*n0-a)+b)/(n0-A))+c))-((b*sDiff)/(s2*n0));
	
	if(f*fmid>=0.0 || std::isnan(f) || std::isnan(fmid))
		return(RayEstimate());
	
	double ds;
	double rtb=(f<0.0?(ds=s2-s1,s1):(ds=s1-s2,s2));
	const unsigned int maxIter=40;
	unsigned int i;
	for(i=0; i<maxIter; i++){
		ds*=0.5;
		s=rtb+ds;
		a=s*s*n0*n0;
		b=sqrt(A*A-a);
		c=A/b;
		fmid=log((((sqrt(n*n-a)+b)/(n-A))+c)/(((sqrt(n0*n0-a)+b)/(n0-A))+c))-((b*sDiff)/(s*n0));
		if(fmid<=0.0)
			rtb=s;
		if(std::abs(fmid)<1.0e-4)
			break;
	}
	if(i==maxIter)
		return(RayEstimate());
	if(swap){
		if(sourceDepth>receiverDepth)
			return(RayEstimate(SOLUTION,RayTrace::pi-asin((n0/n)*s)));
		return(RayEstimate(SOLUTION,asin((n0/n)*s)));
	}
	if(sourceDepth<receiverDepth)
		return(RayEstimate(SOLUTION,asin(s)));
	return(RayEstimate(SOLUTION,RayTrace::pi-asin(s)));
}

simpleExponentialRefractiveIndex::simpleExponentialRefractiveIndex(double a, double b):
A(a),B(b){}

double simpleExponentialRefractiveIndex::indexOfRefraction(double z) const{
	if(z>0.0)
		return(1.0);
	return(A*exp(B*z));
}
double simpleExponentialRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z>0.0)
		return(0.0);
	return(A*B*exp(B*z));
}
void simpleExponentialRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z>0.0){
		n=1.0;
		dndz=0.0;
	}
	else{
		n=A*exp(B*z);
		dndz=B*n;
	}
}

quadraticRefractiveIndex::quadraticRefractiveIndex(double a, double b, double c):
A(a),B(b),C(c),maxPoint(B/(2.*C)),maxVal(A+(B*B)/(4.*C)){}

double quadraticRefractiveIndex::indexOfRefraction(double z) const{
	if(z>0.0)
		return(1.0);
	if(z<maxPoint)
		return(maxVal);
	return(A+(B-C*z)*z);
}
double quadraticRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z>0.0 || z<maxPoint)
		return(0.0);
	return(B-2.*C*z);
}
void quadraticRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z>0.0){
		n=1.0;
		dndz=0.0;
	}
	else if(z<maxPoint){
		n=maxVal;
		dndz=0.0;
	}
	else{
		n=(A+(B-C*z)*z);
		dndz=B-2*C*z;
	}
}


double todorDensity(double z){
	if (z>0)
		throw std::domain_error("todorDensity is defined only for negative z");
	float rho=1.;
	static const float par[7]={0.9283,2.375,-0.0249,-3.095,-0.0386,1.354,-0.0635};
	float f1a=par[2]*z; 
	if(f1a<75.)
		rho-=par[1]/exp(f1a);
	float f2a=par[4]*z; 
	if(f2a<75.)
		rho-=par[3]/exp(f2a);
	float f3a=par[6]*z; 
	if(f3a<75.)
		rho-=par[5]/exp(f3a);
	rho*=par[0];
	return(rho);
}

double approxTodorDensity(double z){
	const double a=-500.0;
	const double b=0.0;
	const double rho_a=0.928291380405426;
	const unsigned int nCoeffs=16;
	const static double coeffs[nCoeffs] = 
	{1.630184163164813,-0.1974340731867372,-0.1301940665484572,-0.0636618804359376,
		-0.02284342020414087,-0.008295837629225138,-0.007815359409374324,-0.01028047190477481,
		-0.01079129490788115,-0.009077891732173929,-0.006456928164827926,-0.004035689213371644,
		-0.002274127986029153,-0.00117587185504285,-0.0005648884951339264,-0.0002544329712926173};
	
	if(z<a)
		return(rho_a);
	double d=0.0, dd=0.0, sv, y, y2;
	y2=2.*(y=(2.*z-a-b)/(b-a));
	for(int i=nCoeffs-1; i>0; i--){
		sv=d;
		d=y2*d-dd+coeffs[i];
		dd=sv;
	}
	return(y*d-dd+.5*coeffs[0]);
}

double approxTodorDensityDerivative(double z){
	const double a=-500.0;
	const double b=0.0;
	const unsigned int nCoeffs=16;
	const static double coeffs[nCoeffs] = 
	{-0.005183430938378898,-0.004694588495944741,-0.003603958352885,-0.002611483431169426,
		-0.002076073222422498,-0.001880493984636918,-0.001744239717253493,-0.00150535673298695,
		-0.001168533290586103,-0.0008147138588825569,-0.0005149250858695804,-0.0002981596056963227,
		-0.0001597844350928757,-7.984331903752403e-05,-3.749376216841926e-05,-1.657580758252428e-05};
	
	if(z<a)
		return(0.0);
	double d=0.0, dd=0.0, sv, y, y2;
	y2=2.*(y=(2.*z-a-b)/(b-a));
	for(int i=nCoeffs-1; i>0; i--){
		sv=d;
		d=y2*d-dd+coeffs[i];
		dd=sv;
	}
	return(y*d-dd+.5*coeffs[0]);
}


todorLinearRefractiveIndex::todorLinearRefractiveIndex(double a, double b):
nref0(a),rho0(b){}

double todorLinearRefractiveIndex::indexOfRefraction(double z) const{
	if(z>0.0)
		return(1.0);
	return(1.+((nref0-1.)/rho0)*approxTodorDensity(z));
}
double todorLinearRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z>0.0 || z<-500.0)
		return(0.0);
	return(approxTodorDensityDerivative(z));
}
void todorLinearRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z>0.0){
		n=1.0;
		dndz=0.0;
	}
	else if(z<-500.0){
		n=1.+((nref0-1.)/rho0)*todorDensity(z);
		dndz=0.0;
	}
	else{
		n=1.+((nref0-1.)/rho0)*todorDensity(z);
		dndz=((nref0-1.)/rho0)*approxTodorDensityDerivative(z);
	}
}

todorChiRefractiveIndex::todorChiRefractiveIndex(double a, double b):
nref0(a),rho0(b){}

double todorChiRefractiveIndex::indexOfRefraction(double z) const{
	if(z>0.0)
		return(1.0);
	double a=(nref0*nref0-1)/rho0;
	return(sqrt(1.+a*approxTodorDensity(z)));
}
double todorChiRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z>0.0 || z<-500.0)
		return(0.0);
	double a=(nref0*nref0-1)/rho0;
	return(0.5*(a*approxTodorDensityDerivative(z))/sqrt(1.+a*approxTodorDensity(z)));
}
void todorChiRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z>0.0){
		n=1.0;
		dndz=0.0;
	}
	double a=(nref0*nref0-1)/rho0;
	if(z<-500.0){
		n=sqrt(1.+a*approxTodorDensity(z));
		dndz=0.0;
	}
	else{
		n=sqrt(1.+a*approxTodorDensity(z));
		dndz=0.5*a*approxTodorDensityDerivative(z)/n;
	}
}


todorLLRefractiveIndex::todorLLRefractiveIndex(double a, double b):
nref0(a),rho0(b){}

double todorLLRefractiveIndex::indexOfRefraction(double z) const{
	if(z>0.0)
		return(1.0);
	double a=(nref0*nref0-1)/(nref0*nref0+2)/rho0;
	double a_rho=a*approxTodorDensity(z);
	return(sqrt((1.+2.*a_rho)/(1.-a_rho)));
}
double todorLLRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z>0.0 || z<-500.0)
		return(0.0);
	double a=(nref0*nref0-1)/(nref0*nref0+2)/rho0;
	double a_rho=a*approxTodorDensity(z);
	return(3*a/(sqrt((1.+2.*a_rho)/(1.-a_rho))*(2.*a_rho*a_rho-4.*a_rho+2))*approxTodorDensityDerivative(z));
}
void todorLLRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z>0.0){
		n=1.0;
		dndz=0.0;
	}
	double a=(nref0*nref0-1)/(nref0*nref0+2)/rho0;
	double a_rho=a*approxTodorDensity(z);
	if(z<-500.0){
		n=sqrt((1.+2.*a_rho)/(1.-a_rho));
		dndz=0.0;
	}
	else{
		n=sqrt((1.+2.*a_rho)/(1.-a_rho));
		dndz=(3*a/(n*(2.*a_rho*a_rho-4.*a_rho+2))*approxTodorDensityDerivative(z));
	}
}


double negligibleAttenuationModel::attenuationLength(double z, double frequency) const{
	return(std::numeric_limits<double>::infinity());
}


double basicAttenuationModel::temperature(double z) const{
	return(-51.5 + z*(-4.5319e-3 + 5.822e-6*z));
}

double basicAttenuationModel::attenuationLength(double z, double frequency) const{
	if(z>0.0)
		return(std::numeric_limits<double>::infinity());
	double t = temperature(z);
	const double f0=0.0001, f2=3.16;
	const double w0=log(f0), w1=0.0, w2=log(f2), w=log(frequency);
	const double b0=-6.74890+t*(0.026709-t*0.000884);
	const double b1=-6.22121-t*(0.070927+t*0.001773);
	const double b2=-4.09468-t*(0.002213+t*0.000332);
	double a,bb;
	if(frequency<1.){
		a=(b1*w0-b0*w1)/(w0-w1);
		bb=(b1-b0)/(w1-w0);
	}
	else{
		a=(b2*w1-b1*w2)/(w1-w2);
		bb=(b2-b1)/(w2-w1);
	}
	return 1./exp(a+bb*w);
}

constantRefractiveIndex::constantRefractiveIndex(double n): n_fixed(n){}

double constantRefractiveIndex::indexOfRefraction(double z) const {
    //n_fixed = 1.5;
    if (z>0.) {
        return 1.0;
    } else
        return n_fixed;
}

    
double constantRefractiveIndex::indexOfRefractionDerivative(double z) const{
	if(z>0.0)
		return(0.0);
	return(0.0);
}

void constantRefractiveIndex::indexOfRefractionWithDerivative(double z, double& n, double& dndz) const{
	if(z>0.0){
		n=1.0;
		dndz=0.0;
	}
	else{
		n=n_fixed;
		dndz=0.0;
	}
}





