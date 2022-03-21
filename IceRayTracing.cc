/*
  This is the IceRayTracing namespace. Author: Uzair Latif 
  released under GPL3.
*/

#include "IceRayTracing.hh"

/* Get the value of the B parameter for the refractive index model */
double IceRayTracing::GetB(double z){
  z=fabs(z);
  double B=0;

  B=-0.43;
  
  // if(z<=IceRayTracing::TransitionBoundary){
  //   B=-0.5019;
  // }else{
  //   B=-0.448023;
  // }
  
  return B;
}

/* Get the value of the C parameter for the refractive index model */
double IceRayTracing::GetC(double z){
  z=fabs(z);
  double C=0;

  C=0.0132;
 
  // if(z<=IceRayTracing::TransitionBoundary){
  //   C=0.03247;
  // }else{
  //   C=0.02469;
  // }
 
  return C;
}

/* Get the value of refractive index model for a given depth  */
double IceRayTracing::Getnz(double z){
  z=fabs(z);
  return IceRayTracing::A_ice+IceRayTracing::GetB(z)*exp(-IceRayTracing::GetC(z)*z);
}

/* E-feild Power Fresnel coefficient for S-polarised wave which is perpendicular to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R */
double IceRayTracing::Refl_S(double thetai){

  double Nair=1;
  double Nice=IceRayTracing::Getnz(0); 
  double n1=Nice;
  double n2=Nair;
  
  double sqterm=sqrt(1-pow((n1/n2)*(sin(thetai)),2));
  double num=n1*cos(thetai)-n2*sqterm;
  double den=n1*cos(thetai)+n2*sqterm;
  double RS=(num*num)/(den*den);

  if(std::isnan(RS)){
    RS=1;
  }
  return (RS);
}

/* E-feild Power Fresnel coefficient for P-polarised wave which is parallel to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R */
double IceRayTracing::Refl_P(double thetai){
   
  double Nair=1;
  double Nice=IceRayTracing::Getnz(0); 
  double n1=Nice;
  double n2=Nair;

  double sqterm=sqrt(1-pow((n1/n2)*(sin(thetai)),2));
  double num=n1*sqterm-n2*cos(thetai);
  double den=n1*sqterm+n2*cos(thetai);
  double RP=(num*num)/(den*den);
  if(std::isnan(RP)){
    RP=1;
  }
  return (RP);
}

/* The temperature and attenuation model has been taken from AraSim which also took it from here http://icecube.wisc.edu/~araproject/radio/ . This is basically Matt Newcomb's icecube directory which has alot of information, plots and codes about South Pole Ice activities. Please read it if you find it interesting. */

/* Temperature model:The model takes in value of depth z in m and returns the value of temperature in Celsius.*/
double IceRayTracing::GetIceTemperature(double z){
  double depth=fabs(z);
  double t = 1.83415e-09*pow(depth,3) + (-1.59061e-08*pow(depth,2)) + 0.00267687*depth + (-51.0696 );
  return t;
}

/* Ice Attenuation Length model: Takes in value of frequency in Ghz and depth z and returns you the value of attenuation length in m */
double IceRayTracing::GetIceAttenuationLength(double z, double frequency){

  double t =IceRayTracing::GetIceTemperature(z);
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
  double Lval=1./exp(a+bb*w);
  return Lval;
}

/* Setup the integrand to calculate the attenuation */
double IceRayTracing::AttenuationIntegrand (double x, void * params) {

  double *p=(double*)params;

  double A0=p[0];
  double Frequency=p[1];
  double L=p[2];

  double Integrand=(A0/IceRayTracing::GetIceAttenuationLength(x,Frequency))*sqrt(1+pow(tan(asin(L/IceRayTracing::Getnz(x))) ,2));
  return Integrand;
}

/* Integrate over the integrand to calculate the attenuation */
double IceRayTracing::IntegrateOverLAttn (double A0, double Frequency, double z0, double z1, double Lvalue) {
  gsl_integration_workspace * w= gsl_integration_workspace_alloc (1000);

  double result, error;
  double zlimit[2] = {z0,z1};
  double param[3] = {A0,Frequency,Lvalue};
  
  gsl_function F;
  F.function = &AttenuationIntegrand;
  F.params = &param;

  gsl_integration_qags (&F, zlimit[0], zlimit[1], 0, 1e-7, 1000,
                        w, &result, &error);

  // printf ("result          = % .18f\n", result);
  // printf ("estimated error = % .18f\n", error);
  // printf ("intervals       = %zu\n", w->size);

  gsl_integration_workspace_free (w);  

  return fabs(result);
}

/* Calculate the total attenuation for each type of ray */
double IceRayTracing::GetTotalAttenuationDirect (double A0, double frequency, double z0, double z1, double Lvalue) {
  z0=fabs(z0);
  z1=fabs(z1);
  return IceRayTracing::IntegrateOverLAttn(A0,frequency,z0,z1,Lvalue);
}

double IceRayTracing::GetTotalAttenuationReflected (double A0, double frequency, double z0, double z1, double Lvalue) {
  z0=fabs(z0);
  z1=fabs(z1);
  return IceRayTracing::IntegrateOverLAttn(A0,frequency,z0,0.000001,Lvalue) + IceRayTracing::IntegrateOverLAttn(A0,frequency,z1,0.000001,Lvalue);
}

double IceRayTracing::GetTotalAttenuationRefracted (double A0, double frequency, double z0, double z1, double zmax, double Lvalue) {
  z0=fabs(z0);
  z1=fabs(z1);
  return IceRayTracing::IntegrateOverLAttn(A0,frequency,z0,zmax,Lvalue) + IceRayTracing::IntegrateOverLAttn(A0,frequency,z1,zmax,Lvalue);
}

/* Use GSL minimiser which relies on calculating function deriavtives. This function uses GSL's Newton's algorithm to find root for a given function. */
double IceRayTracing::FindFunctionRootFDF(gsl_function_fdf FDF,double x_lo, double x_hi){
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double x0 ,x = (x_lo+x_hi)/2;
  
  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, x);

  // printf ("using %s method\n",
  //         gsl_root_fdfsolver_name (s));

  // printf ("%-5s %10s %10s %10s\n",
  //         "iter", "root", "err", "err(est)");
  do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 1e-6);
      
      // if (status == GSL_SUCCESS)
      //   printf ("Converged:\n");
      
      // printf ("%5d %10.7f %10.7f\n",
      //         iter, x, x - x0);
      
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fdfsolver_free (s);
  
  return x;
}

/* Use GSL minimiser which uses GSL's false position algorithm to find root for a given function. */
double IceRayTracing::FindFunctionRoot(gsl_function F,double x_lo, double x_hi)
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
 
  T = gsl_root_fsolver_falsepos;
  s = gsl_root_fsolver_alloc (T);
  gsl_set_error_handler_off();
  status = gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      
      if(x_lo>x_hi){
      	swap(x_lo,x_hi);
      }
      double checkval=(*((F).function))(r,(F).params);
      status = gsl_root_test_residual (checkval,1e-6);
      //status = gsl_root_test_interval (x_lo, x_hi,0, 0.000001);
      
      if (status == GSL_SUCCESS){
	// printf ("Converged:");
	// printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi,r,x_hi - x_lo);
      }
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  //printf ("status = %s\n", gsl_strerror (status));
  gsl_root_fsolver_free (s);

  return r;
}

/* Use GSL minimiser which uses GSL's false position algorithm to find root for a given function. */
double IceRayTracing::FindFunctionRootZmax(gsl_function F,double x_lo, double x_hi)
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
 
  T = gsl_root_fsolver_falsepos;
  s = gsl_root_fsolver_alloc (T);
  gsl_set_error_handler_off();
  status = gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,1e-6, 1e-6);
	
      if (status == GSL_SUCCESS){
	// printf ("Converged:");
	// printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi,r,x_hi - x_lo);
      }
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  gsl_root_fsolver_free (s);

  return r;
}

/* Define the function that will be minimized to get the value of the depth of the turning point for a given refracted ray. This function basically requires the value of the L parameter to find the depth of the turning point.  This comes from the part of the fDnfR function where sqrt( n(z) - L ). This imposes the constraint then n(z)=L at the turning point of the ray from which we can find zmax. */
double IceRayTracing::GetMinnz(double x,void *params){
  struct IceRayTracing::Minnz_params *p= (struct IceRayTracing::Minnz_params *) params;
  double A = p->a;
  double L = p->l;
  return A+IceRayTracing::GetB(x)*exp(-IceRayTracing::GetC(x)*x)-L;
}

/* Get the value of the depth of the turning point for the refracted ray */
double IceRayTracing::GetZmax(double A, double L){
  gsl_function F1;
  struct IceRayTracing::Minnz_params params1= {A,L};
  F1.function = &GetMinnz;
  F1.params = &params1;
  double zmax=IceRayTracing::FindFunctionRootZmax(F1,0.0,5000);
  return zmax;
}

/* Analytical solution describing ray paths in ice as function of depth */
double IceRayTracing::fDnfR(double x,void *params){
  
  struct IceRayTracing::fDnfR_params *p= (struct IceRayTracing::fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*IceRayTracing::Getnz(x)-L*L+sqrt(A*A-L*L)*sqrt(pow(IceRayTracing::Getnz(x),2)-L*L)));;
}

/* Analytical solution describing the ray path in ice as a function of the L parameter */
double IceRayTracing::fDnfR_L(double x,void *params){
  
  struct IceRayTracing::fDnfR_L_params *p= (struct IceRayTracing::fDnfR_L_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Z = p->z;
  
  double result=(x/C)*(1.0/sqrt(A*A-x*x))*(C*Z-log(A*IceRayTracing::Getnz(Z)-x*x+sqrt(A*A-x*x)*sqrt(pow(IceRayTracing::Getnz(Z),2)-x*x)));
  
  return result;
}

/* The function used to calculate ray propogation time in ice */
double IceRayTracing::ftimeD(double x,void *params){

  struct IceRayTracing::ftimeD_params *p= (struct IceRayTracing::ftimeD_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Speedc = p->speedc;
  double L = p->l;
  
  return (1.0/(Speedc*C*sqrt(pow(IceRayTracing::Getnz(x),2)-L*L)))*(pow(IceRayTracing::Getnz(x),2)-L*L+(C*x-log(A*IceRayTracing::Getnz(x)-L*L+sqrt(A*A-L*L)*sqrt(pow(IceRayTracing::Getnz(x),2)-L*L)))*(A*A*sqrt(pow(IceRayTracing::Getnz(x),2)-L*L))/sqrt(A*A-L*L) +A*sqrt(pow(IceRayTracing::Getnz(x),2)-L*L)*log(IceRayTracing::Getnz(x)+sqrt(pow(IceRayTracing::Getnz(x),2)-L*L)) );
}

/* The function is used to calculate ray geometric path in ice */
double IceRayTracing::fpathD(double x,void *params){

  struct IceRayTracing::ftimeD_params *p= (struct IceRayTracing::ftimeD_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Speedc = p->speedc;
  double L = p->l;

  //integral sec(sin^(-1)(L/(A + B e^(C x)))) dx = (log((A + B e^(C x)) (sqrt((A^2 + 2 A B e^(C x) + B^2 e^(2 C x) - L^2)/(A + B e^(C x))^2) + 1)) - (A log(A sqrt(A^2 - L^2) sqrt((A^2 + 2 A B e^(C x) + B^2 e^(2 C x) - L^2)/(A + B e^(C x))^2) + B sqrt(A^2 - L^2) e^(C x) sqrt((A^2 + 2 A B e^(C x) + B^2 e^(2 C x) - L^2)/(A + B e^(C x))^2) + A^2 + A B e^(C x) - L^2))/sqrt(A^2 - L^2) + (A C x)/sqrt(A^2 - L^2))/C;
  
  return (log((A + B*exp(C*x))*(sqrt((A*A + 2*A*B*exp(C*x) + B*B*exp(2*C*x) - L*L)/((A + B*exp(C*x))*(A + B*exp(C*x))) ) + 1)) - (A*log(A*sqrt(A*A - L*L)*sqrt((A*A + 2*A*B*exp(C*x) + B*B* exp(2*C*x) - L*L)/(( A + B*exp(C*x))*(A + B*exp(C*x)))) + B*sqrt(A*A - L*L)*exp(C*x)*sqrt((A*A + 2*A*B*exp(C*x) + B*B* exp(2*C*x) - L*L)/((A + B*exp(C*x))*(A + B*exp(C*x)))) + A*A + A*B*exp(C*x) - L*L))/sqrt(A*A - L*L) + (A*C*x)/sqrt(A*A - L*L))/C ;

}

/* The set of functions starting with the name "fDa" are used in the minimisation procedure to find the launch angle (or the L parameter) for the direct ray */
double IceRayTracing::fDa(double x,void *params){
  struct IceRayTracing::fDanfRa_params *p= (struct IceRayTracing::fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  struct IceRayTracing::fDnfR_L_params params1a = {A, IceRayTracing::GetB(z1), IceRayTracing::GetC(z1), z1};
  struct IceRayTracing::fDnfR_L_params params1b = {A, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), z0};
  struct IceRayTracing::fDnfR_L_params params1c = {A, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), IceRayTracing::GetC(IceRayTracing::TransitionBoundary), -IceRayTracing::TransitionBoundary};
  struct IceRayTracing::fDnfR_L_params params1d = {A, IceRayTracing::GetB(-(IceRayTracing::TransitionBoundary+0.000001)), IceRayTracing::GetC(-(IceRayTracing::TransitionBoundary+0.000001)), -(IceRayTracing::TransitionBoundary+0.000001)};

  double distancez0z1=0;
  if(IceRayTracing::TransitionBoundary!=0){
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)>IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1c) + IceRayTracing::fDnfR_L(x,&params1d) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1c) + IceRayTracing::fDnfR_L(x,&params1d) - IceRayTracing::fDnfR_L(x,&params1b);
    }
  }else{
    distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
  }
  // if(std::isnan(distancez0z1)){
  //   distancez0z1=1e9;
  // }
  double output=distancez0z1-x1;

  return output;
}

double IceRayTracing::fDa_df(double x,void *params){
  gsl_function F;
  F.function = &IceRayTracing::fDa;
  F.params = params;
 
  double result,abserr;
  gsl_deriv_central (&F, x, 1e-8, &result, &abserr);

  return result;
}

void IceRayTracing::fDa_fdf (double x, void *params,double *y, double *dy){ 
  *y = IceRayTracing::fDa(x,params);
  *dy = IceRayTracing::fDa_df(x,params);
}

/* The set of functions starting with the name "fRa" are used in the minimisation procedure to find the launch angle (or the L parameter) for the reflected ray */
double IceRayTracing::fRa(double x,void *params){
  struct IceRayTracing::fDanfRa_params *p= (struct IceRayTracing::fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  struct IceRayTracing::fDnfR_L_params params1a = {A, IceRayTracing::GetB(z1), -IceRayTracing::GetC(z1), -z1};
  struct IceRayTracing::fDnfR_L_params params1b = {A, IceRayTracing::GetB(z0), -IceRayTracing::GetC(z0), -z0};
  struct IceRayTracing::fDnfR_L_params params1c = {A, IceRayTracing::GetB(1e-7), -IceRayTracing::GetC(1e-7), 1e-7};
  struct IceRayTracing::fDnfR_L_params params1d = {A, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary), IceRayTracing::TransitionBoundary};
  struct IceRayTracing::fDnfR_L_params params1f = {A, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+0.000001), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary+0.000001), IceRayTracing::TransitionBoundary+0.000001};

  double distancez0z1=0;
  double distancez0surface=0;
  if(IceRayTracing::TransitionBoundary!=0){
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)>IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
      distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1d) + IceRayTracing::fDnfR_L(x,&params1f) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1d) + IceRayTracing::fDnfR_L(x,&params1f) - IceRayTracing::fDnfR_L(x,&params1b);
      distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1d) + IceRayTracing::fDnfR_L(x,&params1f) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
      distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1b); 
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
      distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1b); 
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
      distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1b); 
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1d) + IceRayTracing::fDnfR_L(x,&params1f) - IceRayTracing::fDnfR_L(x,&params1b);
      distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1d) + IceRayTracing::fDnfR_L(x,&params1f) - IceRayTracing::fDnfR_L(x,&params1b); 
    }
  }else{
    distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
    distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1b);
  }
  // if(std::isnan(distancez0z1)){
  //   distancez0z1=1e9;
  // }
  // if(std::isnan(distancez0surface)){
  //   distancez0surface=1e9;
  // }
  double output= distancez0z1 - 2*(distancez0surface) - x1;
  
  return output;
}

double IceRayTracing::fRa_df(double x,void *params){
  gsl_function F;
  F.function = &IceRayTracing::fRa;
  F.params = params;
 
  double result,abserr;
  gsl_deriv_central (&F, x, 1e-8, &result, &abserr);

  return result;
}

void IceRayTracing::fRa_fdf (double x, void *params,double *y, double *dy){ 
  *y = IceRayTracing::fRa(x,params);
  *dy = IceRayTracing::fRa_df(x,params);
}

/* This function is minimised to find the launch angle (or the L parameter) for the refracted ray */
double IceRayTracing::fRaa(double x,void *params){
  struct IceRayTracing::fDanfRa_params *p= (struct IceRayTracing::fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;
  
  double zmax= IceRayTracing::GetZmax(A,x)+1e-7;
  double output=0;
  if(zmax>0){
    struct IceRayTracing::fDnfR_L_params params1a = {A, IceRayTracing::GetB(z1), -IceRayTracing::GetC(z1), -z1};
    struct IceRayTracing::fDnfR_L_params params1b = {A, IceRayTracing::GetB(z0), -IceRayTracing::GetC(z0), -z0};
    struct IceRayTracing::fDnfR_L_params params1c = {A, IceRayTracing::GetB(zmax), -IceRayTracing::GetC(zmax), zmax};
    struct IceRayTracing::fDnfR_L_params params1d = {A, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary), IceRayTracing::TransitionBoundary};
    struct IceRayTracing::fDnfR_L_params params1f = {A, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+1e-7), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary+1e-7), IceRayTracing::TransitionBoundary+1e-7};

    double distancez0z1=0;
    double distancez0surface=0;
    if(IceRayTracing::TransitionBoundary!=0){
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)>IceRayTracing::TransitionBoundary){
	if(zmax<=IceRayTracing::TransitionBoundary){
	  distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
	  distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1d) + IceRayTracing::fDnfR_L(x,&params1f) - IceRayTracing::fDnfR_L(x,&params1b);
	}else{
	  distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
	  distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1b);
	}
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1d) + IceRayTracing::fDnfR_L(x,&params1f) - IceRayTracing::fDnfR_L(x,&params1b);
	distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1d) + IceRayTracing::fDnfR_L(x,&params1f) - IceRayTracing::fDnfR_L(x,&params1b);
      }
      if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
	distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1b); 
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
	distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1b); 
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
	distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1b); 
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1d) + IceRayTracing::fDnfR_L(x,&params1f) - IceRayTracing::fDnfR_L(x,&params1b);
	distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1d) + IceRayTracing::fDnfR_L(x,&params1f) - IceRayTracing::fDnfR_L(x,&params1b); 
      }
    }else{
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
      distancez0surface=IceRayTracing::fDnfR_L(x,&params1c) - IceRayTracing::fDnfR_L(x,&params1b);
    }

    if(std::isnan(distancez0z1)){
      distancez0z1=1e9;
    }
    if(std::isnan(distancez0surface)){
      distancez0surface=1e9;
    }
    output= distancez0z1 - 2*(distancez0surface) - x1;
  }else{
    output=1e9;
  }
  return output;
}

double IceRayTracing::fRaa_df(double x,void *params){
  gsl_function F;
  F.function = &IceRayTracing::fRaa;
  F.params = params;
 
  double result,abserr;
  gsl_deriv_central (&F, x, 1e-8, &result, &abserr);

  return result;
}

void IceRayTracing::fRaa_fdf (double x, void *params,double *y, double *dy){ 
  *y = IceRayTracing::fRaa(x,params);
  *dy = IceRayTracing::fRaa_df(x,params);
}

/* This functions works for the Direct ray and gives you back the launch angle, receive angle and propagation time of the ray together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
double* IceRayTracing::GetDirectRayPar(double z0, double x1, double z1){

  double *output=new double[6];
  
  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* First we setup the fDa function that will be minimised to get the launch angle (or the L parameter) for the direct ray. */
  gsl_function F1;
  struct IceRayTracing::fDanfRa_params params1= {IceRayTracing::A_ice, z0, x1, z1};
  F1.function = &IceRayTracing::fDa;
  F1.params = &params1;

  // gsl_function_fdf F1;
  // struct IceRayTracing::fDanfRa_params params1= {IceRayTracing::A_ice, z0, x1, z1};
  // F1.f = &fDa;
  // F1.df = &fDa_df;
  // F1.fdf = &fDa_fdf;
  // F1.params = &params1;
  
  /* In my raytracing solution given in the function IceRayTracing::fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0 and z1 and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */ 
  double UpLimnz[]={IceRayTracing::Getnz(z1),IceRayTracing::Getnz(z0)};
  double* UpperLimitL=min_element(UpLimnz,UpLimnz+2);

  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function. */
  double lvalueD=IceRayTracing::FindFunctionRoot(F1,1e-7,UpperLimitL[0]);
  double LangD=asin(lvalueD/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
  double checkzeroD=IceRayTracing::fDa(lvalueD,&params1);

  /* Get the propagation time for the direct ray using the ftimeD function after we have gotten the value of the L parameter. */
  struct IceRayTracing::ftimeD_params params2a = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), -IceRayTracing::GetC(z0), IceRayTracing::c_light_ms,lvalueD};
  struct IceRayTracing::ftimeD_params params2b = {IceRayTracing::A_ice, IceRayTracing::GetB(z1), -IceRayTracing::GetC(z1), IceRayTracing::c_light_ms,lvalueD};
  struct IceRayTracing::ftimeD_params params2c = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary), IceRayTracing::c_light_ms, lvalueD};
  struct IceRayTracing::ftimeD_params params2d = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+1e-7), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary+1e-7), IceRayTracing::c_light_ms, lvalueD};

  /* we do the subtraction because we are measuring the time taken between the Tx and Rx positions */
  double timeD=0;
  double pathD=0;
  if(IceRayTracing::TransitionBoundary!=0){
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)>IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(-z1,&params2b);
      pathD=IceRayTracing::fpathD(-z0,&params2a) - IceRayTracing::fpathD(-z1,&params2b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(IceRayTracing::TransitionBoundary+1e-7,&params2d) + IceRayTracing::ftimeD(IceRayTracing::TransitionBoundary,&params2c) - IceRayTracing::ftimeD(-z1,&params2b);
      pathD=IceRayTracing::fpathD(-z0,&params2a) - IceRayTracing::fpathD(IceRayTracing::TransitionBoundary+1e-7,&params2d) + IceRayTracing::fpathD(IceRayTracing::TransitionBoundary,&params2c) - IceRayTracing::fpathD(-z1,&params2b);
    }
    if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(-z1,&params2b);
      pathD=IceRayTracing::fpathD(-z0,&params2a) - IceRayTracing::fpathD(-z1,&params2b);
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(-z1,&params2b);
      pathD=IceRayTracing::fpathD(-z0,&params2a) - IceRayTracing::fpathD(-z1,&params2b);
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(-z1,&params2b);
      pathD=IceRayTracing::fpathD(-z0,&params2a) - IceRayTracing::fpathD(-z1,&params2b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(IceRayTracing::TransitionBoundary+1e-7,&params2d) + IceRayTracing::ftimeD(IceRayTracing::TransitionBoundary,&params2c) - IceRayTracing::ftimeD(-z1,&params2b);
      pathD=IceRayTracing::fpathD(-z0,&params2a) - IceRayTracing::fpathD(IceRayTracing::TransitionBoundary+1e-7,&params2d) + IceRayTracing::fpathD(IceRayTracing::TransitionBoundary,&params2c) - IceRayTracing::fpathD(-z1,&params2b);
    }
  }else{
    timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(-z1,&params2b);
    pathD=IceRayTracing::fpathD(-z0,&params2a) - IceRayTracing::fpathD(-z1,&params2b);
  }

  /* Setup the function that will be used to calculate the angle of reception for all the rays */
  gsl_function F5;
  struct IceRayTracing::fDnfR_params params5a = {IceRayTracing::A_ice, IceRayTracing::GetB(z1), -IceRayTracing::GetC(z1), lvalueD};
  double result, abserr;
  F5.function = &IceRayTracing::fDnfR;

  /* Calculate the recieve angle for direc rays by calculating the derivative of the function at the Rx position */
  F5.params = &params5a;
  gsl_deriv_central (&F5, -z1, 1e-8, &result, &abserr);
  double RangD=atan(result)*(180.0/IceRayTracing::pi);

  /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
  if(z1==z0 && std::isnan(RangD)==true){
    RangD=180-LangD;
  }
  
  /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
  if(z1!=z0 && std::isnan(RangD)==true){
    RangD=90;
  }
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangD;
  output[1]=LangD;
  output[2]=timeD;
  output[3]=lvalueD;
  output[4]=checkzeroD;
  output[5]=pathD;
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangD;
    output[1]=180-RangD;
  }
  
  return output;
}

/* This functions works for the Reflected ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
double *IceRayTracing::GetReflectedRayPar(double z0, double x1 ,double z1){

  double *output=new double[11];

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the reflected ray. */
  gsl_function F3;
  struct IceRayTracing::fDanfRa_params params3= {IceRayTracing::A_ice, z0, x1, z1};
  F3.function = &IceRayTracing::fRa;
  F3.params = &params3;

  // gsl_function_fdf F3;
  // struct IceRayTracing::fDanfRa_params params3= {IceRayTracing::A_ice, z0, x1, z1};
  // F3.f = &IceRayTracing::fRa;
  // F3.df = &IceRayTracing::fRa_df;
  // F3.fdf = &IceRayTracing::fRa_fdf;
  // F3.params = &params3;
  
  /* In my raytracing solution given in the function IceRayTracing::fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0, z1 and also 1e-7 m and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */
  double UpLimnz[]={IceRayTracing::Getnz(z1),IceRayTracing::Getnz(z0),IceRayTracing::Getnz(1e-7)};
  double *UpperLimitL=min_element(UpLimnz,UpLimnz+3);
  
  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fRa function. */
  double lvalueR=IceRayTracing::FindFunctionRoot(F3,1e-7,UpperLimitL[0]);
  double LangR=asin(lvalueR/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
  double checkzeroR=IceRayTracing::fRa(lvalueR,&params3); 

  /* Get the propagation time for the reflected ray using the ftimeD function after we have gotten the value of the L parameter. */
  struct IceRayTracing::ftimeD_params params3a = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), IceRayTracing::c_light_ms,lvalueR};
  struct IceRayTracing::ftimeD_params params3b = {IceRayTracing::A_ice, IceRayTracing::GetB(z1), IceRayTracing::GetC(z1), IceRayTracing::c_light_ms,lvalueR};
  struct IceRayTracing::ftimeD_params params3c = {IceRayTracing::A_ice, IceRayTracing::GetB(1e-7), IceRayTracing::GetC(1e-7), IceRayTracing::c_light_ms,lvalueR};
  struct IceRayTracing::ftimeD_params params3d = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), IceRayTracing::GetC(IceRayTracing::TransitionBoundary), IceRayTracing::c_light_ms, lvalueR};
  struct IceRayTracing::ftimeD_params params3f = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+1e-7), IceRayTracing::GetC(IceRayTracing::TransitionBoundary+1e-7), IceRayTracing::c_light_ms, lvalueR};
  /* We do the subtraction because we are measuring the time taken between the Tx and Rx positions. In the reflected case we basically have two direct rays 1) from Tx to surface 2) from surface to Rx. Also get the time for the two individual direct rays separately */
  double timeR1=0;
  double timeR2=0;
  double pathR1=0;
  double pathR2=0;

  if(IceRayTracing::TransitionBoundary!=0){
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)>IceRayTracing::TransitionBoundary){
      timeR1=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(-IceRayTracing::TransitionBoundary,&params3d) + IceRayTracing::ftimeD(-(IceRayTracing::TransitionBoundary+1e-7),&params3f) - IceRayTracing::ftimeD(z0,&params3a);
      timeR2=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(-IceRayTracing::TransitionBoundary,&params3d) + IceRayTracing::ftimeD(-(IceRayTracing::TransitionBoundary+1e-7),&params3f) - IceRayTracing::ftimeD(z1,&params3b);

      pathR1=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(-IceRayTracing::TransitionBoundary,&params3d) + IceRayTracing::fpathD(-(IceRayTracing::TransitionBoundary+1e-7),&params3f) - IceRayTracing::fpathD(z0,&params3a);
      pathR2=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(-IceRayTracing::TransitionBoundary,&params3d) + IceRayTracing::fpathD(-(IceRayTracing::TransitionBoundary+1e-7),&params3f) - IceRayTracing::fpathD(z1,&params3b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      timeR1=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(-IceRayTracing::TransitionBoundary,&params3d) + IceRayTracing::ftimeD(-(IceRayTracing::TransitionBoundary+1e-7),&params3f) - IceRayTracing::ftimeD(z0,&params3a);
      timeR2=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(z1,&params3b);

      pathR1=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(-IceRayTracing::TransitionBoundary,&params3d) + IceRayTracing::fpathD(-(IceRayTracing::TransitionBoundary+1e-7),&params3f) - IceRayTracing::fpathD(z0,&params3a);
      pathR2=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(z1,&params3b);
    }
    if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      timeR1=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(z0,&params3a);
      timeR2=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(z1,&params3b);

      pathR1=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(z0,&params3a);
      pathR2=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(z1,&params3b);
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
      timeR1=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(z0,&params3a);
      timeR2=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(z1,&params3b);

      pathR1=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(z0,&params3a);
      pathR2=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(z1,&params3b);
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
      timeR1=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(z0,&params3a);
      timeR2=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(z1,&params3b);

      pathR1=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(z0,&params3a);
      pathR2=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(z1,&params3b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
      timeR1=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(-IceRayTracing::TransitionBoundary,&params3d) + IceRayTracing::ftimeD(-(IceRayTracing::TransitionBoundary+1e-7),&params3f) - IceRayTracing::ftimeD(z0,&params3a);
      timeR2=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(z1,&params3b);

      pathR1=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(-IceRayTracing::TransitionBoundary,&params3d) + IceRayTracing::fpathD(-(IceRayTracing::TransitionBoundary+1e-7),&params3f) - IceRayTracing::fpathD(z0,&params3a);
      pathR2=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(z1,&params3b);
    } 
  }else{
    timeR1=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(z0,&params3a);
    timeR2=IceRayTracing::ftimeD(-1e-7,&params3c) - IceRayTracing::ftimeD(z1,&params3b);

    pathR1=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(z0,&params3a);
    pathR2=IceRayTracing::fpathD(-1e-7,&params3c) - IceRayTracing::fpathD(z1,&params3b);
  }
  double timeR=timeR1 + timeR2;
  double pathR=pathR1 + pathR2;
  
  /* flip the times back if the original positions were flipped */
  if(Flip==true){
    double dumR=timeR2;
    timeR2=timeR1;
    timeR1=dumR;

    double dumRb=pathR2;
    pathR2=pathR1;
    pathR1=dumRb;
  }
  timeR1=timeR1;
  timeR2=timeR2;

  pathR1=pathR1;
  pathR2=pathR2;
  
  /* Setup the function that will be used to calculate the angle of reception for all the rays */
  gsl_function F5;
  struct IceRayTracing::fDnfR_params params5b = {IceRayTracing::A_ice, IceRayTracing::GetB(z1), IceRayTracing::GetC(z1), lvalueR};
  double result, abserr;
  F5.function = &IceRayTracing::fDnfR;
  
  /* Calculate the recieve angle for reflected ray by calculating the derivative of the function at the Rx position */
  F5.params = &params5b;
  gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
  double RangR=180-atan(result)*(180.0/IceRayTracing::pi);

  /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
  if(z1==z0 && std::isnan(RangR)==true){
    RangR=180-LangR;
  }

  /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
  if(z1!=z0 && std::isnan(RangR)==true){
    RangR=90;
  }

  if(std::isnan(checkzeroR)){
    checkzeroR=-1000;
  }
  
  /* Calculate the angle of incidence of the reflected ray at the surface ice. This will be used to calculate the Fresnel Coefficients. The angle is calculated by calculating the derivative of the ray path fucnction at the surface*/
  struct IceRayTracing::fDnfR_params paramsIAngB = {IceRayTracing::A_ice, IceRayTracing::GetB(1e-7), IceRayTracing::GetC(1e-7), lvalueR};
  F5.function = &IceRayTracing::fDnfR; 
  F5.params = &paramsIAngB;
  gsl_deriv_central (&F5, -1e-7, 1e-8, &result, &abserr);
  double IncidenceAngleInIce=atan(result)*(180.0/IceRayTracing::pi);
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangR;
  output[1]=LangR;
  output[2]=timeR;
  output[3]=lvalueR;
  output[4]=checkzeroR;
  output[5]=timeR1;
  output[6]=timeR2;
  output[7]=IncidenceAngleInIce;
  output[8]=pathR;
  output[9]=pathR1;
  output[10]=pathR2;
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangR;
    output[1]=180-RangR;
  } 
  
  return output;
}

/* This functions works for the Refracted ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. It requires the launch angle of the reflected ray as an input. */
double *IceRayTracing::GetRefractedRayPar(double z0, double x1 ,double z1, double LangR, double RangR, double checkzeroD, double checkzeroR){
  
  double *output=new double[22];

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }

  dsw=0;
  if(Flip==true){
    dsw=180-LangR;
    LangR=180-RangR;
    RangR=dsw;
  }
  
  /* Set up all the variables that will be used to get the parameters for the refracted ray */
  //double lvalueR=sin(LangR*(IceRayTracing::pi/180))*IceRayTracing::Getnz(z0);
  double lvalueRa[2]={0,0};
  double LangRa[2]={0,0};
  double checkzeroRa[2]={-1000,-1000};

  double timeRa[2]={0,0};
  double timeRa1[2]={0,0};
  double timeRa2[2]={0,0};

  double pathRa[2]={0,0};
  double pathRa1[2]={0,0};
  double pathRa2[2]={0,0};
  
  double raytime[2]={0,0};
  double RangRa[2]={0,0};
  double zmax[2]={10,10};
  
  /* In my raytracing solution given in the function IceRayTracing::fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0, z1 and also 1e-7 m and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */
  double UpLimnz[]={IceRayTracing::Getnz(z0),IceRayTracing::Getnz(z1)};
  double *UpperLimitL=min_element(UpLimnz,UpLimnz+2);
  
  /* First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the refracted ray. */
  gsl_function F4;
  struct IceRayTracing::fDanfRa_params params4= {IceRayTracing::A_ice, z0, x1, z1};
  F4.function = &IceRayTracing::fRaa;
  F4.params = &params4;

  gsl_function_fdf F4b;
  F4b.f = &IceRayTracing::fRaa;
  F4b.df = &IceRayTracing::fRaa_df;
  F4b.fdf = &IceRayTracing::fRaa_fdf;
  F4b.params = &params4;
  
  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fRaa function. The thing to note here is the lower limit of the minimisation function is set to the L value corresponding to the reflected ray launch angle. Since we know the refracted ray always has bigger launch angle the reflected ray this reduces our range and makes the function more efficient at finding the refracted ray launch angle. */
  double LowerLimit=0;
  LowerLimit=IceRayTracing::Getnz(z0)*sin((64.0*(IceRayTracing::pi/180.0)));
  if(LowerLimit>UpperLimitL[0]){
    LowerLimit=IceRayTracing::Getnz(z0)*sin((LangR*(IceRayTracing::pi/180.0)));
  }

  lvalueRa[0]=IceRayTracing::FindFunctionRoot(F4,LowerLimit,UpperLimitL[0]);
  LangRa[0]=asin(lvalueRa[0]/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
  checkzeroRa[0]=(IceRayTracing::fRaa(lvalueRa[0],&params4));
  zmax[0]=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa[0])+1e-7;

  if(fabs(checkzeroRa[0])>0.5){
    lvalueRa[0]=IceRayTracing::FindFunctionRootFDF(F4b,LowerLimit,UpperLimitL[0]);
    LangRa[0]=asin(lvalueRa[0]/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
    checkzeroRa[0]=IceRayTracing::fRaa(lvalueRa[0],&params4);
    zmax[0]=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa[0])+1e-7;
  }  

  if(fabs(checkzeroRa[0])<0.5 && fabs(checkzeroD)>0.5 && fabs(checkzeroR)>0.5){
    lvalueRa[1]=IceRayTracing::FindFunctionRoot(F4,lvalueRa[0]-0.23,lvalueRa[0]-0.023);
    LangRa[1]=asin(lvalueRa[1]/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
    checkzeroRa[1]=IceRayTracing::fRaa(lvalueRa[1],&params4);
    zmax[1]=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa[1])+1e-7;
  
    if(fabs(checkzeroRa[1])>0.5 || std::isnan(checkzeroRa[1])==true || fabs(lvalueRa[1]-lvalueRa[0])<1e-4 ){  
      lvalueRa[1]=IceRayTracing::FindFunctionRoot(F4,lvalueRa[0]-0.15,lvalueRa[0]-0.023);
      LangRa[1]=asin(lvalueRa[1]/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
      checkzeroRa[1]=IceRayTracing::fRaa(lvalueRa[1],&params4);
      zmax[1]=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa[1])+1e-7;
    }

    if(fabs(checkzeroRa[1])>0.5 || std::isnan(checkzeroRa[1])==true || fabs(lvalueRa[1]-lvalueRa[0])<1e-4 ){
      if(lvalueRa[0]+0.005<UpperLimitL[0]){
	lvalueRa[1]=IceRayTracing::FindFunctionRoot(F4,lvalueRa[0]+0.005,UpperLimitL[0]);
	LangRa[1]=asin(lvalueRa[1]/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
	checkzeroRa[1]=IceRayTracing::fRaa(lvalueRa[1],&params4);
	zmax[1]=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa[1])+1e-7;
      }else{
	lvalueRa[1]=IceRayTracing::FindFunctionRoot(F4,lvalueRa[0]-0.1,lvalueRa[0]-0.01);
	LangRa[1]=asin(lvalueRa[1]/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
	checkzeroRa[1]=IceRayTracing::fRaa(lvalueRa[1],&params4);
	zmax[1]=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa[1])+1e-7;
      }
    }

    if(fabs(checkzeroRa[1])>0.5 || std::isnan(checkzeroRa[1])==true || fabs(lvalueRa[1]-lvalueRa[0])<1e-4 ){
      double lvalueRatmp=IceRayTracing::FindFunctionRootFDF(F4b,lvalueRa[0]-0.23,lvalueRa[0]-0.023);
      if(fabs(lvalueRatmp)<IceRayTracing::A_ice){
	lvalueRa[1]=IceRayTracing::FindFunctionRootFDF(F4b,lvalueRa[0]-0.23,lvalueRa[0]-0.023);
	LangRa[1]=asin(lvalueRa[1]/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
	checkzeroRa[1]=IceRayTracing::fRaa(lvalueRa[1],&params4);
	zmax[1]=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa[1])+1e-7;
      }
    }

    if(fabs(checkzeroRa[1])>0.5 || std::isnan(checkzeroRa[1])==true || fabs(lvalueRa[1]-lvalueRa[0])<1e-4 ){
      double lvalueRatmp=IceRayTracing::FindFunctionRootFDF(F4b,lvalueRa[0]-0.1,lvalueRa[0]-0.023);
      if(fabs(lvalueRatmp)<IceRayTracing::A_ice){
	lvalueRa[1]=lvalueRatmp;
	LangRa[1]=asin(lvalueRa[1]/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
	checkzeroRa[1]=IceRayTracing::fRaa(lvalueRa[1],&params4);
	zmax[1]=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa[1])+1e-7;
      }
    }

    if(fabs(checkzeroRa[1])<0.5 && fabs(checkzeroRa[0])<0.5 && fabs(lvalueRa[1]-lvalueRa[0])<1e-4 ){
      checkzeroRa[1]=1000;
    }
    if(std::isnan(LangRa[0])==true){
      LangRa[0]=0;
    }
    if(std::isnan(LangRa[1])==true){
      LangRa[1]=0;
    }

    if(LangRa[1]<LangRa[0] && fabs(checkzeroRa[0])<0.5 && fabs(checkzeroRa[1])<0.5){
      swap(lvalueRa[1],lvalueRa[0]);
      swap(LangRa[1],LangRa[0]);
      swap(checkzeroRa[1],checkzeroRa[0]);
      swap(zmax[1],zmax[0]);
    }

  }else{
    lvalueRa[1]=0;
    LangRa[1]=0;
    checkzeroRa[1]=-1000;
    zmax[1]=-1000;
  }
 
  for(int i=0;i<2;i++){////loop over the two refracted solutions
  
    /* If we still did not find a refracted ray then set the check zero parameter to 1000 to make sure my code does not output this as a possible solution */
    if(std::isnan(checkzeroRa[i])==true){
      checkzeroRa[i]=-1000;
    }

    /* If the turning point depth also came out to be zero then now we are sure that there is no refracted ray */
    if(zmax[i]==1e-7 || zmax[i]<=0){
      checkzeroRa[i]=-1000;
    }

    /* Set parameters for ftimeD function to get the propagation time for the refracted ray */
    struct IceRayTracing::ftimeD_params params4a = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), IceRayTracing::c_light_ms,lvalueRa[i]};
    struct IceRayTracing::ftimeD_params params4b = {IceRayTracing::A_ice, IceRayTracing::GetB(z1), IceRayTracing::GetC(z1), IceRayTracing::c_light_ms,lvalueRa[i]};
    struct IceRayTracing::ftimeD_params params4c = {IceRayTracing::A_ice, IceRayTracing::GetB(zmax[i]), IceRayTracing::GetC(zmax[i]), IceRayTracing::c_light_ms,lvalueRa[i]};
    struct IceRayTracing::ftimeD_params params4d = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), IceRayTracing::GetC(IceRayTracing::TransitionBoundary), IceRayTracing::c_light_ms, lvalueRa[i]};
    struct IceRayTracing::ftimeD_params params4f = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+1e-7), IceRayTracing::GetC(IceRayTracing::TransitionBoundary+1e-7), IceRayTracing::c_light_ms, lvalueRa[i]};
    
    /* This if condition checks if the function has not gone crazy and given us a turning point of the ray which is lower than both Tx and Rx and is shallower in depth than both */
    if((z0<-zmax[i] || zmax[i]<-z1)){
      /* We do the subtraction because we are measuring the time taken between the Tx and Rx positions. In the refracted case we basically have two direct rays 1) from Tx to turning point 2) from turning point to Rx. Also get the time for the two individual direct rays separately */
     
      if(IceRayTracing::TransitionBoundary!=0){
	if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)>IceRayTracing::TransitionBoundary){
	  if(zmax[i]<=IceRayTracing::TransitionBoundary){
	    timeRa1[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(-IceRayTracing::TransitionBoundary,&params4d) + IceRayTracing::ftimeD(-(IceRayTracing::TransitionBoundary+1e-7),&params4f) - IceRayTracing::ftimeD(z0,&params4a);
	    timeRa2[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(-IceRayTracing::TransitionBoundary,&params4d) + IceRayTracing::ftimeD(-(IceRayTracing::TransitionBoundary+1e-7),&params4f) - IceRayTracing::ftimeD(z1,&params4b);

	    pathRa1[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(-IceRayTracing::TransitionBoundary,&params4d) + IceRayTracing::fpathD(-(IceRayTracing::TransitionBoundary+1e-7),&params4f) - IceRayTracing::fpathD(z0,&params4a);
	    pathRa2[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(-IceRayTracing::TransitionBoundary,&params4d) + IceRayTracing::fpathD(-(IceRayTracing::TransitionBoundary+1e-7),&params4f) - IceRayTracing::fpathD(z1,&params4b);
	    
	  }else{
	    timeRa1[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z0,&params4a);
	    timeRa2[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z1,&params4b);

	    pathRa1[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z0,&params4a);
	    pathRa2[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z1,&params4b);
	  }  
	}
	if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
	  timeRa1[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(-IceRayTracing::TransitionBoundary,&params4d) + IceRayTracing::ftimeD(-(IceRayTracing::TransitionBoundary+1e-7),&params4f) - IceRayTracing::ftimeD(z0,&params4a);
	  timeRa2[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z1,&params4b);

	  pathRa1[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(-IceRayTracing::TransitionBoundary,&params4d) + IceRayTracing::fpathD(-(IceRayTracing::TransitionBoundary+1e-7),&params4f) - IceRayTracing::fpathD(z0,&params4a);
	  pathRa2[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z1,&params4b);

	}
	if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
	  timeRa1[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z0,&params4a);
	  timeRa2[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z1,&params4b);

	  pathRa1[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z0,&params4a);
	  pathRa2[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z1,&params4b);
	} 
	if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)<IceRayTracing::TransitionBoundary){
	  timeRa1[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z0,&params4a);
	  timeRa2[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z1,&params4b);

	  pathRa1[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z0,&params4a);
	  pathRa2[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z1,&params4b);
	}
	if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
	  timeRa1[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z0,&params4a);
	  timeRa2[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z1,&params4b);

	  pathRa1[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z0,&params4a);
	  pathRa2[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z1,&params4b);
	}
	if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(z1)==IceRayTracing::TransitionBoundary){
	  timeRa1[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(-IceRayTracing::TransitionBoundary,&params4d) + IceRayTracing::ftimeD(-(IceRayTracing::TransitionBoundary+1e-7),&params4f) - IceRayTracing::ftimeD(z0,&params4a);
	  timeRa2[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z1,&params4b);

	  pathRa1[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(-IceRayTracing::TransitionBoundary,&params4d) + IceRayTracing::fpathD(-(IceRayTracing::TransitionBoundary+1e-7),&params4f) - IceRayTracing::fpathD(z0,&params4a);
	  pathRa2[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z1,&params4b);
	}
      }else{
	timeRa1[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z0,&params4a);
	timeRa2[i]=IceRayTracing::ftimeD(-zmax[i],&params4c) - IceRayTracing::ftimeD(z1,&params4b);

	pathRa1[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z0,&params4a);
	pathRa2[i]=IceRayTracing::fpathD(-zmax[i],&params4c) - IceRayTracing::fpathD(z1,&params4b);
      }
      
      raytime[i]= timeRa1[i] + timeRa2[i];
      pathRa[i]= pathRa1[i] + pathRa2[i];
      
      if(Flip==true){
	double dumRa=timeRa2[i];
	timeRa2[i]=timeRa1[i];
	timeRa1[i]=dumRa;

	dumRa=pathRa2[i];
	pathRa2[i]=pathRa1[i];
	pathRa1[i]=dumRa;
      }
    }
    timeRa[i]=raytime[i];

    /* Setup the function that will be used to calculate the angle of reception for all the rays */
    gsl_function F5;
    struct IceRayTracing::fDnfR_params params5c = {IceRayTracing::A_ice, IceRayTracing::GetB(z1), IceRayTracing::GetC(z1), lvalueRa[i]};
    double result, abserr;
    F5.function = &IceRayTracing::fDnfR;
    
    /* Calculate the recieve angle for refacted ray by calculating the derivative of the function at the Rx position */
    F5.params = &params5c;
    gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
    RangRa[i]=180-atan(result)*(180.0/IceRayTracing::pi);

    /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
    if(z1==z0 && std::isnan(RangRa[i])==true){
      RangRa[i]=180-LangRa[i];
    }

    /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
    if(z1!=z0 && std::isnan(RangRa[i])==true){
      RangRa[i]=90;
    }

  }//// i loop
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangRa[0];
  output[1]=LangRa[0];
  output[2]=timeRa[0];
  output[3]=lvalueRa[0];
  output[4]=checkzeroRa[0];
  output[5]=timeRa1[0];
  output[6]=timeRa2[0];
  output[7]=zmax[0];
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangRa[0];
    output[1]=180-RangRa[0];
  }

  output[8]=RangRa[1];
  output[9]=LangRa[1];
  output[10]=timeRa[1];
  output[11]=lvalueRa[1];
  output[12]=checkzeroRa[1];
  output[13]=timeRa1[1];
  output[14]=timeRa2[1];
  output[15]=zmax[1];

  output[16]=pathRa[0];
  output[17]=pathRa1[0];
  output[18]=pathRa2[0];

  output[19]=pathRa[1];
  output[20]=pathRa1[1];
  output[21]=pathRa2[1];
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[8]=180-LangRa[1];
    output[9]=180-RangRa[1];
  }
  
  return output;
}


/* This function returns the x and z values for the full Direct ray and prints out the ray path in a text file */
void IceRayTracing::GetFullDirectRayPath(double z0, double x1, double z1,double lvalueD, vector <double> &x, vector <double> &z){
  
  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
   
  /* Set the name of the text files */
  //ofstream aoutD("DirectRay.txt");

  /* Set the step size for plotting */
  double h=0.5;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;

  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;
  struct IceRayTracing::fDnfR_params params6a;
  struct IceRayTracing::fDnfR_params params6b;
  struct IceRayTracing::fDnfR_params params6c;
  struct IceRayTracing::fDnfR_params params6d;

  
  for(int i=0;i<dmax;i++){
    params6a = {IceRayTracing::A_ice, IceRayTracing::GetB(zn), IceRayTracing::GetC(zn), lvalueD};
    params6b = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), lvalueD};
    params6c = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), IceRayTracing::GetC(IceRayTracing::TransitionBoundary), lvalueD};
    params6d = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+1e-7), IceRayTracing::GetC(IceRayTracing::TransitionBoundary+1e-7), lvalueD};

    if(IceRayTracing::TransitionBoundary!=0){
      if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(-IceRayTracing::TransitionBoundary,&params6c) + IceRayTracing::fDnfR(-(IceRayTracing::TransitionBoundary+1e-7),&params6d) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)>IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(zn)==IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)==IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(-IceRayTracing::TransitionBoundary,&params6c) + IceRayTracing::fDnfR(-(IceRayTracing::TransitionBoundary+1e-7),&params6d) - IceRayTracing::fDnfR(z0,&params6b);
      }
    }else{
      xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
    }
    
    checknan=IceRayTracing::fDnfR(zn,&params6a);
    if(std::isnan(checknan)==false && Flip==false){
      x.push_back(xn);
      z.push_back(zn);
      //aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
      npnt++;
    }

    if(std::isnan(checknan)==false && Flip==true){
      x.push_back(x1-xn);
      z.push_back(zn);
      //aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;      
    }  
  }

  params6a = {IceRayTracing::A_ice, IceRayTracing::GetB(zn), IceRayTracing::GetC(zn), lvalueD};
  params6b = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), lvalueD};
  xn=IceRayTracing::fDnfR(zn,&params6a)-IceRayTracing::fDnfR(z0,&params6b);  
  if(Flip==true){
    x.push_back(x1-xn);
    z.push_back(zn);
    //aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
  }else{
    x.push_back(xn);
    z.push_back(zn);
    //aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
  }
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }

}

/* This function returns the x and z values for the full Reflected ray path and prints out the ray path in a text file */
void IceRayTracing::GetFullReflectedRayPath(double z0, double x1, double z1,double lvalueR, vector <double> &x, vector <double> &z){
  
  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* Set the name of the text files */
  //ofstream aoutR("ReflectedRay.txt");
  /* Set the step size for plotting. */
  double h=0.5;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;
  
  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;  
  struct IceRayTracing::fDnfR_params params6a;
  struct IceRayTracing::fDnfR_params params6b;
  struct IceRayTracing::fDnfR_params params6c;
  struct IceRayTracing::fDnfR_params params6f;
  struct IceRayTracing::fDnfR_params params6d;

  /* Map out the 1st part of the reflected ray */
  for(int i=0;i<dmax;i++){
    params6a = {IceRayTracing::A_ice, IceRayTracing::GetB(zn), -IceRayTracing::GetC(zn), lvalueR};
    params6b = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), -IceRayTracing::GetC(z0), lvalueR};
    params6c = {IceRayTracing::A_ice, IceRayTracing::GetB(1e-7), -IceRayTracing::GetC(1e-7), lvalueR};
    params6d = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary), lvalueR};
    params6f = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+1e-7), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary+1e-7), lvalueR};

    double distancez0z1=0;
    double distancez0surface=0;  
    if(IceRayTracing::TransitionBoundary!=0){
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)>IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(1e-7,&params6c) - IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary,&params6d) + IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary+1e-7,&params6f) - IceRayTracing::fDnfR(-z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary,&params6d) + IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary+1e-7,&params6f) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(1e-7,&params6c) - IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary,&params6d) + IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary+1e-7,&params6f) - IceRayTracing::fDnfR(-z0,&params6b);
      }
      if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(1e-7,&params6c) - IceRayTracing::fDnfR(-z0,&params6b); 
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(1e-7,&params6c) - IceRayTracing::fDnfR(-z0,&params6b); 
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(zn)==IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(1e-7,&params6c) - IceRayTracing::fDnfR(-z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)==IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary,&params6d) + IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary+1e-7,&params6f) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(1e-7,&params6c) - IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary,&params6d) + IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary+1e-7,&params6f) - IceRayTracing::fDnfR(-z0,&params6b);
      }
    }else{
      distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
      distancez0surface=IceRayTracing::fDnfR(1e-7,&params6c) - IceRayTracing::fDnfR(-z0,&params6b);
    }

    xn= distancez0z1 - 2*(distancez0surface);
 
    checknan=IceRayTracing::fDnfR(-zn,&params6a);
    if(std::isnan(checknan)==false && zn<=0 && Flip==false){
      x.push_back(xn);
      z.push_back(zn);
      //aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(std::isnan(checknan)==false && zn<=0 && Flip==true){
      x.push_back(x1-xn);
      z.push_back(zn);
      //aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
      
    zn=zn+h;
    if(zn>0){
      i=dmax+2;      
    }
  }
  
  /* Map out the 2nd part of the reflected ray */
  zn=-1e-7;
  for(int i=0;i<dmax;i++){  
    params6a = {IceRayTracing::A_ice, IceRayTracing::GetB(zn), IceRayTracing::GetC(zn), lvalueR};
    params6b = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), lvalueR};
    params6c = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), IceRayTracing::GetC(IceRayTracing::TransitionBoundary), lvalueR};
    params6d = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+1e-7), IceRayTracing::GetC(IceRayTracing::TransitionBoundary+1e-7), lvalueR};

    if(IceRayTracing::TransitionBoundary!=0){
      if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(-IceRayTracing::TransitionBoundary,&params6c) + IceRayTracing::fDnfR(-(IceRayTracing::TransitionBoundary+1e-7),&params6d) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)>IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(zn)==IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)==IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(-IceRayTracing::TransitionBoundary,&params6c) + IceRayTracing::fDnfR(-(IceRayTracing::TransitionBoundary+1e-6),&params6d) - IceRayTracing::fDnfR(z0,&params6b);
      }
    }else{
      xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
    }
    
    checknan=IceRayTracing::fDnfR(zn,&params6a);
    if(std::isnan(checknan)==false && Flip==false){
      x.push_back(xn);
      z.push_back(zn);
      //aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }
      
    if(std::isnan(checknan)==false && Flip==true){
      x.push_back(x1-xn);
      z.push_back(zn); 
      //aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }

  params6a = {IceRayTracing::A_ice, IceRayTracing::GetB(zn), IceRayTracing::GetC(zn), lvalueR};
  params6b = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), lvalueR};
  xn=IceRayTracing::fDnfR(zn,&params6a) -IceRayTracing::fDnfR(z0,&params6b);
  if(Flip==true){
    x.push_back(x1-xn);
    z.push_back(zn); 
    //aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
  }else{
    x.push_back(xn);
    z.push_back(zn); 
    //aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
  }
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }

}

/* This function returns the x and z values for the full Refracted ray and prints out the ray path in a text file */
void IceRayTracing::GetFullRefractedRayPath(double z0, double x1, double z1, double zmax, double lvalueRa, vector <double> &x, vector <double> &z,int raynumber){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* Set the name of the text files */
  //ofstream aoutRa;
  /* Set the name of the text files */
  // if(raynumber==1){
  //   aoutRa.open("RefractedRay1.txt");
  // }
  // if(raynumber==2){
  //   aoutRa.open("RefractedRay2.txt");
  // }
  /* Set the step size for plotting. */
  double h=0.5;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;

  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;
  struct IceRayTracing::fDnfR_params params6a;
  struct IceRayTracing::fDnfR_params params6b;
  struct IceRayTracing::fDnfR_params params6c;  
  struct IceRayTracing::fDnfR_params params6f;
  struct IceRayTracing::fDnfR_params params6d;

  /* Map out the 1st part of the refracted ray */
  for(int i=0;i<dmax;i++){
    params6a = {IceRayTracing::A_ice, IceRayTracing::GetB(zn), -IceRayTracing::GetC(zn), lvalueRa};
    params6b = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), -IceRayTracing::GetC(z0), lvalueRa};
    params6c = {IceRayTracing::A_ice, IceRayTracing::GetB(zmax), -IceRayTracing::GetC(zmax), lvalueRa};
    params6d = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary), lvalueRa};
    params6f = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+1e-6), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary+1e-6), lvalueRa};

    double distancez0z1=0;
    double distancez0surface=0;  
    if(IceRayTracing::TransitionBoundary!=0){
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)>IceRayTracing::TransitionBoundary){
	if(zmax<=IceRayTracing::TransitionBoundary){
	  distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
	  distancez0surface=IceRayTracing::fDnfR(zmax,&params6c) - IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary,&params6d) + IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary+1e-6,&params6f) - IceRayTracing::fDnfR(-z0,&params6b);
	}else{
	  distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
	  distancez0surface=IceRayTracing::fDnfR(zmax,&params6c) - IceRayTracing::fDnfR(-z0,&params6b);
	}
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary,&params6d) + IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary+1e-6,&params6f) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(zmax,&params6c) - IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary,&params6d) + IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary+1e-6,&params6f) - IceRayTracing::fDnfR(-z0,&params6b);
      }
      if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(zmax,&params6c) - IceRayTracing::fDnfR(-z0,&params6b);
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(zmax,&params6c) - IceRayTracing::fDnfR(-z0,&params6b); 
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(zn)==IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(zmax,&params6c) - IceRayTracing::fDnfR(-z0,&params6b); 
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)==IceRayTracing::TransitionBoundary){
	distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary,&params6d) + IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary+1e-6,&params6f) - IceRayTracing::fDnfR(-z0,&params6b);
	distancez0surface=IceRayTracing::fDnfR(zmax,&params6c) - IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary,&params6d) + IceRayTracing::fDnfR(IceRayTracing::TransitionBoundary+1e-6,&params6f) - IceRayTracing::fDnfR(-z0,&params6b); 
      }
    }else{
      distancez0z1=IceRayTracing::fDnfR(-zn,&params6a) - IceRayTracing::fDnfR(-z0,&params6b);
      distancez0surface=IceRayTracing::fDnfR(zmax,&params6c) - IceRayTracing::fDnfR(-z0,&params6b);
    }

    xn= distancez0z1 - 2*(distancez0surface);
    
    checknan=IceRayTracing::fDnfR(-zn,&params6a);
    if(std::isnan(checknan)==false && zn<=0 && Flip==false){
      x.push_back(xn);
      z.push_back(zn);
      //aoutRa<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(std::isnan(checknan)==false && zn<=0 && Flip==true){
      x.push_back(x1-xn);
      z.push_back(zn);
      //aoutRa<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
    
    zn=zn+h;
    if(zn>-zmax){
      i=dmax+2;      
    }
  }

  /* Map out the 2nd part of the refracted ray */
  zn=-zmax;
  for(int i=0;i<dmax;i++){  
    params6a = {IceRayTracing::A_ice, IceRayTracing::GetB(zn), IceRayTracing::GetC(zn), lvalueRa};
    params6b = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), lvalueRa};
    params6c = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), IceRayTracing::GetC(IceRayTracing::TransitionBoundary), lvalueRa};
    params6d = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+1e-6), IceRayTracing::GetC(IceRayTracing::TransitionBoundary+1e-6), lvalueRa}; 
   
    if(IceRayTracing::TransitionBoundary!=0){
      if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(-IceRayTracing::TransitionBoundary,&params6c) + IceRayTracing::fDnfR(-(IceRayTracing::TransitionBoundary+1e-6),&params6d) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)>IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(zn)<IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(zn)==IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
      }
      if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(zn)==IceRayTracing::TransitionBoundary){
	xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(-IceRayTracing::TransitionBoundary,&params6c) + IceRayTracing::fDnfR(-(IceRayTracing::TransitionBoundary+1e-6),&params6d) - IceRayTracing::fDnfR(z0,&params6b);
      }
    }else{
      xn=IceRayTracing::fDnfR(zn,&params6a) - IceRayTracing::fDnfR(z0,&params6b);
    }
    
    checknan=IceRayTracing::fDnfR(zn,&params6a);
    if(std::isnan(checknan)==false && Flip==false){
      x.push_back(xn);
      z.push_back(zn);
      //aoutRa<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(std::isnan(checknan)==false && Flip==true){
      x.push_back(x1-xn);
      z.push_back(zn);
      //aoutRa<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
    
    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }

  params6a = {IceRayTracing::A_ice, IceRayTracing::GetB(zn), IceRayTracing::GetC(zn), lvalueRa};
  params6b = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), lvalueRa};  
  xn=IceRayTracing::fDnfR(zn,&params6a)-IceRayTracing::fDnfR(z0,&params6b);
  if(Flip==true){
    x.push_back(x1-xn);
    z.push_back(zn);
    //aoutRa<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
  }else{
    x.push_back(xn);
    z.push_back(zn);
    //aoutRa<<npnt<<" "<<xn<<" "<<zn<<endl;
  }  
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
}

/* function for plotting and storing all the rays */
void IceRayTracing::PlotAndStoreRays(double x0,double z0, double z1, double x1, double zmax[2], double lvalues[4], double checkzeroes[4]){

  
  double lvalueD=lvalues[0];
  double lvalueR=lvalues[1];
  double lvalueRa[2]={lvalues[2],lvalues[3]};

  double checkzeroD=checkzeroes[0];
  double checkzeroR=checkzeroes[1];
  double checkzeroRa[2]={checkzeroes[2],checkzeroes[3]}; 

  vector <double> xD,zD;
  vector <double> xR,zR;
  vector <double> xRa[2],zRa[2];
  GetFullDirectRayPath(z0,x1,z1,lvalueD,xD,zD);
  GetFullReflectedRayPath(z0,x1,z1,lvalueR,xR,zR);
  if((fabs(checkzeroR)>0.5 || fabs(checkzeroD)>0.5) && fabs(checkzeroRa[0])<0.5){
    GetFullRefractedRayPath(z0,x1,z1,zmax[0],lvalueRa[0],xRa[0],zRa[0],1);
    if(fabs(checkzeroRa[1])<0.5){
      GetFullRefractedRayPath(z0,x1,z1,zmax[1],lvalueRa[1],xRa[1],zRa[1],2);
    }
  }
  
  /* An example of how the function returns the filled in vector and you can see the printed out values as a sanity check */
  // for (int i=0;i<xD.size();i++){
  //   cout<<i<<" "<<xD[i]<<" "<<zD[i]<<endl;
  // }
  
}

double *IceRayTracing::IceRayTracing(double x0, double z0, double x1, double z1){

  /* define a pointer to give back the output of raytracing */ 
  double *output=new double[29];

  /* Store the ray paths in text files */
  bool PlotRayPaths=false;
  /* calculate the attenuation (not included yet!) */
  bool attcal=false;
  
  double Txcor[2]={x0,z0};/* Tx positions */
  double Rxcor[2]={x1,z1};/* Rx Positions */
  
  /*  ********This part of the code will try to get the Direct ray between Rx and Tx.********** */
  double* GetDirectRay=IceRayTracing::GetDirectRayPar(z0,x1,z1);
  double RangD=GetDirectRay[0];
  double LangD=GetDirectRay[1];
  double timeD=GetDirectRay[2];
  double lvalueD=GetDirectRay[3];
  double checkzeroD=GetDirectRay[4];
  double pathD=GetDirectRay[5];
  delete []GetDirectRay;
  
  /* ********This part of the code will try to get the Reflected ray between Rx and Tx.********** */
  double* GetReflectedRay=IceRayTracing::GetReflectedRayPar(z0,x1,z1);
  double RangR=GetReflectedRay[0];
  double LangR=GetReflectedRay[1];
  double timeR=GetReflectedRay[2];
  double lvalueR=GetReflectedRay[3];
  double checkzeroR=GetReflectedRay[4];
  double timeR1=GetReflectedRay[5];
  double timeR2=GetReflectedRay[6];
  double AngleOfIncidenceInIce=GetReflectedRay[7];
  double pathR=GetReflectedRay[8];
  double pathR1=GetReflectedRay[9];
  double pathR2=GetReflectedRay[10];
  
  delete []GetReflectedRay;

  /* ********This part of the code will try to get the Refracted ray between Rx and Tx.********** */
  double RangRa[2]={0,0};
  double LangRa[2]={0,0};
  double timeRa[2]={0,0};
  double lvalueRa[2]={0,0}; 
  double checkzeroRa[2]={-1000,-1000};
  double timeRa1[2]={0,0};
  double timeRa2[2]={0,0};
  double zmax[2]={0,0};

  double pathRa[2]={0,0};
  double pathRa1[2]={0,0};
  double pathRa2[2]={0,0};
  
  /* This if condition makes sure that we only try to find a refracted ray if we don't get two possible ray paths from the direct and reflected case. This saves us alot of time since we know that between each Tx and Rx position we only expect 2 rays. */
  if(fabs(checkzeroR)>0.5 || fabs(checkzeroD)>0.5){
    double* GetRefractedRay=IceRayTracing::GetRefractedRayPar(z0,x1,z1,LangR,RangR,checkzeroD,checkzeroR);
    RangRa[0]=GetRefractedRay[0];
    LangRa[0]=GetRefractedRay[1];
    timeRa[0]=GetRefractedRay[2];
    lvalueRa[0]=GetRefractedRay[3]; 
    checkzeroRa[0]=GetRefractedRay[4];
    timeRa1[0]=GetRefractedRay[5];
    timeRa2[0]=GetRefractedRay[6];
    zmax[0]=GetRefractedRay[7];

    if(fabs(checkzeroR)>0.5 && fabs(checkzeroD)>0.5){
      RangRa[1]=GetRefractedRay[8];
      LangRa[1]=GetRefractedRay[9];
      timeRa[1]=GetRefractedRay[10];
      lvalueRa[1]=GetRefractedRay[11]; 
      checkzeroRa[1]=GetRefractedRay[12];
      timeRa1[1]=GetRefractedRay[13];
      timeRa2[1]=GetRefractedRay[14];
      zmax[1]=GetRefractedRay[15];
    }

    pathRa[0]=GetRefractedRay[16];
    pathRa1[0]=GetRefractedRay[17];
    pathRa2[0]=GetRefractedRay[18];

    pathRa[1]=GetRefractedRay[19];
    pathRa1[1]=GetRefractedRay[20];
    pathRa2[1]=GetRefractedRay[21];
    
    delete []GetRefractedRay;
  }
  
  /* This part of the code can be used if the user wants to plot the individual ray paths. This part of the code prints out the individual ray paths in text files and also plots them on a canvas */
  if(PlotRayPaths==true){
    double lvalues[4];
    lvalues[0]=lvalueD;
    lvalues[1]=lvalueR;
    lvalues[2]=lvalueRa[0];
    lvalues[3]=lvalueRa[1];

    double checkzeroes[4];
    checkzeroes[0]=checkzeroD;
    checkzeroes[1]=checkzeroR;
    checkzeroes[2]=checkzeroRa[0];
    checkzeroes[3]=checkzeroRa[1];
    
    IceRayTracing::PlotAndStoreRays(x0,z0,z1,x1,zmax,lvalues,checkzeroes);
  }
  
  /* print out all the output from the code */
  //cout<<0<<" ,x0= "<<x0<<" ,z0= "<<z0<<" ,x1= "<<x1<<" ,z1= "<<z1<<" ,langRa= "<<LangRa[0]<<" ,langR= "<<LangR<<" ,langD= "<<LangD<<" ,langD-langR= "<<LangD-LangR<<" ,langD-langRa= "<<LangD-LangRa[0]<<" ,RangRa= "<<RangRa[0]<<" ,RangR= "<<RangR<<" ,RangD= "<<RangD<<" ,RangR-RangD= "<<RangR-RangD<<" ,RangRa-RangD= "<<RangRa[0]-RangD<<" ,timeRa= "<<timeRa[0]<<" ,timeR= "<<timeRa[0]<<" ,timeD= "<<timeD<<" ,timeR-timeD= "<<timeR-timeD<<" ,timeRa-timeD= "<<timeRa[0]-timeD<<" ,lvalueRa "<<lvalueRa[0]<<" ,lvalueR "<<lvalueR<<" "<<" ,lvalueD "<<lvalueD<<" ,checkzeroRa "<<checkzeroRa[0]<<" ,checkzeroR "<<checkzeroR<<" ,checkzeroD "<<checkzeroD<<endl;

  //cout<<0<<" ,x0= "<<x0<<" ,z0= "<<z0<<" ,x1= "<<x1<<" ,z1= "<<z1<<" ,langRa= "<<LangRa[1]<<" ,langR= "<<LangR<<" ,langD= "<<LangD<<" ,langD-langR= "<<LangD-LangR<<" ,langD-langRa= "<<LangD-LangRa[1]<<" ,RangRa= "<<RangRa[1]<<" ,RangR= "<<RangR<<" ,RangD= "<<RangD<<" ,RangR-RangD= "<<RangR-RangD<<" ,RangRa-RangD= "<<RangRa[1]-RangD<<" ,timeRa= "<<timeRa[1]<<" ,timeR= "<<timeRa[1]<<" ,timeD= "<<timeD<<" ,timeR-timeD= "<<timeR-timeD<<" ,timeRa-timeD= "<<timeRa[1]-timeD<<" ,lvalueRa "<<lvalueRa[1]<<" ,lvalueR "<<lvalueR<<" "<<" ,lvalueD "<<lvalueD<<" ,checkzeroRa "<<checkzeroRa[1]<<" ,checkzeroR "<<checkzeroR<<" ,checkzeroD "<<checkzeroD<<endl;

  /* Fill in the output pointer after calculating all the results */
  output[0]=LangD;
  output[1]=LangR;
  output[2]=LangRa[0];
  output[3]=LangRa[1];
  output[4]=timeD;
  output[5]=timeR;
  output[6]=timeRa[0];
  output[7]=timeRa[1];
  output[8]=RangD;
  output[9]=RangR;
  output[10]=RangRa[0];
  output[11]=RangRa[1];
  
  /* fill in the output array part where you fill in the times for the two parts of the reflected or refracted rays */
  if(fabs(checkzeroR)<0.5){
    output[12]=timeR1;
    output[13]=timeR2;
  }
  
  if(fabs(checkzeroRa[0])<0.5){
    output[14]=timeRa1[0];
    output[15]=timeRa2[0];
  }

  if(fabs(checkzeroRa[1])<0.5){
    output[16]=timeRa1[1];
    output[17]=timeRa2[1];
  }
  
  output[18]=AngleOfIncidenceInIce;
  output[19]=lvalueD;
  output[20]=lvalueR;
  output[21]=lvalueRa[0];
  output[22]=lvalueRa[1];
  output[23]=zmax[0];  
  output[24]=zmax[1];

  output[25]=pathD;
  output[26]=pathR;
  output[27]=pathRa[0];
  output[28]=pathRa[1];
  
  /* Set the recieve angle to be zero for a ray which did not give us a possible path between Tx and Rx. I use this as a flag to determine which two rays gave me possible ray paths. */
  if(fabs(checkzeroD)>0.5){
    output[8]=-1000;
  }
  if(fabs(checkzeroR)>0.5){
    output[9]=-1000;
  }
  if(fabs(checkzeroRa[0])>0.5){
    output[10]=-1000;
  }
  if(fabs(checkzeroRa[1])>0.5){
    output[11]=-1000;
  } 
  
  return output;
}

/* Analytical solution describing ray paths in ice as function of depth for constant refractive index*/
double IceRayTracing::fDnfR_Cnz(double x,void *params){
  
  struct IceRayTracing::fDnfR_params *p= (struct IceRayTracing::fDnfR_params *) params;
  double A = p->a;
  double L = p->l;
  
  return (L/sqrt(A*A-L*L))*x;
}

/* Analytical solution describing the ray path in ice as a function of the L parameter for constant refractive index*/
double IceRayTracing::fDnfR_L_Cnz(double x,void *params){
  
  struct IceRayTracing::fDnfR_L_params *p= (struct IceRayTracing::fDnfR_L_params *) params;
  double A = p->a;
  double Z = p->z;
  
  double out=0;
  if(A>x){
    out=(x/sqrt(A*A-x*x))*Z;
  }else{
    out=tan(asin(x/A))*Z;
  }
  return out;
}

/* This function is minimised to find the launch angle (or the L parameter) for the reflected ray for constant refractive index*/
double IceRayTracing::fRa_Cnz(double x,void *params){
  struct IceRayTracing::fDanfRa_params *p= (struct IceRayTracing::fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  struct IceRayTracing::fDnfR_L_params params1a = {A, 0, 0, -z1};
  struct IceRayTracing::fDnfR_L_params params1b = {A, 0, 0, -z0};
  struct IceRayTracing::fDnfR_L_params params1c = {A, 0, 0, 0.0};

  return IceRayTracing::fDnfR_L_Cnz(x,&params1a) - IceRayTracing::fDnfR_L_Cnz(x,&params1b) - 2*( IceRayTracing::fDnfR_L_Cnz(x,&params1c) - IceRayTracing::fDnfR_L_Cnz(x,&params1b) ) - x1;
}

double IceRayTracing::fRa_Cnz_df(double x,void *params){
  gsl_function F;
  F.function = &IceRayTracing::fRa_Cnz;
  F.params = params;
 
  double result,abserr;
  gsl_deriv_central (&F, x, 1e-8, &result, &abserr);

  return result;
}

void IceRayTracing::fRa_Cnz_fdf (double x, void *params,double *y, double *dy){ 
  *y = IceRayTracing::fRa_Cnz(x,params);
  *dy = IceRayTracing::fRa_Cnz_df(x,params);
}

/* This functions works for the Direct ray and gives you back the launch angle, receive angle and propagation time of the ray together with values of the L parameter. This for constant refractive index*/
double* IceRayTracing::GetDirectRayPar_Cnz(double z0, double x1, double z1, double A_ice_Cnz){

  double *output=new double[4];
  
  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* Calculate the launch angle and the value of the L parameter */
  double LangD=(IceRayTracing::pi*0.5-atan(fabs(z1-z0)/x1))*(180.0/IceRayTracing::pi);
  double lvalueD=A_ice_Cnz*sin(LangD*(IceRayTracing::pi/180.0));
  double timeD=(sqrt( pow(x1,2) + pow(z1-z0,2) )/IceRayTracing::c_light_ms)*A_ice_Cnz;
  
  /* Calculate the recieve angle for direct rays by which is the same as the launch angle */
  double RangD=LangD;
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangD;
  output[1]=LangD;
  output[2]=timeD;
  output[3]=lvalueD;

  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangD;
    output[1]=180-RangD;
  }
  
  return output;
}

/* This functions works for the Reflected ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter. This is for constant refractive index*/
double *IceRayTracing::GetReflectedRayPar_Cnz(double z0, double x1 , double z1, double A_ice_Cnz){

  double *output=new double[8];

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the reflected ray. */
  gsl_function F3;
  struct IceRayTracing::fDanfRa_params params3= {A_ice_Cnz, z0, x1, z1};
  F3.function = &IceRayTracing::fRa_Cnz;
  F3.params = &params3;

  // gsl_function_fdf F3;
  // struct IceRayTracing::fDanfRa_params params3= {A_ice_Cnz, z0, x1, z1};
  // F3.f = &fRa_Cnz;
  // F3.df = &fRa_Cnz_df;
  // F3.fdf = &fRa_Cnz_fdf;
  // F3.params = &params3;

  /* In my raytracing solution given in the function IceRayTracing::fDnfR_Cnz the launch angle (or the L parameter) has limit placed on it by this part in the solution that L<A . This sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit to be the angle of the direct ray as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */
  double UpperLimitL=A_ice_Cnz*sin(IceRayTracing::pi*0.5-atan(fabs(z1-z0)/x1));

  /* Do the minimisation and get the value of the L parameter and the launch angle */
  double lvalueR=IceRayTracing::FindFunctionRoot(F3,0.0,UpperLimitL);
  double LangR=asin(lvalueR/A_ice_Cnz)*(180.0/IceRayTracing::pi);
  
  /* In the reflected case we basically have two direct rays 1) from Tx to surface 2) from surface to Rx. . Also get the time for the two individual direct rays separately */
  double z2=0,x2=fabs(z0)*tan(LangR*(IceRayTracing::pi/180));///coordinates of point of incidence in ice at the surface
  double timeR1=(sqrt( pow(x2,2) + pow(z2-z0,2) )/IceRayTracing::c_light_ms)*A_ice_Cnz;
  double timeR2=(sqrt( pow(x2-x1,2) + pow(z2-z1,2) )/IceRayTracing::c_light_ms)*A_ice_Cnz;
  double timeR= timeR1 + timeR2;
  
  /* flip the times back if the original positions were flipped */
  if(Flip==true){
    double dumR=timeR2;
    timeR2=timeR1;
    timeR1=dumR;
  }
  timeR1=timeR1;
  timeR2=timeR2;
  
  /* Calculate the recieve angle for reflected ray using simple geometry*/
  double RangR=180-LangR;

  /* Calculate the angle of incidence of the reflected ray at the surface ice. This will be used to calculate the Fresnel Coefficients.*/
  double IncidenceAngleInIce=LangR;
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangR;
  output[1]=LangR;
  output[2]=timeR;
  output[3]=lvalueR;
  output[4]=0;
  output[5]=timeR1;
  output[6]=timeR2;
  output[7]=IncidenceAngleInIce;
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangR;
    output[1]=180-RangR;
  } 
  
  return output;
}

/* This function returns the x and z values for the full Direct ray and prints out the ray path in a text file. This is for a constant refractive index. */
void IceRayTracing::GetFullDirectRayPath_Cnz(double z0, double x1, double z1, double lvalueD, double A_ice_Cnz, vector <double> &x, vector <double> &z){
  
  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
   
  /* Set the name of the text files */
  //ofstream aoutD("DirectRay_Cnz.txt");
  /* Set the step size for plotting */
  double h=0.5;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;

  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;
  struct IceRayTracing::fDnfR_params params6a;
  struct IceRayTracing::fDnfR_params params6b;
  
  for(int i=0;i<dmax;i++){
    params6a = {A_ice_Cnz, IceRayTracing::GetB(zn), IceRayTracing::GetC(zn), lvalueD};
    params6b = {A_ice_Cnz, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), lvalueD};
    xn=IceRayTracing::fDnfR_Cnz(zn,&params6a)-IceRayTracing::fDnfR_Cnz(z0,&params6b);
    checknan=IceRayTracing::fDnfR(zn,&params6a);
    if(std::isnan(checknan)==false && Flip==false){
      x.push_back(xn);
      z.push_back(zn);
      //aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
      npnt++;
    }

    if(std::isnan(checknan)==false && Flip==true){
      x.push_back(x1-xn);
      z.push_back(zn);
      //aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;      
    }  
  }

  params6a = {A_ice_Cnz, 0, 0, lvalueD};
  params6b = {A_ice_Cnz, 0, 0, lvalueD};
  xn=IceRayTracing::fDnfR_Cnz(zn,&params6a)-IceRayTracing::fDnfR_Cnz(z0,&params6b); 
  if(Flip==true){
    x.push_back(x1-xn);
    z.push_back(zn);
    //aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
  }else{
    x.push_back(xn);
    z.push_back(zn);
    //aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
  }
  npnt++;
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
}

/* This function returns the x and z values for the full Reflected ray and prints out the ray path in a text file. This is for a constant refractive index. */
void IceRayTracing::GetFullReflectedRayPath_Cnz(double z0, double x1, double z1, double lvalueR, double A_ice_Cnz,vector <double> &x, vector <double> &z){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  // /* Set the name of the text files */
  //ofstream aoutR("ReflectedRay_Cnz.txt");
  /* Set the step size for plotting. */
  double h=0.5;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;
  
  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;  
  struct IceRayTracing::fDnfR_params params6a;
  struct IceRayTracing::fDnfR_params params6b;
  struct IceRayTracing::fDnfR_params params6c;

  /* Map out the 1st part of the reflected ray */
  for(int i=0;i<dmax;i++){
    params6a = {A_ice_Cnz, 0, 0, lvalueR};
    params6b = {A_ice_Cnz, 0, 0, lvalueR};
    params6c = {A_ice_Cnz, 0, 0, lvalueR};
    xn=(IceRayTracing::fDnfR_Cnz(-zn,&params6a)-IceRayTracing::fDnfR_Cnz(-z0,&params6b)+2*fabs(IceRayTracing::fDnfR_Cnz(0.0,&params6c)-IceRayTracing::fDnfR_Cnz(-z0,&params6b)));
    checknan=IceRayTracing::fDnfR_Cnz(-zn,&params6a);
    if(std::isnan(checknan)==false && zn<=0 && Flip==false){
      x.push_back(xn);
      z.push_back(zn); 
      //aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(std::isnan(checknan)==false && zn<=0 && Flip==true){
      x.push_back(x1-xn);
      z.push_back(zn);
      //aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
      
    zn=zn+h;
    if(zn>0){
      i=dmax+2;      
    }
  }
  
  /* Map out the 2nd part of the reflected ray */
  zn=0.0;
  for(int i=0;i<dmax;i++){  
    params6a = {A_ice_Cnz, 0, 0, lvalueR};
    params6b = {A_ice_Cnz, 0, 0, lvalueR};
    xn=IceRayTracing::fDnfR_Cnz(zn,&params6a)-IceRayTracing::fDnfR_Cnz(z0,&params6b);
    checknan=IceRayTracing::fDnfR_Cnz(zn,&params6a);
    if(std::isnan(checknan)==false && Flip==false){
      x.push_back(xn);
      z.push_back(zn);
      //aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }
      
    if(std::isnan(checknan)==false && Flip==true){
      x.push_back(x1-xn);
      z.push_back(zn);
      //aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }

  params6a = {A_ice_Cnz, 0, 0, lvalueR};
  params6b = {A_ice_Cnz, 0, 0, lvalueR};
  xn=IceRayTracing::fDnfR_Cnz(zn,&params6a)-IceRayTracing::fDnfR_Cnz(z0,&params6b);
  if(Flip==true){
    x.push_back(x1-xn);
    z.push_back(zn);
    //aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
  }else{
    x.push_back(xn);
    z.push_back(zn);
    //aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
  }
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }

}

/* function for plotting and storing all the rays. This is for constant refractive index. */
void IceRayTracing::PlotAndStoreRays_Cnz(double x0,double z0, double z1, double x1, double lvalues[2], double A_ice_Cnz){

  double lvalueD=lvalues[0];
  double lvalueR=lvalues[1];

  vector <double> xD,zD;
  vector <double> xR,zR;
  GetFullDirectRayPath_Cnz(z0,x1,z1,lvalueD,A_ice_Cnz,xD,zD);
  GetFullReflectedRayPath_Cnz(z0,x1,z1,lvalueR,A_ice_Cnz,xR,zR);
  
}

/* This is the main raytracing function. x0 always has to be zero. z0 is the Tx depth in m and z1 is the depth of the Rx in m. Both depths are negative. x1 is the distance between them. This functions works for a constant refractive index */
double *IceRayTracing::IceRayTracing_Cnz(double x0, double z0, double x1, double z1, double A_ice_Cnz){

  /* define a pointer to give back the output of raytracing */ 
  double *output=new double[9];

  /* Store the ray paths in text files */
  bool PlotRayPaths=false;
  
  /*  ********This part of the code will try to get the Direct ray between Rx and Tx.********** */
  double* GetDirectRay=GetDirectRayPar_Cnz(z0,x1,z1,A_ice_Cnz);
  double RangD=GetDirectRay[0];
  double LangD=GetDirectRay[1];
  double timeD=GetDirectRay[2];
  double lvalueD=GetDirectRay[3];
  double checkzeroD=GetDirectRay[4];
  delete []GetDirectRay;
  
  /* ********This part of the code will try to get the Reflected ray between Rx and Tx.********** */
  double* GetReflectedRay=GetReflectedRayPar_Cnz(z0,x1,z1,A_ice_Cnz);
  double RangR=GetReflectedRay[0];
  double LangR=GetReflectedRay[1];
  double timeR=GetReflectedRay[2];
  double lvalueR=GetReflectedRay[3];
  double checkzeroR=GetReflectedRay[4];
  double timeR1=GetReflectedRay[5];
  double timeR2=GetReflectedRay[6];
  double AngleOfIncidenceInIce=GetReflectedRay[7];
  delete []GetReflectedRay; 

  /* This part of the code can be used if the user wants to plot the individual ray paths. This part of the code prints out the individual ray paths in text files and also plots them on a canvas */
  if(PlotRayPaths==true){
    double lvalues[2];
    lvalues[0]=lvalueD;
    lvalues[1]=lvalueR;
    
    PlotAndStoreRays_Cnz(x0,z0,z1,x1,lvalues,A_ice_Cnz);
  }  

  /* Fill in the output pointer after calculating all the results */
  output[0]=LangD;
  output[1]=LangR;
  output[2]=timeD;
  output[3]=timeR;
  output[4]=RangD;
  output[5]=RangR;
  
  /* fill in the output array part where you fill in the times for the two parts of the reflected or refracted rays */  
  output[6]=timeR1;
  output[7]=timeR2;
  output[8]=AngleOfIncidenceInIce;

  
  return output;
}

/* The set of functions starting with the name "fDa" are used in the minimisation procedure to find the launch angle (or the L parameter) for the direct ray */
double IceRayTracing::fDa_Air(double x,void *params){
  struct IceRayTracing::fDanfRa_params *p= (struct IceRayTracing::fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  double nz_Air=1;
  double AngleInAir=asin(x/nz_Air);
  double x1_Air=z1*tan(AngleInAir);
  
  struct IceRayTracing::fDnfR_L_params params1a = {A, IceRayTracing::GetB(1e-7), IceRayTracing::GetC(1e-7), -1e-7};
  struct IceRayTracing::fDnfR_L_params params1b = {A, IceRayTracing::GetB(z0), IceRayTracing::GetC(z0), z0};
  struct IceRayTracing::fDnfR_L_params params1c = {A, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), IceRayTracing::GetC(IceRayTracing::TransitionBoundary), -IceRayTracing::TransitionBoundary};
  struct IceRayTracing::fDnfR_L_params params1d = {A, IceRayTracing::GetB(-(IceRayTracing::TransitionBoundary+0.000001)), IceRayTracing::GetC(-(IceRayTracing::TransitionBoundary+0.000001)), -(IceRayTracing::TransitionBoundary+0.000001)};

  double distancez0z1=0;
  
  if(IceRayTracing::TransitionBoundary!=0){
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(1e-7)>IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(1e-7)<IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1c) + IceRayTracing::fDnfR_L(x,&params1d) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(1e-7)<IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(1e-7)<IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(1e-7)==IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(1e-7)==IceRayTracing::TransitionBoundary){
      distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1c) + IceRayTracing::fDnfR_L(x,&params1d) - IceRayTracing::fDnfR_L(x,&params1b);
    }
  }else{
    distancez0z1=IceRayTracing::fDnfR_L(x,&params1a) - IceRayTracing::fDnfR_L(x,&params1b);
  }
  // if(std::isnan(distancez0z1)){
  //   distancez0z1=1e9;
  // }
  if(std::isnan(x1_Air)){
    x1_Air=1e9;
  }
 
  double output=distancez0z1+x1_Air-x1;

  return output;
}

/* This functions works for the Direct ray and gives you back the launch angle, receive angle and propagation time of the ray together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
double* IceRayTracing::GetDirectRayPar_Air(double z0, double x1, double z1){

  double *output=new double[5];  
  
  /* First we setup the fDa function that will be minimised to get the launch angle (or the L parameter) for the direct ray. */
  gsl_function F1;
  struct IceRayTracing::fDanfRa_params params1= {IceRayTracing::A_ice, z0, x1, z1};
  F1.function = &IceRayTracing::fDa_Air;
  F1.params = &params1;
  
  /* In my raytracing solution given in the function IceRayTracing::fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0 and 1e-7 and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */ 
  double UpLimnz[]={IceRayTracing::Getnz(1e-7),IceRayTracing::Getnz(z0)};
  double* UpperLimitL=min_element(UpLimnz,UpLimnz+2);

  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function. */
  double lvalueD=IceRayTracing::FindFunctionRoot(F1,1e-7,UpperLimitL[0]);
  double LangD=asin(lvalueD/IceRayTracing::Getnz(z0))*(180.0/IceRayTracing::pi);
  double checkzeroD=IceRayTracing::fDa_Air(lvalueD,&params1);

  /* Get the propagation time for the direct ray using the ftimeD function after we have gotten the value of the L parameter. */
  struct IceRayTracing::ftimeD_params params2a = {IceRayTracing::A_ice, IceRayTracing::GetB(z0), -IceRayTracing::GetC(z0), IceRayTracing::c_light_ms,lvalueD};
  struct IceRayTracing::ftimeD_params params2b = {IceRayTracing::A_ice, IceRayTracing::GetB(1e-7), -IceRayTracing::GetC(1e-7), IceRayTracing::c_light_ms,lvalueD};
  struct IceRayTracing::ftimeD_params params2c = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary), IceRayTracing::c_light_ms, lvalueD};
  struct IceRayTracing::ftimeD_params params2d = {IceRayTracing::A_ice, IceRayTracing::GetB(IceRayTracing::TransitionBoundary+1e-7), -IceRayTracing::GetC(IceRayTracing::TransitionBoundary+1e-7), IceRayTracing::c_light_ms, lvalueD};

  /* we do the subtraction because we are measuring the time taken between the Tx and Rx positions */
  double timeD=0;
  if(IceRayTracing::TransitionBoundary!=0){
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(1e-7)>IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(1e-7,&params2b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(1e-7)<IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(IceRayTracing::TransitionBoundary+1e-7,&params2d) + IceRayTracing::ftimeD(IceRayTracing::TransitionBoundary,&params2c) - IceRayTracing::ftimeD(1e-7,&params2b);
    }
    if (fabs(z0)<IceRayTracing::TransitionBoundary && fabs(1e-7)<IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(1e-7,&params2b);
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(1e-7)<IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(1e-7,&params2b);
    }
    if (fabs(z0)==IceRayTracing::TransitionBoundary && fabs(1e-7)==IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(1e-7,&params2b);
    }
    if (fabs(z0)>IceRayTracing::TransitionBoundary && fabs(1e-7)==IceRayTracing::TransitionBoundary){
      timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(IceRayTracing::TransitionBoundary+1e-7,&params2d) + IceRayTracing::ftimeD(IceRayTracing::TransitionBoundary,&params2c) - IceRayTracing::ftimeD(1e-7,&params2b);
    }
  }else{
    timeD=IceRayTracing::ftimeD(-z0,&params2a) - IceRayTracing::ftimeD(1e-7,&params2b);
  }

  /* Setup the function that will be used to calculate the angle of reception for all the rays */
  gsl_function F5;
  struct IceRayTracing::fDnfR_params params5a = {IceRayTracing::A_ice, IceRayTracing::GetB(1e-7), -IceRayTracing::GetC(1e-7), lvalueD};
  double result, abserr;
  F5.function = &IceRayTracing::fDnfR;

  /* Calculate the recieve angle for direc rays by calculating the derivative of the function at the Rx position */
  F5.params = &params5a;
  gsl_deriv_central (&F5, 1e-7, 1e-8, &result, &abserr);
  double RangD=atan(result);
  
  /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
  if(z1==z0 && std::isnan(RangD)==true){
    RangD=180-LangD;
  }
  
  /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
  if(z1!=z0 && std::isnan(RangD)==true){
    RangD=90;
  }

  double AirAngle=asin(IceRayTracing::Getnz(1e-7)*sin(RangD));
  double AirHorizontalDistance=tan(AirAngle)*z1;
  double AirTime=AirHorizontalDistance/IceRayTracing::c_light_ms;
 
  timeD=timeD+AirTime;
  RangD=AirAngle*(180.0/IceRayTracing::pi);
  
  output[0]=RangD;
  output[1]=LangD;
  output[2]=timeD;
  output[3]=lvalueD;
  output[4]=checkzeroD;

  if(fabs(checkzeroD)>0.5){
    output[0]=-1000;
  }
  
  return output;
}

double *IceRayTracing::GeantRayTracer(double xT, double yT, double zT, double xR, double yR, double zR){
  double TxCor[3]={xT,yT,zT};
  double RxCor[3]={xR,yR,zR};
  
  ////For recording how much time the process took
  //auto t1b = std::chrono::high_resolution_clock::now();  
  
  double x0=0;/////always has to be zero
  double z0=TxCor[2];
  double x1=sqrt(pow(TxCor[0]-RxCor[0],2)+pow(TxCor[1]-RxCor[1],2));
  double z1=RxCor[2];

  double * getresults=IceRayTracing::IceRayTracing(x0,z0,x1,z1);

  // cout<<"*******For the Direct Ray********"<<endl;
  // cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
  // cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
  // cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;
  // cout<<"*******For the Refracted[1] Ray********"<<endl;
  // cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
  // cout<<"Recieve Angle: "<<getresults[10]<<" deg"<<endl;
  // cout<<"Propogation Time: "<<getresults[6]*pow(10,9)<<" ns"<<endl;
  // cout<<"*******For the Refracted[2] Ray********"<<endl;
  // cout<<"Launch Angle: "<<getresults[3]<<" deg"<<endl;
  // cout<<"Recieve Angle: "<<getresults[11]<<" deg"<<endl;
  // cout<<"Propogation Time: "<<getresults[7]*pow(10,9)<<" ns"<<endl;
  // cout<<"*******For the Reflected Ray********"<<endl;
  // cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
  // cout<<"Recieve Angle: "<<getresults[9]<<" deg"<<endl;
  // cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;   
  // cout<<"Incident Angle in Ice on the Surface: "<<getresults[18]<<" deg"<<endl;
  
  //cout<<" "<<endl;
  //cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;

  vector <double> OutputValues[4];
  
  if(getresults[8]!=-1000){ 
    // cout<<"*******For the Direct Ray********"<<endl;
    // cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;
    OutputValues[0].push_back(getresults[0]);
    OutputValues[1].push_back(getresults[8]);
    OutputValues[2].push_back(getresults[25]);
    OutputValues[3].push_back(getresults[4]*IceRayTracing::c_light_ms);
  }

  if(getresults[10]!=-1000){ 
    // cout<<"*******For the Refracted Ray 1********"<<endl;
    // cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<getresults[10]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<getresults[6]*pow(10,9)<<" ns"<<endl;
    OutputValues[0].push_back(getresults[2]);
    OutputValues[1].push_back(getresults[10]);
    OutputValues[2].push_back(getresults[27]);
    OutputValues[3].push_back(getresults[6]*IceRayTracing::c_light_ms);
  }

  if(getresults[11]!=-1000){ 
    // cout<<"*******For the Refracted Ray 2********"<<endl;
    // cout<<"Launch Angle: "<<getresults[3]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<getresults[11]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<getresults[7]*pow(10,9)<<" ns"<<endl;
    OutputValues[0].push_back(getresults[3]);
    OutputValues[1].push_back(getresults[11]);
    OutputValues[2].push_back(getresults[28]);
    OutputValues[3].push_back(getresults[7]*IceRayTracing::c_light_ms);
  }

  double *output=new double[4];
  
  if(OutputValues[3].size()!=0){
    int MinValueBin=0;//=TMath::LocMin(OutputValues[3].size(),OutputValues[3].data());
    double min=1e9;
    for(int i=0;i<OutputValues[3].size();i++){
      if(OutputValues[3][i]<min){
	min=OutputValues[3][i];
	MinValueBin=i;
      }
    }
    
    output[0]=OutputValues[0][MinValueBin];
    output[1]=OutputValues[1][MinValueBin];
    output[2]=OutputValues[2][MinValueBin];
    output[3]=OutputValues[3][MinValueBin];
  }else{
    output[0]=-1000;
    output[1]=-1000;
    output[2]=-1000;
    output[3]=-1000;
  }
  
  delete []getresults;
  
  // auto t2b = std::chrono::high_resolution_clock::now();
  // double durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();

  // double Duration=durationb;
  //cout<<"total time taken by the script: "<<Duration<<" micro s"<<endl;

  //cout<<output[0]<<" "<<output[1]<<" "<<output[2]<<" "<<output[3]<<endl;

  return output;
  
}

void IceRayTracing::MakeTable(double ShowerHitDistance,double zT){ 
  
  IceRayTracing::TotalStepsX_O=(IceRayTracing::GridWidthX/IceRayTracing::GridStepSizeX_O)+1;
  IceRayTracing::TotalStepsZ_O=(IceRayTracing::GridWidthZ/IceRayTracing::GridStepSizeZ_O)+1;

  IceRayTracing::GridPoints=IceRayTracing::TotalStepsX_O*IceRayTracing::TotalStepsZ_O;
  
  IceRayTracing::GridStartX=ShowerHitDistance-(IceRayTracing::GridWidthX/2);
  IceRayTracing::GridStopX=ShowerHitDistance+(IceRayTracing::GridWidthX/2);

  IceRayTracing::GridStartZ=-IceRayTracing::GridWidthZ;
  IceRayTracing::GridStopZ=0;
  
  //////For recording how much time the process took
  auto t1b = std::chrono::high_resolution_clock::now();  
  
  int totalpoints=0;
  for(int ix=0;ix<IceRayTracing::TotalStepsX_O;ix++){
    for(int iz=0;iz<IceRayTracing::TotalStepsZ_O;iz++){

      double xR=IceRayTracing::GridStartX+IceRayTracing::GridStepSizeX_O*ix;
      double zR=IceRayTracing::GridStartZ+IceRayTracing::GridStepSizeZ_O*iz;

      double *RTresults=IceRayTracing::GeantRayTracer(0, 0, zT, xR,0, zR);
      if(RTresults[2]!=-1000){
	IceRayTracing::GridPositionX.push_back(xR);
	IceRayTracing::GridPositionZ.push_back(zR);

	IceRayTracing::GridZValue[0].push_back(RTresults[0]);
	IceRayTracing::GridZValue[1].push_back(RTresults[1]);
	IceRayTracing::GridZValue[2].push_back(RTresults[2]);
	IceRayTracing::GridZValue[3].push_back(RTresults[3]);
	
	totalpoints++;
      }
    }
  }

  auto t2b = std::chrono::high_resolution_clock::now();
  double durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();

  double Duration=durationb/1000;

  cout<<"The table took "<<Duration<<" ms to make"<<endl;
  
  cout<<"total points are "<<totalpoints<<endl;
}

double IceRayTracing::GetInterpolatedValue(double xR, double zR, int rtParameter){

  int MinDistBin[4];
  double MinDist[9];

  double sum1=0;
  double sum2=0;
  double NewZValue=0;
  
  double minXbin=round((xR-IceRayTracing::GridStartX)/IceRayTracing::GridStepSizeX_O);
  double minZbin=round(fabs(zR-IceRayTracing::GridStartZ)/IceRayTracing::GridStepSizeZ_O);
     
  int newXbin=(minXbin/(IceRayTracing::TotalStepsX_O))*IceRayTracing::GridPoints;
  int newZbin=newXbin+minZbin;
      
  int count=0;
  if(minXbin<1){
    minXbin=1;
  }
  if(minZbin<1){
    minZbin=1;
  }

  if(minXbin+2>IceRayTracing::GridPoints){
    minXbin=IceRayTracing::GridPoints-2;
  }
  if(minZbin+2>IceRayTracing::GridPoints){
    minZbin=IceRayTracing::GridPoints-2;
  }  

  int startbinX=minXbin-1;
  int endbinX=minXbin+1;
  int startbinZ=minZbin-1;
  int endbinZ=minZbin+1;
     
  newXbin=((minXbin-1)/IceRayTracing::TotalStepsX_O)*IceRayTracing::GridPoints;
  newZbin=newXbin+(minZbin-1);
  int newich=newZbin;
  double minDist1=fabs(((xR-IceRayTracing::GridPositionX[newich])*(xR-IceRayTracing::GridPositionX[newich])+(zR-IceRayTracing::GridPositionZ[newich])*(zR-IceRayTracing::GridPositionZ[newich])));

  newXbin=((minXbin+1)/IceRayTracing::TotalStepsX_O)*IceRayTracing::GridPoints;
  newZbin=newXbin+(minZbin+1);
  newich=newZbin;
  double minDist2=fabs(((xR-IceRayTracing::GridPositionX[newich])*(xR-IceRayTracing::GridPositionX[newich])+(zR-IceRayTracing::GridPositionZ[newich])*(zR-IceRayTracing::GridPositionZ[newich])));

  if(minDist1<minDist2){
    startbinX=minXbin-1;
    endbinX=minXbin+1;
    startbinZ=minZbin-1;
    endbinZ=minZbin+1;
  }

  if(minDist1>minDist2){
    startbinX=minXbin;
    endbinX=minXbin+2;
    startbinZ=minZbin;
    endbinZ=minZbin+2;
  }
   
  sum1=0;
  sum2=0;
  NewZValue=0;
    
  for(int ixn=startbinX;ixn<endbinX;ixn++){
    for(int izn=startbinZ;izn<endbinZ;izn++){
      newXbin=((double)ixn/IceRayTracing::TotalStepsX_O)*IceRayTracing::GridPoints;
      newZbin=newXbin+izn;
	  
      newich=newZbin;
      if(newich<IceRayTracing::GridPoints){
	MinDist[count]=fabs(((xR-IceRayTracing::GridPositionX[newich])*(xR-IceRayTracing::GridPositionX[newich])+(zR-IceRayTracing::GridPositionZ[newich])*(zR-IceRayTracing::GridPositionZ[newich])));
	MinDistBin[count]=newich;
	    
	sum1+=(1.0/MinDist[count])*IceRayTracing::GridZValue[rtParameter][MinDistBin[count]];
	sum2+=(1.0/MinDist[count]);
      
	NewZValue=sum1/sum2;
	    
	if(MinDist[count]==0){
	  NewZValue=IceRayTracing::GridZValue[rtParameter][MinDistBin[count]];
	  izn=minZbin+3;
	  ixn=minXbin+3;
	}
	count++;
      }
    }
  }

  return NewZValue;
  
}

void IceRayTracing::GetRayTracingSolutions(double RxDepth, double Distance, double TxDepth, double TimeRay[2], double PathRay[2], double LaunchAngle[2], double RecieveAngle[2], int IgnoreCh[2], double IncidenceAngleInIce[2],vector <double> xRay[2], vector <double> zRay[2]){
  
  int DHits=0,RHits=0;
  double timeD, timeR, timeRa[2];
  double pathD, pathR, pathRa[2];
  double RangD, RangR, RangRa[2];
  double LangD, LangR, LangRa[2];

  int RayType[2];
  
  // double Distance=sqrt(pow(RxCor[0]-TxCor[0],2) + pow(RxCor[1]-TxCor[1],2));
  // double RxDepth=AntennaCoordRx[2]+AvgAntennaCoordRx[2];
  // double TxDepth=AntennaCoordTx[2];    
  //cout<<"rt parameters are "<<0<<" "<<TxDepth<<" "<<Distance<<" "<<RxDepth<<endl;

  double *RTresults=IceRayTracing::IceRayTracing(0, TxDepth, Distance,RxDepth);
   
  // cout<<"*******For the Direct Ray********"<<endl;
  // cout<<"Launch Angle: "<<RTresults[0]<<" deg"<<endl;
  // cout<<"Recieve Angle: "<<RTresults[8]<<" deg"<<endl;
  // cout<<"Propogation Time: "<<RTresults[4]*pow(10,9)<<" ns"<<endl;
  // cout<<"*******For the Reflected Ray********"<<endl;
  // cout<<"Launch Angle: "<<RTresults[1]<<" deg"<<endl;
  // cout<<"Recieve Angle: "<<RTresults[9]<<" deg"<<endl;
  // cout<<"Propogation Time: "<<RTresults[5]*pow(10,9)<<" ns"<<endl;   
  // cout<<"Incident Angle in Ice on the Surface: "<<RTresults[18]<<" deg"<<endl;
  // cout<<"*******For the Refracted[1] Ray********"<<endl;
  // cout<<"Launch Angle: "<<RTresults[2]<<" deg"<<endl;
  // cout<<"Recieve Angle: "<<RTresults[10]<<" deg"<<endl;
  // cout<<"Propogation Time: "<<RTresults[6]*pow(10,9)<<" ns"<<endl;
  // cout<<"*******For the Refracted[2] Ray********"<<endl;
  // cout<<"Launch Angle: "<<RTresults[3]<<" deg"<<endl;
  // cout<<"Recieve Angle: "<<RTresults[11]<<" deg"<<endl;
  // cout<<"Propogation Time: "<<RTresults[7]*pow(10,9)<<" ns"<<endl;
  
  timeD=RTresults[4]*pow(10,9);
  timeR=RTresults[5]*pow(10,9);
  timeRa[0]=RTresults[6]*pow(10,9);
  timeRa[1]=RTresults[7]*pow(10,9);

  pathD=RTresults[25];
  pathR=RTresults[26];
  pathRa[0]=RTresults[27];
  pathRa[1]=RTresults[28];
  
  RangD=RTresults[8];
  RangR=RTresults[9];
  RangRa[0]=RTresults[10];
  RangRa[1]=RTresults[11];

  LangD=RTresults[0];
  LangR=RTresults[1];
  LangRa[0]=RTresults[2];
  LangRa[1]=RTresults[3];
  
  TimeRay[0]=timeD;
  TimeRay[1]=timeR;
  
  PathRay[0]=pathD;
  PathRay[1]=pathR;

  RecieveAngle[0]=RangD;
  RecieveAngle[1]=RangR;

  LaunchAngle[0]=LangD;
  LaunchAngle[1]=LangR;

  RayType[0]=1;
  RayType[1]=2;

  IncidenceAngleInIce[0]=100;
  IncidenceAngleInIce[1]=RTresults[18];
  if(RangR==-1000){
    IncidenceAngleInIce[0]=100;
    IncidenceAngleInIce[1]=100;
  }

  double lvalueD=RTresults[19];
  double lvalueR=RTresults[20];
  double lvalueRa[2]={RTresults[21],RTresults[22]};

  double zmax[2]={0,0};
  zmax[0]=RTresults[23];  
  zmax[1]=RTresults[24];
  
  if(RangD!=-1000){
    TimeRay[0]=timeD;
    PathRay[0]=pathD;
    RecieveAngle[0]=RangD;
    LaunchAngle[0]=LangD;
    RayType[0]=1;
  }

  if(RangR!=-1000){
    TimeRay[1]=timeR;
    PathRay[1]=pathR;
    RecieveAngle[1]=RangR;
    LaunchAngle[1]=LangR;
    RayType[1]=2;
  }
    
  if(RangRa[0]!=-1000 && RangD!=-1000){
    TimeRay[0]=timeD;
    PathRay[0]=pathD;
    RecieveAngle[0]=RangD;
    LaunchAngle[0]=LangD;
    TimeRay[1]=timeRa[0];
    PathRay[1]=pathRa[0];
    RecieveAngle[1]=RangRa[0];
    LaunchAngle[1]=LangRa[0];
    RayType[0]=1;
    RayType[1]=3;
  }
      
  if(RangRa[0]!=-1000 && RangR!=-1000){
    TimeRay[1]=timeR;
    PathRay[1]=pathR;
    RecieveAngle[1]=RangR;
    LaunchAngle[1]=LangR;
    TimeRay[0]=timeRa[0];
    PathRay[0]=pathRa[0];
    RecieveAngle[0]=RangRa[0];
    LaunchAngle[0]=LangRa[0];
    RayType[0]=3;
    RayType[1]=2;
  }

  if(RangRa[1]!=-1000 && RangD!=-1000){
    TimeRay[0]=timeD;
    PathRay[0]=pathD;
    RecieveAngle[0]=RangD;
    LaunchAngle[0]=LangD;
    TimeRay[1]=timeRa[1];
    PathRay[1]=pathRa[1];
    RecieveAngle[1]=RangRa[1];
    LaunchAngle[1]=LangRa[1];
    RayType[0]=1;
    RayType[1]=4;
  }
      
  if(RangRa[1]!=-1000 && RangR!=-1000){
    TimeRay[1]=timeR;
    PathRay[1]=pathR;
    RecieveAngle[1]=RangR;
    LaunchAngle[1]=LangR;
    TimeRay[0]=timeRa[1];
    PathRay[0]=pathRa[1];
    RecieveAngle[0]=RangRa[1];
    LaunchAngle[0]=LangRa[1];
    RayType[0]=4;
    RayType[1]=2;
  }

  if(RangRa[1]!=-1000 && RangRa[0]!=-1000){
    TimeRay[1]=timeRa[1];
    PathRay[1]=pathRa[1];
    RecieveAngle[1]=RangRa[1];
    LaunchAngle[1]=LangRa[1];
    TimeRay[0]=timeRa[0];
    PathRay[0]=pathRa[0];
    RecieveAngle[0]=RangRa[0];
    LaunchAngle[0]=LangRa[0];
    RayType[0]=3;
    RayType[1]=4;
  }

  if(RecieveAngle[1]==-1000 && RecieveAngle[0]==-1000 && RangRa[0]!=-1000){
    TimeRay[0]=timeRa[0];
    PathRay[0]=pathRa[0];
    RecieveAngle[0]=RangRa[0];
    LaunchAngle[0]=LangRa[0];
    RayType[0]=3;
  }

  if(RecieveAngle[1]==-1000 && RecieveAngle[0]==-1000 && RangRa[1]!=-1000){
    TimeRay[1]=timeRa[1];
    PathRay[1]=pathRa[1];
    RecieveAngle[1]=RangRa[1];
    LaunchAngle[1]=LangRa[1];
    RayType[1]=4;
  }

  IgnoreCh[0]=1;
  IgnoreCh[1]=1;
  
  if(RecieveAngle[0]==-1000){
    IgnoreCh[0]=0;
  }

  if(RecieveAngle[1]==-1000){
    IgnoreCh[1]=0;
  }

  if(TimeRay[0]>TimeRay[1] && RecieveAngle[0]!=-1000 && RecieveAngle[1]!=-1000){
    swap(LaunchAngle[0],LaunchAngle[1]);
    swap(RecieveAngle[0],RecieveAngle[1]);
    swap(TimeRay[0],TimeRay[1]);
    swap(RayType[0],RayType[1]);
    swap(PathRay[0],PathRay[1]);
  }
  
  for(int iray=0;iray<2;iray++){
    if(RayType[iray]==1 && IgnoreCh[iray]!=0){
      GetFullDirectRayPath(TxDepth,Distance,RxDepth,lvalueD,xRay[iray],zRay[iray]);
    }
    if(RayType[iray]==2 && IgnoreCh[iray]!=0){
      GetFullReflectedRayPath(TxDepth,Distance,RxDepth,lvalueR,xRay[iray],zRay[iray]);
    }
    if(RayType[iray]==3 && IgnoreCh[iray]!=0){
      GetFullRefractedRayPath(TxDepth,Distance,RxDepth,zmax[0],lvalueRa[0],xRay[iray],zRay[iray],1);
    }
    if(RayType[iray]==4 && IgnoreCh[iray]!=0){
      GetFullRefractedRayPath(TxDepth,Distance,RxDepth,zmax[1],lvalueRa[1],xRay[iray],zRay[iray],2);
    }
    
  }
    
  // cout<<"inside Ray 1 "<<TimeRay[0]<<" "<<LaunchAngle[0]<<" "<<RecieveAngle[0]<<" "<<IgnoreCh[0]<<endl;
  // cout<<"inside Ray 2 "<<TimeRay[1]<<" "<<LaunchAngle[1]<<" "<<RecieveAngle[1]<<" "<<IgnoreCh[1]<<endl;
  
  if(RecieveAngle[0]!=-1000){
    DHits++;
  }

  if(RecieveAngle[1]!=-1000){
    RHits++;
  }
  delete [] RTresults;
}
