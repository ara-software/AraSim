#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include "RayTrace.h"
#include "RayTrace_IceModels.h"

struct underline{
	static char esc;
	const std::string& str;
	underline(const std::string& s):str(s){}
	friend std::ostream& operator<<(std::ostream& os, const underline& u);
};

char underline::esc=0x1B;

std::ostream& operator<<(std::ostream& os, const underline& u){
	return(os << underline::esc << "[4m" << u.str << underline::esc << "[0m");
}

template<typename positionType>
class pathPrinter{
public:
	pathPrinter(){}
	void operator()(const positionType& p, RayTrace::RKStepType stepType){
		std::cout << p.x << ' ' << p.z << '\n';
	}
};

int main(int argc, char* argv[]){
	double ns,nd,nc;
	double src_x, src_y, src_z, trg_x, trg_y, trg_z;
	Vector src,trg;
	bool showLabels=true;
	bool dumpPaths=false;
	bool surface_reflect;
	bool bedrock_reflect;
	double requiredAccuracy;
	double frequency;
	double polarization;
	boost::shared_ptr<RayTrace::indexOfRefractionModel> refractionModel;
	std::string refractionName;
	boost::shared_ptr<RayTrace::attenuationModel> attenuationModel;
	std::string attenuationName;
	
	boost::program_options::options_description desc("Allowed options");
	desc.add_options()
		("help,h","print usage information")
		("src_x",boost::program_options::value<double>(&src_x),"source x coordinate")
		("src_y",boost::program_options::value<double>(&src_y),"source y coordinate")
		("src_z",boost::program_options::value<double>(&src_z),"source z coordinate")
		("trg_x",boost::program_options::value<double>(&trg_x),"target x coordinate")
		("trg_y",boost::program_options::value<double>(&trg_y),"target y coordinate")
		("trg_z",boost::program_options::value<double>(&trg_z),"target z coordinate")
		("quiet,q","quiet, hide column labels")
		("show_paths,p","show paths; print out x and z coordinates of points visited along each path. Suppresses ordinary output")
		("n_s",boost::program_options::value<double>(&ns)->default_value(1.35),"surface index of refraction")
		("n_d",boost::program_options::value<double>(&nd)->default_value(1.78),"deep index of refraction")
		("n_c",boost::program_options::value<double>(&nc)->default_value(.0132),"index of refraction transition coefficient")
		("reflect_surface",boost::program_options::value<bool>(&surface_reflect)->default_value(true),"whether to search for surface\nreflected solutions")
		("reflect_bedrock",boost::program_options::value<bool>(&bedrock_reflect)->default_value(false),"whether to search for bedrock\nreflected solutions")
		("accuracy",boost::program_options::value<double>(&requiredAccuracy)->default_value(0.1),"the maximum acceptable vertical miss distance in meters")
		("frequency",boost::program_options::value<double>(&frequency)->default_value(300),"the frequency of the signal, in MHz")
		("polarization",boost::program_options::value<double>(&polarization)->default_value(RayTrace::pi/2),"the angle of the signal polarization, relative to the plane of propagation, in radians")
		("index_of_refraction",boost::program_options::value<std::string>(&refractionName)->default_value("exponential"),
		 "the index of refraction function of the ice, recognized values are 'exponential', 'inverse_exponential', 'quadratic', 'todor_linear', 'todor_chi', and 'todor_LL'")
		("attenuation",boost::program_options::value<std::string>(&attenuationName)->default_value("besson"),"the attenuation function of the ice, recognized values are 'negligible' and 'besson'")
	;
	boost::program_options::positional_options_description p;
	p.add("src_x", 1);
	p.add("src_y", 1);
	p.add("src_z", 1);
	p.add("trg_x", 1);
	p.add("trg_y", 1);
	p.add("trg_z", 1);
	
	boost::program_options::variables_map vm;
	try{
		boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		boost::program_options::notify(vm);
	}catch(std::exception& except){
		std::cerr << "Caught an exception during argument parsing: " << except.what() << std::endl;
		return(1);
	}catch(...){
		std::cerr << "An unknown exception was caught during argument parsing" << std::endl;
		return(1);
	}
	
	if (vm.count ("help") || argc < 2) {
		std::cout << "Usage: ray_solver [OPTION]... " << underline("src.x") << ' ' << underline("src.y") << ' ' << underline("src.z") << ' '
		<< underline("trg.x") << ' ' << underline("trg.y") << ' ' << underline("trg.z") << std::endl;
		std::cout << desc << std::endl;
		std::cout << "Angles are measured in radians from the downward vertical" << std::endl;
		std::cout << "The reported attenuation includes the effects of the ice attenuation model and reflections." << std::endl;
		std::cout << "The (electric field) amplitude value reported combines the attenuation with the divergence of rays." << std::endl;
		return(0);
	}
	if(vm.count("quiet"))
	   showLabels=false;
	if(vm.count("show_paths")){
		dumpPaths=true;
		showLabels=false;
	}
	if(refractionName=="exponential")
		refractionModel=boost::shared_ptr<exponentialRefractiveIndex>(new exponentialRefractiveIndex(ns,nd,nc));
		//refractionModel=boost::make_shared<exponentialRefractiveIndex>(ns,nd,nc);
	else if(refractionName=="inverse_exponential")
		refractionModel=boost::shared_ptr<inverseExponentialRefractiveIndex>(new inverseExponentialRefractiveIndex(ns,nd,nc));
		//refractionModel=boost::make_shared<inverseExponentialRefractiveIndex>(ns,nd,nc);
	else if(refractionName=="quadratic")
		refractionModel=boost::shared_ptr<quadraticRefractiveIndex>(new quadraticRefractiveIndex(ns,nd,nc));
		//refractionModel=boost::make_shared<quadraticRefractiveIndex>(ns,nd,nc);
	else if(refractionName=="todor_linear")
		refractionModel=boost::shared_ptr<todorLinearRefractiveIndex>(new todorLinearRefractiveIndex(ns,nd));
		//refractionModel=boost::make_shared<todorLinearRefractiveIndex>(ns,nd);
	else if(refractionName=="todor_chi")
		refractionModel=boost::shared_ptr<todorChiRefractiveIndex>(new todorChiRefractiveIndex(ns,nd));
		//refractionModel=boost::make_shared<todorChiRefractiveIndex>(ns,nd);
	else if(refractionName=="todor_LL")
		refractionModel=boost::shared_ptr<todorLLRefractiveIndex>(new todorLLRefractiveIndex(ns,nd));
		//refractionModel=boost::make_shared<todorLLRefractiveIndex>(ns,nd);
	else{
		std::cerr << "Unrecognized refraction model name '" << refractionName << "'" << std::endl;
		return(1);
	}
	if(attenuationName=="besson")
		attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);
		//attenuationModel=boost::make_shared<basicAttenuationModel>();
	else if(attenuationName=="negligible")
		attenuationModel=boost::shared_ptr<negligibleAttenuationModel>(new negligibleAttenuationModel);
		//attenuationModel=boost::make_shared<negligibleAttenuationModel>();
	else{
		std::cerr << "Unrecognized attenuation model name '" << attenuationName << "'" << std::endl;
		return(1);
	}
	
	unsigned short refl = RayTrace::NoReflection;
	if(surface_reflect)
		refl|=RayTrace::SurfaceReflection;
	if(bedrock_reflect)
		refl|=RayTrace::BedrockReflection;
	
	src.SetX(src_x);
	src.SetY(src_y);
	src.SetZ(src_z);
	trg.SetX(trg_x);
	trg.SetY(trg_y);
	trg.SetZ(trg_z);

	std::vector<RayTrace::TraceRecord> paths;
	
	RayTrace::TraceFinder tf(refractionModel,attenuationModel);
	
	paths=tf.findPaths(src,trg,frequency/1.0e3,polarization,refl,requiredAccuracy);
	if(showLabels){
		if(!paths.empty()){
			std::cout << std::fixed 
			<< "path length (m) "
			<< "path time (ns) "
			<< "launch angle "
			<< "recipt angle "
			<< "reflect angle "
			<< "miss dist. "
			<< "attenuation "
			<< "amplitude"
			<< std::endl;
		}
		else
			std::cout << "No solutions" << std::endl;
	}
	if(!dumpPaths){
		for(std::vector<RayTrace::TraceRecord>::const_iterator it=paths.begin(); it!=paths.end(); ++it){
			double signal = tf.signalStrength(*it,src,trg,refl);
			std::cout << std::left << std::fixed 
			<< std::setprecision(2) << std::setw(15) << it->pathLen << ' '
			<< std::setprecision(2) << std::setw(14) << 1e9*it->pathTime << ' '
			<< std::setprecision(4) << std::setw(12) << it->launchAngle << ' '
			<< std::setprecision(4) << std::setw(12) << it->receiptAngle << ' '
			<< std::setprecision(3) << std::setw(13) << it->reflectionAngle << ' '
			<< std::setprecision(2) << std::setw(10) << it->miss << ' ' 
			<< std::scientific << std::setprecision(4) << std::setw(11) << it->attenuation << ' '
			//amplitude calculation, ignoring frequency response at both ends, angular response of receiver
			<< std::setw(10) << (it->attenuation*signal)
			<< std::endl;
		}
	}
	else{ //do write out path data
		pathPrinter<RayTrace::minimalRayPosition> print;
		for(std::vector<RayTrace::TraceRecord>::const_iterator it=paths.begin(); it!=paths.end(); ++it){
			tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, &print);
			std::cout << "\n\n";
		}
	}
	return(0);
}
