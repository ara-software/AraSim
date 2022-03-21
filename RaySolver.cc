#include "RaySolver.h"
#include "Position.h"
#include "IceModel.h"
#include <string>
#include <cstring>
#include <sstream>
#include <iostream>
#include "Settings.h"
#include "IceRayTracing.hh"

RaySolver::RaySolver() {
  //default constructor
  //
}


void RaySolver::Earth_to_Flat_same_depth (Position &source, Position &target, IceModel *antarctica) {   // set both target, source has same depth from the surface : angle should be somewhat difficult to change
  double D = target.Distance( source );
  double target_depth = antarctica->Surface(target.Lon(), target.Lat()) - target.R();
  double source_depth = antarctica->Surface(source.Lon(), source.Lat()) - source.R();

  //--------------------------------------------------
  //     std::cout<<"source R : "<<source.R()<<" source Ice surface : "<<antarctica->Surface(source.Lon(), source.Lat())<<"\n";
  // 
  //     std::cout<<"target_depth : "<<target_depth<<" and source_depth : "<<source_depth<<"\n";
  //     std::cout<<"Surface R at source : "<<antarctica->Surface(source.Lon(), source.Lat())<<"\n";
  //-------------------------------------------------- 
  target.SetX( 0. );
  target.SetY( 0. );
  target.SetZ( -target_depth );

  source.SetX( pow( D*D - pow(target_depth - source_depth, 2), 0.5) );
  source.SetY( 0. );
  source.SetZ( -source_depth );

}


void RaySolver::Earth_to_Flat_same_angle (Position &source, Position &target, IceModel *antarctica) {   // set target depth as a same (but source depth changed). at target, any angle will be not changed... But there are a chance for source be out of surface (which mostly don't have solution for ray_solver).
  double D = target.Distance( source );
  double target_depth = antarctica->Surface(target.Lon(), target.Lat()) - target.R();
  double ang_diff = target.Angle(source);
  double depth_diff = target.R() - source.R()*cos(ang_diff);

  //--------------------------------------------------
  //     std::cout<<"source R : "<<source.R()<<" source Ice surface : "<<antarctica->Surface(source.Lon(), source.Lat())<<"\n";
  // 
  //     std::cout<<"target_depth : "<<target_depth<<" and depth_diff : "<<depth_diff<<"\n";
  //     std::cout<<"Surface R at source : "<<antarctica->Surface(source.Lon(), source.Lat())<<"\n";
  //-------------------------------------------------- 
  target.SetX( 0. );
  target.SetY( 0. );
  target.SetZ( -target_depth );

  source.SetX( pow( D*D - depth_diff*depth_diff, 0.5) );
  source.SetY( 0. );
  source.SetZ( -target_depth - depth_diff );

}




// 
//
//  Solve_Ray_org : don't change any positions from earth to flat.
//
//  input positions are values for flat surface where z = 0 is at ice surface
//
//
void RaySolver::Solve_Ray_org (Position &source, Position &target, std::vector < std::vector <double> > &outputs, Settings *settings1) {
                        
  outputs.clear();

  //--------------------------------------------------
  //     int argc = 7;
  //     int arg_length = 40;
  //     char argv_tmp[argc][arg_length];
  //     char *argv[argc];
  //-------------------------------------------------- 

  int sol_no = 0;                 // solution number (for vector solutions)

  Position source_tmp = source;
  Position target_tmp = target;

  int test;

  //--------------------------------------------------
  //     Earth_to_Flat_same_depth (source_tmp, target_tmp, antarctica);
  //-------------------------------------------------- 
  //--------------------------------------------------
  //     Earth_to_Flat_same_angle (source_tmp, target_tmp, antarctica);
  //-------------------------------------------------- 
    
  // set error message in case posnu is above Surface
  if (source_tmp.GetZ() > 0.) {
    source_over_surface = 1;
  }
  else {
    source_over_surface = 0;
  }

  /*
    test = sprintf(argv_tmp[1], "--src_x=%f", source_tmp.GetX() );
    test = sprintf(argv_tmp[2], "--src_y=%f", source_tmp.GetY() );
    test = sprintf(argv_tmp[3], "--src_z=%f", source_tmp.GetZ() );
    test = sprintf(argv_tmp[4], "--trg_x=%f", target_tmp.GetX() );
    test = sprintf(argv_tmp[5], "--trg_y=%f", target_tmp.GetY() );
    test = sprintf(argv_tmp[6], "--trg_z=%f", target_tmp.GetZ() );


    argv[1] = &argv_tmp[1][0];   //source x
    argv[2] = &argv_tmp[2][0];   //source y
    argv[3] = &argv_tmp[3][0];   //source z
    argv[4] = &argv_tmp[4][0];   //target x
    argv[5] = &argv_tmp[5][0];   //target y
    argv[6] = &argv_tmp[6][0];   //target z


    std::cout<<"argc : "<<argc<<"\n";
    std::cout<<"argv[1] : "<<argv[1]<<"\n";
    std::cout<<"argv[2] : "<<argv[2]<<"\n";
    std::cout<<"argv[3] : "<<argv[3]<<"\n";
    std::cout<<"argv[4] : "<<argv[4]<<"\n";
    std::cout<<"argv[5] : "<<argv[5]<<"\n";
    std::cout<<"argv[6] : "<<argv[6]<<"\n";
  */

  //--------------------------------------------------
  //     std::cout<<"src_x : "<<source_tmp.GetX()<<"\n";
  //     std::cout<<"src_y : "<<source_tmp.GetY()<<"\n";
  //     std::cout<<"src_z : "<<source_tmp.GetZ()<<"\n";
  //     std::cout<<"trg_x : "<<target_tmp.GetX()<<"\n";
  //     std::cout<<"trg_y : "<<target_tmp.GetY()<<"\n";
  //     std::cout<<"trg_z : "<<target_tmp.GetZ()<<"\n";
  //-------------------------------------------------- 




  //--------------------------------------------------
  // int main(int argc, char* argv[]){
  //-------------------------------------------------- 
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
	
  /*
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
    //--------------------------------------------------
    // return(1);
    //-------------------------------------------------- 
    }catch(...){
    std::cerr << "An unknown exception was caught during argument parsing" << std::endl;
    //--------------------------------------------------
    // return(1);
    //-------------------------------------------------- 
    }
	
    if (vm.count ("help") || argc < 2) {
    //--------------------------------------------------
    // std::cout << "Usage: ray_solver [OPTION]... " << underline("src.x") << ' ' << underline("src.y") << ' ' << underline("src.z") << ' '
    // << underline("trg.x") << ' ' << underline("trg.y") << ' ' << underline("trg.z") << std::endl;
    // std::cout << desc << std::endl;
    // std::cout << "Angles are measured in radians from the downward vertical" << std::endl;
    // std::cout << "The reported attenuation includes the effects of the ice attenuation model and reflections." << std::endl;
    // std::cout << "The (electric field) amplitude value reported combines the attenuation with the divergence of rays." << std::endl;
    //-------------------------------------------------- 
    //--------------------------------------------------
    // return(0);
    //-------------------------------------------------- 
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
    //--------------------------------------------------
    // return(1);
    //-------------------------------------------------- 
    }
    if(attenuationName=="besson")
    attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);
    //attenuationModel=boost::make_shared<basicAttenuationModel>();
    else if(attenuationName=="negligible")
    attenuationModel=boost::shared_ptr<negligibleAttenuationModel>(new negligibleAttenuationModel);
    //attenuationModel=boost::make_shared<negligibleAttenuationModel>();
    else{
    std::cerr << "Unrecognized attenuation model name '" << attenuationName << "'" << std::endl;
    //--------------------------------------------------
    // return(1);
    //-------------------------------------------------- 
    }
  */

  // Defaults!!!!
  src = source_tmp;
  trg = target_tmp;

  if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 0){
    // South Pole Values (AraSim original default, based on RICE)
    ns = 1.35;
    nd = 1.78;
    nc = 0.0132;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 1){
    // South Pole Values (RICE (2004))
    ns = 1.35;
    nd = 1.78;
    nc = 0.014;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 2){
    // South Pole Values (Eisen (2003))
    ns = 1.30;
    nd = 1.78;
    nc = 0.02;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 3){
    // South Pole values (Gow)
    ns = 1.345;
    nd = 1.78;
    nc = 0.016;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 10){
    // Moore's Bay values (MB #1)
    ns = 1.32;
    nd = 1.78;
    nc = 0.029;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 11){
    // Moore's Bay values (MB #2)
    ns = 1.299;
    nd = 1.78;
    nc = 0.027;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 20){
    // Byrd values (Ebimuna (1983))
    ns = 1.316;
    nd = 1.78;
    nc = 0.0244;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 30){
    // Mizuho values (Ebimuna (1983))
    ns = 1.357;
    nd = 1.78;
    nc = 0.027;
  }
  else {
    // South Pole Values (AraSim original default, based on RICE)
    ns = 1.35;
    nd = 1.78;
    nc = 0.0132;
  } 




    

  //    cout << "ns, nd, nc: " << ns << ", " << nd << ", " << nc << endl;

  surface_reflect = true;
  bedrock_reflect = false;
  //requiredAccuracy = 0.1;
  //requiredAccuracy = 0.5;
  requiredAccuracy = 0.2;
  frequency = 300;
  polarization = RayTrace::pi/2;
    
    
  if (settings1->NOFZ == 1){
    
    refractionModel=boost::shared_ptr<exponentialRefractiveIndex>(new exponentialRefractiveIndex(ns,nd,nc));
    attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);

  } else if (settings1->NOFZ == 0){
    refractionModel=boost::shared_ptr<constantRefractiveIndex>(new constantRefractiveIndex(1.48));
    attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);
  }
  else if (settings1->NOFZ == 2){
    refractionModel=boost::shared_ptr<inverseExponentialRefractiveIndex>(new inverseExponentialRefractiveIndex(ns,nd,nc));
    attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);
  }

	
  unsigned short refl = RayTrace::NoReflection;
  if(surface_reflect)
    refl|=RayTrace::SurfaceReflection;
  if(bedrock_reflect)
    refl|=RayTrace::BedrockReflection;
  /*	
	src.SetX(src_x);
	src.SetY(src_y);
	src.SetZ(src_z);
	trg.SetX(trg_x);
	trg.SetY(trg_y);
	trg.SetZ(trg_z);
  */


  std::vector<RayTrace::TraceRecord> paths;
  std::vector<RayTrace::TraceRecord> paths_exc;	

  int sol_cnt;
  int sol_error;

  int sol_cnt_exc;
  int sol_error_exc;

	
  RayTrace::TraceFinder tf(refractionModel,attenuationModel);
                            
  std::vector < std::vector < std::vector <double> > > testvector;

  int SrcTrgExc = 0;  // 0 not exchanged, 1 exchanged

	
  paths=tf.findPaths(src,trg,frequency/1.0e3,polarization, sol_cnt, sol_error, refl,requiredAccuracy);

  // have to fix src trg exchanged case angle outputs
  //
  if ( sol_cnt > 0 && sol_error > 0 ) {   // this case do findPahts again with src, trg exchanged
    paths_exc=tf.findPaths(trg,src,frequency/1.0e3,polarization, sol_cnt_exc, sol_error_exc, refl,requiredAccuracy);

    if ( sol_error > sol_error_exc ) {  // exchanged one is better (less sol_error)
      paths = paths_exc;
      SrcTrgExc = 1;
    }
    else if (sol_error == sol_error_exc) {  // sol_error are same for both paths. Then use the better miss dist.
      std::vector<RayTrace::TraceRecord>::const_iterator it_tmp=paths.begin();
      std::vector<RayTrace::TraceRecord>::const_iterator it_tmp_exc=paths_exc.begin();
      if ( it_tmp->miss > it_tmp_exc->miss ) {    // if fist one miss dist is worse
	paths = paths_exc;
	SrcTrgExc = 1;
      }
    }
  }



  if(showLabels){
    if(!paths.empty()){
      //--------------------------------------------------
      // std::cout << std::fixed 
      // << "path length (m) "
      // << "path time (ns) "
      // << "launch angle "
      // << "recipt angle "
      // << "reflect angle "
      // << "miss dist. "
      // << "attenuation "
      // << "amplitude"
      // << std::endl;
      //-------------------------------------------------- 
      solution_toggle = 1;
    }
    else {
      //std::cout << "No solutions" << std::endl;
      solution_toggle = 0;
    }
  }
  if(!dumpPaths){
    for(std::vector<RayTrace::TraceRecord>::const_iterator it=paths.begin(); it!=paths.end(); ++it){
      if ( it->sol_error == 0 ) {
	//double signal = tf.signalStrength(*it,src,trg,refl, sol_error );
	//--------------------------------------------------
	// std::cout << std::left << std::fixed 
	// << std::setprecision(2) << std::setw(15) << it->pathLen << ' '
	// << std::setprecision(2) << std::setw(14) << 1e9*it->pathTime << ' '
	// << std::setprecision(4) << std::setw(12) << it->launchAngle << ' '
	// << std::setprecision(4) << std::setw(12) << it->receiptAngle << ' '
	// << std::setprecision(3) << std::setw(13) << it->reflectionAngle << ' '
	// << std::setprecision(2) << std::setw(10) << it->miss << ' ' 
	// << std::scientific << std::setprecision(4) << std::setw(11) << it->attenuation << ' '
	// //amplitude calculation, ignoring frequency response at both ends, angular response of receiver
	// << std::setw(10) << (it->attenuation*signal)
	// << std::endl;
	//-------------------------------------------------- 

	if ( SrcTrgExc == 0 ) { // src, trg as original
	  outputs.resize(3);

	  outputs[0].push_back(it->pathLen);
	  outputs[1].push_back(it->launchAngle);
	  outputs[2].push_back(it->receiptAngle);
	  //std::cout<<"outputs[0]["<<sol_no<<"] : "<<outputs[0][sol_no]<<"pathLen : "<<it->pathLen<<"\n";


	  std::string pathfilename;
	  if ( sol_no == 0 ) {
	    pathfilename = "./pathfile_0.txt";
	  }
	  else if ( sol_no == 1 ) {
	    pathfilename = "./pathfile_1.txt";
	  }
	  else {
	    pathfilename = "./pathfile.txt";
	  }

	  testvector.resize(sol_no+1); // x, z

	  // construct
	  pathStore_vector<RayTrace::minimalRayPosition> pathsave_test;
	  //pathStore_vector_2<RayTrace::minimalRayPosition> pathsave_test (testvector);
	  //pathStore_test<RayTrace::minimalRayPosition> pathsave_test ("./pathfile.txt");
	  //pathStore_test<RayTrace::minimalRayPosition> pathsave_test (pathfilename);
                
	  //tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, sol_error, &pathsave_test);
	  tf.doTrace<RayTrace::minimalRayPosition>(src.GetZ(), it->launchAngle, RayTrace::rayTargetRecord(trg.GetZ(),sqrt((trg.GetX()-src.GetX())*(trg.GetX()-src.GetX())+(trg.GetY()-src.GetY())*(trg.GetY()-src.GetY()))), refl, 0.0, 0.0, sol_error, &pathsave_test);


	  pathsave_test.CopyVector( testvector, sol_no );
	  pathsave_test.DelVector();

	  /*
	    double totalpath = 0.;
	    double dx, dz;

	    //cout<<"\nRayStep : "<<(int)testvector[0].size()<<endl;
	    for (int step=0; step<(int)testvector[sol_no][0].size(); step++ ) {

	    if ( step > 0 ) {
	    dx = fabs(testvector[sol_no][0][step-1] - testvector[sol_no][0][step]);
	    dz = fabs(testvector[sol_no][1][step-1] - testvector[sol_no][1][step]);
	    totalpath += sqrt( (dx*dx) + (dz*dz) );
	    }
	    }
	  */

	  //cout<<"pathLen : "<<it->pathLen<<", pathSum : "<<totalpath<<endl;


	  //testvector.clear();



	  sol_no++;
	}
	else if ( SrcTrgExc == 1 ) { // src, trg exchanged
	  outputs.resize(3);

	  outputs[0].push_back(it->pathLen);
	  outputs[1].push_back(PI - it->receiptAngle);
	  outputs[2].push_back(PI - it->launchAngle);
	  //std::cout<<"outputs[0]["<<sol_no<<"] : "<<outputs[0][sol_no]<<"pathLen : "<<it->pathLen<<"\n";


	  std::string pathfilename;
	  if ( sol_no == 0 ) {
	    pathfilename = "./pathfile_0.txt";
	  }
	  else if ( sol_no == 1 ) {
	    pathfilename = "./pathfile_1.txt";
	  }
	  else {
	    pathfilename = "./pathfile.txt";
	  }

	  testvector.resize(sol_no+1); // x, z

	  // construct
	  pathStore_vector<RayTrace::minimalRayPosition> pathsave_test;
	  //pathStore_vector_2<RayTrace::minimalRayPosition> pathsave_test (testvector);
	  //pathStore_test<RayTrace::minimalRayPosition> pathsave_test ("./pathfile.txt");
	  //pathStore_test<RayTrace::minimalRayPosition> pathsave_test (pathfilename);
                
	  //tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, sol_error, &pathsave_test);
	  tf.doTrace<RayTrace::minimalRayPosition>(src.GetZ(), PI - it->receiptAngle, RayTrace::rayTargetRecord(trg.GetZ(),sqrt((trg.GetX()-src.GetX())*(trg.GetX()-src.GetX())+(trg.GetY()-src.GetY())*(trg.GetY()-src.GetY()))), refl, 0.0, 0.0, sol_error, &pathsave_test);

	  pathsave_test.CopyVector( testvector, sol_no );
	  pathsave_test.DelVector();

	  /*
	    double totalpath = 0.;
	    double dx, dz;

	    //cout<<"\nRayStep : "<<(int)testvector[0].size()<<endl;
	    for (int step=0; step<(int)testvector[sol_no][0].size(); step++ ) {

	    if ( step > 0 ) {
	    dx = fabs(testvector[sol_no][0][step-1] - testvector[sol_no][0][step]);
	    dz = fabs(testvector[sol_no][1][step-1] - testvector[sol_no][1][step]);
	    totalpath += sqrt( (dx*dx) + (dz*dz) );
	    }
	    }
	  */

	  //cout<<"pathLen : "<<it->pathLen<<", pathSum : "<<totalpath<<endl;



	  //testvector.clear();



	  sol_no++;
	}

	


      }

    }
  }
  else{ //do write out path data
    // I didn't fix this part yet!
    int sol_error;

    pathPrinter<RayTrace::minimalRayPosition> print;
    for(std::vector<RayTrace::TraceRecord>::const_iterator it=paths.begin(); it!=paths.end(); ++it){
      //tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, &print);
      tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, sol_error, &print);
      //std::cout << "\n\n";
    }
  }
  //--------------------------------------------------
  // return(0);
  //-------------------------------------------------- 
}




// 
//
//  Solve_Ray_org : don't change any positions from earth to flat. (also this version read x, y, z conponents and output travel time and travel distance only)
//
//  input positions are values for flat surface where z = 0 is at ice surface
//
//
//void RaySolver::Solve_Ray_org (double source_x, double source_y, double source_z, double target_x, double target_y, double target_z ) {
void RaySolver::Solve_Ray_org (double source_x, double source_y, double source_z, double target_x, double target_y, double target_z, bool print_path ) {
                        
  //outputs.clear();

  int sol_no = 0;                 // solution number (for vector solutions)

  //Position source_tmp = source;
  //Position target_tmp = target;

  int test;

  // set error message in case posnu is above Surface
  if (source_z > 0.) {
    source_over_surface = 1;
  }
  else {
    source_over_surface = 0;
  }



  //--------------------------------------------------
  // int main(int argc, char* argv[]){
  //-------------------------------------------------- 
  double ns,nd,nc;
  double src_x, src_y, src_z, trg_x, trg_y, trg_z;
  Vector src,trg;
  bool showLabels=true;
  //bool dumpPaths=false;
  bool dumpPaths=print_path;
  bool surface_reflect;
  bool bedrock_reflect;
  double requiredAccuracy;
  double frequency;
  double polarization;
  boost::shared_ptr<RayTrace::indexOfRefractionModel> refractionModel;
  std::string refractionName;
  boost::shared_ptr<RayTrace::attenuationModel> attenuationModel;
  std::string attenuationName;
	

  // Defaults!!!!
  //src = source_tmp;
  //trg = target_tmp;
    
    
  //    if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 0){
  // South Pole Values    
  ns = 1.35;
  nd = 1.78;
  nc = 0.0132;

  //    } else {
  // Moore's Bay values (MB #1)
  //ns = 1.32;
  // nd = 1.78;
  //nc = 0.029;
  //    }

  //    cout << "ns, nd, nc: " << ns << ", " << nd << ", " << nc << endl;


  surface_reflect = true;
  bedrock_reflect = false;
  //requiredAccuracy = 0.1;
  //requiredAccuracy = 0.5;
  requiredAccuracy = 0.2;
  frequency = 300;
  polarization = RayTrace::pi/2;
    
  //if (settings1->NOFZ == 1){
    
  refractionModel=boost::shared_ptr<exponentialRefractiveIndex>(new exponentialRefractiveIndex(ns,nd,nc));
  attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);

  /*
    } else if (settings1->NOFZ == 0){
    refractionModel=boost::shared_ptr<constantRefractiveIndex>(new constantRefractiveIndex(1.48));
    attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);
    }
  */

	
  unsigned short refl = RayTrace::NoReflection;
  if(surface_reflect)
    refl|=RayTrace::SurfaceReflection;
  if(bedrock_reflect)
    refl|=RayTrace::BedrockReflection;

  src.SetX(source_x);
  src.SetY(source_y);
  src.SetZ(source_z);
  trg.SetX(target_x);
  trg.SetY(target_y);
  trg.SetZ(target_z);





  std::vector<RayTrace::TraceRecord> paths;
  std::vector<RayTrace::TraceRecord> paths_exc;


  RayTrace::TraceFinder tf(refractionModel,attenuationModel);

                        
  std::vector < std::vector < std::vector <double> > > testvector;

  int sol_cnt;
  int sol_error;

  int sol_cnt_exc;
  int sol_error_exc;

  int SrcTrgExc = 0;  // 0 not exchanged, 1 exchanged

	
  paths=tf.findPaths(src,trg,frequency/1.0e3,polarization, sol_cnt, sol_error, refl,requiredAccuracy);

  // have to fix src trg exchanged case angle outputs
  //
  if ( sol_cnt > 0 && sol_error > 0 ) {   // this case do findPahts again with src, trg exchanged
    paths_exc=tf.findPaths(trg,src,frequency/1.0e3,polarization, sol_cnt_exc, sol_error_exc, refl,requiredAccuracy);

    if ( sol_error > sol_error_exc ) {  // exchanged one is better (less sol_error)
      paths = paths_exc;
      SrcTrgExc = 1;
    }
    else if (sol_error == sol_error_exc) {  // sol_error are same for both paths. Then use the better miss dist.
      std::vector<RayTrace::TraceRecord>::const_iterator it_tmp=paths.begin();
      std::vector<RayTrace::TraceRecord>::const_iterator it_tmp_exc=paths_exc.begin();
      if ( it_tmp->miss > it_tmp_exc->miss ) {    // if fist one miss dist is worse
	paths = paths_exc;
	SrcTrgExc = 1;
      }
    }
  }






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
      solution_toggle = 1;
    }
    else {
      std::cout << "No solutions" << std::endl;
      solution_toggle = 0;
    }
  }
  if(!dumpPaths){
    for(std::vector<RayTrace::TraceRecord>::const_iterator it=paths.begin(); it!=paths.end(); ++it){
      if ( it->sol_error == 0 ) {
	//double signal = tf.signalStrength(*it,src,trg,refl);
	double signal = tf.signalStrength(*it,src,trg,refl, sol_error );
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


        

	std::string pathfilename;
	if ( sol_no == 0 ) {
	  pathfilename = "./pathfile_0.txt";
	}
	else if ( sol_no == 1 ) {
	  pathfilename = "./pathfile_1.txt";
	}
	else {
	  pathfilename = "./pathfile.txt";
	}

	testvector.resize(sol_no+1); // x, z

	// construct
	pathStore_vector<RayTrace::minimalRayPosition> pathsave_test;
	//pathStore_vector_2<RayTrace::minimalRayPosition> pathsave_test (testvector);
	//pathStore_test<RayTrace::minimalRayPosition> pathsave_test ("./pathfile.txt");
	//pathStore_test<RayTrace::minimalRayPosition> pathsave_test (pathfilename);


            
	//tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, sol_error, &pathsave_test);
	tf.doTrace<RayTrace::minimalRayPosition>(src.GetZ(), it->launchAngle, RayTrace::rayTargetRecord(trg.GetZ(),sqrt((trg.GetX()-src.GetX())*(trg.GetX()-src.GetX())+(trg.GetY()-src.GetY())*(trg.GetY()-src.GetY()))), refl, 0.0, 0.0, sol_error, &pathsave_test);

	pathsave_test.CopyVector( testvector, sol_no );
	pathsave_test.DelVector();

	/*
	  double totalpath = 0.;
	  double dx, dz;

	  //cout<<"\nRayStep : "<<(int)testvector[0].size()<<endl;
	  for (int step=0; step<(int)testvector[sol_no][0].size(); step++ ) {

	  if ( step > 0 ) {
	  dx = fabs(testvector[sol_no][0][step-1] - testvector[sol_no][0][step]);
	  dz = fabs(testvector[sol_no][1][step-1] - testvector[sol_no][1][step]);
	  totalpath += sqrt( (dx*dx) + (dz*dz) );
	  }
	  }
	*/

	//cout<<"pathLen : "<<it->pathLen<<", pathSum : "<<totalpath<<endl;


	//testvector.clear();


	sol_no++;
	if (sol_no == 1) { // save only first sol
	  //travel_dist = it->pathLen;
	  //travel_time = it->pathTime; // in s
	}
      }

    }
  }
  else{ //do write out path data
    // I didn't fix this part yet!
    int sol_error;

    pathPrinter<RayTrace::minimalRayPosition> print;
    for(std::vector<RayTrace::TraceRecord>::const_iterator it=paths.begin(); it!=paths.end(); ++it){
      //tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, &print);
      tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, sol_error, &print);
      std::cout << "\n\n";
    }
  }
  //--------------------------------------------------
  // return(0);
  //-------------------------------------------------- 

  //no_sol = sol_no;
}













//--------------------------------------------------
// void RaySolver::Solve_Ray (Position &source, Position &target, IceModel *antarctica, std::vector <double> &outputs) {
//-------------------------------------------------- 
//void RaySolver::Solve_Ray (Position &source, Position &target, IceModel *antarctica, std::vector < std::vector <double> > &outputs, Settings *settings1) {
void RaySolver::Solve_Ray (Position &source, Position &target, IceModel *antarctica, std::vector < std::vector <double> > &outputs, Settings *settings1, std::vector < std::vector < std::vector <double> > > &RayStep ) {
                        
  outputs.clear();

  //--------------------------------------------------
  //     int argc = 7;
  //     int arg_length = 40;
  //     char argv_tmp[argc][arg_length];
  //     char *argv[argc];
  //-------------------------------------------------- 

  int sol_no = 0;                 // solution number (for vector solutions)

  Position source_tmp = source;
  Position target_tmp = target;

  double distance_org = source.Distance( target );

  int test;

  //--------------------------------------------------
  //     Earth_to_Flat_same_depth (source_tmp, target_tmp, antarctica);
  //-------------------------------------------------- 
  Earth_to_Flat_same_angle (source_tmp, target_tmp, antarctica);
    
  double distance_flat = source_tmp.Distance( target_tmp );
    
  //cout<<"source X Y Z are "<<source_tmp.GetX()<<" "<<source_tmp.GetY()<<" "<<source_tmp.GetZ()<<" target X Y Z are "<<target_tmp.GetX()<<" "<<target_tmp.GetY()<<" "<<target_tmp.GetZ()<<" ANALYTIC_RAYTRACE_MODE is "<<settings1->ANALYTIC_RAYTRACE_MODE<<endl;

  if ( fabs(distance_org - distance_flat) > 1. ) // if more than 1m difference
    cout<<"source, target distance, org: "<<distance_org<<", flat: "<<distance_flat<<endl;


  // set error message in case posnu is above Surface
  if (source_tmp.GetZ() > 0.) {
    source_over_surface = 1;
  }
  else {
    source_over_surface = 0;
  }

  /*
    test = sprintf(argv_tmp[1], "--src_x=%f", source_tmp.GetX() );
    test = sprintf(argv_tmp[2], "--src_y=%f", source_tmp.GetY() );
    test = sprintf(argv_tmp[3], "--src_z=%f", source_tmp.GetZ() );
    test = sprintf(argv_tmp[4], "--trg_x=%f", target_tmp.GetX() );
    test = sprintf(argv_tmp[5], "--trg_y=%f", target_tmp.GetY() );
    test = sprintf(argv_tmp[6], "--trg_z=%f", target_tmp.GetZ() );


    argv[1] = &argv_tmp[1][0];   //source x
    argv[2] = &argv_tmp[2][0];   //source y
    argv[3] = &argv_tmp[3][0];   //source z
    argv[4] = &argv_tmp[4][0];   //target x
    argv[5] = &argv_tmp[5][0];   //target y
    argv[6] = &argv_tmp[6][0];   //target z


    std::cout<<"argc : "<<argc<<"\n";
    std::cout<<"argv[1] : "<<argv[1]<<"\n";
    std::cout<<"argv[2] : "<<argv[2]<<"\n";
    std::cout<<"argv[3] : "<<argv[3]<<"\n";
    std::cout<<"argv[4] : "<<argv[4]<<"\n";
    std::cout<<"argv[5] : "<<argv[5]<<"\n";
    std::cout<<"argv[6] : "<<argv[6]<<"\n";
  */


  //--------------------------------------------------
  //     std::cout<<"src_x : "<<source_tmp.GetX()<<"\n";
  //     std::cout<<"src_y : "<<source_tmp.GetY()<<"\n";
  //     std::cout<<"src_z : "<<source_tmp.GetZ()<<"\n";
  //     std::cout<<"trg_x : "<<target_tmp.GetX()<<"\n";
  //     std::cout<<"trg_y : "<<target_tmp.GetY()<<"\n";
  //     std::cout<<"trg_z : "<<target_tmp.GetZ()<<"\n";
  //-------------------------------------------------- 





  //--------------------------------------------------
  // int main(int argc, char* argv[]){
  //-------------------------------------------------- 
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

  /*
	
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
    //--------------------------------------------------
    // return(1);
    //-------------------------------------------------- 
    }catch(...){
    std::cerr << "An unknown exception was caught during argument parsing" << std::endl;
    //--------------------------------------------------
    // return(1);
    //-------------------------------------------------- 
    }
	
    if (vm.count ("help") || argc < 2) {
    //--------------------------------------------------
    // std::cout << "Usage: ray_solver [OPTION]... " << underline("src.x") << ' ' << underline("src.y") << ' ' << underline("src.z") << ' '
    // << underline("trg.x") << ' ' << underline("trg.y") << ' ' << underline("trg.z") << std::endl;
    // std::cout << desc << std::endl;
    // std::cout << "Angles are measured in radians from the downward vertical" << std::endl;
    // std::cout << "The reported attenuation includes the effects of the ice attenuation model and reflections." << std::endl;
    // std::cout << "The (electric field) amplitude value reported combines the attenuation with the divergence of rays." << std::endl;
    //-------------------------------------------------- 
    //--------------------------------------------------
    // return(0);
    //-------------------------------------------------- 
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
    //--------------------------------------------------
    // return(1);
    //-------------------------------------------------- 
    }
    if(attenuationName=="besson")
    attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);
    //attenuationModel=boost::make_shared<basicAttenuationModel>();
    else if(attenuationName=="negligible")
    attenuationModel=boost::shared_ptr<negligibleAttenuationModel>(new negligibleAttenuationModel);
    //attenuationModel=boost::make_shared<negligibleAttenuationModel>();
    else{
    std::cerr << "Unrecognized attenuation model name '" << attenuationName << "'" << std::endl;
    //--------------------------------------------------
    // return(1);
    //-------------------------------------------------- 
    }

  */


  // Defaults!!!!
  src = source_tmp;
  trg = target_tmp;
    
  if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 0){
    // South Pole Values (AraSim original default, based on RICE)
    ns = 1.35;
    nd = 1.78;
    nc = 0.0132;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 1){
    // South Pole Values (RICE (2004))
    ns = 1.35;
    nd = 1.78;
    nc = 0.014;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 2){
    // South Pole Values (Eisen (2003))
    ns = 1.30;
    nd = 1.78;
    nc = 0.02;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 3){
    // South Pole values (Gow)
    ns = 1.345;
    nd = 1.78;
    nc = 0.016;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 10){
    // Moore's Bay values (MB #1)
    ns = 1.32;
    nd = 1.78;
    nc = 0.029;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 11){
    // Moore's Bay values (MB #2)
    ns = 1.299;
    nd = 1.78;
    nc = 0.027;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 20){
    // Byrd values (Ebimuna (1983))
    ns = 1.316;
    nd = 1.78;
    nc = 0.0244;
  } else if (settings1->RAY_TRACE_ICE_MODEL_PARAMS == 30){
    // Mizuho values (Ebimuna (1983))
    ns = 1.357;
    nd = 1.78;
    nc = 0.027;
  }
  else {
    // South Pole Values (AraSim original default, based on RICE)
    ns = 1.35;
    nd = 1.78;
    nc = 0.0132;
  } 


  //    cout << "ns, nd, nc: " << ns << ", " << nd << ", " << nc << endl;


  surface_reflect = true;
  bedrock_reflect = false;
  //requiredAccuracy = 0.1;
  //requiredAccuracy = 0.5;
  //requiredAccuracy = 0.2;
  frequency = 300;
  polarization = RayTrace::pi/2;


  // change accuracy value as a functin of phys_dist
  //if ( distance_org / 1000. > requiredAccuracy ) 
  //    requiredAccuracy = distance_org / 1000.; // loosen accuracy cut value

  if (settings1->Z_THIS_TOLERANCE == 0){
    requiredAccuracy = 0.2;

    if ( distance_org / 1000. > requiredAccuracy ) // change accuracy value as a functin of phys_dist
      requiredAccuracy = distance_org / 1000.; // loosen accuracy cut value
  }
  if (settings1->Z_THIS_TOLERANCE == 1){

    requiredAccuracy = settings1->Z_TOLERANCE;
    //cout<<"Z tolerance is "<<requiredAccuracy<<"m"<<endl;
  }


  if (settings1->NOFZ == 1){
        
    refractionModel=boost::shared_ptr<exponentialRefractiveIndex>(new exponentialRefractiveIndex(ns,nd,nc));
    attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);
        
  } else if (settings1->NOFZ == 0){
    refractionModel=boost::shared_ptr<constantRefractiveIndex>(new constantRefractiveIndex(1.5));
    attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);
  }
  else if (settings1->NOFZ == 2){
    refractionModel=boost::shared_ptr<inverseExponentialRefractiveIndex>(new inverseExponentialRefractiveIndex(ns,nd,nc));
    attenuationModel=boost::shared_ptr<basicAttenuationModel>(new basicAttenuationModel);
  }
    
  unsigned short refl = RayTrace::NoReflection;
  if(surface_reflect)
    refl|=RayTrace::SurfaceReflection;
  if(bedrock_reflect)
    refl|=RayTrace::BedrockReflection;

  /*	
	src.SetX(src_x);
	src.SetY(src_y);
	src.SetZ(src_z);
	trg.SetX(trg_x);
	trg.SetY(trg_y);
	trg.SetZ(trg_z);
  */


  // test print out location
  //std::cout<<"--src_x="<<src.GetX()<<" --src_y="<<src.GetY()<<" --src_z="<<src.GetZ()<<std::endl;
  //std::cout<<"--trg_x="<<trg.GetX()<<" --trg_y="<<trg.GetY()<<" --trg_z="<<trg.GetZ()<<std::endl;

  int sol_cnt;
  int sol_error;

  int sol_cnt_exc;
  int sol_error_exc;
	
  if(settings1->ANALYTIC_RAYTRACE_MODE==1){

    double TimeRay[2];
    double PathRay[2];
    double LaunchAngle[2];
    double RecieveAngle[2];
    int IgnoreCh[2];
    double IncidenceAngleInIce[2];
    vector <double> xRay[2];
    vector <double> zRay[2];
      
    double Distance= sqrt( (source_tmp.GetX()-target_tmp.GetX())*(source_tmp.GetX()-target_tmp.GetX()) + (source_tmp.GetY()-target_tmp.GetY())*(source_tmp.GetY()-target_tmp.GetY()) ) ;
      
    IceRayTracing::GetRayTracingSolutions(source_tmp.GetZ(), Distance, target_tmp.GetZ(), TimeRay, PathRay, LaunchAngle, RecieveAngle, IgnoreCh, IncidenceAngleInIce,xRay,zRay);

    swap(LaunchAngle[0],RecieveAngle[0]);

    TimeRay[0]=TimeRay[0]*pow(10,-9);
    TimeRay[1]=TimeRay[1]*pow(10,-9);

    LaunchAngle[0]=(180-LaunchAngle[0])*(PI/180);
    RecieveAngle[0]=(180-RecieveAngle[0])*(PI/180);

    LaunchAngle[1]=(LaunchAngle[1])*(PI/180);
    RecieveAngle[1]=(RecieveAngle[1])*(PI/180);

    if(IncidenceAngleInIce[0]!=100){
      IncidenceAngleInIce[0]=IncidenceAngleInIce[0]*(PI/180);
    }
    if(IncidenceAngleInIce[1]!=100){
      IncidenceAngleInIce[1]=IncidenceAngleInIce[1]*(PI/180);
    }
    // PathRay[0]=TimeRay[0]*IceRayTracing::c_light_ms;
    // PathRay[1]=TimeRay[1]*IceRayTracing::c_light_ms;
      
    solution_toggle = 0;
    for(int iray=0;iray<2;iray++){
      if(IgnoreCh[iray]==1){
	outputs.resize(5);

	outputs[0].push_back(PathRay[iray]);
	outputs[1].push_back(LaunchAngle[iray]);
	outputs[2].push_back(RecieveAngle[iray]);
	outputs[3].push_back(IncidenceAngleInIce[iray]);
	outputs[4].push_back(TimeRay[iray]);   // time in s (not ns)
      
	std::string pathfilename;
	if ( sol_no == 0 ) {
	  pathfilename = "./pathfile_0.txt";
	}
	else if ( sol_no == 1 ) {
	  pathfilename = "./pathfile_1.txt";
	}

	//cout<<"Ray "<<iray<<" "<<PathRay[iray]<<" "<<TimeRay[iray]<<" "<<LaunchAngle[iray]<<" "<<RecieveAngle[iray]<<" "<<IncidenceAngleInIce[iray]<<" "<<IgnoreCh[iray]<<endl;

	solution_toggle = 1;
	RayStep.resize(sol_no+1);

	RayStep[sol_no].resize(2); // x and z for the direct
	for(int step=0; step<xRay[iray].size(); step++){
	  RayStep[sol_no][0].push_back(xRay[iray][step]); // x component
	  RayStep[sol_no][1].push_back(zRay[iray][step]); // z component
	}   
	      
	sol_no++;
      }
    }	  
      
  }
	
  if(settings1->ANALYTIC_RAYTRACE_MODE==0){////do AraSim Numerical Raytracing

	
    RayTrace::TraceFinder tf(refractionModel,attenuationModel);

    std::vector<RayTrace::TraceRecord> paths;
    std::vector<RayTrace::TraceRecord> paths_exc;
	
    int SrcTrgExc = 0;  // 0 not exchanged, 1 exchanged

	
    paths=tf.findPaths(src,trg,frequency/1.0e3,polarization, sol_cnt, sol_error, refl,requiredAccuracy);

    //std::vector < std::vector < std::vector <double> > > testvector;

    //std::cout<<"sol cnt : "<<sol_cnt<<" sol error : "<<sol_error<<std::endl;
    // have to fix src trg exchanged case angle outputs
    //
    if ( sol_cnt > 0 && sol_error > 0 ) {   // this case do findPahts again with src, trg exchanged
      paths_exc=tf.findPaths(trg,src,frequency/1.0e3,polarization, sol_cnt_exc, sol_error_exc, refl,requiredAccuracy);

      std::cout<<"sol cnt exc : "<<sol_cnt<<" sol error exc : "<<sol_error<<std::endl;

      if ( sol_error > sol_error_exc ) {  // exchanged one is better (less sol_error)
	paths = paths_exc;
	SrcTrgExc = 1;
      }
      else if (sol_error == sol_error_exc) {  // sol_error are same for both paths. Then use the better miss dist.


	std::vector<RayTrace::TraceRecord>::const_iterator it_tmp=paths.begin();
	std::vector<RayTrace::TraceRecord>::const_iterator it_tmp_exc=paths_exc.begin();
	if ( it_tmp->miss > it_tmp_exc->miss ) {    // if fist one miss dist is worse
	  paths = paths_exc;
	  SrcTrgExc = 1;
	}
      }
    }
	
	
    if(showLabels){
      if(!paths.empty()){
	//--------------------------------------------------
	// std::cout << std::fixed 
	// << "path length (m) "
	// << "path time (ns) "
	// << "launch angle "
	// << "recipt angle "
	// << "reflect angle "
	// << "miss dist. "
	// << "attenuation "
	// << "amplitude"
	// << std::endl;
	//--------------------------------------------------
	solution_toggle = 1;
      }
      else {
	std::cout << "No solutions" << std::endl;
	solution_toggle = 0;
      }
    }
	
    if(!dumpPaths){
      for(std::vector<RayTrace::TraceRecord>::const_iterator it=paths.begin(); it!=paths.end(); ++it){
	if ( it->sol_error == 0 ) {
	  //double signal = tf.signalStrength(*it,src,trg,refl);
	  //--------------------------------------------------
	  // std::cout << std::left << std::fixed 
	  // << std::setprecision(2) << std::setw(15) << it->pathLen << ' '
	  // << std::setprecision(2) << std::setw(14) << 1e9*it->pathTime << ' '
	  // << std::setprecision(4) << std::setw(12) << it->launchAngle << ' '
	  // << std::setprecision(4) << std::setw(12) << it->receiptAngle << ' '
	  // << std::setprecision(3) << std::setw(13) << it->reflectionAngle << ' '
	  // << std::setprecision(2) << std::setw(10) << it->miss << ' ' 
	  // << std::scientific << std::setprecision(4) << std::setw(11) << it->attenuation << ' '
	  // //amplitude calculation, ignoring frequency response at both ends, angular response of receiver
	  // << std::setw(10) << (it->attenuation*signal)
	  // << std::endl;
	  //-------------------------------------------------- 

	  if ( SrcTrgExc == 0 ) { // src, trg as original
	    outputs.resize(5);

	    outputs[0].push_back(it->pathLen);
	    outputs[1].push_back(it->launchAngle);
	    outputs[2].push_back(it->receiptAngle);
	    outputs[3].push_back(it->reflectionAngle);
	    outputs[4].push_back( it->pathTime );   // time in s (not ns)
	    //std::cout<<"outputs[0]["<<sol_no<<"] : "<<outputs[0][sol_no]<<"pathLen : "<<it->pathLen<<"\n";
	    //std::cout<<"here outputs[0]["<<sol_no<<"] : "<<outputs[0][sol_no]<<" pathLen : "<<it->pathLen<<" "<<it->launchAngle<<" "<<it->receiptAngle<<" "<<it->reflectionAngle<<" "<<it->pathTime<<endl;

	    // test if raysol dist is shorter than physical distance
	    if ( distance_flat - it->pathLen >= 10. ) // if more than 10m difference
	      {
		//cout<<"source, target physical distance: "<<distance_flat<<", RayTrace: "<<it->pathLen<<", orgsol"<<endl;
		//cout<<"source(z:"<<src.GetZ()<<"), target(z:"<<trg.GetZ()<<") physical distance: "<<distance_flat<<", RayTrace: "<<it->pathLen<<", orgsol"<<endl;
		cout<<"source(x:"<<src.GetX()<<" y:"<<src.GetY()<<" z:"<<src.GetZ()<<"), target(x:"<<trg.GetX()<<" y:"<<trg.GetY()<<" z:"<<trg.GetZ()<<") physical distance: "<<distance_flat<<", RayTrace: "<<it->pathLen<<", orgsol"<<endl;
	      }


	    std::string pathfilename;
	    if ( sol_no == 0 ) {
	      pathfilename = "./pathfile_0.txt";
	    }
	    else if ( sol_no == 1 ) {
	      pathfilename = "./pathfile_1.txt";
	    }
	    else {
	      pathfilename = "./pathfile.txt";
	    }

	    //testvector.resize(sol_no+1);
	    RayStep.resize(sol_no+1);

	    // construct
	    pathStore_vector<RayTrace::minimalRayPosition> pathsave_test;
	    //pathStore_vector_2<RayTrace::minimalRayPosition> pathsave_test (testvector);
	    //pathStore_test<RayTrace::minimalRayPosition> pathsave_test ("./pathfile.txt");
	    //pathStore_test<RayTrace::minimalRayPosition> pathsave_test (pathfilename);
                
	    //tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, sol_error, &pathsave_test);
	    tf.doTrace<RayTrace::minimalRayPosition>(src.GetZ(), it->launchAngle, RayTrace::rayTargetRecord(trg.GetZ(),sqrt((trg.GetX()-src.GetX())*(trg.GetX()-src.GetX())+(trg.GetY()-src.GetY())*(trg.GetY()-src.GetY()))), refl, 0.0, 0.0, sol_error, &pathsave_test);

	    //pathsave_test.CopyVector( testvector, sol_no );
	    pathsave_test.CopyVector( RayStep, sol_no );
	    pathsave_test.DelVector();

	    /*
	      double totalpath = 0.;
	      double dx, dz;

	      //cout<<"\nRayStep : "<<(int)testvector[0].size()<<endl;
	      //for (int step=0; step<(int)testvector[sol_no][0].size(); step++ ) {
	      for (int step=0; step<(int)RayStep[sol_no][0].size(); step++ ) {

	      if ( step > 0 ) {
	      //dx = fabs(testvector[sol_no][0][step-1] - testvector[sol_no][0][step]);
	      //dz = fabs(testvector[sol_no][1][step-1] - testvector[sol_no][1][step]);
	      dx = fabs(RayStep[sol_no][0][step-1] - RayStep[sol_no][0][step]);
	      dz = fabs(RayStep[sol_no][1][step-1] - RayStep[sol_no][1][step]);
	      totalpath += sqrt( (dx*dx) + (dz*dz) );
	      }
	      }
	    */

	    //cout<<"pathLen : "<<it->pathLen<<", pathSum : "<<totalpath<<endl;



	    //testvector.clear();


	    sol_no++;
	  }
	  else if ( SrcTrgExc == 1 ) { // src, trg exchanged
	    outputs.resize(5);

	    outputs[0].push_back(it->pathLen);
	    outputs[1].push_back(PI - it->receiptAngle);
	    outputs[2].push_back(PI - it->launchAngle);
	    outputs[3].push_back(it->reflectionAngle);
	    outputs[4].push_back( it->pathTime );   // time in s (not ns)
	    //std::cout<<"outputs[0]["<<sol_no<<"] : "<<outputs[0][sol_no]<<"pathLen : "<<it->pathLen<<"\n";

	    // test if raysol dist is shorter than physical distance
	    if ( distance_flat - it->pathLen >= 10. ) // if more than 10m difference
	      {
		//cout<<"source, target physical distance: "<<distance_flat<<", RayTrace: "<<it->pathLen<<", excsol"<<endl;
		//cout<<"source(z:"<<src.GetZ()<<"), target(z:"<<trg.GetZ()<<") physical distance: "<<distance_flat<<", RayTrace: "<<it->pathLen<<", orgsol"<<endl;
		cout<<"source(x:"<<src.GetX()<<" y:"<<src.GetY()<<" z:"<<src.GetZ()<<"), target(x:"<<trg.GetX()<<" y:"<<trg.GetY()<<" z:"<<trg.GetZ()<<") physical distance: "<<distance_flat<<", RayTrace: "<<it->pathLen<<", orgsol"<<endl;
	      }



	    std::string pathfilename;
	    if ( sol_no == 0 ) {
	      pathfilename = "./pathfile_0.txt";
	    }
	    else if ( sol_no == 1 ) {
	      pathfilename = "./pathfile_1.txt";
	    }
	    else {
	      pathfilename = "./pathfile.txt";
	    }


	    //testvector.resize(sol_no+1);
	    RayStep.resize(sol_no+1);

	    // construct
	    pathStore_vector<RayTrace::minimalRayPosition> pathsave_test;
	    //pathStore_vector_2<RayTrace::minimalRayPosition> pathsave_test (testvector);
	    //pathStore_test<RayTrace::minimalRayPosition> pathsave_test ("./pathfile.txt");
	    //pathStore_test<RayTrace::minimalRayPosition> pathsave_test (pathfilename);
                
	    //tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, sol_error, &pathsave_test);
	    tf.doTrace<RayTrace::minimalRayPosition>(src.GetZ(), PI - it->receiptAngle, RayTrace::rayTargetRecord(trg.GetZ(),sqrt((trg.GetX()-src.GetX())*(trg.GetX()-src.GetX())+(trg.GetY()-src.GetY())*(trg.GetY()-src.GetY()))), refl, 0.0, 0.0, sol_error, &pathsave_test);

	    //pathsave_test.CopyVector( testvector, sol_no );
	    pathsave_test.CopyVector( RayStep, sol_no );
	    pathsave_test.DelVector();

	    /*
	      double totalpath = 0.;
	      double dx, dz;

	      //cout<<"\nRayStep : "<<(int)testvector[0].size()<<endl;
	      //for (int step=0; step<(int)testvector[sol_no][0].size(); step++ ) {
	      for (int step=0; step<(int)RayStep[sol_no][0].size(); step++ ) {

	      if ( step > 0 ) {
	      //dx = fabs(testvector[sol_no][0][step-1] - testvector[sol_no][0][step]);
	      //dz = fabs(testvector[sol_no][1][step-1] - testvector[sol_no][1][step]);
	      dx = fabs(RayStep[sol_no][0][step-1] - RayStep[sol_no][0][step]);
	      dz = fabs(RayStep[sol_no][1][step-1] - RayStep[sol_no][1][step]);
	      totalpath += sqrt( (dx*dx) + (dz*dz) );
	      }
	      }
	    */

	    //cout<<"pathLen : "<<it->pathLen<<", pathSum : "<<totalpath<<endl;


	    //testvector.clear();


	    sol_no++;
	  }





	}

      }


		
      // reorder the solutions (BAC Jan 2021)
      // it is useful at this point to re-order the solutions so that the direct solution 
      // (if it exists) is reported first
      // this copying method is rather inefficient in terms of memory allocation
      // but also, it's not obvious that it ever gets called 
      // (Brian couldn't find a case in testing on a few hundred events)
      // so probably it doesn't matter
      if (outputs.size()>0){
	int num_solutions = outputs[0].size();
	if(num_solutions>1){
	  // we only need to do something in the case of two solutions
	  double path_len_first_sol = outputs[0][0];
	  double path_len_second_sol = outputs[0][1];
	  if(path_len_second_sol < path_len_first_sol){
	    printf("The second path length (%.2f) is shorter than the first (%.2f)\n",path_len_second_sol, path_len_first_sol);
	    printf("Flipping the order in the output vectors...");

	    // if the shorter path length (which must be the direct solution)
	    // has ended up second, then we need to make a swap
	    // so that it's first in the output

	    // for the outputs vector, we flip the first and second columns
	    // that is, the direct and refracted/reflected solution
	    // first, make a copy of the original
	    vector< vector <double> > tmp_outputs(outputs);
	    // then make the flip
	    for(int item=0; item<outputs.size(); item++){
	      outputs[item][0] = tmp_outputs[item][1];
	      outputs[item][1] = tmp_outputs[item][0];
	    }

	    // for the RayStep, it's [solution][dimension (x or z)][step]
	    // so we need to be a little smarter
	    // first make a copy of the original
	    vector< vector< vector <double> > > tmp_RayStep(RayStep);

	    // now, totally empty the older vector (the one we need to return to the user at the end)
	    RayStep.clear();
	    RayStep.resize(2); // the RayStep needs a direct and refr/refl solution

	    // first, the direct solution (which is currently in tmp_RayStep[1])
	    RayStep[0].resize(2); // x and z for the direct
	    for(int step=0; step<tmp_RayStep[1][0].size(); step++){
	      RayStep[0][0].push_back(tmp_RayStep[1][0][step]); // x component
	      RayStep[0][1].push_back(tmp_RayStep[1][1][step]); // z component
	    }
	    // next, the reflected/refracted solution (which is currently in tmp_RayStep[0])
	    RayStep[1].resize(2); // x and z for the refracted/reflected
	    for(int step=0; step<tmp_RayStep[0][0].size(); step++){
	      RayStep[1][0].push_back(tmp_RayStep[0][0][step]); // x component
	      RayStep[1][1].push_back(tmp_RayStep[0][1][step]); // y component
	    }
	  }
	}		
      }
    }
    else{ //do write out path data
      // I didn't fix this part yet!
      int sol_error;

      pathPrinter<RayTrace::minimalRayPosition> print;
      for(std::vector<RayTrace::TraceRecord>::const_iterator it=paths.begin(); it!=paths.end(); ++it){
	tf.doTrace<RayTrace::minimalRayPosition>(src_z, it->launchAngle, RayTrace::rayTargetRecord(trg_z,sqrt((trg_x-src_x)*(trg_x-src_x)+(trg_y-src_y)*(trg_y-src_y))), refl, 0.0, 0.0, sol_error, &print);
	//std::cout << "\n\n";
      }
    }
    //--------------------------------------------------
    // return(0);
    //-------------------------------------------------- 


  }////AraSim numerical raytracing
	
}
