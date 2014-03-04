// ------------  IMPORTANT!  ------------ //
//  Do not include this header directly!  //
// (It is for internal use by RayTrace.h) //
// -------------------------------------- //

//theta is the angle off of vertical (0 = straight up)
template<typename positionType
#ifdef DO_TRACE_WITH_CALLBACK
, typename callbackType
#endif
>
TraceRecord TraceFinder::doTrace(const double depth, const double theta, const rayTargetRecord& target, unsigned short allowedReflections, double frequency, double polarization
//Eugene added sol_error
, int &sol_error
#ifdef DO_TRACE_WITH_CALLBACK
								 , callbackType callback
#endif
								 ) const{

        //sol_error = 0; // begin with no solution error

	const unsigned long maxIter = 10000;
	//const unsigned long maxIter = 10;
	//const unsigned long maxIter = 100000;
	const double eps=1e-10;
	const double targetTol = 0.01; //try to get the final value of x this close to target.distance
	
	//std::cout << "  TraceFinder::doTrace with inital angle " << theta << std::endl;
	
	TraceRecord result;
	result.launchAngle=theta;
	
	positionType pos, scale, tiny;
	typename positionType::derivativeType derivatives;
	pos.z=depth;
	pos.theta=theta;
	tiny.makeTiny(1.0e-30);
	double length=0.0;
	
#ifdef DO_TRACE_WITH_CALLBACK
	callCallback(callback,pos,RK_FIRST_STEP);
#endif
	
	//intro for the case of downward going rays above the ice
	if(pos.z>0.0 && pos.theta>(pi/2.)){
		if(target.distance > pos.x-pos.z*tan(pos.theta)){
			//cos(theta)<0, so use subtraction here
			length-=pos.z/cos(pos.theta);
			pos.addTime(-pos.z/(speedOfLight*cos(pos.theta)));
			pos.x-=pos.z*tan(pos.theta);
			pos.z=0.0;
			result.reflectionAngle=-pos.theta; //store negative of value to distinguish from reflection
			correctAmplitudeTransmit(pos,polarization,1.0,rModel->indexOfRefraction(0.0));
			pos.theta=pi-asin(sin(pi-pos.theta)/rModel->indexOfRefraction(0.0));
#ifdef DO_TRACE_WITH_CALLBACK
			callCallback(callback,pos,RK_AIR_STEP);
#endif
		}
		else{
			pos.giveData(result);
			result.pathLen=(target.distance-pos.x)/sin(pos.theta);
			result.pathTime=result.pathLen/speedOfLight;
			result.miss = pos.z-target.depth+(target.distance-pos.x)/tan(pos.theta); //theta>pi/2 => cos(theta)<0
#ifdef DO_TRACE_WITH_CALLBACK
			callCallback(callback,pos,RK_AIR_STEP);
#endif
			return(result);
		}
	}
	
	//choose starting stepsize using a heuristic:
	//the greater the depth the smoother the ice, so the larger the steps we can use
	//however, make sure that the proposed step size is smaller than the total distance
	double h=std::max(0.1,std::min(0.2*-depth,0.1*target.distance));
	double hdid, hnext;
	unsigned long stepCount;
	for(stepCount=0; stepCount<maxIter; stepCount++){
		computeRayDerivatives(pos, derivatives, frequency);
		scale = abs(pos)+abs(h*derivatives)+tiny;
		
		rkStepControl(length, frequency, pos, derivatives, scale, h, eps, hdid, hnext);
#ifdef DO_TRACE_WITH_CALLBACK
		callCallback(callback,pos,RK_STEP);
#endif
		
		//std::cout << "Did step of size " << hdid << " to depth " << pos.z << ", angle now " << pos.theta << ", step " << stepCount << std::endl;
		//std::cout << "Elapsed time is now " << pos.time << std::endl;
		
		//check for completion
		if(pos.x >= (target.distance-targetTol))
                        {
                        //std::cout<<"doTrace break1"<<std::endl;
			break;
                        }
		//prevent overshooting
		if((pos.x+hnext*derivatives.x) > (target.distance+targetTol)){
			//std::cout << "Expect to go to x=" << (pos.x+hnext*derivatives.x) << std::endl;
			hnext = /*0.9**/(target.distance-pos.x)/derivatives.x;
			//std::cout << "Limiting step size, now expect to go to x=" << (pos.x+hnext*derivatives.x) << std::endl;
		}
		
		//the next step would cross the ice surface, going up
		if(pos.theta<(pi/2.) && (-pos.z<=targetTol || pos.z+hnext*derivatives.z>0.0)){
			if(-pos.z<=targetTol){ //the ray is already about touching the surface, send it on through
				if(allowedReflections & SurfaceReflection){ //just reflect the ray back down, recording the angle of incidence
					result.reflectionAngle=pos.theta;
					//do this _before_ changing theta:
					correctAmplitudeReflect(pos,polarization,rModel->indexOfRefraction(0.0),1.0);
					pos.theta=pi-pos.theta;
#ifdef DO_TRACE_WITH_CALLBACK
					callCallback(callback,pos,RK_REFLECT_STEP);
#endif
				}
				else //stop the ray here
                                        {
                                        //std::cout<<"doTrace break2"<<std::endl;
					break;
                                        }
			}
			else{ //the ray isn't close enough yet, so just stop it from overshooting
				//std::cout << "Expect to go to z=" << pos.z+hnext*derivatives.z << std::endl;
				hnext = /*0.9**/(-pos.z/derivatives.z); //dz_ds > 0, z<0, so htemp > 0
				//std::cout << "Limiting step size, now expect to go to z=" << pos.z+hnext*derivatives.z << std::endl;
			}
		}
		else if(pos.theta>(pi/2.) && (pos.z<=(maximum_ice_depth+targetTol) || pos.z+hnext*derivatives.z<maximum_ice_depth)){
			if(pos.z<=(maximum_ice_depth+targetTol)){ //the ray is suitably close to the rock
				if(allowedReflections & BedrockReflection){ //reflect the ray back up, recording the angle of incidence
					result.reflectionAngle=pos.theta;
					pos.theta=pi-pos.theta;
					//do this _after_ changing theta
					correctAmplitudeReflect(pos,polarization,rModel->indexOfRefraction(pos.z),/*TODO: correct IoR for bedrock?*/1.0);
#ifdef DO_TRACE_WITH_CALLBACK
					callCallback(callback,pos,RK_REFLECT_STEP);
#endif
				}
				else //stop the ray here
                                        {
                                        //std::cout<<"doTrace break3"<<std::endl;
					break;
                                        }
			}
			else{ //the ray isn't close enough yet, so just stop it from overshooting
				//std::cout << "Expect to go to z=" << pos.z+hnext*derivatives.z << std::endl;
				hnext = /*0.9**/((maximum_ice_depth-pos.z)/derivatives.z);
				//std::cout << "Limiting step size, now expect to go to z=" << pos.z+hnext*derivatives.z << std::endl;
			}
		}
		
		h = std::max(targetTol, hnext); //don't let step sizes become too small
                // test out log
                
                // test 
                //if(stepCount==maxIter-1){
                if(stepCount==maxIter-1 || (pos.z != pos.z || pos.theta != pos.theta) ){ // maybe make to do same thing with src, trg switched?
                        //std::cout << "Step arrived maxIter!" << std::endl;
                        //std::cout << "solution error!" << std::endl;
                        //sol_error = 1; // solution error
                        sol_error++; // solution error (error occurs -> sol_error > 0)
                        result.sol_error=sol_error;
                        break;
                        //std::cout << "Was attempting to trace from depth of " << depth << " at an angle of " << theta
                        //<< " to target " << target.distance << " meters away at a depth " << target.depth << std::endl;
                        //throw std::runtime_error("TraceFinder::doTrace: exceeded maximum allowed number of steps.");
                }


	}
        /*
	if(stepCount==maxIter){
		//std::cout << "Was attempting to trace from depth of " << depth << " at an angle of " << theta
		//<< " to target " << target.distance << " meters away at a depth " << target.depth << std::endl;
		throw std::runtime_error("TraceFinder::doTrace: exceeded maximum allowed number of steps.");
	}
        */
        result.sol_error=sol_error;


	//std::cout << "TraceFinder::doTrace: final position is (" << pos.x << ',' << pos.z << ')' << std::endl;
	result.pathLen=length;
	pos.giveData(result);
	//result.miss=(target.distance-pos.x)*tan(pi/2-result.receiptAngle) - target.depth + pos.z;
	//std::cout << result.receiptAngle << " \t";
	//std::cout << "  " << result.receiptAngle << ' ' << (target.distance-pos.x) << ' ' << target.depth << ' ' << pos.z << std::endl;
	result.miss=(target.distance-pos.x)/tan(result.receiptAngle) - target.depth + pos.z;
	
	return(result);
}

//theta is the angle off of vertical (0 = straight _down_, due to the direction of the z axis)
template<typename positionType
#ifdef DO_TRACE_WITH_CALLBACK
, typename callbackType
#endif
>
TraceRecord TraceFinder::doVerticalTrace(const double depth, const double theta, const rayTargetRecord& target, unsigned short allowedReflections, double frequency, double polarization
#ifdef DO_TRACE_WITH_CALLBACK
										 , callbackType callback
#endif
										 ) const{
	const unsigned long maxIter = 10000;
	const double eps=1e-8;
	const double targetTol = 0.01; //try to get the final value of x this close to target.distance
	
	//std::cout << "TraceFinder::doVerticalTrace with inital angle " << theta << std::endl;
	
	bool willReflect=false;
	if((allowedReflections&SurfaceReflection && theta<pi/2.) ||
	   (allowedReflections&BedrockReflection && theta>pi/2.))
		willReflect=true;
	
	TraceRecord result;
	result.launchAngle=theta;
	
	positionType pos,scale,tiny;
	pos.z=depth;
	pos.theta=theta;
	typename positionType::derivativeType derivatives;
	tiny.makeTiny(1.0e-30);
	double length=0.0;
	
#ifdef DO_TRACE_WITH_CALLBACK
	callCallback(callback,pos,RK_FIRST_STEP);
#endif
	
	//intro for the case of downward going rays above the ice
	if(pos.z>0.0 && pos.theta>(pi/2.)){
		//cos(theta)<0, so use subtraction here
		length-=pos.z/cos(pos.theta);
		pos.time-=pos.z/(speedOfLight*cos(pos.theta));
		pos.x-=pos.z*tan(pos.theta);
		pos.z=0.0;
		result.reflectionAngle=pos.theta;
		correctAmplitudeTransmit(pos,polarization,1.0,rModel->indexOfRefraction(0.0));
		pos.theta=pi-asin(sin(pi-pos.theta)/rModel->indexOfRefraction(0.0));
#ifdef DO_TRACE_WITH_CALLBACK
		callCallback(callback,pos,RK_FIRST_STEP);
#endif
	}
	
	//choose starting stepsize using a heuristic:
	//the greater the depth the smoother the ice, so the larger the steps we can use
	//however, make sure that the proposed step size is smaller than the total distance
	double h=std::max(0.1,std::min(0.2*-depth,0.1*std::abs(depth-target.depth)));
	double hdid, hnext;
	unsigned long stepCount;
	for(stepCount=0; stepCount<maxIter; stepCount++){
		computeRayDerivatives(pos, derivatives, frequency);
		scale = abs(pos)+abs(h*derivatives)+tiny;
		
		rkStepControl(length, frequency, pos, derivatives, scale, h, eps, hdid, hnext);
#ifdef DO_TRACE_WITH_CALLBACK
		callCallback(callback,pos,RK_STEP);
#endif
		//std::cout << "Did step of size " << hdid << " to depth " << pos.z << " at distance " << pos.x << ", angle now " << pos.theta << std::endl;
		//std::cout << "Elapsed time is now " << pos.time << std::endl;
		
		if(!willReflect){
			if(pos.theta>pi/2.){ //down-going
				//check for completion
				if(pos.z <= target.depth)
					break;
				//prevent overshooting
				if((pos.z+hnext*derivatives.z) < (target.depth-targetTol))
					hnext = /*0.9**/(target.depth-pos.z)/derivatives.z;
			}
			else{ //up-going
				//check for completion
				if(pos.z >= target.depth)
					break;
				//prevent overshooting
				if((pos.z+hnext*derivatives.z) > (target.depth+targetTol))
					hnext = /*0.9**/(pos.z-target.depth)/derivatives.z;
			}
		}
		
		//the next step would cross the ice surface, going up
		if(pos.theta<(pi/2.) && (-pos.z<=targetTol || pos.z+hnext*derivatives.z>0.0)){
			if(-pos.z<=targetTol){ //the ray is already about touching the surface, send it on through
				if(allowedReflections & SurfaceReflection){ //just reflect the ray back down, recording the angle of incidence
					result.reflectionAngle=pos.theta;
					//do this _before_ changing theta:
					correctAmplitudeReflect(pos,polarization,rModel->indexOfRefraction(0.0),1.0);
					pos.theta=pi-pos.theta;
					willReflect=false;
#ifdef DO_TRACE_WITH_CALLBACK
					callCallback(callback,pos,RK_REFLECT_STEP);
#endif
				}
				else //stop the ray here
					break;
			}
			else{ //the ray isn't close enough yet, so just stop it from overshooting
				//std::cout << "Expect to go to z=" << pos.z+hnext*derivatives.z << std::endl;
				hnext = /*0.9**/(-pos.z/derivatives.z); //dz_ds < 0, so htemp > 0
				//std::cout << "Limiting step size, now expect to go to z=" << pos.z+hnext*derivatives.z << std::endl;
			}
		}
		else if(pos.theta>(pi/2.) && (pos.z<=(maximum_ice_depth+targetTol) || pos.z+hnext*derivatives.z<maximum_ice_depth)){
			if(pos.z<=(maximum_ice_depth+targetTol)){ //the ray is suitably close to the rock
				if(allowedReflections & BedrockReflection){ //reflect the ray back up, recording the angle of incidence
					result.reflectionAngle=pos.theta;
					pos.theta=pi-pos.theta;
					//do this _after_ changing theta
					correctAmplitudeReflect(pos,polarization,rModel->indexOfRefraction(pos.z),/*TODO: correct IoR for bedrock?*/1.0);
					willReflect=false;
#ifdef DO_TRACE_WITH_CALLBACK
					callCallback(callback,pos,RK_REFLECT_STEP);
#endif
				}
				else //stop the ray here
					break;
			}
			else{ //the ray isn't close enough yet, so just stop it from overshooting
				//std::cout << "Expect to go to z=" << pos.z+hnext*derivatives.z << std::endl;
				hnext = /*0.9**/((maximum_ice_depth-pos.z)/derivatives.z);
				//std::cout << "Limiting step size, now expect to go to z=" << pos.z+hnext*derivatives.z << std::endl;
			}
		}
		
		h = std::max(targetTol, hnext); //don't let step sizes become too small
	}
	if(stepCount==maxIter){
		//std::cout << "Was attempting to trace from depth of " << depth << " at an angle of " << theta
		//<< " to target " << target.distance << " meters away at a depth " << target.depth << std::endl;
		throw std::runtime_error("TraceFinder::doTrace: exceeded maximum allowed number of steps.");
	}
	//std::cout << "TraceFinder::doTrace: final position is (" << pos.x << ',' << pos.z << ')' << std::endl;
	result.pathLen=length;
	pos.giveData(result);
	result.miss=sqrt((target.depth-pos.z)*(target.depth-pos.z) + (target.distance-pos.x)*(target.distance-pos.x));
	
	return(result);
}
