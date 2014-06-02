#ifndef RAYSOLVER_H
#define RAYSOLVER_H


#include <vector>
#include <iostream>
#include "RayTrace.h"
#include "RayTrace_IceModels.h"
#include "Vector.h"
#include <fstream>
#include <iomanip>
#include "Settings.h"
#ifndef __CINT__
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#endif

class Position;
class IceModel;

//--------------------------------------------------
// struct underline{
// 	static char esc;
// 	const std::string& str;
// 	underline(const std::string& s):str(s){}
// 	friend std::ostream& operator<<(std::ostream& os, const underline& u);
// };
// 
// char underline::esc=0x1B;
// 
// std::ostream& operator<<(std::ostream& os, const underline& u){
// 	return(os << underline::esc << "[4m" << u.str << underline::esc << "[0m");
// }
//-------------------------------------------------- 

template<typename positionType>

class pathPrinter{
public:
	pathPrinter(){}
	void operator()(const positionType& p, RayTrace::RKStepType stepType){
		std::cout << p.x << ' ' << p.z << '\n';
	}
};


template<typename positionType>

class pathStore_test{
public:

        ofstream pathfile;

        // constructor
	pathStore_test(std::string filename){
            pathfile.open( filename.c_str() );
        }

	void operator()(const positionType& p, RayTrace::RKStepType stepType){
            pathfile << p.x << ' ' << p.z << '\n';
	}
};




template<typename positionType>

class pathStore_vector{
public:

        std::vector < std::vector <double> > vectorout;

        // constructor
	pathStore_vector(){
            vectorout.resize(2); // x, z
        }


	void CopyVector( std::vector < std::vector < std::vector<double> > > &vectorname, int solnum ) {
            //cout<<"testCopyVector"<<endl;
            //cout<<"size : "<<(int)vectorname.size()<<endl;

            vectorname[solnum].resize(2); // x, z

            for (int step=0; step<(int)vectorout[0].size(); step++) {
            
                vectorname[solnum][0].push_back( vectorout[0][step] );
                vectorname[solnum][1].push_back( vectorout[1][step] );
                //cout<<"1 x : "<<vectorout[0][step]<<", y : "<<vectorout[1][step]<<endl;
            }
        }

        // delete
	void DelVector(){
        
            vectorout.clear();
        }

        /*
        // copy 
	void CopyVector()(std::vector < std::vector<double> > &vectorname){
        
            for (int step=0; step<(int)vectorout.size(); step++) {
                vectorname[0].push_back( vectorout[0][step] );
                vectorname[1].push_back( vectorout[1][step] );
            }
        }

        // delete
	void DelVector()(){
        
            vectorout.clear();
        }
        */

	void operator()(const positionType& p, RayTrace::RKStepType stepType){
            //pathfile << p.x << ' ' << p.z << '\n';
            vectorout[0].push_back( p.x );
            vectorout[1].push_back( p.z );
            //cout<<"0 x : "<<p.x<<", y : "<<p.z<<endl;
            //cout<<"0 x : "<<vectorout[0][(int)vectorout[0].size()-1]<<", y : "<<vectorout[1][(int)vectorout[1].size()-1]<<endl;
	}
};

/*
template<typename positionType>

class pathStore_vector_2{
public:

        std::vector < std::vector <double> > vectorout;

        // constructor
	pathStore_vector_2(std::vector < std::vector<double> > &vectorname){
            vectorout.resize(2); // x, z
            vectorout = &vectorname;
        }

	void operator()(const positionType& p, RayTrace::RKStepType stepType){
            //pathfile << p.x << ' ' << p.z << '\n';
            vectorout[0].push_back( p.x );
            vectorout[1].push_back( p.z );
	}
};
*/





class RaySolver {

    private:
        void Earth_to_Flat_same_depth(Position &source, Position &target, IceModel* antarctica);
        void Earth_to_Flat_same_angle(Position &source, Position &target, IceModel* antarctica);

    public:
        RaySolver();
//        RaySolver(int argc, char* argv[]);

        void test1();
        void Solve_Ray_org(Position &source, Position &target, std::vector < std::vector <double> > &outputs, Settings *settings1);

        //void Solve_Ray_org (double source_x, double source_y, double source_z, double target_x, double target_y, double target_z, double &travel_time, double &travel_dist, int &no_sol, Settings *settings1);
        void Solve_Ray_org (double source_x, double source_y, double source_z, double target_x, double target_y, double target_z, bool print_path=false );

        //void Solve_Ray(Position &source, Position &target, IceModel *antarctica, std::vector < std::vector <double> > &outputs, Settings *settings1);
        void Solve_Ray(Position &source, Position &target, IceModel *antarctica, std::vector < std::vector <double> > &outputs, Settings *settings1, std::vector < std::vector < std::vector <double> > > &RayStep );

        int source_over_surface;
        int solution_toggle;    // no solution : 0  solution exist : 1
};


#endif //RAYSOLVER_H
