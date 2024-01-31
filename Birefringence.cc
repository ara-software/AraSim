#include "Settings.h"
#include "Detector.h"
#include "Position.h"
#include "Constants.h"
#include "TVector3.h"
#include "TGraph.h"
#include "Vector.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#include "Birefringence.h"

Birefringence::Birefringence(Settings *settings1) {

/*
 *	This birefringence model is based on AC paper: https://arxiv.org/pdf/2110.09015.pdf
 *	Parts of the code from cpol.cc: https://github.com/osu-particle-astrophysics/birefringence
 */

	//Reading values for principal axes and their asociated depths 
	string sn1file="./data/birefringence/n1.txt";
	string sn2file="./data/birefringence/n2.txt";
	string sn3file="./data/birefringence/n3.txt";

	Read_Indicatrix_Par(sn1file, sn2file, sn3file, settings1); //loading principal axes n1 (alpha in AC's paper), n2 (beta in AC's paper), n3 (gamma in AC's paper) and depths into nvec1, nvec2, nvec3, vdepths_n1, vdepths_n2, vdepths_n3
 	
	Smooth_Indicatrix_Par(); //smoothing values for principal axes nvec1, nvec2, nvec3
}


Birefringence::~Birefringence() {
	//default destructor
}


double Birefringence::getDeltaN(int BIAXIAL, vector<double> nvec, TVector3 rhat, double angle_iceflow, double n_e1, double n_e2, TVector3 &p_e1, TVector3 &p_e2) {

/*
 *
 *	- ASG 11/7/2023
 *
 *	This function finds the intersection elipse given by the incoming plane wave into the indicatrix. 
 *
 *	Here, nvec contains the principal axes for the indicatrix,
 *	rhat defines the direction of the incoming wave, 
 *	angle_iceflow is the angle of the iceflow from surveyor's coordinates, 
 *	n_e1 and n_e2 are the indices of refraction experienced by eigenstates p_e1 and p_e2 respectively. 
 *	p_e1 and p_e2 are found here and stored globally for further use.
 *
 *	Note, p_e1 and p_e2 are defined using surveyor's coordinates. The are rotated
 *	to local station-centric coordinates in Birefringence::Time_Diff_TwoRays() 
 *
 *	Something like n_e2 - n_e1 is the output of the function. 
 *
 *	AC would know further details.
 *
 */                   
                                                                                                                                                                   
	int FLIPPED=0;                                                                                                                                                   
  
 	if (rhat[2]<0.) {                                                                                                                                                
    		FLIPPED=1;                                                          
		rhat[2]=-1.*rhat[2];                                                                                                                                           	     }                 
  
  	TVector3 myy;
  	myy[0]=0.;
  	myy[1]=-1.;
  	myy[2]=0.;
  
  	TVector3 myz;
  	myz[0]=0.;
  	myz[1]=0.;
  	myz[2]=1.;

  	double rotate_toxalongiceflow[3][3]={{cos(angle_iceflow) , 1.*sin(angle_iceflow),0. },
                                       {-1.*sin(angle_iceflow), cos(angle_iceflow),0.},
                                       {0.,0.,1.}};

  	TVector3 rhat_iceflowalongx;
  	TVector3 myy_iceflowalongx;
  	TVector3 myz_iceflowalongx;
  	TVector3 x_iceflowalongx(1.,0.,0.);
  	TVector3 y_iceflowalongx(0.,1.,0.);
  	TVector3 z_iceflowalongx(0.,0.,1.);


  	for (int i=0;i<3;i++) {
    		double sum=0.;
    		for (int j=0;j<3;j++) {
      			sum+=rotate_toxalongiceflow[i][j]*rhat[j];
    		}
    		rhat_iceflowalongx[i]=sum;
  	}

  	TVector3 nominal_pe1=myz.Cross(rhat_iceflowalongx);

  	if (nominal_pe1.Mag()<HOWSMALLISTOOSMALL) {
    		cout << "myz is " << myz[0] << "\t" << myz[1] << "\t" << myz[2] << "\n";
    		cout << "rhat is " << rhat[0] << "\t" << rhat[1] << "\t" << rhat[2] << "\n";
    		cout << "rhat_iceflowalongx is " << rhat_iceflowalongx[0] << "\t" << rhat_iceflowalongx[1] << "\t" << rhat_iceflowalongx[2] << "\n";
    		cout << "cross of them is " << nominal_pe1[0] << "\t" << nominal_pe1[1] << "\t" << nominal_pe1[2] << "\n";
    		cout << "nominal_pe1 mag is " << nominal_pe1.Mag() << "\n";
  	}
  	nominal_pe1.SetMag(1.);
  	TVector3 nominal_pe2=rhat_iceflowalongx.Cross(nominal_pe1);

  	for (int i=0;i<3;i++) {
    		double sum=0.;
    		double sum2=0.;
    		for (int j=0;j<3;j++) {
      			sum+=rotate_toxalongiceflow[i][j]*myy[j];
      			sum2+=rotate_toxalongiceflow[i][j]*myz[j];
    		}
    		myy_iceflowalongx[i]=sum;
    		myz_iceflowalongx[i]=sum2;
  	}


  	double a=nvec[0];
  	double b=nvec[1];
  	double c=nvec[2];


  	double Ax=rhat_iceflowalongx[0];
  	double Ay=rhat_iceflowalongx[1];
  	double Az=rhat_iceflowalongx[2];


  	double A=1/(a*a)+(Ax*Ax)/(Az*Az*c*c);
  	double B=(2.*Ax*Ay)/(Az*Az*c*c);
  	double C=1/(b*b)+(Ay*Ay)/(Az*Az*c*c);

  	double theta_initial=atan2( B , A - C )/2.; // not sure this is rotated in the right direction - check this.                                                     

  	double Psi=PI/2.-atan2(abs(Az),sqrt(Ax*Ax+Ay*Ay));
  	double omega=-1.*(PI - atan2(Ay,Ax));

  	double rotate[3][3]={{cos(Psi)*cos(omega),cos(Psi)*sin(omega),sin(Psi)},
                       {-1.*sin(omega),cos(omega),0.},
                       {-1.*sin(Psi)*cos(omega),-1.*sin(Psi)*sin(omega),cos(Psi)}};


  	double myy_rotate[3];
  	double myz_rotate[3];

  	TVector3 nominal_pe1_rotate;
  	TVector3 nominal_pe2_rotate;

  	TVector3 x_iceflowalongx_rotate;
  	TVector3 y_iceflowalongx_rotate;

  	TVector3 rhat_rotate;
  	TVector3 tmpvec;
  	TVector3 tmpvec2;

  	for (int i=0;i<3;i++) {
    		double sum=0.;
    		double sum1=0.;
    		double sum2=0.;
    		double sum3=0.;
    		double sum4=0.;
   
		for (int j=0;j<3;j++) {
      			sum+=rotate[i][j]*rhat_iceflowalongx[j];
      			sum1+=rotate[i][j]*nominal_pe1[j];
      			sum2+=rotate[i][j]*nominal_pe2[j];
      			sum3+=rotate[i][j]*x_iceflowalongx[j];
      			sum4+=rotate[i][j]*y_iceflowalongx[j];

    		}
    		rhat_rotate[i]=sum;
    		nominal_pe1_rotate[i]=sum1;
    		nominal_pe2_rotate[i]=sum2;
    		x_iceflowalongx_rotate[i]=sum3;
    		y_iceflowalongx_rotate[i]=sum4;
  	}

  	for (int i=0;i<3;i++) {
    		double sum=0.;
    		double sum2=0.;
    		for (int j=0;j<3;j++) {
      			sum+=rotate[i][j]*myy_iceflowalongx[j];
      			sum2+=rotate[i][j]*myz_iceflowalongx[j];
    		}
    		myy_rotate[i]=sum;
    		myz_rotate[i]=sum2;
  	}


  	double rotate_T[3][3]={{0.}};
  	for (int i=0;i<3;i++) {
    		for (int j=0;j<3;j++) {
      			rotate_T[i][j]=rotate[j][i];
    		}
  	}

  	double epsilon_T=-1.*atan2(rotate_T[1][2],rotate_T[2][2]);
  	double Psi_T= asin(rotate_T[0][2]);
  	double omega_T=-1.*atan2(rotate_T[0][1],rotate_T[0][0]);

	TVector3 rhat_rotateback;
	TVector3 myy_rotateback;
  	TVector3 myz_rotateback;

	for (int i=0;i<3;i++) {
    		double sum=0.;
    		double sum2=0.;
    		double sum3=0.;
    		for (int j=0;j<3;j++) {
      			sum+=rotate_T[i][j]*rhat_rotate[j];
      			sum2+=rotate_T[i][j]*myy_rotate[j];
      			sum3+=rotate_T[i][j]*myz_rotate[j];
    		}
    		rhat_rotateback[i]=sum;
    		myy_rotateback[i]=sum2;
    		myz_rotateback[i]=sum3;
  	}

  	double Anew=1/(a*a)*(cos(Psi_T)*cos(Psi_T)*cos(omega_T)*cos(omega_T)) +
    		1/(b*b)*pow(cos(epsilon_T)*sin(omega_T)+sin(epsilon_T)*sin(Psi_T)*cos(omega_T),2) +
    		1/(c*c)*pow(sin(epsilon_T)*sin(omega_T)-cos(epsilon_T)*sin(Psi_T)*cos(omega_T),2);

  	double Bnew=1/(a*a)*(-2.*cos(Psi_T)*cos(Psi_T)*cos(omega_T)*sin(omega_T)) +
    		1/(b*b)*2.*(cos(epsilon_T)*sin(omega_T)+sin(epsilon_T)*sin(Psi_T)*cos(omega_T))*(cos(epsilon_T)*cos(omega_T)-sin(epsilon_T)*sin(Psi_T)*sin(omega_T)) +
    		1/(c*c)*2.*(sin(epsilon_T)*sin(omega_T)-cos(epsilon_T)*sin(Psi_T)*cos(omega_T))*(sin(epsilon_T)*cos(omega_T)+cos(epsilon_T)*sin(Psi_T)*sin(omega_T));

  	double Cnew=1/(a*a)*(cos(Psi_T)*cos(Psi_T)*sin(omega_T)*sin(omega_T)) +
    		1/(b*b)*pow(cos(epsilon_T)*cos(omega_T)-sin(epsilon_T)*sin(Psi_T)*sin(omega_T),2) +
    		1/(c*c)*pow(sin(epsilon_T)*cos(omega_T)+cos(epsilon_T)*sin(Psi_T)*sin(omega_T),2);

  	double Mnew=Anew;
  	double Pnew=Bnew/2.;
  	double Qnew=Bnew/2.;
  	double Rnew=Cnew;

	double lambda2_new=(1.*(Mnew+Rnew)+sqrt((Mnew-Rnew)*(Mnew-Rnew)+4*Pnew*Qnew))/2.;
  	double lambda1_new=(1.*(Mnew+Rnew)-sqrt((Mnew-Rnew)*(Mnew-Rnew)+4*Pnew*Qnew))/2.;

  	double ne1new=sqrt(1./lambda2_new);
  	double ne2new=sqrt(1./lambda1_new);

  	p_e1[0]=1.;
  	p_e1[1]=0.;
  	p_e1[2]=0.;

  	TVector3 rhat_rotatetheta=rhat_rotate;

  	double theta2=atan2( Bnew , Anew - Cnew )/2.; 


  	TVector3 findthefreakingaxis(1.,0.,0.);
  	TVector3 findthefreakingaxis_perp=findthefreakingaxis;

  	findthefreakingaxis.RotateZ(theta2);
  	findthefreakingaxis_perp.RotateZ(theta2+PI/2.);


  	TVector3 rhat_unrotate=rhat_rotate;

  	TVector3 tmpvec3;

  	for (int i=0;i<3;i++) {
    		double sum1=0.;
    		double sum2=0.;
    		double sum3=0.;
    		for (int j=0;j<3;j++) {
      			sum1+=rotate_T[i][j]*findthefreakingaxis[j];
      			sum2+=rotate_T[i][j]*findthefreakingaxis_perp[j];
      			sum3+=rotate_T[i][j]*rhat_rotate[j];
    		}
    		tmpvec[i]=sum1;
    		tmpvec2[i]=sum2;
    		tmpvec3[i]=sum3;
  	}

  	findthefreakingaxis=tmpvec;
  	findthefreakingaxis_perp=tmpvec2;
  	rhat_unrotate=tmpvec3;

  	TVector3 findthefreakingaxis_projecttoXY(findthefreakingaxis[0],findthefreakingaxis[1],0.);
  	TVector3 findthefreakingaxis_perp_projecttoXY(findthefreakingaxis_perp[0],findthefreakingaxis_perp[1],0.);
  	TVector3 yaxis(0.,1.,0.);

  	double anglebetweenthem=findthefreakingaxis_projecttoXY.Angle(yaxis);


  	double diffangle=theta_initial-anglebetweenthem;

  	diffangle=0.;
  	findthefreakingaxis_projecttoXY.RotateZ(diffangle);
  	findthefreakingaxis_perp_projecttoXY.RotateZ(diffangle);
  	rhat_unrotate.RotateZ(diffangle);

  	findthefreakingaxis[0]=findthefreakingaxis_projecttoXY[0];
  	findthefreakingaxis[1]=findthefreakingaxis_projecttoXY[1];

  	findthefreakingaxis_perp[0]=findthefreakingaxis_perp_projecttoXY[0];
  	findthefreakingaxis_perp[1]=findthefreakingaxis_perp_projecttoXY[1];

  	TVector3 rhat_rotateawayfromiceflow;


  	double rotate_backtonormal[3][3];

  	for (int i=0;i<3;i++) {
    		for (int j=0;j<3;j++) {
      			rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
    		}
  	}

  	p_e1=findthefreakingaxis;
  	p_e2=findthefreakingaxis_perp;

  	// p_e1 is in 1st or 4th quadrant in coordinate system where                                                                                                       
  	// ice is along the x axis                                                                                                                                         
  	if (!(p_e1.Phi()>-1*PI/2. && p_e1.Phi()<PI/2.))
    		p_e1=-1.*p_e1;
  	// p_e1 crossed with p_e2 should be in the vertical z direction                                                                                                  
  	if ((p_e1.Cross(p_e2)).Dot(myz)<0.)
    		p_e2=-1.*p_e2;



  	for (int i=0;i<3;i++) {
    		double sum3=0.;
    		double sum4=0.;
    		double sum5=0.;
    		for (int j=0;j<3;j++) {
      			sum3+=rotate_backtonormal[i][j]*rhat_unrotate[j];
      			sum4+=rotate_backtonormal[i][j]*findthefreakingaxis[j];
      			sum5+=rotate_backtonormal[i][j]*findthefreakingaxis_perp[j];
    		}
    		rhat_rotateawayfromiceflow[i]=sum3;
    		tmpvec[i]=sum4;
    		tmpvec2[i]=sum5;
  	}

  	findthefreakingaxis=tmpvec;
  	findthefreakingaxis_perp=tmpvec2;

  	p_e1=findthefreakingaxis;
  	p_e2=findthefreakingaxis_perp;


  	if (FLIPPED==1) {
    		p_e1[2]=-1.*p_e1[2];
    		p_e2[2]=-1.*p_e2[2];
    		rhat_rotateawayfromiceflow[2]=-1.*rhat_rotateawayfromiceflow[2];
  	}

  	n_e1=ne1new;
  	n_e2=ne2new;

  	// ray 1 is the one with the shortest index of refraction.                                                                                                       
  	if (n_e2<n_e1) {

    		double n_e1_temp=n_e1;
    		n_e1=n_e2;
    		n_e2=n_e1_temp;

    		TVector3 p_e1_temp=p_e1;
    		p_e1=p_e2;
    		p_e2=p_e1_temp;

  	}

  // for an isotropic medium it can tend to pick the orientation of the axes                                                                                       
  // somewhat randomly.                                                                                                                                            
  // this is why when an isotropic medium is chosen, I pick the n1                                                                                                 
  // principal axis to be very slightly longer than the other two.                                                                                                 
  // here is make sure that the 2st eigenvector is the one at 12 o'clock 

//POSSIBLE IMPROVEMENT: I don't think we want the isotropic case in AraSim. Just turn off birefringence for a better isotropic treatment. I will make BIAXIAL==-1 an inconsitency in the Settings.cc.
                                                                                          
  	if (BIAXIAL==-1) {
    		TVector3 temp1=myz.Cross(rhat);
    		TVector3 twelveoclock=rhat.Cross(temp1);
    		twelveoclock.SetMag(1.);
    		double mindotproduct=1000.;
    		TVector3 threeoclock=-1.*temp1;
    		threeoclock.SetMag(1.);
    		mindotproduct=fabs(p_e2.Dot(twelveoclock));
    		if (fabs(p_e1.Dot(twelveoclock)>mindotproduct)) {

      			double n_e1_temp=n_e1;
      			n_e1=n_e2;
      			n_e2=n_e1_temp;
      
      			TVector3 p_e1_temp=p_e1;
      			p_e1=p_e2;
      			p_e2=p_e1_temp;
      
    		}
    		p_e1=threeoclock;
    		p_e2=twelveoclock;
    
  	}
  
  	double deltan=0.;
  	if (BIAXIAL==-1)
    		deltan=0.;
  	else
    		deltan=ne2new-ne1new;

  	return deltan;
  
} // End of getDeltaN()


void Birefringence::Read_Indicatrix_Par(string sn1file, string sn2file, string sn3file, Settings *settings1){ //reads in data from n1file, n2file, n3file. These define the principal axes of the indicatrix


	int BIAXIAL = settings1->BIAXIAL; //BIAXIAL = 1 and = 0 are the only ones allowed now (by settings inconsistencies) for modeling biaxial (1) and uniaxial (0) birefringence, but I'll leave the code and capability for the isotropic case (BIAXIAL = -1). -ASG 11/7/2023 

	double thisn;
	double thisdepth;
	double firstn1;
	double firstn2;
	double firstn3;
	string stemp;

	ifstream n1file(sn1file.c_str());
	ifstream n2file(sn2file.c_str());
	ifstream n3file(sn3file.c_str());
	string line;
	char errorMessage[400];

	//counting lines
	int lineCount1 = 0;
        if(n1file.is_open()){
                while(n1file.peek()!=EOF){
                        getline(n1file, line);
                        lineCount1++;
                }
        }

	int lineCount2 = 0;
        if(n2file.is_open()){
                while(n2file.peek()!=EOF){
                        getline(n2file, line);
                        lineCount2++;
                }
        }

	int lineCount3 = 0;
        if(n3file.is_open()){
                while(n3file.peek()!=EOF){
                        getline(n3file, line);
                        lineCount3++;
                }
        }

	if ((lineCount1 != lineCount2) || (lineCount2 != lineCount3)){
		sprintf(errorMessage, "Files with principal axis for birefringence do not have the same number of lines");
                throw std::runtime_error(errorMessage);
	}
	
	int NDEPTHS_NS = lineCount1 - 1;

        n1file.clear(); // back to the beginning of the file again
        n1file.seekg(0, ios::beg);

        n2file.clear(); // back to the beginning of the file again
        n2file.seekg(0, ios::beg);
    	
        n3file.clear(); // back to the beginning of the file again
        n3file.seekg(0, ios::beg);

		
	//Reading in the principal axes
	n1file >> stemp;
    	for (int i=0;i<NDEPTHS_NS;i++) {
      		n1file >> thisdepth >> thisn;
      		vdepths_n1.push_back(-1.*thisdepth); 
      		n1vec.push_back(thisn);
     	}

    	n2file >> stemp;
    	for (int i=0;i<NDEPTHS_NS;i++) {//loops through this data
      		n2file >> thisdepth >> thisn;//piping into thisn
      		vdepths_n2.push_back(-1.*thisdepth);
      		if (BIAXIAL==1)//
			n2vec.push_back(thisn);//adds our data into thisn for certain properties
      		else if (BIAXIAL==0 || BIAXIAL==-1)
			n2vec.push_back(n1vec[i]);//adds data into n1vec? little confused on this part 

    	}
    
	n3file >> stemp; //n3file data into our stemp file!
    	for (int i=0;i<NDEPTHS_NS;i++) {//same loop as for n1 and n2
      		n3file >> thisdepth >> thisn;//from here below, same stuff as the last one for different biaxial values
      		vdepths_n3.push_back(-1.*thisdepth);
      		if (BIAXIAL==0 || BIAXIAL==1)
			n3vec.push_back(thisn);
      		else if (BIAXIAL==-1)
			n3vec.push_back(n1vec[i]+1.E-5); // the 1.E-5 is so the eigenvectors don't just go in completely random directions                    
    	}

    	if (CONSTANTINDICATRIX==1) {//Use first principal axes values always to define the indicatrix. I will set it to 0 in Birefringence.h but leave this code here. -ASG 11/7/2023
      		firstn1=n1vec[0];
      		firstn2=n2vec[0];
      		firstn3=n3vec[0];

      		int thissize=(int)n1vec.size();//length of vector

      		//empties the vectors
      		n1vec.clear();
      		n2vec.clear();
      		n3vec.clear();
      
      		//adds the first values to the vectors
      		for (int i=0;i<thissize;i++) {
			n1vec.push_back(firstn1);
			n2vec.push_back(firstn2);
			n3vec.push_back(firstn3);
      		}

    	}
} //End of Read_Indicatrix_Par 

void Birefringence::Smooth_Indicatrix_Par(){ //Smoothing the principal axes

    	vector<double> tmp;
    	tmp.clear();
    	tmp.resize(n1vec.size());

    	int NSMOOTH=5;
    	int min=(int)(((double)NSMOOTH)/2.);
    	for (int i=0;i<min;i++) {
      		tmp[i]=n1vec[i];
    	}
    	for (int i=n1vec.size()-(NSMOOTH-min);i<n1vec.size();i++) {
      		tmp[i]=n1vec[i];
    	}
    	for (int i=min;i<n1vec.size()-(NSMOOTH-min);i++) {
      		double tmpdouble=0.;
      		for (int j=i-min;j<i+(NSMOOTH-min);j++) {
			tmpdouble+=n1vec[j];
      		}
      		tmpdouble=tmpdouble/(double)NSMOOTH;
      		tmp[i]=tmpdouble;
    	}
    	n1vec=tmp;
    	tmp.clear();
    	tmp.resize(n2vec.size());
 
    	min=(int)(((double)NSMOOTH)/2.);
    	for (int i=0;i<min;i++) {
      		tmp[i]=n2vec[i];
    	}
    	for (int i=n2vec.size()-(NSMOOTH-min);i<n2vec.size();i++) {
     	 	tmp[i]=n2vec[i];
    	}
    	for (int i=min;i<n2vec.size()-(NSMOOTH-min);i++) {
      		double tmpdouble=0.;
      		for (int j=i-min;j<i+(NSMOOTH-min);j++) {
			tmpdouble+=n2vec[j];
      		}
      		tmpdouble=tmpdouble/(double)NSMOOTH;
      		tmp[i]=tmpdouble;
    	}
    	n2vec=tmp;

    	tmp.clear();
    	tmp.resize(n3vec.size());
 
    	min=(int)(((double)NSMOOTH)/2.);
    	for (int i=0;i<min;i++) {
      		tmp[i]=n3vec[i];
    	}
    	for (int i=n3vec.size()-(NSMOOTH-min);i<n3vec.size();i++) {
      		tmp[i]=n3vec[i];
    	}
    	for (int i=min;i<n3vec.size()-(NSMOOTH-min);i++) {
      		double tmpdouble=0.;	
      		for (int j=i-min;j<i+(NSMOOTH-min);j++) {
			tmpdouble+=n3vec[j];
      		}
      		tmpdouble=tmpdouble/(double)NSMOOTH;
      		tmp[i]=tmpdouble;
    	}
    	n3vec=tmp;

}// End of Smooth_Indicatrix_Par()


double Birefringence::Time_Diff_TwoRays(vector <double> res, vector <double> zs, double refl_angle, Position interaction_vertex, Settings *settings1){

        int stationID = settings1->DETECTOR_STATION;
        int BIREFRINGENCE = settings1->BIREFRINGENCE;

	if (BIREFRINGENCE==1 && refl_angle >= PI/2.0 ){ //Applying birefringence if BIREFRINGENCE is on and if we have no relected ray	

     		vector<double> nvec_thisstep;
        	nvec_thisstep.resize(3);

		TGraph *gn1=new TGraph(n1vec.size(),&vdepths_n1[0],&n1vec[0]);
 		TGraph *gn2=new TGraph(n2vec.size(),&vdepths_n2[0],&n2vec[0]);
 		TGraph *gn3=new TGraph(n3vec.size(),&vdepths_n3[0],&n3vec[0]);
	
		TVector3 rhat_thisstep;

		double n_e1;
		double n_e2;	
	
  		double rotate_toxalongiceflow[3][3]={{cos(angle_iceflow) , 1.*sin(angle_iceflow),0. },
                                       {-1.*sin(angle_iceflow), cos(angle_iceflow),0.},
                                       {0.,0.,1.}}; //For rotation from surveyor's to local station coordinates (same as with x along iceflow)
  	
		double rotate_toxalongglobal[3][3];

  		for (int i=0;i<3;i++) {
    			for (int j=0;j<3;j++) {
      				rotate_toxalongglobal[i][j]=rotate_toxalongiceflow[j][i];
    			}
  		} //For rotations from local station coordinates (same as with x along iceflow) to surveyor's 		

		TVector3 yhat(0.-interaction_vertex.GetX(),
                	      0.-interaction_vertex.GetY(),
                	      0.); // yhat points from interaction vertex to station in local coordinates	

		if (yhat.Mag()<HOWSMALLISTOOSMALL){
        		cout << "yhat mag is " << yhat.Mag() << "\n";
        	}

		yhat.SetMag(1.); //setting yhat to unity
	
		double deltantimeslength_alongpath=0.;

		for (int istep=0;istep<res.size();istep++) {
		
			nvec_thisstep.resize(3);

			// Getting values of principal axes at zs[istep] depth
         		nvec_thisstep[0]=gn1->Eval(zs[istep]);
        		nvec_thisstep[1]=gn2->Eval(zs[istep]);
        		nvec_thisstep[2]=gn3->Eval(zs[istep]);

			if (istep>0) {

				//The next three lines differ from AC's code by "-1.*" because she uses Uzair's ray tracer that returns positions on the ray solutions as moving down from target to source
       				rhat_thisstep[0]=(res[istep]-res[istep-1])*yhat[0];
           			rhat_thisstep[1]=(res[istep]-res[istep-1])*yhat[1];
           			rhat_thisstep[2]=(zs[istep]-zs[istep-1]);

				double length=rhat_thisstep.Mag();
           		
				if (rhat_thisstep.Mag()<HOWSMALLISTOOSMALL)
             				cout << "rhat_thisstep mag is " << rhat_thisstep.Mag() << "\n";

           			rhat_thisstep.SetMag(1.);

				if (rhat_thisstep.Mag()<1.E-8){
              				cout << "before calling getDeltaN at place 2, rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
				}
	

				TVector3 rhat_unrotate(rhat_thisstep[0], rhat_thisstep[1], rhat_thisstep[2]); //direction of propagation of the ray in local ccordinates

				for (int i=0;i<3;i++) {
    					double sum=0.;
    					for (int j=0;j<3;j++) {
      						sum+=rotate_toxalongglobal[i][j]*rhat_unrotate[j];
    					}
    					rhat_thisstep[i]=sum;
  				} //Overwrite rhat_thisstep so that it's given in surveyor's coordinates

	  
				double deltan_alongpath=getDeltaN(settings1->BIAXIAL,nvec_thisstep,rhat_thisstep,angle_iceflow,n_e1,n_e2,p_e1,p_e2); //In this function, we need to work with vectors in surveyor's coordinates; p_e1 and p_e2 are returned back in surveyor's coordinates as well
		
				TVector3 p_e1_toxalongiceflow(p_e1[0], p_e1[1], p_e1[2]);
				TVector3 p_e2_toxalongiceflow(p_e2[0], p_e2[1], p_e2[2]);

				for (int i=0;i<3;i++) {
    					double sum=0.;
    					double sum2=0.;
    					for (int j=0;j<3;j++) {
      						sum+=rotate_toxalongiceflow[i][j]*p_e1_toxalongiceflow[j];
      						sum2+=rotate_toxalongiceflow[i][j]*p_e2_toxalongiceflow[j];
    					}
    					p_e1[i]=sum;
    					p_e2[i]=sum2;
  				} //Overwrite p_e1 and p_e2 so that they are given in local coordinaes


				if (istep==1) {
					p_e1_src = p_e1.Unit();
					p_e2_src = p_e2.Unit();				
				} //store eigenstate polarizations at the location of the source
	
	
				if (p_e2.Mag()<HOWSMALLISTOOSMALL){
              				cout << "2, p_e2 is " << p_e2.Mag() << "\n";
				}

				deltantimeslength_alongpath+=deltan_alongpath*length;
	
			} //end if(istep>0)

		} //end for for(int istep...) loop

		double vtimediff = deltantimeslength_alongpath/CLIGHT*1.E9; //time difference in nanoseconds
	
		return vtimediff;
	}
     	else{
		return 0.;	
     	}

} // End of Time_Diff_TwoRays()


TVector3 Birefringence::Get_p_e1(){ 
	return p_e1; //Note this is p_e1 at the target, since we already iterated from source to target 
}


TVector3 Birefringence::Get_p_e2(){ 
	return p_e2; //Note this is p_e2 at the target, since we already iterated from source to target  
}


double Birefringence::Power_split_factor(Vector Pol_vector, int bire_ray_cnt, double refl_angle, Settings *settings1){ //split the power of ray between the two initial eigenstates (p_e1_src and p_e2_src) calculated in Time_Diff_TwoRays() 

/*
 *
 *	The initial ray (with polarization Pol_vector) splits into two (with polarizations p_e1_src and p_e2_src). The magnitude of their waveforms
 *	are scaled by Pol_Vec.Dot(p_e1_src) and Pol_Vec.Dot(p_e2_src). You can think of these factors like cos(theta) and sin(theta) where theta
 *	is the angle between Pol_vector and p_e1_src.
 *
 *	This way, (Vcos(theta))^2 + (Vsin(theta))^2 = V^2 ~ Power in the initial ray. So we have split the power appropiately.
 *
 *	Splitting the power only at the source (or vertex) ensures that the power deposited in eigenstate 1 and eigenstate 2
 *	stays withing the corresponding eigenstate during propagation.
 *
 *	-ASG 11/8/2023  
 *
 */

	int BIREFRINGENCE = settings1->BIREFRINGENCE; 

    	double split_factor = 0; 

    	if (BIREFRINGENCE==1){ // split power only if birefringence is on
 
		TVector3 Pol_Vec; 

		for(int i=0; i<3; i++){
			Pol_Vec[i] = Pol_vector[i];
		} // perhaps unnecesary but p_e1 and p_e2 are TVector3 so I converted Pol_vector to TVector3

		Pol_Vec = Pol_Vec.Unit(); 

		if (refl_angle < PI/2.){
			return split_factor = 1; // don't split if the ray is reflected
		}
		else if (bire_ray_cnt == 0 ){ // first ray 
			split_factor = Pol_Vec.Dot(p_e1_src); //p_e1 comes first 
		}
		else if (bire_ray_cnt == 1 ){ //second ray
			split_factor = Pol_Vec.Dot(p_e2_src); //p_e2 comes second
		}

     	}
     	else if (BIREFRINGENCE==0){
		split_factor = 1; // no birefringence, there's only one ray so no split is needed
     	}

	return split_factor;
} //End of Power_split_factor


void Birefringence::Principal_axes_polarization(Vector &Pol_vector, int bire_ray_cnt, int max_bire_ray_cnt, Settings *settings1){ //function to redefine the polarization of a ray

/*
 *	We'll need to make ray 1 have the polarization of eigenstate 1 at the target (p_e1). Similarly for ray 2.
 */

     	int BIREFRINGENCE = settings1->BIREFRINGENCE;

     	if(BIREFRINGENCE==1 && max_bire_ray_cnt == 2){ //only modify polarization if birefirngence is on and have two rays (we default to 1 ray if our ray solution reflects off the surface even if birefringence is on)

		if (bire_ray_cnt == 0){ //first ray
			Pol_vector = Vector(p_e1[0], p_e1[1], p_e1[2]); 
		}
		else if(bire_ray_cnt == 1){ //second ray
			Pol_vector = Vector(p_e2[0], p_e2[1], p_e2[2]);
		}
     	}
} //End of Principal_axes_polarization()


void Birefringence::Time_shift_and_power_split(double *V_forfft, int size, int T_shift, double split_factor, int bire_ray_cnt, int max_bire_ray_cnt, Settings *settings1){ //This function applies power split and time shift into the waveforms for eigenstates 1 and 2 


	int BIREFRINGENCE = settings1->BIREFRINGENCE;
     
     	if(BIREFRINGENCE==1 && max_bire_ray_cnt ==2){ //only apply changes if we have birefringence on and have two rays (we default to one ray if the ray solution is reflected off the surface)

		for (int n = 0; n < size; n++){
			if (bire_ray_cnt == 1){ //second ray
				V_forfft[n] *= split_factor; //only apply power split
                	}
                	else if (bire_ray_cnt == 0){ //first ray
				if ( n + T_shift < size){ //used defined bins
					V_forfft[n] = V_forfft[n + T_shift] * split_factor; //apply split and shift the waveform to the left (earlier in time)
                        	}     
                        	else{
                               		V_forfft[n] = 0.; //don't do anything with undefined bins      
                        	}
                	}
         	} //end of for-loop
  	}
} //End of Time_shift_and_power_split()


void Birefringence::Store_V_forfft_for_interference(double *V_forfft, double *V_forfft_bire, int size){ //store waveforms of each eigenstate after detector response to apply interference 
	
	for (int n = 0; n < size ;n++){
		V_forfft_bire[n] = V_forfft[n];
	}

}

void Birefringence::Two_rays_interference(double *V_forfft, double *V_forfft_bire_1, double *V_forfft_bire_2, int size, int max_bire_ray_cnt, Settings *settings1){ //apply interference of the waveforms from both eigenstates

     int BIREFRINGENCE = settings1->BIREFRINGENCE;
     
     if(BIREFRINGENCE==1 && max_bire_ray_cnt == 2){//overwrite V_forfft only if we had two rays in birefringence	
	for ( int n = 0; n < size; n++ ){
		V_forfft[n] = V_forfft_bire_1[n] + V_forfft_bire_2[n];
	}
     }
}

int Birefringence::Reflected_ray_remove_bire(double refl_angle, int max_bire_ray_cnt){ //default back to one ray if our ray solution is reflected

	if (refl_angle < PI/2.) {   // the ray is reflected at the surface
		return 1;
	}
	return max_bire_ray_cnt;
}
