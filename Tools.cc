#include <cmath>
#include <fstream> 
#include <iostream>
#include <vector>

#include "Tools.h"
#include "TSpline.h"
#include "TH2F.h"
#include "Constants.h"
#include <boost/math/interpolators/whittaker_shannon.hpp>

using std::cout;

void  Tools::MakeGraph(const int n,double *time,double *volts,TGraph *&mygraph,TH2F *&h2, double scalex,double scaley,string xaxistitle,string yaxistitle) {

    double maxtime=-1.E20;
    double maxv=-1.E20;
    double mintime=1E20;
    double minv=1.E20;

    double timecopy[n];
    double voltscopy[n];

    for (int i=0;i<n;i++) {
        timecopy[i]=time[i];
        voltscopy[i]=volts[i];
        timecopy[i]*=scalex;
        voltscopy[i]*=scaley;
    }
    

    for (int i=0;i<n;i++) {
        
        if (timecopy[i]>maxtime)
            maxtime=timecopy[i];
        if (timecopy[i]<mintime)
            mintime=timecopy[i];
        if (voltscopy[i]>maxv)
            maxv=voltscopy[i];
        if (voltscopy[i]<minv)
            minv=voltscopy[i];
    }

    mygraph=new TGraph(n,timecopy,voltscopy);

    h2=new TH2F("h2","",10*n,mintime*1.1,maxtime*1.1,100,minv*1.1,maxv*1.1);
    h2->SetLineWidth(3);
        h2->SetXTitle(xaxistitle.c_str());
            h2->SetYTitle(yaxistitle.c_str());
}

int Tools::iSum(int* thisarray,int n) {

    int sum=0;
    for (int i=0;i<n;i++) {
        sum+=thisarray[i];
    } //for
    return sum;
} //iSum

double Tools::getMaxMagnitude(vector<double> v) {
    double mag=0.;
    for (int i=0;i<(int)v.size();i++) {
        if (v[i]>mag)
            mag=v[i];

    }
    return mag;

}

void Tools::ShiftLeft(double *x,const int n,int ishift) {

    double x_temp[n];
    // shift the x array to the left by ishift bins and fill the gap with zeroes
    for (int i=0;i<n;i++) {
        x_temp[i]=x[i];
    }
    for (int i=0;i<n-ishift;i++) {
        x[i]=x_temp[i+ishift];
    }
    for (int i=n-ishift;i<n;i++) {
        x[i]=0.;
    }

}

void Tools::ShiftRight(double *x,const int n,int ishift) {

    double x_temp[n];
    // shift the x array to the right by ishift bins and fill the gap with zeroes
    for (int i=0;i<n;i++) {
        x_temp[i]=x[i];
    }
    for (int i=ishift;i<n;i++) {
        x[i]=x_temp[i-ishift];
    }
    for (int i=0;i<ishift;i++) {
        x[i]=0.;
    }

}

void Tools::realft(double *data, const int isign, int nsize){
    int i, i1, i2, i3, i4;
    double c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;
    //theta=3.141592653589793238/(nsize>>1);
    theta=3.141592653589793238/(double)(nsize>>1);
    if (isign == 1) {
        c2 = -0.5;
        four1(data,1,nsize);
    } else {
        c2=0.5;
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    for (i=1;i<(nsize>>2);i++) {
        i2=1+(i1=i+i);
        i4=1+(i3=nsize-i1);
        h1r=c1*(data[i1]+data[i3]);
        h1i=c1*(data[i2]-data[i4]);
        h2r= -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i;
        data[i2]=h1i+wr*h2i+wi*h2r;
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4]= -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
        data[0] = (h1r=data[0])+data[1];
        data[1] = h1r-data[1];
    } else {
        data[0]=c1*((h1r=data[0])+data[1]);
        data[1]=c1*(h1r-data[1]);
        four1(data,-1,nsize);
    }
}

void Tools::four1(double *data, const int isign,int nsize) {
    //int n,mmax,m,j,istep,i;
    int nn,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

    //int nn=nsize/2;
    int n=nsize/2;

    nn=n << 1;
    j=1;
    for (i=1;i<nn;i+=2) {
        if (j > i) {
            //SWAP(data[j-1],data[i-1]);
            //SWAP(data[j],data[i]);
            Exchange(data[j-1],data[i-1]);
            Exchange(data[j],data[i]);
        }
        m=n;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (nn > mmax) {
        istep=mmax << 1;
        theta=isign*(6.28318530717959/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=nn;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j-1]-wi*data[j];
                tempi=wr*data[j]+wi*data[j-1];
                data[j-1]=data[i-1]-tempr;
                data[j]=data[i]-tempi;
                data[i-1] += tempr;
                data[i] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

double Tools::dSquare(double *p) {
    return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
} //dSquare

int Tools::WhichIsMin(double *x,int n) {
    double min=1.E22;
    int imin=0;
    for (int i=0;i<n;i++) {
        if (x[i]<min) {
            min=x[i];
            imin=i;
        }
    } //for
    return imin;
} //WhichIsMin

void Tools::Print(double *p,int i) {
    for (int j=0;j<i;j++) {
        cout << p[j] << " ";
    } //for
    cout << "\n";
} //Print (double*,int)

void Tools::Print(int *p,int i) {
    for (int j=0;j<i;j++) {
        cout << p[j] << " ";
    } //for
    cout << "\n";
} //Print (double*,int)

void Tools::GetNextNumberAsString(ifstream& fin,ofstream& fout,string& number) {
    string temp;  
    getline(fin,temp); // get next line of the input file 
 
    fout << temp << "\n"; // output this line to the summary file

    int place=0; 
    place=temp.find_first_of(" \t"); // find where the first space or tab is

    number=temp.substr(0,place); // everything up until the first space is what we're looking for
} //GetNextNumberAsString

void Tools::GetNumbersAsStringArray(ifstream& fin, ofstream& fout,vector<string>& vnumbers, int nelements) {
    string temp;
    // getline(fin,temp);
    
    // fout << temp << "\n";

    // int place_previous=0;
    // int place_next;
    vnumbers.clear();
    string s;
    for (int n=0;n<nelements;n++) {
        fin >> s;
        fout << s << "\t"; // output this line to the summary file
        vnumbers.push_back(s);
        //    place_next=temp.find_first_of("\t",place_previous+1); // find where first tab is
        //vnumbers.push_back(temp.substr(place_previous,place_next-place_previous));
                
        //place_previous=place_next;
    }
    getline(fin,temp);
    fout << temp << "\n";
    //  cout << "temp is " << temp << "\n";
}

void Tools::GetNext2NumbersAsString(ifstream& fin,ofstream& fout,string& number1,string& number2, string& stherest) {

    string temp;  
    getline(fin,temp); // get next line of the input file 
 
    fout << temp << "\n"; // output this line to the summary file

    int place=0; 
    place=temp.find_first_of(" \t"); // find where the first space is

    number1=temp.substr(0,place); // everything up until the first space is what we're looking for

    temp=temp.substr(place+1,temp.size());
 

    number2=temp.substr(0,temp.find_first_of(" "));

    stherest=temp.substr(2,temp.size());
} //GetNext2NumbersAsString

double Tools::GetFWHM(TH1 *h1) {
    
    int imax=h1->GetMaximumBin();
    double max=h1->GetMaximum();

    //  cout << "imax, max are " << imax << " " << max << "\n";
    int ibin_plus=0;
    int ibin_minus=0;
    // now step to the right until it's half
    for (int ibin=imax;ibin<=h1->GetNbinsX();ibin++) {

        if (h1->GetBinContent(ibin)<max/2.) {
            ibin_plus=ibin;
            ibin=h1->GetNbinsX()+1;
            //  cout << "ibin_plus is " << ibin_plus << "\n";
        }
    }
    // now step to the left
    for (int ibin=imax;ibin>=1;ibin--) {

        if (h1->GetBinContent(ibin)<max/2.) {
            ibin_minus=ibin;
            ibin=0;
            //  cout << "ibin_minus is " << ibin_minus << "\n";
        }
    }
    if (ibin_plus>0 && ibin_minus==0) {
        ibin_minus=1;
        //cout << "bin_minus is " << ibin_minus << "\n";
    }
    

    if (ibin_plus==0 && ibin_minus==0) {
        cout << "Found 0 FWHM.\n";
        return 0.;
    }

    return (h1->GetBinCenter(ibin_plus)-h1->GetBinCenter(ibin_minus))/2.;

}

void Tools::Zero(int *anarray,int n) {
    for (int i=0;i<n;i++) {
        anarray[i]=0;
    } //for
} //Zero (int*,int)

void Tools::Zero(double *anarray,int n) {
    for (int i=0;i<n;i++) {
        anarray[i]=0.;
    } //for
} //Zero (int*,int)

double Tools::dMinNotZero(const double *x,int n) {
    double min=dMax(x,n);
    if (min==0)
        cout << "max is 0.\n";
    for (int k=1;k<n;k++) {
        if (x[k]<min && x[k]!=0)
            min=x[k];
    }
    return min;
} //dMinNotZero(double*, int)

double Tools::dMin(const double *x,int n) {
    double min=x[0];
    for (int k=1;k<n;k++) {
        if (x[k]<min)
            min=x[k];
    }
    return min;
} //dMin(double*, int)

double Tools::dMin(double x,double y) {
    double min=1.E22;
    if (x<y)
        min=x;
    else
        min=y;
    
    return min;
} //dMin(double,double)


double Tools::dMax(const double *x,int n) {
    
    double max=x[0];
    for (int k=1;k<n;k++) {
        if (x[k]>max)
            max=x[k];
    }
    return max;
} //dMax(double*, int)


double Tools::dvMax(const vector<double> x) {
    
    double max=x[0];
    for (int k=1;k<(int)x.size();k++) {
        if (x[k]>max)
            max=x[k];
    }
    return max;
} //dMax(double*, int)

double Tools::dsMax(TSpline5 *sp) {
    vector<double> y;
    double maxn;
    double blah1,blah2;
    for (int i=0;i<sp->GetNp();i++) {
        sp->GetKnot(i,blah1,blah2);
        y.push_back(blah2);
    }
    maxn=Tools::dvMax(y);
    return maxn;
}

double Tools::dMax(double a,double b) {
    if (a>b)
        return a;
    else if (a<b)
        return b;
    else if (a==b)
        return a;
    return 0;
} //dMax(double,double

int Tools::Getifreq(double freq,double freq_low,double freq_high,int n) {

    if (freq>=freq_high)
        return -1;
    if (freq<freq_low)
        return -1;

    return (int)((freq-freq_low)/(freq_high-freq_low)*(double)n);
} //Getifreq

void Tools::InterpolateComplex(double *array, const int n) {
    // to get rid of the zero bins
    double previous_nonzero=0.;
    double next_nonzero=0.;
    double check;
    int ifirstnonzero=0;
    int ilastnonzero=0;
    int k;
    int m=0;
    int count_nonzero=0;

    // find the first nonzero even element
    while (array[2*m]==0) {
        m++;
    }
    ifirstnonzero=m;
    
    
    // count the nonzero elements
    for (int i=0;i<n;i++) {
        if (array[2*i]!=0)
            count_nonzero++;
    }
    if (count_nonzero!=0) {

        // loop through the elements of the array and replace the zeros with interpolated values
        for (int i=ifirstnonzero;i<n;i++) {
            
            if (array[2*i]!=0.) {
                // set the lower nonzero value that we are interpolating from 
                previous_nonzero=array[2*i];
            
            }
            else {
                check=0.;
                k=i;
                while (check==0. && k<n) {
                    check=array[2*k];
                    k++;
                }
                if (k<n) {
                    next_nonzero=check;
            
                    for (int j=i;j<k;j++) {
                        array[2*j]=previous_nonzero+(next_nonzero-previous_nonzero)*(double)(j-(i-1))/(double)(k-(i-1));
                        array[2*j+1]=array[2*j];
                    }
                    i=k-1;
                
                    previous_nonzero=next_nonzero;
                }
                else {
                    ilastnonzero=i-1;
                    i=n;
                }
            } // end if array=0
        } // end loop over i
        
        //if (inu==49416)
        //cout << "inu, count_nonzero, diff are " << inu << " " << count_nonzero << " " << ilastnonzero << " " << ifirstnonzero << "\n";
        //cout << "factor is " << (double)count_nonzero/(double)(ilastnonzero-ifirstnonzero) << "\n";
        for (int j=0;j<n;j++) {
            array[2*j]*=sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
            array[2*j+1]*=sqrt((double)count_nonzero/(double)(ilastnonzero-ifirstnonzero));
        }
    }
}

void Tools::NormalTimeOrdering(const int n,double *volts) {
    double volts_temp[n];
    for (int i=0;i<n/2;i++) {
        volts_temp[i]=volts[i+n/2];
        volts_temp[i+n/2]=volts[i];
    }
    for (int i=0;i<n;i++) {
        volts[i]=volts_temp[i];
    }
}

void Tools::NormalTimeOrdering_InvT(const int n,double *volts) {
    double volts_temp[n];
    for (int i=0;i<n/2;i++) {
        volts_temp[i]=volts[i+n/2];
        volts_temp[i+n/2]=volts[i];
    }
    for (int i=0;i<n;i++) {
        volts[i]=volts_temp[n-i-1]; // inverse time
    }

}

//! A function to do sinc interpolation from time basis of x1 to the time basis of x2
/*!
    
    The function takes an input array (defined by n1, x1, y1), and interpolates
    that data to a new time base, provided by x2, and puts the values into y2.
    The user must therefore provide the number of input samples (n1)
    and the x and y values of the data to be interpolated (x1, y1).
    The user must also provide the number of points at which they would like
    the function interpolated (n2) and the x-values where the function is to be
    interpolation (x2). The content of y2[i] will be replaced with the interpolated values.

    \param n1 number of points in the input array
    \param x1 array of points representing the x-values of the input array
    \param y1 array of points representing the y-values of the input array
    \param n2 number of points in the output array
    \param x1 array of points representing the x-values of the output array
    \param y1 array of points representing the y-values of the output array
    \return void
*/

void Tools::SincInterpolation(int n1, double *x1, double *y1, int n2, double *x2, double *y2){

    /*
    * The Whittaker-Shannon interpolator is useful in the case of band-limited data.
    * Otherwise known as "sinc" interpolation, it protects the fidelity of the frequency spectrum of the signal.
    * Unlike, say, cubic-spline interpolation--which is faster, but can leave artifacts.
    * See https://en.wikipedia.org/wiki/Whittaker–Shannon_interpolation_formula for information,
    * and https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/whittaker_shannon.html
    * for implementation details from the boost documentation.
    * This method is slower than linear or spline interpolation--so its use is probably 
    * probably not ideal/necessary in cases where preserving spectral shape is not important.
    */
    
    // the whittaker-shannon method likes the data to be in a vector
    size_t num_input_samps = n1;
    std::vector<double> input_y(num_input_samps);
    for(size_t samp=0; samp<num_input_samps; samp++){
        input_y[samp] = y1[samp];
    }
    double t0 = x1[0];
    double dT = x1[1]-x1[0];
    double first_input_sample = x1[0];
    double last_input_sample = x1[n1-1];

    auto interpolator = boost::math::interpolators::whittaker_shannon<std::vector<double>>(std::move(input_y), t0, dT);

    for(int samp=0; samp<n2; samp++){
        // check if the sample comes before the first sample of the input array (x1[0])
        // or after the last sample of the input array (x1[n1-1])
        // if so, then we are asking for the function to *extrapolate*, not *interpolate*
        // just use the first/last sample, which replicates the behavior in SimpleLinearInterpolation_OutZero
        
        if(x2[samp]<first_input_sample){
            // before first sample, use first sample y1
            y2[samp] = y1[0];
        }
        else if(x2[samp]>last_input_sample){
            // after last sample, use last sample of y1
            y2[samp] = y1[n1-1];
        }
        else{
            // in the range of support, do interpolation
            y2[samp] = interpolator(x2[samp]);
        }
    }
}

void Tools::SimpleLinearInterpolation(int n1, double *x1, double *y1, int n2, double *x2, double *y2 ) {    // reads n1 array values x1, y1 and do simple linear interpolation and return n2 array with values x2, y2.
    // if interploated array has wider range (x values) than original array, it will use the first original value
    //

    int first=0;    // first x2 array which is bigger than x1[0]
    int last=0;       // last x2 array which is smaller than x1[n1-1]

    int cnt = 0;
    
    for (int i=0; i<n2; i++) {
            if (x2[i] < x1[0] ) {
                    first++;
            }
            else if (x2[i] > x1[n1-1]) {
                    last++;
            }
    }


    for (int i=0; i<n2; i++) {

        if (i<first) {  // if x2 has smaller x values than x1, just use x1[0] value
            y2[i] = y1[0];
        }
        else if (i>n2-last-1) {   // if x2 has bigger x values than x1, just use x1[n1-1] value
            y2[i] = y1[n1-1];
        }

        else {
            cnt=-1;
            for (int j=0; j<n1; j++) {
                //if (x2[i] < x1[j] && cnt==-1) {
                if (x2[i] <= x1[j] && cnt==-1) {
                    cnt = j;
                }
            }

            y2[i] = y1[cnt-1] + (x2[i]-x1[cnt-1])*(y1[cnt]-y1[cnt-1])/(x1[cnt]-x1[cnt-1]);
        }
    }
}



void Tools::SimpleLinearInterpolation_OutZero(int n1, double *x1, double *y1, int n2, double *x2, double *y2 ) {    // reads n1 array values x1, y1 and do simple linear interpolation and return n2 array with values x2, y2.
    // if interploated array has wider range (x values) than original array, it will use the first original value
    //

    int first=0;    // first x2 array which is bigger than x1[0]
    int last=0;       // last x2 array which is smaller than x1[n1-1]

    int cnt = 0;
    
    for (int i=0; i<n2; i++) {
        if (x2[i] < x1[0] ) {
                first++;
        }
        else if (x2[i] > x1[n1-1]) {
                last++;
        }
    }


    for (int i=0; i<n2; i++) {

        if (i<first) {  // if x2 has smaller x values than x1, just use x1[0] value
                //y2[i] = y1[0];
                y2[i] = 0.;
        }
        else if (i>n2-last-1) {   // if x2 has bigger x values than x1, just use x1[n1-1] value
                //y2[i] = y1[n1-1];
                y2[i] = 0.;
        }

        else {
            cnt=-1;
            //for (int j=0; j<n1; j++) {
            for (int j=1; j<n1; j++) {
                //if (x2[i] < x1[j] && cnt==-1) {
                if (x2[i] <= x1[j] && cnt==-1) {
                    cnt = j;
                }
            }

            y2[i] = y1[cnt-1] + (x2[i]-x1[cnt-1])*(y1[cnt]-y1[cnt-1])/(x1[cnt]-x1[cnt-1]);

            // test
            // if y2[i] goes outside y1[cnt] and y1[cnt-1] range (which should not)
            if ( ( y2[i] > y1[cnt] && y2[i] > y1[cnt-1] ) || 
                ( y2[i] < y1[cnt] && y2[i] < y1[cnt-1] ) ) {
                    cout<<"SimpleLinearInterpolation bug?! y2["<<i<<"] : "<<y2[i]<<" where y1["<<cnt<<"] : "<<y1[cnt]<<", y1["<<cnt-1<<"] : "<<y1[cnt-1]<<std::endl;
            }
        }
    }
}

double Tools::SimpleLinearInterpolation_extend_Single(int n1, double *x1, double *y1, double x2 ) {    // reads n1 array values x1, y1 and do simple linear interpolation and return y2 value

    int first=0;    // first x2 array which is bigger than x1[0]
    int last=0;       // last x2 array which is smaller than x1[n1-1]

    int cnt = 0;
    
    double slope_1;
    double slope_2;

    slope_1 = (y1[1] - y1[0]) / (x1[1] - x1[0]);
    slope_2 = (y1[n1-1] - y1[n1-2]) / (x1[n1-1] - x1[n1-2]);

    double y2;
            
    if ( x2 < x1[0] ) {  // if x2 has smaller x values      
        y2 = slope_1 * (x2 - x1[0]) + y1[0];
    }
    else if ( x2 > x1[n1-1] ) {         
        y2 = slope_2 * (x2 - x1[n1-1]) + y1[n1-1];
    }
    else { // in between

        cnt=-1;
        for (int j=0; j<n1; j++) {
            //if (x2[i] < x1[j] && cnt==-1) {
            if (x2 <= x1[j] && cnt==-1) {
                cnt = j;
            }
        }

        y2 = y1[cnt-1] + (x2-x1[cnt-1])*(y1[cnt]-y1[cnt-1])/(x1[cnt]-x1[cnt-1]);            
    }

    return y2;

}

void Tools::get_random_rician(double signal_amplitude, double signal_phase, double sigma, double &amplitude, double &phase){
    double rand_gauss_a, rand_gauss_b;
    get_circular_bivariate_normal_random_variable(rand_gauss_a, rand_gauss_b);

    // check the value
    //cout<<"Factor from random rician : "<<sqrt(rand_gauss_a * rand_gauss_a +rand_gauss_b *rand_gauss_b) * sqrt(2./M_PI)<<"\n";
    
    // Gives the gaussian-distributed random variables a standard deviation of sigma
    rand_gauss_a *= sigma;
    rand_gauss_b *= sigma;
    
    // Gives the gaussian-distributed random variables a mean of (v*cos(theta), v*sin(theta)) when v is the mean of the desired rician distribution
    rand_gauss_a += signal_amplitude * cos(signal_phase);
    rand_gauss_b += signal_amplitude * sin(signal_phase);
    
    // The Rician Distribution produces the probability of the the absolute value (radius) of a circular bivariate normal random variable:
    amplitude = sqrt(rand_gauss_a * rand_gauss_a + rand_gauss_b * rand_gauss_b);
    // Thus, the descriptor other than amplitude for the circular bivariate is given by a phase:
    phase = atan2(rand_gauss_b, rand_gauss_a);
    return;
}

void Tools::get_circular_bivariate_normal_random_variable(double& rand_gauss_a, double& rand_gauss_b){
    double rand_uni_a = gRandom->Rndm(); //gRandom->Rndm() produces uniformly-distributed floating points in ]0,1]
    double rand_uni_b = gRandom->Rndm();
        
    // Box-Muller transform from a bivariate uniform distribution from 0 to 1 to a gaussian with mean = 0 and sigma = 1
    rand_gauss_a = sqrt(-2. * log(rand_uni_a)) * cos(2. * M_PI * rand_uni_b);
    rand_gauss_b = sqrt(-2. * log(rand_uni_a)) * sin(2. * M_PI * rand_uni_b);
    return;
}


void Tools::Exchange( double &a, double &b ) {
    double tmp = a;
    a = b;
    b = tmp;
}

