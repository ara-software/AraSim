# Updated date: 2020/12/22
#
# This script is designed to make the Chiba in-situ antenna model.
# The model is based on two polarization (Vpol and Hpol)
# Each polarization model is a 2D array (Theta angle, frequency range).
# The unit of resulted gain is [dB].
# Additional script parts that makes the gain into AraSim antenna format will be added here soon.
#
# The related talk:http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1838
# The related wiki page:http://ara.icecube.wisc.edu/wiki/index.php/Antenna_model
#
# Contacts info.: Myoungchul Kim, mkim@icecube.wisc.edu

import os, sys
import numpy as np
from matplotlib import pyplot as plt
import h5py

def VPol_Model(freq, angle): # from Simon Archambault
    A = 3.97624997308*pow(freq, 0)  - 315.233620075*pow(freq, 1) + 5051.87095934*pow(freq, 2) - 35817.642524*pow(freq, 3) + 140907.450596*pow(freq, 4) - 329031.929939*pow(freq, 5) + 453097.902201*pow(freq, 6) - 339172.365085*pow(freq, 7) + 106162.954882*pow(freq, 8)
    #A = np.array(A, dtype=np.double)
    B = -0.980440058135*pow(freq, 0) + 19.708692819*pow(freq, 1) - 58.9056647312*pow(freq, 2) + 52.76218623105*pow(freq, 3) + 0.996765601448*pow(freq, 4)
    #B = np.array(B, dtype=np.double)
    C = 0.79019478805*pow(freq, 0) + 1.48169370259*pow(freq, 1) - 55.7820913455*pow(freq, 2) + 171.131052418*pow(freq, 3) - 133.907945563*pow(freq, 4)
    #C = np.array(C, dtype=np.double)
    D = 1.98220468619*pow(freq, 0) + 0.23451969725*pow(freq, 1) - 0.944977363565*pow(freq, 2) + 1.52844938585*pow(freq,3) - 0.862249036112*pow(freq,4)
    #D = np.array(D, dtype=np.double)
    rgain = A*pow(np.cos(B*angle), 2) + C*pow(np.sin(D*angle), 2)

    rgain[freq<0.150] = 1.e-10
    rgain[freq>0.700] = 1.e-10
    rgain[rgain<=0] = 1.e-10
    return rgain

def HPol_Model(freq, angle): # from Simon Archambault
    A = 73.7257*pow(freq, 0) - 1156.5*pow(freq,1) + 7178.79*pow(freq, 2) - 22526.7*pow(freq, 3) + 37982.6*pow(freq, 4) - 32721.5*pow(freq, 5) + 11269.1*pow(freq, 6)
    B = 1.415*pow(freq, 0) + 1.1625*pow(freq, 1) - 16.436*pow(freq, 2) + 24.945*pow(freq, 3)
    C = -22.5365*pow(freq, 0) + 366.68*pow(freq, 1) - 2360.94*pow(freq, 2) + 7715.72*pow(freq, 3) - 13517.22*pow(freq, 4) + 12046.69*pow(freq, 5) - 4277.36*pow(freq, 6)
    D = 1.9527*pow(freq, 0) + 0.3485*pow(freq, 1) - 0.7970*pow(freq, 2) + 0.5743*pow(freq, 3)
    rgain = A*pow(np.cos(B*angle), 2) + C*pow(np.sin(D*angle), 2)
    rgain[freq<=0.200] = 1.e-10
    rgain[freq>0.700] = 1.e-10
    rgain[rgain<=0] = 1.e-10

    if len(angle) > 1:
        for ii in range(len(freq)):
            if (-0.3*(1000*freq[ii]-800)) < np.abs(180.*angle[ii]/np.pi):
                rgain[ii] = 1.e-10
    else:
        for ii in range(len(freq)):
            if (-0.3*(1000*freq[ii]-800)) < np.abs(180.*angle/np.pi):
                rgain[ii] = 1.e-10

    return rgain

def model_2d_plot(gains, title, file_name, freq_range, freq_bin, theta_range, theta_bin):
    
    fig, ax = plt.subplots(figsize=(9, 6))
    plt.xlabel(r'Elevation Angle [ $deg.$ ]', fontsize=25)
    plt.ylabel(r'Frequency [ $MHz$ ]', fontsize=25)
    plt.grid(color='k',linestyle=':')
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)
    plt.title(title, y=1.02,fontsize=25)
    #plt.ylim(freq_range[0]-freq_bin/2,freq_range[-1]+freq_bin/2)
    plt.ylim(150,800)
    plt.xlim(theta_range[0]-theta_bin/2,theta_range[-1]+theta_bin/2)

    cc = plt.imshow(gains, interpolation='none',cmap='viridis',extent=[theta_range.min()-theta_bin/2, theta_range.max()+theta_bin/2, freq_range.max()+freq_bin/2, freq_range.min()-freq_bin/2], aspect='auto')

    cbar1 = plt.colorbar(cc, ax=ax)
    cbar1.ax.tick_params(axis='y', labelsize=15)
    cbar1.ax.set_ylabel(r'Gain [ $dB$ ]', fontsize=15)#, color= 'red')

    plt.tight_layout()
    #r_path = '/path/to/somewhere/'
    #os.chdir(r_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close()

def main():

    # frequency range. For now, it is following AraSim table frequency range 
    start_freq_0=0
    end_freq=1066
    interval_freq=100/6
    freq = np.arange(start_freq_0,end_freq+1,interval_freq)[5:]

    # theta range
    theta_bin = 5
    theta = np.arange(-90,90.001,theta_bin)
    theta_rad = np.radians(theta)

    # array for gain results
    Gain = np.full((2,len(theta),len(freq)),np.nan)

    # making the model
    for angle in range(len(theta)):

        chosen_theta_rad = np.array([theta_rad[angle]]) # need to cooperate with the antenna model function.

        Gain[0,angle] = VPol_Model(freq/1000, chosen_theta_rad)
        Gain[1,angle] = HPol_Model(freq/1000, chosen_theta_rad)
    
    print('The model making is done!')

    #saving
    hf = h5py.File('In_situ_Vpol_Model.h5', 'w')
    hf.create_dataset('theta_range_deg', data=theta, compression="gzip", compression_opts=9)
    hf.create_dataset('freq_range_MHz', data=freq, compression="gzip", compression_opts=9)
    hf.create_dataset('gain_hpol_dB', data=Gain[0], compression="gzip", compression_opts=9)
    hf.close()

    hf = h5py.File('In_situ_Hpol_Model.h5', 'w')
    hf.create_dataset('theta_range_deg', data=theta, compression="gzip", compression_opts=9)
    hf.create_dataset('freq_range_MHz', data=freq, compression="gzip", compression_opts=9)
    hf.create_dataset('gain_hpol_dB', data=Gain[1], compression="gzip", compression_opts=9)
    hf.close()

    #vpol plot
    model_2d_plot(np.copy(Gain[0]).T, r'Chiba In-situ Vpol Model', r'In_situ_Vpol_Model.png', freq, interval_freq, theta, theta_bin)

    #hpol plot
    model_2d_plot(np.copy(Gain[1]).T, r'Chiba In-situ Hpol Model', r'In_situ_Hpol_Model.png', freq, interval_freq, theta, theta_bin)

    print('Done! The output h5 and png are in the same path.')

if __name__ == "__main__":

    main()
