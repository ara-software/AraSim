#!/usr/bin/env python

#EXTREMELY HACKY example of using AraSim raytracer

accuracy = 0.05

import ctypes
import sys


import ROOT



import numpy as np
import matplotlib.pylab as plt


ROOT.gSystem.Load("libAra.so") 

def R(cmd):
  return ROOT.gInterpreter.ProcessLine(cmd)


R('#include "RayTrace.h"')
R('#include "RayTrace_IceModels.h"')



# attenuation model 
R('auto atten_model = boost::shared_ptr<basicAttenuationModel>( new basicAttenuationModel );'); 

#exponential model 
R('auto refr_model = boost::shared_ptr<exponentialRefractiveIndex>(new exponentialRefractiveIndex(1.35,1.78,0.0132));')


R('RayTrace::TraceFinder tf(refr_model, atten_model);')
R('Vector src; Vector rec; std::vector<RayTrace::TraceRecord> paths;')


#add additional parameters to dictionary if needed
def trace(r_src, z_src, z_rec=-200): 

  ROOT.src.SetXYZ(r_src,0,z_src); 
  ROOT.rec.SetXYZ(0,0,z_rec); 
  sol_cnt = ctypes.c_int()
  sol_err = ctypes.c_int()
  ROOT.paths = ROOT.tf.findPaths(ROOT.src, ROOT.rec, 0.3, ROOT.TMath.Pi()/2, sol_cnt, sol_err, ROOT.RayTrace.SurfaceReflection,accuracy)

  times = []
  for path in ROOT.paths: 
    times.append({'tof':path.pathTime*1e9,'dist': path.pathLen, 'atten': path.attenuation})
  return times 


nr = 501; 
rmax = 5000
rmin = 0

nz = 251
zmax = -2500
zmin = 0

t0s = np.zeros((nr,nz)) 
t1s = np.zeros((nr,nz))


jj = 0

for z in np.linspace(zmin,zmax,nz):
  ii = 0
  for r in np.linspace(rmin,rmax,nr): 
    results = trace(r,z)
    t0s[ii][jj] = results[0]["tof"] if len(results) else -1
    t1s[ii][jj] = results[1]["tof"] if len(results) > 1 else -1
    ii+=1
  jj+=1

plt.figure("Direct") 
plt.imshow(t0s.T, extent=(rmin,rmax,zmax,zmin), origin='upper')
plt.figure("Reflected") 
plt.imshow(t1s.T, extent=(rmin,rmax,zmax,zmin), origin='upper')
plt.show()






