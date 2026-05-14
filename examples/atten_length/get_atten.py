#!/usr/bin/env python
"""Compute AraSim ice attenuation length vs depth at a few frequencies."""

import numpy as np
import ROOT
from ROOT import gSystem

gSystem.Load('/home/baclark/scratch/ARA/AraSim/libAra.so')

icemodel = ROOT.IceModel(10, 0, 0)
settings = ROOT.Settings()

depths_m = np.linspace(1.0, 3000.0, 301)
frequencies_mhz = [100.0, 300.0, 700.0]

attens = {f: np.empty_like(depths_m) for f in frequencies_mhz}
for f_mhz in frequencies_mhz:
    f_ghz = f_mhz / 1000.0
    for i, z in enumerate(depths_m):
        attens[f_mhz][i] = icemodel.GetFreqDepIceAttenuLength(
            float(z), float(f_ghz), settings
        )

temperatures_C = np.array([icemodel.temperature(float(z), False) for z in depths_m])

csv_path = "atten_vs_depth.csv"
header = "depth_m,temperature_C," + ",".join(
    f"L_{int(f)}MHz_m" for f in frequencies_mhz
)
stacked = np.column_stack(
    [depths_m, temperatures_C] + [attens[f] for f in frequencies_mhz]
)
np.savetxt(csv_path, stacked, delimiter=",", header=header, comments="", fmt="%.4f")
print(f"Wrote {csv_path}")
