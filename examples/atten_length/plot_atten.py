#!/usr/bin/env python
"""Plot ice attenuation length and temperature vs depth from CSV."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("atten_vs_depth.csv")
depths_m = df["depth_m"].values
temperatures_C = df["temperature_C"].values

freq_cols = [c for c in df.columns if c.startswith("L_") and c.endswith("MHz_m")]
frequencies_mhz = [float(c.split("_")[1].replace("MHz", "")) for c in freq_cols]

linestyles = ["-", "--", ":", "-."]
line_color = "tab:blue"
temp_color = "tab:red"

fig, ax = plt.subplots(figsize=(7.0, 6.5))

for col, f_mhz, ls in zip(freq_cols, frequencies_mhz, linestyles):
    ax.plot(df[col].values, depths_m,
            color=line_color, ls=ls, lw=1.8, label=f"{int(f_mhz)} MHz")

ax.set_ylabel("Depth (m)")
ax.set_xlabel("Attenuation length (m)", color=line_color)
ax.tick_params(axis="x", colors=line_color)
ax.spines["bottom"].set_color(line_color)
ax.invert_yaxis()
ax.grid(True, alpha=0.3)

south_pole_ice_depth_m = 2850.0
ax.axhline(south_pole_ice_depth_m, color="k", lw=1.0, ls="-.", alpha=0.7)
ax.text(ax.get_xlim()[1], south_pole_ice_depth_m,
        f" South Pole bedrock ({south_pole_ice_depth_m:.0f} m)",
        color="k", fontsize=8, va="bottom", ha="right")

ax_T = ax.twiny()
ax_T.plot(temperatures_C, depths_m, color=temp_color, lw=1.8, label="Temperature")
ax_T.set_xlabel("Temperature (\N{DEGREE SIGN}C)", color=temp_color)
ax_T.tick_params(axis="x", colors=temp_color)
ax_T.spines["top"].set_color(temp_color)
ax_T.spines["bottom"].set_visible(False)

lines_L, labels_L = ax.get_legend_handles_labels()
lines_T, labels_T = ax_T.get_legend_handles_labels()
ax.legend(lines_L + lines_T, labels_L + labels_T, loc="center right")

ax.set_title("AraSim ice attenuation length and temperature vs depth", pad=28)

fig.tight_layout()
fig.savefig("atten_vs_depth.png", dpi=150)
print("Wrote atten_vs_depth.png")
