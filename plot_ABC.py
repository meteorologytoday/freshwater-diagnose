import matplotlib as mplt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib import rc



default_linewidth = 1.0;
default_ticksize = 10.0;

mplt.rcParams['lines.linewidth'] =   default_linewidth;
mplt.rcParams['axes.linewidth'] =    default_linewidth;
mplt.rcParams['xtick.major.size'] =  default_ticksize;
mplt.rcParams['xtick.major.width'] = default_linewidth;
mplt.rcParams['ytick.major.size'] =  default_ticksize;
mplt.rcParams['ytick.major.width'] = default_linewidth;

rc('font', **{'size': 15.0});
rc('axes', **{'labelsize': 15.0});
rc('mathtext', **{'fontset':'stixsans'});

import matplotlib.pyplot as plt

import os,sys, argparse
from netCDF4 import Dataset
import numpy as np
from pprint import pprint
import sys, argparse

with Dataset("output.nc", "r") as f:
    y_T = f.variables["y_T"][:]
    z_T = f.variables["z_T"][:]
    y_V = f.variables["y_V"][:]
    z_W = f.variables["z_W"][:]
    
    psi = f.variables["psi"][0, :, :].transpose()
    A_W = f.variables["A_W"][0, :, :].transpose()
    theta_T = f.variables["theta_T"][0, :, :].transpose()
    B_W = f.variables["B_W"][0, :, :].transpose()
    B_V = f.variables["B_V"][0, :, :].transpose()
    C_V = f.variables["C_V"][0, :, :].transpose()


cmap1 = cm.get_cmap("gnuplot")
cmap2 = cm.get_cmap("bwr")

clevs1 = np.linspace(240, 350, 11)

factor=1e-4
A_levs = np.linspace(-5, 5, 11)



xticks = [-90, -60, -30, 0, 30, 60, 90]

fig = plt.figure(figsize=(12, 6))

hr = [1,1,1]
wr = [1,0.1]
gs = fig.add_gridspec(nrows=len(hr), ncols=len(wr), width_ratios=wr, height_ratios=hr, hspace=0.3)

axes = []
for i in range(len(hr)):
    axes.append(fig.add_subplot(gs[i, 0]))

cax = fig.add_subplot(gs[:, 1])

for i, _ax in enumerate(axes):

    mappable1 = _ax.contourf(y_T, z_T, theta_T, clevs1, cmap=cmap1, extend="both")

    _ax.set_ylim([0, 15000.0])
    _ax.set_xticks(xticks)
    if i != len(axes)-1:
        _ax.set_xticklabels([""] * len(xticks))

cb = plt.colorbar(mappable1, cax=cax, orientation="vertical")

cl = axes[0].contour(y_T, z_W, A_W / factor, A_levs, colors='k'); axes[0].clabel(cl, fmt="%.1f")
cl = axes[1].contour(y_V, z_T, B_V / factor, A_levs, colors='k'); axes[1].clabel(cl, fmt="%.1f")
cl = axes[2].contour(y_V, z_T, C_V / factor, A_levs, colors='k'); axes[2].clabel(cl, fmt="%.1f")


plt.show()

