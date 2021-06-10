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
    z_T = f.variables["z_T"][:] / 1e3
    y_V = f.variables["y_V"][:]
    z_W = f.variables["z_W"][:] / 1e3
    
    psi1   = f.variables["PSI1"][0, :, :].transpose()
    psi2F  = f.variables["PSI2F"][0, :, :].transpose()
    psi2QE = f.variables["PSI2QE"][0, :, :].transpose()
    Q1   = f.variables["Q1"][0, :, :].transpose()
    F2   = f.variables["F2"][0, :, :].transpose()
    QE2  = f.variables["QE2"][0, :, :].transpose()





cmap1 = cm.get_cmap("bwr")
cflevs1 = np.linspace(-10, 10, 21)
clevs1 = np.linspace(-100, 100, 101)

cmap2 = cm.get_cmap("bwr")
cflevs2 = np.linspace(-100, 100, 21)
clevs2 = np.linspace(-10, 10, 41)

cmap3 = cm.get_cmap("bwr")
cflevs3 = np.linspace(-200, 200, 21)
clevs3 = np.linspace(-100, 100, 41)


xticks = [-90, -60, -30, 0, 30, 60, 90]

fig, ax = plt.subplots(3, 1, sharex=True, figsize=(12, 12), constrained_layout=True)

for i, _ax in enumerate(ax):

    _ax.set_ylim([0, 15])
    _ax.set_xticks(xticks)
    _ax.grid(True)
    _ax.set_ylabel("Height [ km ]")

    
mappable1 = ax[0].contourf(y_T, z_T, Q1 * 86400, cflevs1, cmap=cmap1, extend="both")
cb = plt.colorbar(mappable1, ax=ax[0], orientation="vertical")
cb.ax.set_ylabel(r"Heating rate [$\mathrm{K} / \mathrm{day}$]")
cl = ax[0].contour(y_V, z_W, psi1 / 1e10, clevs1, colors='k'); ax[0].clabel(cl, fmt="%.1f")
ax[0].set_title(r"Imaginary adiabatic heating and diagnosed streamfunction $\psi$")



mappable2 = ax[1].contourf(y_V, z_T, F2, cflevs2, cmap=cmap2, extend="both")
cb = plt.colorbar(mappable2, ax=ax[1], orientation="vertical")
cb.ax.set_ylabel(r"$\overline{v'u'}$ [$\mathrm{m}^2 / \mathrm{s}^2$]")
cl = ax[1].contour(y_V, z_W, psi2F / 1e10, clevs2, colors='k'); ax[1].clabel(cl, fmt="%.1f")
ax[1].set_title(r"Simulated eddy momentum flux $\overline{v'u'}$ and diagnosed streamfunction $\psi$")

mappable3 = ax[2].contourf(y_V, z_T, QE2, cflevs3, cmap=cmap3, extend="both")
cb = plt.colorbar(mappable3, ax=ax[2], orientation="vertical")
cb.ax.set_ylabel(r"$\overline{v'\theta'} $[$\mathrm{K} \, \mathrm{m} / \mathrm{s}$]")
cl = ax[2].contour(y_V, z_W, psi2QE / 1e10, clevs3, colors='k'); ax[2].clabel(cl, fmt="%.1f")
ax[2].set_title(r"Simulated eddy heat flux $\overline{v'\theta'}$ and diagnosed streamfunction $\psi$")



ax[2].set_xlabel(r"Latitude [ degree north ]")

plt.show()

fig.savefig("Example_invert.png", dpi=200)
