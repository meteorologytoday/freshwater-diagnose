#!/bin/bash

ncra -O -v U,V,T,lev,ilev,lat,slat atm/hist/paper2021_CTL_POP2.cam.h0.00{41..45}-{01..12}.nc avg.nc
ncwa -O -a lon,time avg.nc avgzm.nc
