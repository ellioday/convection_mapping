#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 13:56:40 2021

@author: elliott
"""

import numpy as np
import pydatadarn
import datetime
import potential_functions as pf
import matplotlib.pyplot as plt

import fpipy

stations=["ade", "adw", "bks", "cly", "cve", "cvw", "fhe", "fhw", "gbr", "han",
		   "hkw", "hok", "inv", "kap", "kod", "ksr", "lyr", "pgr", "pyk", 
		   "rkn", "sas", "sto", "wal"]
FPI_stations=["ann", "uao"]

##########################
### Get SuperDARN data ###
##########################

vectors = pydatadarn.GridData(is_map=True)
vectors.add_data("all", "2013/10/02 07:00:00", "2013/10/02 08:00:00")

time_i = "2013/10/02 07:57:00"
time_index = np.where(vectors.times == time_i)

#get heppner-maynard boundary information
hmb_mlat = vectors.boundary_mlats[time_i]
hmb_mcolat = 90-hmb_mlat
hmb_mlon = vectors.boundary_mlons[time_i]

#get spherical fitting coefficients
N = vectors.N[time_i]
N1 = vectors.N1[time_i]
N2 = vectors.N2[time_i]
N3 = vectors.N3[time_i]
solution = np.zeros([4, len(N)])
solution[0,] = N
solution[1,] = N1
solution[2,] = N2
solution[3,] = N3

#get fitting parameters
latmin = vectors.latmin[time_i]
lon_shft = vectors.lon_shft[time_i]
lat_shft = vectors.lat_shft[time_i]
order = 8

#get data and model lat/lon
mlats = vectors.mlats[time_index]
mlons = vectors.mlons[time_index]
mod_mlats = vectors.mod_mlats[time_i]
mod_mlons = vectors.mod_mlons[time_i]
#append model to data
mlats = np.append(mlats, mod_mlats)
mlons = np.append(mlons, mod_mlons)
mcolats = 90-mlats
#affix lats/lons into one array for spherical fitting input
pos = np.array([mlats, mlons])

#get velocity and direction from spherical fit
mag, az = pf.find_gradV(pos, solution, latmin, lon_shft, lat_shft, order)


################
# Get FPI data #
################

fpi_path = "/home/elliott/Documents/madrigalWeb-3.2/madrigalWeb/Data/"
fname1 = fpi_path+"minime05_uao_20131001.cedar.010.hdf5"
fname2 = fpi_path+"minime05_uao_20131002.cedar.008.hdf5"

FPI = fpipy.FPIData(fname1)

FPI.add_HDF5(fname2)
FPI.get_table()

time_min = "2013/10/02 00:30:00"
time_max = "2013/10/02 08:30:00"
dtime_min = pydatadarn.tools.time_to_dtime(time_min)
dtime_max = pydatadarn.tools.time_to_dtime(time_max)

#convert time_i into seconds
dtime_i = pydatadarn.tools.time_to_dtime(time_i)
time_s = (dtime_i - dtime_min).seconds

N, E, S, W, zen, dN, dE, dS, dW, dzen, N_times, E_times, S_times, W_times, zen_times = FPI.get_azm_vels(
	dtime_min = dtime_min, dtime_max = dtime_max)

oN = N
oE = E
oS = S
oW = W

#get dtimes
N_dtimes = pydatadarn.tools.time_to_dtime(N_times)
E_dtimes = pydatadarn.tools.time_to_dtime(E_times)
S_dtimes = pydatadarn.tools.time_to_dtime(S_times)
W_dtimes = pydatadarn.tools.time_to_dtime(W_times)
zen_dtimes = pydatadarn.tools.time_to_dtime(zen_times)

#calculate number of seconds between start and end time
dtime_min = pydatadarn.tools.time_to_dtime(time_min)
dtime_max = pydatadarn.tools.time_to_dtime(time_max)
num_seconds = (dtime_max - dtime_min).seconds

#if we have an array of all seconds between the start and end time
#get the indexes where measurement times match seconds
N_time_indexes = np.empty(len(N), dtype="int")
E_time_indexes = np.empty(len(E), dtype="int")
S_time_indexes = np.empty(len(S), dtype="int")
W_time_indexes = np.empty(len(W), dtype="int")
zen_time_indexes = np.empty(len(zen), dtype="int")

for i in range(len(N)):
	N_time_indexes[i] = (N_dtimes[i]-dtime_min).total_seconds()
for i in range(len(E)):
	E_time_indexes[i] = (E_dtimes[i]-dtime_min).total_seconds()
for i in range(len(S)):
	S_time_indexes[i] = (S_dtimes[i]-dtime_min).total_seconds()
for i in range(len(W)):
	W_time_indexes[i] = (W_dtimes[i]-dtime_min).total_seconds()
for i in range(len(zen)):
	zen_time_indexes[i] = (zen_dtimes[i]-dtime_min).total_seconds()

N = fpipy.hor_vel_calc(N, N_time_indexes, zen, zen_time_indexes)
E = fpipy.hor_vel_calc(E, E_time_indexes, zen, zen_time_indexes)
S = fpipy.hor_vel_calc(S, S_time_indexes, zen, zen_time_indexes)
W = fpipy.hor_vel_calc(W, W_time_indexes, zen, zen_time_indexes)

#get data for specific time_i
N_time_i = fpipy.interpolate(N, N_time_indexes, time_s)
E_time_i = fpipy.interpolate(E, E_time_indexes, time_s)
S_time_i = fpipy.interpolate(S, S_time_indexes, time_s)
W_time_i = fpipy.interpolate(W, W_time_indexes, time_s)

#get East-West and North-South flow across the FPI station
EW = (E_time_i + W_time_i)/2
NS = (N_time_i + S_time_i)/2

#compute full 2D vector
FPI_vel, FPI_az = fpipy.merge_vecs(NS, EW)

tick_labels = np.array([])

time_tick = dtime_min
minute_interval = 60
while time_tick <= dtime_max:
	hour = time_tick.hour
	minute = time_tick.minute
	tick_labels = np.append(tick_labels, "{:02d}:{:02d}".format(hour, minute))
	time_tick = time_tick + datetime.timedelta(minutes=60)

xlabels = np.array([])
for i in range(0, 11):
	xlables = np.append(xlabels, "00:{:02d}".format(i))

fig, ax = plt.subplots(3, 1, sharex=True)
ax0 = ax[0]
ax1 = ax[1]
ax2 = ax[2]

ax0.errorbar(W_time_indexes, -W, yerr = dW, marker="o", markersize=2.5, label = "West Look")
ax0.errorbar(E_time_indexes, E, yerr = dE, marker="o", markersize=2.5, label = "East Look")
ax0.set_xticks(np.arange(0, num_seconds+minute_interval*60, minute_interval*60))
ax0.set_xticklabels(tick_labels)
#ax0.set_ylim([-150, 100])
ax0.axhline(y=0, color = "k", linestyle = "-.")
ax0.grid(linestyle = "--")
ax0.legend()

ax1.errorbar(S_time_indexes, -S, yerr = dS, marker="o", markersize=2.5, label = "South Look")
ax1.errorbar(N_time_indexes, N, yerr = dN, marker="o", markersize=2.5, label = "North Look")
ax1.axhline(y=0, color = "k", linestyle = "-.")
#ax1.set_ylim([-600, 150])
ax1.grid(linestyle = "--")
ax1.legend()

ax2.errorbar(zen_time_indexes, zen, yerr=dzen, marker="o", markersize=2.5, label = "Zenith")
ax2.axhline(y=0, color = "k", linestyle = "-.")
ax2.set_yticklabels([-600, -450, -300, -150, 0, 150])
ax2.grid(linestyle = "--")
ax2.legend()

plt.show()

#testing
UAO = fpipy.FPIStation("uao")
ANN = fpipy.FPIStation("ann")

uao_mlat, uao_mlon = UAO.get_aacgm(dtime_i)
uao_mcolat = 90-uao_mlat

from pydatadarn.utils import tools as tools
FPI_dr, FPI_dtheta = tools.vector_change(uao_mcolat, uao_mlon, FPI_vel, FPI_az)

#plot
dr, dtheta = pydatadarn.plotting.vector_plot(mcolats, mlons, az, mag, time=time_i, 
								station_names=stations, FPI_names=FPI_stations, 
								FPI_kvecs=FPI_az, FPI_vels=FPI_vel, mlt=True, 
								theta_min=0, theta_max=90, cbar_min=min(mag), 
								cbar_max=max(mag))

