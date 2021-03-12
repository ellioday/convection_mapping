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
import aacgmv2

stations=["ade", "adw", "bks", "cly", "cve", "cvw", "fhe", "fhw", "gbr", "han",
		   "hkw", "hok", "inv", "kap", "kod", "ksr", "lyr", "pgr", "pyk", 
		   "rkn", "sas", "sto", "wal"]
FPI_stations=["ann", "uao"]

##########################
### Get SuperDARN data ###
##########################

vectors = pydatadarn.GridData(is_map=True)
vectors.add_data("all", "2013/10/02 00:00:00", "2013/10/02 09:00:00")

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

#get data 
N, E, S, W, zen, dN, dE, dS, dW, dzen, N_times, E_times, S_times, W_times, zen_times = FPI.get_azm_vels(
	dtime_min = dtime_min, dtime_max = dtime_max)

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

#get times in seconds from start time instead of as dtime objects
N_time_indexes = fpipy.times_to_secs(N_dtimes, dtime_min)
E_time_indexes = fpipy.times_to_secs(E_dtimes, dtime_min)
S_time_indexes = fpipy.times_to_secs(S_dtimes, dtime_min)
W_time_indexes = fpipy.times_to_secs(W_dtimes, dtime_min)
zen_time_indexes = fpipy.times_to_secs(zen_dtimes, dtime_min)

#calculate horizontal component of line of sight velocities
hN = fpipy.hor_vel_calc(N, N_time_indexes, zen, zen_time_indexes)
hE = fpipy.hor_vel_calc(E, E_time_indexes, zen, zen_time_indexes)
hS = fpipy.hor_vel_calc(S, S_time_indexes, zen, zen_time_indexes)
hW = fpipy.hor_vel_calc(W, W_time_indexes, zen, zen_time_indexes)
#South and west need to flip signs so that neutral wind direction is consistent
#Either side of FPI (so -ve = Westwards/Southwards not -ve=towards FPI)
hS = -hS
hW = -hW

#######################
### Get Time sample ###
#######################

neutral_vel_interval = np.array([])
neutral_az_interval = np.array([])
plasma_vel_interval = np.array([])
plasma_az_interval = np.array([])
latmins = np.array([])
times = np.array([])

hmin=3
hmax=8
min_interval=10

loop=True
if loop == True:
	for hour in range(hmin, hmax+1):
		for minute in range(0, 60, min_interval):
			time_i = "2013/10/02 {:02d}:{:02d}:00".format(hour, minute)
			print(time_i)
			time_index = np.where(vectors.times == time_i)
			
			#convert time_i into seconds
			dtime_i = pydatadarn.tools.time_to_dtime(time_i)
			time_s = (dtime_i - dtime_min).seconds
			
			##########################
			### Get superDARN info ###
			##########################
			
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
			plasma_mlats = vectors.mlats[time_index]
			plasma_mcolats = 90-plasma_mlats
			plasma_mlons = vectors.mlons[time_index]
			mod_mlats = vectors.mod_mlats[time_i]
			mod_mlons = vectors.mod_mlons[time_i]
			#append model to data
			all_mlats = np.append(plasma_mlats, mod_mlats)
			all_mlons = np.append(plasma_mlons, mod_mlons)
			all_mcolats = 90-all_mlats
			#affix lats/lons into one array for spherical fitting input
			pos = np.array([plasma_mlats, plasma_mlons])
			
			#get velocity and direction from spherical fit
			plasma_vels, plasma_az = pf.find_gradV(pos, solution, latmin, lon_shft, lat_shft, order, dtime_i)
			
			####################
			### Get FPI Data ###
			####################
			
			#get data for specific time_i
			N_time_i = fpipy.interpolate(hN, N_time_indexes, time_s)
			E_time_i = fpipy.interpolate(hE, E_time_indexes, time_s)
			S_time_i = fpipy.interpolate(hS, S_time_indexes, time_s)
			W_time_i = fpipy.interpolate(hW, W_time_indexes, time_s)
			
			#get East-West and North-South flow across the FPI station
			hEW = (E_time_i + W_time_i)/2
			hNS = (N_time_i + S_time_i)/2
			
			#compute full 2D vector for neutral wind
			uao_neutral_vel, uao_neutral_az = fpipy.merge_vecs(hNS, hEW)
			uao_neutral_vel = np.array(uao_neutral_vel)
			uao_neutral_kvec = np.array(uao_neutral_az)
			
			#get magnetic coordinates of the station
			UAO = fpipy.FPIStation("uao")
			ANN = fpipy.FPIStation("ann")
	
			#get plasma velocity at the station
			uao_mlat, uao_mlon = UAO.get_coords(dtime_i, aacgm=True)
			uao_mcolat = 90-uao_mlat
			uao_pos = np.atleast_2d(np.array([uao_mlat, uao_mlon])).T
			uao_plasma_vel, uao_plasma_az = pf.find_gradV(uao_pos, solution, latmin, lon_shft, lat_shft, order, dtime_i)
			
			if not isinstance(uao_plasma_vel, np.ndarray):
				uao_plasma_vel= np.array(uao_plasma_vel)
			if not isinstance(uao_plasma_az, np.ndarray):
				uao_plasma_az = np.array(uao_plasma_az)
			if not isinstance(uao_mcolat, np.ndarray):
				uao_mcolat = np.array(uao_mcolat)
			if not isinstance(uao_mlon, np.ndarray):
				uao_mlon = np.array(uao_mlon)
			
			#append uao station to plasma data locations
			plasma_mcolats = np.append(plasma_mcolats, uao_mcolat)
			plasma_mlons = np.append(plasma_mlons, uao_mlon)
			plasma_vels = np.append(plasma_vels, uao_plasma_vel)
			plasma_az = np.append(plasma_az, uao_plasma_az)
			
			#get heppner-maynard boundary
			boundary_mlats = vectors.boundary_mlats[time_i]
			boundary_mlons = vectors.boundary_mlons[time_i]
			
			uao_dr, uao_dtheta = pydatadarn.plotting.vector_plot(plasma_mcolats, plasma_mlons, plasma_az, plasma_vels, time=time_i, 
											station_names=stations, FPI_names=["uao"], 
											FPI_kvecs=uao_neutral_az, FPI_vels=uao_neutral_vel, boundary_mlats=boundary_mlats, boundary_mlons=boundary_mlons, mlt=False, cart=True,
											mcolat_min=0, mcolat_max=45, theta_min=0, theta_max=360, save=True)
			
			#append to array
			neutral_vel_interval = np.append(neutral_vel_interval, uao_neutral_vel)
			neutral_az_interval = np.append(neutral_az_interval, uao_neutral_az)
			plasma_vel_interval = np.append(plasma_vel_interval, uao_plasma_vel)
			plasma_az_interval = np.append(plasma_az_interval, uao_plasma_az)
			latmins = np.append(latmins, latmin)
			times = np.append(times, time_i)
			
			print("\n")
			
#get x array
x = np.arange(0, int(60/min_interval)*(hmax+1-hmin))
tick_labels = np.array([])

time_tick = pydatadarn.tools.time_to_dtime(times[0])
minute_interval = 60
while time_tick <= pydatadarn.tools.time_to_dtime(times[-1]):
 	hour = time_tick.hour
 	minute = time_tick.minute
 	tick_labels = np.append(tick_labels, "{:02d}:{:02d}".format(hour, minute))
 	time_tick = time_tick + datetime.timedelta(minutes=60)	
	 
#add -ve sign for westwards flowing currents
# neutral_vel_interval = neutral_vel_interval*np.sign(neutral_az_interval)
# plasma_vel_interval = plasma_vel_interval*np.sign(plasma_az_interval)

plt.figure()
plt.plot(x, neutral_vel_interval, label="neutral")
plt.plot(x, plasma_vel_interval, label="plasma")
plt.axhline(0, color="r", linestyle="--")
plt.legend()
plt.xticks(np.arange(0, len(x), 6), tick_labels)
plt.grid(linestyle = "--")
plt.show()

#get east and north look directions at uao site of plasma
plasma_EW = abs(plasma_vel_interval*np.cos(np.deg2rad(plasma_az_interval)))
plasma_NS = abs(plasma_vel_interval*np.sin(np.deg2rad(plasma_az_interval)))

plasma_N = np.empty(len(plasma_NS))
plasma_E = np.empty(len(plasma_EW))

#determine if flow is positive/negative northwards/eastwards
plasma_E = plasma_EW*np.sign(plasma_az_interval)
plasma_N = plasma_NS*np.sign(np.cos(np.deg2rad(plasma_az_interval)))

plt.figure()
plt.plot(plasma_N, label = "North Look")
plt.plot(plasma_E, label = "East Look")
plt.axhline(0, color="k", linestyle = "--")
plt.legend()
plt.xticks(np.arange(0, len(x), 6), tick_labels)
plt.grid(linestyle = "--")
plt.show()
	
	

time_i = "2013/10/02 08:30:00"
time_index = np.where(vectors.times == time_i)

#convert time_i into seconds
dtime_i = pydatadarn.tools.time_to_dtime(time_i)
time_s = (dtime_i - dtime_min).seconds

##########################
### Get superDARN info ###
##########################

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

sample = vectors.get_time_data(time_i)

#get data and model lat/lon
plasma_mlats = sample["mlats"]
plasma_mlons = sample["mlons"]
mod_mlats = vectors.mod_mlats[time_i]
mod_mlons = vectors.mod_mlons[time_i]
#append model to data
all_mlats = np.append(plasma_mlats, mod_mlats)
all_mlons = np.append(plasma_mlons, mod_mlons)
plasma_mcolats = 90-plasma_mlats
#affix lats/lons into one array for spherical fitting input
pos = np.array([all_mlats, all_mlons])

#get velocity and direction from spherical fit
all_plasma_vels, all_plasma_az = pf.find_gradV(pos, solution, latmin, lon_shft, lat_shft, order, dtime_i)
plasma_vels = all_plasma_vels[0:len(plasma_mlats)]
plasma_az = all_plasma_az[0:len(plasma_mlats)]

boundary_mlats = vectors.boundary_mlats[time_i]
boundary_mlons = vectors.boundary_mlons[time_i]

####################
### Get FPI info ###
####################

#get data for specific time_i
N_time_i = fpipy.interpolate(hN, N_time_indexes, time_s)
E_time_i = fpipy.interpolate(hE, E_time_indexes, time_s)
S_time_i = fpipy.interpolate(hS, S_time_indexes, time_s)
W_time_i = fpipy.interpolate(hW, W_time_indexes, time_s)

#get East-West and North-South flow across the FPI station
hEW = (E_time_i + W_time_i)/2
hNS = (N_time_i + S_time_i)/2

#compute full 2D vector for neutral wind
uao_neutral_vel, uao_neutral_az = fpipy.merge_vecs(hNS, hEW)
uao_neutral_vel = np.array(uao_neutral_vel)
uao_neutral_kvec = np.array(uao_neutral_az)

#compute full 2D vector for plasma wind
UAO = fpipy.FPIStation("uao")
ANN = fpipy.FPIStation("ann")

#obtain pasma flow at station
uao_mlat, uao_mlon = UAO.get_coords(dtime_i, aacgm=True)
uao_mcolat = 90-uao_mlat
uao_pos = np.atleast_2d(np.array([uao_mlat, uao_mlon])).T
uao_plasma_vel, uao_plasma_az = pf.find_gradV(uao_pos, solution, latmin, lon_shft, lat_shft, order, dtime_i)

######################
### Make the Plots ###
######################

tick_labels = np.array([])

time_tick = dtime_min
minute_interval = 60
while time_tick <= dtime_max:
 	hour = time_tick.hour
 	minute = time_tick.minute
 	tick_labels = np.append(tick_labels, "{:02d}:{:02d}".format(hour, minute))
 	time_tick = time_tick + datetime.timedelta(minutes=60)

fig, ax = plt.subplots(3, 1, sharex=True)
ax0 = ax[0]
ax1 = ax[1]
ax2 = ax[2]

ax0.errorbar(W_time_indexes, hW, yerr = dW, marker="o", markersize=2.5, label = "West Look")
ax0.errorbar(E_time_indexes, hE, yerr = dE, marker="o", markersize=2.5, label = "East Look")
ax0.set_xticks(np.arange(0, num_seconds+minute_interval*60, minute_interval*60))
ax0.set_xticklabels(tick_labels)
ax0.axhline(y=0, color = "k", linestyle = "-.")
ax0.grid(linestyle = "--")
ax0.legend()

ax1.errorbar(S_time_indexes, hS, yerr = dS, marker="o", markersize=2.5, label = "South Look")
ax1.errorbar(N_time_indexes, hN, yerr = dN, marker="o", markersize=2.5, label = "North Look")
ax1.axhline(y=0, color = "k", linestyle = "-.")
ax1.grid(linestyle = "--")
ax1.legend()

ax2.errorbar(zen_time_indexes, zen, yerr=dzen, marker="o", markersize=2.5, label = "Zenith")
ax2.axhline(y=0, color = "k", linestyle = "-.")
ax2.set_yticklabels([-600, -450, -300, -150, 0, 150])
ax2.grid(linestyle = "--")
ax2.legend()

plt.show()

#plot full convection map
#convert plasma vectors from magnetic to geographic
#plasma_mlats, plasma_mlons, alt = aacgmv2.convert_latlon_arr(mlats, mlons, 0, dtime_i, "A2G")
plasma_mcolats = 90-plasma_mlats

#convert neutral wind vector from geographic to aacgm
#find how much look direction is different
#get coordinates of geographic north in aacgm
glat, glon, galt = aacgmv2.convert_latlon(90, 0, 0, dtime_i)
#get coordinates of UAO station in mlat
uao_lat, uao_lon = UAO.get_coords(dtime_i, aacgm=True)
#calculate angle difference
look_change = pydatadarn.tools.cosine_rule([glat, glon], [uao_lat, uao_lon], [0, 0], polar=True)
#figure out if geographic north is east or west of magnetic to determine which way the look changes
if pydatadarn.tools.lon_look(0, glon) == "E":
	look_change = -look_change
#add our look difference
uao_neutral_az += look_change

#plot only negative plasma kvectors
# test_indexes = np.where(plasma_az < 0)
# plasma_mcolats = plasma_mcolats[test_indexes]
# plasma_mlons = plasma_mlons[test_indexes]
# plasma_az = plasma_az[test_indexes]
# plasma_vels = plasma_vels[test_indexes]

dr, dtheta = pydatadarn.plotting.vector_plot(plasma_mcolats, plasma_mlons, plasma_az, plasma_vels, time=time_i, 
								station_names=stations, FPI_names=["uao"], 
								FPI_kvecs=uao_neutral_az, FPI_vels=uao_neutral_vel, 
								boundary_mlats=boundary_mlats, boundary_mlons=boundary_mlons, 
								mlt=False, cart=True, 
								mcolat_min=0, mcolat_max=45, theta_min=0, theta_max=360)

#plot only neutral wind and ion vector at FPI station
#get velocity and direction from spherical fit
uao_pos = np.atleast_2d(np.array([uao_mlat, uao_mlon])).T
uao_plasma_vel, uao_plasma_az = pf.find_gradV(uao_pos, solution, latmin, lon_shft, lat_shft, order, dtime_i)

if not isinstance(uao_plasma_vel, np.ndarray):
	uao_plasma_vel= np.array(uao_plasma_vel)
if not isinstance(uao_plasma_az, np.ndarray):
	uao_plasma_az = np.array(uao_plasma_az)
if not isinstance(uao_mcolat, np.ndarray):
	uao_mcolat = np.array(uao_mcolat)
if not isinstance(uao_mlon, np.ndarray):
	uao_mlon = np.array(uao_mlon)

uao_dr, uao_dtheta = pydatadarn.plotting.vector_plot(uao_mcolat, uao_mlon, uao_plasma_az, uao_plasma_vel, time=time_i, 
								station_names=stations, FPI_names=["uao"], 
								FPI_kvecs=uao_neutral_az, FPI_vels=uao_neutral_vel, mlt=False, cart=True,
								mcolat_min=0, mcolat_max=45, theta_min=0, theta_max=360, save=True)

area_lats = np.arange(40, 90, 1)
area_lons = np.arange(-20, 20, 1)
areas = np.meshgrid(area_lons, area_lats)
test_lons = areas[0].flatten()
test_lats = areas[1].flatten()
test_mcolats = 90-test_lats

test_pos = np.atleast_2d(np.array([test_lats, test_lons]))

test_vels, test_az = pf.find_gradV(test_pos, solution, latmin, lon_shft, lat_shft, order)

dr, dtheta = pydatadarn.plotting.vector_plot(test_mcolats, test_lons, test_az, test_vels, time=time_i, 
								station_names=stations, FPI_names=["uao"], 
								FPI_kvecs=uao_neutral_az, FPI_vels=uao_neutral_vel, mlt=False, 
								mcolat_min=0, mcolat_max=45, theta_min=0, theta_max=90, cbar_min=min(plasma_vels), 
								cbar_max=max(plasma_vels))