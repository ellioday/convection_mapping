#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 17:21:52 2021

@author: elliott
"""

import numpy as np
import pydatadarn
import datetime
import elliotools

import matplotlib.pyplot as plt

import fpipy
import aacgmv2
import json

stations=["ade", "adw", "bks", "cly", "cve", "cvw", "fhe", "fhw", "gbr", "han",
		   "hkw", "hok", "inv", "kap", "kod", "ksr", "lyr", "pgr", "pyk", 
		   "rkn", "sas", "sto", "wal"]
FPI_stations=["ann", "uao"]

##########################
### Get SuperDARN data ###
##########################

vectors = pydatadarn.GridData(is_map=True)
vectors.add_data("all", "2013/10/02 00:00:00", "2013/10/02 11:00:00")

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
time_max = "2013/10/02 11:00:00"
dtime_min = elliotools.time_to_dtime(time_min)
dtime_max = elliotools.time_to_dtime(time_max)

#get data 
N, E, S, W, zen, dN, dE, dS, dW, dzen, N_times, E_times, S_times, W_times, zen_times = FPI.get_azm_vels(
	dtime_min = dtime_min, dtime_max = dtime_max)

#get dtimes
N_dtimes = elliotools.time_to_dtime(N_times)
E_dtimes = elliotools.time_to_dtime(E_times)
S_dtimes = elliotools.time_to_dtime(S_times)
W_dtimes = elliotools.time_to_dtime(W_times)
zen_dtimes = elliotools.time_to_dtime(zen_times)

#calculate number of seconds between start and end time
dtime_min = elliotools.time_to_dtime(time_min)
dtime_max = elliotools.time_to_dtime(time_max)
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

hmin=1
hmax=11
min_interval=10

#get x array
x = np.arange(0, int(60/min_interval)*(hmax+1-hmin))
class Lat_Lon_Storage():
	def __init__(self, lats, lons):
		self.lats=lats
		self.lons=lons
lat_lon_dict = dict()
HMB_storage = np.empty([len(x), 2, 73])
HMB_storage[:] = np.nan

count = 0

loop=True
if loop == True:
	for hour in range(hmin, hmax+1):
		for minute in range(0, 60, min_interval):
			time_i = "2013/10/02 {:02d}:{:02d}:00".format(hour, minute)
			print(time_i)
			time_index = np.where(vectors.times == time_i)
			
			#convert time_i into seconds
			dtime_i = elliotools.time_to_dtime(time_i)
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
			plasma_vels, plasma_az = pydatadarn.find_gradV(pos, solution, latmin, lon_shft, lat_shft, order, dtime_i)
			
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
			uao_plasma_vel, uao_plasma_az = pydatadarn.find_gradV(uao_pos, solution, latmin, lon_shft, lat_shft, order, dtime_i)
			
			lat_lon_dict[time_i] = Lat_Lon_Storage(plasma_mlats, plasma_mlons)
			
			if not isinstance(uao_plasma_vel, np.ndarray):
				uao_plasma_vel= np.array(uao_plasma_vel)
			if not isinstance(uao_plasma_az, np.ndarray):
				uao_plasma_az = np.array(uao_plasma_az)
			if not isinstance(uao_mcolat, np.ndarray):
				uao_mcolat = np.array(uao_mcolat)
			if not isinstance(uao_mlon, np.ndarray):
				uao_mlon = np.array(uao_mlon)
			
			#get heppner-maynard boundary
			boundary_mlats = vectors.boundary_mlats[time_i]
			boundary_mlons = vectors.boundary_mlons[time_i]
			
			los_vs_index = np.where(vectors.times == time_i)
			los_vs = vectors.los_vs[los_vs_index]
			mcolats = vectors.mcolats[los_vs_index]
			mlons = vectors.mlons[los_vs_index]
			kvecs = vectors.kvecs[los_vs_index]
			
			#create dictionary from all relevant information
			data_dict = dict()
			data_dict["time"] = time_i
			data_dict["uao_hEW"] = hEW.tolist()
			data_dict["uao_hNS"] = hNS.tolist()
			data_dict["uao_neutral_vel"] = uao_neutral_vel.tolist()
			data_dict["uao_neutral_kvec"] = uao_neutral_kvec.tolist()
			data_dict["uao_mcolat"] = uao_mcolat.tolist()
			data_dict["uao_mlon"] = uao_mlon.tolist()
			data_dict["uao_plasma_vel"] = uao_plasma_vel.tolist()
			data_dict["uao_plasma_az"] = uao_plasma_az.tolist()
			data_dict["plasma_mcolats"] = plasma_mcolats.tolist()
			data_dict["plasma_mlons"] = plasma_mlons.tolist()
			data_dict["plasma_vels"] = plasma_vels.tolist()
			data_dict["plasma_az"] = plasma_az.tolist()
			data_dict["hmb_mlats"] = boundary_mlats.tolist()
			data_dict["hmb_mlons"] = boundary_mlons.tolist()
			data_dict["N"] = N.tolist()
			data_dict["N1"] = N1.tolist()
			data_dict["N2"] = N2.tolist()
			data_dict["N3"] = N3.tolist()
			data_dict["latmin"] = latmin
			data_dict["lat_shft"] = lat_shft
			data_dict["lon_shft"] = lon_shft
			data_dict["order"] = order
			data_dict["los_vs"] = los_vs.tolist()
			data_dict["los_mcolats"] = mcolats.tolist()
			data_dict["los_mlons"] = mlons.tolist()
			data_dict["los_kvecs"] = kvecs.tolist()
			
			with open("data/{}{}{}_{}{}{}.json".format(time_i[0:4], time_i[5:7], 
										 time_i[8:10], time_i[11:13], time_i[14:16], 
										 time_i[17:19]), "w") as outfile:
				json.dump(data_dict, outfile)
			
			count += 1