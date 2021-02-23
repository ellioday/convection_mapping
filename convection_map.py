#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 13:56:40 2021

@author: elliott
"""

import numpy as np
import pydatadarn

stations=["ade", "adw", "bks", "cly", "cve", "cvw", "fhe", "fhw", "gbr", "han",
		   "hkw", "hok", "inv", "kap", "kod", "ksr", "lyr", "pgr", "pyk", 
		   "rkn", "sas", "sto", "wal"]

#in map_plot.c line 1080 (1191) -> 1604


vectors = pydatadarn.GridData(is_map=True)
vectors.add_data("all", "2013/10/02 06:00:00", "2013/10/02 07:00:00")
sample = vectors.sample

hmb_mlat = sample.boundary.mlat
hmb_mcolat = 90-hmb_mlat
hmb_mlon = sample.boundary.mlon








time_i = "2013/10/02 06:41:00"
time_index = np.where(vectors.times == time_i)
mcolats = vectors.mcolats[time_index]
mlons = vectors.mlons[time_index]
kvecs = vectors.kvecs[time_index]
los_vs = vectors.los_vs[time_index]

#plot
dr, dtheta = pydatadarn.plotting.vector_plot(mcolats, mlons, kvecs, los_vs, hmb_mcolat, hmb_mlon, time=time_i, 
								station_names = stations, mlt=True, 
								theta_min=0, theta_max=360, 
								cbar_min=min(los_vs), cbar_max=max(los_vs))

mod_time_index = np.where(vectors.mod_times == time_i)
mod_mcolats = vectors.mod_mcolats[mod_time_index]
mod_mlons = vectors.mod_mlons[mod_time_index]
mod_kvecs = vectors.mod_kvecs[mod_time_index]
mod_los_vs = vectors.mod_los_vs[mod_time_index]

#plot
mod_dr, mod_dtheta = pydatadarn.plotting.vector_plot(mod_mcolats, mod_mlons, mod_kvecs, mod_los_vs, hmb_mcolat, hmb_mlon,
								time=time_i, 
								station_names = stations, mlt=True, 
								theta_min=0, theta_max=360, 
								cbar_min=min(mod_los_vs), cbar_max=max(mod_los_vs))

