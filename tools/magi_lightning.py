# create the environment with all packages in conda with
# conda env create -f lightning.yml
import sys
import netCDF4 as ncdf
import numpy as np
import math
from calendar import monthrange
import time

# TODO:
# make work with standard CMIP mc
n_args = len(sys.argv)

pl_input = sys.argv[1] #'Pressure_aij_E213f10aF40oQ40_climatology_1961-1990.nc'
mc_input = sys.argv[2] #'MCamFX_aijl_E213f10aF40oQ40_climatology_1961-1990.nc'
land_input = sys.argv[3] #'GISS_E2_144_90_landmask.nc'
if n_args==6:
    lgtng_present_file = sys.argv[4] #'lightning_GISS-E2-R_RCP26_all.nc'
else:
    lgtng_present_file = None
lgtng_file = sys.argv[n_args-1] #'Magi_lightning_aij_E213f10aF40oQ40_climatology_1961-1990.nc'
#########################################################
# Find pressure level
pl_nc = ncdf.Dataset(pl_input, 'r')
pl_var = pl_nc.variables['Pressure'][:]
lat_var = pl_nc.variables['lat'][:]
lon_var = pl_nc.variables['lon'][:]
lev_var = pl_nc.variables['levels'][:]
time_var = pl_nc.variables['time'][:]

# Calculate weights
weight_arr = np.ones([time_var.size, lat_var.size, lon_var.size])
weight = np.cos(lat_var * math.pi / 180)
weight = abs(weight)
weight = weight.reshape(weight.size, 1)
weight = weight_arr * weight

# Calculate weighted mean for each layer
pl_mean = np.empty(lev_var.size)
for l in range(lev_var.size):
    pl_mean[l] = np.ma.average(pl_var[:,l,:], weights=weight)

# Closest mean pressure to 427 hPa
idx = (np.abs(pl_mean - 427)).argmin()

pl_nc.close()

# Land mask
land_nc = ncdf.Dataset(land_input, 'r')
land_mask = land_nc.variables['land'][:]
sea_mask = 1 - land_mask
masks = [land_mask, sea_mask]
land_nc.close()

# extract weighted-mean Mass Convective Flux at p level
mc_nc = ncdf.Dataset(mc_input, 'r')
mc_name = [var for var in mc_nc.variables if 'mc' in var or 'MCamFX' in var]
mc_var = mc_nc.variables[mc_name[0]][:]
mc_var = mc_var[:,idx,:,:]

# fifth polynomial fit from Magi (2015)
# L = a1*M + a2*M^{2} + a3*M^{3} + a4*M^{4} + a5*M^{5}; M = kg m^{-2} h^{-1}
# convert from 1e-4 kg m2 s1 if neeed
mc_var = mc_var * 1 / 3.6
mc_var.min()
a1 = [1.31e-1, 1.06e-2]
a2 = [-1.93e-3, 6.71e-4]
a3 = [8.87e-5, -8.58e-5]
a4 = [-6.84e-6, 2.41e-6]
a5 = [1.17e-7, -1.97e-8]

# calculate lightning for mask
lgtng_var = np.zeros(mc_var.shape)
for i in range(1):
    mc_masked = mc_var * masks[i]
    lgtng_var = lgtng_var + a1[i] * mc_masked + a2[i] * mc_masked**2 + a3[i] * mc_masked**3 + a4[i] * mc_masked**4 + a5[i] * mc_masked**5

# convert from flashes km-2 month-1 to flashes km-2 day-1
year_unit = mc_nc.variables.get('time').calendar

days_month = list()
if year_unit=='365_day' or year_unit=='noleap':
    for i in range(1,13): days_month.append(monthrange(2002, i)[1])

for i in range(12):
    lgtng_var[i,:,:] = lgtng_var[i,:,:] / days_month[i]

# set min-max to present-day
if lgtng_present_file is None:
    L_max = 0.06375045 # max value from lightning_GISS-E2-R_RCP26 climatology
else:
    lgtng_present_nc = ncdf.Dataset(lgtng_present_file, 'r')
    L_present_var = lgtng_present_nc.variables['L'][:]
    L_max = np.max(L_present_var)

lgtng_var[lgtng_var > L_max] = L_max
lgtng_var[lgtng_var < 0] = 0

# create netcdf with lightning data
lgtng_nc = ncdf.Dataset(lgtng_file, 'w')

# dimensions
lgtng_nc.createDimension('time', size = None)
lgtng_nc.createDimension('lon', size = 144)
lgtng_nc.createDimension('lat', size = 90)

# variables
month = lgtng_nc.createVariable(varname = 'time', datatype = 'i', dimensions = ('time',))
lon = lgtng_nc.createVariable(varname = 'lon', datatype = 'f4', dimensions = ('lon',))
lat = lgtng_nc.createVariable(varname = 'lat', datatype = 'f4', dimensions = ('lat',))
lgtng = lgtng_nc.createVariable(varname = 'L', datatype = 'f8', dimensions = ('time', 'lat', 'lon'), fill_value = 1e+20)

lon[:] = lon_var
lat[:] = lat_var
month[:] = mc_nc.variables['time'][:]
lgtng[:] = lgtng_var

# attributes
lgtng_nc.description = 'Lightning density calculated from convective mass flux following Magi (2015)'
lgtng_nc.history = 'Created ' + time.ctime(time.time())
lgtng_nc.source = 'Alexander Koch [akoch@hku.hk]'
lgtng.standard_name = 'lighting density'
lgtng.long_name = 'Parameterized lightning density'
lgtng.units = 'flashes km^-2 day^-1'
lgtng_nc['time'].setncatts(mc_nc['time'].__dict__)
month.standard_name = 'time'
lgtng_nc['lon'].setncatts(mc_nc['lon'].__dict__)
lgtng_nc['lat'].setncatts(mc_nc['lat'].__dict__)

# done
lgtng_nc.close()
mc_nc.close()
