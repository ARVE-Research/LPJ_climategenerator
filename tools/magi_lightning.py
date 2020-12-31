# create the environment with all packages in conda with
# conda env create -f lightning.yml
import netCDF4 as ncdf
import numpy as np
import math
from calendar import monthrange
import time
import sys

## INPUT ################################################
mc_input = str(sys.argv[1]) #'/home/akoch/group_data/CMIP6/scenarioMIP/ssp370/Amon/mc/mc_Amon_BCC-CSM2-MR_ssp370_r1i1p1f1_gn_201501-210012.nc'
land_input = str(sys.argv[2]) #'/home/akoch/group_data/CMIP6/land_mask/sftlf_fx_BCC-CSM2-MR_hist-resIPO_r1i1p1f1_gn.nc'
lgtng_file = str(sys.argv[3]) #'/home/akoch/data/future_lightning/lght_Amon_BCC-CSM2-MR_ssp370_r1i1p1f1_gn_201501-210012.nc' 
# set upper boundary
if len(sys.argv)==5:
   lgtng_present_file = str(sys.argv[4])
   if len(lgtng_present_file) < 4: # catch numbers instead of filenames
      lgtng_present_file = float(lgtng_present_file)
else:
   lgtng_present_file = None
     

#########################################################
# Land mask
land_nc = ncdf.Dataset(land_input, 'r')
land_name = [var for var in land_nc.variables 
             if 'land' in var or 'landmask' in var or 'sftlf' in var]
land_mask = land_nc.variables[land_name[0]][:]
if land_mask.max() > 1:
    land_mask  = np.ma.where(land_mask > 0, 1, 0) 

sea_mask = land_mask.max() - land_mask
if land_mask.ndim==2:
    masks = [land_mask[:], sea_mask[:]]
else:
    masks = [land_mask[0,:,:], sea_mask[0,:,:]]

# extract Mass Convective Flux
mc_nc = ncdf.Dataset(mc_input, 'r')
mc_name = [var for var in mc_nc.variables if 'mc' in var or 'MCamFX' in var]
mc_var = mc_nc.variables[mc_name[0]][:]
lon_var = mc_nc.variables['lon'][:]
lat_var = mc_nc.variables['lat'][:]
time_var = mc_nc.variables['time'][:]

# fifth polynomial fit from Magi (2015)
# L = a1*M + a2*M^{2} + a3*M^{3} + a4*M^{4} + a5*M^{5}; M = kg m^{-2} h^{-1}
mc_var = mc_var * 3600

a1 = [1.31e-1, 1.06e-2]
a2 = [-1.93e-3, 6.71e-4]
a3 = [8.87e-5, -8.58e-5]
a4 = [-6.84e-6, 2.41e-6]
a5 = [1.17e-7, -1.97e-8]

# calculate lightning for mask
lgtng_var = np.zeros(mc_var.shape)
for i in range(2):
    mc_masked = mc_var * masks[i]
    lgtng_var = (lgtng_var 
                + a1[i] 
                * mc_masked 
                + a2[i] 
                * mc_masked**2 
                + a3[i] 
                * mc_masked**3 
                + a4[i] 
                * mc_masked**4 
                + a5[i] 
                * mc_masked**5)


# convert from flashes km-2 month-1 to flashes km-2 day-1
year_unit = mc_nc.variables.get('time').calendar

days_month = list()
if year_unit=='365_day' or year_unit=='noleap':
    for i in range(1,13): days_month.append(monthrange(2002, i)[1])

days_month = days_month * int((lgtng_var.shape[0] / 12))

for i in range(0,len(days_month)):
    lgtng_var[i,:,:] = lgtng_var[i,:,:] / days_month[i]

# set min-max to present-day
if lgtng_present_file is None:
	print('Unbound lightning used. Check output for extremes.')
elif type(lgtng_present_file) is float:
    L_max = lgtng_present_file
elif lgtng_present_file=='observations':
    L_max = 10 # max value from observations
else:
    lgtng_present_nc = ncdf.Dataset(lgtng_present_file, 'r')
    L_present_var = lgtng_present_nc.variables['L'][:]
    L_max = np.max(L_present_var)

if lgtng_present_file is not None:
    lgtng_var[lgtng_var > L_max] = L_max
    
lgtng_var[lgtng_var < 0] = 0

# create netcdf with lightning data
lgtng_nc = ncdf.Dataset(lgtng_file, 'w')

# dimensions
lgtng_nc.createDimension('time', size = None)
lgtng_nc.createDimension('lon', size = len(lon_var))
lgtng_nc.createDimension('lat', size = len(lat_var))

# variables
month = lgtng_nc.createVariable(varname = 'time', 
                                datatype = 'i', 
                                dimensions = ('time',))
lon = lgtng_nc.createVariable(varname = 'lon', 
                              datatype = 'f4', 
                              dimensions = ('lon',))
lat = lgtng_nc.createVariable(varname = 'lat', 
                              datatype = 'f4', 
                              dimensions = ('lat',))
lgtng = lgtng_nc.createVariable(varname = 'L', 
                                datatype = 'f8', 
                                dimensions = ('time', 'lat', 'lon'), 
                                fill_value = 1e+20)

lon[:] = lon_var
lat[:] = lat_var
month[:] = time_var
lgtng[:] = lgtng_var

# attributes
lgtng_nc.description = 'Lightning density calculated from \
                        convective mass flux following Magi (2015)'
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
land_nc.close()
print(lgtng_file)
