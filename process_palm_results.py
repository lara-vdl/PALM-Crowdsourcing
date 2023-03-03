'''
Script to process the PALM results
- first: potential temperature is extracted and converted to air temperature, then buildings are clipped out of the data
- second: urban rural temperatures are calculated for both domains
- third: palm results at the location of the LMSS stations are extracted for a comparison 
- fourth: palm energy balance results at the location of the LMSS stations are extracted 
- fifth: palm results at the location of the Netatmo stations are extracted for statistical evaluation

Results from step 1 and 2 are saved as netcdf, coordinates are defined in those files to allow for display in QGIS at the correct georeferende location
Results from step 3 and 4 are saved as csv files for analysis in R
'''

## load required libraries
import sys
import xarray as xr
import numpy as np
import geopandas
import rioxarray as rio
from shapely.geometry import mapping
from pyproj import Proj
from xarray.core.common import T
import pandas as pd

## load data
# palm results
palm_parent = xr.open_dataset("/media/lara/2TB SSD/Data_MA/Model_output/city/isotropic/v9/MA_v9_av_xy.nc")
palm_child = xr.open_dataset("/media/lara/2TB SSD/Data_MA/Model_output/city/isotropic/v9/MA_v9_av_xy_N02.nc")
child_3d = xr.open_dataset("/media/lara/2TB SSD/Data_MA/Model_output/city/isotropic/v9/MA_v9_av_3d_N02.nc")

# shape files
buildings_p = geopandas.read_file("/media/lara/2TB SSD/Data_MA/Model_input/10m_city/buildings_10m_city.shp", crs = "epsg:25832")
buildings_c = geopandas.read_file("/media/lara/2TB SSD/Data_MA/Model_input/city/buildings_city_new2.shp", crs = "epsg:25832")
lcz = geopandas.read_file("~/Dokumente/Data_processing/data_in/LCZ_14.shp",  crs = "epsg:25832")

# netatmo data
netatmo = pd.read_csv("~/Dokumente/Data_processing/data_in/netatmo_utm.csv")

##### part 1: convert potential temperature to air temperature for the parent and child domain ###############################################
# extract theta at 2m
theta_2m_par = palm_parent.variables["theta_2m*_xy"]
theta_2m_child = palm_child.variables["theta_2m*_xy"]

## convert to absolute temperature
# necessary variables for conversion
# gas constant of dry air
r = 287.05
# specific heat capacity of air
cp = 1004.67
# reference pressure
p0 = 1000
# current pressure
p = 1006

# calculate ta
divider = (p0/p)**(r/cp)

ta_2m_p = theta_2m_par/divider
ta_2m_c = theta_2m_child/divider

# convert to °C
ta_2m_p = np.round(ta_2m_p - 273.15, 1)
ta_2m_p = np.squeeze(ta_2m_p)

ta_2m_c = np.round(ta_2m_c - 273.15, 1)
ta_2m_c = np.squeeze(ta_2m_c)

## extract latitude and longitude to define crs for display in QGIS
# utm coordinates
utm_y_p = palm_parent.variables["N_UTM"]
utm_x_p = palm_parent.variables["E_UTM"]

utm_y_c = palm_child.variables["N_UTM"]
utm_x_c = palm_child.variables["E_UTM"]

# define utm projection as proj object
utm_proj = Proj("+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# coordinates as arrays
utm_x_arr_p = np.array(utm_x_p)
utm_y_arr_p = np.array(utm_y_p)

utm_x_arr_c = np.array(utm_x_c)
utm_y_arr_c = np.array(utm_y_c)

# meshgrid of coordinates
x_mesh_p, y_mesh_p = np.meshgrid(utm_x_arr_p, utm_y_arr_p)
x_mesh_c, y_mesh_c = np.meshgrid(utm_x_arr_c, utm_y_arr_c)

# convert utm coordinates to lat/lon
lon_p, lat_p = utm_proj(x_mesh_p, y_mesh_p, inverse=True)
lon_c, lat_c = utm_proj(x_mesh_c, y_mesh_c, inverse=True)

## extract time
time = palm_parent.coords["time"]

## convert to xarray DataArray, x and y coordinates must contain UTM coordinates
ta_2m_parent =xr.DataArray(data=ta_2m_p, dims=["time","y","x"], \
    coords=dict(x=(["x"], utm_x_p), y=(["y"], utm_y_p), time = time, lat = (["y","x"], lat_p), lon = (["y","x"], lon_p)))

ta_2m_child =xr.DataArray(data=ta_2m_c, dims=["time","y","x"], \
    coords=dict(x=(["x"], utm_x_c), y=(["y"], utm_y_c), time = time, lat = (["y","x"], lat_c), lon = (["y","x"], lon_c)))

## save as tiff
ta_2m_parent = ta_2m_parent.rio.set_spatial_dims(x_dim="x",y_dim="y")
ta_2m_parent.rio.write_crs("epsg:25832", inplace=True)
ta_2m_parent.rio.to_raster("data_paper/ta_2m_parent.tiff")

ta_2m_child = ta_2m_child.rio.set_spatial_dims(x_dim="x",y_dim="y")
ta_2m_child.rio.write_crs("epsg:25832", inplace=True)
ta_2m_child.rio.to_raster("data_paper/ta_2m_child.tiff")

## clip out values at buildings locations
# use an invert = True in clip to set values to NA at building locations
ta_2m_parent_clip = ta_2m_parent.rio.clip(geometries = buildings_p.geometry.apply(mapping), crs=25832, invert = True)
ta_2m_child_clip = ta_2m_child.rio.clip(geometries = buildings_c.geometry.apply(mapping), crs=25832, invert = True)

## save as tiff
ta_2m_parent_clip = ta_2m_parent_clip.rio.set_spatial_dims(x_dim="x",y_dim="y")
ta_2m_parent_clip.rio.write_crs("epsg:25832", inplace=True)
ta_2m_parent_clip.rio.to_raster("data_paper/ta_2m_parent_clip.tiff")

ta_2m_child_clip = ta_2m_child_clip.rio.set_spatial_dims(x_dim="x",y_dim="y")
ta_2m_child_clip.rio.write_crs("epsg:25832", inplace=True)
ta_2m_child_clip.rio.to_raster("data_paper/ta_2m_child_clip.tiff")

#### part 2: calculate the urban rural temperature differences ##############################################################################
## clip out results pixel within LCZ 14, shape with LCZ preprocessed in QGIS
ta_rural_parent = ta_2m_parent_clip.rio.clip(geometries = lcz.geometry.apply(mapping), crs = 25832)

## calculate mean "rural" ta for city and la
ta_rural_mean = ta_rural_parent.reduce(np.nanmean, dim = ("x","y"))

## calculate uhi
uhi_p = np.round(ta_2m_p - ta_rural_mean, 1)
uhi_c = np.round(ta_2m_c - ta_rural_mean, 1)

## combine as DataArray
uhi_parent =xr.DataArray(data=uhi_p, dims=["time","y","x"], \
    coords=dict(x=(["x"], utm_x_p), y=(["y"], utm_y_p), time = time, lat = (["y","x"], lat_p), lon = (["y","x"], lon_p)))

uhi_child =xr.DataArray(data=uhi_c, dims=["time","y","x"], \
    coords=dict(x=(["x"], utm_x_c), y=(["y"], utm_y_c), time = time, lat = (["y","x"], lat_c), lon = (["y","x"], lon_c)))

## save as tiff
uhi_parent = uhi_parent.rio.set_spatial_dims(x_dim="x",y_dim="y")
uhi_parent.rio.write_crs("epsg:25832", inplace=True)
uhi_parent.rio.to_raster("data_paper/uhi_parent.tiff")

uhi_child = uhi_child.rio.set_spatial_dims(x_dim="x",y_dim="y")
uhi_child.rio.write_crs("epsg:25832", inplace=True)
uhi_child.rio.to_raster("data_paper/uhi_child.tiff")

##### part 3: combine the palm results with the lmss data ######################################################################################
## get the index of the required coordinates
# garden
i = np.argmin([np.abs(utm_x_c - 376130.12)])
j = np.argmin([np.abs(utm_y_c - 5706005.70)])

# roof
k = np.argmin([np.abs(utm_x_c - 375959.80)])
l = np.argmin([np.abs(utm_y_c - 5705596.18)])

# extract wind and rh variables from dataset
wind = palm_child.variables["wspeed_10m*_xy"]
rh = child_3d.variables["rh"]

# extract variable for desired location via indexing
ta_lmss = ta_2m_c[:,j-1:j+2,i-1:i+2]
wind_lmss = wind[:,:,l-1:l+2,k-1:k+2]
rh_lmss = rh[:,:,j-1:j+2,i-1:i+2]

# calculate mean
ta_lmss_mean = ta_lmss.reduce(np.nanmean, dim=("x","y"))
wind_lmss_mean = wind_lmss.reduce(np.nanmean, dim=("x","y"))
rh_lmss_mean = rh_lmss.reduce(np.nanmean, dim=("x","y"))

# change into array and round values
ta_lmss_mean = np.array(ta_lmss_mean).round(1)
wind_lmss_mean = np.array(wind_lmss_mean).round(1)
rh_lmss_mean = np.array(rh_lmss_mean).round(1)

# remove unneccessary dimensions
ta_lmss_squz = np.squeeze(ta_lmss_mean)
wind_lmss_sqz = np.squeeze(wind_lmss_mean)

# get index of first non nan value in rh array and slice array
rh_index = np.where(np.isnan(rh_lmss_mean[0,:]))[0][-1] +1
rh_lmss_mean = rh_lmss_mean[:,rh_index]

## get time array
# add timedelta to origin time to get each time step as datetime
date = np.datetime64('2020-08-08T00:00:00') + np.array(time)

## combine variable and datetime
# join in pd dataframe and print to check
ta_rh_time = pd.DataFrame({"date":date,"palm_ta":ta_lmss_squz,"palm_rh":rh_lmss_mean})
wind_time = pd.DataFrame({"date":date, "palm_ws":wind_lmss_sqz})

# save as csv
ta_rh_time.to_csv("~/Dokumente/Data_processing/data_paper/palm_lmss.csv")
wind_time.to_csv("~/Dokumente/Data_processing/data_paper/palm_lmss_ws.csv")

##### part 4: extract palm energy balance at lmss location #########################################################################################
# get energy balance variables
sens_hf = palm_child.variables["shf*_xy"]
lat_hf = palm_child.variables["qsws*_xy"]
ground_hf = palm_child.variables["ghf*_xy"]

# extract energy balance at lmss location
sens_hf_lmss = sens_hf[:,:,j-1:j+2,i-1:i+2]
lat_hf_lmss = lat_hf[:,:,j-1:j+2,i-1:i+2]
ground_hf_lmss = ground_hf[:,:,j-1:j+2,i-1:i+2]

# calculate means
sens_hf_mean = sens_hf_lmss.reduce(np.nanmean, dim=("x","y"))
lat_hf_mean = lat_hf_lmss.reduce(np.nanmean, dim=("x","y"))
ground_hf_mean = ground_hf_lmss.reduce(np.nanmean, dim=("x","y"))

# round and change into np array
sens_hf_mean = np.array(sens_hf_mean).round(1)
lat_hf_mean = np.array(lat_hf_mean).round(1)
ground_hf_mean = np.array(ground_hf_mean).round(1)

# remove empty dimensions
sens_hf_squz = np.squeeze(sens_hf_mean)
lat_hf_squz = np.squeeze(lat_hf_mean)
ground_hf_squz = np.squeeze(ground_hf_mean)

# combine variables in dataframe
energy = pd.DataFrame({"date":date,"sens_hf":sens_hf_squz,"lat_hf":lat_hf_squz,"ground_hf":ground_hf_squz})

# save
energy.to_csv("~/Dokumente/Data_processing/data_out/palm_lmss_energy.csv")

##### part 5: combine palm results with netatmo data ################################################################################################
## reduce netatmo dataframe to individual stations and only utm coordinates
stations = netatmo.drop_duplicates("intern_id")
stations = stations[["intern_id","utm_x","utm_y"]]

# filter for study areas by coordinates
# parent
stations = stations.loc[(stations["utm_x"] >= 372200.0) & (stations["utm_x"] <=381000.0)]
stations = stations.loc[(stations["utm_y"] >= 5702550.0) & (stations["utm_y"] <=5711050.0)]

# child
stations2 = stations.loc[(stations["utm_x"] >= 375480.0) & (stations["utm_x"] <=377280.0)]
stations2 = stations2.loc[(stations2["utm_y"] >= 5705110.0) & (stations2["utm_y"] <=5706610.0)]

## extract palm cells for each station
# change coordinates to numpy array
x = np.array(utm_x_p)
y = np.array(utm_y_p)

x2 = np.array(utm_x_c)
y2 = np.array(utm_y_c)

# parent
# initialise empty numpy arrays, specify data type for time array
palm_data = np.zeros(0)
time_arr = np.zeros(0, dtype=np.datetime64)
ids = np.zeros(0)

# use for loop to get data at netatmo locations and average data from radius
for i in range(0,len(stations.index)):
    # get index off cell with station
    index_x = np.abs(x - stations.iloc[i,1]).argmin()
    index_y = np.abs(y - stations.iloc[i,2]).argmin()

    # slice out cell with station and surrounding stations, calculate mean
    # parent: 9 cells
    cells = ta_2m_parent_clip[:,index_y-1:index_y+2,index_x-1:index_x+2]
    cell_mean = cells.reduce(np.nanmean, dim=("x","y"))

    # convert to numpy array and round values
    cell_mean = np.array(cell_mean).round(1)

    # create numpy array with station id
    station_id = np.full(len(cell_mean), stations.iloc[i,0])

    # append to numpy arrays
    palm_data = np.r_[palm_data, cell_mean]
    ids = np.r_[ids, station_id]
    time_arr = np.r_[time_arr, date]

# child
# initialise empty numpy arrays, specify data type for time array
palm_data2 = np.zeros(0)
time_arr2 = np.zeros(0, dtype=np.datetime64)
ids2 = np.zeros(0)

# use for loop to get data at netatmo locations and average data from radius
for i in range(0,len(stations2.index)):
    # get index off cell with station
    index_x = np.abs(x2 - stations2.iloc[i,1]).argmin()
    index_y = np.abs(y2 - stations2.iloc[i,2]).argmin()

    # slice out cell with station and surrounding stations, calculate mean
    # child: 25 cells
    cells = ta_2m_child_clip[:,index_y-2:index_y+3,index_x-2:index_x+3]
    cell_mean = cells.reduce(np.nanmean, dim=("x","y"))

    # convert to numpy array and round values
    cell_mean = np.array(cell_mean).round(1)

    # create numpy array with station id
    station_id = np.full(len(cell_mean), stations2.iloc[i,0])

    # append to numpy arrays
    palm_data2 = np.r_[palm_data2, cell_mean]
    ids2 = np.r_[ids2, station_id]
    time_arr2 = np.r_[time_arr2, date]

# merge all arrays in dataframe
palm_df = pd.DataFrame({"station_id":ids, "time":time_arr, "palm_ta":palm_data}, columns=["station_id","time","palm_ta"])
palm_df2 = pd.DataFrame({"station_id":ids2, "time":time_arr2, "palm_ta":palm_data2}, columns=["station_id","time","palm_ta"])

# save as csv
palm_df.to_csv("~/Dokumente/Data_processing/data_paper/palm_netatmo_parent.csv")
palm_df2.to_csv("~/Dokumente/Data_processing/data_paper/palm_netatmo_child.csv")