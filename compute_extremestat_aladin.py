# Description: 
# 
# Step 1: Compute precipitation maxima at different temporal and spatial scales
#
#
# Step 2: Use Christian Steger's script to compute high percentiles from NetCDF
#              climate data with high spatiotemporal resolution. Statistics are
#              computed from yearly blocks and subsequently merged. For per-
#              centiles, this is achieved by keeping the largest occurring
#              precipitation values per grid cell (and updating them during the
#              iteration through the yearly blocks). All specified percentiles
#              can then be computed from the maximal values kept in memory.
#
# Author: Guillaume Evin, Dec 2025
# partly based on the script written by Christian R. Steger, March 2023

# Load modules
import os
import numpy as np
import time
import xarray as xr
import textwrap
import warnings
from scipy import interpolate
import numba as nb

            
# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# Input/output names
files_pattern = "D:/CNRM_ERA5/pr_EUR-12_ERA5_evaluation_r1i1p1f1_CNRM-MF_CNRM-ALADIN64E1_v1-r1_1hr_" \
                + "YYYY01010030-YYYY12312330.nc"
# input file pattern -> important: label year(s) with 'YYYY'
path_out = "D:/CNRM_ERA5/"  # output directory
file_out_fp = "CNRM_ERA5_ALADIN_"
# engine used to read the file, can be "netcdf4", "h5netcdf", "zarr", "cfgrib"
myengine = 'netcdf4'
# name of the variable corresponding to precipitation
namevar = "pr"
# factor to convert input units to precipitation flux per temporal input
# frequency.
unit_con_fac = 3600.0

# Other settings
time_freq_in = "1hr"  # temporal input frequency ("1hr" or "day")

# define period of analysis
years = np.arange(2011, 2016+1)
nyear = len(years)

### SETTINGS FOR MAXIMA AT DIFFERENT TEMPORAL AND SPATIAL SCALES ###
# define duration
if time_freq_in=="1hr":
    timescale_dim = [1,3,6,12,24,72]
else:
    timescale_dim = [24,72]

ntimescale = len(timescale_dim)

# define spatial scales: we aggregate several times 3x3 cells, e.g. 1x1, 3x3, 9x9
# for 0.11deg resolution, 3x3 ~ 0.33deg, 9x9 ~ 1deg
# len_spatialscale defines how many spatial scales we consider
# for EUROCORDEX at 0.11deg, we can consider up to 3 spatial scales: 1x1, 3x3, 9x9
# CERRA: 3
# ERA5-Land: 3
# COMEPHORE: 5 
nspatialscale = 3

# -----------------------------------------------------------------------------

@nb.jit((nb.float32[:, :, :], nb.float32[:, :, :], nb.int64, nb.int64,
         nb.int64), nopython=True, parallel=True)
def update_max_values_all_day(prec_keep, prec, len_y, len_x, num_keep):
    """Update maximal precipitation values for all day/hour percentile
    calculation.

    Parameters
    ----------
    prec_keep : ndarray of float
        Array (three-dimensional) with retained precipitation data (y, x, time)
    prec : ndarray of float
        Array (three-dimensional) with precipitation data (time, y, x)
    len_y : int
        Dimension length in y-direction
    len_x : int
        Dimension length in x-direction
    num_keep : int
        Number of elements to keep"""

    for i in nb.prange(len_y):
        for j in range(len_x):
            mask = (prec[:, i, j] > prec_keep[i, j, 0])
            prec_keep[i, j, :] \
                = np.sort(np.append(prec_keep[i, j, :],
                                    prec[mask, i, j]))[-num_keep:]

############### SETTINGS FOR PERCENTILES ######################
percentile_method = "all"
# percentile method according to Schär et al. (2016) ("all", "wet")
qs = np.array([99.0, 99.9])  # percentiles
prec_thresh = {"1hr": 0.1, "day": 1.0}  # [mm]
# threshold for wet day/hour according to Ban et al. (2021)
interp_kind = "linear"
# interpolation method used to compute percentiles. Either "linear" (linear),
# "previous" (lower) or "next" (higher). The bracket values indicates the
# corresponding 'method' in np.percentile(). In Chinita et al. (2021), the
# method 'higher' is applied (-> see supplementary material)

# -----------------------------------------------------------------------------
# Preprocessing steps
# -----------------------------------------------------------------------------
print(" Preprocessing steps " .center(79, "#"))

# Check input settings
if time_freq_in not in ("1hr", "day"):
    raise ValueError("Unknown selected temporal granularity")
if percentile_method not in ("all", "wet"):
    raise ValueError("Unknown value for 'percentile_method'")
if (qs.min() < 75.0) or (qs.max() > 100.0):
    raise ValueError("Allowed range for qs of [85.0, 100.0] is exceeded")
if interp_kind not in ("linear", "previous", "next"):
    raise ValueError("Unknown value for 'interp_kind'")

# Check if all files exist
print("\n".join(textwrap.wrap(("Process files "
                               + files_pattern.split("/")[-1]), 79)))
files = [files_pattern.replace("YYYY", str(i)) for i in years]
if not all([os.path.isfile(i) for i in files]):
    raise ValueError("Not all input files exist")

# Load metadata of NetCDF files
ds = xr.open_dataset(files[0], engine=myengine)
if "calendar" in list(ds["time"].attrs):
    mod_cal = ds["time"].calendar
elif "calender" in list(ds["time"].attrs):
    mod_cal = ds["time"].calender
else:
    mod_cal = "365_day"
if mod_cal not in ("standard", "gregorian", "proleptic_gregorian", "360_day","365_day"):
    raise ValueError("Unknown calendar")
if "rlon" in list(ds.coords):
    len_x = ds.coords["rlon"].size
    len_y = ds.coords["rlat"].size
    out_dim = ("rlat", "rlon")
elif "x" in list(ds.coords):
    len_x = ds.coords["x"].size
    len_y = ds.coords["y"].size
    out_dim = ("y", "x")
else:
    raise ValueError("Unknown spatial coordinates")
ds.close()

# Compute total number of time steps
map_freq = {"1hr": "H", "day": "D"}
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    time_axis = xr.date_range(start=str(years[0]), end=str(years[-1] + 1),
                              freq=map_freq[time_freq_in], calendar=mod_cal,
                              inclusive ="left")
da = xr.DataArray(time_axis, [("time", time_axis)])

# temporal subselection for year ("yearly", "JJA", "SON", "DJF", "MAM")
for year_subsel in ["JJA", "SON", "DJF", "MAM"]:
    print("Start computation for season: " + year_subsel)

    # Compute total number of time steps
    da_sel = da.sel(time=da["time.season"] == year_subsel)

    if time_freq_in=="1hr":
        num_ts_tot_1hr = da_sel["time"].size
        num_ts_tot_day = int(num_ts_tot_1hr/24)
    else:
        num_ts_tot_day = da_sel["time"].size

    t_beg_tot = time.time()

    # Dictionary for maxima at different temporal and spatial scales (e.g., 10km, 25km, 50km, 100km)
    # Thet do not have the same spatial resolution for different spatial scales
    Xagg = {}
    xmax = {}
    coords_x = {}
    coords_y = {}
    coords_lat = {}
    coords_lon = {}

    # Allocate arrays for percentile
    num_keep_day = (np.ceil(num_ts_tot_day * (1.0 - qs.min() / 100.0)) + 1.0) \
            .astype(np.int32)
    prec_keep_day = np.empty((len_y, len_x, num_keep_day), dtype=np.float32)
    prec_keep_day.fill(-999.0)
    
    if time_freq_in=="1hr":
        num_keep_1hr = (np.ceil(num_ts_tot_1hr * (1.0 - qs.min() / 100.0)) + 1.0) \
            .astype(np.int32)
        prec_keep_1hr = np.empty((len_y, len_x, num_keep_1hr), dtype=np.float32)
        prec_keep_1hr.fill(-999.0)
    
    # Loop through years
    out_unit = {"1hr": "mm h-1", "day": "mm day-1"}  # output units
    for ind, year in enumerate(years):

        print((" Process year " + str(year) + " ").center(79, "-"))

        # Load data
        t_beg = time.time()
        file_in = files[ind]
        print(textwrap.fill("Load file " + file_in.split("/")[-1], 79))
        
        # -------------------------------------------------------------------------
        # Load data in one block:
        # Could be modified as in Christian's script to avoid loading too much data
        # however, it seems ok to load hourly precipitation data for one season
        # -------------------------------------------------------------------------
        ds = xr.open_dataset(file_in, engine=myengine)
        ds = ds.sel(time=ds["time.season"] == year_subsel)
        len_t = ds.coords["time"].size
        ds_prec = ds[namevar]*unit_con_fac  # conversion to [mm h-1] or [mm day-1]
        prec = ds_prec.values
        ds.close()
        
        print("Data loaded (" + "%.1f" % (time.time() - t_beg) + " s)")

        # -----------------------------------------------------------------------------
        # Compute maxima at different temporal and spatial scales
        # -----------------------------------------------------------------------------
        for i_timescale in range(ntimescale):
            #print('Processing temporal scale '+str(timescale_dim[i_timescale])+'h')

            # temporal scale, e.g. 6h
            timescale = timescale_dim[i_timescale]

            # tag temporal scale
            tagTS = 'T'+str(timescale)
            Xagg[tagTS] = {}
            for i_spatialscale in range(nspatialscale):

                # tag spatial scale, S1 is the original data at 1x1 resolution
                # S2 is the spatial scale at the aggregated resolution 3x3
                # S3 is the spatial scale at the aggregated resolution 9x9
                # etc.
                tagSS = 'S'+str(i_spatialscale+1)
                #print('   Processing spatial scale '+tagSS)

                # raw precipitation data
                if(i_timescale==0 and i_spatialscale==0):
                    Xagg[tagTS][tagSS] = ds_prec
                # if this is the finer time scale, we aggregate at different spatial resolution
                elif i_timescale==0:
                    # when i=1, S1 is the finer scale 1x1, we want to aggregate to 3x3, when i=2, S2 is 3x3, we want to aggregate to 9x9, etc. 
                    fine_spatialscale = 'S'+str(i_spatialscale)
                    if(out_dim == ("y", "x")):
                        Xagg[tagTS][tagSS] = Xagg['T1'][fine_spatialscale].coarsen(x=3, y=3, boundary='pad').mean()
                    else:
                        Xagg[tagTS][tagSS] = Xagg['T1'][fine_spatialscale].coarsen(rlon=3, rlat=3, boundary='pad').mean()
                # for larger time scales, we aggregate temporally from the finer time scale
                else:
                    # temporal aggregation: get sum in mm
                    temporalAgg = str(timescale)+'h'
                    dsagg = Xagg['T1'][tagSS].resample(time=temporalAgg).sum()
                    Xagg[tagTS][tagSS] = dsagg

                    # keep daily prec if not already in prec
                    if(i_spatialscale==0 and timescale==24):
                        precday = dsagg.values

                if(i_timescale==0 and ind==0):
                    # get dimensions to save in a xarray Dataset later
                    if(out_dim == ("y", "x")):
                        x_dim = Xagg[tagTS][tagSS].x
                        y_dim = Xagg[tagTS][tagSS].y
                    else:
                        x_dim = Xagg[tagTS][tagSS].rlon
                        y_dim = Xagg[tagTS][tagSS].rlat

                    nx = len(x_dim)
                    ny = len(y_dim)
                    
                    # define arrays to store aggregated data
                    # modif CC : nx and ny reversed
                    xmax[tagSS] = np.empty((nyear, ntimescale, ny, nx))

                    # store coordinates
                    coords_x[tagSS] = x_dim
                    coords_y[tagSS] = y_dim
                    coords_lat[tagSS] = Xagg[tagTS][tagSS].lat
                    coords_lon[tagSS] = Xagg[tagTS][tagSS].lon


                # get maxima for this temporal and spatial scale
                xmax[tagSS][ind, i_timescale, :, :] = Xagg[tagTS][tagSS].max(dim='time')

        # Update maximal values for percentile computation
        if time_freq_in=="1hr":
            update_max_values_all_day(prec_keep_1hr, prec, len_y, len_x, num_keep_1hr)
            update_max_values_all_day(prec_keep_day, precday, len_y, len_x, num_keep_day)
        else:                
            update_max_values_all_day(prec_keep_day, prec, len_y, len_x, num_keep_day)

    # Compute percentiles for entire period
    t_beg = time.time()
    if time_freq_in=="1hr":
        prec_per_1hr = np.empty((len(qs), len_y, len_x), dtype=np.float32)
        prec_per_1hr.fill(np.nan)
        x = np.linspace(0.0, 100.0, num_ts_tot_1hr, dtype=np.float32)
        if qs.min() < x[-num_keep_1hr]:
            raise ValueError("x-position for interpolation is out of range")
        for i in range(len_y):
            for j in range(len_x):
                f_ip = interpolate.interp1d(x[-num_keep_1hr:], prec_keep_1hr[i, j, :],
                                            bounds_error=True, kind=interp_kind,
                                            assume_sorted=True)
                prec_per_1hr[:, i, j] = f_ip(qs)
        
        prec_per_day = np.empty((len(qs), len_y, len_x), dtype=np.float32)
        prec_per_day.fill(np.nan)
        x = np.linspace(0.0, 100.0, num_ts_tot_day, dtype=np.float32)
        if qs.min() < x[-num_keep_day]:
            raise ValueError("x-position for interpolation is out of range")
        for i in range(len_y):
            for j in range(len_x):
                f_ip = interpolate.interp1d(x[-num_keep_day:], prec_keep_day[i, j, :],
                                            bounds_error=True, kind=interp_kind,
                                            assume_sorted=True)
                prec_per_day[:, i, j] = f_ip(qs)
    else:
        prec_per_day = np.empty((len(qs), len_y, len_x), dtype=np.float32)
        prec_per_day.fill(np.nan)
        print("Compute all day precipitation percentiles")
        x = np.linspace(0.0, 100.0, num_ts_tot_day, dtype=np.float32)
        if qs.min() < x[-num_keep_day]:
            raise ValueError("x-position for interpolation is out of range")
        for i in range(len_y):
            for j in range(len_x):
                f_ip = interpolate.interp1d(x[-num_keep_day:], prec_keep_day[i, j, :],
                                            bounds_error=True, kind=interp_kind,
                                            assume_sorted=True)
                prec_per_day[:, i, j] = f_ip(qs)
    print("Compute percentiles (" + "%.1f" % (time.time() - t_beg) + " s)")


    # -----------------------------------------------------------------------------
    # Save precipitation indices to NetCDF file
    # -----------------------------------------------------------------------------
    print(" Save statistics to NetCDF file ".center(79, "#"))

    # 1. ___________precipitation maxima at different temporal and spatial scales________
    for i_spatialscale in range(nspatialscale):
        tagSS = 'S'+str(i_spatialscale+1)
        print('Saving maxima for spatial scale '+tagSS)

        # create xarray dataset to save maxima
        # modif CC : nx and ny reversed
        max_SS_dict = {
            'maxpr': (('year','timescale','y', 'x'), xmax[tagSS]),
            'year': years,
            'x': coords_x[tagSS],
            'y': coords_y[tagSS],
            'lat': coords_lat[tagSS],
            'lon': coords_lon[tagSS],
            'timescale': timescale_dim
        }

        ds = xr.Dataset(
            data_vars={
                'maxpr': max_SS_dict['maxpr']
            },
            coords={
                'year': max_SS_dict['year'],
                'timescale': max_SS_dict['timescale'],
                'x': max_SS_dict['x'],
                'y': max_SS_dict['y'],
                'lat': max_SS_dict['lat'],
                'lon': max_SS_dict['lon']
            }
        )

        # Add attributes
        ds['timescale'].attrs.update({
            'long_name': 'Temporal aggregation scale',
            'units': 'hours'
        })
        
        ds['maxpr'].attrs.update({
            'long_name': 'Maximum precipitation at different temporal and spatial scales',
            'units': 'mm/hour'
        })

        # Save to NetCDF
        # modif CC : years added in the file name
        ds.to_netcdf(path_out + file_out_fp + str(year[0]) + "-" + str(year[-1]) + "_" + str(year_subsel) + "_" + 'max_scale_'+tagSS+'.nc')


    #_________________ 2. High percentiles ______________

    # Processing information and addition to output file name
    info = "Considered years: " + str(years[0]) + " - " + str(years[-1]) \
                + ", sub-yearly period: " + str(year_subsel) + ", threshold for " \
                    + "wet day: %.2f" % prec_thresh[time_freq_in] + " mm"
    
    # Save to NetCDF file
    nan_val = -999.0
    ds = xr.open_dataset(files[0], engine=myengine)
    ds = ds.drop_dims("time")
    ds.attrs["precipitation_indices"] = info

    # add encoding nan
    encoding_nan = {i: {"_FillValue": nan_val} for i in ["perc_%.2f" % i for i in qs]}
    encoding_no_nan = {i: {"_FillValue": None} for i in list({"rlon", "rlat", "x", "y", "lon", "lat"}
                            .intersection(set(ds.variables)))}

    if time_freq_in=="1hr":
        # 1hr
        for ind, q in enumerate(qs):
            name = "perc_%.2f" % q
            ds[name] = (out_dim, np.nan_to_num(prec_per_1hr[ind, :, :], nan=nan_val))
            ds[name].attrs["units"] = 'mm h-1'
        
        fn_add = str(years[0]) + "-" + str(years[-1]) + "_" + str(year_subsel) \
        + "_" + percentile_method + "_1hr_perc"
        ds.to_netcdf(path_out + file_out_fp + fn_add + ".nc",
                    encoding={**encoding_nan, **encoding_no_nan})

    # in all cases, daily scale
    for ind, q in enumerate(qs):
        name = "perc_%.2f" % q
        ds[name] = (out_dim, np.nan_to_num(prec_per_day[ind, :, :], nan=nan_val))
        ds[name].attrs["units"] = 'mm day-1'

    fn_add = str(years[0]) + "-" + str(years[-1]) + "_" + str(year_subsel) \
    + "_" + percentile_method + "_day_perc"
    ds.to_netcdf(path_out + file_out_fp + fn_add + ".nc",
                encoding={**encoding_nan, **encoding_no_nan})


    print("Total elapsed time: %.1f" % (time.time() - t_beg_tot) + " s")