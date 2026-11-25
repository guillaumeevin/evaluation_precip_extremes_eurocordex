import xarray as xr
import cartopy.crs as ccrs
import os
os.environ['ESMFMKFILE'] = '/exec/gevin/.conda/envs/myxscen/lib/esmf.mk'
import xesmf as xe
import numpy as np
np.float_ = np.float64
import geopandas as gpd
import rioxarray
import cartopy
import cartopy.feature as cf
from shapely.geometry import mapping

from utils import unrot_lon, unrot_lat

# folder_10yr with the 10-yr RL
folder_10yr = "/jarre/gevin/EUROCORDEX/RL/10YR/"

hersh_correction = [1,1.1,1.1,1.11,1.14]

vec_time_scale = ['1','3','6','12','24']

# https://ec.europa.eu/eurostat/web/gisco/geodata/administrative-units/countries
folder_main = "/exec/gevin/evaluation_precip_extremes_eurocordex/"
df_border = gpd.read_file(folder_main+"DATA_FIG/CNTR_RG_20M_2024_4326.geojson")

#________________________________________________________________
# Ref dataset
#________________________________________________________________

if(True):
    for time_scale in vec_time_scale:

        fRefRaw = f"return_levels_10a_{time_scale}h.nc"
        dsRefRaw = xr.load_dataset(folder_10yr+fRefRaw, engine="netcdf4")
        var  = dsRefRaw[f"return_levels_10a_{time_scale}h"]
        rlat = dsRefRaw.rlat
        rlon = dsRefRaw.rlon
        pole_lon = -162
        pole_lat = 39.25

        lon = unrot_lon(rlat, rlon, pole_lat, pole_lon)
        lat = unrot_lat(rlat, rlon, pole_lat, pole_lon)

        dsRef = xr.Dataset(
            data_vars={"RL": (("y", "x"), var.values)},
            coords={
                "x": rlon.values,
                "y": rlat.values,
                "lon": (("y", "x"), lon),
                "lat": (("y", "x"), lat),
            }
        )
        dsRef = dsRef.cf.add_bounds(["lon","lat"])

        f10 = f"DS10RL_REF_{time_scale}h.nc"
        dsRef.to_netcdf(folder_10yr+f10, engine="netcdf4")

#________________________________________________________________
# Obs dataset
#________________________________________________________________
# Step 1: regrid all datasets to the reference
vec_highres_obs = ["OBS_CEH-GEAR_011EUi_1990-2016", "OBS_COMEPHORE_011EUi_1997-2025",
           "OBS_CPC_011EUi_2005-2025", "OBS_GRIPHO_011EUi_2001-2016",
           "OBS_KLIMGRID_011EUi_2011-2019", "OBS_RADKLIM_011EUi_2001-2025",
           "OBS_KNMI_011EUi_1998-2025","OBS_INCA_011EUi_2013-2025"]
namecountry_highres_obs = ["United Kingdom","France","Switzerland","Italy","Denmark",
                   "Germany","Netherlands","Austria"]

if(True):
    vec_obs = vec_highres_obs.copy()
    vec_obs.append("OBS_ERA5_011EUi_1980-2020")
    for time_scale in vec_time_scale:
        # if greater or equal to 6h, add EMO
        if(time_scale in ['6','12','24']):
            vec_obs.append("OBS_EMO_011EUi_1990-2025")

        # if equal to 24h, add CERRA
        if(time_scale in ['24']):
            vec_obs.append("OBS_CERRALAND_011EUi_1986-2020")

        # load dsRef
        fRef = f"DS10RL_REF_{time_scale}h.nc"
        dsRef = xr.load_dataset(folder_10yr+fRef, engine="netcdf4")

        # Hershfield correction
        HershCorr = hersh_correction[vec_time_scale.index(time_scale)]

        if(time_scale in ['6','12','24']):
            vec_obs = ["OBS_EMO_011EUi_1990-2025"]
        else:
            vec_obs = []

        for obs in vec_obs:
            filecsv_obs = f"matrix_10yrRL_{time_scale}hr_{obs}.csv"
            mat_obs = np.genfromtxt(folder_10yr+filecsv_obs, delimiter=',', missing_values='NA', filling_values=np.nan)*HershCorr

            fObsRaw = f'/jarre/gevin/EUROCORDEX/DATA/OBS/{obs}_JJA_all_24h_perc.nc'
            dsObsRaw = xr.load_dataset(fObsRaw, engine="netcdf4")            

            # 1. Define your 1D coordinate vectors
            lat_vec = dsObsRaw.lat.values
            lon_vec = dsObsRaw.lon.values

            # 2. Create the 2D coordinate matrices
            # Use indexing='ij' for matrix/geospatial indexing (Rows=Lat, Cols=Lon)
            lon_mat, lat_mat = np.meshgrid(lon_vec, lat_vec, indexing='xy')

            # 3; buil clean dataset
            dsObs2regrid = xr.Dataset(
                data_vars={"RL": (("lat", "lon"), mat_obs),
                           "x": (("lat", "lon"), lon_mat),
                           "y": (("lat", "lon"), lat_mat)},
                coords={
                    "lon": lon_vec,
                    "lat": lat_vec,
                }
            )
            dsObs2regrid = dsObs2regrid.cf.add_bounds(["lon","lat"])

            # 4: if high-resolution dataset, clip to the country
            if (obs in vec_highres_obs):
                dsObs2regrid.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
                dsObs2regrid.rio.write_crs("epsg:4326", inplace=True)
                
                # border for the corresponding country
                country = namecountry_highres_obs[vec_highres_obs.index(obs)]
                border_country = df_border[df_border["NAME_ENGL"] == country]

                # clip
                dsObs2regrid = dsObs2regrid.rio.clip(border_country.geometry.apply(mapping), border_country.crs, drop=False)

            # 5. Regrid to the reference dataset
            # regrid obs
            regridder = xe.Regridder(
                dsObs2regrid, dsRef, "conservative", unmapped_to_nan=True, ignore_degenerate=True
            )
            dsObs = regridder(dsObs2regrid, keep_attrs=True)

            fobs = f"DS10RL_{obs}_withHershCorr_{time_scale}h.nc"
            dsObs.to_netcdf(folder_10yr+fobs, engine="netcdf4")

# Step 2: merge high-resolution datasets
for time_scale in vec_time_scale:
    for i_obs,obs in enumerate(vec_highres_obs):
        # read raw data and format to be clipped
        fobs = f"DS10RL_{obs}_withHershCorr_{time_scale}h.nc"
        dsObs = xr.load_dataset(folder_10yr+fobs, engine="netcdf4")
        
        # for the first product we retrieve he raw data
        if i_obs==0:
            RL_merge = dsObs.RL.values
        else:
            RL_country = dsObs.RL.values
            mask2add = ~np.isnan(RL_country)
            RL_merge[mask2add] = RL_country[mask2add]

    dsObs_cor = dsObs.cf.add_bounds(["lon","lat"])
    dsObs_cor.RL.values = RL_merge

    fobs = f"DS10RL_HIGHRES_withHershCorr_{time_scale}h.nc"
    dsObs_cor.to_netcdf(folder_10yr+fobs, engine="netcdf4")

#________________________________________________________________
# Model dataset
#________________________________________________________________
# list models
vec_model = ["ERA5_ALARO1-SFX","ERA5_CCLM6-0-1","ERA5_CCLM6-0-1-URB", "ERA5_CCLM6-0-1-URB-ESG", 
           "ERA5_CNRM-ALADIN64E1","ERA5_GCOAST-AHOIB1-1", "ERA5_HCLIM43-ALADIN",
           "ERA5_ICON-CLM-202407-1-1","ERA5_RACMO23E","ERA5_RegCM5-0",
           "ERA5_REMO2020-2-2","ERA5_REMO2020-2-2-MR2","ERA5_REMO2020-2-2-TEB", 
           "ERA5_REMO2020-2-2-iMOVE", "ERA5_REMO2020-2-2-iMOVE-LUC",
           "ERA5_ROAM-NBS","ERA5_TSMP1-140-E",
           "ERA5_WRF451Q-r1","ERA5_WRF451Q-r2","ERA5_WRF451Q"]

if(False):
    # reference file for the regridding
    fMod = 'ERA5_WRF451Q1980-2020_SON_all_24h_perc.nc'
    folder_10yr_model = '/jarre/gevin/EUROCORDEX/DATA/SIM/'
    dsMod = xr.load_dataset(folder_10yr_model+fMod, engine="netcdf4")
    rlat_424_412 = dsMod.rlat.values
    rlon_424_412 = dsMod.rlon.values

    for time_scale in vec_time_scale:

        # load dsRef
        fRef = f"DS10RL_REF_{time_scale}h.nc"
        dsRef = xr.load_dataset(folder_10yr+fRef, engine="netcdf4")

        # Hershfield correction
        HershCorr = hersh_correction[vec_time_scale.index(time_scale)]

        for model in vec_model:
            filecsv_model = f"matrix_10yrRL_{time_scale}hr_{model}.csv"
            mat_model = np.genfromtxt(folder_10yr+filecsv_model, delimiter=',', missing_values='NA', filling_values=np.nan)*HershCorr

            fModelRaw = f"{folder_10yr_model}{model}1980-2020_JJA_all_24h_perc.nc"
            dsModelRaw = xr.load_dataset(fModelRaw, engine="netcdf4")

            if model == 'ERA5_RegCM5-0':
                dsModelRaw = dsModelRaw.assign_coords(rlat=rlat_424_412, rlon=rlon_424_412)

            if 'rlat' in list(dsModelRaw.coords._names):
                x = dsModelRaw.rlon.values
                y = dsModelRaw.rlat.values
            else:
                x = dsModelRaw.x.values
                y = dsModelRaw.y.values

            # buil clean dataset
            dsModel2regrid = xr.Dataset(
                data_vars={"RL": (("y", "x"), mat_model)},
                coords={
                    "x": x,
                    "y": y,
                    "lon": (("y", "x"), dsModelRaw.lon.values),
                    "lat": (("y", "x"), dsModelRaw.lat.values),
                }
            )
            dsModel2regrid = dsModel2regrid.cf.add_bounds(["lon","lat"])

            # 4. Regrid to the reference dataset
            # regrid obs
            regridder = xe.Regridder(
                dsModel2regrid, dsRef, "conservative", unmapped_to_nan=True, ignore_degenerate=True
            )
            dsModel = regridder(dsModel2regrid, keep_attrs=True)

            fmodel = f"DS10RL_{model}_withHershCorr_{time_scale}h.nc"
            dsModel.to_netcdf(folder_10yr+fmodel, engine="netcdf4")