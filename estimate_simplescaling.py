import xarray as xr
import xesmf as xe
from utils import getScalingParameter

folder_obs = "/jarre/gevin/EUROCORDEX/DATA/OBS/"
folder_sim = "/jarre/gevin/EUROCORDEX/DATA/SIM/"
folder_SS = "/jarre/gevin/EUROCORDEX/SIMPLESCALING/"

#______________________________________________________________________________________________________________
# estimate scaling parameters from observations
#______________________________________________________________________________________________________________
name_obs = ["OBS_CEH-GEAR_011EUi_1990-2016", "OBS_COMEPHORE_011EUi_1997-2025",
           "OBS_CPC_011EUi_2005-2025", "OBS_GRIPHO_011EUi_2001-2016",
           "OBS_KLIMGRID_011EUi_2011-2019", "OBS_RADKLIM_011EUi_2001-2025",
           "OBS_KNMI_011EUi_1998-2025","OBS_INCA_011EUi_2013-2025"]

if(True):
    for season in ["JJA","SON","DJF","MAM"]:
        for obs in name_obs:
        
            fObs = f'{obs}_{season}_max_scale_S1.nc'
            dsObs = xr.load_dataset(folder_obs+fObs, engine="netcdf4",decode_timedelta=True)
            dsObs = dsObs.cf.add_bounds(["lon","lat"])

            print(f'start estimation for season {season} and dataset {obs}')
            scalingPar = getScalingParameter(dsObs.maxpr.values)

            dsObs_cor = dsObs.drop_vars("maxpr")

            dsObs_cor["scalingPar"] = (("y", "x"), scalingPar)

            dsObs_cor.to_netcdf(f'{folder_SS}simplescaling_{obs}_{season}.nc')


#______________________________________________________________________________________________________________
# estimate scaling parameters from simulations
#______________________________________________________________________________________________________________

# reference file for the regridding
fMod = 'ERA5_WRF451Q1980-2020_SON_all_24h_perc.nc'
dsMod = xr.load_dataset(folder_sim+fMod, engine="netcdf4")
rlat_424_412 = dsMod.rlat.values
rlon_424_412 = dsMod.rlon.values

# list models and obs
name_models = ["ERA5_ALARO1-SFX","ERA5_CNRM-ALADIN64E1","ERA5_HCLIM43-ALADIN",
             "ERA5_CCLM6-0-1","ERA5_GCOAST-AHOIB1-1", "ERA5_ICON-CLM-202407-1-1",
             "ERA5_RACMO23E","ERA5_RegCM5-0","ERA5_REMO2020-2-2",
             "ERA5_ROAM-NBS","ERA5_TSMP1-140-E","ERA5_WRF451Q-r1"]

if(True):
    for season in ["JJA","SON","DJF","MAM"]:
    
        # loop over model
        for model in name_models:
            print("process data for model " + model + " and season " + season)
        
            fMod = f'{model}1980-2020_{season}_max_scale_S1.nc'
            dsMod = xr.load_dataset(folder_sim+fMod, engine="netcdf4",decode_timedelta=True)

            if model == 'ERA5_RegCM5-0':
                dsMod = dsMod.assign_coords(rlat=rlat_424_412, rlon=rlon_424_412)

            print(f'start estimation for season {season} and model {model}')
            scalingPar = getScalingParameter(dsMod.maxpr.values)

            dsMod = dsMod.cf.add_bounds(["lon","lat"])

            dsMod_cor = dsMod.drop_vars("maxpr")

            dsMod_cor["scalingPar"] = (("y", "x"), scalingPar)

            dsMod_cor.to_netcdf(f'{folder_SS}simplescaling_{model}_{season}.nc')