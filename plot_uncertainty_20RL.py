import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import os
os.environ['ESMFMKFILE'] = '/exec/gevin/.conda/envs/myxscen/lib/esmf.mk'
import xesmf as xe
import numpy as np
np.float_ = np.float64
import matplotlib.colors as mcolors
from clisops.core import subset
from metpy.plots import ctables
import matplotlib.ticker as ticker
import cmaps
import geopandas as gpd
import seaborn as sns
import rioxarray
from math import nan
from scipy import stats
import statsmodels.api as sm
import cartopy
from shapely.geometry import mapping
cartopy.config['data_dir'] = '/exec/gevin/'
os.environ["XDG_CACHE_HOME"] = "/exec/gevin/"
from pyproj import CRS
import copy

# folder_obs and folder_sim contain the maxima aggregated at different spatial scales
# folder_RL20 contains the 20-year return level estimates (obs and sim)
from utils import get_lon_lat

folder_obs = "/jarre/gevin/EUROCORDEX/DATA/OBS/"
folder_sim = "/jarre/gevin/EUROCORDEX/DATA/SIM/"
folder_RL20 = "/jarre/gevin/EUROCORDEX/RL/20YR/"
folder_main = "/exec/gevin/evaluation_precip_extremes_eurocordex/"

# regions PRUDENCE
gdf = gpd.read_file("regions.geojson")

# list models and labels
vec_model = ["ERA5_ALARO1-SFX","ERA5_CNRM-ALADIN64E1","ERA5_HCLIM43-ALADIN","ERA5_CCLM6-0-1",
             "ERA5_GCOAST-AHOIB1-1", "ERA5_ICON-CLM-202407-1-1",
             "ERA5_RACMO23E","ERA5_RegCM5-0","ERA5_REMO2020-2-2",
             "ERA5_ROAM-NBS","ERA5_TSMP1-140-E","ERA5_WRF451Q-r1"]
vec_label = ["ALARO1-SFX","ALADIN64E1","HCLIM43","CCLM6-0-1",
             "GC-AHOIB1", "ICON-CLM", "RACMO23E","RegCM5-0","REMO2020-2-2",
             "ROAM-NBS","TSMP1-140-E","WRF451Q-r1"]

# vector of regions PRUDENCE


vec_region_PRUDENCE = ["France","Mid-Europe","Mediterranean","British Isles"]
vec_tag_region = ["France","Germany","Italy","Great Britain"]
vec_prod = ["OBS_COMEPHORE_011EUi_1997-2025","OBS_RADKLIM_011EUi_2001-2025",
            "OBS_GRIPHO_011EUi_2001-2016","OBS_CEH-GEAR_011EUi_1990-2016"]

##############################################################################
# STEP 1: Compute spatial medians for different regions and different bootstrap 
# estimates of the 20-year return level (RL20) from observations and simulations
##############################################################################
if(False):
    for iregion in [2,3]: # range(len(vec_region_PRUDENCE)):
        region_prudence = vec_region_PRUDENCE[iregion]
        tag_region = vec_tag_region[iregion]
        prod = vec_prod[iregion]

        # select one region
        region_gdf = gdf.loc[gdf['names'] == region_prudence]
        reg_minx, reg_miny, reg_maxx, reg_maxy = region_gdf.geometry.total_bounds

        # initialise list
        dict_median = {}
            
        # initialise vector time/space
        vec_temporal_scale = ["1hr","3hr","6hr","12hr","24hr","72hr"]
        vec_spatial_scale = ["S1","S2","S3","S4"]

        for model, label in zip(vec_model, vec_label):
            for season in ["JJA","SON","DJF","MAM"]:

                df = []

                for itime, temporal_scale in enumerate(vec_temporal_scale):
                    for ispace, spatial_scale in enumerate(vec_spatial_scale):
                        # RL20 from obs
                        file_RL20_obs = f"matrix_20yrRL_{temporal_scale}_{spatial_scale}_{season}_{prod}.nc"
                        ds_RL20_obs = xr.load_dataset(folder_RL20+file_RL20_obs, engine="netcdf4", decode_timedelta=False)
                        RL20_obs = ds_RL20_obs.RL20_boot.values
                        RL20_obs[RL20_obs<0] = np.nan

                        # RL20 from sim
                        file_RL20_sim = f"matrix_10yrRL_{temporal_scale}_{spatial_scale}_{season}_{model}_011EUi.csv"
                        RL20_sim = np.loadtxt(fname=folder_RL20+file_RL20_sim, delimiter=',')
                        RL20_sim[RL20_sim<0] = np.nan

                        # get lon/lat reg
                        file_obs = f"{prod}_JJA_max_scale_{spatial_scale}.nc"
                        lon, lat = get_lon_lat(folder_obs+file_obs)

                        # Create a 2D boolean mask for the bounding box
                        box_mask = (lat >= reg_miny) & (lat <= reg_maxy) & (lon >= reg_minx) & (lon <= reg_maxx)

                        # first loop over boot estimates
                        for iboot in range(RL20_obs.shape[0]):
                            Z = (RL20_sim / RL20_obs[iboot,:,:] - 1) * 100
                            zvec = Z[box_mask].flatten()
                            zmed = np.nanmedian(zvec)

                            # append to dataframe
                            df.append(pd.DataFrame({"Timescale": [temporal_scale],
                                                    "Spatialscale": [spatial_scale],
                                                    "Bias": [zmed]
                                                    }))

                df_concat = pd.concat(df, ignore_index=True)
                df_concat.to_csv(f"{folder_main}DATA_FIG/RL20_BOOT/uncertainty_RL20diff_{tag_region}_{season}_{label}.csv", index=False)


##############################################################################
# STEP 2: Produce plots
##############################################################################

nmodel = len(vec_label)
nseason = 4

hue_order = ["S4","S3","S2","S1"]

vec_season = ["JJA","SON","DJF","MAM"]

for iregion in range(len(vec_tag_region)):
    tag_region = vec_tag_region[iregion]

    fig, axes = plt.subplots(nrows=nseason,ncols=nmodel,figsize=(15, 5), constrained_layout=True)

    for i, season in enumerate(vec_season):
        for j, label in enumerate(vec_label):
            df_concat = pd.read_csv(f"{folder_main}DATA_FIG/RL20_BOOT/uncertainty_RL20diff_{tag_region}_{season}_{label}.csv")
            ax = axes[i,j]
            
            sns_plot = sns.lineplot(data=df_concat, x='Timescale', y='Bias', hue_order=hue_order, hue='Spatialscale', ax=ax, 
                                    estimator='median', errorbar=("pi", 90), marker='',linestyle='-', linewidth = 0.5, palette='tab10')
            sns_plot.grid(visible=True, which='major', axis='y', color='lightgrey', linestyle='-', linewidth=0.5)
            ax.get_legend().remove()
            sns_plot.set_ylim([-65, 65])
            ax.set_yticks([-50,-25,0,25,50])
            ax.set(xlabel=None, ylabel=None)
            ax.axhline(0, linewidth=1, color='black', linestyle='--')

            # Get the confidence interval polygons and adjust transparency
            for artist in ax.collections:
                artist.set_alpha(0.6)  # Adjust alpha for CI only

            if (i==nseason-1):
                ax.set_xticklabels(['T1','T2','T3','T4','T5','T6'], fontsize=8)
            else:
                ax.set_xticklabels([])

            if (j==0):
                ax.set(xlabel=None)
                ax.set_ylabel(season, fontsize=10)
                ax.set_yticklabels([-50,-25,0,25,50],fontsize=8)
            else:
                ax.set(xlabel=None, ylabel=None)
                ax.set_yticklabels([])

    handles, labels = ax.get_legend_handles_labels()
    # 2. Copier les handles pour ne pas modifier les tracés principaux
    handles = [copy.copy(h) for h in handles]
    order = [3, 2, 1, 0]  # Ordre souhaité des étiquettes dans la légende
    for h in handles:
        h.set_linewidth(5.0)  # Changez 5.0 par l'épaisseur voulue
    fig.legend([handles[idx] for idx in order],[labels[idx] for idx in order], 
    loc='outside lower center', ncol=4, bbox_to_anchor=(0.5, -0.05), frameon=False)

    for ax, col in zip(axes[0], vec_label):
        ax.set_title(col, fontsize=10)

    plt.savefig(f"/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/uncertainty_RL20diff_{tag_region}.pdf", 
                bbox_inches="tight")