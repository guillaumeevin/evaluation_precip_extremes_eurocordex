import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import os
import numpy as np
np.float_ = np.float64
from metpy.units import units
import matplotlib.colors as mcolors
from clisops.core import subset
from metpy.plots import ctables
import matplotlib.ticker as ticker
import cmaps
import geopandas as gpd
import seaborn as sns
from shapely.geometry import Point, Polygon
import rioxarray
from math import nan
from scipy import stats
import statsmodels.api as sm
from collections import OrderedDict
import sys
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm


from utils import quatromatrix, get_freq_max


folder_obs = '/jarre/gevin/EUROCORDEX/DATA/OBS/'
folder_model = '/jarre/gevin/EUROCORDEX/DATA/SIM/'

# list models and obs
vec_model = ["ERA5_ALARO1-SFX","ERA5_CNRM-ALADIN64E1", "ERA5_HCLIM43-ALADIN",
            "ERA5_CCLM6-0-1","ERA5_CCLM6-0-1-URB", "ERA5_CCLM6-0-1-URB-ESG", 
            "ERA5_GCOAST-AHOIB1-1", "ERA5_ICON-CLM-202407-1-1","ERA5_RACMO23E","ERA5_RegCM5-0",
            "ERA5_REMO2020-2-2","ERA5_REMO2020-2-2-MR2","ERA5_REMO2020-2-2-TEB", 
            "ERA5_REMO2020-2-2-iMOVE", "ERA5_REMO2020-2-2-iMOVE-LUC",
            "ERA5_ROAM-NBS","ERA5_TSMP1-140-E",
            "ERA5_WRF451Q-r1","ERA5_WRF451Q-r2","ERA5_WRF451Q"]
vec_tag = ["AL","CN","HC","C1","C2", "C3","GC", "IC","RA","Re",
           "R1","R2","R3","R4", "R5","RO","TS","W1","W2","W3"]
vec_label = ["ALARO1-SFX","CNRM-ALADIN64E1","HCLIM43-ALADIN",
             "CCLM6-0-1","CCLM6-0-1-URB", "CCLM6-0-1-URB-ESG", 
           "GCOAST-AHOIB1-1", "ICON-CLM-202407-1-1","RACMO23E","RegCM5-0",
           "REMO2020-2-2","REMO2020-2-2-MR2","REMO2020-2-2-TEB", 
           "REMO2020-2-2-iMOVE", "REMO2020-2-2-iMOVE-LUC",
           "ROAM-NBS","TSMP1-140-E",
           "WRF451Q-r1","WRF451Q-r2","WRF451Q-r3"]


# regions prudence
gdf = gpd.read_file("regions.geojson")

# only three regions have high-resolution data so far
vecregion = ["France","Mid-Europe","Mediterranean","British Isles"]

if True:
    df = []

    # at an hourly scale
    for agg in [1, 3, 6, 12, 24]:

        for region in vecregion:
            for season in ["JJA","SON","DJF","MAM"]:

                if region == "France":
                    name_obs = "OBS_COMEPHORE_011EUi_1997-2025"
                elif region == "Mid-Europe":
                    name_obs = "OBS_RADKLIM_011EUi_2001-2025"
                elif region == "Mediterranean":
                    name_obs = "OBS_GRIPHO_011EUi_2001-2016"
                elif region == "British Isles":
                    name_obs = "OBS_CEH-GEAR_011EUi_1990-2016"

                gdfregion = gdf.loc[gdf['names'] == region]
                reg_minx, reg_miny, reg_maxx, reg_maxy = gdfregion.geometry.total_bounds

                # 1 for the reference
                fRef = f'{folder_obs}{name_obs}_{season}_max_scale_S1.nc'
                dsRef = xr.load_dataset(fRef, engine="netcdf4",decode_timedelta=True)

                # get lon, lat in matrices
                lat_vec = dsRef.lat.values
                lon_vec = dsRef.lon.values
                lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')

                freq_ref = get_freq_max(dsRef,agg,reg_minx, reg_miny, reg_maxx, reg_maxy, lon, lat)

                # 2 for ERA5
                name_era = "OBS_ERA5_011EUi_1980-2020"
                fEra = f'{folder_obs}{name_era}_{season}_max_scale_S1.nc'
                dsEra = xr.load_dataset(fEra, engine="netcdf4",decode_timedelta=True)

                # get lon, lat in matrices
                lat_vec = dsEra.lat.values
                lon_vec = dsEra.lon.values
                lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')

                freq_era = get_freq_max(dsEra,agg,reg_minx, reg_miny, reg_maxx, reg_maxy, lon, lat)


                # 3 for a simulation
                for name_sim in vec_model:
                    lab_sim = vec_label[vec_model.index(name_sim)]

                    fSim = f'{folder_model}{name_sim}1980-2020_{season}_max_scale_S1.nc'
                    dsSim = xr.load_dataset(fSim, engine="netcdf4",decode_timedelta=True)

                    # get lon, lat in matrices
                    lat = dsSim.lat.values
                    lon = dsSim.lon.values

                    freq_sim = get_freq_max(dsSim,agg,reg_minx, reg_miny, reg_maxx, reg_maxy, lon, lat)

                    # DAV
                    Ssim = np.sum(np.minimum(freq_ref, freq_sim))
                    Sera = np.sum(np.minimum(freq_ref, freq_era))
                    DistAddedValue = 100*(Ssim-Sera)/Sera

                    # append to dataframe
                    df.append(pd.DataFrame({"Region": (region), 
                                            "Model": [lab_sim], 
                                            "DAV": [DistAddedValue],
                                            "Season": (season),
                                            "Agg": (str(agg) + "h")
                                            }))

    df_DAV = pd.concat(df, ignore_index=True)
    df_DAV.to_csv("/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/DAV.csv", index=False)

#======================================================================
# FIGURE HEATMAP MEDIAN
#======================================================================

# read csv file
df_DAV = pd.read_csv('/exec/gevin/evaluation_precip_extremes_eurocordex/DATA_FIG/DAV.csv')
vecagg = ["1h", "3h", "6h", "12h", "24h"]

label_region = {"France":"France / COMEPHORE",
                "Mid-Europe": "Germany / RADKLIM",
                "Mediterranean": "Italy / GRIPHO",
                "British Isles": "Great Britain / CEH-GEAR"}

fig, axes = plt.subplots(nrows=2,ncols=3, width_ratios=[3, 3, 1],figsize=(6.5, 3.5), constrained_layout=True)

i_ax = [0,0,1,1]
j_ax = [0,1,0,1]

# Define discrete boundaries
bounds = [-100,-90,-80,-70,-60,-50,-40,-30,-20,-5,0,5,20,30,40,50,60,70,80,90,100]
norm = BoundaryNorm(boundaries=bounds, ncolors=cmaps.precip4_diff_19lev.N)

for iregion, region in enumerate(vecregion):
    i = i_ax[iregion]
    j = j_ax[iregion]
    ax = axes[i, j]

    dict_plot = {}
    for season in ['JJA','SON','DJF','MAM']:
        dict_plot[season] = np.zeros((len(vecagg), len(vec_label)))
        for iagg, agg in enumerate(vecagg):
            for im, model in enumerate(vec_label):
                dict_plot[season][iagg, im] = df_DAV[(df_DAV['Region'] == region) & 
                                    (df_DAV['Season'] == season) & 
                                        (df_DAV['Agg'] == agg) & 
                                        (df_DAV['Model'] == model)]['DAV'].values[0]

    
    p = quatromatrix(dict_plot['JJA'], dict_plot['MAM'], dict_plot['DJF'], dict_plot['SON'], ax=ax,
                triplotkw={"color":"k", "lw":1},
                tripcolorkw={"cmap": cmaps.precip4_diff_19lev, "norm": norm, "ec": "k"}) 

    ax.margins(0)
    ax.set_aspect("equal")

    # add label model
    ax.set_title(label_region[region], fontweight='bold', fontsize=10)

    # change the ticks and labels
    ax.set_xticks(np.arange(len(vec_tag))+0.5)
    ax.set_xticklabels(vec_tag, rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(vecagg))+0.5)
    ax.set_yticklabels(vecagg, fontsize=8)

fig.delaxes(axes[0,2])
axes[1, 2].set_axis_off()
vec_leg_model = [vec_tag[i] + ": " + vec_label[i] for i in range(len(vec_tag))]
vec_leg_col = vec_leg_model[0:]
vec_leg_model_withsep = "\n".join(vec_leg_col)
axes[1, 2].text(0, 0.2, """
{}
""".format(vec_leg_model_withsep), fontsize=4, ha='left', va='center', transform=axes[1, 2].transAxes)

# add colorbar
cbar = fig.colorbar(p, ax=axes, aspect=30, orientation="horizontal", shrink=0.8, 
                    extend='both', extendfrac=0.1, spacing='proportional')
cbar.set_label('DAV (%)', fontsize=10)
cbar.ax.tick_params(labelsize=10)

# Coordinates [left, bottom, width, height] as fractions of the figure (0 to 1)
ax2 = fig.add_axes([0.87, 0.8, 0.1, 0.1])

# 2. Cross the square with diagonal lines to form 4 triangles
ax2.plot([0, 1], [0, 1], color='black', lw=1.5) # Diagonal from bottom-left to top-right
ax2.plot([0, 1], [1, 0], color='black', lw=1.5) # Diagonal from top-left to bottom-right

# add rectangle to the axes
rect = plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor='black', lw=1.5)
ax2.add_patch(rect)

# 3. Place seasonal text in the perfect center coordinates of each triangle quadrant
ax2.text(0.5, 0.85, 'SON', fontsize=8, ha='center', va='center', weight='bold')
ax2.text(0.5, 0.15, 'MAM', fontsize=8, ha='center', va='center', weight='bold')
ax2.text(0.2, 0.5, 'JJA', fontsize=8, ha='center', va='center', weight='bold')
ax2.text(0.8, 0.5, 'DJF', fontsize=8, ha='center', va='center', weight='bold')

# Remove outer tick marks and coordinate grid axes for a clean graphic
ax2.set_xlim(-0.05, 1.05)
ax2.set_ylim(-0.05, 1.05)
ax2.axis('off')

# save the figure
plt.savefig('/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/DAV.png', dpi=300, bbox_inches='tight')