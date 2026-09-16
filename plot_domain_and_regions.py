import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import os
os.environ['ESMFMKFILE'] = 'C:/Users/eving/AppData/Local/miniconda3/envs/evalprecextremes/Library/lib/esmf.mk'
import xesmf as xe
import numpy as np
np.float_ = np.float64
import matplotlib.colors as mcolors
from clisops.core import subset
from metpy.plots import ctables
import matplotlib.ticker as ticker
import cmaps
import geopandas as gpd
from scipy import stats
import cartopy
import matplotlib.patches as mpatches

# Original Projection
proj_name  = "PlateCarree"
proj_dict  = dict()
proj       = ccrs.__dict__[proj_name].__call__(**proj_dict)

# Define the rotated pole coordinate system used by EURO-CORDEX EUR-11
rot_pole = ccrs.RotatedPole(pole_longitude=162.0, pole_latitude=39.25)

fig, ax = plt.subplots(
    subplot_kw={'projection': proj}, figsize=(10, 8)
)

# Plot the outer boundary edges of the 2D grid
fMod = f'D:/EUROCORDEX_extremes/DATA/SIM/ERA5_WRF451Q1980-2020_SON_all_1h_perc.nc'
dsMod = xr.load_dataset(fMod, engine="netcdf4")
lon = dsMod.lon
lat = dsMod.lat

ax.plot(lon[0, :], lat[0, :], color="red", transform=proj, zorder=1)  # Bottom edge
ax.plot(lon[-1, :], lat[-1, :], color="red", transform=proj, zorder=1)  # Top edge
ax.plot(lon[:, 0], lat[:, 0], color="red", transform=proj, zorder=1)  # Left edge
ax.plot(lon[:, -1], lat[:, -1], color="red", transform=proj, zorder=1)  # Right edge

ax.add_feature(cartopy.feature.OCEAN, facecolor="white", zorder=-13)
ax.add_feature(cartopy.feature.BORDERS, linewidth=0.3, color='black',zorder=0)
ax.gridlines(draw_labels=True, linewidth=0.5, color='gray',
             xlocs=range(-180,180,15), ylocs=range(-90,90,15),x_inline=False,y_inline=False)
ax.coastlines(resolution='50m', linewidth=0.3, color='black', zorder=0)


vec_obs = ["OBS_CEH-GEAR_011EUi_1990-2016", "OBS_COMEPHORE_011EUi_1997-2022",
           "OBS_CPC_011EUi_2005-2022", "OBS_GRIPHO_011EUi_2001-2016",
           "OBS_KLIMGRID_011EUi_2011-2019", "OBS_RADKLIM_011EUi_2001-2022",
           "OBS_KNMI_011EUi_1998-2022"]
label_obs = ["CEH-GEAR","COMEPHORE","CPC","GRIPHO","Klimagrid","RADKLIM","RAD_NL21"]
col_obs = {
    "#DB5856",
    "#FFAB0F",
    "#FCFC81",
    "#6ECB3C",
    "#53FCA1",
    "#8B88F8",
    "#F075E6"
}
zorder_obs = [0, 6, 1, 3, 4, 5, 2]  # Define the desired order of the legend entries

# build legend
vec_handles = []
for col, label in zip(col_obs, label_obs):
    vec_handles.append(mpatches.Patch(color=col, label=label))
leg = plt.legend(handles=vec_handles, loc='lower left', fontsize=10)
leg.set_zorder(100)

# 3. Plot the background layer (fully opaque)

for obs, col in zip(vec_obs, col_obs):
    fMod = f'D:/EUROCORDEX_extremes/DATA/OBS/{obs}_DJF_all_1h_perc.nc'
    dsMod = xr.load_dataset(fMod, engine="netcdf4")
    prec = dsMod["perc_99.90"]
    precPlot = prec.where(prec.isnull(), 1)
    single_color_cmap = mcolors.ListedColormap([col])

    mesh1 = ax.pcolormesh(precPlot.lon, precPlot.lat, precPlot.values, zorder=zorder_obs[vec_obs.index(obs)]-20,
                      vmin=0,vmax=2,cmap=single_color_cmap, shading='auto', transform=proj)

# regions prudence
gdf = gpd.read_file("regions.geojson")

gdf.plot(ax=ax, facecolor="none", edgecolor='black', alpha=0.7, transform=proj, zorder=-11)

# 4. Iterate and add names at the center of each polygon
for idx, row in gdf.iterrows():
    # Get a point guaranteed to be inside the polygon
    rep_point = row['geometry'].representative_point()
    x, y = rep_point.x, rep_point.y
    
    # Annotate the map
    labRegion = row['names']
    labRegion2lines = labRegion.replace(" ", "\n")

    ax.annotate(
        text=labRegion2lines, 
        xy=(x, y), 
        horizontalalignment='center', 
        fontsize=9, 
        weight='bold',
        color='black',
        zorder=-10
    )


# save the figure
plt.savefig(f'D:/EUROCORDEX_extremes/PAPER/Domain_and_Regions.png', 
            dpi=300, bbox_inches='tight')
