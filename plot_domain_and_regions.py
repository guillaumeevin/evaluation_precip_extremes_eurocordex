import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import os
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
from shapely.geometry import mapping
import rioxarray
import regionmask
import matplotlib.patheffects as pe
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
cartopy.config['data_dir'] = '/exec/gevin/'

# read orography of CRCM simulations
ds_oro = xr.open_dataset("/jarre/gevin/EUROCORDEX/DATA/orog_EUR-12_ERA5_evaluation_r1i1p1f1_CNRM-MF_CNRM-ALADIN64E1_v1-r1_fx.nc")
level_oro = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500]
norm_oro = mcolors.BoundaryNorm(boundaries=level_oro, ncolors=256, clip=True)

# https://ec.europa.eu/eurostat/web/gisco/geodata/administrative-units/countries
folder_obs = "/jarre/gevin/EUROCORDEX/DATA/OBS/"
folder_sim = "/jarre/gevin/EUROCORDEX/DATA/SIM/"
folder_main = "/exec/gevin/evaluation_precip_extremes_eurocordex/"
df_border = gpd.read_file(folder_main+"DATA_FIG/CNTR_RG_20M_2024_4326.geojson")

# high-res products
vec_obs = ["OBS_CEH-GEAR_011EUi_1990-2016", "OBS_COMEPHORE_011EUi_1997-2025",
           "OBS_CPC_011EUi_2005-2025", "OBS_GRIPHO_011EUi_2001-2016",
           "OBS_INCA_011EUi_2013-2025", "OBS_KLIMGRID_011EUi_2011-2019",
           "OBS_RADKLIM_011EUi_2001-2025", "OBS_KNMI_011EUi_1998-2025"]
namecountry_highres_obs = ["United Kingdom","France","Switzerland","Italy","Austria","Denmark",
                   "Germany","Netherlands"]
label_obs = ["CEH-GEAR","COMEPHORE","CPC","GRIPHO","INCA","Klimagrid","RADKLIM","RAD_NL21"]

# colors for products
col_obs = {
    "#DB5856",
    "#FFAB0F",
    "#FCFC81",
    "#6ECB3C",
    "#53FCA1",
    "#64CDFA",
    "#8B88F8",
    "#F075E6"
}

# Plot the outer boundary edges of the 2D grid
fMod = f'{folder_sim}ERA5_WRF451Q1980-2020_SON_all_1h_perc.nc'
dsMod = xr.load_dataset(fMod, engine="netcdf4")
lon = dsMod.lon
lat = dsMod.lat
extent = [lon.min(), lon.max(), lat.min(), lat.max()]

# Original Projection
projorig_name  = "PlateCarree"
projorig_dict  = dict()
projorig       = ccrs.__dict__[projorig_name].__call__(**projorig_dict)

# Graphic projection
projplot_name       = "LambertConformal"
projplot_dict  = dict(central_longitude=10.5, central_latitude=49.5, standard_parallels=(33,45))
projplot       = ccrs.__dict__[projplot_name].__call__(**projplot_dict)

fig, axes = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': projplot}, figsize=(8.5,5))


#####################################################################
# Second plot, domain + regions
ax = axes[0]

# Add cartographic background visual anchors
ax.add_feature(cfeature.LAND, facecolor="#f7f7f7", zorder=1)
ax.add_feature(cfeature.OCEAN, facecolor="#e0f3f8", zorder=1)
ax.add_feature(cfeature.COASTLINE, edgecolor="#555555", linewidth=0.8, zorder=2)
ax.add_feature(cfeature.BORDERS, edgecolor="#aaaaaa", linewidth=0.5, zorder=3)

# build legend
vec_handles = []
for col, label in zip(col_obs, label_obs):
    vec_handles.append(mpatches.Patch(color=col, label=label))
leg = ax.legend(handles=vec_handles, loc='upper left', fontsize=6,bbox_to_anchor=(0.1,1))
leg.set_zorder(100)

# 3. Plot the background layer (fully opaque)
for obs, col in zip(vec_obs, col_obs):
    fMod = f'{folder_obs}{obs}_DJF_all_1h_perc.nc'
    dsMod = xr.load_dataset(fMod, engine="netcdf4")
    dsMod.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    dsMod.rio.write_crs("epsg:4326", inplace=True)
    
    # border for the corresponding country
    country = namecountry_highres_obs[vec_obs.index(obs)]
    border_country = df_border[df_border["NAME_ENGL"] == country]

    # clip
    dsClip = dsMod.rio.clip(border_country.geometry.apply(mapping), border_country.crs, drop=False)
    
    prec = dsClip["perc_99.90"]
    precPlot = prec.where(prec.isnull(), 1)
    single_color_cmap = mcolors.ListedColormap([col])

    mesh1 = ax.pcolormesh(precPlot.lon, precPlot.lat, precPlot.values, zorder=4,
                      vmin=0,vmax=2,cmap=single_color_cmap, shading='auto', transform=projorig)

# 4. Tracer les contours des régions PRUDENCE
# Définition du style du texte des labels pour assurer leur lisibilité
text_kws = dict(
    bbox=dict(color="white", alpha=0.7, boxstyle="round,pad=0.2"),
    path_effects=[pe.withStroke(linewidth=2, foreground="white")],
    color="k",
    fontsize=10,
    weight='bold', 
    zorder=5
)

# 1. Charger les régions officielles du projet PRUDENCE
prudence_regions = regionmask.defined_regions.prudence

prudence_regions.plot(
    ax=ax,
    projection=projplot,
    add_ocean=False,
    resolution="50m",
    label="abbrev",          # Affiche les abréviations (BI, IP, FR, etc.)
    text_kws=text_kws,
    line_kws=dict(color='k', linewidth=1.5, linestyle='--', zorder=5)
)

# Display grid lines and cosmetics
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, color='#cccccc', linestyle='--',  zorder=3)
gl.xlabel_style = {'size': 6, 'color': 'k'}  # Taille 6 points
gl.ylabel_style = {'size': 6, 'color': 'k'}  # Taille 6 points


# Add the "a" label to the top-left corner
ax.text(0.02, 0.98, 'a)', 
    transform=ax.transAxes, 
    fontsize=12, 
    fontweight='bold', 
    va='top', 
    ha='left',
    zorder=10)

ax.set_extent([-10, 30, 35, 72], crs=ccrs.PlateCarree())

#####################################################################
# First plot, domain + regions
ax = axes[1]
p0 = ax.pcolormesh(ds_oro['lon'], ds_oro['lat'], ds_oro["orog"], cmap="terrain", norm=norm_oro, transform=projorig, zorder=0)
ax.add_feature(cfeature.BORDERS, linewidths=0.5, zorder=9)
ax.add_feature(cfeature.LAKES,edgecolor='k', facecolor='#e0f3f8', zorder=9, linewidths=0.5)
ax.add_feature(cfeature.OCEAN, facecolor='#e0f3f8', zorder=9)

# Add the "a" label to the top-left corner
ax.text(0.02, 0.98, 'b)', 
        transform=ax.transAxes, 
        fontsize=12, 
        fontweight='bold', 
        va='top', 
        ha='left', zorder=10)
ax.set_extent([-10, 30, 35, 72], crs=ccrs.PlateCarree())

# 2. Create the inset axes inside the map frame
# width/height can be % of the parent map or absolute inches
cax = inset_axes(
    ax,
    width="3%",       # Width of the colorbar
    height="30%",       # Height of the colorbar
    loc="upper left",  # Location inside the map
    bbox_to_anchor=(0.2,-0.03,1,1), 
    bbox_transform=ax.transAxes
)

# 3. Draw the colorbar in the inset axis
cbar = fig.colorbar(p0, cax=cax, orientation="vertical", ticks=level_oro, format="%d")
cbar.ax.tick_params(labelsize=6)  # Set to a smaller number like 8
cbar.ax.set_title('Elev. (m)', loc = "left", pad = 5, fontsize=6)

# save the figure
plt.savefig(f'/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/Domain_and_Regions.png', 
            dpi=300, bbox_inches='tight')
