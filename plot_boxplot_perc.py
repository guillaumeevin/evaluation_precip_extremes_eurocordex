import xarray as xr
import pandas as pd
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
import seaborn as sns
import matplotlib.pyplot as plt
from shapely.geometry import mapping

from utils import unrot_lon, unrot_lat

# regions prudence
gdf = gpd.read_file("regions.geojson")
vecregion = gdf['names']
vectagregion = gdf['abbrevs']

# folder with the 10-yr RL
folder_obs = "/jarre/gevin/EUROCORDEX/DATA/OBS/"
folder_sim = "/jarre/gevin/EUROCORDEX/DATA/SIM/"
folder_main = "/exec/gevin/evaluation_precip_extremes_eurocordex/"

#________________________________________________________________
# Obs dataset
#________________________________________________________________
# https://ec.europa.eu/eurostat/web/gisco/geodata/administrative-units/countries
df_border = gpd.read_file(folder_main+"DATA_FIG/CNTR_RG_20M_2024_4326.geojson")

# vector of time scales
vec_time_scale = [1,6,24]

# list of high-res products
vec_highres_obs = ["OBS_CEH-GEAR_011EUi_1990-2016", "OBS_COMEPHORE_011EUi_1997-2025",
           "OBS_CPC_011EUi_2005-2025", "OBS_GRIPHO_011EUi_2001-2016",
           "OBS_KLIMGRID_011EUi_2011-2019", "OBS_RADKLIM_011EUi_2001-2025",
           "OBS_KNMI_011EUi_1998-2025","OBS_INCA_011EUi_2013-2025"]
namecountry_highres_obs = ["United Kingdom","France","Switzerland","Italy","Denmark",
                   "Germany","Netherlands","Austria"]

# list models
vec_model = ["ERA5_ALARO1-SFX","ERA5_CNRM-ALADIN64E1","ERA5_HCLIM43-ALADIN",
             "ERA5_CCLM6-0-1","ERA5_GCOAST-AHOIB1-1", "ERA5_ICON-CLM-202407-1-1",
             "ERA5_RACMO23E","ERA5_RegCM5-0","ERA5_REMO2020-2-2",
             "ERA5_ROAM-NBS","ERA5_TSMP1-140-E","ERA5_WRF451Q-r1","OBS_ERA5_011EUi_"]
tag_model = ["AL","CN","HC","C1","GC", "IC","RA","Re","R1","RO","TS","W1","ERA"]

# reference file for the regridding
fMod = 'ERA5_WRF451Q1980-2020_SON_all_24h_perc.nc'
dsMod = xr.load_dataset(folder_sim+fMod, engine="netcdf4")
rlat_424_412 = dsMod.rlat.values
rlon_424_412 = dsMod.rlon.values

# Associer chaque valeur de la colonne 'sex' à une couleur
palette_dict = {"EMO-1": "#009E73", "High-res. products": "#56B4E9", "CERRA-Land": "#CC79A7"}
hue_order = ["EMO-1", "High-res. products", "CERRA-Land"]

# vector of seasons
vec_season = ['JJA','SON','DJF','MAM']

for season in vec_season:
    if True:
        vec_obs = ["HIGHRES"]
        tag_obs = ["High-res. products"]
        
        ##################################################################################
        # Start compiling the dataframe here
        ##################################################################################

        df = []

        for time_scale in vec_time_scale:
            ts_str = str(time_scale)
            
            if time_scale < 24:
                quantile = 'perc_99.90'
            else: 
                quantile = 'perc_99.00' 
        
            # if greater or equal to 6h, add EMO
            if(time_scale in [6,24]):
                vec_obs.append("OBS_EMO_011EUi_1990-2022")
                tag_obs.append("EMO-1")

            # if equal to 24h, add CERRA
            if(time_scale in [24]):
                vec_obs.append("OBS_CERRALAND_011EUi_1986-2020")
                tag_obs.append("CERRA-Land")

            for obs in vec_obs:
                print("=============================================================")
                print(ts_str + " / " + season + " / " + obs)        
                print("=============================================================")
                
                if obs == "HIGHRES":
                    for iHR in range(len(vec_highres_obs)):
                        # read raw data and format to be clipped
                        fHR = f'{folder_obs}{vec_highres_obs[iHR]}_{season}_all_{ts_str}h_perc.nc'
                        dsHR = xr.load_dataset(fHR, engine="netcdf4")
                        dsHR.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
                        dsHR.rio.write_crs("epsg:4326", inplace=True)

                        # for the first product we retrieve he raw data
                        if iHR==0:
                            precmerge = dsHR[quantile]
                        
                        # border for the corresponding country
                        border_country = df_border[df_border["NAME_ENGL"] == namecountry_highres_obs[iHR]]

                        # clip
                        dsClip = dsHR.rio.clip(border_country.geometry.apply(mapping), border_country.crs, drop=False)

                        prec_country = dsClip[quantile]
                        mask2add = ~np.isnan(prec_country.values)
                        precmerge.values[mask2add] = prec_country.values[mask2add]
                    
                    dsObs_cor = dsHR.cf.add_bounds(["lon","lat"])
                    dsObs_cor[quantile].values = precmerge.values
                else:
                    fObs = f'{obs}_{season}_all_{ts_str}h_perc.nc'
                    dsObs = xr.load_dataset(folder_obs+fObs, engine="netcdf4")
                    dsObs_cor = dsObs.cf.add_bounds(["lon","lat"])

                tago = tag_obs[vec_obs.index(obs)]

                for imodel, model in enumerate(vec_model):
                    # corresponding label for the model
                    tagm = tag_model[vec_model.index(model)]
                    
                    # reference file for the regridding
                    fMod = f'{model}1980-2020_{season}_all_{ts_str}h_perc.nc'
                    
                    if(tagm=="ERA"):
                        # ERA5 is in the obs folder
                        dsMod = xr.load_dataset(folder_obs+fMod, engine="netcdf4")
                    else:
                        dsMod = xr.load_dataset(folder_sim+fMod, engine="netcdf4")

                    if model == 'ERA5_RegCM5-0':
                        dsMod = dsMod.assign_coords(rlat=rlat_424_412, rlon=rlon_424_412)

                    dsMod = dsMod.cf.add_bounds(["lon","lat"])

                    # regrid obs
                    regridder = xe.Regridder(
                        dsMod, dsObs_cor, "conservative", unmapped_to_nan=True, ignore_degenerate=True
                    )
                    dsMod_regrid = regridder(dsMod, keep_attrs=True)

                    ZOBS = dsObs_cor[quantile].values
                    ZMOD = dsMod_regrid[quantile].values
                    ZOBS[ZOBS==0] = np.nan
                    Z = (ZMOD / ZOBS - 1) * 100

                    # Use indexing='ij' for matrix/geospatial indexing (Rows=Lat, Cols=Lon)
                    lon_vec = dsObs_cor.lon.values
                    lat_vec = dsObs_cor.lat.values
                    lon, lat = np.meshgrid(lon_vec, lat_vec, indexing='xy')

                    for iregion, region in enumerate(vecregion):
                        gdfregion = gdf.loc[gdf['names'] == region]
                        reg_minx, reg_miny, reg_maxx, reg_maxy = gdfregion.geometry.total_bounds

                        # Create a 2D boolean mask for the bounding box
                        box_mask = (lat >= reg_miny) & (lat <= reg_maxy) & (lon >= reg_minx) & (lon <= reg_maxx)

                        # get difference for this region in a vector
                        zvec = Z[box_mask].flatten()
                        zvec = zvec[~np.isnan(zvec)]
                        nz = len(zvec)

                        if(nz>0):
                            # append to dataframe
                            df.append(pd.DataFrame({"Timescale": [time_scale]*nz,
                                                    "Region": [region]*nz,
                                                    "Model": [tagm]*nz, 
                                                    "Obs": [tago]*nz,
                                                    "Bias": (zvec)
                                                    }))

            df_concat = pd.concat(df, ignore_index=True)
            for iregion, region in enumerate(vecregion):
                df_subplot = df_concat[(df_concat["Region"] == region) & (df_concat["Timescale"] == time_scale)]
                df_subplot.to_csv(f"{folder_main}DATA_FIG/BOXPLOT_PREC/boxplot_diff_perc_{season}_{region}_{ts_str}h.csv", index=False)

    ##################################################################################
    # Start figure here
    ##################################################################################

    # Initialiser la figure et les axes
    fig = plt.figure(figsize=(10,8.5))

    # split into three groups
    gs = fig.add_gridspec(3, 1, hspace=0.3)  # spacing between the two groups

    gs0 = gs[0].subgridspec(2, 4, hspace=0.35, wspace=0.05)  # spacing within the groups
    gs1 = gs[1].subgridspec(2, 4, hspace=0.35, wspace=0.05)
    gs2 = gs[2].subgridspec(2, 4, hspace=0.35, wspace=0.05)

    for it, time_scale in enumerate(vec_time_scale):
        ts_str = str(time_scale)
        if(it == 0):
            subaxes = gs0.subplots()
            ax_title = fig.add_subplot(gs0[:])
        elif(it == 1):
            subaxes = gs1.subplots()
            ax_title = fig.add_subplot(gs1[:])
        else:
            subaxes = gs2.subplots()
            ax_title = fig.add_subplot(gs2[:])

        # Add ghost axes and titles on gs1 and gs2
        ax_title.axis('off')
        ax_title.set_title(f"{ts_str}-hour time scale", fontsize=12, pad=5.0, fontweight='bold')

        for iregion, region in enumerate(vecregion):
            ax = subaxes[iregion//4, iregion%4]

            if ((region in ["Eastern Europe", "Iberian Peninsula", "Scandinavia"]) and (time_scale == 1)):
                ax.set_facecolor("lightgrey")
                ax.set_xticks([])
                ax.set_ylim([-100, 100])
                ax.set_yticks([-50, 0, 50])
            else:    
                # read csv file
                df_subplot = pd.read_csv(f"{folder_main}DATA_FIG/BOXPLOT_PREC/boxplot_diff_perc_{season}_{region}_{ts_str}h.csv")
                if (region in ["Eastern Europe", "Iberian Peninsula", "Scandinavia"]):
                    df_subplot = df_subplot[df_subplot["Obs"] != "High-res. products"]
                sns_plot = sns.lineplot(data=df_subplot, x='Model', y='Bias', hue_order=hue_order, hue='Obs', palette=palette_dict, ax=ax, 
                                        estimator='median', errorbar=("pi", 90), marker='o', linestyle='', err_style='bars')
                ax.get_legend().remove()
                sns_plot.set_ylim([-100, 100])
                ax.set_yticks([-50, 0, 50])
                sns_plot.grid()
                sns_plot.set_xticklabels(sns_plot.get_xticklabels(),rotation=90,fontsize=6)
            
            if (iregion in [0,4]):
                ax.set(xlabel=None, ylabel="Rel. diff. (%)")
            else:
                ax.set(xlabel=None, ylabel=None)
                ax.set_yticklabels([])

            ax.annotate(
                vectagregion[iregion],
                xy=(0, 1), xycoords='axes fraction',
                xytext=(+0.3, -0.2), textcoords='offset fontsize',
                fontsize='small', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='0.95', edgecolor='none', pad=3.0))

    handles, labels = subaxes[0, 2].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=4, bbox_to_anchor=(0.5, 0.02), frameon=False)

    plt.savefig(f"/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/boxplot_diff_perc_{season}.pdf", 
                bbox_inches="tight")