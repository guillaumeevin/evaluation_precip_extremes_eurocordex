from matplotlib import patches
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import numpy as np
import cmaps
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# folder with the 10-yr RL
folder = "/jarre/gevin/EUROCORDEX/RL/10YR/"

# vector of time scales
vec_time_scale = ['1','3','6','12','24']

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
vec_region = gdf['names']

#_________________________________________________________________
# First step: build a dataframe with the median bias per region, model and time scale
#_________________________________________________________________

if True:
    df = []

    # at an hourly scale
    for time_scale in vec_time_scale:
        # load dsRef
        fRef = f"DS10RL_REF_{time_scale}h.nc"
        dsRef = xr.load_dataset(folder+fRef, engine="netcdf4")

        for model in vec_model:
            lab_sim = vec_label[vec_model.index(model)]

            fmodel = f"DS10RL_{model}_withHershCorr_{time_scale}h.nc"
            dsModel = xr.load_dataset(folder+fmodel, engine="netcdf4")

            # compute relative difference
            lon = dsModel.lon.values
            lat = dsModel.lat.values
            Zref = dsRef.RL.values
            Zmodel = dsModel.RL.values
            Zdiff = (Zmodel - Zref)/Zref * 100

            # loop over regions to calculate a median bias per region
            for region in vec_region:
                # retrieve limits of the region
                gdfregion = gdf.loc[gdf['names'] == region]
                reg_minx, reg_miny, reg_maxx, reg_maxy = gdfregion.geometry.total_bounds

                # Create a 2D boolean mask for the bounding box
                box_mask = (lat >= reg_miny) & (lat <= reg_maxy) & (lon >= reg_minx) & (lon <= reg_maxx)

                zvec = Zdiff[box_mask].flatten()
                zvec = zvec[~np.isnan(zvec)]

                # compute median bias
                medianBias = np.median(zvec)
                    
                df.append(pd.DataFrame({"Region": (region), 
                                        "Model": [lab_sim], 
                                        "Bias": [medianBias],
                                        "Agg": (str(time_scale) + "h")
                                        }))

    df_diff = pd.concat(df, ignore_index=True)
    df_diff.to_csv("/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/df_heatmap_10yrRL.csv", index=False)


#_________________________________________________________________
# Second step: plot this median bias with one axis per time scale
#_________________________________________________________________

# read csv file
df_diff = pd.read_csv('/exec/gevin/evaluation_precip_extremes_eurocordex/DATA_FIG/df_heatmap_10yrRL.csv')

fig, axes = plt.subplots(nrows=2,ncols=3, figsize=(6.5, 3.5), constrained_layout=True)

for i, ax in enumerate(axes.flat):
    if(i<len(vec_time_scale)):
        time_scale = vec_time_scale[i]

        mat_plot = np.zeros((len(vec_region), len(vec_label)))
        mat_mask = mat_plot < 0
        for ir, region in enumerate(vec_region):
            if (region in ["Eastern Europe", "Iberian Peninsula"]) and (time_scale in ["1","3","6","12"]):
                mat_mask[ir, :] = True
            for im, model in enumerate(vec_label):
                mat_plot[ir, im] = df_diff[(df_diff['Region'] == region) &
                                        (df_diff['Agg'] == (str(time_scale) + "h")) & 
                                        (df_diff['Model'] == model)]['Bias'].values[0]
            
        im = sns.heatmap(mat_plot, 
                    cmap = cmaps.precip_diff_12lev,
                    vmin=-65, vmax=65,
                    mask=mat_mask, 
                    ax=ax,
                    cbar=False,
                    linecolor='black', linewidths=0.1,
                    xticklabels=False, yticklabels=False)
        # Set the background color to gray so masked cells appear gray
        ax.set_facecolor("lightgray")
        
        ax.margins(0)
        ax.set_aspect("equal")
        ax.tick_params(left=False, bottom=False) ## other options are right and top

        # add label model
        ax.set_title((str(time_scale) + "h"), fontweight='bold', fontsize=10)

        #add label region for the first column
        if(i%3==0):
            ax.set_yticks(np.arange(len(vec_region))+0.5)
            ax.set_yticklabels(vec_region, rotation=0, fontsize=6)

        # add label model for the last row
        ax.set_xticks(np.arange(len(vec_label))+0.5)
        ax.set_xticklabels(vec_tag, rotation=90, fontsize=6)

axes[1, 2].set_axis_off()
vec_leg_model = [vec_tag[i] + ": " + vec_label[i] for i in range(len(vec_tag))]
vec_leg_col = vec_leg_model[0:]
vec_leg_model_withsep = "\n".join(vec_leg_col)
axes[1, 2].text(0.1, 0, """
{}
""".format(vec_leg_model_withsep), fontsize=7, ha='left', va='center', transform=axes[1, 2].transAxes)

# add colorbar
mappable = im.get_children()[0]
cbar_ax = fig.add_axes([0.1, 0.2, 0.6, 0.05])
cbar = fig.colorbar(mappable, cax = cbar_ax, aspect=50, orientation="horizontal", shrink=0.8, 
                    extend='both', extendfrac=0.1)
cbar.set_label('Median relative difference (%)', fontsize=10)
cbar.ax.tick_params(labelsize=10)

# save the figure
plt.savefig(f'/exec/gevin/evaluation_precip_extremes_eurocordex/PAPER/heatmap_10yrRL.png', dpi=300, bbox_inches='tight')