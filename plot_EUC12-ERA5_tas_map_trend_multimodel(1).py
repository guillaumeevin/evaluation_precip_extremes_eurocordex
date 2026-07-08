#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Script to plot EURO-CORDEX-CMIP6 data.
"""

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import argparse
import json

# from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import datetime
import xarray as xr
import pandas as pd
from cartopy.util import add_cyclic_point
from collections import OrderedDict
from cmcrameri import cm

###############################
####   Input Data to plot   ###
###############################

# General characteristics
var          = 'tas'
freq         = 'mon'
domain       = 'EUR-12'
expid        = 'evaluation'
driver       = 'ERA5'
driver_mb    = 'r1i1p1f1'
name_simu    = driver+"_"+expid+"_"+driver_mb
y0           = 1950
y1           = 2023
y0P          = 1981
y1P          = 2019
nban         = y1-y0+1
season       = 'DJF'
months       = [0,1,11]
#season       = 'JJA'
#months       = [5,6,7]
#season       = 'ANN'
#months       = [0,1,2,3,4,5,6,7,8,9,10,11]
nbm          = len(months)

dir_simu     = "/cnrm/mosca/USERS/nabat/DATA/EURO-CORDEX-CMIP6/"
dir_save     = dir_simu+"save_calcul_python/"
#dir_save     = "/home/nabat/ALADIN/EURO-CORDEX/save_calcul_python/"


# Original Projection
projorig_name  = "PlateCarree"
projorig_dict  = dict()
projorig       = ccrs.__dict__[projorig_name].__call__(**projorig_dict)

# Input data
# EURO-CORDEX RCMs
RCMs = OrderedDict() 
RCMs['CNRM-ALADIN64E1'] = {'id':'CNRM-ALADIN64E1',   'inst':'CNRM-MF','vr':'v1-r1',
                                                'projplot' : 'LambertConformal',
                                                'projplot_dict' : dict(central_longitude=10.5, central_latitude=49.5, standard_parallels=(33,45), false_easting = 2925000., false_northing = 2925000.),
                                                'pers':[(1959,1960),(1961,1970),(1971,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020),(2021,2022)]}
RCMs['HCLIM43-ALADIN'] = {'id':'HCLIM43-ALADIN',   'inst':'HCLIMcom-SMHI','vr':'v1-r1',
                                                'projplot' : 'LambertConformal',
                                                'projplot_dict' : dict(central_longitude=10.5, central_latitude=49.5, standard_parallels=(33,45), false_easting = 2925000., false_northing = 2925000.),
                                                'pers':[(1980,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020),(2021,2023)]}
RCMs['ALARO1-SFX'] = {'id':'ALARO1-SFX',   'inst':'RMIB-UGent','vr':'v1-r1',
                                                'projplot' : 'LambertConformal',
                                                'projplot_dict' : dict(central_longitude=9.9, central_latitude=49.0, standard_parallels=(33,45), false_easting = 3012.5, false_northing = 3012.5),
                                                'pers':[(1980,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020),(2021,2022)]}
RCMs['RACMO23E']   = {'id':'RACMO23E',          'inst':'KNMI','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1980,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020)]}
RCMs['ICON-CLM-202407-1-1'] = {'id':'ICON-CLM-202407-1-1',   'inst':'CLMcom-Hereon','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1950,1950),(1951,1960),(1961,1970),(1971,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020),(2021,2023)]}
RCMs['CCLM6-0-1'] = {'id':'CCLM6-0-1',         'inst':'CLMcom-Hereon','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1960,1960),(1961,1970),(1971,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020),(2021,2021)]}
RCMs['CCLM6-0-1-URB-ESG'] = {'id':'CCLM6-0-1-URB-ESG',   'inst':'CLMcom-KUL','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1979,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020),(2021,2023)]}
RCMs['CCLM6-0-1-URB'] = {'id':'CCLM6-0-1-URB',   'inst':'CLMcom-CMCC','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1980,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020),(2021,2021)]}
RCMs['GCOAST-AHOIB1-1'] = {'id':'GCOAST-AHOIB1-1',         'inst':'CLMcom-Hereon','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1959,1960),(1961,1970),(1971,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2018)]}
RCMs['RegCM5-0'] = {'id':'RegCM5-0',   'inst':'ICTP','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=198, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1970,1970),(1971,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2019)]}
RCMs['REMO2020'] = {'id':'REMO2020',   'inst':'GERICS','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1979,1988),(1989,1998),(1999,2008),(2009,2018),(2019,2020)]}
RCMs['REMO2020-2-2'] = {'id':'REMO2020-2-2',   'inst':'GERICS','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1979,1988),(1989,1998),(1999,2008),(2009,2018),(2019,2020)]}
RCMs['REMO2020-2-2-MR2'] = {'id':'REMO2020-2-2-MR2',   'inst':'GERICS','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1979,1988),(1989,1998),(1999,2008),(2009,2018),(2019,2020)]}
RCMs['REMO2020-2-2-iMOVE'] = {'id':'REMO2020-2-2-iMOVE',   'inst':'GERICS','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1979,1988),(1989,1998),(1999,2008),(2009,2018),(2019,2020)]}
RCMs['REMO2020-2-2-iMOVE-LUC'] = {'id':'REMO2020-2-2-iMOVE-LUC',   'inst':'GERICS','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1979,1988),(1989,1998),(1999,2008),(2009,2018),(2019,2020)]}
RCMs['REMO2020-2-2-TEB'] = {'id':'REMO2020-2-2-TEB',   'inst':'GERICS','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1979,1988),(1989,1998),(1999,2008),(2009,2018),(2019,2020)]}
RCMs['WRF451Q-r1'] = {'id':'WRF451Q',   'inst':'IDL-FCUL','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1980,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020)]}
RCMs['WRF451Q-r2'] = {'id':'WRF451Q',   'inst':'CESAM-UA','vr':'v1-r2',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1980,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020)]}
#RCMs['WRF451Q-r3'] = {'id':'WRF451Q',   'inst':'AUTH','vr':'v1-r3',
#                                                'projplot' : 'RotatedPole',
#                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
#                                                'pers':[(1980,1980),(1981,1990),(1991,2000)]}
RCMs['ROAM-NBS'] = {'id':'ROAM-NBS',   'inst':'DWD-BSH','vr':'v1-r1',
                                                'projplot' : 'RotatedPole',
                                                'projplot_dict' : dict(pole_longitude=-162, pole_latitude=39.25, central_rotated_longitude=-162),
                                                'pers':[(1979,1980),(1981,1990),(1991,2000),(2001,2010),(2011,2020),(2021,2021)]}
    
# Get the data
for model in RCMs:
    # First file
    tper0 = RCMs[model]['pers'][0]
    filein0    = os.path.abspath(dir_simu+"/"+RCMs[model]['id']+"/"+var+"_"+domain+"_"+name_simu+"_"+RCMs[model]['inst']+'_'+RCMs[model]['id']+'_'+RCMs[model]['vr']+'_'+freq+"_"+str(tper0[0])+"01-"+str(tper0[1])+"12.nc")
    with xr.open_dataset(filein0,decode_times=False) as dataset0:
        RCMs[model]['lat'] = dataset0.variables["lat"][:,:]
        RCMs[model]['lon'] = dataset0.variables["lon"][:,:]
    # SAve file
    filesave = os.path.abspath(dir_save+"/"+"save_EUC12-ERA5_"+model+"_"+var+"_"+season+"_"+str(y0P)+"-"+str(y1P)+"_reglin.data")
    if os.path.isfile(filesave+".npy"):
        # Data already computed
        print(model+" already computed")
        RCMs[model]['reglin'] = np.load(filesave+".npy",allow_pickle=False)
    else:
        # Data not computed yet
        print(model+" to compute")
        nlat,nlon = RCMs[model]['lat'].shape
        RCMs[model]['reglin'] = np.zeros((nlat,nlon))
        RCMs[model]['data'] = np.zeros((nban,nbm,nlat,nlon))
        RCMs[model]['data'][:,:,:,:] = np.nan
        for tper in RCMs[model]['pers']:
            filein    = os.path.abspath(dir_simu+"/"+RCMs[model]['id']+"/"+var+"_"+domain+"_"+name_simu+"_"+RCMs[model]['inst']+'_'+RCMs[model]['id']+'_'+RCMs[model]['vr']+'_'+freq+"_"+str(tper[0])+"01-"+str(tper[1])+"12.nc")
            print(filein)
            with xr.open_dataset(filein,decode_times=False) as dataset:
                for imm,mm in enumerate(months):
                    RCMs[model]['data'][(tper[0]-y0):(tper[1]-y0+1),imm,:,:] = dataset.variables[var][mm::12, :, :]
        if season=='DJF':
            RCMs[model]['data'][1:nban,2,:,:] = RCMs[model]['data'][0:nban-1,2,:,:] # DJF1980 = D1979+JF1980
            RCMs[model]['data'][0,2,:,:] = np.nan
        for i in range(0,nlat):
            for j in range(0,nlon):
                tmp = np.mean(RCMs[model]['data'][y0P-y0:y1P-y0+1,:,i,j],axis=1)
                RCMs[model]['reglin'][i,j] = 10*np.polyfit(range(y0P,y1P+1),tmp,1)[0]
                #if i==100 and j==50:
                #    print(tmp)
                #    print( RCMs[model]['reglin'][i,j])

        np.save(filesave,RCMs[model]['reglin'],allow_pickle=False)


nbRCM = len(RCMs)

#################################
####   Plot characteristics   ###
#################################

# Output file
#dir_out   = "/home/nabat/ALADIN/EURO-CORDEX/figures/"
dir_out   = "/d0/images/nabat/"
fileout   = dir_out + "plot_EUC12-ERA5"+var+"_"+season+"_"+str(y0P)+"-"+str(y1P)+"_map_trend_multimodel.eps"

if fileout is not None:
    fileout = os.path.abspath(fileout)
    if os.path.isfile(fileout):
        os.remove(fileout)

# Values
min_value    = -0.1
max_value    = 1.0
step         = 0.05
range_values = np.arange(min_value, max_value, step)

# Graphic projection
projplot       = "LambertConformal"
projplot_dict  = dict(central_longitude=10.5, central_latitude=49.5, standard_parallels=(33,45))
projplot       = ccrs.__dict__[projplot].__call__(**projplot_dict)

# Graphic aspects 
colormap     = cm.batlowK  #"CMRmap_r"
#colormap     = "CMRmap_r"


##########################
####   Make the plot   ###
##########################

fig = plt.figure()
fig, axes = plt.subplots(nrows=5,ncols=4,
                        subplot_kw={'projection': projplot})
                        #figsize=(11,8.5))

# axs is a 2 dimensional array of `GeoAxes`.  We will flatten it into a 1-D array
axes=axes.flatten()

# Plot data
for imod,model in enumerate(RCMs):
    print("Plot "+model)
    #if imod == 14:
    projplot       = ccrs.__dict__[RCMs[model]['projplot']].__call__(**RCMs[model]['projplot_dict'])
    map_plot = axes[imod].contourf(RCMs[model]['lon'], RCMs[model]['lat'], RCMs[model]['reglin'], 
                       range_values,   cmap=colormap, extend="both",
                       projection = projplot,
                       transform = projorig)
    axes[imod].set_title(model, fontsize=3.5, pad=1.0)
    axes[imod].add_feature(cf.COASTLINE.with_scale('110m'), linewidth=0.3)
    #axes[imod].add_feature(cf.BORDERS)
    axes[imod].set_extent((-15,38, 25, 70))  # All Europe domain
    #axes[imod].set_extent((-5,13, 38, 54))  # France zoom
    #ax = plt.axes(projection=projplot)
    #map_plot = ax.contourf(RCMs[model]['lon'], RCMs[model]['lat'], RCMs[model]['reglin'], 
    #                   range_values,   cmap=colormap, extend="max",
    #                   transform=projorig)
    #plt.title(model)

# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(left = 0.3, bottom = 0.15, right = 0.7, top = 0.9, wspace = 0.2, hspace = 0.25)

# Add a colorbar axis at the bottom of the graph
cbar_ax = fig.add_axes([0.35, 0.1, 0.3, 0.02]) #[left, bottom, width, height]
cbar = fig.colorbar(map_plot, cax=cbar_ax,orientation='horizontal',
        ticks = [0,0.2,0.4,0.6,0.8],
            shrink=0.85, pad=0.05)
cbar.set_label("Near-surface temperature (K.dec"+r'$^{-1}$'+")", fontsize=7.5)
cbar.ax.tick_params(labelsize=7.5)

# Add a big title at the top
plt.suptitle("Trend for "+var+" "+season+" "+str(y0P)+"-"+str(y1P))


# File output
if fileout is None:
    plt.show()
else:
    plt.savefig(fileout,bbox_inches='tight')

