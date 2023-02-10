# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 04:24:25 2023

@author: matth
"""

#%% Imports
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import glob
from windrose import WindroseAxes
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import matplotlib as mpl
import math
from math import radians, degrees, pi
import csv
from scipy.interpolate import griddata
import addcopyfighandler


def PollutionRose(windDirection, windSpeed, extraMeasure, numBins = 100, rotate = True, CScale = 'Log'):
    #%% Making inputs easier
    wd = windDirection
    ws = windSpeed
    pm = extraMeasure
    
    #%%
    
    (H2d, wd_edges, ws_edges) = np.histogram2d(wd, ws, bins = numBins)
    wd_idx = np.digitize(wd, wd_edges)
    ws_idx = np.digitize(ws, ws_edges)
    
    wd_idx[wd_idx==numBins+1] = numBins
    ws_idx[ws_idx==numBins+1] = numBins
    
    Hnew = np.zeros(H2d.shape)
    Hcount = np.zeros(H2d.shape)
    
    for i, v in enumerate(pm):
        wd_val = wd_idx[i]-1
        ws_val = ws_idx[i]-1
        Hnew[wd_val, ws_val] = Hnew[wd_val, ws_val] + v
        Hcount[wd_val, ws_val] += 1
    
    
    Hnew = np.divide(Hnew,Hcount)
    Hnew = np.nan_to_num(Hnew)
    
    
    #%% stackoverflow bit2  pretty good
    if rotate:
        wd = (wd + 90) %360
    
    wd_rad = np.radians(wd)
    
    WD, WS = np.meshgrid(np.linspace(0, 2*np.pi, 36), np.linspace(min(ws), max(ws), numBins ))
    Z = griddata((wd_rad, ws), pm, (WD, WS), method='linear')
    
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    cmap = plt.get_cmap('plasma')
    cmap.set_under('none')
    img = ax.pcolormesh(WD, WS, Z, cmap=cmap)
    if CScale == 'Log':
        img.norm = mpl.colors.LogNorm()
    elif CScale == 'Lin':
        pass
    plt.colorbar(img)
    plt.show()
    
    return ax

if __name__ == '__main__':
    #%% Pulling data
    ## Met Data
    MetPath = r'D:\ACE-ENA Data\Peripheral\IOP1_Met\*.cdf'
    Metds = xr.concat([xr.open_dataset(f) for f in glob.glob(MetPath)], dim='time')

    # CO and n2o Data
    COPath = r'D:\ACE-ENA Data\Peripheral\IOP1_CO\*.nc'
    COds = xr.concat([xr.open_dataset(f) for f in glob.glob(COPath)], dim='time')



    # Metds.wdir_vec_mean.plot()
    wd = Metds.wdir_vec_mean
    wd = np.array(wd)
    ws = Metds.wspd_vec_mean
    ws = np.array(ws)
    # ax = WindroseAxes.from_ax()
    # ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white')
    # ax.set_legend()

    ## CPC Data
    CPCPath = 'D:\\ACE-ENA Data\\Peripheral\\IOP1_CPC\\*.nc'
    CPCds = xr.concat([xr.open_dataset(f) for f in glob.glob(CPCPath)], dim='time')

    cpcTime = CPCds.time.data
    pm10 = CPCds.concentration
    pm10 = np.array(pm10)
    pm10[pm10 < 0] = 0
    pm10_sparse = []
    cpcTime_sparse = []
    for p in np.arange(0,len(wd)-1):
        pm10_sparse.append(pm10[int(np.floor(p*59.77)):int(np.floor((p+1)*59.77))].mean())
        cpcTime_sparse.append(cpcTime[int(np.floor(p*59.77))])

    pm10_sparse = np.array(pm10_sparse)
    cpcTime_sparse = np.array(cpcTime_sparse)


    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # COds.co.plot()
    start = mdates.date2num(datetime(2017, 6, 25, 0, 0, 0))
    end = mdates.date2num(datetime(2017, 7, 20, 0, 0, 0))
    # plot_microscopy_samples(ax, iop1_df, start, end, 0.5, 1)
    
    wd = wd[:-1]
    ws = ws[:-1]
    
    PollutionRose(wd, ws, pm10_sparse)