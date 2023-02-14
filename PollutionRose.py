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
import time
from tqdm import tqdm


def PollutionRose(windDirection, windSpeed, pollutionVec,
                  wslim:float or None = [0, 10],
                  vlim:float or None = [None, None],
                  pThresh:float or None = 0,
                  numBins = 100,
                  interpolation_method:str = 'nearest',
                  rotate:bool = True,
                  CScale:str = 'Log'):
    #%% Making inputs easier
    wd = windDirection
    ws = windSpeed
    pm = pollutionVec
    
    #%% Input Checking
    if rotate:
        wd = (wd + 90) %360
    
    if isinstance(numBins,int):
        numBins_wd = numBins
        numBins_ws = numBins
        
    elif len(numBins) == 1:
        numBins_wd = numBins
        numBins_ws = numBins
        
    elif len(numBins) == 2:
        numBins_wd = numBins[0]
        numBins_ws = numBins[1]
        
    vmin = vlim[0]
    vmax = vlim[1]
    
    if wslim[0]== None:
        wsmin = min(ws)
        
    if wslim[1]== None:
        wsmax = max(ws)
        
    wsmin = wslim[0]
    wsmax = wslim[1]
    
    #%%
    
    (H2d, wd_edges, ws_edges) = np.histogram2d(wd, ws, bins = numBins)
    WD, WS = np.meshgrid(np.linspace(0, 2*np.pi, numBins_wd), np.linspace(min(ws), max(ws), numBins_ws ))
    wd_rad = np.radians(wd)
    wd_idx = np.digitize(wd_rad, WD[0,:])
    ws_idx = np.digitize(ws, WS[:,0])
    
    wd_idx -= 1
    ws_idx -= 1
    # wd_idx[wd_idx==numBins] = numBins-1
    # ws_idx[ws_idx==numBins] = numBins-1
    
    Hnew = np.zeros(H2d.shape)
    Hcount = np.zeros(H2d.shape)
    
    
    Hnew[ws_idx, wd_idx] = Hnew[ws_idx, wd_idx] + pm
    
    if pThresh != None:
        Hnew[Hnew <= pThresh] = np.nan;
    
    # wd_edges_rad = np.radians(wd_edges)
    # wdmesh, wsmesh = np.meshgrid(wd_edges_rad, ws_edges)
    
    fig, ax = plt.subplots(subplot_kw={'projection':'polar'})
    cmap = plt.get_cmap('plasma')
    cmap.set_under('none')
    img = ax.pcolormesh(WD, WS, Hnew, cmap = cmap)
    plt.ylim([wsmin, wsmax])
    if CScale == 'Log':
        img.norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax)
    elif CScale == 'Lin':
        pass
    plt.colorbar(img)
    plt.show()
    
    return ax
    
    #%% stackoverflow bit2  pretty good
    
    
    # wd_rad = np.radians(wd)
    
    # WD, WS = np.meshgrid(np.linspace(0, 2*np.pi, numBins_wd), np.linspace(min(ws), max(ws), numBins_ws ))
    # Z = griddata((wd_rad, ws), pm, (WD, WS), method=interpolation_method)
    # Z[Z==0] = np.nan
    
    # fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    # cmap = plt.get_cmap('plasma')
    # cmap.set_under('none')
    # img = ax.pcolormesh(WD, WS, Z, cmap=cmap)
    # plt.ylim([wsmin, wsmax])
    # if CScale == 'Log':
    #     img.norm = mpl.colors.LogNorm(vmin = vmin, vmax = vmax)
    # elif CScale == 'Lin':
    #     pass
    # plt.colorbar(img)
    # plt.show()
    
    # return ax

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
    
    wd_test = np.array([0, 90, 180, 350])
    ws_test = np.array([2, 6, 8, 10])
    pm102 = np.array([5,20,40,50])
    
    PollutionRose(wd, ws, pm10_sparse, rotate = True, CScale = 'Log')