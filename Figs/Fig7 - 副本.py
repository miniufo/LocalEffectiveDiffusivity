# -*- coding: utf-8 -*-
'''
Created on 2020.06.19

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
'''

#%% load data
import xarray as xr
import xmitgcm
import numpy as np
from GeoApps.GridUtils import add_MITgcm_missing_metrics
from GeoApps.DiagnosticMethods import Dynamics


path = 'I:/cartRL_advSchemes/Leith1_k0/model/'


ds0 = xmitgcm.open_mdsdataset(path, grid_dir=path, delta_t=300, prefix=['Stat'])
ds1 = xr.open_dataset(path + 'diags/nKeffXY2.nc',
                      chunks={'XC':560, 'Y':400, 'time':1})
ds2 = xr.open_dataset(path + 'diags/binning/binLag_dec_gt.nc',
                      chunks={'x_bin':280, 'y_bin':200, 'time':5})
ds3 = xr.open_dataset(path + 'diags/Keff0Interp_dec_gt2.nc')

ds0['time'] = np.arange(0, 86400*3600, 86400)
ds1['time'] = np.arange(0, 86400*1441, 86400)
ds2['time'] = np.arange(2*86400, 86400*1440, 86400)
ds3['time'] = np.arange(0, 86400*1441, 86400)

#%%
dset, grid = add_MITgcm_missing_metrics(ds0, boundary={'X':'extend','Y':'extend'})

dyn = Dynamics(dset, grid)

tr = dset.TRAC09.where(dset.TRAC09!=0).astype(np.float64)
grdS = dyn.cal_squared_gradient(tr, dims=['X', 'Y'])

trave   = grid.integrate(tr     , ['X', 'Y']).load()
tr2ave  = grid.integrate(tr**2  , ['X', 'Y']).load()
grdSave = grid.integrate(grdS   , ['X', 'Y']).load()

a = (tr2ave.diff('time') / 2 / 86400).rolling(time=1, center=True,
                                            min_periods=1).mean()
b = (grdSave.rolling(time=1, center=True, min_periods=1).mean())[1:]
Knum = -a/b
KnumSmth = Knum.rolling(time=301, center=True, min_periods=1).mean()

# change unit to km
dset['XC'] = dset['XC'] / 1000
dset['YC'] = dset['YC'] / 1000
grdS['XC'] = grdS['XC'] / 1000
grdS['YC'] = grdS['YC'] / 1000
ds1 = ds1.rename({'Y':'YC'})
ds1['XC'] = ds1['XC'] / 1000
ds1['YC'] = ds1['YC'] / 1000
ds2['x_bin'] = ds2['x_bin'] / 1000
ds2['y_bin'] = ds2['y_bin'] / 1000


delY = ds2['sumD']/ds2['cntD']
delY2 = delY / KnumSmth



#%% plot of different time steps
import numpy as np
import proplot as pplt
from utils.XarrayUtils import coarsen


times = [120, 240, 480, 960]

fontsize = 16

fig, axes = pplt.subplots(nrows=2, ncols=2, figsize=(10,9),
                          sharex=True, sharey=True)

for i, time in enumerate(times):
    KLL = np.log10(delY2[time]).rename({'y_bin':'YC', 'x_bin':'XC'}
                                        ).drop_vars(['time','iter']).load()
    
    KeL = np.log10(grdS[time]/(ds1.grdXY[time])**2.0
                   ).interp_like(KLL).drop_vars(['time', 'iter']).load()
    
    KLL_c = coarsen(KLL[1:-1,1:-1], dims=['YC','XC'], smooth_data=False, ratio=1)
    KeL_c = coarsen(KeL[1:-1,1:-1], dims=['YC','XC'], smooth_data=False, ratio=1)
    
    ax = axes[i]
    ax.scatter(KeL_c.values, KLL_c.values, s=1, c='r')
    
    H, x_edges, y_edges = np.histogram2d(KeL_c.values.ravel(),KLL_c.values.ravel(),
                                         bins=(np.linspace(0,5,21),
                                               np.linspace(0,5,21)),
                                         density=True)
    m1 = ax.contour(x_edges[:-1], y_edges[:-1], H.T, levels=[0.1,0.15,0.2,0.25],
                    colors='k')
    ax.set_ylabel('$log10(\\tilde{K}_L)$', fontsize=fontsize)
    ax.set_xlabel('$log10(\\tilde{K}_{eff})$', fontsize=fontsize)
    ax.set_xticklabels([0,1,2,3,4,5], fontsize=fontsize)
    ax.set_yticklabels([0,1,2,3,4,5], fontsize=fontsize)


axes.format(abc='(a)', xlim=[0, 5], ylim=[0,5])


#%% plot of difference resolutions
import numpy as np
import proplot as pplt
from utils.XarrayUtils import coarsen


resos = [1, 2, 4]
titles = ['original', 'one-half resolution','one-fourth resolution']

time = 940

fontsize = 16

fig, axes = pplt.subplots(nrows=1, ncols=3, figsize=(12,4),
                          sharex=True, sharey=True)

for i, (reso, title) in enumerate(zip(resos, titles)):
    KLL = np.log10(delY2[time]).rename({'y_bin':'YC', 'x_bin':'XC'}
                                        ).drop_vars(['time','iter']).load()
    
    KeL = np.log10(grdS[time]/(ds1.grdXY[time])**2.0
                   ).interp_like(KLL).drop_vars(['time', 'iter']).load()
    
    KLL_c = coarsen(KLL[1:-1,1:-1], dims=['YC','XC'], smooth_data=True, ratio=reso)
    KeL_c = coarsen(KeL[1:-1,1:-1], dims=['YC','XC'], smooth_data=True, ratio=reso)
    
    ax = axes[i]
    ax.scatter(KeL_c.values, KLL_c.values, s=1, c='r')
    
    H, x_edges, y_edges = np.histogram2d(KeL_c.values.ravel(),KLL_c.values.ravel(),
                                         bins=(np.linspace(0,5,21),
                                               np.linspace(0,5,21)),
                                         density=True)
    m1 = ax.contour(x_edges[:-1], y_edges[:-1], H.T, levels=[0.1,0.15,0.2,0.25],
                    colors='k')
    ax.set_ylabel('$log10(\\tilde{K}_L)$', fontsize=fontsize)
    ax.set_xlabel('$log10(\\tilde{K}_{eff})$', fontsize=fontsize)
    ax.set_xticklabels([0,1,2,3,4,5], fontsize=fontsize)
    ax.set_yticklabels([0,1,2,3,4,5], fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)


axes.format(abc='(a)', xlim=[0, 5], ylim=[0,5])


