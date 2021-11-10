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


dset, grid = add_MITgcm_missing_metrics(ds0)

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

ds1 = ds1.rename({'Y':'YC'})


delY = ds2['sumD']/ds2['cntD']/KnumSmth


#%% calculate time mean
import numpy as np

rng = slice(400, 1439)

dispm = delY[rng].mean('time').load()
grdrm = (grdS/(ds1.grdXY**2.0))[rng].mean('time').load()
grd2Ym= (ds1.grdXY**2.0)[rng].mean('time').load()

tracm = dset.TRAC09[rng].mean('time').load()

# calculate Koc
grdmS = dyn.cal_squared_gradient(tracm, dims=['X', 'Y']).load()
grdSm  = grdS[rng].mean('time').load()

Koc = grdSm / grdmS


#%%
#change unit to km
dset['XC'] = dset['XC'] / 1000
dset['YC'] = dset['YC'] / 1000
grdS['XC'] = grdS['XC'] / 1000
grdS['YC'] = grdS['YC'] / 1000

ds1['XC'] = ds1['XC'] / 1000
ds1['YC'] = ds1['YC'] / 1000
ds2['x_bin'] = ds2['x_bin'] / 1000
ds2['y_bin'] = ds2['y_bin'] / 1000

delY['y_bin'] = delY['y_bin'] / 1000
delY['x_bin'] = delY['x_bin'] / 1000

dispm['x_bin'] = dispm['x_bin'] / 1000
dispm['y_bin'] = dispm['y_bin'] / 1000

grdrm['XC'] = grdrm['XC'] / 1000
grdrm['YC'] = grdrm['YC'] / 1000

grd2Ym['XC'] = grd2Ym['XC'] / 1000
grd2Ym['YC'] = grd2Ym['YC'] / 1000

grdmS['XC'] = grdmS['XC'] / 1000
grdmS['YC'] = grdmS['YC'] / 1000

grdSm['XC'] = grdSm['XC'] / 1000
grdSm['YC'] = grdSm['YC'] / 1000

Koc['XC'] = Koc['XC'] / 1000
Koc['YC'] = Koc['YC'] / 1000


#%% plotting KL Keff Koc and gradients in 6 panels
import proplot as pplt
import matplotlib.pyplot as plt

fig, axes = pplt.subplots(nrows=2, ncols=3, figsize=(13,9), wspace=(0.2))

fontsize = 16

ax = axes[0,0]
m1 = ax.pcolormesh(np.log10(dispm), cmap=plt.cm.jet,
                   vmin=0, vmax=4, levels=21)
ax.set_title('time mean $\\tilde{K}_L$', fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar1 = ax.colorbar(m1, loc='b', width=0.15)
# bar1.set_ticks([1, 1.5, 2, 2.5, 3])

ax = axes[0,1]
m2 = ax.pcolormesh(np.log10(grdrm), cmap=plt.cm.jet,
                   vmin=0, vmax=4, levels=21)
ax.set_title('time-mean $\\tilde{K}_{eff}$', fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar2 = ax.colorbar(m2, loc='b', width=0.15)

ax = axes[0,2]
m1 = ax.pcolormesh(np.log10(Koc[:-2,:-2]), cmap=plt.cm.jet,
                   vmin=0, vmax=4, levels=21)
ax.set_title('$K_{OC}$', fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar1 = ax.colorbar(m1, loc='b', width=0.15)

ax = axes[1,0]
m2 = ax.pcolormesh(np.log10(grdSm[:-2,:-2]), cmap=plt.cm.jet,
                   vmin=-14, vmax=-9, levels=21)
ax.set_title('time mean $|\\nabla_{xy} q|^2$', fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar2 = ax.colorbar(m2, loc='b', width=0.15)

ax = axes[1,1]
m2 = ax.pcolormesh(np.log10(grd2Ym[:-2,:-2]), cmap=plt.cm.jet,
                   vmin=-14, vmax=-9, levels=21)
ax.set_title('time mean $|\\nabla_{XY} q|^2$', fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar2 = ax.colorbar(m2, loc='b', width=0.15)

ax = axes[1,2]
m2 = ax.pcolormesh(np.log10(grdmS[:-2,:-2]), cmap=plt.cm.jet,
                    vmin=-14, vmax=-9, levels=21)
# m2 = ax.pcolormesh(np.log10((grdSm/grd2Ym)[:-2,:-2]), cmap=plt.cm.jet,
#                    vmin=0, vmax=4.2, levels=21)
ax.set_title('$|\\nabla_{xy} \overline{q}|^2$', fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar2 = ax.colorbar(m2, loc='b', width=0.15)


axes.format(abc=True, abcloc='l', abcstyle='(a)',
            xlabel='x-coordinate (km)', ylabel='y-coordinate (km)',
            xticks=[0, 500, 1000, 1500, 2000, 2500, 3000],
            yticks=[0, 500, 1000, 1500, 2000])



