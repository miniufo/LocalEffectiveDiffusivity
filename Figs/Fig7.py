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


#%% plotting snapshots
import proplot as pplt
import numpy as np


times = [0, 120, 240, 480, 960]
# times = [240]

array = [  # the "picture" (0 == nothing, 1 == subplot A, 2 == subplot B, etc.)
    [1, 6, 11],
    [2, 7, 12],
    [3, 8, 13],
    [4, 9, 14],
    [5,10, 15],
]

fig, axes = pplt.subplots(array, figsize=(12,14),
                          hspace=(0.9,0.9,0.9,0.9), wspace=(1.2, 1.2))


for i in range(len(times)):
    ###### column 1 ######
    v = dset.TRAC09[times[i]].load()
    
    v[:, -2:] = np.nan
    v[-2:, :] = np.nan
    
    ax = axes[i]
    m1 = ax.pcolormesh(v, cmap='jet', levels=np.linspace(1,3,41))
    ax.format(urtitle='day '+str(times[i])+'')
    # if i == 4:
    #     bar1 = ax.colorbar(m1, loc='b', width=0.15)
    #     bar1.set_ticks([1, 1.5, 2, 2.5, 3])
    
    
    ###### column 2 ######
    ax = axes[i+5]
    m2 = ax.pcolormesh(np.log10(delY2[times[i]]/0.6), cmap='jet',
                       levels=np.linspace(1,4.5,36))
    # ax.format(urtitle='(day '+str(times[i])+')')
    # if i == 4:
    #     bar2 = ax.colorbar(m2, loc='b', width=0.15)
    #     bar2.set_ticks([1, 2, 3, 4, 5])
    #     bar2.set_ticklabels(['$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'])
    
    
    ###### column 3 ######
    v = np.log10(grdS[times[i]]/(ds1.grdXY[times[i]])**2.0).load()
    
    v[:, -2:] = np.nan
    v[-2:, :] = np.nan
    
    ax = axes[i+10]
    m4 = ax.pcolormesh(v, cmap='jet', levels=np.linspace(1,4.5,36))
    # ax.format(urtitle='(day '+str(times[i])+')')
    # if i == 4:
    #     bar4 = ax.colorbar(m4, loc='b', width=0.15)
    #     bar4.set_ticks([1, 2, 3, 4, 5])
    #     bar4.set_ticklabels(['$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'])
    
    
    ###### column 4 ######
    # v = np.log10(grdS[times[i]].load() / grdm2Ysmth[times[i]])
    
    # v[:, -2:] = np.nan
    # v[-2:, :] = np.nan
    
    # ax = axes[i*4+3]
    # m3 = ax.pcolormesh(v, cmap=plt.cm.jet,
    #                  vmin=1, vmax=4.5, levels=36)
    # ax.format(urtitle='(day '+str(times[i])+')')
    # if i == 4:
    #     bar3 = ax.colorbar(m3, loc='b', width=0.15)
    #     bar3.set_ticks([1, 2, 3, 4, 5])
    #     bar3.set_ticklabels(['$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'])

bar1 = fig.colorbar(m1, loc='b', width=0.15, cols=1, ticks=[1, 1.5, 2, 2.5, 3])
bar2 = fig.colorbar(m4, ticks=[1,2,3,4,4.5], loc='b', cols=(2, 3),
                   ticklabels=['$10^1$','$10^2$','$10^3$','$10^4$','$10^{4.5}$'])

axes.format(abc='(a)', abcloc='ul',
            collabels=['$q$',
                       '$\\tilde{K}_L$',
                       '$\\tilde{K}_{eff}$',
                       # '$(\\nabla_{xy}q/\\overline{\\nabla_{XY}q})^2$'
                       ],
            xlabel='', ylabel='',
            xticks=[0, 500, 1000, 1500, 2000, 2500],
            yticks=[0, 500, 1000, 1500, 2000])

