# -*- coding: utf-8 -*-
'''
Created on 2020.11.13

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
'''
#%% read gridded data
import xmitgcm
import numpy as np

path = 'I:/cartRL_advSchemes/Leith1_k0/model/'

dset = xmitgcm.open_mdsdataset(data_dir=path,
                               grid_dir=path,
                               prefix=['Stat'], delta_t=300)

dset = dset.isel({'time':slice(0,1441)})
dset['time'] = np.arange(0, 1441, 1) * 86400

#%% reading contour data
import xarray as xr

increase = False
lt = False

tcount = len(dset.time)
dt = 86400

dsKeff = xr.open_dataset(path + 'diags/Keff0Interp_dec_gt.nc')

dsKeff.coords['time'] = np.linspace(0, dt*(tcount-1), tcount).astype('int')
dsKeff = dsKeff.rename({'new':'Y'})

#%% sample along track
import numpy as np
from GeoApps.GridUtils import add_MITgcm_missing_metrics
from GeoApps.ContourMethods import ContourAnalysis
from utils.IOUtils import read_flt_bcolz

tstr = 300
tend = 300 + 50
rng  = dict(time=slice(tstr*dt, tend*dt))
    
# read_flt
fset = read_flt_bcolz('I:/cartRL_advSchemes/genFLT4/.bcolz',
                cond="(time>="+str(tstr*dt)+") & (time<="+str(tend*dt)+")",
                # cond="(npart>10000000) & (npart<=10001000) & (time>0)",
                fields=['npart', 'time', 'x', 'y'],
                out='xarray').dropna('npart')

xpos = fset.x.loc[rng]
ypos = fset.y.loc[rng]
time = fset.time.loc[rng]

# cal_Yeq
ds = dset.loc[rng]
ds, grid = add_MITgcm_missing_metrics(ds, boundary={'X':'extend', 'Y':'extend'})

mask   = ds['TRAC09'] != 0
tracer = ds['TRAC09'].where(mask).load().drop_vars(['rA','Depth','hFacC',
                                                    'maskC','iter','dxF','dyF'])

cm = ContourAnalysis(ds, tracer,
                     dims={'X':'XC','Y':'YC'},
                     dimEq={'Y':'YC'}, grid=grid,
                     increase=increase, lt=lt)
    
ctr   = dsKeff.TRAC09
Yeq   = dsKeff.Yeq

#%%
xpos, time = xr.broadcast(xpos, time)

trajT = tracer.interp(coords={'XC':xpos, 'YC':ypos, 'time':time}) \
                     .dropna('npart')
trajY = cm.interp_to_coords(trajT,
                            ctr.loc[rng], Yeq.loc[rng],
                            interpDim='Y') \
          .transpose('npart','time') \
          .reset_coords(drop=True)

tracerPlt = tracer
tracerPlt.coords['YC'] = ds['YC'] / 1000.0
tracerPlt.coords['XC'] = ds['XC'] / 1000.0

#%%
disT = np.abs(trajT.diff('time'))
disY = np.abs(trajY.diff('time'))

maxID = disY.mean('time').argmax()

maxT = disT.sel(npart=10185620)
maxY = disY.sel(npart=10185620)

#%% sel IDs where (x0 y0) is close to a particle
tstep = 32

x0 = xpos.sel(npart=10004165)[tstep]
y0 = ypos.sel(npart=10004165)[tstep]

IDs = xpos.where(np.hypot(xpos[:,tstep]-x0.values,
                          ypos[:,tstep]-y0.values)<5500*2,
                 drop=True).npart

#%%
large = [10185620, 10372793]
median = [10004165, 10004128]
test = [ 10026421.,  10405323.,
       ]
trajY.sel(npart=test).plot(ylim=[0, 2200000],hue='npart',
                                               add_legend=False)

#%% plot for selected particles
import proplot as plot
import numpy as np


nparts = [10185620, 10372793]
tsteps = [29376000, 29376000]
xlimstr= [690, 690]
xlimend= [1140, 1140]
ylimstr= [950, 950]
ylimend= [1250, 1250]
vmaxs  = [2.8, 2.8]
vmins  = [1.2, 1.2]
levels = [33, 33]
colors = [(0,0.7,0),(0.9,0.2,0.9)]

trajTplt = trajT
trajYplt = trajY

trajTplt.coords['time'] = time.coords['time'] / 86400
trajYplt.coords['time'] = time.coords['time'] / 86400

def plotPanel(ax, tstep, vmax, vmin, levels):
    m=ax.contourf(tracerPlt.sel(time=tstep), levels=levels,
                vmax=vmax, vmin=vmin, cmap='seismic', extend='both')
    ax.plot(trajSample.x[0]/1000.0, trajSample.y[0]/1000.0,
            color=colors[0], linewidth=3)
    ax.plot(trajSample.x[1]/1000.0, trajSample.y[1]/1000.0,
            color=colors[1], linewidth=3)
    ax.scatter(trajSample.x[0].sel(time=tstep)/1000.0,
               trajSample.y[0].sel(time=tstep)/1000.0, s=180, color=colors[0], zorder=10)
    ax.scatter(trajSample.x[1].sel(time=tstep)/1000.0,
               trajSample.y[1].sel(time=tstep)/1000.0, s=180, color=colors[1], zorder=10)
    ax.set_title('tracer and tracks at t={0:d}'.format(int(tstep/86400)))
    return m

array = [
    [1, 1, 2, 2, 3],
    [4, 4, 5, 5, 6]
]

fig, axes = plot.subplots(array, figsize=(12.3,7),
                          sharex=1, sharey=1, wspace=(0.8, 0.8, 0.8, 4.5))

trajSample = fset.sel(npart=nparts)
values = trajT.sel(npart=nparts)
m = plotPanel(axes[0], tsteps[0]+(-1-1)*dt*2,
              vmax=vmaxs[0], vmin=vmins[0], levels=levels[0])
m = plotPanel(axes[1], tsteps[0]+(0-1)*dt*2,
              vmax=vmaxs[0], vmin=vmins[0], levels=levels[0])
m = plotPanel(axes[3], tsteps[0]+(1-1)*dt*2,
              vmax=vmaxs[0], vmin=vmins[0], levels=levels[0])
m = plotPanel(axes[4], tsteps[0]+(2-1)*dt*2,
              vmax=vmaxs[0], vmin=vmins[0], levels=levels[0])
# trajTplt.sel(npart=nparts).plot(hue='npart', linewidth=2, ax=axes[0, 2])
axes[2].plot(trajTplt.sel(npart=nparts[0]),
                linewidth=2, color=colors[0], label='P2')
axes[2].plot(trajTplt.sel(npart=nparts[1]),
                linewidth=2, color=colors[1], label='P1')
axes[5].plot((trajYplt/1000.0).sel(npart=nparts[0]),
                linewidth=2, color=colors[0], label='P1')
axes[5].plot((trajYplt/1000.0).sel(npart=nparts[1]),
                linewidth=2, color=colors[1], label='P2')
# (trajYplt/1000.0).sel(npart=nparts).plot(hue='npart', linewidth=2, ax=axes[1, 2])

fig.colorbar(m, length=0.55, loc='b', label='')

axes.format(abc='(a)',
    xlim=[xlimstr[0], xlimend[0]], ylim=[ylimstr[0], ylimend[0]],
    xlabel='', ylabel=''
)
axes[2].set_title('along-track tracer')
axes[5].set_title('along-track Y-coords')

axes[1].set_yticklabels('')
axes[4].set_yticklabels('')
axes[2].set_ylim([1, 3])
axes[2].set_xlim([330, trajTplt.time[-1]])
axes[5].set_ylim([0, 2200])
axes[5].set_xlim([330, trajTplt.time[-1]])
axes[2].axvline(336, color=(0.7,0.7,0.7), linestyle='--')
axes[2].axvline(338, color=(0.7,0.7,0.7), linestyle='--')
axes[2].axvline(340, color=(0.7,0.7,0.7), linestyle='--')
axes[2].axvline(342, color=(0.7,0.7,0.7), linestyle='--')
axes[5].axvline(336, color=(0.7,0.7,0.7), linestyle='--')
axes[5].axvline(338, color=(0.7,0.7,0.7), linestyle='--')
axes[5].axvline(340, color=(0.7,0.7,0.7), linestyle='--')
axes[5].axvline(342, color=(0.7,0.7,0.7), linestyle='--')


axes[0].text(828 , 1190, 'P1', fontsize=14)
axes[0].text(920 ,  990, 'P2', fontsize=14)
axes[1].text(890 , 1180, 'P1', fontsize=14)
axes[1].text(960 ,  990, 'P2', fontsize=14)
axes[3].text(920 , 1180, 'P1', fontsize=14)
axes[3].text(1050,  980, 'P2', fontsize=14)
axes[4].text(960 , 1200, 'P1', fontsize=14)
axes[4].text(1100,  980, 'P2', fontsize=14)
axes[2].text(345 ,  2.4, 'P2', fontsize=14)
axes[2].text(345 ,  1.6, 'P1', fontsize=14)
axes[5].text(345 , 1580, 'P1', fontsize=14)
axes[5].text(345 ,  340, 'P2', fontsize=14)



#%% plot for selected particles 2
import proplot as plot
import numpy as np


nparts = [10026421.,  10405323.]
# tsteps = [29376000, 29376000]
tsteps = [29376000-86400*9, 29376000-86400*9]
xlimstr= [1800, 1800]
xlimend= [2250, 2250]
ylimstr= [920, 920]
ylimend= [1220, 1220]
vmaxs  = [2.8, 2.8]
vmins  = [1.2, 1.2]
levels = [33, 33]
colors = [(0,0.7,0),(0.9,0.2,0.9)]

trajTplt = trajT
trajYplt = trajY

trajTplt.coords['time'] = time.coords['time'] / 86400
trajYplt.coords['time'] = time.coords['time'] / 86400

def plotPanel(ax, tstep, vmax, vmin, levels):
    m=ax.contourf(tracerPlt.sel(time=tstep), levels=levels,
                vmax=vmax, vmin=vmin, cmap='seismic', extend='both')
    ax.plot(trajSample.x[0]/1000.0, trajSample.y[0]/1000.0,
            color=colors[0], linewidth=3)
    ax.plot(trajSample.x[1]/1000.0, trajSample.y[1]/1000.0,
            color=colors[1], linewidth=3)
    ax.scatter(trajSample.x[0].sel(time=tstep)/1000.0,
               trajSample.y[0].sel(time=tstep)/1000.0, s=100, color=colors[0], zorder=10)
    ax.scatter(trajSample.x[1].sel(time=tstep)/1000.0,
               trajSample.y[1].sel(time=tstep)/1000.0, s=100, color=colors[1], zorder=10)
    ax.set_title('tracer and tracks at t={0:d}'.format(int(tstep/86400)))
    return m

array = [
    [1, 1, 2, 2, 3],
    [4, 4, 5, 5, 6]
]

fig, axes = plot.subplots(array, figsize=(12.3,7),
                          sharex=1, sharey=1, wspace=(0.8, 0.8, 0.8, 4.5))

trajSample = fset.sel(npart=nparts)
values = trajT.sel(npart=nparts)
m = plotPanel(axes[0], tsteps[0]+(-1-1)*dt*2,
              vmax=vmaxs[0], vmin=vmins[0], levels=levels[0])
m = plotPanel(axes[1], tsteps[0]+(0-1)*dt*2,
              vmax=vmaxs[0], vmin=vmins[0], levels=levels[0])
m = plotPanel(axes[3], tsteps[0]+(1-1)*dt*2,
              vmax=vmaxs[0], vmin=vmins[0], levels=levels[0])
m = plotPanel(axes[4], tsteps[0]+(2-1)*dt*2,
              vmax=vmaxs[0], vmin=vmins[0], levels=levels[0])
# trajTplt.sel(npart=nparts).plot(hue='npart', linewidth=2, ax=axes[0, 2])
axes[2].plot(trajTplt.sel(npart=nparts[0]),
                linewidth=2, color=colors[0], label='P2')
axes[2].plot(trajTplt.sel(npart=nparts[1]),
                linewidth=2, color=colors[1], label='P1')
axes[5].plot((trajYplt/1000.0).sel(npart=nparts[0]),
                linewidth=2, color=colors[0], label='P1')
axes[5].plot((trajYplt/1000.0).sel(npart=nparts[1]),
                linewidth=2, color=colors[1], label='P2')
# (trajYplt/1000.0).sel(npart=nparts).plot(hue='npart', linewidth=2, ax=axes[1, 2])

fig.colorbar(m, length=0.55, loc='b', label='')

axes.format(abc='(a)',
    xlim=[xlimstr[0], xlimend[0]], ylim=[ylimstr[0], ylimend[0]],
    xlabel='', ylabel=''
)
axes[2].set_title('along-track tracer')
axes[5].set_title('along-track Y-coords')

axes[1].set_yticklabels('')
axes[4].set_yticklabels('')
axes[2].set_ylim([1, 3])
axes[2].set_xlim([330-10, trajTplt.time[-1]-10])
axes[5].set_ylim([0, 2200])
axes[5].set_xlim([330-10, trajTplt.time[-1]-10])
axes[2].axvline(327, color=(0.7,0.7,0.7), linestyle='--')
axes[2].axvline(329, color=(0.7,0.7,0.7), linestyle='--')
axes[2].axvline(331, color=(0.7,0.7,0.7), linestyle='--')
axes[2].axvline(333, color=(0.7,0.7,0.7), linestyle='--')
axes[5].axvline(327, color=(0.7,0.7,0.7), linestyle='--')
axes[5].axvline(329, color=(0.7,0.7,0.7), linestyle='--')
axes[5].axvline(331, color=(0.7,0.7,0.7), linestyle='--')
axes[5].axvline(333, color=(0.7,0.7,0.7), linestyle='--')


axes[0].text(1965 , 1000, 'P1', fontsize=14)
axes[0].text(2015 , 1050, 'P2', fontsize=14)
axes[1].text(2005 , 1005, 'P1', fontsize=14)
axes[1].text(2051 , 1067, 'P2', fontsize=14)
axes[3].text(2027 , 1010, 'P1', fontsize=14)
axes[3].text(2047,  1085, 'P2', fontsize=14)
axes[4].text(2060 , 1026, 'P1', fontsize=14)
axes[4].text(2061,  1097, 'P2', fontsize=14)
axes[2].text(323 ,  2.4, 'P2', fontsize=14)
axes[2].text(323 ,  1.6, 'P1', fontsize=14)
axes[5].text(323 , 1530, 'P1', fontsize=14)
axes[5].text(323 ,  340, 'P2', fontsize=14)

