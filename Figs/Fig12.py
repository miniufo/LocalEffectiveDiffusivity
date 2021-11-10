# -*- coding: utf-8 -*-
'''
Created on 2020.06.19

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
'''

#%% load data
import xarray as xr
import xmitgcm
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

grdS = dyn.cal_squared_gradient(dset.TRAC09, dims=['X', 'Y'])

# change unit to km
# dset['XC'] = dset['XC'] / 1000
# dset['YC'] = dset['YC'] / 1000
# grdS['XC'] = grdS['XC'] / 1000
# grdS['YC'] = grdS['YC'] / 1000
ds1 = ds1.rename({'Y':'YC'})
# ds1['XC'] = ds1['XC'] / 1000
# ds1['YC'] = ds1['YC'] / 1000
# ds2['x_bin'] = ds2['x_bin'] / 1000
# ds2['y_bin'] = ds2['y_bin'] / 1000


delY = ds2['sumD']/ds2['cntD']

grdm2Y = ((ds3.TRAC09[:,0]-ds3.TRAC09[:,-1])/(ds3.new[0]-ds3.new[-1]))**2
grdm2Ysmth = grdm2Y.rolling(time=31, center=True, min_periods=1).mean()
grdm2Y['time'] = dset['time'][:len(grdm2Y)]
grdm2Ysmth['time'] = dset['time'][:len(grdm2Y)]

# dYsmth = delY.rolling(time=5, center=True, min_periods=1).mean() \
#              .rolling(x_bin=3,center=True, min_periods=1).mean() \
#              .rolling(y_bin=3,center=True, min_periods=1).mean()


#%% calculate time mean
import numpy as np

rng = slice(400, 1439)

dispm = delY[rng].mean('time').load()
grdrm = (grdS/(ds1.grdXY**2.0))[rng].mean('time').load()
grd2m = (grdS/grdm2Y)[rng].mean('time').load()
grd2Ym= (ds1.grdXY**2.0)[rng].mean('time').load()

tracm = dset.TRAC09[rng].mean('time').load()
trvar = ((dset.TRAC09-tracm)**2.0)[rng].mean('time').load()
uvelm = dset.UE_VEL_C[rng].mean('time').load()
vvelm = dset.VN_VEL_C[rng].mean('time').load()
uvar = ((dset.UE_VEL_C-uvelm) ** 2)[rng].mean('time').load()
vvar = ((dset.VN_VEL_C-vvelm) ** 2)[rng].mean('time').load()
cvar = ((dset.UE_VEL_C-uvelm) * (dset.VN_VEL_C-vvelm))[rng].mean('time').load()
EKE  = (((dset.UE_VEL_C-uvelm) ** 2)/2 +
        ((dset.VN_VEL_C-vvelm) ** 2)/2)[rng].mean('time').load()

# calculate Koc
grdmS = dyn.cal_squared_gradient(tracm, dims=['X', 'Y']).load()
grdSm  = grdS[rng].mean('time').load()

Koc = grdSm / grdmS

# calculate EKE ellipse
major = (uvar + vvar + np.sqrt((uvar-vvar)**2.0 + 4.0*cvar**2.0)) / 2.0
minor = uvar + vvar - major
angle = np.arctan2(major - uvar, cvar)

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

tracm['XC'] = tracm['XC'] / 1000
tracm['YC'] = tracm['YC'] / 1000

delY['y_bin'] = delY['y_bin'] / 1000
delY['x_bin'] = delY['x_bin'] / 1000

dispm['x_bin'] = dispm['x_bin'] / 1000
dispm['y_bin'] = dispm['y_bin'] / 1000

grdrm['XC'] = grdrm['XC'] / 1000
grdrm['YC'] = grdrm['YC'] / 1000

grd2Ym['XC'] = grd2Ym['XC'] / 1000
grd2Ym['YC'] = grd2Ym['YC'] / 1000

grd2m['XC'] = grd2m['XC'] / 1000
grd2m['YC'] = grd2m['YC'] / 1000

EKE['XC'] = EKE['XC'] / 1000
EKE['YC'] = EKE['YC'] / 1000

grdmS['XC'] = grdmS['XC'] / 1000
grdmS['YC'] = grdmS['YC'] / 1000

grdSm['XC'] = grdSm['XC'] / 1000
grdSm['YC'] = grdSm['YC'] / 1000

Koc['XC'] = Koc['XC'] / 1000
Koc['YC'] = Koc['YC'] / 1000

major['XC'] = major['XC'] / 1000
major['YC'] = major['YC'] / 1000

minor['XC'] = minor['XC'] / 1000
minor['YC'] = minor['YC'] / 1000

angle['XC'] = angle['XC'] / 1000
angle['YC'] = angle['YC'] / 1000


#%% plotting time mean and EKE
import proplot as pplt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Rectangle


def _make_ellipse(xpos, ypos, major, minor, angle, ax, lw=1, color='k'):
    """Support function for scatter_ellipse."""
    
    for x, y, ma, mj, agl in zip(xpos, ypos, major, minor, angle):
        if not np.isnan(ma):
            ell = Ellipse((x, y), ma, mj, agl,
                          facecolor='none', edgecolor=color,  lw=lw)
            ell.set_clip_box(ax.bbox)
            # ell.set_alpha(0.5)
            ax.add_artist(ell)
            ell.set(clip_box=ax.bbox)


fig, axes = pplt.subplots(nrows=1, ncols=2, figsize=(11.6,5.8), wspace=(0.16))

fontsize = 16

ax = axes[0]
m1 = ax.pcolormesh(tracm.where(tracm!=0), cmap=plt.cm.jet,
                   vmin=1.2, vmax=2.8, levels=33)
ax.set_title('time-mean $q$', fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar1 = ax.colorbar(m1, loc='b', width=0.15)
# bar1.set_ticks([1, 1.5, 2, 2.5, 3])

ax = axes[1]
ypos, xpos = xr.broadcast(major.YC, major.XC)
skip = 8
mag = np.sqrt(major + minor)

xposL = xpos.where(mag>=0.15)
yposL = ypos.where(mag>=0.15)
majorL = major.where(mag>=0.15)
minorL = minor.where(mag>=0.15)
angleL = angle.where(mag>=0.15)

xposS = xpos.where(np.logical_and(mag>=0.08, mag<0.15))
yposS = ypos.where(np.logical_and(mag>=0.08, mag<0.15))
majorS = major.where(np.logical_and(mag>=0.08, mag<0.15))
minorS = minor.where(np.logical_and(mag>=0.08, mag<0.15))
angleS = angle.where(np.logical_and(mag>=0.08, mag<0.15))

m1 = ax.pcolormesh(np.log10(EKE[:-2,:-2]), cmap=plt.cm.jet,
                    vmin=-4.5, vmax=-1, levels=36, zorder=-2)
_make_ellipse(xposL.values[::skip, ::skip].ravel(),
              yposL.values[::skip, ::skip].ravel(),
              majorL.values[::skip, ::skip].ravel()*1200,
              minorL.values[::skip, ::skip].ravel()*1200,
              angleL.values[::skip, ::skip].ravel(),
              ax, lw=0.8, color='k')
_make_ellipse(xposS.values[::skip, ::skip].ravel(),
              yposS.values[::skip, ::skip].ravel(),
              majorS.values[::skip, ::skip].ravel()*6000,
              minorS.values[::skip, ::skip].ravel()*6000,
              angleS.values[::skip, ::skip].ravel(),
              ax, lw=0.8, color='gray')

# make ellipse legends
ox, oy = 2560, 1950
ax.add_patch(Rectangle((ox, oy), 460, 220,
             ec = 'k', fc = 'w', fill=True, lw=0.6))
_make_ellipse([ox+100], [oy+60], [10e-2*1200], [5e-2*1200], [0],
              ax, lw=0.9, color='k')
_make_ellipse([ox+100], [oy+160], [2e-2*6000], [1e-2*6000], [0],
              ax, lw=0.9, color='gray')
ax.text(ox+200, oy+90, '10, 5', fontsize=11, verticalalignment='top')
ax.text(ox+250, oy+190, '2, 1', fontsize=11, verticalalignment='top')

ax.set_title('time-mean EKE & variance ellipses', fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar1 = ax.colorbar(m1, loc='b', width=0.15)


axes.format(abc=True, abcloc='l', abcstyle='(a)', ylabel='y-coordinate (km)',
            xticks=[0, 500, 1000, 1500, 2000, 2500, 3000],
            yticks=[0, 500, 1000, 1500, 2000])



#%% cal local Keff for time-mean q
import numpy  as np
import xmitgcm
import xarray as xr
from GeoApps.DiagnosticMethods import Dynamics
from GeoApps.ContourMethods import ContourAnalysis
from GeoApps.GridUtils import add_MITgcm_missing_metrics


trName = 'TRAC09'
increase = False
lt = False
rng = slice(400, 1439)

print('{0:s}, increase={1:s}, lt={2:s}'.format(trName, str(increase), str(lt)))

path = 'I:/cartRL_advSchemes/Leith1_k0/model/'


dset = xmitgcm.open_mdsdataset(path, grid_dir=path, delta_t=300, prefix=['Stat'])
dset['time'] = np.arange(3600)

dset, grid = add_MITgcm_missing_metrics(dset,
                                        boundary={'Y':'extend','X':'extend'})

tracer = (dset[trName].where(dset[trName]!=0))[rng].mean('time').load()

grdS = Dynamics(dset, grid).cal_squared_gradient(tracer, dims=['X','Y'])

# test linear contour levels
cm = ContourAnalysis(dset, tracer,
                     dims={'X':'XC','Y':'YC'},
                     dimEq={'Y':'YC'}, grid=grid,
                     increase=increase, lt=lt)

mask = tracer.copy()
mask += 1 - mask

#%%
N = 201

width = dset.XG[-1].values

preY = np.linspace(0, 5500*399, 400)

table = cm.cal_area_eqCoord_table(mask)

# ctr     = cm.cal_contours(N)
ctr     = cm.cal_contours_at(preY, table)
area    = cm.cal_integral_within_contours(ctr, out_name='intArea')
intgrdS = cm.cal_integral_within_contours(ctr, integrand=grdS, out_name='intgrdS')

Yeq   = table.lookup_coordinates(area).rename('Yeq')
Lmin  = (area / Yeq).rename('Lmin')

# Yeq2  = (area / width).rename('Yeq').load()
# Lmin2 = cm.cal_minimum_possible_length(lambda x:x-x+width, Yeq2).load()


dintgrdSdA = cm.cal_gradient_wrt_area(intgrdS, area)
dqdA    = cm.cal_gradient_wrt_area(ctr, area)
Leq2    = cm.cal_sqared_equivalent_length(dintgrdSdA, dqdA)
nkeff   = cm.cal_normalized_Keff(Leq2, Lmin)

vs = [ctr, area, intgrdS, Yeq, dintgrdSdA, dqdA, Leq2, Lmin, nkeff]
# vs = [ctr, nkeff]

origin1 = xr.merge(vs)
interp1 = cm.interp_to_dataset(preY, Yeq, origin1).rename({'new':'YC'})

#%% cal local Keff
# mapping back to geographic space
from GeoApps.ArrayUtils import interp1d

def interp_coords(trajThe, TRAC01):
    interpDim = 'YC'
    coord = TRAC01[interpDim]
    
    trajThe = trajThe.rename({'YC':'tmp'})
    
    increasing = False
    
    coord, TRAC01 = xr.broadcast(coord, TRAC01)
    
    varIntp = xr.apply_ufunc(interp1d, trajThe, TRAC01, coord,
              kwargs={'inc': increasing},
              dask='allowed',
              input_core_dims =[[], [interpDim], [interpDim]],
              # output_core_dims=[[interpDim]],
              exclude_dims=set((interpDim,)),
              vectorize=True
              ).rename('yEq')

    return varIntp.rename({'tmp':'YC'})

tr_xy = tracer
tr_XY = interp1[trName]
ke_XY = interp1.nkeff
gd_XY = (interp1.dTRAC09dA*interp1.Lmin)
ke_XY[np.isnan(ke_XY)] = 1
gd_XY[np.isnan(gd_XY)] = 0

ypos  = interp_coords(tr_xy, tr_XY)
mask = np.isnan(ypos)
ypos.values[mask.values] = 1

ke_xy = ke_XY.interp(coords={'YC':ypos})
ke_xy.values[mask.values] = np.nan

gd_xy = gd_XY.interp(coords={'YC':ypos})
gd_xy.values[mask.values] = np.nan

#%% load in time mean local Keff
ds1 = xr.open_dataset(path + 'diags/nKeffXY2.nc',
                      chunks={'XC':560, 'Y':400, 'time':1})

ds1['time'] = ds1['time'] / 86400

grdrm = ds1.nkeff.rename({'Y':'YC'})[rng].mean('time').load()

#%%
tracer['XC'] = tracer['XC'] / 1000
tracer['YC'] = tracer['YC'] / 1000
gd_xy['XC'] = gd_xy['XC'] / 1000
gd_xy['YC'] = gd_xy['YC'] / 1000
grdS['XC'] = grdS['XC'] / 1000
grdS['YC'] = grdS['YC'] / 1000
grdrm['XC'] = grdrm['XC'] / 1000
grdrm['YC'] = grdrm['YC'] / 1000

#%% plotting time mean and local Keff
import proplot as pplt
import matplotlib.pyplot as plt


fig, axes = pplt.subplots(nrows=1, ncols=3, figsize=(12.6,4.6), wspace=(1, 1))

fontsize = 16

ax = axes[0]
m1 = ax.pcolormesh(tracer, cmap='jet',
                   levels=np.linspace(1.2, 2.8, 33))
ax.set_title('time-mean $q$', fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar1 = ax.colorbar(m1, loc='b', width=0.15, label='')
bar1.set_ticks([1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6])

Keff_qm = grdS/gd_xy**2

ax = axes[1]
m1 = ax.pcolormesh(np.log10(Keff_qm), cmap='jet',
                   levels=np.linspace(0, 4, 21))
ax.set_title('$\\tilde{K}_{eff}^{\\overline{q}}$ for $\\overline{q}$', fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar1 = ax.colorbar(m1, loc='b', width=0.15)
bar1.set_ticks([0, 0.8, 1.6, 2.4, 3.2, 4])

ax = axes[2]
m1 = ax.pcolormesh(np.log10(grdrm/Keff_qm), cmap='jet',
                   levels=np.linspace(0, 4, 21))
ax.set_title('$\\overline{\\tilde{K}_{eff}}$ / $\\tilde{K}_{eff}^{\\overline{q}}\\approx K_{OC}$',
             fontsize=fontsize)
ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
bar1 = ax.colorbar(m1, loc='b', width=0.15)
bar1.set_ticks([0, 0.8, 1.6, 2.4, 3.2, 4])


axes.format(abc='(a)', ylabel='y-coordinate (km)',
            xticks=[0, 500, 1000, 1500, 2000, 2500, 3000],
            yticks=[0, 500, 1000, 1500, 2000])



#%% plot time mean q* and (time-mean q)*
import proplot as pplt
import matplotlib.pyplot as plt


ds3 = xr.open_dataset(path + 'diags/Keff0Interp_dec_gt2.nc')

fig, axes = pplt.subplots(nrows=1, ncols=2, figsize=(11,6),
                          sharex=0, sharey=2)

fontsize = 16

rng = slice(400, 1439)

a1 = ds3.TRAC09[rng].mean('time').rename({'new':'YC'})
a2 = interp1.TRAC09
a1['YC'] = a1['YC'] / 1000
a2['YC'] = a2['YC'] / 1000

ax = axes[0]
m1 = ax.plot(a1, a1.YC, linewidth=3)
m2 = ax.plot(a2, a2.YC, linewidth=3)
ax.set_title('$\\overline{q*(Y)}$ v.s. $[\\overline{q}(x,y)]*$', fontsize=fontsize)
ax.set_xlabel('tracer', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
ax.legend([m1,m2], labels=['$\\overline{q*(Y)}$', '$[\\overline{q}(x,y)]*$'],
          fontsize=fontsize-1)


ax = axes[1]
m3 = ax.plot((a2/a1)**2, a2.YC, linewidth=3)
ax.set_title('squared $[\\overline{q}(x,y)]*$ / $\\overline{q*(Y)}$', fontsize=fontsize)
ax.set_xlabel('squared ratio', fontsize=fontsize-2)
ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)

# ax = axes[2]
# m1 = ax.pcolormesh(np.log10(grdrm/Keff_qm), cmap='jet',
#                    levels=np.linspace(0, 4, 21))
# ax.set_title('$\\overline{\\tilde{K}_{eff}}$ / $\\tilde{K}_{eff}^{\\overline{q}}\\approx K_{OC}$',
#              fontsize=fontsize)
# ax.set_xlabel('x-coordinate (km)', fontsize=fontsize-2)
# ax.set_ylabel('y-coordinate (km)', fontsize=fontsize-2)
# bar1 = ax.colorbar(m1, loc='b', width=0.15)
# bar1.set_ticks([0, 0.8, 1.6, 2.4, 3.2, 4])


axes.format(abc='(a)', ylabel='y-coordinate (km)',
            yticks=[0, 500, 1000, 1500, 2000], ylim=[0, 2200])

