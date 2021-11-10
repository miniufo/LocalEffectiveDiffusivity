# -*- coding: utf-8 -*-
'''
Created on 2020.10.19

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
'''

#%% load data
import xmitgcm


trName = 'TRAC09'
path = 'I:/cartRL_advSchemes/Leith1_k0/model/'

ds = xmitgcm.open_mdsdataset(path, grid_dir=path,
                             delta_t=300, prefix=['Stat'])

#%%
tracer = ds[trName].where(ds[trName]!=0)

tracer['XC'] = tracer['XC'] / 1000
tracer['YC'] = tracer['YC'] / 1000

#%% plotting
import proplot as pplt
import matplotlib.pyplot as plt
import numpy as np


tsteps = [0, 150, 300, 600, 1200, 2000]

fontsize = 16

fig, axes = pplt.subplots(nrows=2, ncols=3, figsize=(12,8.5), sharey=3, sharex=3,
                          wspace=(0.16,0.16), hspace=(0.3))


for t, ax in zip(tsteps, np.array(axes).flatten()):
    m = ax.pcolormesh(tracer[t], cmap=plt.cm.jet, vmin=1, vmax=3, levels=41)

    ax.set_title('tracer at t={0:s} days'.format(str(t)), fontsize=fontsize)
    # ax.set_xlabel('', fontsize=fontsize-1)
    # ax.set_ylabel('', fontsize=fontsize-1)
    ax.set_xticklabels([0, 500, 1000, 1500, 2000, 2500, 3000], fontsize=fontsize-1)
    ax.set_yticklabels([0, 500, 1000, 1500, 2000], fontsize=fontsize-1)

fig.colorbar(m, ticks=0.1, loc='b', minorticks=0.05)

axes.format(abc=True, abcloc='l', abcstyle='(a)',
            xlabel='', ylabel='',
            xticks=[0, 500, 1000, 1500, 2000, 2500, 3000],
            yticks=[0, 500, 1000, 1500, 2000])

