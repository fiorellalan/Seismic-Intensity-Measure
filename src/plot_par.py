#!/usr/bin/env python -W warning
# plot library 
# graphical represantation of results obteined with comp_par.py
# (c) 2019 Maria Lancieri <maria.lancieri@irsn.fr>


from __future__ import division
from math import *
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.transforms as transforms
import matplotlib.ticker as mtick
from matplotlib import rc
from matplotlib.colors import LogNorm



def int_surf(ax,lmin,lmax,per,rsp):
    per_shdw = per[(per>=lmin) * (per<=lmax)]
    rsp_shdw = rsp[(per>=lmin) * (per<=lmax)]
    per_shdw = [lmin] + list(per_shdw) 
    per_shdw = per_shdw + [lmax]

    rsp_shdw = [0] + list(rsp_shdw)
    rsp_shdw = rsp_shdw + [0]
    ax.fill(per_shdw,rsp_shdw,facecolor='blue',alpha=0.5)

    
def plot_wave (ax,xvar,yvar, **kwargs):
    col = 'k'
    yscl = 'lin'
    for key in kwargs:
        if key == "color":  col = kwargs[key]
        if key == "yscale":  yscl = kwargs[key] 
    ax.plot(xvar,yvar,color = col, lw = 1)
    ax.set_xlim(xvar[0],xvar[-1])
    ax.tick_params(axis='both',  labelsize=18)

    if yscl == 'lin':
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%1.2f'))
    
        if max(abs(yvar)) <= 0.02 and max(abs(yvar)) != 0 :
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%1.1e'))
    
    return ax


def add_im (ax, imlabel, imvalue, imunit, xtesto,ytesto):
    if imvalue == 0:
        txtim  = " %s  0 $%5s$"  % (imlabel,imunit)
    else:
        txtim  = " %s  %3.2e $%5s$"  % (imlabel,imvalue,imunit)
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transAxes)
    txt = ax.text(xtesto,ytesto,txtim,transform=trans)
    txt.set_bbox(dict(boxstyle="round,pad=0.5",edgecolor='k',fill=False,alpha=1))
    txt.set_fontsize(15)
    return ax

def ccmap (colori): 
    cmap = colori
    cmaplist = [cmap(i) for i in range(cmap.N)] #Get all the colors of "hot" colorbar
    cmaplist[0] = (1,1,1,1.0)   #Fix the first level to white
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom', cmaplist, 131)   #Customize our cmap

    return (cmap)


            
