#!/usr/bin/env python -W ignore::Warning
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:20:51 2019

@author: Maria Lancieri and LAURENDEAU-AUR
"""
import os

from obspy import read, read_inventory, UTCDateTime
from matplotlib.patches import Rectangle
from obspy.signal import invsim
from scipy.signal import welch
from obspy import Stream
from glob import glob


import signal_analysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pylab as p
import plot_par as ppr
import pickle
import sys
from matplotlib.cbook import get_sample_data

from pyproj import Geod 


def calcola_dist (coord):
    evla = 39.76
    evlo = 22.21
    evdp = 8000
    g = Geod(ellps='WGS84')
    az, baz, dist = g.inv(evlo, evla, coord['longitude'],coord['latitude'])

    hypo = np.sqrt(evdp**2+dist**2)
    hypo = hypo/1000.
    return hypo


t = UTCDateTime("2021-03-03")
evname ="20210303_0000071"


tool = signal_analysis.mytools()
st = Stream()
for file in glob(os.path.join("Input","MSEED","DATA","*.mseed")):
    try:
        st += read(file)
    except Exception as m:
        print(m)


filename = os.path.join("Output","Parameters_20210303_0000071")
fout = open(filename, 'w+b')
Parameters = {}

cut = {}

for ll, tr in enumerate(st):

    trace_name = ("%s.%s.%s") % (tr.stats.network, tr.stats.station, tr.stats.location)
    print (trace_name)
    Parameters[trace_name] = {}

for ll, tr in enumerate(st):
    inv = read_inventory(os.path.join("Input","MSEED","inventory.xml"))
    tr.remove_response(inventory=inv,zero_mean='True',output="ACC")  
    trace_name = ("%s.%s.%s") % (tr.stats.network, tr.stats.station, tr.stats.location)
    component = tr.stats.channel
    Parameters[trace_name]['station'] = {}
    Parameters[trace_name][component] = {}

    
    figname =os.path.join("FIGURES", evname, "FIG_IM_" + tr.stats.network +'_'+ tr.stats.station + '_' + tr.stats.channel + ".pdf")
    print (figname)

    sta_coord = {}
    sta_coord = inv.get_coordinates(("%s.%s.%s.%s") % (tr.stats.network, tr.stats.station, tr.stats.location,tr.stats.channel), t)
    dist = calcola_dist(sta_coord )

    Parameters[trace_name]['station']['stla'] = sta_coord['latitude']
    Parameters[trace_name]['station']['stlo'] = sta_coord['longitude']
    Parameters[trace_name]['station']['distance'] = dist
    print (trace_name,dist, sta_coord['latitude'], sta_coord['longitude'])

    delta = tr.stats.delta
    npts  = tr.stats.npts
    time = np.arange(npts)*delta
    
    Parameters[trace_name]['delta'] = delta 
    # const facteur de conversion cm/s2
    print ("ACC")
    time,acc,vel,spo, fminu, fmaxu = tool.baseline(tr.data,time,const=0.01,fmin=0.1,npoles=4,coswin=0.05) 

    pga = tool.pga()
    pgv = tool.pgv()
    pgd = tool.pgd()

    Parameters[trace_name][component]['time'] = time 
    Parameters[trace_name][component]['acceleration'] = acc
    Parameters[trace_name][component]['velocity'] = vel 
    Parameters[trace_name][component]['displacement'] = spo 
    freq_welch, trace_welch = welch(acc, fs=(1/delta),nperseg=512)

    Parameters[trace_name][component]['freq_psd'] = freq_welch
    Parameters[trace_name][component]['psd'] = trace_welch
    
    scav,scavtime, cav,cavtime, bcav, bcavtime = tool.cav(unit='cm/s2')
    ari,aritime,husid = tool.arias(unit='cm/s2')

    # when Mseed are blindly processed is suitable to avoid Stockwell because it is time consuming
    #stock,X,Y,extent =  tool.stockwell(acc,time)

    Parameters[trace_name][component]['cav']   = cav
    Parameters[trace_name][component]['scav']  = scav
    Parameters[trace_name][component]['bcav']  = bcav
    Parameters[trace_name][component]['arias'] = ari
    Parameters[trace_name][component]['husid'] = husid

    freq =  np.logspace(-1,np.log10(50),100)
    psa, psv, freq = tool.response(freq=freq)
    spi  = tool.specint('sa')
    hous = tool.specint('psv')

    Parameters[trace_name][component]['freq_response'] = freq
    Parameters[trace_name][component]['sa']  = psa
    Parameters[trace_name][component]['psv'] = psv
    Parameters[trace_name][component]['spi'] = spi
    Parameters[trace_name][component]['hous'] = hous


    fig = plt.figure(1,figsize=(25,16))
    gs = gridspec.GridSpec(4, 4)
    gs.update(left=0.1, right=0.95, bottom=0.05, top=0.95, hspace=0.3, wspace=0.4)
    fsz =18

    ax0 = fig.add_subplot(gs[0,0:2])
    ppr.plot_wave(ax0, time,acc)
    ppr.add_im(ax0,"PGA",pga,"cm/s^2",0.65,0.8)
    ax0.set_ylabel(r'$Acceleration [cm/s^2]$', fontsize=fsz)
    ax0.set_xlabel(r'$Time [s]$',fontsize=fsz)
    ax0.set_title("%s %2.2f km" % (trace_name,dist),fontsize=20 )
    ylim_acc = 1.2*np.max(np.abs(acc))
    ax0.set_ylim(-1*ylim_acc, ylim_acc)


    """
    ax00 = fig.add_subplot(gs[1,0:2])
    ax00.set_ylabel(r'$Frequency [Hz]$', fontsize=fsz)
    ax00.set_xlabel(r'$Time [s]$',fontsize=fsz)
    cmap = ppr.ccmap(plt.cm.plasma_r)
    ax00.imshow(np.abs(stock), interpolation="nearest", extent=extent, cmap=cmap)
    print ("Stock ext",extent)
    ax00.axis('tight')
    ax00.set_yscale('log')
    ax00.tick_params(axis='both',  labelsize=18)
    ax00.set_xlim(extent[0],extent[1])
    ax00.set_ylim(0.1,extent[3])
    """

    ax1 = fig.add_subplot(gs[2,0:2])
    ppr.plot_wave(ax1, time,vel)
    ppr.add_im(ax1,"PGV",pgv,"cm/s",0.65,0.8)
    ax1.set_ylabel(r'$Velocity [cm/s]$',fontsize=fsz)
    ax1.set_xlabel(r'$Time [s]$',fontsize=fsz)
    ylim_vel = 1.2*np.max(np.abs(vel))
    ax1.set_ylim(-1*ylim_vel, ylim_vel)

    ax2 = fig.add_subplot(gs[3,0:2])
    ppr.plot_wave(ax2, time,spo)
    ppr.add_im(ax2,"PGD",pgd,"cm",0.65,0.8)
    ax2.set_ylabel(r'$Displacement [cm]$', fontsize=fsz)
    ax2.set_xlabel(r'$Time [s]$', fontsize=fsz)
    ylim_spo = 1.2*np.max(np.abs(spo))
    ax2.set_ylim(-1*ylim_spo, ylim_spo)

    ax3 = fig.add_subplot(gs[0,2])

    ppr.plot_wave(ax3, time,bcavtime, color="r")
    ppr.add_im(ax3,"BCAV",bcav,"gs",0.05,0.75)

    ppr.plot_wave(ax3,time[:len(scavtime)],scavtime, color="b")
    ppr.add_im(ax3,"SCAV",scav,"gs",0.05,0.90)

    ppr.plot_wave(ax3, time[:len(cavtime)],cavtime)
    ppr.add_im(ax3,"CAV",cav,"gs",0.05,0.60)

    ax3.set_ylabel(r'$CAV [gs]$', fontsize=fsz)
    ax3.set_xlabel(r'$Time [s]$', fontsize=fsz)
    ax3.set_ylim(bottom=0)

    ax4 = fig.add_subplot(gs[1,2])
    ax4.set_xscale('log')
    ppr.plot_wave(ax4, freq_welch,trace_welch)
    ax4.axvspan(0.05, fminu, facecolor='grey', alpha=0.5)
    ax4.axvspan(fmaxu,50, facecolor='grey', alpha=0.5)
    ax4.set_ylabel(r'$PSD Amplitude [gs]$',fontsize=fsz)
    ax4.set_xlabel(r'$Frequence [Hz]$', fontsize=fsz)
    ax4.set_ylim(bottom=0)
    ax4.set_xlim(0.05,50)



    ax5 = fig.add_subplot(gs[2,2])
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ppr.plot_wave(ax5, freq, psa, yscale='log')
    ax5.axvspan(0.05, fminu, facecolor='grey', alpha=0.5)
    ax5.axvspan(fmaxu,50, facecolor='grey', alpha=0.5)
    ax5.set_ylabel(r'$SA [cm/s^2]$', fontsize=fsz)
    ax5.set_xlabel(r'$Frequency [Hz]$', fontsize=fsz)
    ax5.set_ylim(bottom=0)
    ax5.grid(which='both')
    ax5.set_xlim(0.05,50)

    ax6 = fig.add_subplot(gs[0,3])
    ppr.plot_wave(ax6, time[:len(aritime)],aritime)
    ppr.add_im(ax6,"Arias",ari,"gs",0.55,0.1)
    ppr.add_im(ax6,"Husid",husid,"s",0.55,0.3)
    ax6.set_ylabel(r'$Arias [m/s]$', fontsize=fsz)
    ax6.set_xlabel(r'$Time [s]$', fontsize=fsz)
    ax6.set_ylim(bottom=0)

    ax7 = fig.add_subplot(gs[1,3])
    ax7.set_xscale('log')
    ppr.plot_wave(ax7, 1/freq[::-1], psa[::-1])
    ppr.int_surf(ax7,0.1,2.5, 1/freq[::-1],psa[::-1])
    ppr.add_im(ax7,"SPI",spi,"cm/s",0.55,0.9)
    ax7.set_ylabel(r'$SA [cm/s^2]$', fontsize=fsz)
    ax7.set_xlabel(r'$Period [s]$', fontsize=fsz)
    ax7.set_ylim(bottom=0)

    
    ax8 = fig.add_subplot(gs[2,3])
    ax8.set_xscale('log')
    ppr.plot_wave(ax8, 1/freq[::-1], psv[::-1])
    ppr.int_surf(ax8,0.1,2.5, 1/freq[::-1],psv[::-1])
    ppr.add_im(ax8,"HI",hous,"cm",0.55,0.9)
    ax8.set_ylabel(r'$PSV  [cm/s]$', fontsize=fsz)
    ax8.set_xlabel(r'$Period [s]$', fontsize=fsz)
    ax8.set_ylim(bottom=0)

    

    plt.savefig(figname)
    plt.close()

pickle.dump(Parameters, fout)

"""
par  = open('Parameters','rb')
dic = pickle.Unpickler(par).load()
"""
