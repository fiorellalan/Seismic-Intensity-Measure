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
tool = signal_analysis.mytools()



##### GMPE Input  Parameters ######
# The enduser can use pandas to read Excell or any other tabular format 
# To avoid a further dependency I define a dictionary with station parameters

rake = -85
dip = 50
mag = 6.0
ztor = 0.4
depth = 8.1

"""
station_code	st_nation_code	vs30_m_sec	epi_dist	epi_az	JB_dist	rup_dist	Rx_dist	
AMT	                IT	            670	        8.5	    329.5	1.38	4.62	    5.43	
NRC	                IT	            498	        15.3	133	    2.01	9.7	        12.33
GSA	                IT	            492	        38.8	322.8	30.74	30.74	    -2.49	
"""

Rx = {}
Rrup = {}
Rjb = {}
Epi = {}
Vs30 = {}

Rx['IT.AMT'] = 5.43
Rx['IT.NRC'] = 12.33
Rx['IT.GSA'] = -2.49

Rrup['IT.AMT'] = 4.62
Rrup['IT.NRC'] = 9.7
Rrup['IT.GSA'] = 20.74

Rjb['IT.AMT'] = 1.38
Rjb['IT.NRC'] = 2.0
Rjb['IT.GSA'] = 30.74

Epi['IT.AMT'] = 8.5
Epi['IT.NRC'] = 15.3
Epi['IT.GSA'] = 38.8

Vs30['IT.AMT'] = 670
Vs30['IT.NRC'] = 498
Vs30['IT.GSA'] = 492


filename = os.path.join("Output","IM_Amatrice_Corrected")
fout = open(filename, 'w+b')
Parameters = {}
cut = {}

def calcola_dist (evcoord,coord):
    g = Geod(ellps='WGS84')
    az, baz, dist = g.inv(evcoord['longitude'], evcoord['latitude'], coord['longitude'],coord['latitude'])

    hypo = np.sqrt(evcoord['depth']**2+dist**2)
    hypo = hypo/1000.
    #hypo = dist/1000
    return hypo

dire = [os.path.join("Input","CORR_AMA","*")]
for ll, doss in enumerate(dire):
 print (doss)
 for asc in glob(doss):
    tr = tool.read_orfeus(asc)
    trace_name = ("%s.%s.%s") % (tr.stats.network, tr.stats.station, tr.stats.location)
    component = tr.stats.channel
    Parameters[trace_name] = {}

 for asc in glob(doss):
    tr = tool.read_orfeus(asc)
    trace_name = ("%s.%s.%s") % (tr.stats.network, tr.stats.station, tr.stats.location)
    component = tr.stats.channel
    Parameters[trace_name]['station'] = {}
    Parameters[trace_name][component] = {}
    cname = ("%s.%s") % (tr.stats.network, tr.stats.station)
    print ("TEST", trace_name, component, cname)
    
    figname =os.path.join("FIGURES","FIG_IM_" + tr.stats.network +'_'+ tr.stats.station + '_' + tr.stats.channel + "_Corr.pdf")
    print (figname)


    ## Uncomment if you have coordinates but you don't have distances
    sta_coord = {}
    sta_coord['latitude']= tr.stats.sac.stla
    sta_coord['longitude']= tr.stats.sac.stlo
    
    ev_coord = {}
    ev_coord['latitude']= tr.hdr.evla
    ev_coord['longitude']= tr.hdr.evlo
    ev_coord['depth']= tr.hdr.evdp

    dist = calcola_dist(ev_coord,sta_coord)
     
    #dist = np.sqrt(tr.stats.epi**2 + tr.hdr.evdp**2)
    Parameters[trace_name]['station']['stla'] = tr.stats.sac.stla
    Parameters[trace_name]['station']['stlo'] = tr.stats.sac.stlo
    Parameters[trace_name]['station']['distance'] = dist

    delta = tr.stats.delta
    npts  = tr.stats.npts
    time = np.arange(npts)*delta
    
    Parameters[trace_name][component]['delta'] = delta 

    try:
        cutmin, cutmax = cut[cname]
        time,acc,vel,spo, fminu, fmaxu = tool.baseline(tr.data,time,cutmin=cut[cname][0],cutmax=cut[cutname][1],const=1,base="Non" )  
        print(cutmin, cutmax)
    except KeyError:
        time,acc,vel,spo, fminu, fmaxu = tool.baseline(tr.data,time,const=1,base="Non" )  

    pga = tool.pga()
    pgv = tool.pgv()
    pgd = tool.pgd()




    freq_gmpe = np.logspace(np.log10(0.5),np.log10(33),100)
    output_CA2014, output_sigma1m, output_sigma1p = tool.computegmpe( 'CauzziEtAl2014FixedVs30', freq_gmpe, mag, depth, ztor, dip, rake, Rx[cname], Rrup[cname], Rjb[cname], Epi[cname], Vs30[cname])
    output_KO2020, output_sigma1m, output_sigma1p = tool.computegmpe( 'KothaEtAl2020Site', freq_gmpe, mag, depth, ztor, dip, rake, Rx[cname], Rrup[cname], Rjb[cname], Epi[cname], Vs30[cname])
    output_BD2014, output_sigma1m, output_sigma1p = tool.computegmpe( 'BindiEtAl2014Rjb', freq_gmpe, mag, depth, ztor, dip, rake, Rx[cname], Rrup[cname], Rjb[cname], Epi[cname], Vs30[cname])
    output_BO2014, output_sigma1m, output_sigma1p = tool.computegmpe( 'BooreEtAl2014', freq_gmpe, mag, depth, ztor, dip, rake, Rx[cname], Rrup[cname], Rjb[cname], Epi[cname], Vs30[cname])
    output_CY2014, output_sigma1m, output_sigma1p = tool.computegmpe( 'ChiouYoungs2014', freq_gmpe, mag, depth, ztor, dip, rake, Rx[cname], Rrup[cname], Rjb[cname], Epi[cname], Vs30[cname])
    
    Parameters[trace_name][component]['time'] = time 
    Parameters[trace_name][component]['acceleration'] = acc
    Parameters[trace_name][component]['velocity'] = vel 
    Parameters[trace_name][component]['displacement'] = spo 
    
    freq_welch, trace_welch = welch(acc, fs=(1/delta),nperseg=1024)

    Parameters[trace_name][component]['freq_psd'] = freq_welch
    Parameters[trace_name][component]['psd'] = trace_welch
    
    scav,scavtime, cav,cavtime, bcav, bcavtime = tool.cav(unit='cm/s2')
    ari,aritime,husid = tool.arias(unit='cm/s2')

    stock,X,Y,extent =  tool.stockwell(acc,time)

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

    #ax3b = ax3.twinx()
    #ax3b.set_ylim(-1*ylim_acc, ylim_acc)
    #ppr.plot_wave(ax3b, time[:-1],acc, color='grey')
    ax3.set_ylabel(r'$CAV [gs]$', fontsize=fsz)
    ax3.set_xlabel(r'$Time [s]$', fontsize=fsz)
    ax3.set_ylim(bottom=0)

    ax4 = fig.add_subplot(gs[1,2])
    ax4.set_xscale('log')
    ppr.plot_wave(ax4, freq_welch,trace_welch)
    ax4.axvspan(0.05, fminu, facecolor='grey', alpha=0.5)
    ax4.axvspan(fmaxu,50, facecolor='grey', alpha=0.5)
    ax4.set_ylabel(r'$PSD Amplitude [cm/s]$',fontsize=fsz)
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




    ax50 = fig.add_subplot(gs[3,2])
    ax50.set_xscale('log')
    ax50.set_yscale('log')
    ppr.plot_wave(ax50, freq, psa, yscale='log')
    ax50.axvspan(freq[0], fminu, facecolor='grey', alpha=0.5)
    ax50.axvspan(fmaxu,freq[-1], facecolor='grey', alpha=0.5)
    ax50.set_ylabel(r'$SA [cm/s^2]$', fontsize=fsz)
    ax50.set_xlabel(r'$Frequency [Hz]$', fontsize=fsz)
    ax50.set_ylim(np.min(psa),np.max(np.array([np.max(output_CA2014*980), np.max(psa), np.max(output_KO2020*980), np.max(output_BD2014*980), np.max(output_BO2014*980), np.max(output_CY2014*980)])))
    ax50.set_ylim(bottom=0)

    ax50.plot(freq_gmpe, output_CA2014*980,  color = 'blue', label="Cauzzi 2014")
    ax50.plot(freq_gmpe, output_KO2020*980,  color = 'darkorange', label="Kotha 2020")
    ax50.plot(freq_gmpe, output_BD2014*980,  color = 'red', label="Bindi 2014")
    ax50.plot(freq_gmpe, output_BO2014*980,  color = 'green', label="Boore 2014")
    ax50.plot(freq_gmpe, output_CY2014*980,  color = 'yellow', label="ChiouYoungs 2014")

    ax50.grid(which='both')
    ax50.legend()

    ax6 = fig.add_subplot(gs[0,3])
    ppr.plot_wave(ax6, time[:len(aritime)],aritime)
    ppr.add_im(ax6,"Arias",ari,"m/s",0.55,0.1)
    ppr.add_im(ax6,"Husid",husid,"s",0.55,0.3)
    ax6.set_ylabel(r'$Arias [gs]$', fontsize=fsz)
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
sys.exit()

"""
par  = open('Parameters','rb')
dic = pickle.Unpickler(par).load()
"""
