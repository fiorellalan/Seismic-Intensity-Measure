#!/usr/bin/env python -W ignore::Warning
# signal_analysis.py 
# Class for seismic acceleration analysis traitment 

# Copyright (C) 2010-2022 Maria LANCIERI

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from scipy import signal
from obspy.signal.invsim import cosine_taper
from obspy.signal.filter import bandpass
from stockwell import st

import scipy as sp
import numpy as np
import sys

class Object(object):
    pass

class mytools:
    
    def read_knet(self,filename):
        tr = Object()
        tr.stats = Object()
        tr.stats.sac = Object()
        tr.hdr = Object()
        build_array = False
        amp = []

        for line in open(filename, "r"):
            word = line.rstrip("\n").split()
            if build_array == True:
                for w in word:
                    val = float(w)
                    #transform to gal
                    val *= scale0
                    val /= scale1
                    #val /=980.
                    amp.append(val)

            if word[0] == "Origin":
                datad = word[2]
                year  =  word[2].split("/")[0]
                month =  word[2].split("/")[1]
                day   =  word[2].split("/")[2]
                hour  = word[3].split(":")[0]
                mnt   = word[3].split(":")[1]
                sec   = word[3].split(":")[2]
                orario = word[3]

            if word[0] == "Lat.":
                tr.hdr.evla = float(word[1])
            if word[0] == "Long.":
                tr.hdr.evlo = float(word[1])
            if word[0] == "Depth.":
                tr.hdr.evdp = float(word[2])
            if word[0] == "Mag.":
                mag = float(word[1])
            if word[0] == "Station":
                if word[1] == "Code":
                    tr.stats.stnm = word[2]
                    tr.stats.station = word[2]
                if word[1] == "Lat.":
                    tr.stats.sac.stla = float(word[2])
                if word[1] == "Long.":
                    tr.stats.sac.stlo = float(word[2])
                if word[1] == "Height(m)":
                    tr.stats.stelv = float(word[2])
            """
            if word[0] == "Record":
                tmptimestr = "%s %s" % (word[2], word[3])
                timestr = tmptimestr.replace("/", "-")
                start_time = UTCDateTime(timestr)
                #start_time must be decreased of 15 secs
                #see http://www.k-net.bosai.go.jp/k-net/man/knetform_en.html
                start_time -= 15
                """    
            if word[0] == "Sampling":
                sampl = float(word[2].replace("Hz", ""))
                delta = 1/sampl
                tr.stats.delta = delta
            if word[0] == "Duration":
                duration = float(word[2])
            if word[0] == "Dir.":
                cmpnm = word[1].replace("-","")
                #KiK-net format:
                if cmpnm == "1": cmpnm="NS1"
                if cmpnm == "2": cmpnm="EW1"
                if cmpnm == "3": cmpnm="UD1"
                if cmpnm == "4": cmpnm="NS2"
                if cmpnm == "5": cmpnm="EW2"
                if cmpnm == "6": cmpnm="UD2"
                tr.stats.channel = cmpnm 
            if word[0] == "Scale":
                scale0, scale1 = word[2].split("/")
                scale0 = float(scale0.replace("(gal)", ""))
                scale1 = float(scale1)
            if word[0] == "Max.":
                max_acc = float(word[3])
            if word[0] == "Memo.":
                build_array = True

             
        tr.hdr.ml =  mag
        evname =  "Undef"
        comp = filename.split('.')[1]
        print ("comp", comp)
        if comp == "NS1" or comp == "EW1" or comp == "UD1":
            tr.stats.network = "KiK-Net"
            tr.stats.location = "Borehole"
            site = "2"
        if comp == "NS2" or comp == "EW2" or comp == "UD2":
            tr.stats.network = "KiK-Net"
            tr.stats.location = "Surface"
            site = "1"
        if comp == "NS" or comp == "EW" or comp == "UD":
             tr.stats.network = "K-Net"
             tr.stats.location = " "
             site = "2"
        tr.stats.npts = len(np.array(amp))
        tr.data = amp
        return tr
        




    def read_orfeus(self,filename):
        tr = Object()
        tr.stats = Object()
        tr.stats.sac = Object()
        tr.hdr = Object()
        
        amp = []
        with open(filename) as ff:
            for ii,line in enumerate(ff):
                if ii==0:
                    tr.hdr.evname = line.split()[1]
                if ii==1:
                    tr.hdr.evid = line.split()[1]
                elif ii==2:
                    tr.hdr.evdate = line.split()[1]
                elif ii==3:
                    tr.hdr.evtime = line.split()[1]
                elif ii==4:
                    tr.hdr.evla = line.split()[1]
                elif ii==5:
                    tr.hdr.evlo = line.split()[1]
                elif ii==6:
                    tr.hdr.evdp = np.float(line.split()[1])
                elif ii==8:
                    try:
                        tr.hdr.mw = line.split()[1]
                    except:
                        tr.hdr.mw = -1234
                elif ii==10:
                    try:
                        tr.hdr.ml = line.split()[1]
                    except:
                        tr.hdr.ml = -1234
                elif ii==12:
                    try:
                        tr.hdr.foca = line.split()[1]
                    except:
                        tr.hdr.foca = -1234
                elif ii==13:
                    tr.stats.network = line.split()[1]
                elif ii==14:
                    tr.stats.station = line.split()[1]
                elif ii==15:
                     tr.stats.stnm = line.split()[1]
                elif ii==16:
                    tr.stats.sac.stla = line.split()[1]
                elif ii==17:
                     tr.stats.sac.stlo = line.split()[1]
                elif ii==18:
                    try:
                         tr.stats.sac.stelv = line.split()[1]
                    except:
                         tr.stats.sac.stelv = "-12345"
                elif ii==19:
                    try:
                        tr.stats.location = line.split()[1]
                    except:
                        tr.stats.location = ""
                elif ii==21:
                    try:
                         tr.stats.vs30 =  line.split()[1]
                    except:
                         tr.stats.vs30 = "-12345"
                elif ii==22:
                     tr.stats.sclass =  line.split()[1]
                elif ii==24:
                     tr.stats.epi =  np.float(line.split()[1])
                elif ii==28:
                    tr.stats.delta =  np.float(line.split()[1])
                elif ii==29:
                    tr.stats.npts =  np.float(line.split()[1])
                elif ii==31:
                    tr.stats.channel = line.split()[1]
                if ii>=64:
                    amp.append(np.float(line))
        amp = np.array(amp)
        tr.data = amp
        return tr 
        

    def baseline (self,acc,time,**kwargs):
        """
        acc  : array 
        time : array
        A bandpass Butterworth filter zero phase shift is applied, the kwargs are:
        fmin   : the high-pass cutoff filter frequency, default 0.1
        fmax   : the low-pass  cutoff filter frequency, default Nyquist frequency
        npoles :number of poles, default 4
        
        coswin : cosine windowing lenght in percentage of signal length, default is 0.1
        
        const : normalization value (i.e. 9.8), default is 1 it is applied only to acceleration
    
        cutmin : cut window beg time, default time[0]
        cutmax : cut window beg time, default time[-1]
    
        """
        delta = time[1]-time[0]
        df = 1/delta
        ny = 1/(2*delta)
        fmin = 0.1
        fmax = ny
        npoles = 4
        const = 1
        cutmin = time[0]
        cutmax = time [-1]
        coswin = 0.1
        for key in kwargs:
            if key == "fmin":    fmin   = kwargs[key]
            if key == "fmax":    fmax   = kwargs[key]
            if key == "npoles":  npoles = kwargs[key]
            if key == "coswin":  coswin = kwargs[key]
            if key == "const":   const  = kwargs[key]
            if key == "cutmin":  cutmin = kwargs[key]
            if key == "cutmax":  cutmax = kwargs[key]
            if key == "base":    option = kwargs[key]

        # cut
        cut_beg = np.int(np.ceil(cutmin/delta))
        cut_end = np.int(np.floor(cutmax/delta))
        #print (cutmin,cutmax,cut_beg,cut_end)
        
        acc = acc[cut_beg:cut_end]
        tempo = np.arange(len(acc)) * delta
    
        # detrend
        acc = signal.detrend(acc)
        # Tapering
        nlen = len(tempo)
        win = cosine_taper(nlen,coswin)
        acc = acc * win
        famp = bandpass(acc,fmin,fmax,df,corners=npoles,zerophase=True)
    
        vel = sp.cumsum(famp)*delta
     
        vel = signal.detrend(vel)
        #win = cosine_taper(len(vel),coswin)
        #vel = vel * win
        fvel = vel
        #fvel = bandpass(vel,fmin,fmax,df,corners=npoles,zerophase=True)
        fvel = signal.detrend(fvel)
        
        spo = sp.cumsum(fvel)*delta
        fspo = spo
        #fspo = bandpass(spo,fmin,fmax,df,corners=npoles,zerophase=True)
        
        fspo = signal.detrend(fspo)
        win = cosine_taper(len(fspo),coswin)
        fspo = fspo * win
    
        velyy = np.gradient(fspo,delta)
        accyy = np.gradient(velyy, delta)
        accyy = accyy/const
    
    
        vely = np.cumsum(accyy)*delta
        vely = vely
    
        spoy = np.cumsum(vely)*delta
    
        self.acc = accyy
        self.vel = vely
        self.spo = spoy 
        self.ny = ny
        self.delta = delta
        self.time = tempo 

        # this equation is true if a zerophase true butterwoth is applied once
        gain = -0.2
        # tranformation en decibel
        eta = 10**(gain/20)
        ratio = (eta/(1-eta))**(1/(2*npoles))
       
        fminu = fmin*ratio
        fmaxu = fmax/ratio 
        
        try:
            if option == "Non":
                print ("No baseline")
                acc = acc/const
                self.acc = acc
                self.vel = sp.cumsum(acc)*delta
                self.spo = sp.cumsum(self.vel)*delta
                self.ny = ny
                self.delta = delta
                self.time = tempo 


                accyy = acc
                vely = sp.cumsum(acc)*delta
                spoy = sp.cumsum(self.vel)*delta
                fminu = 0.001
                fmaxu = ny
                 
        except:
            pass


        return tempo, accyy, vely, spoy, fminu, fmaxu

    def baseline_vel (self,vel,time,**kwargs):
        """
        acc  : array 
        time : array
        A bandpass Butterworth filter zero phase shift is applied, the kwargs are:
        fmin   : the high-pass cutoff filter frequency, default 0.1
        fmax   : the low-pass  cutoff filter frequency, default Nyquist frequency
        npoles :number of poles, default 4
        
        coswin : cosine windowing lenght in percentage of signal length, default is 0.1
        
        const : normalization value (i.e. 9.8), default is 1 it is applied only to acceleration
    
        cutmin : cut window beg time, default time[0]
        cutmax : cut window beg time, default time[-1]
    
        """
        delta = time[1]-time[0]
        df = 1/delta
        ny = 1/(2*delta)
        fmin = 0.1
        fmax = ny
        npoles = 4
        const = 1
        cutmin = time[0]
        cutmax = time [-1]
        coswin = 0.1
        for key in kwargs:
            if key == "fmin":    fmin   = kwargs[key]
            if key == "fmax":    fmax   = kwargs[key]
            if key == "npoles":  npoles = kwargs[key]
            if key == "coswin":  coswin = kwargs[key]
            if key == "const":   const  = kwargs[key]
            if key == "cutmin":  cutmin = kwargs[key]
            if key == "cutmax":  cutmax = kwargs[key]
        
        # cut
        cut_beg = np.int(np.ceil(cutmin/delta))
        cut_end = np.int(np.floor(cutmax/delta))
        print (cutmin,cutmax,cut_beg,cut_end)
        
        vel = vel[cut_beg:cut_end]
        vel = vel/const
        tempo = np.arange(len(vel)) * delta
        nlen = len(tempo)
        
        
        vel = signal.detrend(vel)
        win = cosine_taper(nlen,coswin)
        vel *= win
        vel = bandpass(vel,fmin,fmax,df,corners=npoles,zerophase=True)

        acc = np.gradient(vel, delta)
    
        # detrend
        acc = signal.detrend(acc)
        # Tapering
        win = cosine_taper(nlen,coswin)
        acc = acc * win
    
        #famp = bandpass(acc,fmin,fmax,df,corners=npoles,zerophase=True)
        #acc = famp 
        #amp = signal.detrend(tamp_pad)
        #famp = bandpass(tamp_pad,fmin,fmax,df,corners=npoles,zerophase=True)
    
        #beg = np.int(np.ceil(40./delta))
        #end = np.int()
        """
        vel = sp.cumsum(famp)*delta
     
        vel = signal.detrend(vel)
        win = cosine_taper(len(vel),coswin)
        #vel = vel * win
        """
        #fvel = vel
        #fvel = bandpass(vel,fmin,fmax,df,corners=npoles,zerophase=True)
        #fvel = signal.detrend(fvel)
        
        spo = sp.cumsum(vel)*delta
        #fspo = spo
        #fspo = bandpass(spo,fmin,fmax,df,corners=npoles,zerophase=True)
        """
        fspo = signal.detrend(fspo)
        win = cosine_taper(len(fspo),coswin)
        fspo = fspo * win
        """
        #velyy = np.gradient(fspo,delta)
        #accyy = np.gradient(velyy, delta)
        #accyy = accyy/const
    
    
        #vely = np.cumsum(accyy)*delta
        #vely = vely
    
        #spoy = np.cumsum(vely)*delta
    
        self.acc = acc
        self.vel = vel
        self.spo = spo 
        self.ny = ny
        self.delta = delta
        self.time = tempo 

        # this equation is true if a zerophase true butterwoth is applied once
        gain = -0.2
        # tranformation en decibel
        eta = 10**(gain/20)
        ratio = (eta/(1-eta))**(1/(2*npoles))
       
        fminu = fmin*ratio
        fmaxu = fmax/ratio 
        return tempo, acc, vel, spo, fminu, fmaxu
        
    def read_wave(self,acc,time):
        """
        for signals to be just checked
        just to put in memory the acc,time values
        and compute velocity and displacement
        """
        delta = time[1]-time[0]
        df = 1/delta
        ny = 1/(2*delta)
        vel = sp.cumsum(acc)*delta
        spo = sp.cumsum(vel)*delta
        self.acc = acc
        self.vel = vel
        self.spo = spo 
        self.ny = ny
        self.delta = delta
        self.time = time 


        return time, acc, vel, spo,delta

    def pga (self):
        pga = np.max(np.abs(self.acc))
        return pga

    def pgv (self):
        pgv = np.max(np.abs(self.vel))
        return pgv
    
    def pgd (self):
        pgd = np.max(np.abs(self.spo))
        return pgd


    """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       CAV calculation    
       EPRI 1988 -> CAV
       RFS V.II  -> Bracketed CAV
       EPRI 2006 -> Standardizd CAV
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """

    def cav (self, **kwargs):
    
        """
        the kwargs unit :  string g,cm/s2 or m/s2 
        since the baseline normalize the acceleration the default value is g

        the 3 CAV are computed
        scav :  standardized CAV 
        cav  :  classical cav
        bcav : bracketed duration cav (5% of g)
        """    

        unit = "g"
        for key in kwargs:
            if key == "unit":    unit = kwargs[key] 

        if unit == "g"       :  const = 1
        elif unit == "cm/s2" :  const = 980
        elif unit == "m/s2"  : const = 9.8

        
        acc = self.acc/const

        # classic CAV
        cav  = np.trapz(abs(acc))*self.delta
        cavtime = sp.integrate.cumtrapz(abs(acc))*self.delta

        # bracketed duration CAV
        bdur_beg = np.where(acc>= 0.05)[0]
        bcavtime = np.zeros(len(acc))
        if len(bdur_beg) > 2:
            bdur_beg = np.where(acc>= 0.05)[0][0]
            bdur_end = np.where(acc<= 0.05)[0][-1]
            bcav = np.trapz(abs(acc[bdur_beg:bdur_end]))*self.delta
            bcavtime[bdur_beg:bdur_end-1] = sp.integrate.cumtrapz(abs(acc[bdur_beg:bdur_end]))*self.delta
        else:
            bcav = 0

        # standardized CAV 
        t0   = 1.0                                        # time window for check
        nw   = int(np.floor(len(acc)*self.delta/t0));     # number of windows
        ns   = int(np.ceil(t0/self.delta));               # npts in each time window

        scav = 0;
        scavtime  = np.zeros(int(nw*ns));
        for ii in range(1,nw):
            w    = [];
            imin = int((ii-1)*ns+1);
            imax = int(ii*ns+1);
            if imax < len(acc):
                maxw = np.max(np.abs(acc[imin:imax]))
                if maxw >= 0.025:
                    scav =scav+np.trapz(abs(acc[imin:imax]))*self.delta
                    scavtime[imin:imax]=scav
                else:
                    scav = scav
                    scavtime[imin:imax]=scav
        return scav,scavtime, cav,cavtime, bcav, bcavtime

    """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Arias Intensity 
     Defined as: pi/2g int[a(t)**2 dt]
     Defined in m/s the accelartion has to be in m/s2
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    def arias(self, **kwargs):
        """
        the kwargs unit :  string g,cm/s2 or m/s2 
        since the baseline normalize the acceleration the default value is g
        """    
        unit = "g"
        for key in kwargs:
            if key == "unit":    unit = kwargs[key] 

        if   unit == "g"     :  const = (1/9.8)
        elif unit == "cm/s2" :  const = 100
        elif unit == "m/s2"  : const = 1
     
        
        ncoef  = np.pi/(2*9.8)
        acc = self.acc/const * np.sqrt(ncoef)
        aritime = sp.integrate.cumtrapz(acc**2)*self.delta
        arias = np.trapz(acc**2)*self.delta

        hus_to = np.where(aritime >= 0.05 * arias)[0][0]
        hus_tf = np.where(aritime >= 0.95 * arias)[0][0]
        
        husid = self.time[hus_tf]-self.time[hus_to ]
        return arias,aritime,husid

    """
    Response spectra
    Acceleration and pseudo-velocity
    """
    def response (self, **kwargs):
        """
        Compute acceleration and pseudo-velocity spectra
        kwargs 
            damp : integer damping value, default 5
            freq : array frequency, default logspace -1,25,100 points
        """
    
        damp = 5
        freq =  np.logspace(-1,np.log10(25),100)
        for key in kwargs: 
            if key == "damp" : damp = kwargs[key]
            if key == "freq" : freq = kwargs[key]
    
        
        amp = self.acc
        damp = damp/100.
        per = 1/freq
        self.per = per
        npts_amp = len(self.acc)
        npts_out = len(per)
    
        sqd = np.sqrt(1-damp**2)
        dp = per/10.0
        sa = np.zeros(npts_out)
        sv = np.zeros(npts_out)
    
        #compteur=0.0
        for j in range(0,npts_out):
            amax = 0.0
            vmax = 0.0
            elle = 1.0
            if dp[j]<self.delta: 
                elle = int(self.delta/dp[j]+1.0 - 0.00001)
            dtl = self.delta/elle
            omega = 2*np.pi/per[j]
            w2 = omega**2
            w2d = w2*dtl
            wsqd = omega*sqd
            zeta = np.exp(-damp*omega*dtl)
            xt = wsqd*dtl
            sxt = np.sin(xt)
            dsxt = damp*sxt/sqd
            cxt = np.cos(xt)
            a11 = zeta*(dsxt+cxt)
            a12 = zeta*sxt/wsqd
            a21 = -a12*w2
            a22 = zeta*(cxt-dsxt)
            ga1 = amp[0]
            v0 = 0.0
            a0 = -ga1
            for i in range(0,npts_amp-2):
                dx = (amp[i+1]-amp[i])/elle
                dxwd = dx/w2d
                elle=int(elle)
                if elle==1.0:
                    ga1 = amp[i+1]
                    vdxwd = v0+dxwd
                    v1 = a11*vdxwd+a12*a0-dxwd
                    a1 = a21*vdxwd+a22*a0
                    xa1 = ga1+a1
                    sa1 = abs(xa1)
                    sv1 = abs(v1)
                    if sa1>amax: amax=sa1                    
                    if sv1>vmax: vmax=sv1                    
                    v0 = v1
                    a0 = a1
                else:                    
                    for k in range(1,int(elle)+1):
                        ga1 = ga1+dx
                        vdxwd = v0+dxwd
                        v1 = a11*vdxwd+a12*a0-dxwd
                        a1 = a21*vdxwd+a22*a0
                        xa1 = ga1+a1
                        sa1 = abs(xa1)
                        if sa1>amax: amax=sa1                    
                        if sv1>vmax: vmax=sv1
                        v0 = v1
                        a0 = a1
    
            sa[j] = amax
            sv[j] = vmax
            self.sa = sa 
            self.psv  = sv
        return sa, sv, freq
    
        
        """
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         Spectral Intensities
         Housner Instenisty defined as the integral of PSV between 0.1s and 2.5s
         Spectral Intensity defined as the integral of PSA between 0.1 and 0.5s
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
    
    def specint (self,pst):
        """
        :param psds : pseudo spectrum (frequency, amplitude) , type numpy.array 
        :param pst : pseudospectrum type,  type boolean (psa=1 psv=0)     
        """
        if pst=="sa":
            pmin = 0.1
            pmax = 0.5
            rsp = self.sa
        elif pst=="psv":    
            pmin = 0.1
            pmax = 2.5
            rsp = self.psv
           
        # it seems stupid but delta is not uniform because of log spacing 
        somma = 0 
        for i,per in enumerate(self.per) :
            if per <= pmax and per >= pmin:
                delta  = self.per[i]-self.per[i+1]
                somma = somma+(rsp[i]+rsp[i+1])*delta/2
        somma=abs(somma)
        return somma



    """
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 Stockwell Transform
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    def stockwell(self,data,tempo):

        """
        :param data : array of amplitude, type=np.array
        :parm tipe : array of time samples, type=np.array
        """

        delta = self.delta
        fny = 1/(2*delta)
        df = 1/(tempo[-1]-tempo[0])
        nfreq = int(fny/df)
        
        print ("Nyquist", fny)
        while fny > 25:
            data,tempo=sp.signal.resample(data,int(len(data)/2),t=tempo)
            delta =  tempo[1]-tempo[0]
            fny = (1./(2*delta))
            nfreq = int(fny/df)
            df = 1/(tempo[-1]-tempo[0])
            print ("Nyquist", fny)

        fmin = df
        low = int(fmin/df) 
            
        stock = st.st(data,low,nfreq)
        stock = np.flipud(stock)

        time = []
        for i_time in range(len(data)):
            time.append(i_time*delta)

        freq=[]
        for i_freq in range(1,nfreq+1):
            freq.append(i_freq*df)

        X,Y = np.meshgrid(time,freq)
        halfbin_time = delta/2.
        halfbin_freq = df/2.

        extent = (time[0] , time[-1], freq[0], freq[-1])
        return stock,X,Y,extent


    #### GMPE computing ####


    def read_gmpeinput (self):
        from openquake.hazardlib import gsim, imt

        
        
        if self.vs30<180:
            z1pt0 = np.exp(6.745)*0.001
        elif self.vs30<500:
            z1pt0 = np.exp(6.745-1.35*np.log(self.vs30/180))*0.001
        else:
            z1pt0 = np.exp(5.394-4.48*np.log(self.vs30/500))*0.001

        if z1pt0 >0:
            z2pt5 = (519 + 3.595*(z1pt0*1000.))*0.001
        else: 
            z2pt5 = (519 + 3.595*(self.vs30*1000.))*0.001
            
        Rhypo = np.sqrt(self.Epi**2+self.depth**2)
        
        sites = gsim.base.SitesContext()
        rupture = gsim.base.RuptureContext()
        distances = gsim.base.DistancesContext()   
                
        setattr(sites,'vs30', np.array([np.float(self.vs30)]))
        setattr(sites,'z2pt5',np.array([np.float(z2pt5)]))
        setattr(sites,'z1pt0',np.array([np.float(z1pt0)]))
        setattr(sites,'vs30measured', np.array([0]) )
                
        setattr(rupture,'rake',self.rake) 
        setattr(rupture,'dip',np.float(self.dip))
                
        setattr(rupture,'mag',       np.array([np.float(self.mag)]))
        setattr(rupture,'ztor',      np.array([np.float(self.ztor)]))
        setattr(rupture,'hypo_depth',np.array([np.float(self.depth)]))

        setattr(distances,'rhypo',np.array([np.float(Rhypo)]))
        setattr(distances,'rjb',  np.array([np.float(self.Rjb)]))
        setattr(distances,'rrup', np.array([np.float(self.Rrup)]))
        setattr(distances,'rx',   np.array([np.float(self.Rx)]))

        #decomment to check in input
        #print ("Vs30",sites.vs30, "Z2500",sites.z2pt5, "Z1000", sites.z1pt0, "FlagVs30",sites.vs30measured)
        #print ("rake",rupture.rake, "dip", rupture.dip,"mag", rupture.mag, "ztor", rupture.ztor, "depth", rupture.hypo_depth)
        #print ("rhypo", distances.rhypo, "rjb", distances.rjb, "rrup", distances.rrup, "rx", distances.rx)

        return sites,rupture,distances


    def computegmpe(self, gmpe, freq, mag, depth, ztor, dip, rake, Rx, Rrup, Rjb, Epi, vs30):  


        """
        :param gmpe : name of the GMPE on the bases of OQ GSIM list, type=str
        :param freq_gmp : frequencies, type=np.array
        :param mag : event magnitude, type=float
        :param depth : event depth, type=float
        :param ztor : ztor paramaters, type = float
        :param dip : fault dip, type=float
        :param rake : fault rake, type = float
        :param  Rx : Rx parameters, type =float
        :param Rrup : distance from the rupture, type = float 
        :param Rjb : Joyner and Boore distance, type = float
        :param Epi : epicentral distance, type = float
        :param Vs30: velocity of S waves in first 30 meters, type =float

        """


        from openquake.hazardlib import gsim, imt
        self.mag = mag
        self.depth = depth
        self.ztor = ztor
        self.dip = dip
        self.rake = rake
        self.Rx =  Rx
        self.Rrup = Rrup
        self.Rjb = Rjb
        self.Epi = Epi
        self.vs30 = vs30
        AVAILABLE_GSIMS = gsim.get_available_gsims()
        #print (AVAILABLE_GSIMS)
        
        gmpe = AVAILABLE_GSIMS[gmpe]()  
        print (gmpe)
        #print (gmpe.REQUIRES_SITES_PARAMETERS)
        #print (gmpe.REQUIRES_RUPTURE_PARAMETERS)
        #print (gmpe.REQUIRES_DISTANCES)


        sites,rupture,distances = self.read_gmpeinput()  
        periods = np.array(1./freq)        
        periods_entry = ["SA(%s)" % per for per in periods]
        stddev = ['Total']

        output_mean = []
        output_sigma1m = []
        output_sigma1p = []
        output_sigma2m = []
        output_sigma2p = []
        sss = []
        for i_m in periods_entry:
            means,sigma = gmpe.get_mean_and_stddevs(sites,rupture,distances,imt.from_string(i_m),stddev)
            output_mean.append(np.exp(means[0]))
            output_sigma1p.append(np.exp(means[0]+sigma[0][0]))
            output_sigma1m.append(np.exp(means[0]-sigma[0][0]))
            output_sigma2p.append(np.exp(means[0]+2*sigma[0][0]))
            output_sigma2m.append(np.exp(means[0]-2*sigma[0][0]))
            sss.append(sigma[0][0])
        output_mean = np.array(output_mean) 
        output_sigma1p = np.array(output_sigma1p)
        output_sigma1m = np.array(output_sigma1m)
        output_sigma2p = np.array(output_sigma2p)
        output_sigma2m = np.array(output_sigma2m)
        sss = np.array(sss)

        return output_mean, output_sigma1m, output_sigma1p #sss,sites,rupture,distances
