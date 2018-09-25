import os
import numpy as np
from scipy.io import loadmat, savemat
import matplotlib.pyplot as plt
from tqdm import tqdm

import general
import stfft
import eegFilter

class EEGAnalysis(object):

    def __init__(self, datadir, resultdir, patientname, expname, fs, roi):

        # meta data
        self.name = expname
        self.resultdir = os.path.join(resultdir, patientname)
        self.datadir = os.path.join(datadir, patientname, "EEG")

        general.dircheck(self.resultdir, patientname)

        # analysis parameters
        self.fs = fs
        self.roi = roi

        # analysis constant
        self.bandname = ["delta", "theta", "alpha", "beta",
            "gamma", "low gamma", "high gamma",
            "lowpass", "highpass"]

        # import data
        data = loadmat(os.path.join(self.datadir, "Compact", self.name+".mat"))
        self.channels = data["channels"]
        #self.cue_onset = data["cue_onset"][0,:]
        self.times = data["times"][0,:]
        self.markers = data["markers"]

        # results
        self.tfdata = None
        self.latency = None
        self.rho = None


    def tfdomain_analysis(self, frange, markername, baseroi=(1,2), rho=5, 
                          needsavemat=False, matsuffix="_tf", averaged=False,
                          needpreview=False):
        cue_onset = self.markers[markername][0][0][0,:]
        ch_split = general.split_datawithmarker(
                self.channels[1, :], cue_onset, self.roi, self.fs)
        power_channels = np.zeros((np.size(self.channels,0),
                                   np.size(frange), np.size(ch_split, 1)))

        print("processing time-frequency domain analysis")
        for chidx in tqdm(range(np.size(self.channels, 0))):
            ch_split = general.split_datawithmarker(
                    self.channels[chidx, :], cue_onset, self.roi, self.fs)
            pxx = stfft.dwt_tf(ch_split, self.fs, frange, baseroi, reflection=True)

            power_channels[chidx, :, :] = pxx

        self.tfdata = power_channels
        tspec = np.linspace(self.roi[0], self.roi[1], np.size(self.tfdata, 2))

        if needsavemat:
            savemat(os.path.join(self.resultdir, "tf_domain", self.name + matsuffix + ".mat"),
                {"name": self.name,
                 "method": "morlet wavelet",
                 "power_channels": self.tfdata[:,:,::rho],
                 "times": tspec[::rho],
                 "frange": frange,
                 "rho": rho,
                 })
        if needpreview:
            if not os.path.isdir(os.path.join(self.resultdir, 'preview', 'tf_domain', self.name)):
                os.mkdir(os.path.join(self.resultdir, 'preview', 'tf_domain', self.name))
            
            print("start rendering contourf preview of tf_domain ...")
            for chidx in tqdm(range(np.size(self.tfdata, 0))):
                plt.figure(figsize=(5,3))
                plt.contourf(tspec, frange, self.tfdata[chidx, :, :], 30, cmap=plt.get_cmap("jet"))
                plt.savefig(os.path.join(self.resultdir, 'preview', 'tf_domain', self.name, "ch%03d"%chidx+matsuffix+'.png'), bbox_inches='tight')
                plt.close()


    def fast_latency_detection(self, markername, cutoffband="gamma",
                               n=3, gapsec = 0.1, rho=0.5, baseroi=(1.5,2), validlen=0.1,
                               needpreview=False):
        cutoff = eegFilter.getbandrange(cutoffband)
        frange = np.logspace(np.log10(cutoff[0]), np.log10(cutoff[1]), 20)
        cue_onset = self.markers[markername][0][0][0,:]
        
        if needpreview and not os.path.isdir(os.path.join(self.resultdir, 'preview', 'bandpower', self.name)):
            os.mkdir(os.path.join(self.resultdir, 'preview', 'bandpower', self.name))

        datarise = []
        print("processing fast latency analysis for %s"%markername)
        for chidx in tqdm(range(np.size(self.channels, 0))):
            ch_split = general.split_datawithmarker(
                        self.channels[chidx, :], cue_onset, self.roi, self.fs)
            pxx = stfft.dwt_tf(ch_split, self.fs, frange, baseroi, reflection=True)

            power_channels = np.mean(pxx, 0)
            tfcurve = power_channels
            tspec = np.linspace(self.roi[0], self.roi[1], (self.roi[1]-self.roi[0])*self.fs)

            tfcurve_mu = np.mean(tfcurve[np.where(tspec < 0)])
            tfcurve_shift = tfcurve - tfcurve_mu
            tfcurve_sig = np.std(tfcurve[np.where(tspec < 0)])
            thresh = n*tfcurve_sig
            
            if needpreview:
                Ntrial = np.size(pxx, 0)
                plt.figure(figsize=(10,Ntrial))
                for trialn in range(Ntrial):
                    plt.subplot(Ntrial,1,trialn+1)
                    plt.plot(tspec, pxx[trialn, :])
                    plt.vlines(0, np.min(pxx[trialn, :]), np.max(pxx[trialn, :]), linewidth=2)
                    plt.xlim([tspec[0], tspec[-1]])
                plt.savefig(os.path.join(self.resultdir, 'preview', 'bandpower', self.name, "ch%03d_%s_%s.png"%(chidx, markername, cutoffband)), bbox_inches='tight')
                plt.close()

            datarise_item = general.group_consecutive(np.where(tfcurve_shift > thresh)[0], gap=gapsec*self.fs)
            if len(datarise_item[0]) < validlen*self.fs:
                datarise.append([np.nan])
            else:
                points = []
                for each in datarise_item:
                    x0, y0 = tspec[each[0]], tfcurve_shift[each[0]]
                    xm1, ym1 = tspec[each[0]-1], tfcurve_shift[each[0]-1]
                    points.append(x0 - y0 * (x0-xm1)/(y0-ym1))
                datarise.append(points)

        ## latency
        latency = np.array([])
        for idx, item in enumerate(datarise):
            point = np.nan
            for each in item:
                if each>0:
                    point = each
                    break
            latency = np.append(latency, point)
        self.latency = latency

    def entrain_latency_detection(self, markername="entrain", cutoffband="gamma",
                                  n=3, gapsec=0.1, validlen=0.1, rho=2, iti=5):
        cutoff = eegFilter.getbandrange(cutoffband)
        frange = frange = np.logspace(np.log10(cutoff[0]), np.log10(cutoff[1]), 20)
        cue_onset = self.markers[markername][0][0][0, :]
        tspec = np.linspace(self.roi[0], self.roi[1], (self.roi[1]-self.roi[0])*self.fs)
        baseline_cue_onset = np.array([cue_onset[0]-iti*1, 
                                      cue_onset[0]-iti*2,
                                      cue_onset[0]-iti*3,
                                      cue_onset[0]-iti*4,
                                      cue_onset[0]-iti*5])
        
        datarise = []
        
        print("processing entrain latency analysis for %s"%markername)
        for chidx in tqdm(range(np.size(self.channels, 0))):
            baseline_ch_split = general.split_datawithmarker(
                        self.channels[chidx, :], baseline_cue_onset, self.roi, self.fs)
            baseline_pxx = stfft.dwt_tf(baseline_ch_split, self.fs, frange, reflection=True)
            baseline_mu = np.mean(baseline_pxx[:, np.where(tspec < 0)])
            baseline_sig = np.std(baseline_pxx[:, np.where(tspec < 0)])
            thresh = baseline_mu + n*baseline_sig
            
            ch_split = general.split_datawithmarker(
                        self.channels[chidx, :], cue_onset, self.roi, self.fs)
            pxx = stfft.dwt_tf(ch_split, self.fs, frange, reflection=True)
            ch_power = np.mean(pxx, 0)
            
            datarise_item = general.group_consecutive(np.where(ch_power > thresh)[0], gap=gapsec*self.fs)
            if len(datarise_item[0]) < validlen*self.fs:
                datarise.append([np.nan])
            else:
                points = []
                for each in datarise_item:
                    x0, y0 = tspec[each[0]], ch_power[each[0]]
                    xm1, ym1 = tspec[each[0]-1], ch_power[each[0]-1]
                    points.append(x0 - y0 * (x0-xm1)/(y0-ym1))
                datarise.append(points)
            
        ## latency
#         latency = np.array([])
#         for idx, item in enumerate(datarise):
#             point = np.nan
#             for each in item:
#                 if each>0:
#                     point = each
#                     break
#             latency = np.append(latency, point)
        self.latency = np.array(datarise)

    
    def global_events_detection(self, markername="grating", cutoffband="gamma",
                                  n=3, gapsec=0.1, validlen=0.1, rho=2, iti=5):
        cutoff = eegFilter.getbandrange(cutoffband)
        frange = frange = np.logspace(np.log10(cutoff[0]), np.log10(cutoff[1]), 20)
        cue_onset = self.markers[markername][0][0][0, :]
        tspec = np.linspace(self.roi[0], self.roi[1], (self.roi[1]-self.roi[0])*self.fs)
        
        
        datarise = []
        for chidx in tqdm(range(np.size(self.channels, 0))):
            baseline = general.split_datawithmarker(self.channels[chidx, :], 
                                                cue_onset, self.roi, self.fs)
            baseline_mu = np.mean(baseline[:,np.where(tspec<0)])
            baseline_sig = np.std(baseline[:,np.where(tspec<0)])
            thresh = baseline_mu + n*baseline_sig
        
            pxx = stfft.dwt_tf(self.channels[chidx:chidx, :], self.fs, frange, reflection=False)
            pxx[:500,:] = 0
            pxx[-1000:,:] = 0
            datarise_item = general.group_consecutive(np.where(pxx > thresh)[0], gap=gapsec*self.fs)
            if len(datarise_item[0]) < validlen*self.fs:
                datarise.append([np.nan])
            else:
                points = []
                for each in datarise_item:
                    x0, y0 = tspec[each[0]], tfcurve_shift[each[0]]
                    xm1, ym1 = tspec[each[0]-1], tfcurve_shift[each[0]-1]
                    points.append(x0 - y0 * (x0-xm1)/(y0-ym1))
                datarise.append(points)
                
        return datarise
        
        