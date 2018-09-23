import os
import numpy as np
from scipy.io import loadmat
from tqdm import tqdm

import general
import stfft
import eegFilter

class EEGAnalysis(object):
    
    def __init__(self, datadir, resultdir, filename, fs, roi):
        #TODO: data clean
        
        # meta data
        self.name = os.path.splitext(filename)[0]
        self.resultdir = resultdir
        general.dircheck(resultdir, self.name)
        
        # analysis parameters
        self.fs = fs
        self.roi = roi
        
        # analysis constant
        self.bandname = ["delta", "theta", "alpha", "beta",
            "gamma", "low gamma", "high gamma",
            "lowpass", "highpass"]
        
        # import data
        data = loadmat(os.path.join(datadir, filename))
        self.channels = data["channels"]
        self.cue_onset = data["cue_onset"][0,:]
        self.times = data["times"][0,:]
        self.markers = data["markers"]
        
        # results
        self.tfdata = None
        self.latency = None
        
    
    
    def tfdomain_analysis(self, frange, markername, baseroi=(1,2), rho=5, savemat=False, matsuffix="_tf"):
        cue_onset = self.marker[markername][0][0][0,:]
        ch_split = general.split_datawithmarker(
                self.channels[1, :], cue_onset, self.roi, self.fs)
        power_channels = np.zeros((np.size(self.channels,0), 
                                   np.size(frange), np.size(ch_split, 1)))
        
        print("processing time-frequency domain analysis")
        for chidx in tqdm(range(np.size(self.channels, 0))):
            ch_split = general.split_datawithmarker(
                    self.channels[chidx, :], cue_onset, self.roi, self.fs)
            pxx = stfft.dwt_tf(ch_split, self.fs, frange, baseroi)
            
            power_channels[chidx, :, :] = np.mean(pxx, 0)
            
        self.tfdata = power_channels
        self.rho = rho
        
        if savemat:
            savemat(os.path.join(self.resultdir, "tfERP",
                self.name, self.name + matsuffix + ".mat"),
                {"name": self.name,
                 "method": "morlet wavelet",
                 "power_channels": power_channels,
                 "frange": frange,
                 "rho": rho,
                 })
    
    
    def fast_latency_detection(self, markername, cutoffband="gamma", 
                               n=3, gapsec = 0.1, rho=0.5, baseroi=(1.5,2), validlen=0.1):
        cutoff = eegFilter.getbandrange(cutoffband)
        frange = np.linspace(cutoff[0], cutoff[1], int((cutoff[1]-cutoff[0])*rho))
        cue_onset = self.markers[markername][0][0][0,:]
        
        datarise = []
        print("processing fast latency analysis for %s"%markername)
        for chidx in tqdm(range(np.size(self.channels, 0))):
            ch_split = general.split_datawithmarker(
                        self.channels[chidx, :], cue_onset, self.roi, self.fs)
            pxx = stfft.dwt_tf(ch_split, self.fs, frange, baseroi)
            
            power_channels = np.mean(pxx, 0)
            tfcurve = power_channels[100:-200]  # trade off of edge effect
            tspec = np.linspace(self.roi[0], self.roi[1], (self.roi[1]-self.roi[0])*self.fs)[100:-200]
            
            tfcurve_mu = np.mean(tfcurve[np.where(tspec < 0)])
            tfcurve_shift = tfcurve - tfcurve_mu
            tfcurve_sig = np.std(tfcurve[np.where(tspec < 0)])
            thresh = n*tfcurve_sig
            
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
        