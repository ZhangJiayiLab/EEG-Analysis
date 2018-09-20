"""
data analysis for project Entrain.

Further plan:
- generalize EntrainAnalysis class for future use.
- documentation
- methodology check

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Sept 20 2018
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat
from tqdm import tqdm
import os

import general
import eegFilter as eegfilter
import stfft
import datacleaner

class EntrainAnalysis(object):
    def __init__(self, datadir, resultdir, filename,
                fs, fftwindow, fftoverlap, roi, entrain_cue_onset):
        datacleaner.processCompact(datadir, filename)
        data = loadmat(os.path.join(datadir, "compact"+filename))
        self.channels = data["channels"]
        self.cue_onset = data["cue_onset"][0,:]
        self.entrain_cue_onset = entrain_cue_onset
        self.times = data["times"][0,:]
        self.name = os.path.splitext(filename)[0]
        self.resultdir = resultdir
        
        self.fs = fs
        self.fftwindow = fftwindow
        self.fftoverlap = fftoverlap
        self.roi = roi
        self.bandname = ["delta", "theta", "alpha", "beta",
            "gamma", "low gamma", "high gamma",
            "lowpass", "highpass"]
        
        general.dircheck(resultdir, os.path.splitext(filename)[0])
        
    def tfdomain_analysis(self, plot_latency=False, rho=2):
        ch_split = general.split_datawithmarker(
                self.channels[1, :], self.cue_onset, self.roi, self.fs)
        pxx, tspec = stfft.stfft(ch_split, self.fftwindow, self.fftoverlap, self.fs)
        power_channels = np.zeros((np.size(self.channels,0), np.size(pxx,1),
                                   np.size(pxx, 2)))
        power_entrain = np.zeros((np.size(self.channels,0), np.size(pxx,1),
                                  np.size(pxx, 2)))
        
        print("processing time-frequency domain analysis")
        for chidx in tqdm(range(np.size(self.channels, 0))):
            ch_split = general.split_datawithmarker(
                    self.channels[chidx, :], self.cue_onset, self.roi, self.fs)
            pxx, tspec = stfft.stfft(ch_split, self.fftwindow, self.fftoverlap, 
                                     self.fs, rho=rho)
            
            ch_split_entrain = general.split_datawithmarker(
                    self.channels[chidx, :], self.entrain_cue_onset, self.roi, self.fs)
            pxx_entrain, tspec_entrain = stfft.stfft(ch_split_entrain, self.fftwindow, self.fftoverlap, 
                                     self.fs, rho=rho)
            
            if plot_latency:
                plt.figure(figsize=(10,4))
                plt.imshow(np.mean(pxx, 0), aspect='auto')
                plt.savefig(os.path.join(self.resultdir, "tfERP",
                    self.name, "channel %d-tf.png"%chidx))
                plt.close()

            power_channels[chidx, :, :] = np.mean(pxx, 0)
            power_entrain[chidx, :, :]  = np.mean(pxx_entrain, 0)
            
        self.tfdata = power_channels
        self.tfdata_entrain = power_entrain
        self.tspec = tspec + self.roi[0]
        self.rho = rho
        
        savemat(os.path.join(self.resultdir, "tfERP",
            self.name, self.name + "tf.mat"),
            {"name":self.name,
             "power_channels": power_channels,
             "entrain_channels": power_entrain,
             "tspec": tspec + self.roi[0],
             "rho": rho,
             })
        
    def latency_analysis(self, cutoffbans="gamma", n=3 , fromfile="", 
                         tfdataname="entrain_channels", savedir="", plot_latency=False):
        if fromfile != "":
            data = loadmat(fromfile)
            self.tfdata = data[tfdataname]
            self.tspec = data["tspec"][0,:]
            self.rho = data["rho"][0,0]

        cutoff = eegfilter.getbandrange(cutoffbans)
        
        datarise = []
        for chidx in range(np.size(self.tfdata,0)):
            tfcurve = np.mean(self.tfdata[chidx, self.rho*cutoff[0]:self.rho*cutoff[1], :], 0)
            tfcurve_mu = np.mean(tfcurve[np.where(self.tspec<0)])
            tfcurve_shift = tfcurve - tfcurve_mu
            tfcurve_sig = np.std(tfcurve[np.where(self.tspec<0)])
            thresh = n*tfcurve_sig

            datarise_item = general.group_consecutive(np.where(tfcurve_shift > thresh)[0], gap=5)
            if len(datarise_item[0]) == 0:
                datarise.append([np.nan])
            else:
                points = []
                for each in datarise_item:
                    x0, y0 = self.tspec[each[0]], tfcurve_shift[each[0]]
                    xm1, ym1 = self.tspec[each[0]-1], tfcurve_shift[each[0]-1]
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
        
        

##### ARCHIVE #####
## import data
def import_entrain_data(datadir, filename):
    data = loadmat(os.path.join(datadir, "compact"+filename))
    channels = data["channels"]
    cue_onset = data["cue_onset"][0,:]
    times = data["times"][0,:]
    return (channels, cue_onset, times)

## preview filted wave in each band
def previewbandpass(filename, resultdir, channels, cue_onset, times, bandname):
    print("rendering bandpass preview")
    for chidx in tqdm(range(np.size(channels, 0))):
        for targetname in bandname:
            ch_bp = eegfilter.butter_bandpass(channels[chidx, :],
                                              targetname, fs, 5)
            ch_bp_split = general.split_datawithmarker(
                    ch_bp, cue_onset, plot_bandpass_roi, fs)

            plt.figure(figsize=(25, np.size(cue_onset)))
            for idx in range(np.size(cue_onset)):
                plt.subplot(np.size(cue_onset), 1, idx+1)
                plt.plot(ch_bp_split[idx, :])
                plt.xticks([])
                plt.xlim([0, len(ch_bp_split[idx, :])])
            plt.savefig(os.path.join(resultdir, "bandpass",
                os.path.splitext(filename)[0],
                "channel %d-%s.png"%(chidx, targetname)))
            plt.close()

## time-frequency domain analysis
def tfdomain(filename, resultdir, channels, cue_onset, times):
    ch_split = general.split_datawithmarker(
                channels[1, :], cue_onset, plot_bandpass_roi, fs)
    pxx, tspec = stfft.stfft(ch_split, stfft.hantaper,
                             fftwindow, fftoverlap, fs)
    power_channels = np.zeros((np.size(channels,0), np.size(pxx,1),
                               np.size(pxx, 2)))

    print("processing time-frequency domain analysis")
    for chidx in tqdm(range(np.size(channels, 0))):
        ch_split = general.split_datawithmarker(
                channels[chidx, :], cue_onset, plot_bandpass_roi, fs)
        pxx, tspec = stfft.stfft(ch_split, stfft.hantaper,
                                 fftwindow, fftoverlap, fs)

        if plot_latency:
            plt.figure(figsize=(10,4))
            plt.imshow(np.mean(pxx, 0), aspect='auto')
            plt.savefig(os.path.join(resultdir, "tfERP",
                os.path.splitext(filename)[0],
                "channel %d-tf.png"%chidx))
            plt.close()

        power_channels[chidx, :, :] = np.mean(pxx, 0)

    ## save data
    savemat(os.path.join(resultdir, "tfERP",
        os.path.splitext(filename)[0],
        os.path.splitext(filename)[0] + "tf.mat"),
        {"title":os.path.splitext(filename)[0],
         "power_channels":power_channels,
         "tspec": tspec,
         "rho":5,
         })


