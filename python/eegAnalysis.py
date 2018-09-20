import numpy as np
from scipy.io import loadmat
from tqdm import tqdm

import general
import stfft

class EEGAnalysis(object):
    
    def __init__(self, datadir, resultdir, filename,
                fs, fftwindow, fftoverlap, roi,
                compact_prefix="compact"):
        #TODO: data clean
        
        # meta data
        self.name = os.path.splitext(filename)[0]
        self.resultdir = resultdir
        general.dircheck(resultdir, self.name)
        
        # analysis parameters
        self.fs = fs
        self.fftwindow = fftwindow
        self.fftoverlap = fftoverlap
        self.roi = roi
        
        # analysis constant
        self.bandname = ["delta", "theta", "alpha", "beta",
            "gamma", "low gamma", "high gamma",
            "lowpass", "highpass"]
        
        # import data
        data = loadmat(os.path.join(datadir, compact_prefix+filename))
        self.channels = data["channels"]
        self.cue_onset = data["cue_onset"][0,:]
        self.times = data["times"][0,:]
        
        # results
        self.tfdata = None
        self.tspec = None
        
    
    def tfdomain_analysis(self, rho=5):
        
        # create result container
        ch_split = general.split_datawithmarker(
                self.channels[1, :], self.cue_onset, self.roi, self.fs)
        pxx, tspec = stfft.stfft(ch_split, self.fftwindow, self.fftoverlap, self.fs)
        power_channels = np.zeros((np.size(self.channels,0), 
                                   np.size(pxx,1), np.size(pxx, 2)))
        
        print("processing time-frequency domain analysis")
        for chidx in tqdm(range(np.size(self.channels, 0))):
            ch_split = general.split_datawithmarker(
                    self.channels[chidx, :], self.cue_onset, self.roi, self.fs)
            pxx, tspec = stfft.stfft(ch_split, self.fftwindow,
                                     self.fftoverlap, self.fs, rho=rho)
            
            power_channels[chidx, :, :] = np.mean(pxx, 0)  # average cross trials
            
        # store data
        self.tfdata = power_channels
        self.tspec = tspec - self.roi[0]
        
        savemat(os.path.join(self.resultdir, "tfERP",
            self.name, self.name + "tf.mat"),
            {"name":self.name,
             "power_channels": power_channels,
             "tspec": tspec - self.roi[0],
             "rho": rho,
             })