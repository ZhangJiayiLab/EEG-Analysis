from scipy.fftpack import fft
import numpy as np


def hantaper(timepoints):
    return 0.5-0.5*np.cos(2*np.pi*np.linspace(0,1,timepoints))

def stfft(data, taper, window, noverlap, fs):
    step = window - noverlap
    start = window // 2
    nstep = (np.size(data,-1) - window) // step

    Tspec = np.linspace(start, step*nstep, nstep) / fs
    Pxx = np.zeros((np.size(data,0), 1000, nstep))

    taper_list = taper(window)
    for idx in range(nstep):
        temp = data[(slice(None), slice(idx*step,idx*step+window))] * taper_list
        power = np.abs(fft(temp, n=2*fs)/(1.2*fs)) ** 2
        Pxx[(slice(None), slice(None), idx)] = np.log10(power[(slice(None),slice(0,1000))])

    return Pxx, Tspec
