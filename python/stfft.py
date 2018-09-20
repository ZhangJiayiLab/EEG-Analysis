"""
short time fast fourier transformation

basic stfft function, based on scipy.fftpack.fft,
with a taper function.

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Sept 20 2018
"""

from scipy.fftpack import fft
import numpy as np


def hantaper(timepoints):
    """Han Taper function. i.e. shifted half cosine."""
    return 0.5-0.5*np.cos(2*np.pi*np.linspace(0,1,timepoints))


def stfft(data, window, noverlap, fs, taper=hantaper, rho=2):
    """Perform Short Time Fast Fourier Transformation
    
    Syntax: (Pxx, Tspec) = stfft(data, taper, window, noverlap, fs, rho=5)
    
    Keyword arguments:
    data     -- a m*n numpy.ndarray. 
                columns as observations (or channels), 
                rows as raw data.
    window   -- fft window, as time points
    noverlap -- overlap points between two windows
    fs       -- sampling frequency
    taper    -- taper function, take timepoints as input
    rho      -- density of fft readout frequency, relative to sampling frequency
    """
    step = window - noverlap
    start = window // 2
    nstep = (np.size(data,-1) - window) // step

    Tspec = np.linspace(start, step*nstep, nstep) / fs
    Pxx = np.zeros((np.size(data,0), rho*500, nstep))

    taper_list = taper(window)
    for idx in range(nstep):
        temp = data[(slice(None), slice(idx*step,idx*step+window))] * taper_list
        power = np.abs(fft(temp, n=rho*fs)/(1.2*fs)) ** 2
        Pxx[(slice(None), slice(None), idx)] = np.log10(power[(slice(None),slice(0,rho*500))])

    return Pxx, Tspec
