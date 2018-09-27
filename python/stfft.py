"""
time-frequency decomposition analysis using
short time fast fourier transform and
morlet wavelet transform

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Sept 27 2018
"""

from scipy.fftpack import fft, ifft
import numpy as np
# from numba import jit  # numba acceleration


# @jit(nopython=True)
def hantaper(timepoints):
    """Han Taper function. i.e. shifted half cosine."""
    return 0.5-0.5*np.cos(2*np.pi*np.linspace(0,1,timepoints))


# @jit(nopython=True)
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
        Pxx[(slice(None), slice(None), idx)] = np.log10(power[(slice(None), slice(0,rho*500))])

    return Pxx, Tspec


# @jit(nopython=True)
def dwt_tf(eeg_data, fs, frange, baseroi=None, reflection=False, zscore=False, averageall=True):
    """time frequency decomposition analysis with wavelet tranform
    
    Syntax: Pxx = dwt_tf(eeg_data, fs, frange, baseroi)
    
    Keyword arguments:
    data       -- a m*n numpy.ndarray. 
                  columns as observations (or channels), 
                  rows as raw data.
    fs         -- sampling frequency
    frange     -- frequency array for calculation
    baseroi    -- range for normalization baseline
    reflection -- perform data reflection, compensate the edge effect
    zscore     -- using zscore normalization
    averageall -- return averaged total power
    """
    if reflection:
        eeg_data_flip = np.fliplr(eeg_data)
        fft_eeg_data = np.hstack((eeg_data_flip, eeg_data))
        fft_eeg_data = np.hstack((fft_eeg_data, eeg_data_flip))
        pass
    else:
        fft_eeg_data = eeg_data
        
    # wavelet parameters
    wtime = np.linspace(-1,1,2*fs)
    nConv = np.size(fft_eeg_data, 1) + 2*fs
    fft_eeg = np.fft.fft(fft_eeg_data, nConv)
    
    if averageall:
        Pxx = np.zeros((np.size(frange), np.size(fft_eeg_data, 1)))
    else:
        Pxx = np.zeros((np.size(frange), np.size(eeg_data, 0), np.size(fft_eeg_data, 1)))
    
    for idx, F in enumerate(frange):
        s = 6 / (2 * np.pi * F)
        wavelet = np.exp(2*1j*np.pi*wtime*F) * np.exp(- wtime**2/(2*s**2))  # morlet wavelet
        fft_wavelet = np.fft.fft(wavelet, nConv)
        
        conv_wave = np.fft.ifft(fft_wavelet*fft_eeg, nConv)
        conv_wave = conv_wave[:, fs:-fs]
        
        if averageall:
            temppow = np.mean(np.abs(conv_wave)**2,0)
        else:
            Pxx[idx, :, :] = np.abs(conv_wave)**2
            continue
        
        if zscore:
            temppow_cal = 10*np.log10(temppow)  # db-calibration
            temppow_cal = (temppow_cal-np.mean(temppow_cal))/np.std(temppow_cal)
                
        else:
            temppow_cal = 10*np.log10(temppow/np.mean(temppow))  # db-calibration
#         temppow_cal = 10*np.log10(temppow)
        Pxx[idx, :] = temppow_cal

    if reflection:
        if averageall:
            return Pxx[:, np.size(eeg_data,1):-np.size(eeg_data,1)]
        else:
            return Pxx[:,:,np.size(eeg_data,1):-np.size(eeg_data,1)]
    else:
        return Pxx