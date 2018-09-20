"""
EEG bandpass filters
based on scipy.signal

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Sept 20 2018
"""

from scipy.signal import butter, lfilter, filtfilt
import numpy as np


def getbandrange(bandname):
    """Return the cutoff frequecy of each band."""
    bandrange = {
            "delta":(  2,  4),
            "theta":(  4,  8),
            "alpha":(  8, 12),
            "beta": ( 15, 30),
            "gamma":( 30,150),
            "low gamma": (30, 80),
            "high gamma":(80,150),
            "lowpass":(1,200),
            "highpass":(200,500)}
    return bandrange[bandname]


def butter_bandpass_filter(lowcut, highcut, fs, order):
    """Butterworth bandpass filter designer"""
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype="band")
    return b,a


def butter_bandpass(data, bandname, fs, order=5):
    """Butterworth bandpass filter
    
    Syntax: y = butter_bandpass(data, bandname, fs, order=5)
    
    Keyword arguments:
    data     -- array of data to be filtered.
    bandname -- desired band name, one of ["delta","theta","alpha",
                "beta","gamma","low gamma","high gamma","lowpas","highpass"]
    fs       -- sampling frequenc
    order    -- butterworth filter order
    """
    lowcut, highcut = getbandrange(bandname)
    b, a = butter_bandpass_filter(lowcut, highcut, fs, order)
    # y = lfilter(b, a, data)  # Apply a digital filter only forward
    y = filtfilt(b, a, data)  # Apply a digital filter forward and backward
    return y


