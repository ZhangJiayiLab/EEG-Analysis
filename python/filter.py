from scipy.signal import butter, lfilter
import numpy as np

def getbandrange(bandname):
    bandrange = {
            "delta":(  2,  4),
            "theta":(  4,  8),
            "alpha":(  8, 12),
            "beta": ( 15, 30),
            "gamma":( 30,150),
            "low gamma": (30, 80),
            "high gamma":(80,150),
            "lowpass":(1,200),
            "highpass"(200,500)}
    return bandrange[bandname]

def butter_bandpass_filter(lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = higcut / nyq
    b, a = butter(order, [low, high], btype="band")
    return b,a

def butter_bandpass(data, bandname, fs, order=5):
    lowcut, highcut = getbandrange(bandname)
    b, a = butter_bandpass_filter(lowcut, highcut, fs, order)
    y = lfilter(b, a, data)
    return y


