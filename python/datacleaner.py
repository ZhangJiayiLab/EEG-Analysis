"""
Re-organize the data structure for python analysis

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Sept 23 2018
"""

# %matplotlib inline
import sys, os, re
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import savemat
from tqdm import tqdm

def previewChannel(previewdir, chn, save=False):
    plt.figure(figsize=(20,3))
    plt.plot(times, channels[chn,:])
    setoff = np.max(channels[chn,:])
    plt.vlines(cueonset, setoff, setoff+100, colors = "red")
    plt.xlim((0, np.max(times)))
    plt.title("channel %d"%chn)

    plt.savefig(os.path.join(previewdir,"%d.png"%chn))
    plt.close()


def checkCompact(datadir, filename):
    """Check if compact data file exists, return bool."""
    return os.path.isfile(os.path.join(datadir, filename))


def processCompact(datadir, subdir, filename, prefix="Compact", fs=2000, overwrite=False, renderpreview=False):
    """Convert structure of spike2 raw data for python analysis"""
    if not checkCompact(os.path.join(datadir, prefix), filename) or overwrite:
        with h5py.File(os.path.join(datadir, subdir, filename), 'r') as f:
            channels = np.zeros((len(f), len(f["Chan__1"]["times"][0])))
            times = np.array(f["Chan__1"]["times"])[0,:]
            cueonset = np.array(f["Memory"]["times"])[0,:]

            for i in f.keys():
                if i == "Memory" or i == "file":
                    continue
                idx = int(np.array(f[i]["title"], dtype="uint8")[4:].tostring())
                channels[idx-1, :] = np.array(f[i]["values"])[0,:]

        length = np.size(channels, 1) / fs
        residual = length - cueonset[-1]
        entrain_cue = [cueonset[-1]+i*5+5 for i in range(min(2, int(residual//5-1)))]
        savemat(os.path.join(datadir, prefix, filename),
                {"times":times,
                 "cue_onset":cueonset,
                 "channels":channels,
                 "markers":{
                     "grating":cueonset,
                     "entrain":entrain_cue}
                })
        print(entrain_cue)


if __name__ == "__main__":
    datadir = sys.argv[1]
    patientName = sys.argv[2]
    datadir = os.path.join(datadir, patientName, "EEG")

    ignorefiles = [".DS_Store"]
    match_pattern = r"\d{6}-.*?\.mat"

    rawfiles = []
    for item in os.listdir(os.path.join(datadir,"Spike2")):
        if item in ignorefiles:
            continue

        if re.match(match_pattern, item):
            if not checkCompact(os.path.join(datadir, "Compact"), item):
                rawfiles.append(item)

    print(rawfiles)
    input("press any key to start ... ")
    for idx, item in enumerate(rawfiles):
        print("%d/%d: %s"%(idx, len(rawfiles), item))
        processCompact(datadir, "Spike2", item)

