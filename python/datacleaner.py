"""
Re-organize the data structure for python analysis

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Sept 20 2018
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
    return os.path.isfile(os.path.join(datadir, "compact"+filename))


def processCompact(datadir, filename, prefix="compact", overwrite=False, renderpreview=False):
    """Convert structure of raw data for python analysis"""
    if not checkCompact(datadir, filename) or overwrite:
        with h5py.File(os.path.join(datadir,filename), 'r') as f:
            channels = np.zeros((len(f), len(f["Chan__1"]["times"][0])))
            times = np.array(f["Chan__1"]["times"])[0,:]
            cueonset = np.array(f["Memory"]["times"])[0,:]
            
            for i in f.keys():
                if i == "Memory" or i == "file":
                    continue
                idx = int(np.array(f[i]["title"], dtype="uint8")[4:].tostring())
                channels[idx-1, :] = np.array(f[i]["values"])[0,:]
        
        savemat(os.path.join(datadir, prefix+filename), 
                {"times":times, 
                 "cue_onset":cueonset, 
                 "channels":channels, 
                 "markers":{
                     "grating":cueonset, 
                     "entrain":[cueonset[-1]+5, cueonset[-1]+10]}
                })
            
    
if __name__ == "__main__":
    datadir = sys.argv[1]
    
    ignorefiles = [".DS_Store"]
    ignore_pattern = "^compact"
    match_pattern = r".*?\.mat"
    
    rawfiles = []
    for item in os.listdir(datadir):
        if item in ignorefiles:
            continue
        if re.match(ignore_pattern, item):
            continue
        
        if re.match(match_pattern, item):
            rawfiles.append(item)
    
    print(rawfiles)
    for idx, item in enumerate(rawfiles):
        print("%d/%d: %s"%(idx, len(rawfiles), item))
        if not checkCompact(datadir, item):
            processCompact(datadir, item)