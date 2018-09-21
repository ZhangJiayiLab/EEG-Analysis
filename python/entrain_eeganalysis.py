import os, re
import numpy as np
from scipy.io import savemat
from eegAnalysis import EEGAnalysis
import general
import time

datadir = "../../Data"
resultdir = "../../Result/"
fs = 2000
latency_cutoffband = ["gamma"]

##### set target files #####
files = []
matching_pattern = r"^compact"
for item in os.listdir(datadir):
    if re.match(matching_pattern, item):
        files.append(item)
print(files)

# files = [
#     "0903-1-5-rawdata.mat"
# ]


##### protocol #####
head_t = time.time()

for fidx, eachfile in enumerate(files):
    start_t = time.time()
    
    analysis = EEGAnalysis(
        datadir = datadir,
        resultdir = resultdir,
        filename = eachfile,
        fs=fs,
        roi=(-2,5))
    
    print("starting processing %d/%d: %s"%(fidx, len(files), analysis.name))
    
    for targetband in latency_cutoffband:
        analysis.fast_latency_detection(markername="grating", cutoffband=targetband)
        latency = analysis.latency
        
        analysis.fast_latency_detection(markername="entrain", cutoffband=targetband)
        latency_entrain = analysis.latency
        
        latencydatafile = os.path.join(analysis.resultdir, "latency", 
                                       analysis.name+"_"+targetband+".mat")
        savemat(latencydatafile, {"name":analysis.name, 
                                  "targetband":targetband,
                                  "grating":latency,
                                  "entrain":latency_entrain})
    elapsedtime = time.time() - start_t
    remaintime = (len(files)-fidx-1)*elapsedtime
    print("estimated remaining time: %.2f sec"%remaintime)
    
print("all finished! elapsed time: %.2f sec"%(time.time()-head_t))