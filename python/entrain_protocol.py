import os, re
import numpy as np
from scipy.io import savemat
import entrain
import general
import time


datadir = "../../Data"
resultdir = "../../Result/"
fs = 2000

files = []
ignore_files = [".DS_Store"]
ignore_pattern = r"^compact"
for item in os.listdir(datadir):
    if item in ignore_files:
        continue
    if re.match(ignore_pattern, item):
        continue
    
    files.append(item)
print(files)

# files = [
#     "0903-1-5-rawdata.mat"
# ]

latency_cutoffband = ["gamma"]

head_t = time.time()

for fidx, eachfile in enumerate(files):
    start_t = time.time()
    
    analysis = entrain.EntrainAnalysis(
        datadir = datadir,
        resultdir = resultdir,
        filename = eachfile,
        fs=fs,
        fftwindow  = 200,  #points
        fftoverlap = 150,  #points
        roi=(-2,5)
    )
    
    print("starting processing %d/%d: %s"%(fidx, len(files), analysis.name))
    analysis.tfdomain_analysis()
    
    for targetband in latency_cutoffband:
        analysis.latency_analysis(cutoffbans=targetband)
        latencydatafile = os.path.join(analysis.resultdir, "latency", 
                                       analysis.name+"_"+targetband+".mat")
        savemat(latencydatafile, {"name":analysis.name, 
                                  "targetband":targetband,
                                  "latency":analysis.latency})
    
    elapsedtime = time.time() - start_t
    remaintime = (len(files)-fidx-1)*elapsedtime
    print("estimated remaining time: %.2f sec"%remaintime)
    
print("all finished! elapsed time: %.2f sec"%(time.time()-head_t))