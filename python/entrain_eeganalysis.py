import os, re
import numpy as np
import pandas as pd
from scipy.io import savemat
from eegAnalysis import EEGAnalysis
import general
import time

datadir = "../../Data"
resultdir = "../../Result/"
patientName = "Chen Zhou"
fs = 2000
latency_cutoffband = ["gamma"]

##### set target files #####
files = []
matching_pattern = r"\d{6}-.*?\.mat"
for item in os.listdir(os.path.join(datadir, patientName, "EEG", "Compact")):
    if re.match(matching_pattern, item):
        possibleResultFile = os.path.splitext(item)[0] + ".csv"
        if not os.path.exists(os.path.join(resultdir, patientName, possibleResultFile)):
            files.append(os.path.splitext(item)[0])
print(files)
input("press any key to start ...")

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
        patientname = patientName,
        expname = eachfile,
        fs=fs,
        roi=(-2,5))

    print("starting processing %d/%d: %s"%(fidx, len(files), analysis.name))

    for targetband in latency_cutoffband:
        analysis.fast_latency_detection(markername="grating", cutoffband=targetband)
        latency = analysis.latency

        analysis.fast_latency_detection(markername="entrain", cutoffband=targetband)
        latency_entrain = analysis.latency

        #latencydatafile = os.path.join(analysis.resultdir, "latency",
        #                               analysis.name+"_"+targetband+".mat")
        #savemat(latencydatafile, {"name":analysis.name,
        #                          "targetband":targetband,
        #                          "grating":latency,
        #                          "entrain":latency_entrain})

        layout = pd.read_csv(os.path.join(datadir, patientName, "EEG", "Layout", "%s-layout.csv"%patientName))
        chind = np.argsort(layout.channel.values)
        channel = layout.channel.values[chind]
        position = layout.position.values[chind]
        result = pd.DataFrame({
            "channel":channel,
            "position":position,
            "grating": latency[channel-1],
            "entrain": latency_entrain[channel-1]})
        result.to_csv(os.path.join(resultdir, patientName, eachfile+".csv"),
                index=False, sep=',', encoding='utf-8')

    elapsedtime = time.time() - start_t
    remaintime = (len(files)-fidx-1)*elapsedtime
    print("estimated remaining time: %.2f sec"%remaintime)
    exit(0)

print("all finished! elapsed time: %.2f sec"%(time.time()-head_t))
