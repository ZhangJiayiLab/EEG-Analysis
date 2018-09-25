import os, re
import numpy as np
import pandas as pd
from scipy.io import savemat
from eegAnalysis import EEGAnalysis
import eegFilter
import general
import time

datadir = "../../Data"
resultdir = "../../Result/"
patientName = "Chen Zhou"
fs = 2000
cutoffband = ["gamma"]

plot_tfdomain = True
plot_tfdomain_entrain = False
plot_bandpower = False  #TODO
plot_rawdata = False    #TODO
detect_latency = True
detect_auc = False      #TODO
detect_global_events = False  #TODO

##### set target files #####
files = []
matching_pattern = r"\d{6}-.*?\.mat"
for item in os.listdir(os.path.join(datadir, patientName, "EEG", "Compact")):
    if re.match(matching_pattern, item):
        possibleResultFile = os.path.splitext(item)[0] + ".csv"
        if not os.path.exists(os.path.join(resultdir, patientName, possibleResultFile)):
            files.append(os.path.splitext(item)[0])

files = [
    "180904-1-5"
]

print(files)
input("press any key to start ...")

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
    
    ##### time-frequency domain analysis #####
    if plot_tfdomain:
        analysis.tfdomain_analysis(np.logspace(np.log10(1), np.log10(200), 10),
                "grating", needsavemat=True, matsuffix="_grating", averaged=True, rho=5, needpreview=True)

        if plot_tfdomain_entrain:
            analysis.tfdomain_analysis(np.logspace(np.log10(1), np.log10(200), 10),
                   "entrain", needsavemat=True, matsuffix="_entrain", averaged=True, rho=5, needpreview=True)

    for targetband in cutoffband:
        
        ##### latency detection #####
        if detect_latency:
            analysis.fast_latency_detection(markername="grating", cutoffband=targetband)
            latency = analysis.latency

            analysis.entrain_latency_detection(markername="entrain", cutoffband=targetband)
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
                #"entrain": latency_entrain[channel-1]
            })
            result.to_csv(os.path.join(resultdir, patientName, "latency", eachfile+".csv"),
                    index=False, sep=',', encoding='utf-8')
            savemat(os.path.join(resultdir, patientName, "latency", eachfile+".mat"), 
                    {"channel":channel,
                     "position":position,
                     "grating": latency[channel-1],
                     "entrain": latency_entrain[channel-1]})
            
        if detect_auc:
            pass
        
        if detect_global_events:
            analysis.global_events_detection()

    elapsedtime = time.time() - start_t
    remaintime = (len(files)-fidx-1)*elapsedtime
    print("estimated remaining time: %.2f sec"%remaintime)

print("all finished! elapsed time: %.2f sec"%(time.time()-head_t))
