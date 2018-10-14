import os, re
import EEGAnalysis as ea
import general
import numpy as np
import matplotlib.pyplot as plt

datadir = "../../Data"
resultdir = "../../Result"
patientName = "Chen Zhou"
Fs = 2000
cutoffband = ["gamma"]

##### set target files #####
files = []
ignorefiles = []
matching_pattern = r"(\d{6})-(.*?)\.mat"
processedfiles = general.getprocessedfiles(resultdir, patientName)
for item in os.listdir(os.path.join(datadir, patientName, "EEG", "Compact")):
    if re.match(matching_pattern, item):
        possiblefile = os.path.splitext(item)[0]
        if possiblefile not in processedfiles["processed"] and possiblefile not in ignorefiles:
            files.append(possiblefile)

print(files)
input("press any key to start ...")

##### import layout #####
layout = pd.read_csv(os.path.join(datadir, patientName, "EEG", "Layout", "%s-layout.csv"%patientName))
chind = np.argsort(layout.channel.values)
layout_channel = layout.channel.values[chind]
layout_position = layout.position.values[chind]

##### protocol #####
head_t = time.time()

for fidx, eachfile in enumerate(files):
    start_t = time.time()

    container = ea.container.CompactDataContainer(
        datadir, resultdir, patientName, eachfile, Fs)

    print("starting processing %d/%d: %s"%(fidx, len(files), eachfile))
    print("ITI = %.1f"%container.iti)

    ##### decomposition analysis #####
    frange = np.logspace(np.log10(1), np.log(200), 30)
    tf_decmp_grating = {}
    tf_decmp_entrain = {}
    for chidx in range(np.size(container.channels, 0)):
        grating_split = container.group_channel_by_marker(chidx, "grating")
        entrain_split = container.group_channel_by_marker(chidx, "entrain")

        tf_decmp_grating[chidx] = ea.decomposition.dwt(grating_split, Fs, frange)
        tf_decmp_entrain[chidx] = ea.decomposition.dwt(entrain_split, Fs, frange)

    ##### total power #####
    for chidx in range(np.size(container.channels, 0)):
        grating_pwr = ea.power.dwt_power(tf_decmp_grating[chidx], Fs)
        entrain_pwr = ea.power.dwt_power(tf_decmp_entrain[chidx], Fs)

        pass
        #TODO:
