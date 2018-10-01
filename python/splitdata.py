import os, re
from scipy.io import loadmat, savemat
from fastprogress import master_bar, progress_bar
# import numpy as np

datadir = "../../Data/Chen Zhou/EEG/Compact/"
exportdir = "../../Data/Chen Zhou/EEG/SgCh/"

chrange = range(123)

files = []
matching_pattern = r"\d{6}-.*?\.mat"
for item in os.listdir(datadir):
    if re.match(matching_pattern, item):
        possibleResultFile = os.path.splitext(item)[0]
        files.append(possibleResultFile)
print(len(files))

mb = master_bar(chrange)
for chidx in mb:
    exportfiledir = os.path.join(exportdir, "ch%03d"%(chidx+1))
    if not os.path.isdir(exportfiledir):
        os.mkdir(exportfiledir)

    for each_source in progress_bar(files, parent=mb):
        exportmat = os.path.join(exportfiledir, each_source+'_ch%03d.mat'%chidx)
        sourcemat = os.path.join(datadir, each_source+'.mat')
        if os.path.isfile(exportmat):
            continue

        rawmat = loadmat(sourcemat)
        savemat(exportmat, {
            "values": rawmat['channels'][chidx, :],
            "times":  rawmat["times"][0,:],
            "markers": rawmat["markers"]
        })
