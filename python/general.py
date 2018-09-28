"""
General functions and tools

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Sept 20 2018
"""

import numpy as np
import os
import json


def split_datawithmarker(data, marker, roi, fs):
    """Splite data array by markers and roi."""
    groupdata = np.zeros((len(marker), int((roi[1]-roi[0])*fs)))
    gap = np.size(groupdata, 1)
    for idx, each in enumerate(marker):
        groupdata[idx, :] = data[int(np.floor((each+roi[0])*fs)):int(np.floor((each+roi[0])*fs))+gap]

    return groupdata

def dircheck(resultdir, title):
    """Check and create working directory"""
    checklist = [
        resultdir,
        os.path.join(resultdir, "tf_domain"),
        os.path.join(resultdir, "latency"),
        os.path.join(resultdir, "event_times"),
        os.path.join(resultdir, "auc"),

        os.path.join(resultdir, "preview"),
        os.path.join(resultdir, "preview", "tf_domain"),
        os.path.join(resultdir, "preview", "tf_domain", title),
        os.path.join(resultdir, "preview", "bandpower"),
        os.path.join(resultdir, "preview", "bandpower", title),
        os.path.join(resultdir, "preview", "raw"),
        os.path.join(resultdir, "preview", "raw", title),
#         os.path.join(resultdir, "preview", "bandpass"),
    ]

    for item in checklist:
        if not os.path.exists(item) or not os.path.isdir(item):
            os.mkdir(item)

    return


def group_consecutive(a, gap=1):
    ''' group consecutive numbers in an array
        modified from https://zhuanlan.zhihu.com/p/29558169'''
    return np.split(a, np.where(np.diff(a) > gap)[0] + 1)


def getprocessedfiles(resultdir, patientName):
    if os.path.isfile(os.path.join(resultdir, patientName, "processed.json")):
        with open(os.path.join(resultdir, patientName, "processed.json"), 'r') as plist:
            processedfiles = json.loads(plist.read())
        return processedfiles
    else:
        with open(os.path.join(resultdir, patientName, "processed.json"), 'w') as plist:
            plist.write("{\"processed\":[]}")
        return {"processed":[]}


def writeprocessedfile(resultdir, patientName, filename):
    processedfiles = getprocessedfiles(resultdir, patientName)
    if filename in processedfiles:
        pass
    else:
        processedfiles['processed'].append(filename)
        with open(os.path.join(resultdir, patientName, "processed.json"), 'w') as plist:
            plist.write(json.dumps(processedfiles))
        return
