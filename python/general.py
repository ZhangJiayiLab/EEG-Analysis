"""
General functions and tools

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Sept 20 2018
"""

import numpy as np
import os


def split_datawithmarker(data, marker, roi, fs):
    """Splite data array by markers and roi."""
    groupdata = np.zeros((len(marker), int(roi[1]-roi[0])*fs))
#     print(roi)
    for idx, each in enumerate(marker):
#         print(np.shape(data[int((each+roi[0])*fs):int((each+roi[1])*fs)]))
        groupdata[idx, :] = data[np.floor((each+roi[0])*fs):np.floor((each+roi[1])*fs)]

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
        os.path.join(resultdir, "preview", "bandpower"),
        os.path.join(resultdir, "preview", "raw"),
        os.path.join(resultdir, "preview", "bandpass"),
    ]

    for item in checklist:
        if not os.path.exists(item) or not os.path.isdir(item):
            os.mkdir(item)

    return

def group_consecutive(a, gap=1):
    ''' group consecutive numbers in an array
        modified from https://zhuanlan.zhihu.com/p/29558169'''
    return np.split(a, np.where(np.diff(a) > gap)[0] + 1)
