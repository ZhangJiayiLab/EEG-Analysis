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
    groupdata = np.zeros((len(marker), (roi[1]-roi[0])*fs))
    
    for idx, each in enumerate(marker):
        groupdata[idx, :] = data[(int(each)+roi[0])*fs:(int(each)+roi[1])*fs]
    
    return groupdata

def dircheck(resultdir, title):
    """Check and create working directory"""
    checklist = [
        resultdir,
        os.path.join(resultdir, "bandpass"),
        os.path.join(resultdir, "bandpass", title),
        os.path.join(resultdir, "tfERP"),
        os.path.join(resultdir, "tfERP", title),
        os.path.join(resultdir, "latency")
    ]
    
    for item in checklist:
        if not os.path.exists(item) or not os.path.isdir(item):
            os.mkdir(item)
    
    return

def group_consecutive(a, gap=1):
    ''' group consecutive numbers in an array
        modified from https://zhuanlan.zhihu.com/p/29558169'''
    return np.split(a, np.where(np.diff(a) > gap)[0] + 1)