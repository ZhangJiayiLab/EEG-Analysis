"""
EEG analysis data container

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Oct 1 2018
"""

import os
from scipy.io import loadmat, savemat
import numpy as np


def dircheck(resultdir, expname):
    pass

class CompactDataContainer(object):
    """data container of compact data"""

    def __init__(self, datadir, resultdir, patientname, expname, fs, roi):

        # meta
        self.name = expname
        self.resultdir = os.path.join(resultdir, patientname)
        self.datadir = os.path.join(datadir, patientname, "EEG")

        dircheck(self.resultdir, self.name)

        # import data
        data = loadmat(os.path.join(self.datadir, "Compact", self.name+".mat"))
        self.channels = data["channels"]
        self.times = data["times"][0, :]
        self.markers = data["markers"]  #TODO: convert to dict

        # constants
        self.fs = fs
        #self.iti = iti  #TODO: calculate iti from markers
        self.roi = roi  # TODO: generate roi from iti

    def group_channel_by_marker(self, chidx, markername):
        gap = int((self.roi[1]-self.roi[0])*self.fs)
        marker = self.markers[markername][0][0][0, :]  #TODO conver to dict
        groupdata = np.zeros((len(marker), int((self.roi[1]-self.roi[0])*self.fs)))
        for idx, each in enumerate(marker):
            groupdata[idx, :] = self.channels[chidx, int(np.floor((each+self.roi[0])*self.fs)):int(np.floor((each+self.roi[0])*self.fs))+gap]

        return groupdata

class SplitDataContainer(object):
    pass