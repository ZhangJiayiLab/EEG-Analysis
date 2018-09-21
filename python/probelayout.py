"""
Probe layout information

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Sept 21 2018
"""

import os
import pandas as pd

class ProbeLayout(object):
    def __init__(self, datadir, layout_dic=):
        self.layout = pd.read_csv(os.path.join(datadir, layout_dic))

    def findchannel(self, area):
        return self.layout[self.layout.position == area].channel.values

    def findarea(self, channel):
        return self.layout[self.layout.channel == channel].position.values
