"""
Probe layout information

author: Yizhan Miao
email: yzmiao@protonmail.com
last update: Sept 20 2018
"""

import pandas as pd

class ProbeLayout(object):
    def __init__(self, layout_dic="dic.csv"):
        self.layout = pd.read_csv(layout_dic)
        
    def findchannel(self, area):
        return self.layout[self.layout.position == area].channel.values
    
    def findarea(self, channel):
        return self.layout[self.layout.channel == channel].position.values