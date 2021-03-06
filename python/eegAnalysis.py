import os
import numpy as np
from scipy.io import loadmat, savemat
import matplotlib.pyplot as plt
from tqdm import tqdm

import general
import stfft
import eegFilter

class EEGAnalysis(object):

    def __init__(self, datadir, resultdir, patientname, expname, fs, roi):

        # meta data
        self.name = expname
        self.resultdir = os.path.join(resultdir, patientname)
        self.datadir = os.path.join(datadir, patientname, "EEG")

        general.dircheck(self.resultdir, expname)

        # analysis parameters
        self.fs = fs
        self.roi = roi

        # analysis constant
        self.bandname = ["delta", "theta", "alpha", "beta",
            "gamma", "low gamma", "high gamma",
            "lowpass", "highpass"]

        # import data
        data = loadmat(os.path.join(self.datadir, "Compact", self.name+".mat"))
        self.channels = data["channels"]
        #self.cue_onset = data["cue_onset"][0,:]
        self.times = data["times"][0,:]
        self.markers = data["markers"]

        # results
        self.tfdata = None
        self.latency = None
        self.rho = None


    def tfdomain_analysis(self, frange, markername, baseroi=(1,2), rho=5,
                          needsavemat=False, matsuffix="_tf", averaged=False,
                          needpreview=False, layout=None):

        cue_onset = self.markers[markername][0][0][0,:]
        ch_split = general.split_datawithmarker(
                self.channels[1, :], cue_onset, self.roi, self.fs)
        power_channels = np.zeros((np.size(self.channels,0),
                                   np.size(frange), np.size(ch_split, 1)))

        print("processing time-frequency domain analysis")
        for chidx in tqdm(range(np.size(self.channels, 0)), ascii=True):
            ch_split = general.split_datawithmarker(
                    self.channels[chidx, :], cue_onset, self.roi, self.fs)
            pxx = stfft.dwt_tf(ch_split, self.fs, frange, reflection=True)

            power_channels[chidx, :, :] = pxx

        self.tfdata = power_channels
        tspec = np.linspace(self.roi[0], self.roi[1], np.size(self.tfdata, 2))

        if needsavemat:
            savemat(os.path.join(self.resultdir, "tf_domain", self.name + matsuffix + ".mat"),
                {"name": self.name,
                 "method": "morlet wavelet",
                 "power_channels": self.tfdata[:,:,::rho],
                 "times": tspec[::rho],
                 "frange": frange,
                 "rho": rho,
                 })
        if needpreview:
            if not os.path.isdir(os.path.join(self.resultdir, 'preview', 'tf_domain', self.name)):
                os.mkdir(os.path.join(self.resultdir, 'preview', 'tf_domain', self.name))

            print("start rendering contourf preview of tf_domain ...")
            for chidx in tqdm(range(np.size(self.tfdata, 0)), ascii=True):
                plt.figure(figsize=(5,3))
                plt.contourf(tspec, frange, self.tfdata[chidx, :, :], 30, cmap=plt.get_cmap("jet"))

                if type(layout) != type(None):
                    if len(layout[layout.channel==chidx+1]) == 1:
                        plt.title("ch%d - %s"%(chidx, layout[layout.channel==chidx+1].position.values[0]))
                    else:
                        plt.title("ch%d - n.a."%chidx)

                if markername == "entrain":
#                     plt.clim([-5,5])
                    pass
                plt.savefig(os.path.join(self.resultdir, 'preview', 'tf_domain', self.name, "ch%03d"%chidx+matsuffix+'.png'), bbox_inches='tight')
                plt.close()


    def fast_latency_detection(self, markername, cutoffband="gamma",
                               n=3, gapsec = 0.1, rho=0.5, baseroi=(1.5,2), validlen=0.1,
                               needpreview=False):
        cutoff = eegFilter.getbandrange(cutoffband)
        frange = np.logspace(np.log10(cutoff[0]), np.log10(cutoff[1]), 20)
        cue_onset = self.markers[markername][0][0][0,:]

        auc_on = np.zeros((np.size(self.channels,0), np.size(cue_onset)))

        if needpreview and not os.path.isdir(os.path.join(self.resultdir, 'preview', 'bandpower', self.name)):
            os.mkdir(os.path.join(self.resultdir, 'preview', 'bandpower', self.name))

        datarise = []
        print("processing fast latency analysis for %s"%markername)
        for chidx in tqdm(range(np.size(self.channels, 0)), ascii=True):
            ch_split = general.split_datawithmarker(
                        self.channels[chidx, :], cue_onset, self.roi, self.fs)
            pxx = stfft.dwt_tf(ch_split, self.fs, frange, baseroi, reflection=True)

            power_channels = np.mean(pxx, 0)
            tfcurve = power_channels
            tspec = np.linspace(self.roi[0], self.roi[1], (self.roi[1]-self.roi[0])*self.fs)

            tfcurve_mu = np.mean(tfcurve[np.where(tspec < 0)])
            tfcurve_shift = tfcurve - tfcurve_mu
            tfcurve_sig = np.std(tfcurve[np.where(tspec < 0)])
            thresh = n*tfcurve_sig

            if needpreview:
                Ntrial = np.size(pxx, 0)
                plt.figure(figsize=(10,Ntrial))
                for trialn in range(Ntrial):
                    plt.subplot(Ntrial,1,trialn+1)
                    plt.plot(tspec, pxx[trialn, :])
                    plt.vlines(0, np.min(pxx[trialn, :]), np.max(pxx[trialn, :]), linewidth=2)
                    plt.xlim([tspec[0], tspec[-1]])
                plt.savefig(os.path.join(self.resultdir, 'preview', 'bandpower', self.name, "ch%03d_%s_%s.png"%(chidx, markername, cutoffband)), bbox_inches='tight')
                plt.close()

            datarise_item = general.group_consecutive(np.where(tfcurve_shift > thresh)[0], gap=gapsec*self.fs)
            if len(datarise_item[0]) < validlen*self.fs:
                datarise.append([np.nan])
            else:
                points = []
                for each in datarise_item:
                    x0, y0 = tspec[each[0]], tfcurve_shift[each[0]]
                    xm1, ym1 = tspec[each[0]-1], tfcurve_shift[each[0]-1]
                    points.append(x0 - y0 * (x0-xm1)/(y0-ym1))
                datarise.append(points)

        ## latency
        latency = np.array([])
        for idx, item in enumerate(datarise):
            point = np.nan
            for each in item:
                if each>0:
                    point = each
                    break
            latency = np.append(latency, point)
        self.latency = latency
        self.events = datarise
#         self.auc_on =
#         self.auc_off =

    def entrain_latency_detection(self, markername="entrain", cutoffband="gamma",
                                  n=3, gapsec=0.1, validlen=0.1, rho=2, iti=5):
        cutoff = eegFilter.getbandrange(cutoffband)
        frange = frange = np.logspace(np.log10(cutoff[0]), np.log10(cutoff[1]), 20)
        cue_onset = self.markers[markername][0][0][0, :]
        tspec = np.linspace(self.roi[0], self.roi[1], (self.roi[1]-self.roi[0])*self.fs)
        baseline_cue_onset = np.array([cue_onset[0]-iti*1,
                                      cue_onset[0]-iti*2,
                                      cue_onset[0]-iti*3,
                                      cue_onset[0]-iti*4,
                                      cue_onset[0]-iti*5])

        datarise = []

        print("processing entrain latency analysis for %s"%markername)
        for chidx in tqdm(range(np.size(self.channels, 0)), ascii=True):
            baseline_ch_split = general.split_datawithmarker(
                        self.channels[chidx, :], baseline_cue_onset, self.roi, self.fs)
            baseline_pxx = stfft.dwt_tf(baseline_ch_split, self.fs, frange, reflection=True)
            baseline_mu = np.mean(baseline_pxx[:, np.where(tspec < 0)])
            baseline_sig = np.std(baseline_pxx[:, np.where(tspec < 0)])
            thresh = baseline_mu + n*baseline_sig

            ch_split = general.split_datawithmarker(
                        self.channels[chidx, :], cue_onset, self.roi, self.fs)
            pxx = stfft.dwt_tf(ch_split, self.fs, frange, reflection=True)
            ch_power = np.mean(pxx, 0)

            datarise_item = general.group_consecutive(np.where(ch_power > thresh)[0], gap=gapsec*self.fs)
            if len(datarise_item[0]) < validlen*self.fs:
                datarise.append([np.nan])
            else:
                points = []
                for each in datarise_item:
                    x0, y0 = tspec[each[0]], ch_power[each[0]]
                    xm1, ym1 = tspec[each[0]-1], ch_power[each[0]-1]
                    points.append(x0 - y0 * (x0-xm1)/(y0-ym1))
                datarise.append(points)

        ## latency
#         latency = np.array([])
#         for idx, item in enumerate(datarise):
#             point = np.nan
#             for each in item:
#                 if each>0:
#                     point = each
#                     break
#             latency = np.append(latency, point)
        self.latency = np.array(datarise)
        self.events = np.array(datarise)


    def global_events_detection(self, markername="grating", cutoffband="gamma",
                                  n=3, gapsec=0.1, validlen=0.1, rho=2, iti=5):
        cutoff = eegFilter.getbandrange(cutoffband)
        frange = frange = np.logspace(np.log10(cutoff[0]), np.log10(cutoff[1]), 20)
        cue_onset = self.markers[markername][0][0][0, :]
        tspec = np.linspace(self.roi[0], self.roi[1], (self.roi[1]-self.roi[0])*self.fs)


        datarise = []
        print("global events detections")
        for chidx in tqdm(range(np.size(self.channels, 0)), ascii=True):
            baseline = general.split_datawithmarker(self.channels[chidx, :],
                                                cue_onset, self.roi, self.fs)
            baseline_mu = np.mean(baseline[:,np.where(tspec<0)])
            baseline_sig = np.std(baseline[:,np.where(tspec<0)])
            thresh = baseline_mu + n*baseline_sig

            pxx = stfft.dwt_tf(self.channels[chidx:chidx, :], self.fs, frange, reflection=False)
            pxx[:500,:] = 0
            pxx[-1000:,:] = 0
            datarise_item = general.group_consecutive(np.where(pxx > thresh)[0], gap=gapsec*self.fs)
            if len(datarise_item[0]) < validlen*self.fs:
                datarise.append([np.nan])
            else:
                points = []
                for each in datarise_item:
                    x0, y0 = tspec[each[0]], tfcurve_shift[each[0]]
                    xm1, ym1 = tspec[each[0]-1], tfcurve_shift[each[0]-1]
                    points.append(x0 - y0 * (x0-xm1)/(y0-ym1))
                datarise.append(points)
            if chidx == 6:
                break

        return datarise


    def auc_detection(self, markername, cutoffband="gamma",
                               n=3, gapsec = 0.1, rho=0.5, baseroi=(1.5,2), validlen=0.1,
                               needpreview=False):
        cutoff = eegFilter.getbandrange(cutoffband)
        frange = np.logspace(np.log10(cutoff[0]), np.log10(cutoff[1]), 20)
        cue_onset = self.markers[markername][0][0][0,:]
        tspec = np.linspace(self.roi[0], self.roi[1], (self.roi[1]-self.roi[0])*self.fs)

        auc_on = np.zeros((np.size(self.channels,0), np.size(cue_onset)))
        auc_off = np.zeros((np.size(self.channels,0), np.size(cue_onset)))

        print("auc detection for %s"%markername)
        for chidx in tqdm(range(np.size(self.channels, 0)), ascii=True):
            ch_split = general.split_datawithmarker(
                        self.channels[chidx, :], cue_onset, self.roi, self.fs)
            pxx = stfft.dwt_tf(ch_split, self.fs, frange, reflection=True, needaverage=False)

            auc_on[chidx, :] = np.squeeze(np.sum(np.sum((pxx-np.mean(pxx, 0))/np.std(pxx, 0), 0)[:,np.where((0<tspec)&(tspec<2))],2))
            auc_off[chidx, :] = np.squeeze(np.sum(np.sum((pxx-np.mean(pxx, 0))/np.std(pxx, 0), 0)[:,np.where((2<tspec)&(tspec<2.5))],2))

        self.auc_on = auc_on
        self.auc_off = auc_off


    def bandpower_curve_preview(self, bandname="all", markername="entrain", iti=5):

        cue_onset = self.markers[markername][0][0][0,:]
        tspec = np.linspace(self.roi[0], self.roi[1], (self.roi[1]-self.roi[0])*self.fs)

        if bandname == "all":
            bandnames = self.bandname
        else:
            bandnames = [bandname]

        for chidx in tqdm(range(np.size(self.channels, 0)), ascii=True):
            power_curve = []

            for cutoffband in bandnames:
                cutoff = eegFilter.getbandrange(cutoffband)
                frange = np.logspace(np.log10(cutoff[0]), np.log10(cutoff[1]), 5)

                ch_split = general.split_datawithmarker(self.channels[chidx, :],
                                                        cue_onset, self.roi, self.fs)
                pxx = stfft.dwt_tf(ch_split, self.fs, frange, reflection=True, zscore=True)
                power_curve.append(np.mean(pxx, 0))

            plt.figure(figsize=(3,1*len(bandnames)))
            for pidx in range(len(bandnames)):
                plt.subplot(len(bandnames), 1, pidx+1)
                plt.plot(tspec[::100], power_curve[pidx][::100], 'k')
                plt.title(bandnames[pidx])
                plt.xlim(tspec[0], tspec[-1])

            plt.savefig(os.path.join(self.resultdir, 'preview', 'bandpower', self.name, "ch%03d"%chidx+markername+'.png'), bbox_inches='tight')
            plt.close()



    def bandpower_auc_detection(self, bandname="all", markername="entrain", iti=5, auc_roi=(-2,2)):

        cue_onset = self.markers[markername][0][0][0,:]
        tspec = np.linspace(self.roi[0], self.roi[1], (self.roi[1]-self.roi[0])*self.fs)

        if bandname == "all":
            bandnames = ["delta", "theta", "alpha", "beta", "low gamma", "high gamma", "highpass"]
        else:
            bandnames = [bandname]

        exportfile = os.path.join(self.resultdir, "auc", "bandpower_auc.csv")
        exportfile_raw = os.path.join(self.resultdir, "auc", "bandpower_auc_raw.csv")
        exportentry = "{exp},{channel},{id},{mode},{delta},{theta},{alpha},{beta},{lowgamma},{highgamma},{highpass}\n"

        if not os.path.isfile(exportfile):
            with open(exportfile, 'w') as csvf:
                csvf.write("exp,channel,id,mode,delta,theta,alpha,beta,lowgamma,highgamma,highpass\n")

        if not os.path.isfile(exportfile_raw):
            with open(exportfile_raw, 'w') as csvf:
                csvf.write("exp,channel,id,mode,delta,theta,alpha,beta,lowgamma,highgamma,highpass\n")

        for chidx in tqdm(range(np.size(self.channels, 0)), ascii=True):
            power_curve = {}

            ch_split = general.split_datawithmarker(self.channels[chidx, :],
                                                    cue_onset, self.roi, self.fs)
            for cutoffband in bandnames:
                cutoff = eegFilter.getbandrange(cutoffband)
                frange = np.logspace(np.log10(cutoff[0]), np.log10(cutoff[1]), 5)

                pxx = stfft.dwt_tf(ch_split, self.fs, frange, reflection=True, zscore=True)
                power_curve[cutoffband] = np.sum(np.mean(pxx, 0)[np.where((tspec>auc_roi[0])&(tspec<auc_roi[1]))])


            with open(exportfile_raw, 'a') as csvf:
                csvf.write(exportentry.format(
                    exp = self.name,
                    channel = chidx+1,
                    id = self.name+'-ch%03d'%(chidx+1),
                    mode = markername,
                    delta = power_curve['delta'],
                    theta = power_curve['theta'],
                    alpha = power_curve['alpha'],
                    beta = power_curve['beta'],
                    lowgamma = power_curve['low gamma'],
                    highgamma = power_curve['high gamma'],
                    highpass = power_curve['highpass'],
                ))

            with open(exportfile, 'a') as csvf:
                csvf.write(exportentry.format(
                    exp = self.name,
                    channel = chidx+1,
                    id = self.name+'-ch%03d'%(chidx+1),
                    mode = markername,
                    delta = power_curve['delta']/self.fs,
                    theta = power_curve['theta']/self.fs,
                    alpha = power_curve['alpha']/self.fs,
                    beta = power_curve['beta']/self.fs,
                    lowgamma = power_curve['low gamma']/self.fs,
                    highgamma = power_curve['high gamma']/self.fs,
                    highpass = power_curve['highpass']/self.fs,
                ))
