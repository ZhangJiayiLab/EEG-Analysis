clear; clc; close all
datadir = "/Users/yizhan/Desktop/Observer/Data";
filename = "compact0903-1-5-rawdata.mat";
fs = 2000;  % Hz

load(fullfile(datadir, filename))

%% 

roi = -2*fs:5*fs;
targetband = "gamma";

for idx = 1:size(channels, 1)
    fprintf("processing channel: %d / %d\n", idx, size(channels, 1))
    
    ch_split = split_datawithmarker(channels(idx, :), cue_onset, roi, fs);
    ch_bp = get_bandresponse(channels(idx, :), targetband, fs);
    ch_bp_split = split_datawithmarker(ch_bp, cue_onset, roi, fs);
    
    group_mu = mean(ch_bp_split(:, 1:2*fs), 2);
    group_sig = std(ch_bp_split(:, 1:2*fs), [], 2);
    
    [psd_pxx, psd_f] = periodogram(ch_split', [], size(ch_split,2), fs);
    
    frange = 1:0.5:200;
    average_psdg = zeros(length(frange), (7-2)*10+1);
    for i = 1:size(cue_onset,1)
        [~,~,~,ps] = spectrogram(ch_split(5,:), 2*fs, 1.9*fs, frange, fs, 'yaxis');
        average_psdg = average_psdg + log10(ps);
    end
    
    if idx == 5
        break
    end
end

plot(ch_bp_split(5,:))
hold on
plot([1 size(ch_bp_split,2)], [group_mu(5) + 3*group_sig(5) group_mu(5) + 3*group_sig(5)])
hold off

