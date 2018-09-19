function bandresponse = get_bandresponse(rawdata, bandname, fs)
%GET_BANDRESPONSE Summary of this function goes here
%   Detailed explanation goes here

switch bandname
    case "delta"
        bprange = [1 4];
    case "theta"
        bprange = [4 8];
    case "alpha"
        bprange = [8 12];
    case "beta"
        bprange = [13 30];
    case "gamma"
        bprange = [30 150];
    case "low gamma"
        bprange = [30 80];
    case "high gamma"
        bprange = [80 150];
    case "lowpass"
        bprange = [0.5 200];
    case "highpass"
        bprange = [200 500];
    otherwise
        bprange = [80 150];
end

% bpFilt = designfilt('bandpassfir','FilterOrder',20, ... % which order should we use???
%          'CutoffFrequency1',bprange(1),'CutoffFrequency2',bprange(2), ...
%          'SampleRate',fs);
% bandresponse = filtfilt(bpFilt, rawdata);

bandresponse = bandpass(rawdata, bprange, fs);
end

