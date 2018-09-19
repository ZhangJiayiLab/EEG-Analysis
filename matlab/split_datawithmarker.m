function groupdata = split_datawithmarker(bpdata, marker, roi, fs)
%SPLIT_DATAWITHMARKER Summary of this function goes here
%   Detailed explanation goes here

groupdata = zeros(length(marker), length(roi));

for idx = 1:length(marker)
    groupdata(idx, :) = bpdata(marker(idx)*fs+roi);
end
end

