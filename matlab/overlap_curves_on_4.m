clear all;clc;close all;
format long;

rawdata = load('D:\SYF\0901-1-10-lowgamma-0.01.mat');%更改文件及目录
channel_name = fieldnames(rawdata);
channel_num = length(channel_name)-1;
grating_interval = 10;%这里要记得，10s 改成10
sampleinterval = 0.01;
grating_time = rawdata.Memory.times;
xtime = 0:sampleinterval:10;%这里要记得，10s 改成10！！！！！！！！！

res = 0;
nanres = 0;
response_latency = 0;

for i = 2:channel_num
    rawdataa = getfield(rawdata,channel_name{i});
    power = rawdataa.values;
    
    
    for k = 1:length(grating_time)-1
        powertrace(:,k) = power(grating_time(k,1)/sampleinterval:grating_time(k+1,1)/sampleinterval,1)/...
            max(power(grating_time(k,1)/sampleinterval:grating_time(k+1,1)/sampleinterval,1));
        powertrace_raw(:,k) = power(grating_time(k,1)/sampleinterval:grating_time(k+1,1)/sampleinterval,1);
    end
    figure(1);
    s = subplot(12,12,i); 
    averagepower = mean(powertrace,2);
    averagepowertrace_raw = mean(powertrace_raw,2);
    
    errBar = (std(powertrace,0,2))/sqrt(k);
    p = plot(xtime,mean(powertrace,2));
    
    ShadedErrorBar(xtime,averagepower,errBar);
    set(gca,'xtick',[],'xtickLabel',[]);
    set(gca,'ytick',[],'ytickLabel',[]); 
    title(rawdataa.title);
    xlswrite(strcat('C:\Users\lxd\Desktop\human data\response curve lowgamma\ ',num2str(i)),averagepower);
    xtime_on = 0:sampleinterval:2;
    xtime_off = 0:sampleinterval:grating_interval/2.5;
    grating_on = averagepower(1:(2/sampleinterval)+1);% 光刺激引起的averagepower
    grating_off = averagepower(((grating_interval/2.5)/sampleinterval)+1:((grating_interval/1.25)/sampleinterval)+1); %光刺激结束后2-4秒内的averagepower，
    area_on = trapz(xtime_on,grating_on); 
    area_off = (trapz(xtime_off,grating_off))/(grating_interval/2.5);
    
    %判断该通道是否具有光反应
    if area_on/area_off >= 2.3; %这里需要根据反应的信噪比调整
        res = res+1;
        response{res,1} = rawdataa.title;
    else
        nanres=nanres+1;
        nonresponse{nanres,1} = rawdataa.title;
    end  
    
    %计算光反应的amplitude
    response_amplitude = max(averagepowertrace_raw(1:2/sampleinterval+1));
    xlswrite(strcat('C:\Users\lxd\Desktop\human data\response_amplitude lowgamma\ ',num2str(i)),response_amplitude);
    
    %计算光反应的latency
    baseline = averagepowertrace_raw((grating_interval-2)/sampleinterval:grating_interval/sampleinterval);
    baseline_average = mean(baseline);
    baseline_sd = std(baseline,0,1);
    
    for n = 1:length(averagepowertrace_raw(1:300,1));
        responsevalue = averagepowertrace_raw(n,1);
        if responsevalue > baseline_average+3*baseline_sd;
           responselatency = responsevalue;
           latencytime = n*sampleinterval;
           response_latency = response_latency+1;
           latencyexport{response_latency,1} = latencytime;
           xlswrite(strcat('C:\Users\lxd\Desktop\human data\latency lowgamma\',num2str(i)),latencytime); 
        end
    end
end
clc;
%     for j = 1:i;
%         cellReference = sprintf('j%s',num2str(j));
%         xlswrite('1-5s.xlsx',averagepower,cellReference);
%     end
 
    
% xtime = 0.25:0.25:5;
% averagepower = mean(powertrace,2);
% errBar = (std(powertrace,0,2))/sqrt(length);
% plot(xtime,mean(powertrace,2));
% ShadedErrorBar(xtime,averagepower,errBar);
% xlswrite('1-5s.xlsx',powertrace);



% %读取同一文件件下所有mat文件
% folder = 'E:\matfiles\20180831\1-5s';
% filetype = '*.mat';
% f = fullfile(folder,filetype);
% files = dir(f);

%     
% %对同一文件夹下面的同类型文件按名字进行排序
% nameCell = cell(length(files)-2,1);
% for i = 1:length(files);
%     disp(files(i).name);
%     nameCell{i} = files(i).name;
% end
% newfiles = sort_nat(nameCell);






