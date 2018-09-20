% mikexcohen.com
clear;clc;close all
Rawdata=load('saliency.mat');

Allsame=Rawdata.all_same.times;
Saliency=Rawdata.saliency.times;
Strucname=fieldnames(Rawdata);
Channame=fieldnames(Rawdata);
Channum=length(Channame);
Fs=2000;
freqmax=300;
freqmin=2;
for n=1:Channum-3
    Field=getfield(Rawdata,Channame{n});
    value=Field.values;
    for stinum=1:length(Allsame)
        allsamestitime=Allsame(stinum,1);
        saliencystitime=Saliency(stinum,1);
        Allsame_powervalue(n,stinum)=bandpower(value(round((allsamestitime-0.2)*Fs)+1:round(allsamestitime*Fs)+2*Fs),Fs,[freqmin,freqmax]);
        Saliency_powervalue(n,stinum)=bandpower(value(round((allsamestitime-0.2)*Fs)+1:round(allsamestitime*Fs)+2*Fs),Fs,[freqmin,freqmax]);
        Allsametrace(n,:,stinum)=value(round((allsamestitime-0.2)*Fs)+1:round(allsamestitime*Fs)+2*Fs);
        Saliencytrace(n,:,stinum)=value(round((saliencystitime-0.2)*Fs)+1:round(saliencystitime*Fs)+2*Fs);
    end
end
Allsame_plot=mean(Allsametrace,3);
Saliency_plot=mean(Saliencytrace,3);
Allsame_std=std(mean(Allsametrace,3),0,1)/Channum;
Saliency_std=std(mean(Saliencytrace,3),0,1)/Channum;
aveplotx=-0.2+1/Fs:1/Fs:2.2;
figure(41)
xplot1=zeros(101,1)+400;
yplot=0:100;

imagesc(-1*Allsame_plot);
hold on
plot(xplot1,yplot,'--','linewidth',4)
hold on
plot(mean(Allsame_plot)/2+80);

figure(42)
imagesc(-1*Saliency_plot);
hold on
plot(xplot1,yplot,'--','linewidth',2)
hold on
plot(mean(Saliency_plot)/2+80);
order=1;

%% prepare details and parameters for time-frequency decomposition
xplot=-199.5:0.5:2000;
load sampleEEGdata
for sb=83:83
    Rawtrace(1,:,:)=Allsametrace(sb,:,:);
    Rawtrace(2,:,:)=Allsametrace(sb,:,:);
    % 
    % % plot ERP on top
    % hold on
    % plot(xplot,erpt,'k')

    %% end.


    timewin    = 400; % in ms
    timewinidx = round(timewin/(3000/Fs));
    tapers     = dpss(timewinidx,3); % this line will crash without matlab signal processing toolbox
    channel2plot = 'o1';
    times2save   = -200:25:2000;
    basetime     = [-300 -100];

    % convert time points to indices
    times2saveidx = dsearchn(EEG.times',times2save'); 

    % find baseline time point range
    baseidx = dsearchn(times2save',basetime');

    % define frequencies for FFT
    hz = linspace(0,Fs,timewinidx/2+1);



    % initialize output matrix
    tf = zeros(floor(timewinidx/2)+1,length(times2save));
    tf1 = zeros(floor(timewinidx/2)+1,length(times2save));
    % loop through time bins
    for ti=1:length(times2saveidx)

        % initialize power vector (over tapers)
        taperpow = zeros(floor(timewinidx/2)+1,1);
        taperpow1 = zeros(floor(timewinidx/2)+1,1);
        % loop through tapers
        for tapi = 1:size(tapers,2)-1

            % get data from this time window and taper
            tempEEG  = squeeze(Rawtrace(1,times2saveidx(ti)-floor(timewinidx/2)+1:times2saveidx(ti)+ceil(timewinidx/2),:));
            tempEEG1  = squeeze(Rawtrace(2,times2saveidx(ti)-floor(timewinidx/2)+1:times2saveidx(ti)+ceil(timewinidx/2),:));
            data     = bsxfun(@times,tempEEG,tapers(:,tapi));
            data1     = bsxfun(@times,tempEEG1,tapers(:,tapi));
            % compute FFT and extract power
            pow      = fft(data)/timewinidx;
            pow      = pow(1:floor(timewinidx/2)+1,:);
            taperpow = taperpow + mean(abs(pow).^2,2);
            pow1      = fft(data1)/timewinidx;
            pow1      = pow1(1:floor(timewinidx/2)+1,:);
            taperpow1 = taperpow1 + mean(abs(pow1).^2,2);
        end

        % divide by N tapers for average
        tf(:,ti) = taperpow/tapi;
        tf1(:,ti) = taperpow1/tapi;
    end

    % db-correct
    tf = 10*log10( bsxfun(@rdivide,tf,mean(tf(:,baseidx(1):baseidx(2)),2)) );
    tf1 = 10*log10( bsxfun(@rdivide,tf1,mean(tf1(:,baseidx(1):baseidx(2)),2)));

    % plot!
    figurenum=ceil(sb/5);
    if rem(sb,5)~=0
        subnum=rem(sb,5);
    else
        subnum=rem(sb,5)+5;
    end
    figure(figurenum)
    colormap(hsv)
    set(gcf,'outerposition',get(0,'screensize'));
    subplot(5,2,subnum*2-1)
    contourf(times2save,hz,tf/abs(max(max(max(tf(:,:,:))))),40,'linecolor','none')
    set(gca,'clim',[-2 2],'ylim',[2 150])
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    subplot(5,2,subnum*2)
    contourf(times2save,hz,tf1/abs(max(max(max(tf1(:,:,:))))),40,'linecolor','none')
    set(gca,'clim',[-2 2],'ylim',[2 150])
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title([ 'Power via multitaper from channel ' channel2plot ])
end
    channel2plot = 'o1';

    % wavelet parameters
    min_freq =  2;
    max_freq = 150;
    num_frex = 20;


    % baseline time window
    baseline_time = [ -300 -100 ];
    % convert baseline from ms to indices
    baseidx = dsearchn(xplot',baseline_time');


    % other wavelet parameters
    frex = logspace(log10(min_freq),log10(max_freq),num_frex);
    time = -1:1/Fs:1;
    half_wave = (length(time)-1)/2;

    % FFT parameters
    nKern = length(time);
    nData = 2.2*Fs*stinum;
    nConv(1:2) = nKern+nData-1;
    nConv(3)   = nKern+2.2*Fs-1; % ERP is only one trial-length

    % find sensor index
    %chanidx = find(strcmpi(channel2plot,{EEG.chanlocs.labels}));

    % initialize output time-frequency data
    tf = zeros(4,length(frex),2.2*Fs);

    %% prepare non-phase-locked activity

    % compute ERP
    dsb=1;
    erp = squeeze(mean(Rawtrace(dsb,:,:),3));

    % compute non-phase-locked power by subtracting ERP from each trial
    nonphase_EEG = squeeze( bsxfun(@minus,Rawtrace(dsb,:,:),erp) );


    % figure(1), clf
    % subplot(311)
    % plot(xplot,erp)
    % xlabel('Time (ms)'), ylabel('Voltage (\muV)')
    % set(gca,'xlim',[-200 1000])
    % 
    % subplot(312)
    % plot(xplot,Rawtrace(order,:,10))
    % hold on
    % plot(xplot,squeeze(nonphase_EEG(:,10)),'r')
    % legend({'total';'non-phase-locked'})
    % xlabel('Time (ms)'), ylabel('Voltage (\muV)')
    % set(gca,'xlim',[-200 1000])
    % 
    % subplot(313)
    % plot(xplot,erp)
    % hold on
    % plot(xplot,mean(nonphase_EEG,2),'r')
    % legend({'total ERP';'non-phase-locked ERP'})
    % xlabel('Time (ms)'), ylabel('Voltage (\muV)')
    % set(gca,'xlim',[-200 1000])
    % 
    % %% time-frequency decomposition

    % FFT of total data
    fft_EEG{1} = fft( reshape(Rawtrace(dsb,:,:),1,[]), nConv(1));
    % FFT of nonphase-locked data
    fft_EEG{2} = fft( reshape(nonphase_EEG,1,[]), nConv(2));
    % FFT of ERP (phase-locked data)
    fft_EEG{3} = fft( erp ,nConv(3));


for fi=1:length(frex)

    % create wavelet and get its FFT
    s        = 6/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));

    % run convolution for each of total, induced, and evoked
    for methodi=1:3

        % need separate FFT 
        waveletX = fft(wavelet,nConv(methodi));

        % notice that the fft_EEG cell changes on each iteration
        as = ifft(waveletX.*fft_EEG{methodi},nConv(methodi));
        as = as(half_wave+1:end-half_wave);
        if methodi<3
            as = reshape(as,2.2*Fs,stinum);

            % compute power
            temppow = mean(abs(as).^2,2);
        else
            temppow = abs(as).^2;
        end

        % db correct power
        tf(methodi,fi,:) = 10*log10( temppow ./ mean(temppow(baseidx(1):baseidx(2))) );

        % inter-trial phase consistency on total EEG
        if methodi==1
            tf(4,fi,:) = abs(mean(exp(1i*angle(as)),2));
        end
    end % end loop around total, evoked, induced
end % end frequency loop



    analysis_labels = {'Total';'Non-phase-locked';'ERP power';'ITPC'};
    
    % color limits
    clims = [ -3 3; -3 3; -15 15; 0 .6 ];
    
    % next two lines scale the ERP for plotting on top of the TF plot
    erpt = (erp-min(erp))./max(erp-min(erp));
    erpt = erpt*(frex(end)-frex(1))+frex(1);
    
    figure(31)
    colormap jet
    for methodi=1:3
        
        subplot(2,2,methodi)
        %contourf(xplot,frex,squeeze(tf(methodi,:,:)),40,'linecolor','none')
        imagesc(tf(methodi,:,:))
        set(gca,'clim',clims(methodi,:),'xlim',[-200 2000],'xtick',-200:200:2000)
        xlabel('Time (ms)')
        ylabel('Frequency (Hz)')
        title(analysis_labels{methodi})
        
%         % plot ERP on top
%         hold on
%         plot(xplot,erpt,'k')
    end
    
    subplot(224)
    colormap jet
    % estimate phase-locked part of the signal as the difference 
    % between the phase-locked and non-phase-locked
    phaselockedTF = squeeze(tf(1,:,:) - tf(2,:,:));
    % the next line is equivalent to the previous line, just FYI
    % phaselockedTF = squeeze(diff(tf([2 1],:,:),[],1));
    
    contourf(xplot,frex,phaselockedTF,40,'linecolor','none')
    set(gca,'clim',clims(1,:),'xlim',[-200 1000],'xtick',-200:200:1000)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title('Phase-locked')
    eegX = fft( Allsame_plot ,10000,2)/(1.2*Fs);
    ampl = 2*(abs(eegX));

    eegX1 = fft( Saliency_plot,10000,2)/(1.2*Fs);
    amp2 = 2*abs(eegX1);
    hz = linspace(0,Fs,5*Fs);

    figure(32)
    plot(hz,mean(amp2,1),'g','linewidth',2)
    hold on
    plot(hz,mean(ampl,1),'b','linewidth',2)
    colormap hsv
    set(gca,'xlim',[0 150])
    xlabel('Frequency (Hz)'), ylabel('Amplitude (\muV)')
    title( 'Single-trial and trial-average power spectrum from channel ')