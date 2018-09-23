clear
[data, Fs] = audioread("Woodcreeper.wav");
data = data(1:3*Fs);
xplot = linspace(0,3,3*Fs);
%% wavelet parameters
min_freq = 1;
max_freq = 15000;
num_frex = 50;

% other wavelet parameters
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1:1/Fs:1;
half_wave = (length(time)-1)/2;

% FFT parameters
nKern = length(time);
nData = 3*Fs;
nConv = nKern+nData-1;


% initialize output time-frequency data
tf = zeros(length(frex), 3*Fs);

fft_EEG = fft(reshape(data,1,[]), nConv);

%%
for fi=1:length(frex)
    
    % create wavelet and get its FFT
    s        = 6/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    waveletX = fft(wavelet,nConv);
    % notice that the fft_EEG cell changes on each iteration
    as = ifft(waveletX.*fft_EEG,nConv);
    as = as(half_wave+1:end-half_wave);
    temppow = mean(abs(as).^2,1);
    
    tf(fi,:) = 10*log10( temppow ./ mean(temppow(1:1000),2));
    
end

%%
figure
contourf(xplot,frex,tf,'linecolor','none')