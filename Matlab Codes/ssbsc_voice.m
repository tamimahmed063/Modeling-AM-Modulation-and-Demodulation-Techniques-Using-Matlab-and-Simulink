clc;
clear all;
close all;
%% Input
[y, Fs] = audioread('tamimDSB.wav');
y = y';
n = length(y);
t = 0:1/Fs:(n-1)*(1/Fs);
ts=1/Fs;
m_sig=y;
bw=3400;
ts=1/Fs;
Lm_sig=length(m_sig);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft)) ;
M_sig=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
h=fir1(40,[bw*ts]); 

%% SSB modulation
fc=4000; % carrier frequency
s_dsb=(m_sig).*cos(2*pi*fc*t); 
Lfft=length(t); 
Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft)); 
L_lsb=floor(fc*ts*Lfft);
SSBfilt=ones(1,Lfft); 
%% SSBfilt(Lfft/2-L_lsb+1:Lfft/2+L_lsb)=zeros(1,2*L_lsb); 
S_ssb=S_dsb.*SSBfilt;
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts); 
s_ssb=real(ifft(fftshift(S_ssb)));
s_ssb=s_ssb(1:Lm_sig);
%% Channel
s_noise=awgn(100*s_ssb,20);
fr=[2*4500/Fs 2*7500/Fs];
%% Reciever
s_rcv=bandpass(s_noise,fr);


%% Demodulation begins by using a rectifier 
s_dem=s_rcv.*cos(2*pi*fc*t)*2; 
S_dem=fftshift(fft(s_dem,Lfft)); 
%% Using an ideal low pass filter with bandwidth 150 Hz
s_rec=filter(h,1,s_dem);

S_rec=fftshift(fft(s_rec,Lfft)); 
%% Plot 
figure(1)
subplot(221); plot(t,m_sig,'Linewidth',1.5) 
title('message signal')

subplot(222); plot(t,s_ssb,'Linewidth',1.5)
title('SSB-SC modulated signal')

subplot(223); plot(t, s_dem,'Linewidth',1.5)
title('After multiplying local carrier')

subplot(224); plot(t,s_rec,'Linewidth',1.5)
title('Recovered signal')


figure(2)
subplot(221); plot(freqm,abs(M_sig),'Linewidth',1.5)
title('Message Spectrum')

subplot(222); plot(freqs,abs(S_ssb),'Linewidth',1.5)
title('Upper Sideband SSB-SC spectrum')

subplot(223); plot(freqs, abs(S_dem),'Linewidth',1.5) 
title('Detector spectrum')

subplot(224); plot(freqs,abs(S_rec),'Linewidth',1.5)
title('Recovered signal')
%% Output
sound(s_rec,Fs)