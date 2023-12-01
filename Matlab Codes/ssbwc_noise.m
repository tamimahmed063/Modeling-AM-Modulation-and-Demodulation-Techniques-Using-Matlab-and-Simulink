clear all;
clc;

ts=10^(-4);
Fs=1/ts;
t=-0.04:ts:0.04;
Ta=0.01;
m_sig=triangl((t+0.01)/0.01)-triangl((t-0.01)/0.01);
Lm_sig=length(m_sig);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft)) 
M_sig=fftshift(fft(m_sig,Lfft))

freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts) 
B_m=150;
h=fir1(40,[B_m*ts]) 

%% SSB modulation
fc=300; % carrier frequency
Ac=1;
carrier=Ac*cos(2*pi*fc*t);
s_dsb=carrier+(m_sig).*cos(2*pi*fc*t) 
Lfft=length(t); 
Lfft=2^ceil(log2(Lfft)+1)
S_dsb=fftshift(fft(s_dsb,Lfft)); 
L_lsb=floor(fc*ts*Lfft)
SSBfilt=ones(1,Lfft); 
%% SSBfilt(Lfft/2-L_lsb+1:Lfft/2+L_lsb)=zeros(1,2*L_lsb); 
S_ssb=S_dsb.*SSBfilt;
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts) 
s_ssb=real(ifft(fftshift(S_ssb)));
s_ssb=s_ssb(1:Lm_sig);
s_noise=awgn(100*s_ssb,-3);
fr=[2*100/Fs 2*500/Fs];
s_rcv=bandpass(s_noise,fr);


%% Demodulation begins by using a rectifier 
s_dem=s_rcv.*cos(2*pi*fc*t)*2 -Ac*100*ones(1,length(s_rcv));
S_dem=fftshift(fft(s_dem,Lfft)); 
% Using an ideal low pass filter with bandwidth 150 Hz
s_rec=filter(h,1,s_dem)
S_rec=fftshift(fft(s_rec,Lfft)); 
 


% Trange=[-0.025 0.025 -1 1] 
% Frange=[-700 700 0 200] 
figure(1)
subplot(221); plot(t,m_sig,'Linewidth',1.5)
% axis(Trange) 
title('message signal')

subplot(222); plot(t,s_ssb,'Linewidth',1.5)
% axis(Trange)
title('SSB-SC modulated signal')

subplot(223); plot(t, s_dem,'Linewidth',1.5)
% axis(Trange) 
title('After multiplying local carrier')

subplot(224); plot(t,s_rec,'Linewidth',1.5)
% axis(Trange)
title('Recovered signal')


figure(2)
subplot(221); plot(freqm,abs(M_sig),'Linewidth',1.5)
% axis(Frange)
title('Message Spectrum')

subplot(222); plot(freqs,abs(S_ssb),'Linewidth',1.5)
% axis(Frange)
title('Upper Sideband SSB-SC spectrum')

subplot(223); plot(freqs, abs(S_dem),'Linewidth',1.5)
% axis(Frange) 
title('Detector spectrum')

subplot(224); plot(freqs,abs(S_rec),'Linewidth',1.5)
% axis(Frange)
title('Recovered signal')