clc;
close all;
clear all;

%% generation of message signal and carrier signal

[y_tr, Fs] = audioread('speech_dft.wav');
nBits = 16;
% y_tr = y(0.5*Fs: 1*Fs);
t = linspace(0,length(y_tr)/Fs,length(y_tr));
% figure(1)
% plot(t,y_tr),title 'truncated audio file'
% sound(y_tr,Fs,nBits)


% Fs =16000;
ts=1/Fs;
t = linspace(0,length(y_tr)/Fs,length(y_tr));
BW = 5000;
fc = 5000;

% % t=-0.04:ts:0.04;
% Ta=0.01;
m_sig=y_tr;%triangl((t+0.01)/0.01)-triangl((t-0.01)/0.01);
Lm_sig=length(m_sig);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft)); 
M_sig=fftshift(fft(m_sig,Lfft)); 

freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts) ;
% BW=150;

h=fir1(40,[BW*ts]);

% 
% figure(1);
% subplot(211),plot(t,m_sig)
% % title 'Double sideband'
% % xlim([-1000 1000])
% subplot(212),plot(freqm,abs(M_sig))
% % title 'vestigial sideband'
% % xlim([-1000 1000])

%% VSB Modulation
% fc=300; % carrier frequency
s_dsb=(m_sig)'.*cos(2*pi*fc*t); %dimension milanor jonno signal transpose korsi
Lfft=length(t); 
Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft)); 
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

percentageOfLsb = 0.25;
fv = percentageOfLsb*BW/2;
L_vsb=floor((fc-fv)*ts*Lfft);
hfv = floor((fc+fv)*ts*Lfft);

%setting up the vsb filter
VSBfilt=ones(1,Lfft); 

p = Lfft/2+L_vsb:Lfft/2+hfv;
q = Lfft/2-hfv+1:Lfft/2-L_vsb+1;

x = linspace(fc-fv,fc+fv,length(p));
m = 0.5/(2*fv);
y = m.*x - m*(fc-fv);

VSBfilt(p)=y;
VSBfilt(q)=fliplr(y);

VSBfilt(Lfft/2-L_vsb+1:Lfft/2+L_vsb)=zeros(1,2*L_vsb); 

% figure(1);
% subplot(211),plot(freqs,fftshift(abs(VSBfilt)))
% % title 'Double sideband'
% % xlim([-1000 1000])
% subplot(212),plot(freqs,abs(S_dsb))
% % title 'vestigial sideband'
% % xlim([-1000 1000])

S_vsb=S_dsb.*VSBfilt;
% freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
s_vsb=real(ifft(fftshift(S_vsb)));
s_vsb=s_vsb(1:Lm_sig);


s_noise=awgn(100*s_vsb,-3);
fr=[2*5300/Fs 2*8400/Fs];
s_rcv=bandpass(s_noise,fr);

% figure(1);
% subplot(211),plot(freqs,abs(S_dsb))
% title 'Double sideband'
% xlim([-1000 1000])
% subplot(212),plot(freqs,abs(S_vsb))
% title 'vestigial sideband'
% xlim([-1000 1000])

%% demodulation
%%% Demodulation begins by using a rectifier 
s_dem=s_rcv.*cos(2*pi*fc*t)*2; 
S_dem=fftshift(fft(s_dem,Lfft)); 
% Using an ideal low pass filter with bandwidth 150 Hz
s_rec=filter(h,1,s_dem);%ei filter ta niye pech lagbe

%need to design only this filter
S_rec=fftshift(fft(s_rec,Lfft)); 
 


% Trange=[-0.025 0.025 -1 1] ;
% Frange=[-3*fc 3*fc 0 200] ;%-700 700


figure(2)
subplot(221); plot(t,m_sig,'Linewidth',1.5)
% axis(Trange) 
title('message signal')

subplot(222); plot(t,s_vsb,'Linewidth',1.5)
% axis(Trange)
title('VSB-SC modulated signal')

subplot(223); plot(t, s_dem,'Linewidth',1.5)
% axis(Trange) 
title('After multiplying local carrier')

subplot(224); plot(t,s_rec,'Linewidth',1.5)
% axis(Trange)
title('Recovered signal')

figure(3)
subplot(221); plot(freqm,abs(M_sig),'Linewidth',1.5)
% axis(Frange)
title('Message Spectrum')

subplot(222); plot(freqs,abs(S_vsb),'Linewidth',1.5)
% axis(Frange)
title('Upper Sideband VSB-SC spectrum')

subplot(223); plot(freqs, abs(S_dem),'Linewidth',1.5)
% axis(Frange) 
title('Detector spectrum')

subplot(224); plot(freqs,abs(S_rec),'Linewidth',1.5)
% axis(Frange)
title('Recovered signal')

%sound(y_tr,Fs,nBits)

sound(s_rec,Fs,nBits)