clc;
close all;
clear all;

%% generation of message signal and carrier signal

[y_tr, Fs] = audioread('tamimDSB.wav');
nBits = 16;
t = linspace(0,length(y_tr)/Fs,length(y_tr));
ts=1/Fs;
t = linspace(0,length(y_tr)/Fs,length(y_tr));
bw = 3400;
fc = 4000;
m_sig=y_tr;
Lm_sig=length(m_sig);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft)); 
M_sig=fftshift(fft(m_sig,Lfft)); 
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts) ;
h=fir1(40,[bw*ts]);

%% VSB Modulation
s_dsb=(m_sig)'.*cos(2*pi*fc*t); %Transposed the signal to match dimension
Lfft=length(t); 
Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft)); 
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
percentageOfLsb = 0.25;
fv = percentageOfLsb*bw/2;
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

S_vsb=S_dsb.*VSBfilt;
% freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
s_vsb=real(ifft(fftshift(S_vsb)));
s_vsb=s_vsb(1:Lm_sig);

%% Channel
s_noise=awgn(100*s_vsb,20);
fr=[2*4300/Fs 2*7400/Fs];
%% Reciever
s_rcv=bandpass(s_noise,fr);


%% demodulation
%%% Demodulation begins by using a rectifier 
s_dem=s_rcv.*cos(2*pi*fc*t)*2; 
S_dem=fftshift(fft(s_dem,Lfft)); 
% Using an ideal low pass filter with bandwidth 150 Hz
s_rec=filter(h,1,s_dem);%ei filter ta niye pech lagbe

%need to design only this filter
S_rec=fftshift(fft(s_rec,Lfft)); 
 
%% Plot
figure(2)
subplot(221); plot(t,m_sig,'Linewidth',1.5)
title('message signal')

subplot(222); plot(t,s_vsb,'Linewidth',1.5)
title('VSB-SC modulated signal')

subplot(223); plot(t, s_dem,'Linewidth',1.5)
title('After multiplying local carrier')

subplot(224); plot(t,s_rec,'Linewidth',1.5)
title('Recovered signal')

figure(3)
subplot(221); plot(freqm,abs(M_sig),'Linewidth',1.5)
title('Message Spectrum')

subplot(222); plot(freqs,abs(S_vsb),'Linewidth',1.5)
title('Upper Sideband VSB-SC spectrum')

subplot(223); plot(freqs, abs(S_dem),'Linewidth',1.5)
title('Detector spectrum')

subplot(224); plot(freqs,abs(S_rec),'Linewidth',1.5)
title('Recovered spectrum')
%% Output
%sound(s_rec,Fs,nBits)
y = getaudiodata(s_rec);

audiowrite(sprintf('reconstructedvoice.wav'),y, Fs);