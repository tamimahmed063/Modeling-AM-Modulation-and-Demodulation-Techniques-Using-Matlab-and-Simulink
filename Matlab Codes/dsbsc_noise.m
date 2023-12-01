clc;
clear all;
close all;
%% Input
[y, Fs] = audioread('tamim.wav');
y = y';
y=y(1,:);
n = length(y);
t = 0:1/Fs:(n-1)*(1/Fs);
m_sig=y;
bw=3400;
ts=1/Fs;
Lm_sig=length(m_sig);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
h=fir1(40, [bw*ts]);
%% Modulation
fc=4000;
s_dsb=m_sig.*cos(2*pi*fc*t);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft));
%% Channel
s_ch=awgn(100*s_dsb,20);

%% SIgnal recieve
fr=[2*4300/Fs 2*7400/Fs];
s_rcv=bandpass(s_ch,fr);



%% Demodulation
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
s_dem=s_rcv.*cos(2*pi*fc*t);
S_dem=fftshift(fft(s_dem,Lfft));
s_rec=filter(h,1,s_dem);
S_rec=fftshift(fft(s_rec,Lfft));

%% Plot
figure(1)
subplot(221); td1=plot(t,m_sig);
xlabel('{\it t} (sec)');
ylabel('{\it m}({\it t})')
title('message signal');

subplot(222); td2=plot(t,s_dsb);
xlabel('{it t} (sec)');
ylabel('{\it s}_{\it s}_{\rm DSB}({\it t})')

subplot(223); td3=plot(t,s_dem);
xlabel('{it t} (sec)');
ylabel('{\it e} ({\it t})')
title('{\it e}({\it t})');

subplot(224); td4=plot(t,s_rec);
xlabel('{it t} (sec)');
ylabel('{\it m}_d({\it t})')
title('Recovered signal');


figure(2)
subplot(221); fd1=plot(freqm,abs(M_fre));
xlabel('{\it f} (Hz)');
ylabel('{\it M}({\it f})')
title('message spectrum')

subplot(222);fd2=plot(freqs,abs(S_dsb));
xlabel('{\it f} (Hz)');
ylabel('{\it S}_{rm DSB} ({\it f})')
title('DSB-SC specttrum')

subplot(223); fd3=plot(freqs,abs(S_dem));
xlabel('{\it f} (Hz)');
ylabel('{\it E}({\it f})');
title('spectrum of {\it e}({\it t})')

subplot(224); fd4=plot(freqs,abs(S_rec));
xlabel('{\it f} (Hz)');
ylabel('{\it M}_d({\it f})');
title('recovered spectrum')
%% Output
sound(s_rec,Fs)
audiowrite(sprintf('reconstructedVoice.wav'),s_rec, Fs);