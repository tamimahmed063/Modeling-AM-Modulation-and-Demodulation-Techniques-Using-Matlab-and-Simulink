clc;clear all;close all;
[y, Fs] = audioread('speech_dft.wav');
figure(5)
plot(y)
y = y';
n = length(y);
t = 0:1/Fs:(n-1)*(1/Fs);
ts=1/Fs;
figure(4)
plot(t,y)
x1 = (find(t<=0));
x2 = (find(t>=4));
x1 = x1(end);
x2 = x2(1);
y = y(x1:x2);
t = t(x1:x2);
figure(3)
plot(t,y)
m_sig=y;
ts=1.e-4;
Lm_sig=length(m_sig);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
%B_m=150;
h=fir1(40, [ts]);

% t=-0.04:ts:0.04;
Ta=0.01;
fc=2000;
s_dsb=m_sig.*cos(2*pi*fc*t);

s_ch=awgn(100*s_dsb,-5);
s_chhh=s_ch;
fr=[2*2300/Fs 2*5400/Fs];

%%
%s_recieved=xcorr(s_ch);


s_rcv=bandpass(s_chhh,fr);


%%
Lfft=length(t);
Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft));
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

t1=t;
%t=linspace(t(1),t(end),length(s_recieved));
%s_dem=s_recieved.*cos(2*pi*fc*t)*2;
%s_dem=s_dsb.*cos(2*pi*fc*t)*2;
%s_dem=s_rcv.*cos(2*pi*fc*t);
s_dem=s_rcv.*cos(2*pi*fc*t);
S_dem=fftshift(fft(s_dem,Lfft));

s_rec=filter(h,1,s_dem);
S_rec=fftshift(fft(s_rec,Lfft));

% Trange=[-0.025 0.025 -2 2];
figure(1)
subplot(221); td1=plot(t1,m_sig);
% axis(Trange); set(td1,'Linewidth',2);
xlabel('{\it t} (sec)');
ylabel('{\it m}({\it t})')
title('message signal');
subplot(222); td2=plot(t1,s_dsb);
% axis(Trange);
% set(td2,'Linewidth',2);
xlabel('{it t} (sec)');
ylabel('{\it s}_{\it s}_{\rm DSB}({\it t})')
subplot(223); td3=plot(t,s_dem);
% axis(Trange);
% set(td3,'Linewidth',2);
xlabel('{it t} (sec)');
ylabel('{\it e} ({\it t})')
title('{\it e}({\it t})');
subplot(224); td4=plot(t,s_rec);
% axis(Trange);
% set(td4,'Linewidth',2);
xlabel('{it t} (sec)');
ylabel('{\it m}_d({\it t})')
title('Recovered signal');

Frange =[-700 700 0 200]
figure(2)
subplot(221); fd1=plot(freqm,abs(M_fre));
% axis(Frange); set(fd1,'Linewidth',2);
xlabel('{\it f} (Hz)');
ylabel('{\it M}({\it f})')
title('message spectrum')
subplot(222);fd2=plot(freqs,abs(S_dsb));
% axis(Frange);set(fd2,'Linewidth',2);
xlabel('{\it f} (Hz)');
ylabel('{\it S}_{rm DSB} ({\it f})')
title('DSB-SC specttrum')
subplot(223); fd3=plot(freqs,abs(S_dem));
% axis(Frange); set(fd3,'Linewidth',2);
xlabel('{\it f} (Hz)');
ylabel('{\it E}({\it f})');
title('spectrum of {\it e}({\it t})')
subplot(224); fd4=plot(freqs,abs(S_rec));
% axis(Frange); set(fd4,'Linewidth',2);
xlabel('{\it f} (Hz)');
ylabel('{\it M}_d({\it f})');
title('recovered spectrum')

sound(s_rec,Fs)