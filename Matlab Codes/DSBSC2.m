clc;
clear all;
close all;
ts=1.e-4
t=-0.04:ts:0.04;
Ta=0.01;
m_sig=triangl((t+0.01)/0.01)-triangl((t-0.01)/0.01);
Lm_sig=length(m_sig);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
B_m=150;
%h=fir1(40, [B_m*ts])
h=fir1(40,.015);
t=-0.04:ts:0.04;
Ta=0.01;
fc=300;
s_dsb=m_sig.*cos(2*pi*fc*t);
Lfft=length(t);
Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft));
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
s_dem=s_dsb.*cos(2*pi*fc*t)*2;
S_dem=fftshift(fft(s_dem,Lfft));

s_rec=filter(h,1,s_dem);
S_rec=fftshift(fft(s_rec,Lfft));

Trange=[-0.025 0.025 -2 2];
figure(1)
subplot(221); td1=plot(t,m_sig);
axis(Trange); set(td1,'Linewidth',2);
xlabel('{\it t} (sec)');
ylabel('{\it m}({\it t})')
title('message signal');
subplot(222); td2=plot(t,s_dsb);
axis(Trange);
set(td2,'Linewidth',2);
xlabel('{it t} (sec)');
ylabel('{\it s}_{\it s}_{\rm DSB}({\it t})')
subplot(223); td3=plot(t,s_dem);
axis(Trange);
set(td3,'Linewidth',2);
xlabel('{it t} (sec)');
ylabel('{\it e} ({\it t})')
title('{\it e}({\it t})');
subplot(224); td4=plot(t,s_rec);
axis(Trange);
set(td4,'Linewidth',2);
xlabel('{it t} (sec)');
ylabel('{\it m}_d({\it t})')
title('Recovered signal');

Frange =[-700 700 0 200]
figure(2)
subplot(221); fd1=plot(freqm,abs(M_fre));
axis(Frange); set(fd1,'Linewidth',2);
xlabel('{\it f} (Hz)');
ylabel('{\it M}({\it f})')
title('message spectrum')
subplot(222);fd2=plot(freqs,abs(S_dsb));
axis(Frange);set(fd2,'Linewidth',2);
xlabel('{\it f} (Hz)');
ylabel('{\it S}_{rm DSB} ({\it f})')
title('DSB-SC specttrum')
subplot(223); fd3=plot(freqs,abs(S_dem));
axis(Frange); set(fd3,'Linewidth',2);
xlabel('{\it f} (Hz)');
ylabel('{\it E}({\it f})');
title('spectrum of {\it e}({\it t})')
subplot(224); fd4=plot(freqs,abs(S_rec));
axis(Frange); set(fd4,'Linewidth',2);
xlabel('{\it f} (Hz)');

ylabel('{\it M}_d({\it f})');
title('recovered spectrum')
