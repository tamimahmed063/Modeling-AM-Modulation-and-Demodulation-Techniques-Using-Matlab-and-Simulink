clc;
clear all;
close all;
ts=1.e-4;
t=-0.04:ts:0.04;
Ta=0.01;
m_sig=triplesinc(t,Ta);
figure(1)
plot(t,m_sig)
Lfft=length(t);
Lfft=2^ceil(log2(Lfft))
M_fre=fftshift(fft(m_sig,Lfft));
fregm=(-Lfft/2:Lfft/2-1)/(Lfft*ts)

s_dsb=m_sig.*(2*(cos(2*pi*500*t)));
lfft=length(t);
Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft));
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);



Trange=[-0.03 0.03 -2 2]
figure(2)
subplot(221); td1=plot(t,m_sig);
td1=plot(t,m_sig);
axis(Trange); set(td1,'Linewidth',1);
xlabel('{\it t} (sec)');
ylabel('{\it m}({\it t})')

subplot(223); td2=plot(t,s_dsb);
axis(Trange);
set(td2,'Linewidth',1);
xlabel('{it t} (sec)');
ylabel('{\it s}_{\it s}_{\rm DSB}({\it t})')

Frange =[-600 600 0 200]
subplot(222); fd1=plot(fregm,abs(M_fre));
axis(Frange); set(fd1,'Linewidth',1);
xlabel('{\it f} (Hz)');
ylabel('{\it M}({\it f})')
subplot(224);fd2=plot(freqs,abs(S_dsb));
axis(Frange);set(fd2,'Linewidth',1);
xlabel('{\it f} (Hz)');
ylabel('{\it S}_{rm DSB} ({\it f})')


