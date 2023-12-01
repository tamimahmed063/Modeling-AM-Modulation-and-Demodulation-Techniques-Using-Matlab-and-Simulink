clc;
clear all;
close all;
[m_sig1 Fs]= audioread('tamimDSB.wav');
[m_sig2 Fs]= audioread('tamim.wav');
nBits=16;
BW=3000;
fc=5000;
% [s_rec] = qam_mod_demod(nBits,Fs,m_sig1,BW,fc,m_sig2);
ts = 1/Fs;
t = linspace(0,length(m_sig1)/Fs,length(m_sig1));

Lm_sig = length(m_sig1);
Lfft = length(t);
Lfft = 2^ceil(log2(Lfft));
M1_fre = fftshift(fft(m_sig1,Lfft));
M2_fre = fftshift(fft(m_sig2,Lfft));
freqm = (-Lfft/2:Lfft/2 -1)/(Lfft*ts);

h = fir1(40,[BW*ts]);

%% QAM modulation

% QAM signal generated adding a carrier to DSB-SC
s_qam = m_sig1'.*cos(2*pi*fc*t)+m_sig2'.*sin(2*pi*fc*t);
Lfft = length(t);  Lfft = 2^ceil(log2(Lfft)+1);
S_qam = fftshift(fft(s_qam,Lfft));
freqs = (-Lfft/2:Lfft/2 -1)/(Lfft*ts);

%% QAM demodulation

% demodulation begins by using a rectifier
s_dem1 = s_qam.*cos(2*pi*fc*t)*2;
S_dem1 = fftshift(fft(s_dem1,Lfft));
% demodulate 2nd signal
s_dem2 = s_qam.*sin(2*pi*fc*t)*2;
S_dem2 = fftshift(fft(s_dem2,Lfft));

% using ideal LPF with BW
s_rec1 = filter(h,1,s_dem1);
S_rec1 = fftshift(fft(s_rec1,Lfft));
s_rec2 = filter(h,1,s_dem2);
S_rec2 = fftshift(fft(s_rec2,Lfft));
Trange = [-0.025 0.025 -2 2];
Trange2 = [-0.025 0.025 -2 4];

figure(1)
subplot(221); plot(t,m_sig1,'Linewidth',1.5);  axis(Trange);
xlabel('t,(sec)');  ylabel('m1(t)');  title('message signal 1');  
subplot(222);  plot(t,s_qam,'Linewidth',1.5);  axis(Trange);
xlabel('t,(sec)');  ylabel('s1_DSB(t)');  title('QAM modulated signal');
subplot(223);  plot(t,s_dem1,'Linewidth',1.5);  axis(Trange2);
xlabel('t,(sec)');  ylabel('e1(t)');  title('first demodulator output');
subplot(224);  plot(t,s_rec1,'Linewidth',1.5);  axis(Trange);
xlabel('t,(sec)');  ylabel('m_d1(t)');  title('detected signal 1');

figure(2)
subplot(221);  plot(t,m_sig2,'Linewidth',1.5);  axis(Trange);
xlabel('t(sec)');  ylabel('m2(t)');  title('message signal 2');  
subplot(222);  plot(t,s_qam,'Linewidth',1.5);  axis(Trange);
xlabel('t(sec)');  ylabel('s2_DSB(t)');  title('QAM modulated signal');
subplot(223);  plot(t,s_dem2,'Linewidth',1.5);  axis(Trange2);
xlabel('t(sec)');  ylabel('e2(t)');  title('second demodulator output');
subplot(224); plot(t,s_rec2,'Linewidth',1.5);  axis(Trange);
xlabel('t(sec)');  ylabel('m_d2(t)');  title('detected signal 2');

Frange = [-700 700 0 250];
figure(3)
subplot(221);  plot(freqm,abs(M1_fre),'Linewidth',1.5);  axis(Frange);
xlabel('f(Hz)');  ylabel('M1(f)');  title('message 1 spectrum');  
subplot(222);  plot(freqs,abs(S_qam),'Linewidth',1.5);  axis(Frange);
xlabel('f(Hz)');  ylabel('S1_AM(f)');  title('QAM spectrum magnitude');
subplot(223);  plot(freqs,abs(S_dem1),'Linewidth',1.5);  axis(Frange);
xlabel('f(Hz)');  ylabel('E1(f)');  title('first demodulator spectrum');
subplot(224);  plot(freqs,abs(S_rec1),'Linewidth',1.5);  axis(Frange);
xlabel('f(Hz)');  ylabel('M_d1(f)');  title('recovered spectrum 1');

figure(4)
subplot(221);  plot(freqm,abs(M2_fre),'Linewidth',1.5);  axis(Frange);
xlabel('f(Hz)');  ylabel('M2(f)');  title('message 2 spectrum');  
subplot(222);  plot(freqs,abs(S_qam),'Linewidth',1.5);  axis(Frange);
xlabel('f(Hz)');  ylabel('S2_AM(f)');  title('QAM spectrum magnitude');
subplot(223);  plot(freqs,abs(S_dem2),'Linewidth',1.5);  axis(Frange);
xlabel('f(Hz)');  ylabel('E2(f)');  title('second demodulator spectrum');
subplot(224);  plot(freqs,abs(S_rec2),'Linewidth',1.5);  axis(Frange);
xlabel('f(Hz)');  ylabel('M_d2(f)');  title('recovered spectrum 2');