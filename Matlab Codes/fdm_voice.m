clc;
clear all;
close all;


[y1 Fs1]=audioread('tamim01.wav');
[y2 Fs2]=audioread('tamim02.wav');

fpass=2000;
ssbm1=ssbmod(y1, 3000, 8000);
ssbm2=ssbmod(y2, 3000, 8000);
ssbmodu=ssbm1+ssbm2;
t = linspace(0,length(y1)/Fs1,length(y1));

ssbmod2=ssbmod(ssbmodu, 5000,8000 )
ssbdemod1=bandpass(ssbmod2,[6000 8000],17000);
ssbdemod2=bandpass(ssbmod2,[4000 5900],17000);

r1=ssbdemod1 .* sin(2*pi*3000*t);
r2=ssbdemod2 .* sin(2*pi*3000*t);
recons1= lowpass(r1, fpass, Fs1);
recons2= lowpass(r2, fpass, Fs1);
audiowrite(sprintf('reconstructed1.wav'),recons1, Fs);
audiowrite(sprintf('reconstructed2.wav'),recons2, Fs);
