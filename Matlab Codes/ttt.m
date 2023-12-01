clc;clear all;close all;
[y, Fs] = audioread('speech_dft.wav');
figure(5)
plot(y)
y = y';
y=y(1,:);
n = length(y);
t = 0:1/Fs:(n-1)*(1/Fs);
ts=1/Fs;
figure(4)
plot(t,y)
x1 = (find(t<=0));
x2 = (find(t>=3.5));
x1 = x1(end);
x2 = x2(1);
y = y(x1:x2);
t = t(x1:x2);
figure(3)
plot(t,y)
z=dsbsc(y,t);
sound(z,Fs)
