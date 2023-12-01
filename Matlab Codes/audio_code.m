clc;
clear all;
close all;
T = 3;

Fs = 16000;
bits = 16;
channel = 1;

recObj = audiorecorder(Fs,bits,channel,1);

disp('Start speaking...')
recordblocking(recObj, T);
disp('End of Recording.');
play(recObj);

y = getaudiodata(recObj);

audiowrite(sprintf('tamimDSB.wav'),y, Fs);
