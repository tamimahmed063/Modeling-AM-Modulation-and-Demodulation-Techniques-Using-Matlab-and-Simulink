clc; clear all; close all;
n=0:49;
x=2*n.*0.9.^n;
d=rand(1,50)-0.5;
y=x+d;
M=3;
z=zeros(1,length(n));
for i=0:M-1
shift=zeros(1,length(n));
shift(i+1:length(shift))=y(1:length(y)-i)
z=z+shift
end
z=z/M;
subplot(121), plot(n,x,'k.-')
hold on
plot(n,y)
xlabel('time index')
legend ('x(n)','x(n)+d(n)')
subplot(122), plot(n,x,'k.-')
hold on
plot(n,z)
xlabel('time index')