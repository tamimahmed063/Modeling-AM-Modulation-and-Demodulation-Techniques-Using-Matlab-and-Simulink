function m= triplesinc(t,Ta)

sig_1=sinc(2*t/Ta);
sig_2=sinc(2*t/Ta-1);
sig_3=sinc(2*t/Ta+1);
m=2*sig_1+sig_2+sig_3;
end