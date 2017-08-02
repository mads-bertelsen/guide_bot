% testing a fomula!
clear all;close all;clc;

b=0.03;
t=0.003;
o=0.06;
L1=20:130;
L2=0.8;
L3=150-L1;

phi=atand((b/2+t+o/2)*(L1+L2+L3)./(L3+L2)./L1-t./L1);

plot(L1,phi,'b')
xlabel('L1 / 150-L3')
ylabel('Kink angle in degrees')

hold on

o=0.05;
phi=atand((b/2+t+o/2)*(L1+L2+L3)./(L3+L2)./L1-t./L1);

plot(L1,phi,'r')


o=0.04;
phi=atand((b/2+t+o/2)*(L1+L2+L3)./(L3+L2)./L1-t./L1);

plot(L1,phi,'g')