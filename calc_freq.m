clc
clear all
C = 5e-10;
L = 0.7;
const = 1/sqrt(C*L);
period = 2*pi*sqrt(C*L);
freq = 1/period;
I = @(t) cos(t*const) + sin(t*const);
fplot(I, [-0.0002,0.0002]);
hold on 
plot([-period,0,period], [0,0,0], '*'); 