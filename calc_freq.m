function [period,freq] = calc_freq()
%Analytiskt funnen period för L = L_0, samt plot för denna
C = 5e-7;
L = 0.7;
const = 1/sqrt(C*L);
period = 2*pi*sqrt(C*L);
freq = 1/period;
I = @(t) cos(t*const) + sin(t*const);
fplot(I, [-period,period]);
title('I med L=L_0');
hold on 
plot([-period,0,period], [0,0,0], '*'); 
legend({'I(t)','punkter med 1 periods distans'},'Location','southwest')
hold off
end