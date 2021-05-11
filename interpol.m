clear all;
clc;
%Först få punkterna i I.
% konstanter
C = 5e-7;
L_0 = 0.7;
I_0 = 1;

% Från calc_freq
period = 2*pi*sqrt(C*L_0);

start = 0;
steps = 1000;
stop = 4*period;

maxs = zeros(2, 6);
periods = zeros(1, 6);
mi = 1;
for U_0 = [220,1500,2300]
    for steps = [10000,20000]
        h = (stop -start)/steps;
        
        %find  points to interpolate between
        results = runge_kutta(U_0,start, steps, stop);
        xx = [0:h/100:stop];
        cx = spline([start:h:stop],results(1,:));
        yy = ppval(cx,xx);
        %Hitta max val
        [y,x] = max(yy);
        maxs(:,mi) = [xx(x),y];
        mi = mi +1;
        plot(xx,yy,'-',xx(x),y, 'o');
        hold on;
        
        %räkna ut "roots" (där den korsar x axeln) för att få periodstid
        prev_val = [0,0];
        ii = 1;
        for i = [1:size(xx,2)]
            if prev_val(2) * yy(i) < 0
                root_data(:,ii) = [prev_val(1),xx(i), prev_val(2), yy(i), mean([prev_val(1),xx(i)])];
                ii = ii + 1;
            end
            prev_val = [xx(i),yy(i)];
        end
    end
end

% Topvärden: 220 => 0.18755, 1500 => 1.9971, 2300 => 6.5386
hold off;