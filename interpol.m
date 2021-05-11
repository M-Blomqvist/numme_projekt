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
num_steps = [1000, 2000, 4000, 8000, 16000, 32000, 64000];
stop = 4*period;

maxs = zeros(2, 3*size(num_steps, 2));
periods = zeros(1, 3*size(num_steps, 2));
mi = 1;
for U_0 = [220,1500,2300]
    for steps = num_steps
        h = (stop -start)/steps;
        
        %find  points to interpolate between
        results = runge_kutta(U_0,start, steps, stop);
        xx = [0:h/10:stop];
        cx = spline([start:h:stop],results(1,:));
        yy = ppval(cx,xx);
        %Hitta max val
        [y,x] = max(yy);
        maxs(:,mi) = [xx(x),y];
        
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
        
        % ta diffen mellan "roots", vilket ger T/2
        diffs = zeros(1, size(root_data, 2)-1);
        for i = [2:size(root_data,2)]
            diffs(i) = abs(root_data(5, i)-root_data(5, i-1));
        end
        periods(mi) = mean(diffs)*2;
        
        mi = mi +1;
    end
end
% beräkna diffen för de olika steglängderna
errors = zeros(2, size(periods, 2) - 1);
hs = (stop)*num_steps.^-1;
for i = [1:size(errors,2)]
    errors(1,i) = abs(periods(i+1) - periods(i));
    errors(2,i) = abs(maxs(i+1) - maxs(i));
end
hold off;
pause;
loglog(hs(2:end), errors(1, :));
hold on;
loglog(hs(2:end), errors(2, :));
% Topvärden: 220 => 0.18755, 1500 => 1.9971, 2300 => 6.5386
hold off;

























