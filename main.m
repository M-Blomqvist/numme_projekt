%main.m
clc
clear all;
%Ge uppskattad frekvens & period med konstant L_0
[period,freq] = calc_freq();
pause

%RungeKutta 4 för att ge strömkurvorna I(t) med antalet steg (steps) och
%u-värden U_0s.
hs = [1e-6,5e-7,25e-8,125e-9, 625e-10]; 
start = 0;
periods = [period,1.5*period, 2*period];
U_0s = [220,1500, 2300];
for stop = periods
   for h = hs
      for U_0 = U_0s
          results = runge_kutta(U_0, start, h, stop);
          txt = ['I(t) med U_0 = ',num2str(U_0), ' & h = ', num2str(h)];
          plot([start:h:stop], results(1,:),'DisplayName',txt);
          hold on;
      end
   end
      title(['RungeKutta4 I(t)']);
      %legend('Location','southwest');
      pause;
      hold off;
end

%Visa att E(t) är konstant:
for stop = periods
    for U_0 = U_0s
        for h = hs
          Es = E_const(U_0, start, h, stop);
          txt = ['E(t) med h = ',num2str(h)];
          plot([start:h:stop], Es,'DisplayName',txt);
          hold on;
        end
        %move it here
    end
    %Move this up for more accurate yet silly plots
    title(['RungeKutta4 E(t) med U_0 = ', num2str(U_0)]);
    %legend('Location','southwest');
    pause;
    hold off;
end

%Använd interpolation för att finna I max & perioden T
interpol_steps_per_step = 100;

for stop = periods(2)
   [hs, inter_period, period_errs,inter_max,  max_errs] = interpol_errors(U_0s, start, hs, stop);
end

u_i = 1;
for U_0 = U_0s
    txt = ['U_0 = ',num2str(U_0)];
    plot(hs,period_errs(u_i,:),'DisplayName', txt);
    title('Felvärden interpolation: perioder');
    u_i = u_i +1;
    hold on;
end
pause
hold off;
u_i = 1;
for U_0 = U_0s
    txt = ['U_0 = ',num2str(U_0)];
    plot(hs,max_errs(u_i,:),'DisplayName',txt);
    title('Felvärden interpolation: max värden');
    u_i = u_i +1;
    hold on;
end
pause
hold off;
%%
%Best candidate: 40000 steps = no period errors, max_errors <10^-11
%Hitta U_0* så att 1/period = 400 = 0.0025, U_0 = 220 best candidate
% wanted_period = 400^-1;
% stop = periods(2);
% h = hs(4);
% guess = 220;
% U_h = 1e-9;
% diff = 100;
% while diff > 1e-10 
%     fprintf('\n guesses = %d, h %d, diff = %d\n', guess,U_h,diff);
%     [~,~,~,y_max,i_period] = interpol(guess, start, h, stop,interpol_steps_per_step);
%     guess_h = guess +U_h;
%     [~,~,~,~,i_period_h] = interpol(guess_h, start, h, stop, interpol_steps_per_step);
%     fprintf('\n i_period = %d, i_period_h = %d\n', i_period,i_period_h);
%     guess = guess-U_h*(i_period/(i_period_h-i_period));
%     diff = abs(i_period - wanted_period);
% end
% U_star = guess;
% max_star = y_max;