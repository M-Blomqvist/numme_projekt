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
% Interpolate to find max_i and period
certainty = 1e-12;
for stop = periods(2)
   [hs, inter_period, period_errs,inter_max,  max_errs] = interpol_errors(U_0s, start, hs, stop,certainty);
end

u_i = 1;
for U_0 = U_0s
    txt = ['U_0 = ',num2str(U_0)];
    loglog(hs,period_errs(u_i,:),'DisplayName', txt);
    hold on;
    loglog(hs,hs.^2);
    title('Felvärden interpolation: perioder');
    u_i = u_i +1;
end
pause
hold off;
u_i = 1;
for U_0 = U_0s
    txt = ['U_0 = ',num2str(U_0)];
    loglog(hs,max_errs(u_i,:),'DisplayName',txt);
    hold on;
    loglog(hs,hs.^2);
    title('Felvärden interpolation: max värden');
    u_i = u_i +1;
end
pause
%%
hold off;
%Hitta U_0* så att 1/period = 400 = 0.0025, U_0 = 220 & 1500 best candidate
wanted_period = 400^-1;
stop = periods(2);
guess1 = 220;
guess2 = 1500;
certainty = 1e-14;
diff = 100;
U_stars = zeros(1,size(hs,2));

for i = [1:size(hs,2)]
    h = hs(i);
    f = @(u) f_period(u,start,h,stop,certainty, wanted_period);
    U_stars(i) = sekant(f, guess1,guess2,certainty);
    fprintf('%d', U_stars(i));
end
%use most accurrate h
h = hs(end);

inter_periods = zeros(1,size(U_stars,2));
inter_maxs = zeros(1,size(U_stars,2));
u_i = 1;
for u = U_stars
    [cx,max_x,max_y,period] = interpol(u,start,h, stop, certainty);
    fprintf('\n %d \n', period);
    inter_periods(u_i) = period;
    inter_maxs(u_i) = max_y;
    txt = ['I(t) med interpolerade maxvärden, period = ', num2str(period)];
    x = [start:h:stop];
    plot(x, ppval(cx,x),'DisplayName',txt);
    hold on;
    plot(max_x,max_y, 'o', 'DisplayName',['max värde för ovan']);
    u_i = u_i +1;
end

pause;
hold off;

period_errs = zeros(1,size(U_stars,2));
max_errs = zeros(1,size(U_stars,2));
for u_i = 2:size(U_stars,2)
   max_errs(u_i) = abs(inter_maxs(u_i) - inter_maxs(u_i-1));
   period_errs(u_i) = abs(inter_periods(u_i) - inter_periods(u_i-1));
end  

%%
for i = 2:size(U_stars,2)
   diffs(i-1) = abs(U_stars(i)-U_stars(i-1)); 
end

%%   Utvidgning: v(t)
h = hs(3);
u_i = 1;

% Plotta den omformade v(t) för de tre U_0
for u = U_0s
    for h_i = 1:size(hs,2)
        h = hs(h_i);
        [cx,~,y_max,period] = interpol(u,start,h,stop,certainty);
        vs{u_i,h_i} = @(xs) ppval(cx,xs*period/(2*pi))/y_max;
    end
    xx = [0:h*2*pi/period:2*pi];
    txt = ['V(t) med U_0 = ', num2str(u)];
    vv = vs{u_i,h_i};
    plot(xx,vv(xx),'DisplayName',txt);
    hold on;
    u_i = u_i + 1;
end
pause
hold off;
% Beräkna och plotta fourierutvecklingen av v(t) för de tre U_0
coefs =10; 
for u = 1:size(U_0s,2)
    for h_i = 1:size(hs,2)
        v = vs{u,h_i};
        h = hs(h_i);
        a = @(k,xs) v(xs).*sin(k.*xs);
        xx = [0:h*2*pi/period:2*pi];
        for k = 1:coefs
            as(u,h_i,k) = (1/pi)*trapz(xx, a(k,xx)); % TODO: beräkna för h/2 och jämför för att visa 4 siffrors nogrannhet
        end
    end
end

for u = 1:size(U_0s,2)
    for h_i = 2:size(hs,2)
        for k = 1:coefs
            diff_a(h_i-1,k) = abs(as(u,h_i,k) - as(u,h_i-1,k));
        end
    end
    for k = 1:coefs
        plot(hs(2:end), diff_a(:,k));
        hold on;
    end
    pause
    hold off;
end

%choose h to continue with
h_i = 3;
for u = 1:size(U_0s,2)
    % definera fourierutvecklingarna med beräknade koefficienter
    fourier_funcs{u,1} = @(t) as(u,h_i,1)*sin(t) + as(u,h_i,2)*sin(2*t) + as(u,h_i,3)*sin(3*t);
    fourier_funcs{u,2} = @(t) as(u,h_i,1)*sin(t) + as(u,h_i,2)*sin(2*t) + as(u,h_i,3)*sin(3*t) + as(u,h_i,4)*sin(4*t) + as(u,h_i,5)*sin(5*t) + as(u,h_i,6)*sin(6*t) + as(u,h_i,7)*sin(7*t) + as(u,h_i,8)*sin(8*t) + as(u,h_i,9)*sin(9*t) + as(u,h_i,10)*sin(10*t);
    
    plot(xx, fourier_funcs{u,1}(xx),xx,fourier_funcs{u,2}(xx),xx, vs{u,h_i}(xx));
    legend({'Fourierutveckling med 3 termer','Fourierutveckling med 10 termer',['v(t) för U_0 = ',num2str(U_0s(u))]},'Location','southwest')
    pause
end
% Plotta de udda a-värderna för alla U_0
for u = 1:size(U_0s,2)
    semilogy(abs(squeeze(as(u,h_i,1:2:end)))); 
    title(['semilogy-plot för udda a_k, U_0 = ',num2str(U_0s(u))])
    pause
end

%% sound-part
    h = 0.00001;
for u_i = 1:size(U_0s,2)
    u = U_0s(u_i);
    [cx,~,y_max,period] = interpol(u,start,h,stop,certainty);
    v = @(xs) ppval(cx,xs*period/(2*pi))/y_max;
    xx = [0:h*2*pi/period:2*pi];
    plot(xx, v(xx));
    for i = 0:1:399
        y(i*size(xx,2)+1:(i+1)*size(xx,2)) = v(xx);
    end
    Fs = 400*size(xx,2);
    %sound(y,Fs);
    audiowrite(['audio_',num2str(u),'.ogg'],y,Fs);
    pause;
end

%% Function declaration for U_0* calculations
function period = f_period(u,start,h,stop,certainty,wanted_period)
    [~,~,~,p] = interpol(u, start, h, stop,certainty);
    period = p - wanted_period;
end