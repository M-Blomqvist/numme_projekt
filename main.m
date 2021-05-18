%main.m
clc
clear all;
%Ge uppskattad frekvens & period med konstant L_0
[period,freq] = calc_freq();
pause

%RungeKutta 4 för att ge strömkurvorna I(t) med antalet steg (steps) och
%u-värden U_0s.
% konstanter
C = 5e-7;
L_0 = 0.7;
hs = [1e-6,5e-7,25e-8,125e-9, 625e-10]; 
start = 0;
periods = [period,1.5*period, 2*period];
U_0s = [220,1500, 2300];
for stop = periods
   for h = hs
      for U_0 = U_0s
          results = runge_kutta(U_0, start, h, stop, C, L_0);
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
          Es = E_const(U_0, start, h, stop, C, L_0);
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
   [hs, inter_period, period_errs,inter_max,  max_errs] = interpol_errors(U_0s, start, hs, stop,certainty, C, L_0);
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
    f = @(u) f_period(u,start,h,stop,certainty, wanted_period, C, L_0);
    U_stars(i) = sekant(f, guess1,guess2,certainty);
    fprintf('\n finding U* = %d for h = %d \n', U_stars(i),h);
end
%use most accurrate h
h = hs(end);

inter_periods = zeros(1,size(U_stars,2));
inter_maxs = zeros(1,size(U_stars,2));
u_i = 1;
for u = U_stars
    [cx,max_x,max_y,period] = interpol(u,start,h, stop, certainty, C, L_0);
    fprintf('\n Period = %d for U* = %d \n', period, u);
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

for i = 2:size(U_stars,2)
   diffs(i-1) = abs(U_stars(i)-U_stars(i-1)); 
end

%% störningsanalys
Ls = [L_0*0.95,L_0*1.05];
Cs = [C*0.95,C*1.05];
consts = [Cs(1), Cs(2), Cs(1), Cs(2);
          Ls(1), Ls(1), Ls(2), Ls(2)];
      
U_0_stars = zeros(1,4);
I_maxs = zeros(1,4);
% for every const combiniation
for c_i = 1:size(consts, 4)
    % TODO: beräkna U_0_star och I_max
end

% calc diffs
U_0_diff = abs(max(U_0_stars) - min(U_0_stars)) ;
I_max_diff = abs(max(I_maxs) - min(I_maxs));

% beräkna procentsats av riktiga
h_i = 3; % välj vilket h-värde 
U_0_proc = U_0_diff*100/U_stars(h_i);
I_max_proc = I_max_diff*100/I_maxs(h_i);

fprintf('En störning i L_0 och C på 5% ger ett fel på\nU_0*: %d% I_max*: %d%\n',U_0_proc,I_max_proc); 
%% Utvidgning: v(t)
h = hs(3);
u_i = 1;

% Plotta den omformade v(t) för de tre U_0
for u = U_0s
    for h_i = 1:size(hs,2)
        h = hs(h_i);
        [cx,~,y_max,period] = interpol(u,start,h,stop,certainty, C, L_0);
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
            as(u,h_i,k) = (1/pi)*trapz(xx, a(k,xx)); 
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
        plot(hs(2:end), diff_a(:,k),'DisplayName',['diff: a_', num2str(k)]);
        hold on;
    end
    title(['felvärden/noggranhet a_k för U_0 = ', num2str(U_0s(u))]);
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
    [cx,~,y_max,period] = interpol(u,start,h,stop,certainty, C, L_0);
    v = @(xs) ppval(cx,xs*period/(2*pi))/y_max;
    xx = [0:h*2*pi/period:2*pi];
    y = zeros(1,size(xx,2)*400);
    for i = 0:1:399
        y(i*size(xx,2)+1:(i+1)*size(xx,2)) = v(xx);
    end
    plot(y);
    Fs = 400*size(xx,2);
    %sound(y,Fs);
    audiowrite(['audio_',num2str(u),'.ogg'],y,Fs);
    pause;
end

%% Mystery sound part
load('mysterysound1');
load('mysterysoundlong');
plot(v);
title('short mystery sound')
pause;
plot(y);
title('mysterysound 400')
audiowrite('mysterylong.wav',y,400*size(v,2));
pause;
period = 400;
step_overs = [1,0];

a = @(k,h,step_over) v(1:(step_over+1):end).*sin(k.*[0:h:2*pi]);

for step_i = 1:size(step_overs,2)
    step_over = step_overs(step_i);
    h = (step_over+1)*2*pi/period;
    mys_hs(step_i) = h;
    xx = [0:h:2*pi];
    coefs = 10;
    for k = 1:coefs
        mys_as(step_i,k) = (1/pi)*trapz([0:h:2*pi], a(k,h,step_over)); % TODO: beräkna för h/2 och jämför för att visa 4 siffrors nogrannhet
    end
    mys_kvot(step_i) = mys_as(step_i,3)/mys_as(step_i,1);
    fprintf('\n kvot a3/a1 = %d h = %d \n', mys_kvot(step_i),h);
end

for step_i = 2:size(step_overs,2)
    diff_kvot(step_i-1) = abs(mys_kvot(step_i) - mys_kvot(step_i-1));
end
fprintf('\n Diff a3/a1 = %d \n', diff_kvot(:));
pause

%ta den mest accurate steglängden (alla punkter)
step_i = size(step_overs,2);
% definera fourierutvecklingarna med beräknade koefficienter
mystery_fourier{1} = @(t) mys_as(step_i,1)*sin(t) + mys_as(step_i,2)*sin(2*t) + mys_as(step_i,3)*sin(3*t);
mystery_fourier{2} = @(t) mys_as(step_i,1)*sin(t) + mys_as(step_i,2)*sin(2*t) + mys_as(step_i,3)*sin(3*t) + mys_as(step_i,4)*sin(4*t) + mys_as(step_i,5)*sin(5*t) + mys_as(step_i,6)*sin(6*t) + mys_as(step_i,7)*sin(7*t) + mys_as(step_i,8)*sin(8*t) + mys_as(step_i,9)*sin(9*t) + mys_as(step_i,10)*sin(10*t);

plot(xx, v,xx, mystery_fourier{1}(xx),xx,mystery_fourier{2}(xx));
legend({'mysterysound','Fourierutveckling med 3 termer','Fourierutveckling med 10 termer'},'Location','southwest')
pause
hold off;

%% Bestäm U_0 för mysterysound
h = hs(3);
U_0_2 = 2000;
U_0_1 = 2300;
% definiera funktionen med interpolation
error = @(U_0) norm(mys_as-interpol_eval(U_0, 0, h, stop, certainty, C, L_0, size(v,2)-1));
U_0_mys = sekant(error, U_0_1, U_0_2, certainty*10);

%% Function declaration for U_0* calculations
function period = f_period(u,start,h,stop,certainty,wanted_period, C, L_0)
    [~,~,~,p] = interpol(u, start, h, stop,certainty, C, L_0);
    period = p - wanted_period;
end

%% Function declaration for U_0 mysterysound
function ys = interpol_eval(U_0,start, h,stop,certainty, C, L_0,steps)
    fprintf('Guess: %d\n', U_0);
    [cx,~,y_max,period] = interpol(U_0,start,h,stop,certainty, C, L_0);
    v = @(xs) ppval(cx,xs*period/(2*pi))/y_max;
    xx = [0:2*pi/steps:2*pi];
    ys = v(xx);
end