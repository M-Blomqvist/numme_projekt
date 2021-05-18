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
    %No need to check different hs in this excercise  
    h = hs(2);
      for U_0 = U_0s
          results = runge_kutta(U_0, start, h, stop, C, L_0);
          txt = ['I(t) med U_0 = ',num2str(U_0), ' & h = ', num2str(h)];
          plot([start:h:stop], results(1,:),'DisplayName',txt);
          hold on;
      end
      title(['RungeKutta4 I(t)']);
      xlabel('t');
      ylabel('I');
      legend('Location','southeast');
      pause;
      hold off;
end

%% Visa att E(t) är konstant:
for stop = periods
    for h = hs
        for U_0 = U_0s
              Es = E_const(U_0, start, h, stop, C, L_0);
              txt = ['E(t) med U_0 = ',num2str(U_0)];
              plot([start:h:stop], Es,'DisplayName',txt);
              hold on;
        end
        title(['RungeKutta4 E(t) med h = ', num2str(h)]);
        xlabel('t');
        ylabel('E');
        legend('Location','east');
        pause;
        hold off;
    end
end
%% Interpolate to find max_i and period
certainty = 1e-14;
fprintf('Sekant certainty (break when diff less than) %d \n', certainty);
for stop = periods(2)
   [e_hs, inter_period, period_errs,inter_max,  max_errs,cxs] = interpol_errors(U_0s, start, hs, stop,certainty, C, L_0);
   for u_i = 1:size(U_0s,2)
       %plotta endast de med minsta h
       i_period = inter_period(u_i,end);
       xx = [start:hs(end):stop];
       cx = cxs{u_i,end};
       x_max = inter_max(u_i,end,1);
       y_max = inter_max(u_i,end,2);
       txt = ['I(t) interpolerad U_0 = ',num2str(U_0s(u_i)), ' h = ', num2str(hs(end))];
       plot(xx, ppval(cx,xx),'DisplayName',txt);
       hold on;
       plot(x_max,y_max, 'o', 'DisplayName',['max värde för ovan ', num2str(y_max)]);
       plot([start:i_period:stop],ppval(cx,[start:i_period:stop]), '*', 'DisplayName',['period för ovan ',num2str(i_period)]);
   end
   title('Interpolerat I(t) med maxvärden');
   xlabel('t');
   ylabel('I');
   legend('Location','southeast');
   pause;
   hold off;
end
%plot errors in period
u_i = 1;
for U_0 = U_0s
    txt = ['U_0 = ',num2str(U_0)];
    loglog(e_hs,period_errs(u_i,:),'DisplayName', txt);
    hold on;
    u_i = u_i +1;
end
loglog(e_hs,e_hs.^2,'DisplayName', 'O(h^2)');
title('Felvärden interpolation: perioder');
xlabel('h');
ylabel('T_h-T_{2h}');
legend('Location','northwest');
hold off;
pause

%plot errors in max I
u_i = 1;
for U_0 = U_0s
    txt = ['U_0 = ',num2str(U_0)];
    loglog(e_hs,max_errs(u_i,:),'DisplayName',txt);
    hold on;
    title('Felvärden interpolation: max värden');
    xlabel('h');
    ylabel('{max(I)}_h - {max(I)}_{2h}');
    u_i = u_i +1;
end
loglog(e_hs,e_hs.^2,'DisplayName', 'O(h^2)');
title('Felvärden interpolation: max värden');
xlabel('h');
ylabel('{max(I)}_h - {max(I)}_{2h}');
legend('Location','southeast');
hold off;
for i = 1:size(U_0s,2)
    U_0 = U_0s(i);
    fprintf(['For U_0 = ',num2str(U_0),'\n T = ', num2str(inter_period(i,end)), ' err = ', num2str(period_errs(i,end)),'\n max(I) = ',num2str(inter_max(i,end)), ' err = ', num2str(max_errs(i,end)),' \n']);
end   
pause
%% Hitta U_0* så att 1/period = 400 = 0.0025, U_0 = 220 & 1500 best candidate

wanted_period = 400^-1;
stop = periods(2);
guess1 = 220;
guess2 = 1500;
certainty = 1e-14;
diff = 100;
U_stars = zeros(1,size(hs,2));
star_periods = zeros(1,size(U_stars,2));
star_maxs = zeros(2,size(U_stars,2));
%Find U* for different hs;
for i = [1:size(hs,2)]
    h = hs(i);
    [U_star,cx,x_max,y_max,period] = find_period(wanted_period,guess1,guess2,start,h,stop,certainty,C,L_0);
    U_stars(i) = U_star;
    star_periods(i) = period;
    star_maxs(:,i) = [x_max;y_max];
    star_plots{i} = cx;
end

%use most accurrate h to plot
h = hs(end);
i = 1;
for u = U_stars
    txt = ['I(t),U*= ',num2str(u),' h = ', num2str(hs(i))];
    x = [start:h:stop];
    plot(x, ppval(star_plots{i},x),'DisplayName',txt);
    hold on;
    plot(star_maxs(1,i),star_maxs(2,i), 'o', 'DisplayName',['max värde för ovan ', num2str(star_maxs(2,i))]);
    plot([start:star_periods(i):stop],ppval(star_plots{i},[start:star_periods(i):stop]), '*', 'DisplayName',['period för ovan ',num2str(star_periods(i) )]);
    i = i +1;
end
title('plot I(t) med period 400 och U*');
legend('Location','southeast');
xlabel('t');
ylabel('I');
hold off;

%calculate differences (resulting errors = inaccuracy)
period_errs = zeros(1,size(U_stars,2)-1);
max_errs = zeros(1,size(U_stars,2)-1);
ustar_errs = zeros(1,size(U_stars,2)-1);
for u_i = 2:size(U_stars,2)
   max_errs(u_i-1) = abs(star_maxs(2,u_i) - star_maxs(2,u_i-1));
   period_errs(u_i-1) = abs(star_periods(u_i) - star_periods(u_i-1));
   ustar_errs(u_i-1) = abs(U_stars(u_i)-U_stars(u_i-1)); 
end
%plotta skillnader i U*
loglog(hs(2:end),ustar_errs,'DisplayName', 'U*');
hold on;
loglog(hs(2:end),hs(2:end).^2,'DisplayName', 'O(h^2)');
title('Felvärden: U*');
xlabel('h');
ylabel('{U*}_h-{U*}_{2h}');
legend('Location','northwest');
hold off;
pause;

%plotta skillnader i perioder (om det inte finns någon skillnad skippa plotten)
if any(period_errs)
    loglog(hs(2:end),period_errs,'DisplayName', 'T');
    hold on;
    %loglog(hs(2:end),hs(2:end).^2,'DisplayName', 'O(h^2)');
    title('Felvärden: Perioder');
    xlabel('h');
    ylabel('{T}_h-{T}_{2h}');
    legend('Location','northwest');
    hold off;
    pause;
else
    fprintf('nothing to plot! Periods are identical for all found U*s! \n differences are guaranteed below %d \n', certainty);
end

%plotta skillnader i max(I)
loglog(hs(2:end),max_errs,'DisplayName', 'max(I)');
hold on;
loglog(hs(2:end),hs(2:end).^2,'DisplayName', 'O(h^2)');
title('Felvärden: max(I)');
xlabel('h');
ylabel('{max(I)}_h-{max(I)}_{2h}');
legend('Location','northwest');
hold off;
pause;

%% störningsanalys
Ls = [L_0*0.95,L_0*1.05];
Cs = [C*0.95,C*1.05];
consts = [Cs(1), Cs(2), Cs(1), Cs(2);
          Ls(1), Ls(1), Ls(2), Ls(2)];     
U_0_stars = zeros(1,4);
I_maxs = zeros(1,4);
certainty = 1e-12;

%Use most correct h;
h = hs(end);
% for every const combiniation
for c_i = 1:4
    C_ = consts(1,c_i);
    L_ = consts(2,c_i);
    fprintf('Finding U* with c = %d & L_0 = %d \n', C_,L_);
    [U_star,~,~,y_max,~] = find_period(wanted_period,guess1,guess2,start,h,stop,certainty,C_,L_);
    U_0_stars(c_i) = U_star;
    I_maxs(c_i) = y_max;
end

% calc diffs
U_0_diff = abs(max(U_0_stars) - min(U_0_stars)) ;
I_max_diff = abs(max(I_maxs) - min(I_maxs));

% beräkna procentsats av riktiga
h_i = 3; % välj vilket h-värde 
U_0_proc = U_0_diff*100/U_stars(h_i);
I_max_proc = I_max_diff*100/I_maxs(h_i);

fprintf('En störning i L_0 och C på 5%% ger ett fel på \n U_0*: %d%% I_max*: %d%% \n',U_0_proc,I_max_proc); 
%% Utvidgning: v(t)
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
xlabel('t');
ylabel('v');
legend('Location','southwest');
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
        loglog(hs(2:end), diff_a(:,k),'DisplayName',['diff: a', num2str(k)]);
        hold on;
    end
    loglog(hs(2:end),hs(2:end).^2,'DisplayName', 'O(h^2)');
    title(['felvärden/noggranhet a_k för U_0 = ', num2str(U_0s(u))]);
    legend('Location','southeast');
    xlabel('h');
    ylabel('{a_n}_{ h} - {a_n}_{ 2h}');
    pause
    hold off;
end

%choose h to continue with 3
h_i = size(hs,2);
for u = 1:size(U_0s,2)
    % definera fourierutvecklingarna med beräknade koefficienter
    fourier_funcs{u,1} = @(t) as(u,h_i,1)*sin(t) + as(u,h_i,2)*sin(2*t) + as(u,h_i,3)*sin(3*t);
    fourier_funcs{u,2} = @(t) as(u,h_i,1)*sin(t) + as(u,h_i,2)*sin(2*t) + as(u,h_i,3)*sin(3*t) + as(u,h_i,4)*sin(4*t) + as(u,h_i,5)*sin(5*t) + as(u,h_i,6)*sin(6*t) + as(u,h_i,7)*sin(7*t) + as(u,h_i,8)*sin(8*t) + as(u,h_i,9)*sin(9*t) + as(u,h_i,10)*sin(10*t);
    
    plot(xx, fourier_funcs{u,1}(xx),xx,fourier_funcs{u,2}(xx),xx, vs{u,h_i}(xx));
    legend({'Fourierutveckling med 3 termer','Fourierutveckling med 10 termer',['v(t) för U_0 = ',num2str(U_0s(u))]},'Location','southwest')
    xlabel('t');
    ylabel('v');
    pause
end
% Plotta de udda a-värderna för alla U_0
for u = 1:size(U_0s,2)
    semilogy(abs(squeeze(as(u,h_i,1:2:end)))); 
    xlabel('k');
    ylabel('a_k');
    title(['semilogy-plot för udda a_k, U_0 = ',num2str(U_0s(u))])
    pause
end

%% sound-part
    h = 0.00001;
for u_i = 1:size(U_0s,2)
    u = U_0s(u_i);
    [cx,~,y_max,i_period] = interpol(u,start,h,stop,certainty, C, L_0);
    v = @(xs,p) ppval(cx,xs*p/(2*pi))/y_max;
    xx = [0:h*2*pi/i_period:2*pi];
    y = zeros(1,size(xx,2)*400);
    xxs = zeros(1,size(xx,2)*400);
    for i = 0:1:399
        xxs(i*size(xx,2)+1:(i+1)*size(xx,2)) = xx+(i*2*pi);
        y(i*size(xx,2)+1:(i+1)*size(xx,2)) = v(xx,i_period);
    end
    plot(xxs,y);
    hold off;
    title(['soundplot for U_0 =' num2str(u)]);
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

h_i = size(hs,2);
for u_i = 1:size(U_0s,2)
   kvot(u_i) = as(u_i,h_i,3)/as(u_i,h_i,1);
end
kvot(size(U_0s,2)+1) = mys_kvot(end);
for i = 1:size(U_0s,2)
   u = U_0s(i);
   kv = kvot(i);
   fprintf('\n a3/a1 = %d for %d', kv,u);
end
fprintf('\n a3/a1 = %d for mystery sound', kvot(end));

%ta den mest accurate steglängden (alla punkter)
step_i = size(step_overs,2);
% definera fourierutvecklingarna med beräknade koefficienter
mystery_fourier{1} = @(t) mys_as(step_i,1)*sin(t) + mys_as(step_i,2)*sin(2*t) + mys_as(step_i,3)*sin(3*t);
mystery_fourier{2} = @(t) mys_as(step_i,1)*sin(t) + mys_as(step_i,2)*sin(2*t) + mys_as(step_i,3)*sin(3*t) + mys_as(step_i,4)*sin(4*t) + mys_as(step_i,5)*sin(5*t) + mys_as(step_i,6)*sin(6*t) + mys_as(step_i,7)*sin(7*t) + mys_as(step_i,8)*sin(8*t) + mys_as(step_i,9)*sin(9*t) + mys_as(step_i,10)*sin(10*t);

plot(xx, v,xx, mystery_fourier{1}(xx),xx,mystery_fourier{2}(xx));
legend({'mysterysound','Fourierutveckling med 3 termer','Fourierutveckling med 10 termer'},'Location','southwest')
pause
hold on;

% Bestäm U_0 för mysterysound
%h = hs(3);
used_hs= hs;

i = 1;
for h = used_hs
    U_0_2 = 1500;
    U_0_1 = 2300;
    coefs = 10;
    % definiera funktionen med interpolation
    error = @(U_0) kvot(end)-interpol_eval(U_0, 0, h, stop, certainty, C, L_0, size(v,2)-1,coefs);
    U_0_mys(i) = sekant(error, U_0_1, U_0_2, certainty);
    [cx,~,y_max,period] = interpol(U_0_mys(i),start,h,stop,certainty,C,L_0);
    xx = [0:h*2*pi/period:2*pi];
    plot(xx,ppval(cx,xx*period/(2*pi))/y_max,'DisplayName',['mystery sound approximation plot med U_0 = ', num2str(U_0_mys(i)), ' h = ',num2str(h)]);
    i = i+1;
end

for i = 2:size(used_hs,2)
   diff_U_mys(i-1) = abs(U_0_mys(i)-U_0_mys(i-1));
   fprintf(['\n diff in U for mystery sound: ', num2str(diff_U_mys(i-1)), ' for h =', num2str(used_hs(i))]);
end
hold off;
%% Function declaration for U_0 mysterysound
function kvot = interpol_eval(U_0,start, h,stop,certainty, C, L_0,steps,coefs)
    fprintf('Guess: %d\n', U_0);
    [cx,~,y_max,period] = interpol(U_0,start,h,stop,certainty, C, L_0);
    v = @(xs) ppval(cx,xs*period/(2*pi))/y_max;
    xx = [0:h*2*pi/period:2*pi];
    a = @(k,xs) v(xs).*sin(k.*xs);
    for k = 1:coefs
        as(k) = (1/pi)*trapz(xx, a(k,xx)); 
    end
    kvot = as(3)/as(1);
end
