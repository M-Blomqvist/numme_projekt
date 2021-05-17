function [e_hs,inter_periods,period_errs,inter_maxs, max_errs] = interpol_errors(U_0s, start, hs,stop, certainty)
    inter_periods = zeros(size(U_0s,2),size(hs,2));
    inter_maxs = zeros(size(U_0s,2),size(hs,2));
    
    u_i = 1;
    for U_0 = U_0s
        s_i = 1;
        for h = hs
          [cx,x_max,y_max,i_period] = interpol(U_0, start, h, stop,certainty);
          txt = ['I(t) med interpolerade maxvärden, period = ', num2str(i_period)];
          inter_periods(u_i, s_i) = i_period;
          inter_maxs(u_i,s_i) = y_max;
          xx = [start:h:stop];
          plot(xx, ppval(cx,xx),'DisplayName',txt);
          hold on;
          %plot(x_max,y_max, 'o', 'DisplayName',['max värde för ovan']);
          plot([start:i_period:stop],ppval(cx,[start:i_period:stop]), 'o', 'DisplayName',['period för ovan']);
          s_i = s_i + 1;
        end
        u_i = u_i + 1;
    end
    title('Interpolerat I(t) med maxvärden');
    %legend('Location','southwest');
    pause;
    hold off;
    
    period_errs = zeros(size(U_0s,2), size(hs,2)-1);
    max_errs = zeros(size(U_0s,2), size(hs,2)-1);
    u_i = 1;
    for U_0 = U_0s
        for s_i = 2:size(hs,2)
           max_errs(u_i,s_i-1) = abs(inter_maxs(u_i,s_i-1) - inter_maxs(u_i,s_i));
           period_errs(u_i,s_i-1) = abs(inter_periods(u_i,s_i-1) - inter_periods(u_i,s_i));
        end
        u_i = u_i +1;
    end  
    e_hs = hs(2:end);
end