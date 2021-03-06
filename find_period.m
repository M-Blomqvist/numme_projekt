function [U_star,cx,x_max,y_max,period]= find_period(wanted_period,guess1,guess2, start,h,stop, certainty, C, L_0)
    
    f = @(u) f_period(u,start,h,stop,certainty, wanted_period, C, L_0);
    U_star = sekant(f, guess1,guess2,certainty);
    
    %give piecewise interpolation of above U
    [cx,x_max,y_max,period] = interpol(U_star,start,h, stop, certainty, C, L_0);
end

%% Function declaration for U_0* calculations
function period = f_period(u,start,h,stop,certainty,wanted_period, C, L_0)
    [~,~,~,p] = interpol(u, start, h, stop,certainty, C, L_0);
    period = p - wanted_period;
end