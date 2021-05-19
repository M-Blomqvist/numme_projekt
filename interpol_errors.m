function [e_hs,inter_periods,period_errs,inter_maxs, max_errs,cxs] = interpol_errors(U_0s, start, hs,stop, certainty, C, L_0)
    %för-allokera
    inter_periods = zeros(size(U_0s,2),size(hs,2));
    inter_maxs = zeros(size(U_0s,2),size(hs,2),2);
    %största h-et har ingen noggranhets beräkning
    e_hs = hs(2:end);
    
    u_i = 1;
    
    %För varje U_0, hitta perioder och max värden för alla h
    for U_0 = U_0s
        h_i = 1;
        for h = hs
          %Använd interpol.m för att hitta perioder, interpolerade
          %funktionen och Imax för varje enskild U_0 & h kommbination.
          [cx,x_max,y_max,i_period] = interpol(U_0, start, h, stop,certainty, C, L_0);
          
          %Spara dessa
          inter_periods(u_i, h_i) = i_period;
          inter_maxs(u_i,h_i,:) = [x_max,y_max];
          cxs{u_i,h_i} = cx;
          
          h_i = h_i + 1;
        end
        u_i = u_i + 1;
    end
    
    period_errs = zeros(size(U_0s,2), size(hs,2)-1);
    max_errs = zeros(size(U_0s,2), size(hs,2)-1);
    
    u_i = 1;
    %Räkna ut errors genom att ta absoluta diffen i Imax & perioden för h/2
    %& h
    for U_0 = U_0s
        for h_i = 2:size(hs,2)
           max_errs(u_i,h_i-1) = abs(inter_maxs(u_i,h_i-1,2) - inter_maxs(u_i,h_i,2));
           period_errs(u_i,h_i-1) = abs(inter_periods(u_i,h_i-1) - inter_periods(u_i,h_i));
        end
        u_i = u_i +1;
    end  
end