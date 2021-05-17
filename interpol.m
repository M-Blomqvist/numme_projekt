function [cy,max_x,max_y,period_i] = interpol(U_0, start, h,stop)
    %find  points to interpolate between
    xx = [start:h:stop];
    results = runge_kutta(U_0,start, h, stop);
    cy = spline([start:h:stop],results(1,:));
    cy_prime = spline([start:h:stop],results(2,:));


    %räkna ut "roots" (där den korsar x axeln) för att få periodstid
    prev_val = [0,0];
    ii = 1;
    for i = [1:size(results(1,:),2)]
        x = xx(i);
        y = results(1,i);
        if y ~= 0 && prev_val(2) ~= 0
            if prev_val(2) * y < 0
                f = @(x) ppval(cy,x);
                roots(ii) = sekant(f,prev_val(1),x,1e-11);
                ii = ii +1;
            end
        elseif y == 0
           roots(ii) = x;
           ii = ii +1;
        end    
        prev_val = [x,y];
    end
    
    
    for i = [1:size(results(2,:),2)]
        x = xx(i);
        y_prime = results(2,i);
        if y_prime == 0 
           max_x = x;
           max_y = results(1,i);
           if max_y > 0
               break;
           end
        elseif prev_val(2) * y_prime < 0
            f = @(x) ppval(cy_prime,x);
            root = sekant(f,prev_val(1),x,1e-15);
            max_x = root;
            max_y = ppval(cy,root);
            if max_y > 0
               break;
            end
        end   
        prev_val = [x,y_prime];
    end
    
    % ta diffen mellan "roots", vilket ger T/2
    diffs = zeros(1, size(roots, 2)-1);
    for i = [2:size(roots,2)]
        diffs(i) = abs(roots(i)-roots(i-1));
    end
    period_i = mean(diffs)*2;
end