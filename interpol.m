function [cy,max_x,max_y,period_i] = interpol(U_0, start, h,stop,certainty, C, L_0)
    %find  points to interpolate between
    xx = [start:h:stop];
    %resultat från runge kutta = [y; y'] för alla stegade t.
    results = runge_kutta(U_0,start, h, stop, C, L_0);
    cy = spline([start:h:stop],results(1,:));
    cy_prime = spline([start:h:stop],results(2,:));


    %räkna ut "roots" (där den korsar x axeln) för att få periodstid
    prev_val = [0,0];
    ii = 1;
    for i = [1:size(results(1,:),2)]
        x = xx(i);
        y = results(1,i);
        
        %Om värdet noll är vi på en rot
        if y ~= 0 && prev_val(2) ~= 0
            %Annars, om nya värdet har ett annat tecken än förra så vet vi
            %att en rot ligger emellan
            if prev_val(2) * y < 0
                %Definiera det sekant metoden ska söka nollpunkten för (den
                %interpolerade I(t) funktionen)
                f = @(x) ppval(cy,x);
                %Hitta var I = 0 med en viss tolerans och spara dessa t
                %värden
                roots(ii) = sekant(f,prev_val(1),x,certainty);
                ii = ii +1;
            end
           
        elseif y == 0
           roots(ii) = x;
           ii = ii +1;
        end    
        prev_val = [x,y];
    end
    
    %Perioden är distansen mellan första och tredje nollpunkten
    period_i = abs(roots(3)-roots(1));
    
    %räkna ut extrempunkter (där y' = 0) för att hitta Imax
    for i = [1:size(results(2,:),2)]
        x = xx(i);
        y_prime = results(2,i);
        %Om y' är 0 är vi på en extrempunkt
        if y_prime == 0 
           max_x = x;
           max_y = results(1,i);
           %Om extrempunkten är positiv har vi hittat vår Imax
           if max_y > 0
               break;
           end
        %Annars, om nya y' har ett annat tecken än förra så vet vi
        %att en extrempunkt ligger emellan
        elseif prev_val(2) * y_prime < 0
            %Definiera det sekant metoden ska söka nollpunkten för (den
            %interpolerade I(t) funktionen)
            f = @(x) ppval(cy_prime,x);
            root = sekant(f,prev_val(1),x,certainty);
            max_x = root;
            max_y = ppval(cy,root);
            %Om extrempunkten är positiv har vi hittat vår Imax
            if max_y > 0
               break;
            end
        end   
        prev_val = [x,y_prime];
    end
end