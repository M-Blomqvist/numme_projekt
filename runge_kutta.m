function results = runge_kutta(U_0,start,steps,stop)
    % konstanter
    C = 5e-7;
    L_0 = 0.7;
    I_0 = 1;

    % Vektorvärd funktion y = [i;i']
    F = @(y) [y(2);((-y(1)*I_0^2-y(1)^3)/(C*L_0*I_0^2) + (2*y(1)*y(2)^2)/(I_0^2 + y(1)^2))];

    % Runge kutta helper funktioner
    s1 = @(y) F(y);
    s2 = @(y, h) F(y + (h/2)*s1(y));
    s3 = @(y, h) F(y + (h/2)*s2(y, h));
    s4 = @(y, h) F(y + h*s3(y, h));

    % Runge kutta 4
    theta = @(y,h) y + (h/6)*(s1(y) + 2*s2(y, h) + 2*s3(y,h ) + s4(y, h));
    
    % startvärden
    I = 0;
    d_I = U_0/L_0;
    y = [I;d_I];
    h = (stop-start)/steps;
    
    ts = [start:h:stop];
    results = zeros(2, size(ts, 2));
    for x = [1:size(ts, 2)]
        results(:,x) = y;
        y = theta(y, h);
    end
end
    
    
    
    
    
    