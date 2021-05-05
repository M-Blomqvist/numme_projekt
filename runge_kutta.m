clc 
clear all;
% konstanter
C = 5e-7;
L_0 = 0.7;
I_0 = 1;
period = 2*pi*sqrt(C*L_0);

L = @(i) L_0*((I_0^2)/(I_0^2 + i^2));
d_L = @(i, d_i) (-2*L_0*I_0^2*i*d_i)/((I_0^2 + i^2)^2);

% % Systemet av ODE
% d_u = @(i) -i/C;
% d_i = @(u, i) u/L(i);
% 
% % Vektorvärda funktionen y = [u;i]
% F = @(y) [d_u(y(1));d_i(y(1), y(2))];

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
U_0 = 2300;
d_I = U_0/L_0;
d_U = -I/C;
y = [I;d_I];

% steglängd
h = period/1000;

ts = [0:h:period*4];
ys = zeros(2, size(ts, 2));
for x = [1:size(ts, 2)]
    ys(:,x) = y;
    y = theta(y, h);
end
   
%plot(ts, ys(1,:), '-r');
%hold on
plot(ts, ys(1,:), '-b');
hold on
    
    
    
    
    
    
    