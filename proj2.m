% konstanter
C = 5e-7;
L_0 = 0.7;
I_0 = 1;

% Från calc_freq
period = 2*pi*sqrt(C*L_0);

results = runge_kutta(220,0,1000,4*period);

plot([0:4*period/1000:4*period], results(1,:), '-b');
hold on

results = runge_kutta(1500,0,1000,4*period);

plot([0:4*period/1000:4*period], results(1,:), '-r');

results = runge_kutta(2300,0,1000,4*period);

plot([0:4*period/1000:4*period], results(1,:), '-g');

hold off;
