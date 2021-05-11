%funktionen
%E(t) = (C/2)*U^2 + (1/2)* L_0*I_0^2*ln(I_0^2+I^2);
% Through chain rule: 
% E'(t) = (C/2)*2U*U' + (1/2)*I'*2I*L_O*I_0^2/(I_0^2+I^2)
%={making it more apparent}= (C*U')*U+I*(I'*L)
%={one more step}= -I*U+I*U= 0 IT IS CONSTANT
clc
clear all
% konstanter
C = 5e-7;
L_0 = 0.7;
I_0 = 1;

% Från calc_freq
period = 2*pi*sqrt(C*L_0);

%funktionen med bara I/I' (U = L* I')
E = @(v) (25e-8)*(0.7/(1+v(1)^2))^2*v(2)^2 + 0.35*log(1+v(1)^2);
%tydligast är att plotta båda och se att de slår ut varandra
E_spole = @(v) (25e-8)*(0.7/(1+v(1)^2))^2*v(2)^2;
E_kondens = @(v) 0.35*log(1+v(1)^2);

% steglängd
start = 0;
stop = period*4;
U_0 = 2300;
%for U_0 = [220,1500,2300]
    % startvärden
    I = 0;
    d_I = U_0/L_0;
    y = [I;d_I];
    %Loopa h, se skumheter för 10^-6/10^-7 ge högre start för att se det som
    %konstant linje (ungefär)
    for i_h = [3:7]
        steps = (10^i_h)*(stop-start);
        h = (stop-start)/steps;
        fprintf('Doing h: %d \n',h);

        %Få I & I' värden
        results = runge_kutta(U_0,start,steps,stop);
        %Initialisera
        Es = zeros(1,size(results,2));
        %E_komb = zeros(1,size(results,2));
        %Räkna E med I& I'
        for i = [1:size(Es, 2)]
            %E_sp(i) = E_spole(results(:,i));
            %E_kond(i) = E_kondens(results(:,i));
            %E_komb(i) = E_kondens(results(:,i)) + E_spole(results(:,i));
            Es(i) = E(results(:,i));
        end

        %plot([start:h:stop], E_komb, '-g');
        hold on
        %plot([start:h:stop], E_komb, '-g');
        %plot([start:h:stop], E_komb, '-g');
        plot([start:h:stop], Es, '-r');
        %hold off
    end
%end
hold off
