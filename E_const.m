%funktionen
%E(t) = (C/2)*U^2 + (1/2)* L_0*I_0^2*ln(I_0^2+I^2);
% Through chain rule: 
% E'(t) = (C/2)*2U*U' + (1/2)*I'*2I*L_O*I_0^2/(I_0^2+I^2)
%={making it more apparent}= (C*U')*U+I*(I'*L)
%={one more step}= -I*U+I*U= 0 IT IS CONSTANT
function Es = E_const(U_0,start,h, stop, C, L_0)
    I_0 = 1;
    I_0_2 = I_0^2;
    %funktionen med bara I/I' (U = L* I'), v = [I,I']
    E = @(v) (C/2)*((I_0_2*L_0)/(I_0_2+v(1)^2))^2*v(2)^2 + (1/2)*L_0*I_0_2*log(I_0_2+v(1)^2);

    fprintf('Calulating E for h: %d \n',h);

    %Få I & I' värden
    results = runge_kutta(U_0,start,h,stop, C, L_0);
    %Initialisera
    Es = zeros(1,size(results,2));
    %Räkna E med I& I'
    for i = [1:size(Es, 2)]
        Es(i) = E(results(:,i));
    end
end