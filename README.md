# 3 Kretsen
### Före den numeriska behandlingen kan det vara bra att bedöma storleksordningen påsvängningstiden så att lämpligt tidssteg och simuleringsintervall kan väljas. Gör det genomatt analytiskt räkna ut frekvensen och svängningstiden för en krets medkonstantL=L0.

Då L = L0 = 0.7, C = 5* 10^-8

`I(t) = Acos(1/sqrt(C*L)*t)+Bsin(1/sqrt(C*L)*t)`

Svängningstid = 2 * pi * sqrt(C*L) ~= 1.1755 * 10^-4

Frekvens = 8.5072 * 10^3

Lämpligt intervall: [-svängtid * 2: svängtid * 2]

### Använd Runge–Kutta 4 för att beräkna strömkurvorna.  Några olika värden på *U0* skaprövas, dels spänningen 220 V då järnkärnans inflytande är nästan försumbart, dels tvåhöga spänningsvärden 1500 V och 2300 V då strömkurvan inte blir särskilt sinuslik längre.Plotta ett par perioder av de tre lösningarna

