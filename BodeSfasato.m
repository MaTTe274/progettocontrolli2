clear; close all; clc; 
% DATI

% SISTEMA MECCANICO

r1 = 0.45 * 10^-2;          % raggio interno puleggia (metri)
r2 = 4.5 * 10^-2;           % raggio esterno puleggia (metri)
J = 2 * 10^-4;              % inerzia (kg*m^2)
m = 10 * 10^-3;             % massa (kilogrammi)
beta = 1.6;                 % coefficiente smorzatore (N * s / m)
g = 9.81;                   % accelerazione di gravità (m / s^2)

% SISTEMA DI RISCALDAMENTO

l = 15 * 10^-2;             % lunghezza filo (metri)
d = 0.2 * 10^-3;            % diametro filo (metri)
Cth = 25.1 * 10^-3;         % capacità termica (J/K)
Re = 4.8;                   % resistenza elettrica (ohm)
Le = 55 * 10^-3;            % induttanza (H)
    
    % Cella di Peltier
    Rth = 16.7;             % resistenza termica (ohm)
    Tp = 15+273.15;                % temperatura cella (gradi)

    % Convezione con l'ambiente
    h = 6.55;               % coefficiente di convezione (W / (m^2 * K))
    S = pi * d * l;         % superficie cilindrica (m^2)
    Tamb = 24+273.15;              % temperatura ambiente (gradi)

% ATTUATORE SMA

Kmax = 3.92 * 10^3;        % (N/m)
Ta = 70+273.15;                  % (gradi)
Td = 20;                 % (gradi)
deltalmax = 0.6 * 10^-2;       % deformazione filo, 96% di l (metri)
c = 6.2;                    % costante

zstar = 4.2 * 10^-2;       % quota desiderata (metri)






Kstar = m*g*r2/(r1*(deltalmax-r1*zstar/r2));
x2eq = Ta+Td*log(Kstar/(Kmax-Kstar))/c;
Kx2eq = Kmax*(1-(1/(1+exp(c*(x2eq-Ta)/Td))));
x1eq = sqrt((h*S*(x2eq-Tamb)/Re)+((x2eq-Tp)/Rth));



x11 = -Re/Le;
x21 = 2*Re*x1eq/Cth;
x22 = -(h*S+1/Rth)/Cth;
x42 = r1*(deltalmax-r1*zstar/r2)/r2*(m+J/r2^2)*(exp(c*(x2eq-Ta)/Td)/(1+exp(c*(x2eq-Ta)/Td))^2);
x43 = -r1^2*Kx2eq/(r2^2*(m+J/r2^2));
x44 = -beta/(m+J/r2^2);

A = [x11  0   0   0  ;
     x21 x22  0   0  ;
      0   0   0   1  ;
      0  x42 x43 x44];

B = [1/Le;
      0  ;
      0  ;
      0  ];

C = [0 0 1 0];

D = 0;

sistema = ss(A,B,C,D);

G = tf(sistema);

figure;
bode(G);

figure;
pzmap(G);

figure
step(G);