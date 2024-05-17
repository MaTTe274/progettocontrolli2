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
    Tp = (15+273.15);                % temperatura cella (gradi)

    % Convezione con l'ambiente
    h = 6.55;               % coefficiente di convezione (W / (m^2 * K))
    S = (pi * d * l);         % superficie cilindrica (m^2)
    Tamb = (24+273.15);              % temperatura ambiente (gradi)

% ATTUATORE SMA

K_max = 3.92 * 10^3;        % (N/m)
Tavg = 70+273.15;                  % (gradi)
Tdiff = 20;                 % (gradi)
dl_max = 0.6 * 10^-2;       % deformazione filo, 96% di l (metri)
c = 6.2;                    % costante

z_star = 4.2 * 10^-2;       % quota desiderata (metri)


% STATI DI EQUILIBRIO

K_star = (m * g * r2) / (r1 * (dl_max - (r1/r2) * z_star));

x2_e = Tavg + ((Tdiff/c) * log(K_star / (K_max - K_star)));
x1_e = sqrt((1 / Re) * (h * S * (x2_e - Tamb) + (1 / Rth) * (x2_e - Tp)));

% K linearizzato in x2_e

K_x2_e = K_max * (1 - 1 / (1 + exp(c * (x2_e - Tavg) / Tdiff)));

% derivata parziale di K

derivata_parziale_K = ((c * K_max) / Tdiff) * ((exp(c * ((x2_e - Tavg) / Tdiff))) / (1 + exp(c * ((x2_e - Tavg) / Tdiff)))^2);

% MATRICI

A_1_1 = -Re / Le;

A_2_1 = ((2 * Re) / Cth) * x1_e;

A_2_2 = (-1 / Cth) * (h * S + 1 / Rth);

A_4_2 = ((r1 * (dl_max - (r1/r2) * z_star)) / (r2 * (m + J / r2^2))) * derivata_parziale_K;

A_4_3 = (-r1^2 /( r2^2 * (m +( J / r2^2)))) * K_x2_e;

A_4_4 = -beta / (m + J / r2^2);

A = [A_1_1     0       0      0;
     A_2_1   A_2_2     0      0;
       0       0       0      1;
       0     A_4_2   A_4_3  A_4_4];


B = [1/Le;
     0;
     0;
     0];


C = [0 0 1 0];


D = 0;

pippo = ss(A,B,C,D);

s=tf('s')
G = tf(pippo);


[P,Z]=pzmap(G);
figure;
bode(G);
grid on;

figure;
pzmap(G);
grid on;

figure
step(G);

Ge= G*(1/s);

figure;
pzmap(Ge);
figure;
rlocus(Ge);


