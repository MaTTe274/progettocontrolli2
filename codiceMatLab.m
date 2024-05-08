m = 10;                     % grammi
beta = 1.6;                   % (N*s)/m
Le = 55 * 10^-3;            % H
Re = 4.8;                   % ohm
Rth = 16.7;                 % K/W
Tp = 15;                    % gradi
Cth = 0.0251;       
d = 0.0002;
l = 0.15;
h = 6.55;                   % W/m^2 * K
S = pi*d*l;                 % superficie m^2
Tamb = 24;                  % gradi
r1 = 0.45 * 10^-2;          % metri
r2 = 4.5 * 10^-2;           % metri
J = 2 * 10^-4;              % kg*m^2
K_max = 3.92 * 10^3;        % N/m
Tavg = 70;                  % gradi
Tdiff = 20;                 % gradi
c = 6.2;                    
dl_max = 0.6 * 10^-2;       % metri
z_star = 4.2 * 10^-2;       % metri
g = 9.81;                   % m*s^2


K_star = (m*g*r2)/(r1*(dl_max - r1/(r2 * z_star)));

x2_e = Tavg + (Tdiff/c)*log(K_star/K_max - K_star);

x1_e = sqrt((1/Re) * h * S * (x2_e -Tamb) + (1/Rth) * (x2_e - Tp));

K = K_max * (1- 1/(1 + exp(c * (x2_e - Tavg)/Tdiff)));

derivata_parziale_K = (c*K_max / Tdiff) * (exp(c*(x2_e - Tavg)/Tdiff)) / (1 + exp(c*(x2_e - Tavg)/Tdiff))^2;


A = [-Re/Le 0 0 0;
     2*Re*x1_e/Cth (-1/Cth)*(h*S + 1/Rth) 0 0;
     0 0 0 1;
     0 (r1*(dl_max -(r1/r2)*z_star) / (r2*(m + J/r2^2)) * derivata_parziale_K) -r1^2/(r2^2*(m+J/r2^2) * K) -beta/(m + J/r2^2)];


B = [1/Le;
     0;
     0;
     0];

C = [0 0 1 0];


D = 0;

pippo = ss(A,B,C,D);

G = tf(pippo);

figure;
bode(G);


