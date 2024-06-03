clear; close all; clc; 

%FUNZIONI

function L = praticoRR(Ge_RR,wc_RR)

s = tf('s');

[m_RR,~] = bode(Ge_RR,wc_RR);

alpha_RR = 10^(-log10(m_RR));

tau_RR = 10/(wc_RR*alpha_RR);

RR = (1+alpha_RR*tau_RR*s)/(1+tau_RR*s);

L = Ge_RR * RR;

end


function L = progettaRR(G,wc,Mf)

s = tf('s');

[m,f] = bode(G,wc);

M = 10^(-log10(m));
phi = -180+Mf-f;

if (M>1 || phi>0 || phi<-90 || cosd(phi)<M)
disp('errore nel tuning della rete ritardatrice'); return;
end

tau = (cosd(phi)-(1/M))/(wc*sind(phi));
alpha = (M-cosd(phi))/(wc*sind(phi)*tau);

RR = (1+alpha*tau*s)/(1+tau*s);

L = G * RR;
end

function L = progettaRA(G,wc,Mf)

s = tf('s');

[m,f] = bode(G,wc);

M = 10^(-log10(m));
 
phi= -180+Mf-f;

if (M<1 || phi<0 || cosd(phi)<1/M)
disp('errore nel tuning della rete anticipatrice'); return;
end

tau = (M-cosd(phi))/(wc*sind(phi));

alpha = (cosd(phi)-(1/M))/(wc*sind(phi)*tau);

RA = (1+tau*s)/(1+alpha*tau*s);

L = G * RA;

end

function RA = progettaRegA(G,wc,Mf)

s = tf('s');

[m,f] = bode(G,wc);

M = 10^(-log10(m));
 
phi= -180+Mf-f;

if (M<1 || phi<0 || cosd(phi)<1/M)
disp('errore nel tuning della rete anticipatrice'); return;
end

tau = (M-cosd(phi))/(wc*sind(phi));

alpha = (cosd(phi)-(1/M))/(wc*sind(phi)*tau);

RA = (1+tau*s)/(1+alpha*tau*s);

end

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

K_max = 3.92 * 10^3;        % (N/m)
Tavg = 70+273.15;                  % (gradi)
Tdiff = 20;                 % (gradi)
dl_max = 0.6 * 10^-2;       % deformazione filo, 96% di l (metri)
c = 6.2;                    % costante

z_star = 4.2 * 10^-2;       % quota desiderata (metri)






K_star = m*g*r2/(r1*(dl_max-r1*z_star/r2));
x2_e = Tavg+Tdiff*log(K_star/(K_max-K_star))/c;
x1_e = sqrt((h*S*(x2_e-Tamb)+((x2_e-Tp)/Rth))/Re);

K_x2_e = K_max*(1-(1/(1+exp(c*(x2_e-Tavg)/Tdiff))));


derivata_parziale_K = (c*K_max/Tdiff)*(exp(c*(x2_e-Tavg)/Tdiff)/(1+exp(c*(x2_e-Tavg)/Tdiff))^2);


A_1_1 = -Re/Le;
A_2_1 = 2*Re*x1_e/Cth;
A_2_2 = -(h*S+1/Rth)/Cth;
A_4_2 = (r1*(dl_max-(r1*z_star/r2))/(r2*(m+J/r2^2)))*derivata_parziale_K;
A_4_3 = -r1^2*K_x2_e/(r2^2*(m+J/r2^2));
A_4_4 = -beta/(m+J/r2^2);

A = [A_1_1  0   0   0  ;
     A_2_1 A_2_2  0   0  ;
      0   0   0   1  ;
      0  A_4_2 A_4_3 A_4_4];

B = [1/Le;
      0  ;
      0  ;
      0  ];

C = [0 0 1 0];

D = 0;

sistema = ss(A,B,C,D);


s = tf('s');


G = tf(sistema);

display(G);

W = logspace(-1,1.477,100);

%PROGETTO INNER LOOP

Gm = 2292/(s+87.2727);
%attenuiamo il disturbo d con il guadagno

[m_Gm_01,f_Gm_01] = bode(Gm/s,0.1);

guadagno_Gm = 3000/m_Gm_01;

%mettiamo un polo nell'origine e uno zero in alta frequenza per avere un
%buon margine di fase (opzionale)

zero_Gm = 1+s/150;

Rm = guadagno_Gm*zero_Gm/s;

display(Rm);

Lm = Gm*Rm;

[m_Lm_01,f_Lm_01] = bode(Lm,0.1);

figure;
margin(Lm);

%stiamo attraversando a 200 quindi non ci dovrebbero essere problemi con il
%disaccoppiamento frequenziale

%PROGETTO OUTER LOOP

Gv = 0.5743/(s+2.4103)/(s+5.3572)/(s+9.3533);

%per rispettare la sepcifica sul Ta dobbiamo attraversare dopo 5.5

[m_Gv_30,f_Gv_30] = bode(Gv/s,30);
[m_Gv_20,f_Gv_20] = bode(Gv/s,20);

%mettiamo guadagno statico per aumentare la pulsazione di attraversamento



Gev = Gv/s;

figure;
margin(Gev,W);

[m_Gev_30,f_Gev_30] = bode(Gev,30);

Ra_1 = ((1+0.32*s)^2)/((1+0.017*s)^2);

taup_Gv = 1/17;

Rv = 700 * Ra_1 / (1+taup_Gv*s);

Lv = Gev * Rv;

display(Rv);

figure;
margin(Lv,W);

[m_Lv_30,f_Lv_30] = bode(Lv,30);



zeri_Lv = zero(Lv);


F = minreal(Lv/(1+Lv));

figure;
step(F);

info = stepinfo (F , SettlingTimeThreshold =0.01);

return
