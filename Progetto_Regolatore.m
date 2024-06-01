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

%mettiamo il polo nell'origine e il guadagno per attenuare di 3000 volte in 0.1

[m_G_01,f_G_01] = bode(G/s,0.1);

guadagno_s = 15.84/m_G_01;


R = guadagno_s/s;

Ge = G * R;

[m_Ge_01,f_Ge_01] = bode(Ge,0.1);


figure;
margin(Ge,W);
%in questo modo abbiamo ampiezza di 3000 alla pulsazione 0.1

%progettiamo il regolatore Rd

[m_Ge_30,f_Ge_30] = bode(Ge,30);

%per attenuare di 100 volte il disturbo a 30 r/s dobbiamo ottenere m_Ge_30 = 0.01

%progettiamo una rete ritardatrice che ci dia wc minore per cercare di ottenerlo


%prima rete anticipatrice

%per attraversare a 6 dobbiamo dare circa 120 di anticipo di fase


Ra_1 = (1+1.42*s)/(1+0.0142*s);

L1 = Ge * Ra_1;

figure;
margin(L1,W);

L2 = progettaRA(L1,5,70);

figure;
margin(L2,W);

poli_L2 = pole(L2);
zeri_L2 = zero(L2);
[m_L2_01,f_L2_01] = bode(L2,0.1);
[m_L2_30,f_L2_30] = bode(L2,30);


L3 = L2 /(1+s/15);


figure;
margin(L3,W);


[m_L3_01,f_L3_01] = bode(L3,0.1);
[m_L3_30,f_L3_30] = bode(L3,30);


F3 = minreal(L3/(1+L3));
figure;
stepplot(F3);


Rff = (1+s/2.4103)*(1+s/5.372)*(1+s/9.3533)/((1+s/24.103)^3);


F4 = minreal(L3/(1+L3)+G*Rff/(1+L3));
figure;
stepplot(F4);

