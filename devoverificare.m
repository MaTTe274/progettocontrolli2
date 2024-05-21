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
A_4_2 = (r1*(dl_max-r1*z_star/r2)/r2*(m+J/r2^2))*derivata_parziale_K;
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

G = tf(sistema);


%mettiamo il polo nell'origine e il guadagno per attenuare di 3000 volte in 0.1

s = tf('s');
R = 203400/s;

Ge = G * R;
figure;
bode(Ge);
hold on;
grid on;

[m01,f01] = bode(Ge,0.1);


%in questo modo abbiamo ampiezza di 3000 alla pulsazione 0.1



%progettiamo il regolatore Rd

%per attenuare di 100 volte il disturbo a 30 r/s

[m30,f30] = bode(Ge,30);

%dobbiamo ottenere m30 = 0.01





%affinchè ci sia disaccopiamento frequenziale dobbiamo avere 0.1 < wc < 30

% c_smorzamento = Mf_star/100 = 0.55  =>  wc > 4.6/Ta*0.55 => wc > 7

% quindi per soddisfare la specifica sul tempo di assestamento bisogna
% avere wc > 7

% verifico che tra 7 e 30 esista un sottointervallo in cui abbiamo Mf = 55

% dato che questo sottointervallo non esiste siamo in uno scenario B

[m7,f7] = bode(Ge,7);

% alla pulsazione 7 r/s la fase è -255 quindi dobbiamo guadagnare 130 di
% fase  => mettiamo 2 zeri 






Wcd = 7;
MFd = 55;

[M,P,W]=bode(Ge);
[V,i]=min(abs(W-Wcd));
GeWcd=M(i);
ArgGeWcd=P(i);
Pd=-180+MFd-ArgGeWcd;
Md=1/GeWcd;
Pd=Pd*pi/180; 

if (Md<1 || Pd<0 || cos(Pd)>1/Md)


disp('errore nel tuning analitico');

end

tau=(Md-cos(Pd))/(Wcd*sin(Pd));
alpha=(cos(Pd)-1/Md)/(Wcd*sin(Pd))/tau; 

disp(tau);

disp(alpha);




