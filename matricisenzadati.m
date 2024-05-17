
m = 10;             %grammi 
beta = 1.6;         %Newton*secondi/metri
Le = 0.055;         %H
Re = 4.8;           %ohm
Rth = 16.7;         %kelvin/watt
Tp = 288;           %kelvin
Cth = 0.0251;       %joule/kelvin
d = 0.0002;         %metri
l = 0.15;           %metri
h = 6.55;           %watt/metri^2*kelvin
S = pi*d*l;
Tamb = 297;         %kelvin


Kx2eq = Kmax*(1-(1/(1+e^((x2eq-Ta)/Td))));
Kstar = m*g*r2/r1*(deltalmax-r1*zstar/r2);
x2eq = Ta+Td*log(Kstar/Kmax-Kstar)/c;
x1eq = sqrt((hS*(Ta+Td*ln(m*g*r/r1*(dlmax-r1*zstar/r2/Cth)-Tamb)/Re))+((x2eq-Tp)/Rth));



x11 = -Re/Le;
x21 = 2*Re*x1eq/Cth;
x22 = -(h*S+1/Rth)/Cth;
x42 = (r1*(deltalmax-r1*zstar/r2)/r2*(m+J/r2^2))*(e^(c*(x2eq-Ta)/Td)/(1+e^(c*(x2eq-Ta)/Td))^2);
x43 = -r1^2*Kx2eq/(r2^2*(m+(J/r2^2)));
x44 = -beta/(m+(J/r2^2));

A = [x11  0   0   0  ;
     x21 x22  0   0  ;
      0   0   0   1  ;
      0  x42 x43 x44];

B = [1/le;
      0  ;
      0  ;
      0  ];

C = [0 0 1 0];

D = 0;

sistema = ss(A,B,C,D);

figure;
