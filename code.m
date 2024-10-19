% SCRIPT MATLAB PROGETTO SISTEMI ELETTRICI DI POTENZA 2024
% REINALDO AVELINO DA SILVA JUNIOR 
% 322110

clear; clc;

% Definizione dei nodi e dei rami 
% 26 rami; 22 nodi
nodo1 = zeros(26);
nodo2 = zeros(26);

%                Stazioni
nodo1(1) = 1; % Trapani-Salemi 
nodo2(1) = 2; % Partanna
nodo1(2) = 2;
nodo2(2) = 3; % Partinico
nodo1(3) = 3;
nodo2(3) = 4; % Bellolampo
nodo1(4) = 4;
nodo2(4) = 5; % Termini
nodo1(5) = 5;
nodo2(5) = 6; % Caracoli
nodo1(6) = 5;
nodo2(6) = 7; % Ciminna
nodo1(7) = 7;
nodo2(7) = 3;
nodo1(8) = 5;
nodo2(8) = 8; % Corriolo
nodo1(9) = 8;
nodo2(9) = 9; % S. F. Mela
nodo1(10) = 9;
nodo2(10) = 11; % Sorgente
nodo1(11) = 11;
nodo2(11) = 8;
nodo1(12) = 11;
nodo2(12) = 12; % Misterbianco
nodo1(13) = 12;
nodo2(13) = 13; % Melili
nodo1(14) = 13;
nodo2(14) = 14; % Priolo
nodo1(15) = 13; 
nodo2(15) = 15; % Anapo
nodo1(16) = 13;
nodo2(16) = 16; % Ragusa
nodo1(17) = 16;
nodo2(17) = 18; % Chiaramonte Gulfi
nodo1(18) = 18;
nodo2(18) = 14; %% errore: non ho considerato il generatore ISAB
nodo1(19) = 21;
nodo2(19) = 19; 
nodo1(20) = 19; % Patemo
nodo2(20) = 22;
nodo1(21) = 18;
nodo2(21) = 20; % Favara
nodo1(22) = 20;
nodo2(22) = 2;
nodo1(23) = 16;
nodo2(23) = 17; 
nodo1(24) = 22; 
nodo2(24) = 10; % Scilla
% Rami trasformatore
nodo1(25) = 18; 
nodo2(25) = 21; % Trafo Chiaramonte 220-280
nodo1(26) = 11;
nodo2(26) = 22; % Trafo Sorgente 380-220

% Distanze misurate con Google Maps in [m]
a = [44; 40; 17; 40; 3; 15; 40; 140; 4; 4; 4; 74;
    39; 7; 7; 50; 11; 40; 62; 70; 97; 80; 115; 39]'*1e3;

% Costanti fisiche
v0 = 299.792458*1e6; % [m/s]
mu0 = 4*pi*1e-7; % [H/m]
eps0 = 1/(mu0*v0^2); % [F/m]
Zv = sqrt(mu0/eps0); % [Ohm]

% Dati delle linee (alluminio)
V1 = 380e3; % [V]
V2 = 220e3; % [V]
n1 = 3;
n2 = 1;
f = 50; % [Hz]
w = 2*pi*f;

Sb = 1000e6; % [VA] potenza di base
Zb1 = V1^2/Sb;
Zb2 = V2^2/Sb;
Yb1 = 1/Zb1;
Yb2 = 1/Zb2;

Dab1 = 7.4; % [m]
Dbc1 = 7.4; % [m]
Dca1 = 14.8; % [m]
Dab2 = 9.5; % [m]
Dbc2 = 8.4; % [m]
Dca2 = 6.1; % [m]
d1 = 31.5e-3; % [m]
d2 = 26.9e-3; % [m]
Df = 0.4; % [m]
S1 = 519.5e-6; % [m2]
S2 = 349.2e-6; % [m2]
gama = 34e6; % [Sm/m2]
rho = 1/gama; % [Ohm/m]
B = 228; % [C]
theta = 50; % [C]
theta0 = 20; % [C]
KL1 = 0.810; % fattore per l'induttanza interna
KL2 = 0.826;

% Dati delle linee in cavo
x = 0.17e-3; % [Ohm/m]
c = 250e-12; % [F/m]
r3 = 0.02e-3; % [Ohm/m] tipo 1 380 kV
r4 = 0.06e-3; % [Ohm/m] tipo 2 220 kV

% Calcoli Zcc trafi (assumendo Zcc = j*Xcc)
Xcc = 14/100*(380e3)^2/400e6; % [Ohm]
Xcc = Xcc/Zb1;

% Calcoli rlc delle linee
Dm1 = (Dab1*Dbc1*Dca1)^(1/3);
Dm2 = (Dab2*Dbc2*Dca2)^(1/3);
Delta = Df/sin(pi/3);
deqL1 = (3*KL1*d1*(Delta)^(3-1))^(1/3);
deqL2 = KL2*d2;
deqC1 = (3*d1*(Df/sin(pi/3))^(3-1))^(1/3);
deqC2 = d2;
l1 = mu0/(2*pi)*log(2*Dm1/deqL1);
l2 = mu0/(2*pi)*log(2*Dm2/deqL2);
c1 = 2*pi*eps0/log(2*Dm1/deqC1);
c2 = 2*pi*eps0/log(2*Dm2/deqC2);
r1 = rho/(n1*S1)*(B+theta)/(B+theta0);
r2 = rho/(n2*S2)*(B+theta)/(B+theta0);

% Ammettenze
Z01 = sqrt((r1+j*w*l1)/(j*w*c1));
gama1 = sqrt((r1+j*w*l1)*(j*w*c1));
Z02 = sqrt((r2+j*w*l2)/(j*w*c2));
gama2 = sqrt((r2+j*w*l2)*(j*w*c2));
Z03 = sqrt((r3+j*x)/(j*w*c));
gama3 = sqrt((r3+j*x)*(j*w*c));
Z04 = sqrt((r4+j*x)/(j*w*c));
gama4 = sqrt((r4+j*x)*(j*w*c));
Z01 = Z01/Zb1;
Z02 = Z02/Zb2;
Z03 = Z03/Zb1;
Z04 = Z04/Zb2;
ZL1 = Z01*sinh(gama1.*a);
YT1 = tanh(gama1.*a/2)/Z01;
ZL2 = Z02*sinh(gama2.*a);
YT2 = tanh(gama2.*a/2)/Z02;
ZL3 = Z03*sinh(gama3.*a);
YT3 = tanh(gama3.*a/2)/Z03;
ZL4 = Z04*sinh(gama4.*a);
YT4 = tanh(gama4.*a/2)/Z04;

yl12 = zeros(26);
yt12 = zeros(26);
yt21 = zeros(26);

for r = 1:16 % tra il ramo 1 e il ramo 16
    % tensione 220 kV - caso 2
    yl12(r) = 1/ZL2(r);
    yt12(r) = YT2(r);
    yt21(r) = YT2(r);
end

yl12(17) = (1/ZL2(17))*2; % doppia terna
yt12(17) = YT2(17)*2; % Chiaramonte - Ragusa
yt21(17) = YT2(17)*2;

for r = 18:20 % tra il ramo 18 e il ramo 20
    % tensione 380 kV - caso 1
    yl12(r) = 1/ZL1(r);
    yt12(r) = YT1(r);
    yt21(r) = YT1(r);
end

for r = 21:22
    % tensione 220 kV - caso 2
    yl12(r) = 1/ZL2(r);
    yt12(r) = YT2(r);
    yt21(r) = YT2(r);
end

% Linee in cavo - 2 cavi 
yl12(23) = (1/ZL4(23))*2;
yt12(23) = YT4(23)*2;
yt21(23) = YT4(23)*2;
yl12(24) = (1/ZL3(24))*2;
yt12(24) = YT3(24)*2;
yt21(24) = YT3(24)*2;

% Trasformatori (circuito equivalente pi greco)
yl12(25) = 1/(j*Xcc);
yt12(25) = 0;
yt21(25) = 0;
yl12(26) = 1/(j*Xcc);
yt12(26) = 0;
yt21(26) = 0;

% Ammettenze trasversali
% Reattori di compensazione agli estremi XR = 1/BT => YR = -YT
nodo0 = [16, 17, 22, 10]; % nodi dei rami 23 e 24 (linee in cavo)
yt0 = [-YT4(23)*2, -YT4(23)*2, -YT3(23)*2, -YT3(23)*2];

% Algoritmo matrice ammettenze nodali Ybus
% 1) Costruzione dei nodi
Ybus = zeros(22, 22); % 22 nodi

% 2) Costruzione dei rami
for r=1:26
    k = nodo1(r);
    i = nodo2(r);
    Ybus(k,k) = Ybus(k,k) + yt12(r) + yl12(r);
    Ybus(i,i) = Ybus(i,i) + yt21(r) + yl12(r);
    Ybus(k,i) = Ybus(k,i) - yl12(r);
    Ybus(i,k) = Ybus(i,k) - yl12(r);
end

% 3) Costruzione delle ammettenze trasversali dei nodi
for t=1:4
    k = nodo0(t);
    Ybus(k,k) = Ybus(k,k) + yt0(t);
end

% Calcolo della matrice delle impedenze
Zbus = inv(Ybus);

% Rappresentazione grafica
spy(Ybus)
%spy(Zbus)

% Fine
disp('SEP24 / Reinaldo Avelino / 322110')
