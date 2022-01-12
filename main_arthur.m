%% Autonomie Automatique Lineaire 
% Arthur SOUTELO ARAUJO
% 2021 - 2022
%% 
clc
clear
close all
%% Grandeurs mécaniques, geométriques et eléctriques 
%Paramètres données
h = 0.95;               % hauteur du barycentre du conducteur
mh = 80;                % masse du conducteur
Ihx = 18;               % inertie du conducteur selon l'axe (A,x3)
ms = 25;                % masse du chariot
Is = 0.8;               % inertie du chariot selon l'axe (A,x3)
mr = 5;                 % masse d'une roue
Ir = 0.28;              % inertie d'une roue selon l'axe (A,x3)
r = 240 * 10^(-3);      % rayon d'une roue
g = 9.81;               % accélération de la pesanteur
kr = 200;               % rapport de réduction

ki = 6.5;               % gain de l'inclinomètre
kem = 0.116;            % constante de couple d'un moteur
gi = 0.1;               % gain de l'amplificateur de courant
kg = 2.29;              % gain du gyromètre

%Formules données
J1 = mh*h^2 + Ihx + Is;
M2 = ms + mh + 2*mr + 2*(Ir/r^2);

%% Fonction de Transfert (1) - Procedé mécanique
%Gain et pulsation naturelle
k1 = (mh*h+r*M2) / (mh*h*g*r*M2);
w1 = sqrt((mh*h*g*M2) / (J1*M2-(mh^2)*(h^2)));

s = tf('s');

TF1 = k1 / ((1/w1^2)*s^2 - 1);

%% Fonction de Transfert (2) - Asservissement
wn = w1 * 1.5;
ksi = 0.7;

kp = 22.708;
kv = 10.197;  
k2 = 0.00979;

s = tf('s');

TF2 = k2 / ((1/wn^2)*s^2+(2*ksi/wn)*s+1);

%% Fonction de Transfert (3) - Pertubation
s = tf('s');

TF3 = (w1^2 - s^2) / (s^2 + s*(2*k1*w1^2*kr*kem*gi*kv*kg) + (2*k1*w1^2*kr*kem*gi*kp*ki-w1^2));

%% Cahier des charges 
d = 0.05;
trep = 0.5;                 % secondes

%% Système de référence
ksiref = ksi ;
wnref = 3/trep; 
kref=1;                     % Erreur stationnaire nulle MREF(s=0) = 1

Mref = tf([1],[1/(wnref^2),2*ksiref/wnref,1]);
poleMref = pole(Mref);

% step(TF2);                % pour tracer la boucle ouverte
% hold
% step(TF2/(1+TF2));        % pour tracer la boucle ferme

%%  C  O  R  R  E  C  T  E  U  R
%% Méthodes Empiriques de Synthèse
% sisotool(TF2)
% ESSAI INDICIEL EN BO   =   Pas possible, présente dépassement
% ESSAI EN BF            =   Pas possible, il n'y a pas de pole avec partie réelle nulle 

%% Correcteur PID - Placement de Pôles

% On choisit un correcteur KPID(s) = Kc(1+1/(Ti*s)+Td*s/(1+Tf*s))
% KPID(s) = Kc(((Ti*Tf+Ti*Td)*s²+(Ti+Tf)*s+1/(Ti*Tf*s²+Ti*s))
% On simplifie en faisant compensation des pôles

Tf = 0.1193;  
Ti = 0.1067; 
Td = 0.13;  
kc = 46.79;

s = tf('s');

Kpid = kc * (1 + 1/(Ti*s) + (Td*s)/(1+Tf*s) );

%% Feed-Foward
s = tf('s');

FF = 45.36;

%% Correcteur RST - Placement de Pôles
Te = trep/40;                   % période d'échantillonnage

%poles discrets
z1 = exp(poleMref(1)*Te);
z2 = exp(poleMref(2)*Te);

Pc = [1 -z1-z2 z1*z2];          % polynôme caractéristique

G = c2d(TF2,Te,'zoh');          % procede discret
B = G.num{1}; 
A = G.den{1};

A0 = conv(A,[1 -1]);            % pour obtenir S0

[S0,R] = bezout(A0,B,0,Pc);     % permet de determiner S0 et R a partir de A0, B et Pc
S = conv([1 -1],S0);

z = tf('z');
S = tf(S,1,Te)/z^2;
R = tf(R,1,Te)/z^2;

T = polyval(Pc,1)/polyval(B,1); % on choisit T comme etant un gain : T=Pc(1)/B(1) 

%% Approche Frequentielle

Pref = Mref/(1-Mref);
% sisotool(Pref);

w0freq = wnref * sqrt(-2*ksi^2+sqrt(1+4*ksi^4));    % 3.89 rad/s - avec sisotool(Pref)
Mphase = 65.2;                                      % avec sisotool(Pref)
% valeurs de reference du systeme
% on doit trouver un correcteur pour arriver a ces valeurs

% --------------------------------------------------------------------
% Corriger la frequence (Correcteur PI)             % Pag: 124 Poly
a = w0freq/10;
% sisotool(TF2);                                    % pour regler le correcteur     % Control System - Turning Methos - Nichols Editor - Plot - C (integrateur et a) - varier le valeur de k jusqu'à w=w0

% On fait nichols de T2 en ajoutant un integrateur et zero en "-a", pour trouver le correcteur PI. On teste plusieurs choix pour 'a'
% On varie k du correcteur jusqu a avoir w=w0, et on note la marge de Phase
Mphase_2 = 118;

% --------------------------------------------------------------------
% Corriger la phase     (Correcteur AP)             % Pag: 126 Poly
phi = deg2rad(Mphase - Mphase_2); 
x = (1+sin(phi))/(1-sin(phi));

c = w0freq * sqrt(x);
b = w0freq / sqrt(x);

% --------------------------------------------------------------------
% Corriger le gain                                  % Pag: 126 Poly
% sisotool(TF2);                                    % Control System - Turning Methos - Nichols Editor - Plot - C (integrateur et a) - varier le valeur de k jusqu'à w=w0 (k=125.50)
k = 36.614;  


NumKpidI = k*conv([1 a],[1 b]);
DenKpidI = [1 c 0];
KPIDI = tf(NumKpidI, DenKpidI);


% Avec perturbation : d > 50% et tr = 4s
% Resultats terribles, on essaie avec autre valeur de 'a'

% --------------------------------------------------------------------
% Corriger la frequence (Correcteur PI) (2)         % Pag: 124 Poly
a = w0freq/1.5;
% sisotool(TF2);
Mphase_2 = 90.2;

% --------------------------------------------------------------------
% Corriger la phase     (Correcteur AP) (2)         % Pag: 126 Poly
phi = deg2rad(Mphase-Mphase_2);
x = (1+sin(phi))/(1-sin(phi));

c = w0freq * sqrt(x);
b = w0freq / sqrt(x);

% --------------------------------------------------------------------
% Corriger le gain                                  % Pag: 126 Poly 
% sisotool(TF2);                                    % Control System - Turning Methos - Nichols Editor - Plot - C (integrateur et a) - varier le valeur de k jusqu'à w=w0 (k=370)
k = 57.941;
  
NumKpidII = k*conv([1 a],[1 b]);
DenKpidII = [1 c 0];
KPIDII = tf(NumKpidII,DenKpidII);

