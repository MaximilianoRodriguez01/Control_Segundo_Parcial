clear all;
close all;
clc;

% Configuración de los bodes
my_bode_options = bodeoptions;
my_bode_options.MagVisible='on';
my_bode_options.PhaseMatching = 'on';
my_bode_options.PhaseMatchingFreq = 1;
my_bode_options.PhaseMatchingValue = -180;
my_bode_options.Grid = 'on';

% Defino las variables del problema

tau = 40e-3;
Ts = 20e-3;
s = tf('s');
p = 0.4;

% Defino la planta, los paps y el pmp

P = zpk(2/tau, [-2/tau, p], 1);
pap1 = zpk(-p, p, 1);
pap2 = zpk(50, -50, 1);
pap = pap1 * pap2;

% Quiero que la digitalización tenga 5° de fase

PD = zpk(4/Ts, -4/Ts, 1);
fase_dig = deg2rad(5);

% Teniendo en cuenta que se debe cumplir que: 
% fase_des = 2*arctg(wgc / z) siendo z = Ts/4

wgc = (4/Ts) * (tan(fase_dig / 2));

% Con este valor, puedo calcular la fase del pap2, tal que:

fase_pap2 = 2 * rad2deg(atan(wgc / (2 / tau)));

% Ya que, este es el pap que está completo, sin considerar p.
% De esta manera, teniendo en cuenta que en fase los paps deben sumar 25,
% pues 5 se lleva el PD.

fase_pap1 = 25 - fase_pap2;

% Con esta fase, puedo obtener el valor de p, siguiendo la siguiente
% ecuacion: fase = 2 atan (p / wgc).

p = wgc * tan(deg2rad(fase_pap1 / 2)); %p = 0.4

% Teniendo esto, ahora voy a hcer el punto B

% Defino pmp

pmp = zpk([], p, 1);

% Defino el controlador teniendo en cuenta que en conjunto con el pmp deben
% tener fase 90° y considerando que debe ser propio. Entonces, cancelo el
% polo en p de pmp y agrego un polo en s = 0, asi se tiene la fase buscada.

k = 1;
C = zpk(p, 0, k);

% Como wgc es 8.7 rad/s, debo buscar la ganancia en esa frecuencia.

figure();bode(minreal(pmp * C), my_bode_options);title("PMP * C");
set(findall(gcf,'type','line'),'linewidth',2);

% Se tiene k = -18.8. Por ende.

k = db2mag(18.8);
C = zpk(p, 0, k);

% De esta manera, se grafica el bode de L para ver si cumple con el margen de
% fase requerido:

L = minreal(pmp * C * pap * PD);

figure();bode(L, my_bode_options);title("L");
set(findall(gcf,'type','line'),'linewidth',2);

L_teo = minreal(P * C * PD);
figure();bode(L, my_bode_options);title("L teórica");
set(findall(gcf,'type','line'),'linewidth',2);

T = feedback(L, 1);

figure();bode(T, my_bode_options);title("T");
set(findall(gcf,'type','line'),'linewidth',2);

A = eig(T);
