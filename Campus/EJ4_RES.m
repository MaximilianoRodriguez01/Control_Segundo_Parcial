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

% Defino la planta dada.

P = zpk([], [10 1], 10);

% Separo en parte pasa todo y parte de fase mínima, tal que:

pap = zpk([-10 -1], [10 1], 1);
pmp = zpk([], [-10 -1], 10);

% Tengo en cuenta que para obtener el mayor sampling rate, voy a usar una
% fase de 15° del PD. De esta manera, teniendo en cuenta que entre el PD y
% el PAP deben atrasar 30°, se tiene que:

fase_dig = deg2rad(15);
fase_pap = deg2rad(30) - fase_dig;

% Ahora, puedo obtener la frecuencia de corte de forma tal que en el bode
% del pap se vea donde resta 15°.

figure();bode(pap, my_bode_options);title("PAP");
set(findall(gcf,'type','line'),'linewidth',2);

% En wgc = 83.6 rad / s, se tiene que el pap resta 15°.

wgc = 83.6;

% Teniendo wgc y la fase del PD, se puede obtener Ts de la siguiente forma:
% fase = 2 * atan(w / z) siendo z = 4 / Ts.

Ts = (4 / wgc) * tan(fase_dig / 2);

% De esta manera, se tiene la expresión del padé.

s = tf('s');
PD = (1 - (Ts/4) * s) / (1 + (Ts/4) * s);

% Como se pide que se debe cumplir que el margen de fase sea de 60°, PAP*PD
% debe resultar en 30° puesto que PMP*C es de 90° dadas las características
% propias de un controlador con acción integral.

% De esta manera, defino el controlador de la siguiente forma:

C = 1 / (s * pmp);

% Es un controlador propio, pues el orden del denominador es mayor al del
% numerador.

% Ahora, se busca la ganancia del compensador en wgc.

figure();bode(minreal(pmp * C), my_bode_options);title("PMP * C");
set(findall(gcf,'type','line'),'linewidth',2);

% En esta frecuencia K = -38.4 dB (de la ganancia de PMP). Entonces:

K = db2mag(38.4);
C = K * (1 / (s * pmp));

% De esta manera, se tiene el compensador diseñado y se verifica que el
% margen de fase a lazo abierto sea el deseado.

L = minreal(PD * pap * pmp * C);

figure();bode(L, my_bode_options);title("L");
set(findall(gcf,'type','line'),'linewidth',2);

L_teo = minreal(P * C * PD);
figure();bode(L_teo, my_bode_options);title("L teorica");
set(findall(gcf,'type','line'),'linewidth',2);


