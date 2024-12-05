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

% Declaro el sistema en espacio de estados

syms x1 x2 y u;

x = [x1; x2];
y = x2;

f1 = -sqrt(x1 - x2) + u;
f2 = sqrt(x1 - x2) - sqrt(x2);
f = [f1; f2];

% Defino los puntos de equilibrio.

x1e = 2;
x2e = 1;
ue = 1;
ye = 1;

% Obtengo la linealización jacobiana y evaluo en los puntos de equilibrio.
 
A = jacobian(f, x);
B = jacobian(f, u);
C = jacobian(y, x);
D = jacobian(y, u);

A = double(subs(A, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye}));
B = double(subs(B, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye}));
C = double(subs(C, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye}));
D = double(subs(D, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye}));

% Tengo la siguiente transferencia.

P = zpk(ss(A, B, C, D));

% Agrego acción integral al controlador con un polo en 0. Además, agrego un
% cero en -0.191 y dejo un polo de P en la multiplicación de L = P * C para
% centrar C con una ganancia determinada en una frecuencia dada para que L
% nos de un margen de fase de 65° (pues 5° se van a restar con el padé). 

C = zpk(-0.191, [-10e3 0], 1);

% wgc = 0.621 rad/s --> k = -85.1 dB
% wgc = 0.087 rad/s --> K = -51.3 dB

k = db2mag(85.1);
C = zpk(-0.191, [-10e3 0], k);

% Ahora, hago el bode de L y busco la frecuencia a la cual se producen 65° para
% cumplir con el margen de fase.

figure();bode(minreal(C * P), my_bode_options);title("P * C");
set(findall(gcf,'type','line'),'linewidth',2);

% Ahora, busco Ts tal que la fase del padé dé 5° con la siguiente fórmula.
% fase = 2 * atan(wgc / z) siendo z = 4 / Ts

fase_dig = deg2rad(5);
wgc = 0.621;

s = tf('s');
Ts = (4 / wgc) * tan(fase_dig / 2);

PD = (1 - (Ts/4) * s) / (1 + (Ts/4) * s);

% Entonces, verifico el margen de ganancia con el bode de L.

L = minreal(P * C * PD);
figure();bode(L, my_bode_options);title("L");
set(findall(gcf,'type','line'),'linewidth',2);

% Discretizo el controlador para el Simulink
C_digital = c2d(C, Ts, 'tustin');
