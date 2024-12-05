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

% Sé que la parte pasa todo debe tener 45°, entonces:

PM_pap = 45;

% Considerando que la fase del pasa todo sigue la siguiente expresion:
% 180° - N * fase_red = 45°, como tengo polos y ceros triples, multiplico
% por 3 y me queda: 540° - 6 * PM = 45°

fase_red = (-PM_pap + 540) / 6;

beta = (tan(deg2rad((fase_red + 90) / 2)))^2;

% Tengo que p = 1 / tau = 1, entonces tau = 1. Por lo tanto, defino z, tal
% que: 

tau = 1;
z = beta / tau;
wgc = sqrt(beta) / tau;

% De esta forma, defino la planta, PAP y PMP.

P = zpk([z z z], [1 1 1], -1);
pap = zpk([-1 -1 -1 z z z], [-z -z -z 1 1 1], 1);
pmp = zpk([-z -z -z], [-1 -1 -1], -1);

% Verifico que en la frecuencia dada el PAP presente un PM_pap de 45°.

figure();bode(pap, my_bode_options);title("PAP");
set(findall(gcf,'type','line'),'linewidth',2);

% Teniendo esto, puedo definir el controlador de forma tal que pmp * C
% tengan fase de atraso de 90°. Para ello, aplico acción integral.

s = tf('s');
C = zpk([-1 -1 -1], [0 -z -z -z], 1);

% Verifico con el bode que tengan 90° de fase.

figure();bode(minreal(pmp * C), my_bode_options);title("PMP * C");
set(findall(gcf,'type','line'),'linewidth',2);

% Ahora busco la ganancia en la frecuencia de corte dada wgc = 15.26 rad/s.
% K = -23.6 dB

K = db2mag(23.6);
C = zpk([-1 -1 -1], [0 -z -z -z], K);

figure();
subplot(1, 2, 1)
rlocus(minreal(pmp*C));
title("k positivo");
subplot(1, 2, 2);
rlocus(-minreal(pmp*C));
title("k negativo");

% De esta forma, se tiene el controlador solicitado, que dará un MP de 45°.

L = minreal(C * pmp * pap);

figure();
bode(L, my_bode_options);
title("L");
set(findall(gcf,'type','line'),'linewidth',2);

L_teo = minreal(P * C);

figure();
bode(L_teo, my_bode_options);
title("L teórica");
set(findall(gcf,'type','line'),'linewidth',2);

