clear all; close all; clc;

% Configuración del Bode
my_bode_options = bodeoptions;
my_bode_options.PhaseMatching = 'on';
my_bode_options.PhaseMatchingFreq = 1;
my_bode_options.PhaseMatchingValue = -180;
my_bode_options.Grid = 'on';
%my_bode_options.XLim = {[1 100]};

Tau = 40e-3;
Ts = 20e-3;

Maxima_Fase_PAP = deg2rad(25);
Maxima_Fase_Pade = deg2rad(5);

wgc = tan((Maxima_Fase_Pade/2))*(4/Ts);

fase_pade = 25 - 2*rad2deg(atan((wgc/(2/Tau))));

P = tan(deg2rad(fase_pade/2))*(wgc);

s = tf('s');

Pmp = zpk([],-P,1);

Pd = zpk(200, -200, 1);

Pap = zpk([-P 50], [-50 P], 1);

figure

C1 = (s+P)/s;

PMPC_sin_ganancia = Pmp*C1;

C = 8.810*(s+P)/s;

Pmpc = Pmp * C;

figure();
bode(Pmp*C, my_bode_options);
title("Pmpc con ganancia");

L = minreal (Pmpc * Pap * Pd);



figure();
bode(L, my_bode_options);
title("Bode de L");


T = feedback(L, 1);

eig(T);


figure();
bode(minreal(T), my_bode_options);
title("Bode de T");



