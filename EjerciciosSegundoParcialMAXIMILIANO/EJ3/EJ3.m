clear all
close all


% Configuración del Bode
my_bode_options = bodeoptions;
my_bode_options.PhaseMatching = 'on';
my_bode_options.PhaseMatchingFreq = 1;
my_bode_options.PhaseMatchingValue = -180;
my_bode_options.Grid = 'on';
%my_bode_options.XLim = {[1 100]};


pap = zpk([-10 -1], [1 10], 1);

figure();
bode(pap, my_bode_options);
title("PAP");

wgc = 84;
fase_pap = deg2rad(-15);
Ts = tan(-fase_pap/2)*4/wgc;

% fase_pade_b = deg2rad(-5);
% Ts_b = tan(-fase_pade_b/2)*4/wgc_b;

s = tf('s');

pade = (1 - (Ts/4) * s)/(1 + (Ts/4) * s);

figure();
bode(pade, my_bode_options);
title("PADE");

pmp = zpk([], [-10 -1], 10);

c = 1/(pmp*s*(s+10000));

pmpc = minreal(pmp*c);

figure();
bode(pmpc, my_bode_options);
title("PMPC SIN AJUSTAR");

K = db2mag(119);

c = c*K;

pmpc = minreal(pmp*c);

figure();
bode(pmpc, my_bode_options);
title("PMPC AJUSTADO");

figure();
bode(minreal(pmpc*pap*pade), my_bode_options);
title("Bode compensado de L");

P = 10/((s-1)*(s-10));



figure();
bode(minreal(P*c*pade), my_bode_options);
title("bode de Planta*C");













