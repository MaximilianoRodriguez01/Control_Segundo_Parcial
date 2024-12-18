clear all
close all

% Configuración de los bodes
my_bode_options = bodeoptions;
my_bode_options.MagVisible='on';
my_bode_options.PhaseMatching = 'on';
my_bode_options.PhaseMatchingFreq = 1;
my_bode_options.PhaseMatchingValue = -180;
my_bode_options.Grid = 'on';

beta = (tan(deg2rad(495/6+90)*(1/2))^2);
tau = 1;
z = beta/tau;
wgc = sqrt(beta)/tau;

s = tf('s');

P = zpk([z z z], [1 1 1], -1);

pap = zpk([z z z -1 -1 -1], [1 1 1 -z -z -z], -1);

pmp = zpk([-z -z -z], [-1 -1 -1], 1);

c = zpk([-1 -1 -1], [0 -z -z -z], 1);

pmpc = minreal(pmp*c);

figure();
bode(pap, my_bode_options);
title("pap");

figure();
bode(pmpc, my_bode_options);
title("pmpc sin ajustar");

k = db2mag(23.6);

c = c*k;
pmpc = minreal(pmp*c);

figure();
subplot(1, 2, 1)
rlocus(pmpc);
title("k positivo");
subplot(1, 2, 2);
rlocus(-pmpc);
title("k negativo");

figure();
bode(pmpc, my_bode_options);
title("pmpc");

L = minreal(pmpc * pap);

figure();
bode(L, my_bode_options);
title("L compensado");

figure();
bode(minreal(P*c), my_bode_options);
title("Planta x C");







