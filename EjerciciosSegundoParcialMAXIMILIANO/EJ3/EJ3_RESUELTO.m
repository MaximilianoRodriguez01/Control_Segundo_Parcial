clear all; close all; clc;

% Separo en parte de fase minima y parte pasatodo
Pap1 = zpk(-1,1,1);
Pap2 = zpk(-10,10,1);
Pap = Pap1*Pap2;
Pmp = zpk([],[-1 -10],10);

% Configuración del Bode
my_bode_options = bodeoptions;
my_bode_options.PhaseMatching = 'on';
my_bode_options.PhaseMatchingFreq = 1;
my_bode_options.PhaseMatchingValue = -180;
my_bode_options.Grid = 'on';

% Busco a partir del bode de Pap la frecuencia en donde la parte pasa todo reste 25° de fase, ya que
% los otros 5° los dejo para la parte correspondiente al control digital
figure();
bode(Pap,my_bode_options);
title('Bode Pap');
wgc = 50;  %cumple las limitaciones fundamentales de Pap1 y Pap2

%Controlador
k = 1;
C1 = zpk([],0,k); %Accion integral
figure();
bode(minreal(C1*Pmp), my_bode_options);
title('Bode Lmp = C1*Pmp');

% - Pongo dos ceros en -5 para cumplir con el criterio de bode a w=50
% - Polo en -10000 para regularizar el controlador pero sin que influya en el margen de fase
C2 = zpk([-5 -5],[-10000 0],k); 
figure();
bode(minreal(C2*Pmp), my_bode_options);
title('Bode Lmp = C2*Pmp');

%Ajusto ganancia para llevar wgc a 50 aprox
k = db2mag(94.1);
C3 = zpk([-5 -5],[-10000 0],k); 
figure();
bode(minreal(C3*Pmp*Pap), my_bode_options);
title('Bode L = C3*Pmp*Pap');

% Aproximacion de Pade para digitalizacion
[Gm,Pm,Wcg,Wcp] = margin(C3*Pmp*Pap);
Phase_dig = 5;    %Utilizo 5° para digitalizacion
Ts = (2*Phase_dig*pi)/(180*Wcp);
s = tf('s');
Pd = (1 - Ts/4 * s)/(1 + Ts/4 * s);
figure();
bode(minreal(C3*Pmp*Pap*Pd), my_bode_options);
title('Bode L = C3*P con C digitalizado');

% Paso el controlador a discreto
C3_digital = c2d(C3, Ts, 'tustin');