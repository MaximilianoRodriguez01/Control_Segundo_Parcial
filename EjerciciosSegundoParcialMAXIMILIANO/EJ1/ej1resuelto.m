clear all; close all; clc;

% Separo en parte de fase minima y parte pasatodo
Tau = 40e-3;
p = -0.4; %valor calculado para los Tau y Ts propuestos 
Pap1 = zpk(p,-p,1);
Pap2 = zpk(2/Tau,-2/Tau,-1);
Pmp = zpk([],p,1);

% Configuración del Bode
my_bode_options = bodeoptions;
my_bode_options.PhaseMatching = 'on';
my_bode_options.PhaseMatchingFreq = 1;
my_bode_options.PhaseMatchingValue = -180;
my_bode_options.Grid = 'on';
%my_bode_options.XLim = {[1 100]};

% Aproximacion de Pade para digitalizacion
Ts = 20e-3;
s = tf('s');
Pd = (1 - Ts/4 * s)/(1 + Ts/4 * s);

% Busco la minima Wgc, tal que el Pap reste 25° de fase, ya que
% los otros 5° los dejo para la parte correspondiente al control digital
Phase_dig = 5;    %Utilizo 5° para digitalizacion
wgc = (4/Ts)*tan(deg2rad(Phase_dig/2)); 
Phase_pap2 = rad2deg(2*atan(wgc/(2/Tau)));
Phase_pap1 = 25 - Phase_pap2;

%Calculo de p
p = -wgc*tan(deg2rad(Phase_pap1/2)); 

%Controlador
k = 1;
C = zpk(p,0,k); %Accion integral
figure();
bode(minreal(C*Pmp), my_bode_options);
title('Bode Lmp = C*Pmp con k=1');

%Ajusto ganancia para llevar wgc a 8.7 rad/seg 
k = db2mag(18.8);
C = zpk(p,0,k); %Accion integral
figure();
bode(minreal(C*Pmp*Pap1*Pap2*Pd), my_bode_options);
title('Bode L = C*P con C digitalizado');
