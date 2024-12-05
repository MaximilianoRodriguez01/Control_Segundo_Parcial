%Ej 1
clear all; 
close all;


% Opciones para el Bode
optionss=bodeoptions;
optionss.MagVisible='on';
optionss.PhaseMatching='on';
optionss.PhaseMatchingValue=-180;
optionss.PhaseMatchingFreq=1;
optionss.Grid='on';



%separamos en Pap y Pmp
s= tf ('s');

P= tf(10*(1/(s-10))*(1/(s-1)));

Pap1= tf((s+10)/(s-10));
Pap2= tf((s+1)/(s-1));

%Pmp * Pap = P (OJO AL CANCELAR)

Pmp = tf(10/((s+10)*(s+1)));


%parte all pass (Quiero que aporte 15°)--> VEMOS EN EL BODE QUE ESTO SE DA
%EN w=84

Pap3=Pap1*Pap2;

% figure()
% bode(Pap3, optionss);



% phi= 2*atan(10/50)+ 2*atan(1/50);
% phi= rad2deg(phi);





%retardo de tiempo discreto (Le voy a dejar 15°)

Ts= (4/84)* tan(deg2rad(7.5)); %me dio 6 ms el maximo Ts

Pap4=(1-(Ts/4)*s)/(1+(Ts/4)*s);

Pap=Pap3*Pap4;
% 
% figure()
% bode(Pap, optionss);



%Pmp= tf((s+10)*(s+1));


%Parte B

%controlador propio con accion integral-> polo en el origen y mayor denom que num
k=1;
C= k/s;
% 
% %quiero que Pmp * C me aporte 90° de fase, sumado a los 30 de Pap3 y Pap4,
% %son los 120°   
% 
% 
 Pmp_C= Pmp*C;
% % 
    % figure();
    % bode(Pmp_C);
% 
% 
% % Veo que en el w de interés, el conjunto Pmp-C difiere 170 grados de lo que necesito (en el Bode, la fase da -260). Necesito que sea de 90
% % Cancelo las singularidades de la Pmp y veo que me quedan igual cantidad
% % de polos que ceros, agrego un polo lejos para que sea propio el
% % controlador
% 
C=C*(1/Pmp)/(s+10000);
% 
Pmp_C= Pmp*C;
% figure();
% bode(Pmp_C);
% 
% %veo en el bode que en w=84 la ganancia es -118 dB por lo que la escalo
% %para aumentar esos 118 dB
k= db2mag(118);

C=C*k;
L= C*P*Pap4;

figure();
bode(L, optionss);






