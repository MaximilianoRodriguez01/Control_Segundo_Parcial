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

% Como la planta es de la forma dada, tengo la forma del controlador tal
% que PMP * C tenga un atraso de fase de aproximadamente 90° (110°). De esta forma, aplicando
% acción integral y un polo doble en 10.000, para que C sea propio, se
% tiene.

s = tf('s');
pmp = zpk([], [0 0 -80], 40);
C = zpk([-80 -80 -80], [0 -10000 -10000], 1/40);

% Le pongo un cero triple en -80 para contrarrestar el efecto del polo
% triple y que me de una fase de aproximadamente 90°.

% Teniendo esto, se tiene:

figure();bode(minreal(pmp * C), my_bode_options);title("PMP * C k = 1");
set(findall(gcf,'type','line'),'linewidth',2);
 
wgc = 900;

% Con K = -219 db, se tiene:

K = db2mag(219);
C = zpk([-80 -80 -80], [0 -10000 -10000], -K/40);

figure();bode(minreal(pmp * C), my_bode_options);title("PMP * C");
set(findall(gcf,'type','line'),'linewidth',2);

figure(); rlocus(minreal(pmp * C)); title("PMP * C, k > 0");

figure(); rlocus(-minreal(pmp * C)); title("PMP * C, k < 0");

% Por lo tanto, elijo k > 0.

% Se busca que pap tenga aproximadamente 10° de retraso de fase, ya que,
% pmp * C tendrán 110° de retraso de fase, el cual es lo mínimo.
% De esta manera, considerando que se trata de un sistema con un cero en el
% semiplano derecho y un polo en el semiplano izquierdo, se tiene:

% fase = 2 * atan (w / z) con z = 4/T;

fase_pap = deg2rad(20);

z = wgc / tan(fase_pap / 2);

T = 4 / z;

% Con esto, genero el PAP de forma tal que:

pap = zpk(4/T, -4/T, 1);

figure();bode(pap, my_bode_options);title("PAP");
set(findall(gcf,'type','line'),'linewidth',2);

% Por lo tanto, se define L y se grafica el bode para ver si se cumple con
% el MF pedido.

P = zpk(4/T, [0 0 -80 -4/T], 40);

L = minreal(pmp * C * pap);
figure();bode(L, my_bode_options);title("L");
set(findall(gcf,'type','line'),'linewidth',2);

L_teo = minreal(P * C);
figure();bode(L_teo, my_bode_options);title("L teorica");
set(findall(gcf,'type','line'),'linewidth',2);
