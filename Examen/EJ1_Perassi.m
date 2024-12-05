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

% Defino la planta P.

P = zpk([], [1 10], 10);

% Se sabe que la planta se puede descomponer en una parte pasa todo (PAP),
% que tendrá un retraso de fase máximo de aproximadamente 30° (dado por
% consigna). Por otro lado, también habrá una parte de fase mínima (PMP),
% cuya fase más la del controlador tendrá un retraso de fase máximo de 90°.

% Además, se debe tener en cuenta que se agrega una parte de digitalización
% que estará dada por un padé de primer orden, que sigue la siguiente 
% expresión: pd = (1 - s * Ts / 4) / (1 + s * Ts / 4). Esta parte del muestreo
% se sumará con la fase correspondiente a la del pasa todo para resultar en un
% retraso de fase de 30°.

% De esta forma, se define la fase de la digitalización como 15° con el fin
% de que se obtenga la Ts máxima.

fase_dig = deg2rad(15);
fase_pap1 = deg2rad(30) - fase_dig;

% Además, se define el pap1 como:

pap1 = zpk([-10 -1], [1 10], 1);

% Se grafica el bode y se busca la frecuencia a la cual el pap1 tendrá una
% fase de 15°.

figure();bode(pap1, my_bode_options);title("PAP1");
set(findall(gcf,'type','line'),'linewidth',2);

% En wgc = 84.2 rad/s se tiene el retraso de 15°.

% De esta manera, teniendo en cuenta que el padé está dado por la ecuación
% vista anteriormente y considerando que es una transferencia con un polo
% en el semiplano izquierdo y un cero en el semiplano derecho, se cumple
% con: fase = 2 * atan(wgc / z) con z = 4 / Ts 

wgc = 84.2;
z = wgc / tan(fase_dig / 2);
Ts = 4 / z;

% ---- RTA DEL A: Ts = 0.0063; ----

% Defino el padé:

s = tf('s');
pd = (1 - s * Ts / 4) / (1 + s * Ts / 4);

% Verifico que en wgc = 84.2 rad/s se cumpla con la fase requerida.

figure();bode(pd, my_bode_options);title("PADE");
set(findall(gcf,'type','line'),'linewidth',2);

% Ahora, defino la parte de fase mínima, tal que:

pmp = zpk([], [-10 -1], 10);

% Considerando que PMP * C debe tener retraso de fase 90° y que el controlador debe
% tener fase integral, cancelo los polos del PMP con ceros en el
% controlador y agrego un polo en el cero (acción integral), el cual me va
% a aportar el retraso de fase que quiero. Además, para que el controlador
% sea propio, es decir, el grado del denominador sea mayor o igual al del
% numerador, agrego un polo en s = -10.000, que no me va a afectar el
% retraso de fase.

C = zpk(1 / (pmp * s * (s + 10000)));

% Ahora, verifico en el bode que el retraso de fase se mantenga en 90° y
% obtengo la ganancia en dB para wgc.

figure();bode(minreal(pmp * C), my_bode_options, {.1,100});title("PMP * C sin k definida");
set(findall(gcf,'type','line'),'linewidth',2);

% Para wgc = 84.2 rad/s, la ganancia es de -119 dB. Entonces:

K = db2mag(119);

C = K * C;

PMPC = minreal(pmp * C);

% Ahora, debo analizar el Root Locus con K > 0 y K < 0, de forma tal que
% sea estable.

figure();rlocus(PMPC);title("PMP * C, k > 0");
figure();rlocus(- PMPC);title("PMP * C, k < 0");

% Elijo K > 0, puesto que las raíces no adoptarán valores inestables.

% Habiendo analizado todo esto, verifico que el PM de L sea igual a 60°.

L = minreal(pd * pap1 * PMPC);

figure();bode(L, my_bode_options);title("L");
set(findall(gcf,'type','line'),'linewidth',2);

% Da el valor esperado de PM.

% Además, se verifica que L a partir de la división de la planta y L con la
% planta dada deben ser iguales.

L_teo = minreal(P * C * pd);
figure();bode(L, my_bode_options);title("L teórica");
set(findall(gcf,'type','line'),'linewidth',2);

% Luego, se tiene el grupo de las 4 transferencias.

S = 1 / (1 + L);
T = 1 - S;
PS = minreal(P*S);
CS = minreal(C*S);

% Por último, simulo respuesta al escalón de referencia y al de
% perturbación.

figure(); step(T); title('T');grid on
figure(); step(PS); title('PS');grid on

% Por completitud, dejo las otras respuestas.

figure(); step(S); title('S');grid on
figure(); step(PS); title('PS');grid on

