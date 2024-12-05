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

% Defino constantes.

a = 0.25 * (6 / pi)^2;
p = 20;

% Defino las variables de estado, la entrada y la salida.
% Se debe tener en cuenta que x1 es theta, x2 es theta punto y x3 es tau.

syms x1 x2 x3 u y;

x = [x1; x2; x3];
y = x1;

% Declaro las ecuaciones del espacio de estados, tal que:

f1 = x2;
f2 = -a * abs(x1) * x1 + sin(x1) - abs(x2) * x2 + x3;
f3 = p * u - p * x3;

% Debo aclarar que voy a trabajar sobre la parte positiva de la funcion
% x*|x| de forma tal que, viendo el gráfico auxiliar lo puedo aproximar
% como una función cuadrática. Por lo tanto, se tiene lo siguiente.

f2 = -a * x1^2 + sin(x1) - x2^2 + x3;
f = [f1; f2; f3];

% Busco los valores en el equilibrio considerando que en el punto de 
% equilibrio f = 0. Además, elijo theta = pi / 6.

x1e = pi / 6;
x2e = 0;
x3e = a * abs(x1e) * x1e - sin(x1e);
ue = x3e;
ye = x1e; 

% RTA A: ue = +-0.25

% Ahora, para linealizar aplico el jacobiano y evaluo en el punto de
% equilibrio dado.

A = jacobian(f, x);
B = jacobian(f, u);
C = jacobian(y, x);
D = jacobian(y, u);

Ass = double(subs(A, str2sym({'x1', 'x2', 'x3', 'u', 'y'}), {x1e, x2e, x3e, ue, ye}));
Bss = double(subs(B, str2sym({'x1', 'x2', 'x3', 'u', 'y'}), {x1e, x2e, x3e, ue, ye}));
Css = double(subs(C, str2sym({'x1', 'x2', 'x3', 'u', 'y'}), {x1e, x2e, x3e, ue, ye}));
Dss = double(subs(D, str2sym({'x1', 'x2', 'x3', 'u', 'y'}), {x1e, x2e, x3e, ue, ye}));

% Tengo la siguiente transferencia.

P = zpk(ss(Ass, Bss, Css, Dss));

% Al igual que en el ejercicio anterior, divido la planta en una parte pasa
% todo y una parte de fase mínima, considerando que la expresión del PAP 
% está dada por consigna.

% De esta manera, se busca un MP de 60° para L de forma tal que PMP * C
% tenga un retraso de fase de 90° y el PAP tenga un retraso de 30°.

% A fin de encontrar un retraso de fase de 90°, diseño el controlador
% con acción integral y cancelando los polos complejos conjugados de la
% planta en su módulo con un cero doble.

% Además, agrego un polo en 10.000 para que el controlador sea propio.

% Busco la frecuencia wgc en la cual se den 90° en el bode por lo 
% mencionado anteriormente.

pmp = zpk([], [0.0889i -0.0889i], 1);
C = zpk([-0.0889 -0.0889], [0 -10000], 1);

figure();bode(minreal(pmp * C), my_bode_options);title("P * C sin k definida");
set(findall(gcf,'type','line'),'linewidth',2);

% wgc = 20 rad / s y K = -138 dB
% wgc = 18 rad/s y K = -105 dB.

K = db2mag(138);
C = K * C;
PMPC = minreal(pmp * C);

figure();bode(PMPC, my_bode_options);title("PMP * C con k definida");
set(findall(gcf,'type','line'),'linewidth',2);

% Analizo el Root Locus para ver si tomo K positiva o negativa.

figure();rlocus(PMPC);title("PMP * C, k > 0");
figure();rlocus(- PMPC);title("PMP * C, k < 0");

% Elijo K > 0, puesto que para esos valores las raíces son estables.

% Por otro lado, para el PAP considerando que es una transferencia con un polo 
% del SPI y un cero del SPD. Se tiene la siguiente ecuación de donde se
% puede obtener el valor del cero z para el retraso de fase de 30°.
% fase = 2 * atan(wgc / z) con z = 4 / Ts.

wgc = 20;
fase_pap = rad2deg(30);
z = wgc / tan(fase_pap / 2);
Ts = 4 / z;

% Defino la función para PAP.

s = tf('s');
pap = (1 - s * Ts / 4) / (1 + s * Ts / 4) * (1 / 40);

% Verifico con el bode que en la frec dada me haya dado el valor esperado.

figure();bode(pap, my_bode_options);title("PAP");
set(findall(gcf,'type','line'),'linewidth',2);

% De esta forma, verifico el valor del PM de L, tal que:

L = minreal(pmp * C * pap);

figure();bode(L, my_bode_options);title("L");
set(findall(gcf,'type','line'),'linewidth',2);

% Para Simulink discretizo el controlador de la siguiente forma.
Ts = 0.93;
P_dig = c2d(C, Ts, 'tustin');
