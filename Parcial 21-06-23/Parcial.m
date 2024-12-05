clear all;
close all;
clc;

% Defino las constantes.

J = 1;
MGL2 = 1;

% Defino x1 como theta y x2 como theta^*

syms x1 x2 y u;

x = [x1; x2];

f1 = x2;
f2 = (u / J) - (MGL2 / J) * sin(x1);
f = [f1; f2];

y = x1;

% Obtengo ue en el equilbrio, ya que se los valores de x1e y x2e.

x1e = deg2rad(150);
x2e = 0;
ue = sin(x1e);
ye = x1e;

% Linealizo.

A = jacobian(f, x);
B = jacobian(f, u);
C = jacobian(y, x);
D = jacobian(y, u);

% Evaluo en el equilibrio.

A = subs(A, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye});
B = subs(B, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye});
C = subs(C, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye});
D = subs(D, str2sym({'x1', 'x2', 'u', 'y'}), {x1e, x2e, ue, ye});

Ass = double(A);
Bss = double(B);
Css = double(C);
Dss = double(D);

% Hago la transferencia:

P = zpk(ss(Ass, Bss, Css, Dss));

% Configuración de los bodes
my_bode_options = bodeoptions;
my_bode_options.MagVisible='on';
my_bode_options.PhaseMatching = 'on';
my_bode_options.PhaseMatchingFreq = 1;
my_bode_options.PhaseMatchingValue = -180;
my_bode_options.Grid = 'on';

% Voy a separar la planta en su parte de fase mínima y la parte de pasa
% todo. Considerando que las fases de retraso del PMP * C = 90° y la parte
% del pasa todo debe ser de 30°, ya que, el PM debe ser de 60°.

% Además, tengo en cuenta que como el PAP tendrá una parte correspondiente
% a la digitalización, a la cual le asigno 5°. Entonces, tengo que:

fase_dig = deg2rad(5);
fase_pap = deg2rad(30) - fase_dig;

pap = zpk(-0.9306, 0.9306, 1);

% Como tengo un polo en la derecha y un cero en el semi plano izquierdo,
% puedo obtener el valor de wgc a partir de la siguiente expresión.

% fase = 2 * atan(p / wgc)

p = 0.9306;
wgc = p / tan(fase_pap / 2);

% Verifico el valor de wgc con el bode.

figure();bode(pap, my_bode_options);title("PAP");
set(findall(gcf,'type','line'),'linewidth',2);

% Hay una paridad entre lo sacado a partir de la expresión y lo visto en el
% bode.

% Ahora, con ese valor de wgc obtengo el valor de Ts para una fase de 5°.
% fase = 2 * atan(wgc / z) siendo z = 4 / Ts.

z = wgc / tan(fase_dig / 2);
Ts = 4 / z;

% Defino la expresión del padé de primer orden.

s = tf('s');
pd = (1 - Ts / 4 * s) / (1 + Ts / 4 * s);

% Por otro lado, defino la parte de fase minima y el controlador con
% acción integral, considerando que la suma de sus fases es 90°.

pmp = zpk([], [-0.9306 -0.9306], 1);
C = zpk([-0.9306 -0.9306], [0 -10000], 1);

% Grafico el bode y analizo que la fase sea de 90° y la magnitud en la frec
% dada.

figure();bode(minreal(pmp * C), my_bode_options);title("PMP * C");
set(findall(gcf,'type','line'),'linewidth',2);

% La magnitud dada es la siguiente K = -92.5 dB. Entonces, el controlador:

K = db2mag(92.5);
C = zpk([-0.9306 -0.9306], [0 -10000], K);

% Ahora debo analizar si el signo será positivo o negativo 
% a partir del Rooot Locus.

figure();rlocus(minreal(pmp * C));title("PMP * C, k > 0");
figure();rlocus(- minreal(pmp * C));title("PMP * C, k < 0");

% K me va a dar estable cuando es mayor a 0. De esta forma, compruebo con
% el bode que el margen de fase sea el adecuado.

figure();bode(minreal(pmp * C), my_bode_options);title("PMP * C");
set(findall(gcf,'type','line'),'linewidth',2);

% Ahora, calculo L con todo el análisis llevado a cabo.

L = minreal(pd * pap * pmp * C);
figure();bode(L, my_bode_options);title("L");
set(findall(gcf,'type','line'),'linewidth',2);

L_teo = minreal(pd * P * C);
figure();bode(L_teo, my_bode_options);title("L teórica");
set(findall(gcf,'type','line'),'linewidth',2);

% Ahora, para discretizarlo para el Simulink utilizo la siguiente función.

C_dig = c2d(C, Ts, 'tustin');

