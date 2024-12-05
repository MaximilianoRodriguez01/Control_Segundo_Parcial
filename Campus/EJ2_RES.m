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

% Voy a declarar las variables de estado de la siguiente forma:
% x1 = theta, x2 = theta^*, x3 = v, x4 = v^*, x5 = v^**

syms x1 x2 x3 x4 x5 u y;
x = [x1; x2; x3; x4; x5];
y = x1;

wn = 2*pi*100;

f1 = x2;
f2 = sin(x1) + x3;
f3 = x4;
f4 = x5;
f5 = -3 * wn * x5 - 3 * wn^2 * x4 - wn^3 * x3 + u;
f = [f1; f2; f3; f4; f5];

% Los valores en el equilibrio:

x1e = pi / 6;
x2e = 0;
x3e = -1 / 2;
x4e = 0;
x5e = 0;
ue = wn^3 * x3e;
ye = x1e;

% Ahora, linealizo.

A = jacobian(f, x);
B = jacobian(f, u);
C = jacobian(y, x);
D = jacobian(y, u);

% Evaluo en los puntos de equilibrio.

A = subs(A, str2sym({'x1', 'x2', 'x3', 'x4', 'x5', 'u', 'y'}), {x1e, x2e, x3e, x4e, x5e, ue, ye});
B = subs(B, str2sym({'x1', 'x2', 'x3', 'x4', 'x5', 'u', 'y'}), {x1e, x2e, x3e, x4e, x5e, ue, ye});
C = subs(C, str2sym({'x1', 'x2', 'x3', 'x4', 'x5', 'u', 'y'}), {x1e, x2e, x3e, x4e, x5e, ue, ye});
D = subs(D, str2sym({'x1', 'x2', 'x3', 'x4', 'x5', 'u', 'y'}), {x1e, x2e, x3e, x4e, x5e, ue, ye});

% Puedo hacerlo asi tambien.
% A = double(subs(A,{x1, x2, x3, x4, x5, u, y},{x1e, x2e, x3e, x4e, x5e, ue, ye}));
% B = double(subs(B,{x1, x2, x3, x4, x5, u, y},{x1e, x2e, x3e, x4e, x5e, ue, ye}));
% C = double(subs(C,{x1, x2, x3, x4, x5, u, y},{x1e, x2e, x3e, x4e, x5e, ue, ye}));
% D = double(subs(D,{x1, x2, x3, x4, x5, u, y},{x1e, x2e, x3e, x4e, x5e, ue, ye}));

Ass = double(A);
Bss = double(B);
Css = double(C);
Dss = double(D);

% Transferencia
P = zpk(ss(Ass, Bss, Css, Dss));

% Ahora realizo la compensacion correspondiente. Para ello, defino el pap y
% el pmp, observando la expresión para la planta.

pap = zpk(-0.9306, 0.9306, 1);
pmp = zpk([], [-628.3, -628.3, -628.3, -0.9306, -0.9306], 1);

% Ahora, considerando que voy a implementar un PD con un Ts determinado, le
% asigno a este 5° (para obtener el menor Ts) y el pap deberá tener 25°, 
% pues se debe tener en cuenta que C * PMP tendrán 90° por la acción 
% integral del controlador y el Phase Margin será de 60° por consigna.

% fase = 2 * atan(p / wgc) ----> polo der y cero izq
% fase = 2 * atan(wgc / z) ----> polo izq y cero der

% Considerando la expresion: fase = 2 * atan(p / wgc), despejo wgc.
fase_pap = deg2rad(25);
wgc = 0.9306 / tan(fase_pap / 2);

% Tambien se puede sacar del bode, haciendo el bode del PAP y verificando
% el w para -25°.

figure();bode(pap, my_bode_options);title("PAP");
set(findall(gcf,'type','line'),'linewidth',2);

% Defino el PD considerando un retraso de media fase, con Ts. De esta
% forma, se puede obtener el valor de Ts con la siguiente ecuación.
% fase = 2 * atan(wgc / z) siendo z = 4/Ts.

fase_dig = deg2rad(5);
Ts = (4 / wgc) * tan(fase_dig / 2);

% Defino el PD (mejor ponerlo sin zpk y con la expresion general).

s = tf('s');
PD = (1 - (Ts/4)*s) / (1 + (Ts/4)*s);

% Defino el controlador considerando la acción integral y teniendo en
% cuenta que debe ser propio (igual cantidad de polos que de ceros).

C = 1 / (pmp * s * (s+10000)^4);

% Grafico el bode de PMP * C con K = 1 y determino la ganancia del controlador en 
% wgc = 4.2 rad/s

figure();bode(minreal(pmp * C), my_bode_options);title("PMP * C");
set(findall(gcf,'type','line'),'linewidth',2);

% Para wgc = 4.2 rad/s tenemos un K = -332 dB, por lo tanto:

K = db2mag(332);
C = K * C; 

% Defino L y verifico que se cumpla el margen de fase pedido.

L = minreal(pmp * C * pap * PD);

figure();bode(L, my_bode_options);title("L");
set(findall(gcf,'type','line'),'linewidth',2);

L_teo = minreal(P * C * PD);

figure();bode(L, my_bode_options);title("L teo");
set(findall(gcf,'type','line'),'linewidth',2);
