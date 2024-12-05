%{
+==========================================
+
+                Fórmulas Matlab
+
+==========================================
%}

% =========================================

P = PMP * PAP;

% PAP se encuentra teniendo en cuenta los polos y ceros en el semiplano derecho. 
% Se espejan sus contrarios en el semiplano izquierdo (polos a ceros y ceros a polos).
% PAP generalmente tiene fase 30° (cuando se pide un PM de 60°).

% PMP se ve a partir de la diferencia entre PAP y PMP.

% =========================================

% Para que C tenga acción integral, debe haber un polo en 0.
% Debe ser propio, es decir, el grado del denominador > al del numerador.

% Generalmente, se cancelan los polos o ceros del PMP, para que la multiplicación PMP * C tenga fase 90°.
% Para que sea propio, si hay más ceros que polos en el C, se agregan polos en frecuencias altas 1e3, 10e3, etc.

% Si el controlador tiene que ser PI, entonces tendrá un cero acompañado de acción integral.

% =========================================

% Cuando se pide un Ts o se pide una discretización del controlador,
% se debe crear un padé. Este tiene fase entre 5° y 15° y se tiene que sumar a la fase de PAP.
% Sigue la siguiente expresión:

pade = (1 - (Ts/4) * s) / (1 + (Ts/4) * s);

% pade = (s - 4 / Ts) / (s + 4 / Ts);
% No hay que hacerlo con zpk, no sé por qué no se puede.

% =========================================

% Si se tiene una planta con retardo, entonces la expresión será similar al del padé.

retardo = (1 - (tau/2) * s) / (1 + (tau/2) * s);

% retardo = (s - 2 / tau) / (s + 2 / tau);
% tau no es lo mismo que la tau de las redes de adleanto y atraso.

% =========================================

% Cuando tengo un cero del lado izquierdo y un polo del lado derecho, puedo utilizar:
% fase = 2 * atan(p / wgc). 

% Cuando tengo un polo del lado izquierdo y un cero del lado derecho, puedo utilizar:
% fase = 2 * atan(wgc / z).

% Esto resulta útil si por ejemplo se necesita Ts del padé y se tiene la fase y el wgc.

% =========================================

% La fase de PAP está dada por:

% PM_PAP = 180° - N * Fase_Red, 
% donde N es la cantidad de polos y ceros en el semiplano izquierdo o derecho (ver uno de los 2).

%La Fase_Red está definida por la red de adelanto (cero más cerca del origen) / atraso (polo más cerca del origen), 
% que es la siguiente:

% Fase_Red = 2 * atan(sqrt(beta)) - 90°.
% donde beta es la distancia del polo al cero.

% Red de adelanto: polo = beta / tau, zero = 1 / tau.
% Red de atraso: polo = 1 / tau, zero = beta / tau.

% =========================================

%Para hacer un controlador discretizado (para simulink):

c2d(C, Ts, 'tustin');

%==========================================

% Grupo de las 4:

T = feedback(L,1); 
S = 1 - T;   
CS = minreal(C*S);
PS = minreal(P*S);

% Simulo respuesta al escalon de referencia:
figure(); step(T);

% Simulo respuesta al escalon de perturbacion:
figure(); step(PS);
