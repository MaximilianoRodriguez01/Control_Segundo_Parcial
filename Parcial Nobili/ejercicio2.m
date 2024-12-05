%{
+==========================================
+
+                Ejercicio 2
+
+==========================================
%}

    close all
    clear all
    clc
    
    
    % Config:
    s = tf('s');
    
    optionss=bodeoptions;
    optionss.MagVisible='on';
    optionss.PhaseMatching='on';
    optionss.PhaseMatchingValue=-180;
    optionss.PhaseMatchingFreq=1;
    optionss.Grid='on';
    
    %Constantes:
    
    v0 = 0.1;
    psi = 0.707;
    wn = 2*pi*100;
    alfa = 32 * sqrt(2)/pi^3;

%-----------------------------------

orden = 4;
x=sym('x',[orden 1],'real');
u=sym('u','real');

% Punto de equlibrio (x'=0)
u_e = v0;
x_e = [pi/4 ; 0 ; v0 ; 0];

%vector de x punto
f1 = x(2);
f2 = -alfa * x(1)^3 + sin(x(1)) - x(2)^3 + (x(3) - v0);
f3 = x(4);
f4 = wn^2 * u - wn^2 *x(3) - 2*psi*wn*x(4);

f = [f1;f2;f3;f4];

%salida
g = x(1);

A = jacobian(f,x);
%la funcion subs cambia las ocurrencias de {x(1),x(2),x(3),x(4),u} por {x_e(1),x_e(2),x_e(3),x_e(4),u_e}
A = double(subs(A,{x(1),x(2),x(3),x(4),u},{x_e(1),x_e(2),x_e(3),x_e(4),u_e}));

B = jacobian(f,u);
B = double(subs(B,{x(1),x(2),x(3),x(4),u},{x_e(1),x_e(2),x_e(3),x_e(4),u_e}));

C = jacobian(g,x);
C = double(subs(C,{x(1),x(2),x(3),x(4),u},{x_e(1),x_e(2),x_e(3),x_e(4),u_e}));

D = jacobian(g,u);
D = double(subs(D,{x(1),x(2),x(3),x(4),u},{x_e(1),x_e(2),x_e(3),x_e(4),u_e}));

% Trasnferencia de la Planta Linealizada
P = tf(ss(A,B,C,D));

Avals=eig(A);

figure(); hold on
bode(P,optionss);

% Comentario sobre simulink:

% Para controlar la planta en el template de simulink, use un controlador
% PID y ajuste las ganancias a ojo tal de acomodar la respuesta a la
% referencia.
