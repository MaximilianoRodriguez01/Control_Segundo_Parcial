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
    
    wn = 2*pi*100;

%-----------------------------------

orden = 5;
x=sym('x',[orden 1],'real');
u=sym('u','real');

% Punto de equlibrio (x'=0)
u_e = 0.5* (wn)^3;
x_e = [pi/6 ; 0 ; -0.5 ; 0; 0]; %V0 lo encontre en la transferencia cuando la hice escrita

%vector de x punto
f1 = x(2);
f2 = sin(x(1))+ x(3);
f3 = x(4);
f4 = x(5);
f5 = u- (3*wn* x(5)+ 3*(wn)^2 * x(4)+ (wn)^3* x(3));


f = [f1;f2;f3;f4;f5];

%salida (podria haberle puesto y para no confundir)
g = x(1);


%Linealizo:


A = jacobian(f,x);
%la funcion subs cambia las ocurrencias de {x(1),x(2),x(3),x(4),u} por {x_e(1),x_e(2),x_e(3),x_e(4),u_e}
A = double(subs(A,{x(1),x(2),x(3),x(4),x(5),u},{x_e(1),x_e(2),x_e(3),x_e(4),x_e(5),u_e}));

B = jacobian(f,u);
B = double(subs(B,{x(1),x(2),x(3),x(4),x(5),u},{x_e(1),x_e(2),x_e(3),x_e(4),x_e(5),u_e}));

C = jacobian(g,x);
C = double(subs(C,{x(1),x(2),x(3),x(4),x(5),u},{x_e(1),x_e(2),x_e(3),x_e(4),x_e(5),u_e}));

D = jacobian(g,u);
D = double(subs(D,{x(1),x(2),x(3),x(4),x(5),u},{x_e(1),x_e(2),x_e(3),x_e(4),x_e(5),u_e}));

% Trasnferencia de la Planta Linealizada

%Puedo ponerla asi y me tira los polos como multiplicacion (mas visual)
P = zpk(ss(A,B,C,D));

%O puedo ponerla asi y me tira en polinomio
%P = tf(ss(A,B,C,D));


Avals=eig(A);

% figure(); hold on
% bode(P,optionss);



%P= 1/( (s-0.9306) (s+0.9306) (s+628.3)^3)

Pap_1= zpk((s+0.9306)/(s-0.9306));

% figure();
% bode(Pap_1, optionss);

%Para la digitalizacion tengo en cuenta el Pade, por lo que le dejo 5°
%Entonces el Pap_1 tiene fase 25° en w= 4.2



%De los 5° que dejé para la digitalización, calculo Ts para w=4.42
Ts= (4/4.2)* tan(deg2rad(2.5));

%OJO NUNCA PONER EL PADE CON ZPK PORQUE TE RETRASA 180°
Pd=(1-(Ts/4)*s)/(1+(Ts/4)*s);




Pmp=zpk(1 /((s+0.9306)^2 *(s+628.3)^3));


C= 1/(s*Pmp*(s+10000)^4);

Pmp_C=minreal( Pmp*C);
% figure();
% bode(Pmp_C, optionss);

K= db2mag(332);

C=C*K;

L= minreal (Pmp*C*Pap_1*Pd);
figure();
bode(L, optionss);




% Comentario sobre simulink:

% Para controlar la planta en el template de simulink, use un controlador
% PID y ajuste las ganancias a ojo tal de acomodar la respuesta a la
% referencia.
