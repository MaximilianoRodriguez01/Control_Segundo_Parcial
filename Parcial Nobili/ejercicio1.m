
%{
+==========================================
+
+                Parte a
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


% Transferencia de la Planta Nominal:

P_0 = 10 /((s-10)*(s-1));

% El efecto del retraso de media muestra del muestro se hace a partir de
% una aproximacion PADE de primer orden:
    
% P = P_0 * (1 - s * (Ts/4)) / (1 - s * (Ts/4));  

%{
    Para analizar cual es el mayor sampling rate de control en tiempo
    discreto para que sea viable una compensacion con accion integral
    utilizando loopshaping, separo la planta en parte de fase minima y
    parte pasa todo. A su vez, separo la parte pasa todo en dos, una que
    tiene en cuenta el efecto del cero de fase no minima impuesto por el
    muestreo y otra que integra los polos inestables de la planta nominal:

    Pap_z = (1 - s * (Ts/4)) / (1 - s * (Ts/4))
    Pap_p = (s+10) / (s-10)   *   (s+1)/(s-1)

    ->Pap = P_ap_z * P_ap_p

    ->Pmp = 10 / ((s+10)*(s+1))


    Para que la planta se pueda estabilizar con un margen de fase de 60
    grados, el retraso de fase admisible para la parte pasa todo es de 30
    grados (como dice el enunciado). Es decir, los retrasos de fase de
    P_ap_p y P_ap_z deben sumar 30 grados entre ambos. 

    El minimo valor de Ts tal que a una dada frecuencia existe retraso de 
    fase menor o igual a 30 grados se dara cuando cada parte pasa todo
    aporte 15 grados de retraso. Analizando el bode de P_ap_p encuentro el
    valor de w donde su fase vale -15deg.
%}

Pap_p = ((s+10)/(s-10)) *((s+1)/(s-1));
optionss.MagVisible='off';
figure(); bode(Pap_p,optionss);title('Pap_p');
optionss.MagVisible='on';


%{
    Se tiene fase -15 deg en Pap_p para:    w = 84 rad/s

    Quiero que a dicha frecuencia la transferencia Pap_z aporte otros 15
    grados de retraso de fase:
    

    phi = 2 * arctan(w/z) ------> z = w / tan(phi/2) 

    ...
    
    z = 4/Ts    y   phi = 15

    ...

    z = 638 --> Ts =(aprox) 0.006

%}
optionss.MagVisible='off';

Ts = 0.006 ; 
Pap_z = (1 - s * (Ts/4)) / (1 + s * (Ts/4));
figure(); bode(Pap_z,optionss); title('Pap_z');


%Entonces:

Pap = minreal(Pap_z * Pap_p);
figure(); bode(Pap,Pap_z,Pap_p,optionss); title('Pap'); legend('Pap','Pap_z','Pap_p');

optionss.MagVisible='on';


%{
    Se observa en el bode de fase de Pap que esta alcanza un maximo en w =
    84 rad/s de -30. Es decir que el menor retraso de fase de la
    transferencia Pap es de 30 grados. 

    En conclusion, el mayor sampling time posible tal que se tenga un  un 
    retraso de fase admisible para la parte pasatodo de no más de 30 
    grados es Ts = 0.0063. 

    Para valores menor de Ts, se pueden obtener menores retrasos de fase ya
    que los polos y cero fase no minima se separan aun man. Recordemos que
    los polos inestables imponian un minimo en la frecuencia de cruce por
    cero y los ceros de fase no minima, un maximo.
%}



%%

%{
+==========================================
+
+                Parte b
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


% Transferencia de la Planta Nominal:

P_0 = 10 /((s-10)*(s-1));

%Termino de definir la planta con la tasa de muestreo encontrada:

Ts = 0.006; 
[num,den] = pade(Ts/2,1);
Pd = tf(num,den);

P = minreal(P_0 * Pd);


% Nuevamente separo en parte pasa todo y de minima fase. Como mi controlador
% debe tener accion integral, agrego un integrador a Pmp:

Pmp = 10 / ((s+10)*(s+1));
Pmp_monio = Pmp/s;

Pap_p =((s+10)/(s-10)) *((s+1)/(s-1));
Pap_z = Pd;
Pap = minreal(Pap_p * Pap_z);

figure(); bode(Pap,Pap_z,Pap_p,optionss); title('Pap'); legend('Pap','Pap_z','Pap_p');
figure(); bode(Pmp,Pmp_monio,optionss); title('Pmp'); legend('Pmp','Pmp_{monio}');

%{
    Como se analizo en el ejercicio anterior, debo utilizar como frecuencia
    de cruce por cero aproximadamente w = 84 ya que para dicha frecuencia la parta pasa
    todo cumple con el retraso de fase admisible. Si el retraso de fase
    fuese mayor, no se podria estabilizar la planta con un margen de fase de
    60 grados.

    Teniendo eso en cuenta planteo el controlador con accion integral. Un
    diseno simple consiste en cancelar la parte de minima fase y agregar un
    polo regularizador alejado con tal de que quede un controlador propio y
    no afecte la fase cerca de la frecuencia de cruce por cero.

    Una vez armado el controlador, ajusto la ganancia tal que la frecuencia
    de cruce por cero sea la deseada.
%}

C_monio = 1/Pmp * (1/(s+10000)); %figure(); bode(minreal(Pmp_monio * C_monio),optionss);
k = db2mag(118.45);
C = k*C_monio/s;

% Transferencia de Lazo Abierto:
L = minreal(P*C); figure(); bode(L,optionss); title('L: transferencia de lazo abierto');

% Planteo las transferencias de interes del sistema y grafico su respuesta
% en frecuencia:

T = feedback(L,1); figure(); bode(T,optionss);title('T');
S = 1 - T;         figure(); bode(S,optionss);title('S');
CS = minreal(C*S); figure(); bode(CS,optionss);title('CS');
PS = minreal(P*S); figure(); bode(PS,optionss);title('PS');

% Simulo respuesta al escalon de referencia:
figure(); step(T);

% Simulo respuesta al escalon de perturbacion:
figure(); step(PS);
