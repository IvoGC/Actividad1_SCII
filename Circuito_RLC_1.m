% Caso 1. Sistema de dos variables de Estado

A_u= 12;
Fs =1000000 ;
f  = 500;
T  =5*(1/f);
t = 0:1/Fs:T-1/Fs;
c = 50;
u = A_u*square(2*pi*f*t,c);
%con esto ya queda definida la entrada
%Caso 1.2
R = 4.7*10^3;%[ohms]
L = 10*10^-6;%[H]
C = 100*10^-9;%[F]

% [i, Vc] este es el vector X
% { X_p=AX+BU
% { Y = CX

A = [ -R/L , -1/L ; 1/C , 0 ];
B = [ 1/L ; 0 ];
C = [ R , 0 ] ;% Variando este yo selecciono la salida 
D = [ 0 ];

sys = ss(A, B, C, D);



lsim(sys, u, t);

% plot(t_s, y_s)




