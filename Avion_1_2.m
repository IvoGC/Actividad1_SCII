close all; clear all ; clc;
% Caso 1. Sistema lineal de cuatro variables de estado
% El sistema linealizado  en el punto de trabajo [0,0,0,0] obtenemos el
% siguiente sistema




%A = [ -a,a,0,0 ; 0,0,0,0 ; w^2,-w^2,0,0 ; c,0,0,0 ];
%B = [ 0; 0 ; w^2*b ; 0 ];
%C = [ 1,0,0,0 ; 0,0,0,1 ];
%D = 0;

%sys = ss(A, B, C, D);

% alfa_p  = a(phi-alfa)
% phi_pp  = -w^2(phi-alfa-b*U)
% h_p     = c*alfa
%equivaent
% x1_p = a(x2-x1)
% x3_P = -w^2(x2-x1-b*u)
% x4_p = c*x1

dt = 1e-3;  % step size
w=2;
a=0.05;
b=5;
c=50;

t_0=0;
t_f=1; %tiempo total de simulacion

t = t_0:dt:(t_f-dt);

y =zeros(t_f/dt,4);

%alfa(0)=0;
%phi(0)=0;
%phi_p(0)=0;
%h(0)=0;

y(1,1) = 0; %alfa
y(1,2) = 0; %phi
y(1,3) = 0; %phi_p
y(1,4) = 0; %h
u = ones(t_f/dt,1);
u = u*0.01;

n = t_f/dt;


for j=(n/2) :1: n
    u(j)=-0.01;
end



for i=1 : 1 :(n-1)
    alfa_p = a*(y(i,2)-y(i,1));
    phi_p  = y(i,3);
    phi_pp = (-w^2)*(y(i,2)-y(i,1)-b*u(i));
    h_p    = c*y(i,1);

    y(i+1 , 1) = y(i,1) + dt *alfa_p;
    y(i+1 , 2) = y(i,2) + dt *phi_p;
    y(i+1 , 3) = y(i,3) + dt *phi_pp;
    y(i+1 , 4) = y(i,4) + dt *c*h_p;
end

h_0=0; %h inicial 
plot (t, y(:,4)+h_0)

