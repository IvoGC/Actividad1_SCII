close all; clear all ; clc;
%Caso de estudio 3. Sistema lineal de cuatro variables de estado
%{
--------------------------------------------------------------------------
1- Obtener el sistema lineal en variables de estado para el equilibrio X=[0,0,0,0]
--------------------------------------------------------------------------
%}
w=2; a=0.05; b=5; c=100;
A = [ -a,a,0,0 ; 0,0,0,0 ; w^2,-w^2,0,0 ; c,0,0,0 ];
B = [ 0; 0 ; w^2*b ; 0 ];
C = [ 1,0,0,0 ; 0,0,0,1 ];
D = 0;

sys = ss(A, B, C, D);

%%

%{
--------------------------------------------------------------------------
2- Obtener la solución numérica del sistema lineal para evaluar cuantitativamente el 
comportamiento con intención de verificar el correcto planteo. Para hacerlo, se le asignan 
los valores siguientes a los parámetros, son w=2; a=0,05; b=5; c=100 m/s, (es decir, 
360Km/h), dt=1e-3
; y el tiempo de simulación ts = 5 segundos.
--------------------------------------------------------------------------
%}

dt = 1e-3;  % step size
w=2; a=0.05; b=5; c=100; t_f=5; %tiempo total de simulacion
t = 0:dt:(t_f-dt);

%para el punto de operacion tenemos
alfa =0; phi =0; phi_p =0; h =0;

y =zeros(t_f/dt,4);
y(1,1) = alfa;  %x1
y(1,2) = phi;   %x2
y(1,3) = phi_p; %x3
y(1,4) = h;     %x4


n = t_f/dt;

%definimos una entrada que valga 0.01 durante la mitad de tiempo de
%simulacion y luego cambie su valor a -0.01
u = ones(t_f/dt,1);
u = u*0.01; 
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

h_0=4000; 
figure
plot (t, y(:,4)+h_0); title('entrada u(t), ts=5 , h_0=4000[m] y c=100[m/s](360[Km/h])'); xlabel('tiempo[s]'); ylabel('altura h[m]')

%%

%{
--------------------------------------------------------------------------
Obtener la solución numérica del sistema lineal para c=50 m/s, (es decir, 180Km/h),
dt=1e-3
; y el tiempo de simulación de 20 segundos
--------------------------------------------------------------------------
%}

dt = 1e-3;  % step size
w=2; a=0.05; b=5; c=50; t_f=20; %tiempo total de simulacion
t = 0:dt:(t_f-dt);

%para el punto de operacion tenemos
alfa =0; phi =0; phi_p =0; h =0;

y =zeros(t_f/dt,4);
y(1,1) = alfa;  %x1
y(1,2) = phi;   %x2
y(1,3) = phi_p; %x3
y(1,4) = h;     %x4
n = t_f/dt;

%definimos una entrada que valga 0.01 durante la mitad de tiempo de
%simulacion y luego cambie su valor a -0.01
u = ones(t_f/dt,1);
u = u*0.01; 
for j=1 :1: n
    if (j>=5000)
    u(j)=-0.02;
    end
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

h_0=5400; 
figure
plot (t, y(:,4)+h_0); title('entrada u(t), ts=20 , h_0=5400[m] y c=50[m/s](180[Km/h])'); xlabel('tiempo[s]'); ylabel('altura h[m]')


