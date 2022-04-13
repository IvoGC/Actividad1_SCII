%https://github.com/IvoGC/Actividad1_SCII/edit/main/pendulo_invertido.m#L262
close all; clear  ; clc;
%Caso de estudio 4. Sistema no lineal de cuatro variables de estado
%{
--------------------------------------------------------------------------
1- Obtener simulaciones del sistema en las condiciones iniciales
Xo=[0, 0, -0.01, 0] y Xo=[0, 0, 3.01, 0]
con dt=10^-4, tiempo de simulacion ts=10 y u=0
--------------------------------------------------------------------------
%}

%parametros para la primera consigna puntos de op X1=[0,0,-0.01,0] y X2=[0,0,3.01,0]
m = 0.1; F = 0.1; l = 0.6; g = 9.8; M = 0.5;
dt = 0.1e-3;
tf = 10; %tiempo final de simulacion
u=0;

%defino condifiones iniciales  a partir del punto de operacion
delta = 0; delta_p = 0; phi = -0.01; phi_p = 0;

t=0:dt:tf-dt;
y = zeros(4,tf/dt) ;
y(1,1) = delta;   % = x1
y(2,1) = delta_p; % = x2
y(3,1) = phi;     % = x3
y(4,1) = phi_p;   % = x4

n = tf/dt;

%Bucle para el sistema NO linealizado
for i=1 :1:n-1
    delta_pp = (-m*g*cos(y(3,i))*sin(y(3,i))+(m*l*(y(4,i))^2*sin(y(3,i)))-(F*y(2,i))+u)/(M+m-m*cos(y(3,i))*cos(y(3,i)));
    phi_pp   = ((g*sin(y(3,i)))-(delta_pp*cos(y(3,i))))/l;
    
    y(1,i+1) = y(1,i) + dt* y(2,i);
    y(2,i+1) = y(2,i) + dt* delta_pp;
    y(3,i+1) = y(3,i) + dt* y(4,i);
    y(4,i+1) = y(4,i) + dt* phi_pp ;
    
end
figure
plot (t , y(1,:)) ; title('distancia delta con phi_in=-0.01 y m=0.1 sistema_NL'); xlabel('tiempo[s]'); ylabel('delta[m]');
figure
plot (t, y(3,:)); title('angulo phi_in=-0.01 y m=0.1 sistema_NL'); xlabel('tiempo[s]'); ylabel('phi[rad]');

%defino condifiones iniciales  a partir del punto de operacion
delta = 0; delta_p = 0; phi = 3.01; phi_p = 0;

t=0:dt:tf-dt;
y = zeros(4,tf/dt) ;
y(1,1) = delta;   % = x1
y(2,1) = delta_p; % = x2
y(3,1) = phi;     % = x3
y(4,1) = phi_p;   % = x4

n = tf/dt;


%Bucle para el sistema NO linealizado
for i=1 :1:n-1
    delta_pp = (-m*g*cos(y(3,i))*sin(y(3,i))+(m*l*(y(4,i))^2*sin(y(3,i)))-(F*y(2,i))+u)/(M+m-m*cos(y(3,i))*cos(y(3,i)));
    phi_pp   = ((g*sin(y(3,i)))-(delta_pp*cos(y(3,i))))/l;
    
    y(1,i+1) = y(1,i) + dt* y(2,i);
    y(2,i+1) = y(2,i) + dt* delta_pp;
    y(3,i+1) = y(3,i) + dt* y(4,i);
    y(4,i+1) = y(4,i) + dt* phi_pp ;
    
end
figure
plot (t , y(1,:)) ; title('distancia delta con phi_in=3.01 y m=0.1 sistema_NL'); xlabel('tiempo[s]'); ylabel('delta[m]');
figure
plot (t, y(3,:)); title('angulo phi_in=3.01 y m=0.1 sistema_NL'); xlabel('tiempo[s]'); ylabel('phi[rad]');




%%
%{
--------------------------------------------------------------------------
2- Modificar la masa m al doble y repetir la operaci�n anterior
--------------------------------------------------------------------------
%}
%parametros para la primera consigna puntos de op X1=[0,0,-0.01,0] y X2=[0,0,3.01,0]
m = 0.2; F = 0.1; l = 0.6; g = 9.8; M = 0.5;
u=0; dt = 0.1e-3; tf = 10; %tiempo final de simulacion


%defino condifiones iniciales  a partir del punto de operacion
delta = 0; delta_p = 0; phi = -0.01; phi_p = 0;

t=0:dt:tf-dt;
y = zeros(4,tf/dt) ;
y(1,1) = delta;   % = x1
y(2,1) = delta_p; % = x2
y(3,1) = phi;     % = x3
y(4,1) = phi_p;   % = x4

n = tf/dt;

%Bucle para el sistema NO linealizado
for i=1 :1:n-1
    delta_pp = (-m*g*cos(y(3,i))*sin(y(3,i))+(m*l*(y(4,i))^2*sin(y(3,i)))-(F*y(2,i))+u)/(M+m-m*cos(y(3,i))*cos(y(3,i)));
    phi_pp   = ((g*sin(y(3,i)))-(delta_pp*cos(y(3,i))))/l;
    
    y(1,i+1) = y(1,i) + dt* y(2,i);
    y(2,i+1) = y(2,i) + dt* delta_pp;
    y(3,i+1) = y(3,i) + dt* y(4,i);
    y(4,i+1) = y(4,i) + dt* phi_pp ;
    
end
figure
plot (t , y(1,:)) ; title('distancia delta con phi_in=-0.01 y m=0.2 sistema_NL'); xlabel('tiempo[s]'); ylabel('delta[m]');
figure
plot (t, y(3,:)); title('angulo phi_in=-0.01 y m=0.2 sistema_NL'); xlabel('tiempo[s]'); ylabel('phi[rad]');

%defino condifiones iniciales  a partir del punto de operacion
delta = 0; delta_p = 0; phi = 3.01; phi_p = 0;

t=0:dt:tf-dt;
y = zeros(4,tf/dt) ;
y(1,1) = delta;   % = x1
y(2,1) = delta_p; % = x2
y(3,1) = phi;     % = x3
y(4,1) = phi_p;   % = x4

n = tf/dt;
%Bucle para el sistema NO linealizado
for i=1 :1:n-1
    delta_pp = (-m*g*cos(y(3,i))*sin(y(3,i))+(m*l*(y(4,i))^2*sin(y(3,i)))-(F*y(2,i))+u)/(M+m-m*cos(y(3,i))*cos(y(3,i)));
    phi_pp   = ((g*sin(y(3,i)))-(delta_pp*cos(y(3,i))))/l;
    
    y(1,i+1) = y(1,i) + dt* y(2,i);
    y(2,i+1) = y(2,i) + dt* delta_pp;
    y(3,i+1) = y(3,i) + dt* y(4,i);
    y(4,i+1) = y(4,i) + dt* phi_pp ;
    
end
figure
plot (t , y(1,:)) ; title('distancia delta con phi_in=3.01 y m=0.2 sistema_NL'); xlabel('tiempo[s]'); ylabel('delta[m]');
figure
plot (t, y(3,:)); title('angulo phi_in=3.01 y m=0.2 sistema_NL'); xlabel('tiempo[s]'); ylabel('phi[rad]');


%%

%{
--------------------------------------------------------------------------
3- Obtener la representaci�n lineal en variables de estado para el
equilibrio inestable X=[0,0,0,0]
--------------------------------------------------------------------------
%}
m = 0.1; F = 0.1; l = 0.6; g = 9.8; M = 0.5;
u=0; tf=5; dt=0.1e-3;
t= 0:dt:tf;

A = [0,1,0,0 ; 0,(-F/M),(-m/M),0 ; 0,0,0,1 ; 0,(F/(l*M)),(((M+m)*g)/(l*M)),0 ];
B = [0 ; (1/M) ; 0 ; (-1/(l*M))];
C = [1,0,0,0 ; 0,0,1,0];
D = 0;

sys=ss(A,B,C,D);
%lsim(sys,u,t);

%%







%%
close all; clear  ; clc;
%{
--------------------------------------------------------------------------
4- Obtener la soluci�n num�rica de los dos sistemas, del lineal y del no lineal para evaluar
cuantitativamente la equivalencia, modificando m de 0,1 a 0,01 y la longitud l a 1,2m.
--------------------------------------------------------------------------
%}
m = 0.01; F = 0.1; l = 1.2; g = 9.8; M = 0.5;
u=0; tf=2; dt=0.1e-3;

%defino condifiones iniciales  a partir del punto de operacion
delta = 0; delta_p = 0; phi = -0.01; phi_p = 0;

t=0:dt:tf-dt;
%Sistema NO Linealizado
yNL = zeros(4,tf/dt) ;
yNL(1,1) = delta;   % = x1
yNL(2,1) = delta_p; % = x2
yNL(3,1) = phi;     % = x3
yNL(4,1) = phi_p;   % = x4


n = tf/dt;

%Bucle para el sistema NO Linealizado
for i=1 :1:n-1
    delta_pp = (-m*g*cos(yNL(3,i))*sin(yNL(3,i))+(m*l*(yNL(4,i))^2*sin(yNL(3,i)))-(F*yNL(2,i))+u)/(M+m-m*cos(yNL(3,i))*cos(yNL(3,i)));
    phi_pp   = ((g*sin(yNL(3,i)))-(delta_pp*cos(yNL(3,i))))/l;
    
    yNL(1,i+1) = yNL(1,i) + dt* yNL(2,i);
    yNL(2,i+1) = yNL(2,i) + dt* delta_pp;
    yNL(3,i+1) = yNL(3,i) + dt* yNL(4,i);
    yNL(4,i+1) = yNL(4,i) + dt* phi_pp ;
    
end

%Sistema Linealizado
y = zeros(4,tf/dt) ;
y(1,1) = delta;   % = x1
y(2,1) = delta_p; % = x2
y(3,1) = phi;     % = x3
y(4,1) = phi_p;   % = x4

% % %Bucle para el sistema Linealizado en X=[0,0,0,0]
for i=1 :1:n-1
    y(1,i+1) = y(1,i) + dt* y(2,i);
    
    y(2,i+1) = y(2,i) + dt* (-((m*g*y(3,i))/M)+(((m*l*(y(4,i))^2*y(3,i))/M))-(((F*y(2,i))/M))+(u/M));
    
    y(3,i+1) = y(3,i) + dt* y(4,i);
    
    y(4,i+1) = y(4,i) + dt* ((y(3,i)*g*(m+M))/(M*l)-(((m*y(4,i))^2*y(3,i))/M)+((F*y(2,i))/(M*l))-(u/(M*l))) ;
    
end

figure
plot (t , yNL(1,:),'G') ;
title('distancia delta con phi=-0.01, l=1.2 y m=0.01 sistemaNL');
xlabel('tiempo[s]');
ylabel('delta[m]');
hold on
plot (t , y(1,:),'R') ; title('distancia delta con phi=-0.01, l=1.2 y m=0.01 sistemaL VS sistemaNL'); xlabel('tiempo[s]'); ylabel('delta[m]');

figure
plot (t, yNL(3,:),'G'); title('angulo phi=-0.01, l=1.2 y m=0.01 sistemaNL'); xlabel('tiempo[s]'); ylabel('phi[rad]');
hold on
plot (t, y(3,:),'R'); title('angulo phi=-0.01, l=1.2 y m=0.01 sistemaL VS sistemaNL'); xlabel('tiempo[s]'); ylabel('phi[rad]');

%%
%{
--------------------------------------------------------------------------
5- Obtener el sistema lineal para el equilibrio estable X=[0,0,pi,0]
--------------------------------------------------------------------------
%}
m = 0.01; F = 0.1; l = 12; g = 9.8; M = 0.5;
u=0; tf=5; dt=0.1e-3;
t= 0:dt:tf;

A = [0,1,0,0 ; 0,(-F/M),(-m*g/M),0 ; 0,0,0,1 ; 0,(-F/(l*M)),(-((M+m)*g)/(l*M)),0 ];
B = [0 ; (1/M) ; 0 ; (-1/(l*M))];
C = [1,0,0,0 ; 0,0,1,0];
D = 0;

sys=ss(A,B,C,D);
%lsim(sys,u,t);



%%

%{
--------------------------------------------------------------------------
6- Obtener la soluci�n num�rica de los dos sistemas, del lineal y del no lineal para evaluar
cuantitativamente la equivalencia, modificando m de 0,1 a 0,5 y la longitud l a 12m.
--------------------------------------------------------------------------
%}
m = 0.5; F = 0.1; l = 12; g = 9.8; M = 0.5;
u=0; tf=1; dt=0.1e-3;

%defino condifiones iniciales  a partir del punto de operacion
delta = 0; delta_p = 0; phi = 0.01; phi_p = 0;
tf=5
t=0:dt:tf-dt;
%Sistema NO Linealizado
yNL = zeros(4,tf/dt) ;
yNL(1,1) = delta;   % = x1
yNL(2,1) = delta_p; % = x2
yNL(3,1) = phi+pi;     % = x3
yNL(4,1) = phi_p;   % = x4

n = tf/dt;

%Bucle para el sistema NO Linealizado
for i=1 :1:n-1
    delta_pp = (-m*g*cos(yNL(3,i))*sin(yNL(3,i))+(m*l*(yNL(4,i))^2*sin(yNL(3,i)))-(F*yNL(2,i))+u)/(M+m-m*cos(yNL(3,i))*cos(yNL(3,i)));
    phi_pp   = ((g*sin(yNL(3,i)))-(delta_pp*cos(yNL(3,i))))/l;
    
    yNL(1,i+1) = yNL(1,i) + dt* yNL(2,i);
    yNL(2,i+1) = yNL(2,i) + dt* delta_pp;
    yNL(3,i+1) = yNL(3,i) + dt* yNL(4,i);
    yNL(4,i+1) = yNL(4,i) + dt* phi_pp ;
    
end

%Sistema Linealizado
y = zeros(4,tf/dt) ;
y(1,1) = delta;   % = x1
y(2,1) = delta_p; % = x2
y(3,1) = phi+pi;     % = x3
y(4,1) = phi_p;   % = x4

% % %Bucle para el sistema Linealizado en X=[0,0,pi,0]
for i=1 :1:n-1
    delta_pp =  (u-F*y(2,i)-m*g*(y(3,i)-pi))/M;
    phi_pp = (-F*y(2,i)-g*(M+m)*(y(3,i)-pi)-u)/(l*M) ;    
    y(1,i+1) = y(1,i) + dt* y(2,i);
    y(2,i+1) = y(2,i) + dt* delta_pp;
    y(3,i+1) = y(3,i) + dt* y(4,i);
    y(4,i+1) = y(4,i) + dt* phi_pp ;
end

figure
plot (t , yNL(1,:),'G') ; title('distancia delta con phi=3.01, l=12 y m=0.5 sistemaNL'); xlabel('tiempo[s]'); ylabel('delta[m]');
hold on
plot (t , y(1,:),'R') ; title('distancia delta con phi=3.01, l=12 y m=0.5 sistemaL VS sistemaNL'); xlabel('tiempo[s]'); ylabel('delta[m]');

figure
plot (t, yNL(3,:),'G'); title('angulo phi=3.01, l=12 y m=0.5 sistemaNL'); xlabel('tiempo[s]'); ylabel('phi[rad]');
hold on
plot (t, y(3,:),'+R'); title('angulo phi=3.01, l=12 y m=0.5 sistemaL VS sistemaNL'); xlabel('tiempo[s]'); ylabel('phi[rad]');







%%

% % % %Bucle para el sistema linealizado en X=[0,0,0,0]
% % for i=1 :1:n-1
% %   y(1,i+1) = y(1,i) + dt* y(2,i);
% %
% %   y(2,i+1) = y(2,i) + dt* (-((m*g*y(3,i))/M)+(((m*l*(y(4,i))^2*y(3,i))/M))-(((F*y(2,i))/M))+(u/M));
% %
% %   y(3,i+1) = y(3,i) + dt* y(4,i);
% %
% %   y(4,i+1) = y(4,i) + dt* ((y(3,i)*g*(m+M))/(M*l)-(((m*y(4,i))^2*y(3,i))/M)+((F*y(2,i))/(M*l))-(u/(M*l))) ;
% %
% % end

%ESTO NO ES PARTE DEL CODIGO
% % % %Bucle para el sistema linealizado
% for i=1 :1:n-1
%   y(1,i+1) = y(1,i) + dt* y(2,i);
%
%   y(2,i+1) = y(2,i) + dt* (-((m*g*y(3,i))/M)+(((m*l*(y(4,i))^2*y(3,i))/M))-(((F*y(2,i))/M))+(u/M));
%
%   y(3,i+1) = y(3,i) + dt* y(4,i);
%
%   y(4,i+1) = y(4,i) + dt* ((y(3,i)*g*(m+M))/(M*l)-(((m*y(4,i))^2*y(3,i))/M)+((F*y(2,i))/(M*l))-(u/(M*l))) ;
%
% end
%
%
%
% %A = [0,1,0,0 ; 0,(-F/M),(-m/M),0 ; 0,0,0,1 ; 0,(F/(l*M)),(((M+m)*g)/(l*M))]
% %B = [0 ; (1/M) ; (-1/(l*M)) ; 0]
% %C = [1,0,0,0 ; 0,0,1,0]
% %D = 0
%
