clc, clear all, close all;

%{
-------------------------------------------------------------------------
A partir de las curvas de mediciones de las variables graficadas en la Fig. 1-3, se 
requiere obtener el modelo del sistema considerando como entrada un escalón de 12V, 
como salida a la velocidad angular, y a partir de 0,1segundo se aplica un TL aproximado 
de 7,5 10-2 Nm. En el archivo Curvas_Medidas_Motor.xls están las mediciones, en la 
primer hoja los valores y en la segunda los nombres. Se requiere obtener el modelo 
dinámico, para establecer las constantes de la corriente.

-------------------------------------------------------------------------
%}

mediciones=xlsread('Curvas_Medidas_Motor.xls');

ts=0.6;
deltat=10^-7;
Va=zeros(1,ts/deltat);
for i=200000:1:6000000
    Va(1,i)=12;
end

%wr_Va

x1=78.2856518;
t1=0.02028;
x2=126.558141;
t2=0.02049;
x3=168.588639;
t3=0.02085;

%ALGORITMO DE CHEN para aproximación de FT de la forma
%G(s)=K*(T3*s+1)/((T1*s+1)*(T2*s+1))
ganancia=198.248802;
k1=(x1/ganancia)-1;
k2=(x2/ganancia)-1; 
k3=(x3/ganancia)-1; 
b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta=(2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b));
T1=-0.001/log(alfa1);
T2=-0.001/log(alfa2);
T3=beta*(T1-T2)+T1; %No hay cero, no se usa

s=tf('s');
G=ganancia/((T1*s+1)*(-T2*s+1))
[num,den]=tfdata(G,'v');
num=num/12;

t=0:deltat:(0.6-deltat);
u=zeros(1,ts/deltat);
for j=0.02/deltat:1:ts/deltat
    u(1,j)=1;
end

figure
[yg tg]=lsim(G,u,t);
plot(tg,yg,'r');title('Respuesta de la función de transferencia aproximada con entrada de tensión');ylabel('wr');xlabel('Tiempo');
hold on;
%plot(mediciones(:,1),mediciones(:,2),'r');

%%



%wr_Tl
ti=0.1;
x1i=198.25-187.834556116561;
t1i=0.10011;
x2i=198.25-153.576494352651;
t2i=0.10023;
x3i=198.25-131.424158099762;
t3i=0.10032;
 
%ALGORITMO DE CHEN para aproximación de FT de la forma
%G(s)=K*(T3*s+1)/((T1*s+1)*(T2*s+1))
gananciai=153.4;
k1i=x1i/gananciai-1;
k2i=x2i/gananciai-1; 
k3i=x3i/gananciai-1; 
bi=4*k1i^3*k3i-3*k1i^2*k2i^2-4*k2i^3+k3i^2+6*k1i*k2i*k3i;
alfa1i=(k1i*k2i+k3i-sqrt(bi))/(2*(k1i^2+k2i));
alfa2i=(k1i*k2i+k3i+sqrt(bi))/(2*(k1i^2+k2i));
betai=(2*k1i^3+3*k1i*k2i+k3i-sqrt(bi))/(sqrt(bi));
T1i=-0.00011/log(alfa1i);
T2i=-0.001/log(alfa2i);
T3i=betai*(T1i-T2i)+T1i;

F=gananciai/((T2i*s+1))
[numf,denf]=tfdata(F,'v');
numf=numf/0.075;

utl=zeros(1,ts/deltat);
for j=(1000000):1:ts/deltat
    utl(1,j)=1;
end    
figure
[ytl, ttl]=lsim(-F,utl,t);
plot(ttl,ytl,'r');title('Respuesta de la función de transferencia aproximada con entrada de torque');xlabel('Tiempo');ylabel('wr');

figure
plot(t,yg+ytl,'r');
hold on;
plot(mediciones(:,1),mediciones(:,2),'b');title('aproximación obtenida VS valores medidos');xlabel('Tiempo');ylabel('wr');legend({'Obtenida','Medida'});



