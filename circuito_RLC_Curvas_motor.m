clc, clear all, close all;

%Definicion del sistema RLC serie en variables de estado
%Variables de estado: x1=i, x2=vc

%{
-------------------------------------------------------------------------
    En el archivo Curvas_Medidas_RLC.xls (datos en la hoja 1 y etiquetas en la hoja 2)
encontrarán las series de datos que deberían emplear para deducir los valores de R, L y 
C del circuito. Emplear el método de la respuesta al escalón, tomando como salida la 
tensión en el capacitor.               
-------------------------------------------------------------------------
%}

%Importo gráficos de mediciones y grafico
mediciones=xlsread('Curvas_Medidas_RLC.xls');
%Columnas: Tiempo Corriente Vcap

%Definicion de la entrada
u=zeros(1,1000);
paso=0.1/1000;
t=0:paso:(0.1-paso);

signo=true;
for i=100:1:1000
    if mod(i,500)==0
       signo=not(signo);
    end
    if signo==1
        u(1,i)=12;
    end
    if signo==0
        u(1,i)=-12;
    end
end
figure
plot(t,u)

%Obtención de R L y C

%ALGORITMO DE CHEN para aproximación de FT de la forma
%G(s)=K*(T3*s+1)/((T1*s+1)*(T2*s+1))

k1=(3.46733881642214/12)-1; %t1=0.012 y(t1)=3.467338816 TENER EN CUENTA EL DELAY!
k2=(6.60704515416013/12)-1; %2t1=0.014 y(2t1)=6.607045154
k3=(8.60264182057171/12)-1; %3t1=0.016 y(3t1)=8.602641820
b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta=(2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b));
T1=-0.002/log(alfa1);
T2=-0.002/log(alfa2);
T3=beta*(T1-T2)+T1; %No hay cero, no se usa

s=tf('s');
G=12/((T1*s+1)*(T2*s+1))

[num,den]=tfdata(G,'v');

%Se sabe que Vc/vin (s) = (1/(L*C)) / (s^2+R*s/L+1/(L*C))
den_norm=den/(den(1))
num_norm=num/(12*den(1))
L=0.1;
R=L*den_norm(2);
Cap=1/(L*den_norm(3));

figure
[yaprox,taprox]=lsim(G,u/12,t);
plot(0.012,3.467338816,'x'); %defino T1=0.002 y le sumo el delay
hold on;
plot(0.014,6.607045154,'x');
plot(0.016,8.602641820,'x');
plot(mediciones(:,1),mediciones(:,3));
%plot(taprox,yaprox);
title('Respuesta del voltaje del capacitor aproximada');
xlabel('Tiempo');
legend({'vc(t1)','vc(2t1)','vc(3t1','vc(t)'},'Location','southeast')

%verifico que el sistema está bien aproximado
figure
plot(mediciones(:,1),mediciones(:,3),'G');
hold on;
plot(taprox,yaprox,'R'); %Son muy parecidos, la aproximación es válida 
legend({'vc(t) real','vc(t) aproxoimada'},'Location','southeast')
title('Comparacion de curvas');
xlabel('Tiempo');

%{
-------------------------------------------------------------------------
 Una vez determinados los parámetros R, L y C, emplear la serie de corriente desde 
0.05seg en adelante para validar el resultado.     
-------------------------------------------------------------------------
%}

%Matrices
A=[-R/L -1/L; 1/Cap 0];
B=[1/L; 0];
C=[1  0];
D=0;

%Definicion de la ecuación de estado y de salida (salida de corriente)
sys1=ss(A,B,C,D)
figure
[yout,yt]=lsim(sys1,u,t);

plot(yt,yout,'R'); %% 
hold on;
plot(mediciones(:,1),mediciones(:,2),'G'); %se verifica con los componentes obtenidos
legend({'i(t) aproximada','i(t) real'},'Location','southeast');
title('Contraste de respuestas');
xlabel('Tiempo');