clc, clear all, close all;

%Motor con carga
%Variables de estado: x1=ia, x2=wr, x3=titat

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas
-------------------------------------------------------------------------
%}

Laa=366e-6;
J=5e-9;
Ra=55.6;
Bm=0;
Ki=6.49e-3;
Km=6.53e-3;

%Simulación por integración de Euler
deltat=10^-6;
ts=5;

Va=zeros(1,ts/deltat);
for i=0.5/deltat:1:ts/deltat
    Va(i)=12;
end

cant_torques=6;
Tl=zeros(ts/deltat,cant_torques);
for j=1:1:cant_torques
    for i=2/deltat:1:ts/deltat
        Tl(i,j)=0+j*0.00035;
    end
end %Genero una matriz de distintos torques

t=0:deltat:(ts-deltat);
variables=zeros(3,ts/deltat,cant_torques);
ia_0=0; wr_0=0; titat_0=0; %Cond. iniciales

for j=1:1:cant_torques
    for i=2:1:(ts/deltat-1)
        variables(1,i,j)=variables(1,i-1,j)+deltat*(-Ra*variables(1,i-1,j)/Laa-Km*variables(2,i-1,j)/Laa+Va(1,i-1)/Laa);
        variables(2,i,j)=variables(2,i-1,j)+deltat*(Ki*variables(1,i-1,j)/J-Bm*variables(2,i-1,j)/J-Tl(i-1,j)/J);
        variables(3,i,j)=variables(3,i-1,j)+deltat*variables(2,i-1,j);
    end
end

figure
plot(t,variables(2,:,1));
hold on;
xlabel('Tiempo');
ylabel('wr');
title('Velocidad angular con una entrada de 12 v con delay ante distintos torques');

for i=2:1:cant_torques
    plot(t,variables(2,:,i));
end


figure
plot(t,variables(1,:,1));
hold on;
xlabel('Tiempo');
ylabel('ia');
title('corriente con una entrada de 12 v con delay ante distitntos torques');
for i=2:1:cant_torques
    plot(t,variables(1,:,i));
end  
