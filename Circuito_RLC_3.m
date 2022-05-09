close all; clear all; clc;

mediciones=xlsread('Curvas_Medidas_Motor.xls');

tiempo = mediciones(:,1);
input = 5;
i=  mediciones(:,2);
vc=mediciones(:,3);


figure
plot(mediciones(:,1),mediciones(:,2)); %Corriente
figure
plot(mediciones(:,1),mediciones(:,3)); %Vc -> escalón de 12v


%SystemIdentification(); Ident
