%% CODE FOR LAT AND CONDUCTION VELOCITY FOR OPTICAL AND ELECTRICAL DATA
%% LAT FOR OPTICAL DATA
clear all; close all;
load("OPT_ELE.mat")
lim1=555; %samples
lim2=1066; %samples
Fs=4000; %sample frequency 
% OPTICO
O_SIN_sinle=O_SIN(:,:,lim1:lim2);
T_LAT_ms = O_SIN_sinle(:,:,1)*0;
for i=1:size(O_SIN_sinle,1)
    i
    for j=1:size(O_SIN_sinle,2)
        if max(max(squeeze(O_SIN_sinle(i,j,:))))~=0
        %i=66;j=139;
        [T_LAT_ms(i,j)]= find_LAT_linearFit_1D(squeeze(O_SIN_sinle(i,j,:)),Fs,15, 200, 0);
        end
    end
end

% plotar os lugares exactos dos eletrodos para figura OJO
T_LAT_ms_RA_S=T_LAT_ms;
d=200;color='black'; tam=11;
f1=figure('color','white','Position', [40 40 400 300]);
C = parula(256);
C(1,1:3) = [1 1 1];
J = T_LAT_ms(:);
J = imrotate(T_LAT_ms,90);
Y = prctile(nonzeros(J),[1 95],'all');
imagesc(J-Y(1), [-1, 30]);colormap(C);%Y(2)-Y(1)
hBar1=colorbar('eastoutside');ylabel(hBar1,'Local Activation Time [ms]','FontSize',14);
box off
set(gca,'fontsize', 18);
title('Cam 2');
ylabel('Pixels');xlabel('Pixels');
% Ajustar los límites de los ejes x e y para enfocar el área de interés
xlim([30, 115]);
ylim([30, 120]);
axis off;

%% CONDUCTION VELOCITY FOR OPTICAL DATA 
T_LAT_ms2 = SpatTemp_Filtering(T_LAT_ms,3,0,'CPU');
    for i=1:size(O_SIN,1)
       for j=1:size(O_SIN,2)
            if max(max(squeeze(O_SIN(i,j,:))))~=0
               % i=66;j=139;
               [AvgSpeed_O_RA(i,j),  StdSpeed_O_RA(i,j), Angle_O_RA(i,j)]= CV_CircleMethod(T_LAT_ms2 , 5, i, j, 0.33); %TEM 10 PIXESL ENTRE ELETRODOS E A DISTANCIA É 3.3 MM 
            end
        end
    end
    AvgSpeed=AvgSpeed_O_RA;
    f1=figure('color','white','Position', [40 40 400 300]);
    C = jet(256);
    C(1,1:3) = [1 1 1];
    J = AvgSpeed(:);
    Y = prctile(nonzeros(J),[5 95],'all');
    J = imrotate(AvgSpeed,90);
    teto=Y(2)-Y(1);
    imagesc(J-Y(1), [0, Y(2)-Y(1)]);colormap(C);
    hBar1=colorbar('eastoutside');ylabel(hBar1,'Conduction velocity [cm/s]','FontSize',18);
    hold on;d=200;color='black'; tam=14;
    set(gca,'fontsize', 18);
    title('Cam 2');
    ylabel('Pixels');xlabel('Pixels');
    %Ajustar los límites de los ejes x e y para enfocar el área de interés
    xlim([30, 100]);
    ylim([0, 50]);
    %axis off
    axis off;

%% DOMINANT FREQUENCY FOR OPTICAL DATA 

[DF_O_RA,Sfft_O_RA,fstep]=f_DF_optico(O_SIN,Fs);
d=200;color='black'; tam=11;
f1=figure('color','white','Position', [40 40 450 300]);
J = DF_O_RA;%atrio
Y = prctile(nonzeros(J),[0 95],'all');
J = imrotate(J,90);
imagesc(J);colormap(C);
box off;hold on;
set(gca,'fontsize', 18);
ylabel('Pixels');xlabel('Pixels');
% Ajustar los límites de los ejes x e y para enfocar el área de interés
xlim([30, 100]);
ylim([40, 100]);
axis off;
hBar1=colorbar('eastoutside');ylabel(hBar1,'Dominant Frequency [Hz]','FontSize',14);

b = length(Sfft_O_RA(62,133,:));
figure;plot((1:b)*fstep, squeeze(Sfft_O_RA(62,133,:)));hold on;
plot((1:b)*fstep, squeeze(Sfft_O_RA(48,115,:)));hold on;
plot((1:b)*fstep, squeeze(Sfft_O_RA(72,118,:)));



%% LAT FOR ELECTRICAL DATA
% lim1 and lim2= take only a little segment of the signal 
% OPT_ELE.mat= has the electrical and optical arrays

load("OPT_ELE.mat")
lim1=555; %samples
lim2=1066; %samples
Fs=4000; %sample frequency 
E_SIN_CUT=E_SIN(:,lim1:lim2);

for i=1:size(E_SIN_CUT,1)
  %i=3
  if max(max(E_SIN_CUT(i,:)))~=0
    xs=E_SIN_CUT(i,:);
    [dvdt(i), position(i)]=max(-diff(xs));
    
        figure ()
        plot(xs);
        hold on
        plot(position(i),xs(position(i)),'o','markersize',10,'color','red');
        hold off

  end
end 
position(1)=256; %error corrregido a mão para este esperimento 
LAT_ms_E=(position)*1000/Fs %changing to time in ms

% Crear una matriz de ejemplo con los tiempos de activación
minimo=min(LAT_ms_E);LAT_ms_E=LAT_ms_E-minimo;%lat(7)=243;
tiempos_activacion = [LAT_ms_E(16) LAT_ms_E(12) LAT_ms_E(8) LAT_ms_E(4);
                      LAT_ms_E(15) LAT_ms_E(11) LAT_ms_E(7) LAT_ms_E(3);
                      LAT_ms_E(14) LAT_ms_E(10) LAT_ms_E(6) LAT_ms_E(2);
                      LAT_ms_E(13) LAT_ms_E(9) LAT_ms_E(5) LAT_ms_E(1)];
LAT_16_E=tiempos_activacion;

% Crear una malla de coordenadas de los tiempos de activación
[X, Y] = meshgrid(1:size(tiempos_activacion, 2), 1:size(tiempos_activacion, 1));

% Crear una malla de coordenadas para la interpolación
[Xq, Yq] = meshgrid(1:0.1:size(tiempos_activacion, 2), 1:0.1:size(tiempos_activacion, 1));

% Realizar la interpolación bilineal
tiempos_interpolados = interp2(X, Y, tiempos_activacion, Xq, Yq, 'linear');
tiempos_interpolados_RA_S=tiempos_interpolados;

% Mostrar el mapa de tiempo de activación
tam=12;
f1=figure('color','white','Position', [40 40 200 200]);
imagesc(Xq(1,:), Yq(:,1), tiempos_interpolados,[0,30]);hBar1=colorbar('eastoutside');colormap("parula");ylabel(hBar1,'Local Activation Time [ms]','FontSize',16);%axis square;
axis off
hold on;
text(1, 1, num2str(16), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(1, 2, num2str(15), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(1, 3, num2str(14), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(1, 4, num2str(13), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(2, 1, num2str(12), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(2, 2, num2str(11), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(2, 3, num2str(10), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(2, 4, num2str(9), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(3, 1, num2str(8), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(3, 2, num2str(7), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(3, 3, num2str(6), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(3, 4, num2str(5), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(4, 1, num2str(4), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(4, 2, num2str(3), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(4, 3, num2str(2), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(4, 4, num2str(1), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
set(gca,'fontsize', 18)
title('MEA3');
ylabel('Pixels');xlabel('Pixels');

%% CONDUCTION VELOCITY FOR ELECTRICAL DATA 
% this code uses the Circle method using the Local activation times
% calculated before.
% r= circle radius

r=3;
for i=(r+1):(size(tiempos_interpolados,1)-(r+1))
        for j=(r+1):(size(tiempos_interpolados,2)-(r+1))
            if max(max(squeeze(tiempos_interpolados(i,j,:))))~=0
               [AvgSpeed_E_RA(i,j),  StdSpeed_E_RA(i,j), Angle_E_RA(i,j)]= CV_CircleMethod(tiempos_interpolados , r, i, j, 0.3);
            end
        end
end
AvgSpeed_E2 = SpatTemp_Filtering(AvgSpeed_E_RA,5,0,'CPU');
f1=figure('color','white','Position', [40 40 230 210]);
C = jet(256);
C(1,1:3) = [1 1 1];
J = AvgSpeed_E2(:);
Y = prctile(nonzeros(J),[5 95],'all');
J = imrotate(AvgSpeed_E2,90);
imagesc(J-Y(1), [0, Y(2)-Y(1)]);colormap(C);
hBar1=colorbar('eastoutside');ylabel(hBar1,'Conduction velocity [cm/s]','FontSize',18);
set(gca,'fontsize', 18);
title('Cam 2');
ylabel('Pixels');xlabel('Pixels');
axis off

%% DOMINANT FREQUENCY FOR ELECTRICAL DATA 

[DF_E_RA,Sfft_E_RA,fstep]=f_DF_electrico(E_SIN,Fs);
% Crear una matriz de ejemplo con los tiempos de activación
frec = [DF_E_RA(16)  DF_E_RA(12) DF_E_RA(8) DF_E_RA(4);
        DF_E_RA(15)  DF_E_RA(11) DF_E_RA(7) DF_E_RA(3);
        DF_E_RA(14)  DF_E_RA(10) DF_E_RA(6) DF_E_RA(2);
        DF_E_RA(13)  DF_E_RA(9)  DF_E_RA(5) DF_E_RA(1)];
frec_16_E=frec;
% Crear una malla de coordenadas de los tiempos de activación
[X, Y] = meshgrid(1:size(frec, 2), 1:size(frec, 1));
% Crear una malla de coordenadas para la interpolación
[Xq, Yq] = meshgrid(1:0.1:size(frec, 2), 1:0.1:size(frec, 1));

% Realizar la interpolación bilineal
frec_interpolados = interp2(X, Y, frec, Xq, Yq, 'linear');
frecv_interpolados_RA=frec_interpolados;

% Mostrar el mapa de Potencial
tam=12;
f1=figure('color','white','Position', [40 40 200 200]);
imagesc(Xq(1,:), Yq(:,1), frec_interpolados);hBar1=colorbar('eastoutside');colormap(C);ylabel(hBar1,'DF [Hz]','FontSize',16);%axis square;
axis off
hold on;
text(1, 1, num2str(16), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(1, 2, num2str(15), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(1, 3, num2str(14), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(1, 4, num2str(13), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(2, 1, num2str(12), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(2, 2, num2str(11), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(2, 3, num2str(10), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(2, 4, num2str(9), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(3, 1, num2str(8), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(3, 2, num2str(7), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(3, 3, num2str(6), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(3, 4, num2str(5), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(4, 1, num2str(4), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(4, 2, num2str(3), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(4, 3, num2str(2), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
text(4, 4, num2str(1), 'Color', 'black', 'FontSize', tam, 'FontWeight', 'bold');hold on;
set(gca,'fontsize', 18)
title('MEA1');
ylabel('Pixels');xlabel('Pixels');
b = length(Sfft_E_RA(1,:));
figure;plot((1:b)*fstep, Sfft_E_RA(5,:));hold on;


