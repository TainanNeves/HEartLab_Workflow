<<<<<<< HEAD:Plot_Mapa_Tanque_3D/estimacao/plot_mapa_tanque.m
clear;
% Sinal
file =  load('sinal60_exp8.mat');
signal = file.sinal60_exp8;

% Eletrodos
dados_ele = load('C:\Users\HEartLab\Documents\Heartlab codes\Cardiac-Electrophysiology\Plot_Mapa_Tanque_3D\index_eletrodos_60.mat');
eletrodos = dados_ele.(subsref(fieldnames(dados_ele),substruct('{}',{1})));

% Tanque
dados_tank = load('C:\Users\HEartLab\Documents\Heartlab codes\Cardiac-Electrophysiology\Plot_Mapa_Tanque_3D\geometrie_tank_simples.mat');
geo_tank = dados_tank.(subsref(fieldnames(dados_tank),substruct('{}',{1})));

% definir instantes de tempo para plot
t1 = 100;
t2 = 62;
t3 = 78;

%conversão para mv
signal_mv = signal*0.00019;

%interpolação
y = signal_mv; 

indices = eletrodos';

[lap,edge] = mesh_laplacian(geo_tank.vertices, geo_tank.faces); 
[int, keepindex, repindex] = mesh_laplacian_interp_Gabriel(lap,indices);

% manter os dados medidos
ECG_interp = int *  y;  
for i=1:60
ECG_interp(indices(1,i),:)= y(i,:);  
end

%% Plotagem do resultado

% organizando geometria
faces = geo_tank.faces;
x = geo_tank.vertices(:,1);
y = geo_tank.vertices(:,2);
z = geo_tank.vertices(:,3);

% organizando sinal
sinal1 = ECG_interp(:,t1); % instante 1
sinal2 = ECG_interp(:,t2); % instante 2
sinal3 = ECG_interp(:,t3); % instante 3

% Plotagem
tlo = tiledlayout(1,3); % grade com 1 linha e 3 colunas
% Top plot
f(1) = nexttile(tlo) ;
trisurf(faces,x,y,z,sinal1,'facecolor','interp');
% med plot
f(2) = nexttile(tlo);
trisurf(faces,x,y,z,sinal2,'facecolor','interp');
% Bottom plot
f(3) = nexttile(tlo);
trisurf(faces,x,y,z,sinal3,'facecolor','interp');

% % propriedade da barra de cor
% set(f, 'Colormap', jet, 'CLim', [-8 1]); % mesma barra e limites para todos os plots
% c = colorbar(f(end)); % assign color bar to one tile 
% c.Layout.Tile = 'east'; % posição da barra
=======
clear;
% Sinal
file =  load('sinal60_exp8.mat');
signal = file.filtered;

% Eletrodos
dados_ele = load('E:\HEartLab\HeartLab Codes\Plot_Mapa_Tanque\index_eletrodos_60.mat');
eletrodos = dados_ele.(subsref(fieldnames(dados_ele),substruct('{}',{1})));

% Tanque
dados_tank = load('E:\HEartLab\HeartLab Codes\Plot_Mapa_Tanque\geometrie_tank_simples.mat');
geo_tank = dados_tank.(subsref(fieldnames(dados_tank),substruct('{}',{1})));

% definir instantes de tempo para plot
t1 = 58;
t2 = 62;
t3 = 78;

%conversão para mv
signal_mv = signal*0.00019;

%interpolação
y = signal_mv; 

indices = eletrodos';

[lap,edge] = mesh_laplacian(geo_tank.vertices, geo_tank.faces); 
[int, keepindex, repindex] = mesh_laplacian_interp_Gabriel(lap,indices);

% manter os dados medidos
ECG_interp = int *  y;  
for i=1:60
ECG_interp(indices(1,i),:)= y(i,:);  
end

%% Plotagem do resultado

% organizando geometria
faces = geo_tank.faces;
x = geo_tank.vertices(:,1);
y = geo_tank.vertices(:,2);
z = geo_tank.vertices(:,3);

% organizando sinal
sinal1 = ECG_interp(:,t1); % instante 1
sinal2 = ECG_interp(:,t2); % instante 2
sinal3 = ECG_interp(:,t3); % instante 3

% Plotagem
tlo = tiledlayout(1,3); % grade com 1 linha e 3 colunas
% Top plot
f(1) = nexttile(tlo) ;
trisurf(faces,x,y,z,sinal1,'facecolor','interp');
% med plot
f(2) = nexttile(tlo);
trisurf(faces,x,y,z,sinal2,'facecolor','interp');
% Bottom plot
f(3) = nexttile(tlo);
trisurf(faces,x,y,z,sinal3,'facecolor','interp');

% propriedade da barra de cor
set(f, 'Colormap', jet, 'CLim', [-8 1]); % mesma barra e limites para todos os plots
c = colorbar(f(end)); % assign color bar to one tile 
c.Layout.Tile = 'east'; % posição da barra
>>>>>>> 3370bdbe07e200cc05b32d2c865a1989f2d1863c:3D_tank_plot_map/estimacao/plot_mapa_tanque.m
