%% Interpolatin Signals
% This code create a struct containing the interpolated signals for
% MEAs and Tank. Using a laplacian interpolation methodology based
% in a 3D surface .
clear; clc;

%% Loading Electric signals
load(""); % Load the filtered electrical signals


%% Configuring
data = D_EL.Data;
Fsampling = D_EL.Header.sample_rate;


%% Check electrodes quality
lim1 = 3*4000;
lim2 = 5*4000;


% Plot MEA1 (electrodes 1:16)
figure('Name', 'MEA1 Electrodes Quality Check', 'Position', [100, 100, 1200, 800]);
for i = 1:16
    subplot(4, 4, i);
    plot(data(i, :));
    xlim([lim1, lim2]);
    title(sprintf('Electrode %d', i));
    ylabel('\muV');
    grid on;
end
sgtitle('MEA1 (Electrodes 1-16)');

% Plot MEA2 (electrodes 17:32)
figure('Name', 'MEA2 Electrodes Quality Check', 'Position', [100, 100, 1200, 800]);
for i = 17:32
    subplot(4, 4, i-16);
    plot(data(i, :));
    xlim([lim1, lim2]);
    title(sprintf('Electrode %d', i));
    ylabel('\muV');
    grid on;
end
sgtitle('MEA2 (Electrodes 17-32)');

% Plot MEA3 (electrodes 81:96)
figure('Name', 'MEA3 Electrodes Quality Check', 'Position', [100, 100, 1200, 800]);
for i = 81:96
    subplot(4, 4, i-80);
    plot(data(i, :));
    xlim([lim1, lim2]);
    title(sprintf('Electrode %d', i));
    ylabel('\muV');
    grid on;
end
sgtitle('MEA3 (Electrodes 81-96)');

% Plot TANK electrodes
figure('Name', 'TANK Electrodes Quality Check', 'Position', [100, 100, 1400, 800]);
tank_electrodes = [129:174, 177:190];
num_tank_elec = length(tank_electrodes);
num_rows = 6;
num_cols = ceil(num_tank_elec / num_rows);

for i = 1:num_tank_elec
    subplot(num_rows, num_cols, i);
    elec_num = tank_electrodes(i);
    plot(data(elec_num, :));
    xlim([lim1, lim2]);
    title(sprintf('Electrode %d', elec_num));
    ylabel('\muV');
    grid on;
end
sgtitle('TANK Electrodes');


%% Substituting values
% Define replacement map in format [target_electrode, source_electrode]
Replace_Map = [4, 8;
                82, 81;
                86, 85;
                91, 92;
                94, 93;
                130, 139;
                131, 133;
                135, 144;
                136, 138;
                145, 146;
                154, 153;
                157, 156;
                163, 164;
                165, 166;
                180, 189;
                182, 183];
data(Replace_Map(:,1), :) = data(Replace_Map(:,2), :);


%% Laplacian Interpolation
interpolated = struct();
for i = 1:4
    % interpolated_signals = electric_interp(data, i);
    interpolated_signals = electric_interp_2025(data, i);
    field_name = sprintf('case%d', i);
    interpolated.(field_name) = fill_matrix(interpolated_signals, i);
end

clear i Fsampling interpolated_signals field_name data


%% Sincrhonization
S = 8; % Multiplication factor
Fsampling = 4000; % Sampling frequency

% MEA 1
D_3D = interpolated.case1;
[numX, numY, ~] = size(D_3D);
D_SYNC.case1 = zeros(numX, numY, round(Fsampling * S)+1);
for x = 1:numX
    for y = 1:numY
        channel = squeeze(D_3D(x, y, :)); % Extract signal for electrode at position (x, y)
        [~, ido] = min(abs(D_EL.Timestamps - D_EL.TTL(1, 1)));
        ido = round(ido + (Fsampling * 1));
        idf = round(ido + (Fsampling * S));
        % Extract the signal segment from ido to idf
        signal_segment = channel(ido:idf)';
        % Store the extracted segment
        D_SYNC.case1(x, y, :) = signal_segment;
    end
end

% MEA 2
D_3D = interpolated.case2;
[numX, numY, ~] = size(D_3D);
D_SYNC.case2 = zeros(numX, numY, round(Fsampling * S)+1);
for x = 1:numX
    for y = 1:numY
        channel = squeeze(D_3D(x, y, :)); % Extract signal for electrode at position (x, y)
        [~, ido] = min(abs(D_EL.Timestamps - D_EL.TTL(1, 1)));
        ido = round(ido + (Fsampling * 1));
        idf = round(ido + (Fsampling * S));
        % Extract the signal segment from ido to idf
        signal_segment = channel(ido:idf)';
        % Store the extracted segment
        D_SYNC.case2(x, y, :) = signal_segment;
    end
end

% MEA 3
D_3D = interpolated.case3;
[numX, numY, ~] = size(D_3D);
D_SYNC.case3 = zeros(numX, numY, round(Fsampling * S)+1);
for x = 1:numX
    for y = 1:numY
        channel = squeeze(D_3D(x, y, :)); % Extract signal for electrode at position (x, y)
        [~, ido] = min(abs(D_EL.Timestamps - D_EL.TTL(1, 1)));
        ido = round(ido + (Fsampling * 1));
        idf = round(ido + (Fsampling * S));
        % Extract the signal segment from ido to idf
        signal_segment = channel(ido:idf)';
        % Store the extracted segment
        D_SYNC.case3(x, y, :) = signal_segment;
    end
end

% TANK
D_3D = interpolated.case4;
[numX, numY, ~] = size(D_3D);
D_SYNC.case4 = zeros(numX, numY, round(Fsampling * S)+1);
for x = 1:numX
    for y = 1:numY
        channel = squeeze(D_3D(x, y, :)); % Extract signal for electrode at position (x, y)
        [~, ido] = min(abs(D_EL.Timestamps - D_EL.TTL(1, 1)));
        ido = round(ido + (Fsampling * 1));
        idf = round(ido + (Fsampling * S));
        % Extract the signal segment from ido to idf
        signal_segment = channel(ido:idf)';
        % Store the extracted segment
        D_SYNC.case4(x, y, :) = signal_segment;
    end
end

% Cleaning (optional)
clear Fsampling S ans i idf ido S x y channel numX numY signal_segment D_3D


%% Electrode positions
% Function that finds electrodes position in the matrix based on its number
[row, col] = getElectrodePosition(5);
% [row, col] = getElectrodePosition_2025(5);


%% Export Struct
% Filename to save
FileName = 'E32_F02_R01_new';

% Fill Struct to export
InterpSignal.TTL = D_EL.TTL;
InterpSignal.Timestamp = D_EL.Timestamps;
InterpSignal.Opticalin = D_EL.opticalin;
InterpSignal.Header = D_EL.Header;
InterpSignal.Data.MEA1 = interpolated.case1;
InterpSignal.Data.MEA2 = interpolated.case2;
InterpSignal.Data.MEA3 = interpolated.case3;
InterpSignal.Data.TANK = interpolated.case4;
InterpSignal.Sync.MEA1 = D_SYNC.case1;
InterpSignal.Sync.MEA2 = D_SYNC.case2;
InterpSignal.Sync.MEA3 = D_SYNC.case3;
InterpSignal.Sync.TANK = D_SYNC.case4;

%Exporting
save(['InterpolatedSignals', FileName, '_filtered'], 'InterpSignal', '-v7.3');


%% 
%% Checking Interpolation Results
clear; clc;

%% Loading
load(); % Optical Filtered
load(); % Electrical Filtered
load(); % Electrical Filtered Interpolated


%% Synchronization Check Plot
% Variables
data_o1 = D_SYNC.CAM1;
data_o2 = D_SYNC.CAM2;
data_o3 = D_SYNC.CAM3;
data_e = D_SYNC.EL;
data_e1 = InterpSignal.Sync.MEA1;
data_e2 = InterpSignal.Sync.MEA2;
data_e3 = InterpSignal.Sync.MEA3;
data_e4 = InterpSignal.Sync.TANK;
% Eletrodes
electrodes = [7 22 90 142 178 172]; 

% --- Definição de Parâmetros para Checagem ---
check_time_seconds = 2.4; 
check_window_seconds = 1;

% Parâmetros Fixos
Fs = 4000;
rows = 4; % 4 Linhas
cols = 3; % 3 Colunas

% --- Preparação dos Traços ---
time_axis = (0:size(data_o1, 3)-1) / Fs;
total_time = size(data_o1, 3) / Fs;

% Selecting Optical points
Background_CAM1 = squeeze(imrotate(D_SYNC.IMG.CAM1(:,:,1), -90));
Background_CAM2 = squeeze(imrotate(D_SYNC.IMG.CAM2(:,:,1), -90));
Background_CAM3 = squeeze(imrotate(D_SYNC.IMG.CAM3(:,:,1), -90));
disp('Select 2 points for CAM1');
[x1, y1] = pick_up_a_trace(Background_CAM1, data_o1, 1);
disp('Select 2 points for CAM2');
[x2, y2] = pick_up_a_trace(Background_CAM2, data_o2, 1);
disp('Select 2 points for CAM3');
[x3, y3] = pick_up_a_trace(Background_CAM3, data_o3, 1);
% Selecting Electrical points
Background_MEA1 = squeeze(InterpSignal.Data.MEA1(:,:,1));
Background_MEA2 = squeeze(InterpSignal.Data.MEA2(:,:,1));
Background_MEA3 = squeeze(InterpSignal.Data.MEA3(:,:,1));
Background_TANK = squeeze(InterpSignal.Data.TANK(:,:,1));
disp('Select 1 point for MEA1');
[ex1, ey1] = pick_up_a_trace(Background_MEA1, data_e1, 1);
disp('Select 1 point for MEA2');
[ex2, ey2] = pick_up_a_trace(Background_MEA2, data_e2, 1);
disp('Select 1 point for MEA3');
[ex3, ey3] = pick_up_a_trace(Background_MEA3, data_e3, 1);
disp('Select 1 point for TANK');
[ex4, ey4] = pick_up_a_trace(Background_TANK, data_e4, 1);

% Gerar traços eletricos
trace_e1_map = squeeze(data_e1(ex1, ey1, 1:32001));
trace_e2_map = squeeze(data_e2(ex2, ey2, 1:32001));
trace_e3_map = squeeze(data_e3(ex3, ey3, 1:32001));
trace_e4_map = squeeze(data_e4(ex4, ey4, 1:32001));

% --- Calcular Limites da Janela ---
t_start_desired = check_time_seconds - (check_window_seconds / 2);
t_end_desired = check_time_seconds + (check_window_seconds / 2);
t_start_plot = max(0, t_start_desired);
t_end_plot = min(total_time, t_end_desired);

% --- Plotagem para Checagem de Sincronia ---
figure('color', 'white', 'Name', 'Synchronization Check');
set(gcf, 'Position', [100 100 1200 800]); 
sgtitle(['Checagem de Sincronia | t = ' num2str(check_time_seconds) ' s'], ...
            'FontSize', 14, 'FontWeight', 'bold');


% --- LINHA 1: Traços Ópticos + 1 Eletrodo Bruto ---
% CAM1 Ponto 1 (1)
subplot(rows, cols, 1);
plot(time_axis, squeeze(data_o1(x1(1), y1(1), :)));
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title('CAM1 (Optical)'); ylabel('%F');

% CAM2 Ponto 1 (2)
subplot(rows, cols, 2);
plot(time_axis, squeeze(data_o2(x2(1), y2(1), :)));
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title('CAM2 (Optical)'); ylabel('%F');

% CAM3 Ponto 1 (3)
subplot(rows, cols, 3);
plot(time_axis, squeeze(data_o3(x3(1), y3(1), :)));
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title('CAM3 (Optical)'); ylabel('%F');


% --- LINHA 2: Mapas Elétricos Interpolados (1 ponto de cada MEA) ---
% MEA1 (4)
subplot(rows, cols, 4);
plot(time_axis, trace_e1_map, 'b');
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title('MEA1 (Interpolated)'); ylabel('\muV');

% MEA2 (5)
subplot(rows, cols, 5);
plot(time_axis, trace_e2_map, 'b');
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title('MEA2 (Interpolated)'); ylabel('\muV');

% MEA3 (6)
subplot(rows, cols, 6);
plot(time_axis, trace_e3_map, 'b');
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title('MEA3 (Interpolated)'); ylabel('\muV');


% --- LINHA 3: Eletrodos Brutos (EL) ---
% EL 1 (7)
subplot(rows, cols, 7);
plot(time_axis, data_e(electrodes(1), :), 'k');
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title(['Electrode ' num2str(electrodes(1))]); ylabel('\muV');

% EL 2 (8)
subplot(rows, cols, 8);
plot(time_axis, data_e(electrodes(2), :), 'k');
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title(['Electrode ' num2str(electrodes(2))]); ylabel('\muV');

% EL 3 (9)
subplot(rows, cols, 9);
plot(time_axis, data_e(electrodes(3), :), 'k');
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title(['Electrode ' num2str(electrodes(3))]); ylabel('\muV');


% --- LINHA 4: Eletrodos Brutos (EL) restantes + TANK ---
% TANK (Mapa Interpolado) (10)
subplot(rows, cols, 10);
plot(time_axis, trace_e4_map, 'b');
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title('TANK (Interpolated)'); ylabel('\muV'); xlabel('Time (s)');

% EL 4 (11)
subplot(rows, cols, 11);
plot(time_axis, data_e(electrodes(4), :), 'k');
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title(['Electrode ' num2str(electrodes(4))]); ylabel('\muV'); xlabel('Time (s)');

% EL 5 (12)
subplot(rows, cols, 12);
plot(time_axis, data_e(electrodes(5), :), 'k');
hold on; xline(check_time_seconds, 'r', 'LineWidth', 1.5); hold off;
xlim([t_start_plot t_end_plot]); title(['Electrode ' num2str(electrodes(5))]); ylabel('\muV'); xlabel('Time (s)');

% Ajustar layout
linkaxes(findobj(gcf, 'type', 'axes'), 'x'); % Sincroniza o zoom horizontalmente



