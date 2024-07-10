%% ECGi and MEAs comparison

%% Loading data

% Load the signals
signal_file = load("C:\Users\HeartLAB\Documents\Documents\Current codes\data_filtered_sync_E14_F3_R4.mat");

% Getting only MEAs signals
meas_signal = signal_file.D_SYNC.EL(1:80,:);
% meas_signal = D_EL.Data(1:80,:);

% Load the MEA indices
files_id_meas = load("C:\Users\HeartLAB\Documents\Documents\CinC 2024\ECGi\Dados\projected_signals_exp14.mat");

% Load the estimated potentials
estimated_signal = load("C:\Users\HeartLAB\Documents\Documents\CinC 2024\ECGi\estimated_signal_E14F3R4_sync.mat");
estimated_signal = estimated_signal.(subsref(fieldnames(estimated_signal),substruct('{}',{1})));
% estimated_signal = x_hat;


%% Calculating the statisctics

% Setting the frequency sampling
Fs = 4000;

% Initialize the metrics matrix
metrics = zeros(80,6);

% Getting the correspondent MEA electrode point on the 3D geometry
for electrode = 1:80
    if electrode < 33 || electrode > 64
        if electrode < 17
            mea = meas_signal(1:16,:);
            id_mea = files_id_meas.MEAS_IDX_20000.MEA1;
            electrode_mea = electrode;
        elseif electrode < 33
            mea = meas_signal(17:32,:);
            id_mea = files_id_meas.MEAS_IDX_20000.MEA2;
            electrode_mea = electrode - 16;
        else
            mea = meas_signal(65:80,:);
            id_mea = files_id_meas.MEAS_IDX_20000.MEA3;
            electrode_mea = electrode - 64;
        end

        % Set estimation time for the measured potential
        init_time = 1; %seconds
        final_time = 2; %seconds

        init_sample = init_time * Fs + 1;
        final_sample = final_time * Fs;

        % check the estimation time to be sure the same time window will be
        % compared

        id_pot_est = id_mea(electrode_mea); % mea(electrode)
        pot_estimated = estimated_signal(id_pot_est, init_sample:final_sample); % signal(vertex, time)
        pot_measured = mea(electrode_mea, init_sample:final_sample); % mea(electrode, time)

        % Calculate metrics

        % Squared error
        sq_error = (pot_measured - pot_estimated).^2;

        % Root Mean Square Error (RMSE)
        rmse = sqrt(mean(sq_error));
        metrics(electrode, 1) = rmse;

        % Mean Squared Error (MSE)
        mse = mean(sq_error);
        metrics(electrode, 2) = mse;

        % Mean Absolute Error (MAE)
        mae = mean(abs(sq_error));
        metrics(electrode, 3) = mae;

        % Standard deviation
        std_dev = std(pot_measured);
        metrics(electrode, 4) = std_dev;

        % Normalized RMSE (NRMSE)
        nrmse = rmse / std_dev;
        metrics(electrode, 5) = nrmse;

        % Correlation
        corr_value = corr2(pot_estimated, pot_measured);
        metrics(electrode, 6) = corr_value;
    end
end

% Find the electrode with the minimum RMSE
min_rmse = min(metrics([1:11, 14:32, 65:79], 1));
electrode_min_rmse = find(metrics(:, 1) == min_rmse);

% Find the electrode with the minimum MSE
min_mse = min(metrics([1:32, 65:end], 2));
electrode_min_mse = find(metrics(:, 2) == min_mse);

% Find the electrode with the minimum MAE
min_mae = min(metrics([1:32, 65:end], 3));
electrode_min_mae = find(metrics(:, 3) == min_mae);

% Find the electrode with the minimum NRMSE
min_nrmse = min(metrics([1:32, 65:end], 5));
electrode_min_nrmse = find(metrics(:, 5) == min_nrmse);

% Calculate average RMSE
mean_rmse = mean(metrics([1:11, 14:16, 17:32, 65:79], 1)); % MEA1, MEA2 and MEA3
mean_rmse_mea1 = mean(metrics([1:11, 14:16], 1)); % MEA1
mean_rmse_mea2 = mean(metrics([17:32], 1)); % MEA2
mean_rmse_mea3 = mean(metrics([65:79], 1)); % MEA3

fprintf('Smallest RMSE: %.1f in electrode %d.\n', min_rmse, electrode_min_rmse);
fprintf('Average RMSE: %.4f\n', mean_rmse);
fprintf('Average RMSE for MEA1 (RA): %.4f\n', mean_rmse_mea1);
fprintf('Average RMSE for MEA2 (V): %.4f\n', mean_rmse_mea2);
fprintf('Average RMSE for MEA3 (LA): %.4f\n', mean_rmse_mea3);

%% Correlation analysis

threshold = 0.5;
count = 0; % Initialize the counter

% Initialize vectors to store values less than -0.5 and greater than 0.5
values_less_than = [];
values_greater_than = [];

for electrode = 1:80
if metrics(electrode, 5) <= threshold
count = count + 1; % Increment the counter if the condition is met
end
% Check and store values less than -0.5
if metrics(electrode, 5) < -0.5
values_less_than = [values_less_than; metrics(electrode, 5)];
end
% Check and store values greater than 0.5
if metrics(electrode, 5) > 0.5
values_greater_than = [values_greater_than; metrics(electrode, 5)];
end
end

fprintf('Number of values less than or equal to 0.5: %d\n', count);
fprintf('Values less than -0.5:\n');
disp(values_less_than);
fprintf('Values greater than 0.5:\n');
disp(values_greater_than);


%% Selecting only one MEA electrode

% Set the electrode numbers to the correspondent MEA
mea1 = 1:16;
mea2 = 17:32;
mea3 = 65:80;

% Choose one electrode
electrode = 10;

% Search the index of the MEA electrode
switch true
    case ismember(electrode, mea1)
        mea = meas_signal(mea1, :);
        id_mea = files_id_meas.MEAS_IDX_20000.MEA1;
        electrode_mea = electrode;
        id_est_mea = id_mea(electrode_mea); % MEA(electrode)
        
    case ismember(electrode, mea2)
        mea = meas_signal(mea2, :);
        id_mea = files_id_meas.MEAS_IDX_20000.MEA2;
        electrode_mea = electrode - 16;
        id_est_mea = id_mea(electrode_mea); % MEA(electrode)
        
    case ismember(electrode, mea3)
        mea = meas_signal(mea3, :);
        id_mea = files_id_meas.MEAS_IDX_20000.MEA3;
        electrode_mea = electrode - 64;
        id_est_mea = id_mea(electrode_mea); % MEA(electrode)
        
    otherwise
        error('Electrode not within the defined MEA ranges');
end

pot_estimated = estimated_signal(id_est_mea, :); % SIGNAL(vertex, time)

%% Comparison plot MEA x ECGi

fs = 4000;

% Define time window for plot (in seconds)
t_start = 1;
t_end = 2;

% Convert time to samples
t_start_sample = t_start * fs + 1;
t_end_sample = t_end * fs;

% Time array
time = linspace(t_start, t_end, length(t_start_sample:t_end_sample));
tline1 = 684 / 4000 + 1;

figure();
subplot(2, 1, 1);
plot(time, meas_signal(electrode, t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'black');
title(['MEA ', num2str(electrode)], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);

subplot(2, 1, 2);
plot(time, pot_estimated(t_start_sample:t_end_sample), 'LineWidth', 1, 'Color', 'black');
title(['Vertex ', num2str(v2)], 'FontSize', 15);
ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex', 'FontSize', 14);






