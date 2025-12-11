function metrics = calculate_metrics(pot_measured, pot_estimated)
<<<<<<< HEAD
    % Initialize the metrics array for a single electrode
    metrics = zeros(1, 7);

    % Calculate error metrics
    error = pot_measured - pot_estimated;

    % Squared error
    sq_error = error.^2;

    % Root Mean Square Error (RMSE)
    rmse = sqrt(mean(sq_error));
    metrics(1) = rmse;

    % Mean Squared Error (MSE)
    mse = mean(sq_error);
    metrics(2) = mse;

    % Mean Absolute Error (MAE)
    mae = mean(abs(error));
    metrics(3) = mae;

    % Standard deviation of measured signal
    std_dev_measured = std(pot_measured);
    metrics(4) = std_dev_measured;

    % Standard deviation of estimated signal
    std_dev_estimated = std(pot_estimated);
    metrics(5) = std_dev_estimated;

    % Correlation coefficient
    corr_value = corr2(pot_estimated, pot_measured);
    metrics(6) = corr_value;

    % Mean relative error
    rel_error = mean(abs((pot_measured - pot_estimated) ./ pot_measured));
    metrics(7) = rel_error;
=======
%
% Author: Angélica Quadros
% Affiliation: HeartLab, UFABC
% Year: 2024
%
% This function computes a set of statistical metrics to evaluate the 
% differences between measured and estimated signals for a given electrode.
%
% Input Arguments:
% - `pot_measured`: A vector representing the measured signal from an MEA electrode.
% - `pot_estimated`: A vector representing the estimated signal for the same electrode.
%
% Output:
% - `metrics`:A row vector where each column corresponds to one of the 
%   following seven metrics:
%   1. **Root Mean Square Error (RMSE)**: Measures the average magnitude of error.
%   2. **Mean Squared Error (MSE)**: The average of squared differences between signals.
%   3. **Mean Absolute Error (MAE)**: The average of absolute differences between signals.
%   4. **Standard Deviation (Measured)**: Standard deviation of the measured signal.
%   5. **Standard Deviation (Estimated)**: Standard deviation of the estimated signal.
%   6. **Correlation Coefficient**: Measures the linear relationship between the signals.
%   7. **Mean Relative Error (MRE)**: Average relative error between measured and estimated signals.
%
% Notes:
% - Each metric is stored in its respective column in the `metrics` vector.
% - Correlation is calculated using the `corr2` function.
% - Ensure `pot_measured` and `pot_estimated` are vectors of equal length.
%
%
% Initialize the metrics array for a single electrode
metrics = zeros(1, 7);

% Calculate error metrics
error = pot_measured - pot_estimated;

% Squared error
sq_error = error.^2;

% Root Mean Square Error (RMSE)
rmse = sqrt(mean(sq_error));
metrics(1) = rmse;

% Mean Squared Error (MSE)
mse = mean(sq_error);
metrics(2) = mse;

% Mean Absolute Error (MAE)
mae = mean(abs(error));
metrics(3) = mae;

% Standard deviation of measured signal
std_dev_measured = std(pot_measured);
metrics(4) = std_dev_measured;

% Standard deviation of estimated signal
std_dev_estimated = std(pot_estimated);
metrics(5) = std_dev_estimated;

%     % Correlation coefficient
corr_value = corr2(pot_estimated, pot_measured);
metrics(6) = corr_value;


% Mean relative error
rel_error = mean(abs((pot_measured - pot_estimated) ./ pot_measured));
metrics(7) = rel_error;
>>>>>>> Development
end
