function metrics = calculate_metrics(pot_measured, pot_estimated)
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
end
