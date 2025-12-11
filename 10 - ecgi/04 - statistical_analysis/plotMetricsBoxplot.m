function plotMetricsBoxplot(metrics_data, metric_type, condition_names)
    % metrics_data: cell array of metrics matrices for each condition
    % metric_type: metric to plot (1=RMSE, 2=MSE, 3=MAE, 4=Standard Deviation of MEAs,
    % 5=Standard Deviation of ECGi, 6=Correlation, 7=Relative Error)
    % condition_names: cell array of condition names (can be rhythm,
    % regularizaton method, filtering technique...)

    if nargin < 3
        error('Please provide at least one condition name.');
    end

    % Validate input sizes
    num_conditions = numel(metrics_data);
    if num_conditions ~= numel(condition_names)
        error('Number of condition names must match the number of metrics matrices provided.');
    end

    % Validate metric_type
    if metric_type < 1 || metric_type > 7
        error('Invalid metric type specified.');
    end

    % Determine the metric name based on metric_type
    switch metric_type
        case 1
            metric_name = 'RMSE';
        case 2
            metric_name = 'MSE';
        case 3
            metric_name = 'MAE';
        case 4
            metric_name = 'Standard Deviation of MEAs';
        case 5
            metric_name = 'Standard Deviation of ECGi';
        case 6
            metric_name = 'Correlation';
        case 7
            metric_name = 'Relative Error';
        otherwise
            error('Invalid metric column specified.');
    end

    % Prepare to store all metric data and groups
    metric_all = [];
    groups = {};

    % Define colors for each condition (you can modify this list as needed)
    colors = ['r', 'g', 'b', 'k', 'm', 'c', 'y', 'w'];
        
    % Loop through each MEA region
    for mea_idx = 1:3
        for condition_idx = 1:num_conditions
            % Determine the indices for each MEA region
            switch mea_idx
                case 1
                    indices = 1:16; % MEA1
                    mea_name = 'Right Atrium';
                    mea_label = 'MEA1';
                case 2
                    indices = 17:32; % MEA2
                    mea_name = 'Ventricle';
                    mea_label = 'MEA2';
                case 3
                    indices = 65:80; % MEA3
                    mea_name = 'Left Atrium';
                    mea_label = 'MEA3';
            end

            % Extract metric data for the current MEA region and condition
            metric_data = metrics_data{condition_idx}(indices, metric_type);

            % Append metric data for this condition to overall list
            metric_all = [metric_all; metric_data];

            % Create groups for boxplot
            groups = [groups; ...
                      repmat({sprintf('%s %s %s', mea_label, mea_name, condition_names{condition_idx})}, size(metric_data))];
        end
    end

    % Create the boxplot
    figure;
    h = boxplot(metric_all, groups, 'Colors', colors(1:num_conditions)); % Use defined colors
    title(sprintf('Boxplot of %s by MEA Region and Condition', metric_name));
    ylabel(sprintf('%s Value', metric_name));
    xlabel('MEA Region and Condition');
    set(gca, 'FontSize', 12);

    % Create legend for each condition
    % Find all boxplot objects and assign legend labels
    legend_entries = cell(num_conditions, 1);
    for i = 1:num_conditions
        legend_entries{i} = sprintf('%s', condition_names{i});
    end

    % Get handles for boxplots
    box_handles = findobj(gca, 'Tag', 'Box');

    % Ensure we have the correct number of handles
    if length(box_handles) ~= num_conditions
        warning('The number of boxplot handles does not match the number of conditions.');
    end

    % Create legend
    legend(box_handles(1:num_conditions), legend_entries, 'Location', 'NorthEast');
end
