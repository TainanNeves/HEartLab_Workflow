function plot_CL_E(CL_E, plotTitle, cLimits, plotMode)
    % Function: plot_CL_E
    % Description:
    %   This function plots the CL_E matrix using a color-coded representation
    %   with the jet(256) colormap and customizable options. It uses pcolor
    %   for the plot and enables interactive data cursor to display values and
    %   coordinates on click.
    %
    % Inputs:
    %   - CL_E: A matrix containing the estimated Cicle Lenght for each (x, y) coordinate.
    %   - plotTitle: Title of the plot.
    %   - cLimits: Colorbar limits as a two-element vector [cmin, cmax].
    %   - plotMode: Mode of the plot, either 'MEA' (square plot) or 'TANK' (non-square plot).

    % Define the coordinates for black circles in 'TANK' mode
    blackCircles_TANK = [
        4, 14;
        4, 16;
        7, 15;
        10, 14;
        10, 16;
        17, 14;
        17, 16;
        20, 15;
        23, 14;
        23, 16;
        4, 18;
        4, 20;
        7, 19;
        10, 18;
        10, 20;
        17, 18;
        4, 2;
        4, 4;
        7, 3;
        10, 2;
        10, 4;
        17, 2;
        17, 4;
        20, 3;
        23, 2;
        23, 4;
        4, 6;
        4, 8;
        7, 7;
        10, 6;
        10, 8;
        17, 6;
        17, 8;
        20, 7;
        23, 6;
        23, 8;
        4, 10;
        4, 12;
        7, 11;
        10, 10;
        10, 12;
        17, 10;
        17, 12;
        20, 11;
        23, 10;
        23, 12;
        17, 20;
        20, 19;
        23, 18;
        23, 20;
        4, 22;
        4, 24;
        7, 23;
        10, 22;
        10, 24;
        17, 22;
        17, 24;
        20, 23;
        23, 22;
        23, 24;
    ];

    % Define the coordinates for black circles in 'MEA' mode
    blackCircles_MEA = [
        9, 3;
        9, 5;
        9, 7;
        9, 9;
        7, 3;
        7, 5;
        7, 7;
        7, 9;
        5, 3;
        5, 5;
        5, 7;
        5, 9;
        3, 3;
        3, 5;
        3, 7;
        3, 9;
    ];

    % Create a figure with a specific size
    fig = figure('Position', [100, 100, 800, 600]); % [left, bottom, width, height]

    % Use pcolor to create the plot
    h = pcolor(CL_E);
    set(h, 'EdgeColor', 'none'); % Remove grid lines between cells

    % Set the colormap to jet(256) with smooth color transitions
    colormap(jet(256));
    shading interp; % Interpolate colors for smooth transitions

    % Add a colorbar with specified limits
    c = colorbar;
    if ~isempty(cLimits)
        caxis(cLimits);
    end
    c.Label.String = 'Cycle Length [ms]';
    c.Label.FontSize = 12;

    % Set axis labels
    xlabel('Electrode Y Coordinate');
    ylabel('Electrode X Coordinate');
    
    % Adjust the aspect ratio and axis based on the plot mode
    switch plotMode
        case 'MEA'
            axis equal tight; % Ensure equal aspect ratio and tightens the axes
            
            % Plot black circles at specified points in 'MEA' mode
            hold on;
            for i = 1:size(blackCircles_MEA, 1)
                row = blackCircles_MEA(i, 1);
                col = blackCircles_MEA(i, 2);
                plot(col, row, 'ko', 'MarkerSize', 10, 'LineWidth', 2); % 'ko' means black circle
            end
            hold off;
            
        case 'TANK'
            axis tight; % Tighten the axes without ensuring equal aspect ratio
            pbaspect([2 1 1]); % Make spacing between columns larger than spacing between rows
            
            % Plot black circles at specified points in 'TANK' mode
            hold on;
            for i = 1:size(blackCircles_TANK, 1)
                row = blackCircles_TANK(i, 1);
                col = blackCircles_TANK(i, 2);
                plot(col, row, 'ko', 'MarkerSize', 8, 'LineWidth', 2); % 'ko' means black circle
            end
            
            % Add black lines at specified x-coordinates in 'TANK' mode
            xLines = [1, 5, 9, 13, 17, 21, 25];
            for x = xLines
                plot([x x], ylim, 'k-', 'LineWidth', 1);
            end
            hold off;
            
        otherwise
            error('Unsupported plot mode "%s". Use "MEA" or "TANK".', plotMode);
    end

    % Invert the y-axis
    set(gca, 'YDir', 'reverse');

    % Set the title
    title(plotTitle);

    % Optional: Add grid lines for better visualization
    grid on;
    set(gca, 'GridColor', [1, 1, 1]); % Set grid color to white
    set(gca, 'GridAlpha', 0.5); % Set grid transparency

    % Enable data cursor mode for interactivity
    dcm = datacursormode(fig);
    set(dcm, 'UpdateFcn', @dataCursorCallback);

    % Data cursor callback function
    function output_txt = dataCursorCallback(~, event_obj)
        % Display the position and value on click
        pos = get(event_obj, 'Position');
        x_coord = pos(1);
        y_coord = pos(2);
        value = CL_E(y_coord, x_coord); % Note: MFFTi uses (y, x) convention

        output_txt = {
            ['X: ', num2str(x_coord)], ...
            ['Y: ', num2str(y_coord)], ...
            ['Cycle Length [ms]: ', num2str(value)]
            };
    end
end
