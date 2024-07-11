function plot_CL_O(CL_O, plotTitle, cLimits)
    % Function to plot the Cycle Length (CL) map from the CL_O matrix
    %
    % Inputs:
    % - CL_O: 3D array with dimensions [x, y, 1] containing CL values
    % - plotTitle: Title for the plot
    % - cLimits: Limits for the colorbar [min, max]
    
    % Define the colormap with white for background and gray for max value
    C = jet(256);
    C(1,1:3) = [1 1 1]; % White for background
    C(256,1:3) = [0.5, 0.5, 0.5]; % Gray for the max value
    
    % Create figure
    figure('Color', 'white', 'Position', [50 50 500 500]);
    
    % Plot Cycle Length map
    imagesc(CL_O);
    colormap(C);
    box off;
    set(gca, 'FontSize', 18);
    title(plotTitle, 'FontSize', 18);
    ylabel('Pixels', 'FontSize', 14);
    xlabel('Pixels', 'FontSize', 14);
    axis off;
    
    % Add colorbar
    hBar1 = colorbar('eastoutside');
    ylabel(hBar1, 'Cycle Length [ms]', 'FontSize', 14);
    
    % Set color axis limits
    caxis(cLimits);
end
