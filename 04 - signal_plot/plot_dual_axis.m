function plot_dual_axis(time, optical_data, electrical_data, optical_label, electrical_label)
    yyaxis left
    plot(time, optical_data, 'b-', 'LineWidth', 1);
    ylabel(optical_label, 'Color', 'b', 'Interpreter', 'latex');
    ylim([min(optical_data) max(optical_data)]);
    
    yyaxis right
    plot(time, electrical_data, 'r-', 'LineWidth', 1);
    ylabel(electrical_label, 'Color', 'r', 'Interpreter', 'latex');
    ylim([min(electrical_data) max(electrical_data)]);
    
    xlim([0 time(end)]);
    set(gca, 'fontsize', 10);
    grid on;
end