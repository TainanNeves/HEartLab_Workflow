function [x_points, y_points, Trace1] = pick_up_a_trace(I, DATA, N)
    button = 0;

    % Create figure with 2x3 subplots
    figure(99);
    clf;
    subplot(3, 4, [1 2 3 4]); % Graphs will be plotted here
    title('Data Trace');

    subplot(3, 4, [5 6 9 10]); % Image display
    imagesc(I);
    colormap('gray');
    axis square;
    title('Click on the Image');
    set(gca, 'NextPlot', 'add'); % Ensure points are added to the image

    subplot(3, 4, [7 8 11 12]); % Instructions
    text(0.1, 0.5, {'Instructions:', ...
                    '1. Left-click to plot signal from a point.', ...
                    '2. Right-click to highlight and save point.', ...
                    '3. Press Spacebar to exit.'}, ...
         'FontSize', 10);
    axis off;

    % Initialize arrays to store x and y coordinates of right-clicks
    x_points = [];
    y_points = [];

    % Set the figure's pointer to crosshair
    set(gcf, 'Pointer', 'crosshair');

    % Main loop
    while true
        % Check for spacebar press or mouse click
        was_a_key = waitforbuttonpress;

        if was_a_key
            key = get(gcf, 'CurrentKey');
            if strcmp(key, 'space')
                break;
            end
        else
            % Mouse click detected
            point = get(gca, 'CurrentPoint');
            x = round(point(1, 2));
            y = round(point(1, 1));
            button = get(gcf, 'SelectionType');

            if x <= 0 || y <= 0
                disp("Incorrect Click: Not allowed region");
                continue;
            end

            if strcmp(button, 'alt') % Right-click
                % Save the coordinates to the arrays
                x_points = [x_points; x];
                y_points = [y_points; y];

                subplot(3, 4, [5 6 9 10]); % Highlighted image display
                I(x, y) = max(max(I));
                imagesc(I);
                colormap('gray');
                axis square;
                title('Click on the Image');
            end

            % Extract and plot data trace
            Trace1 = squeeze(mean(mean(DATA(x-N:x+N, y-N:y+N, :), 1), 2));

            subplot(3, 4, [1 2 3 4]); % Plot in the first row
            plot(Trace1);
            title(['Data Trace (row: ', num2str(x), ' col: ', num2str(y), ')']);
            drawnow;

            disp(['(Max - Min) * 100 = ', num2str((max(Trace1) - min(Trace1)) * 100)]);
        end
    end
end





% function [x,y,Trace1, Trace2] = pick_up_a_trace(I,DATA,N)
% 
% button =0;
%    
% while button ~=32
%     
%    
%     figure(99),imagesc(I),colormap('gray'), axis square
%     [y,x button]=ginput(1);
%     x=round(x);y=round(y);
% 
%     if button == 3 
%         I(x,y) = max(max(I));
%         imagesc(I);
%     end
%     if x<=0 || y<=0, break, end
%     X = squeeze(mean(mean(DATA(x-N:x+N, y-N:y+N, :),1),2));
%     %X2 = squeeze(mean(mean(DATA2(x-N:x+N, y-N:y+N, :),1),2));
%     Trace1 = X;
%     %Trace2 = X2;
%     figure(103),plot(X); title(strcat('row: ', num2str(x), ' col: ', num2str(y)))
%     %figure(104),plot(X2); title(strcat('row: ', num2str(x), ' col: ', num2str(y)))
%     display((max(X)-min(X))*100);
% end