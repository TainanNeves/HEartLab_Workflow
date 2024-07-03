function matrix = fill_matrix(interpolated_signals, i)
% FILL_MATRIX Creates a matrix filled with interpolated signals.
% 
%   matrix = FILL_MATRIX(interpolated_signals, i) takes the interpolated
%   signals and the index i, and returns a matrix filled with the 
%   interpolated signals.
%
%   Input:
%       interpolated_signals - A matrix containing the interpolated signals.
%                              The size of this matrix should be at least
%                              121 for MEA data and 625 for TANK data.
%       i                    - An integer index. If i < 4, the function
%                              assumes MEA data and creates an 11x11 matrix.
%                              If i >= 4, the function assumes TANK data and
%                              creates a 25x25 matrix.
%
%   Output:
%       matrix               - A matrix filled with the interpolated signals.
%                              If the size of interpolated_signals is too short,
%                              the function returns an empty matrix and 
%                              displays an error message.
%
%   Example:
%       interpolated_signals = rand(121, 10); % Example interpolated signals
%       i = 2; % Example index for MEA data
%       matrix = fill_matrix(interpolated_signals, i);
%
%   See also: OTHER_FUNCTIONS

% Check if the input index is less than 4
if i < 4
    % Check if the size of interpolated signals is sufficient for MEA data
    if size(interpolated_signals, 1) < 121
        disp('Error: Interpolated signal is too short');
        matrix = [];
        return
    end
    % Create an 11x11 working matrix
    matrix = zeros(11, 11, size(interpolated_signals, 2));
    % Fill the matrix line by line from left to right
    for row = 1:11
        for col = 1:11
            matrix(row, col, :) = interpolated_signals((row - 1) * 11 + col, :);
        end
    end
else
    % Check if the size of interpolated signals is sufficient for TANK data
    if size(interpolated_signals, 1) < 625
        disp('Error: Interpolated signal is too short');
        matrix = [];
        return
    end
    % Create a 25x25 working matrix
    matrix = zeros(25, 25, size(interpolated_signals, 2));
    % Fill the matrix line by line from left to right
    for row = 1:25
        for col = 1:25
            matrix(row, col, :) = interpolated_signals((row - 1) * 25 + col, :);
        end
    end
end

end
