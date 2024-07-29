function [row, col, source] = getElectrodePosition(electrodeNumber)
%GETELECTRODEPOSITION Returns the row, column position, and source of a given electrode.
%
%   [row, col, source] = GETELECTRODEPOSITION(electrodeNumber) returns the row and
%   column position of the specified electrode, as well as the source number. If the
%   electrode number is invalid or the position is not defined, the function returns 
%   an empty array and displays a message indicating the electrode is not found.
%
%   Input:
%       electrodeNumber - A scalar integer representing the electrode number.
%
%   Output:
%       row - The row position of the electrode.
%       col - The column position of the electrode.
%       source - The source number of the electrode.
%
%   If the electrode number is not found or its position is undefined,
%   the function returns empty arrays for row, col, and source, and displays 
%   the message 'Electrode not found'.

    % Matrix for Electrode and its row and col coordinates
    matrix = [
        1, 9, 3; 2, 9, 5; 3, 9, 7; 4, 9, 9; 5, 7, 3; 6, 7, 5; 7, 7, 7; 8, 7, 9; 
        9, 5, 3; 10, 5, 5; 11, 5, 7; 12, 5, 9; 13, 3, 3; 14, 3, 5; 15, 3, 7; 
        16, 3, 9; 17, 3, 3; 18, 5, 3; 19, 7, 3; 20, 9, 3; 21, 3, 5; 22, 5, 5; 
        23, 7, 5; 24, 9, 5; 25, 3, 7; 26, 5, 7; 27, 7, 7; 28, 9, 7; 29, 3, 9; 
        30, 5, 9; 31, 7, 9; 32, 9, 9; 33, -1, -1; 34, -1, -1; 35, -1, -1; 36, -1, -1; 
        37, -1, -1; 38, -1, -1; 39, -1, -1; 40, -1, -1; 41, -1, -1; 42, -1, -1; 
        43, -1, -1; 44, -1, -1; 45, -1, -1; 46, -1, -1; 47, -1, -1; 48, -1, -1; 
        49, -1, -1; 50, -1, -1; 51, -1, -1; 52, -1, -1; 53, -1, -1; 54, -1, -1; 
        55, -1, -1; 56, -1, -1; 57, -1, -1; 58, -1, -1; 59, -1, -1; 60, -1, -1; 
        61, -1, -1; 62, -1, -1; 63, -1, -1; 64, -1, -1; 65, 9, 3; 66, 9, 5; 67, 9, 7; 
        68, 9, 9; 69, 7, 3; 70, 7, 5; 71, 7, 7; 72, 7, 9; 73, 5, 3; 74, 5, 5; 
        75, 5, 7; 76, 5, 9; 77, 3, 3; 78, 3, 5; 79, 3, 7; 80, 3, 9; 81, -1, -1; 
        82, -1, -1; 83, -1, -1; 84, -1, -1; 85, -1, -1; 86, -1, -1; 87, -1, -1; 
        88, -1, -1; 89, -1, -1; 90, -1, -1; 91, -1, -1; 92, -1, -1; 93, -1, -1; 
        94, -1, -1; 95, -1, -1; 96, -1, -1; 97, -1, -1; 98, -1, -1; 99, -1, -1; 
        100, -1, -1; 101, -1, -1; 102, -1, -1; 103, -1, -1; 104, -1, -1; 105, -1, -1; 
        106, -1, -1; 107, -1, -1; 108, -1, -1; 109, -1, -1; 110, -1, -1; 111, -1, -1; 
        112, -1, -1; 113, -1, -1; 114, -1, -1; 115, -1, -1; 116, -1, -1; 117, -1, -1; 
        118, -1, -1; 119, -1, -1; 120, -1, -1; 121, -1, -1; 122, -1, -1; 123, -1, -1; 
        124, -1, -1; 125, -1, -1; 126, -1, -1; 127, -1, -1; 128, -1, -1; 129, 4, 14; 
        130, 4, 16; 131, 7, 15; 132, 10, 14; 133, 10, 16; 134, 17, 14; 135, 17, 16; 
        136, 20, 15; 137, 23, 14; 138, 23, 16; 139, 4, 18; 140, 4, 20; 141, 7, 19; 
        142, 10, 18; 143, 10, 20; 144, 17, 18; 145, 4, 2; 146, 4, 4; 147, 7, 3; 
        148, 10, 2; 149, 10, 4; 150, 17, 2; 151, 17, 4; 152, 20, 3; 153, 23, 2; 
        154, 23, 4; 155, 4, 6; 156, 4, 8; 157, 7, 7; 158, 10, 6; 159, 10, 8; 
        160, 17, 6; 161, 17, 8; 162, 20, 7; 163, 23, 6; 164, 23, 8; 165, 4, 10; 
        166, 4, 12; 167, 7, 11; 168, 10, 10; 169, 10, 12; 170, 17, 10; 171, 17, 12; 
        172, 20, 11; 173, 23, 10; 174, 23, 12; 175, -1, -1; 176, -1, -1; 177, 17, 20; 
        178, 20, 19; 179, 23, 18; 180, 23, 20; 181, -1, 22; 182, -1, 24; 183, 7, 23; 
        184, 10, 22; 185, 10, 24; 186, 17, 22; 187, 17, 24; 188, 20, 23; 189, 23, 22; 
        190, 23, 24; 191, -1, -1; 192, -1, -1
    ];

    % Check if the electrode number is valid
    if electrodeNumber < 1 || electrodeNumber > size(matrix, 1)
        disp('Electrode number is not valid');
        row = [];
        col = [];
        source = [];
        return;
    end

    % Get the position
    row = matrix(electrodeNumber, 2);
    col = matrix(electrodeNumber, 3);

    % Determine the source based on the electrode number
    if electrodeNumber >= 1 && electrodeNumber <= 16
        source = 1;
    elseif electrodeNumber >= 17 && electrodeNumber <= 32
        source = 2;
    elseif electrodeNumber >= 65 && electrodeNumber <= 80
        source = 3;
    elseif (electrodeNumber >= 129 && electrodeNumber <= 174) || (electrodeNumber >= 177 && electrodeNumber <= 190)
        source = 4;
    else
        source = [];
    end

    % Check if the position is defined
    if row == -1 && col == -1
        disp('Electrode do not have signal');
        row = [];
        col = [];
        source = [];
    end
end
