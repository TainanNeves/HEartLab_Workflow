function interpolated_matrix = interpolateZeros(matrix)
    [Zrow, Zcol] = find(matrix == 0);

    interpolated_matrix = matrix;

    for k = 1:length(Zrow)
        row = Zrow(k);
        col = Zcol(k);

        % Find neighboring known values
        neighbors = matrix(max(row-1, 1):min(row+1, size(matrix, 1)), ...
                           max(col-1, 1):min(col+1, size(matrix, 2)));
        neighbors = neighbors(neighbors ~= 0);

        % Interpolate
        if ~isempty(neighbors)
            interpolated_matrix(row, col) = mean(neighbors);
        end
    end
end