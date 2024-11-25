function distance = geodesicDistance(el1, el2, file_id_meas, resolution, scaleX, scaleY, scaleZ)
    % Extract vertices and faces

    if ~isnumeric(resolution)
        vertices = file_id_meas.geometry_HR.vertices;
        faces = file_id_meas.geometry_HR.faces;
  
    else
        vertices = file_id_meas.(['geometry_' num2str(resolution)]).vertices;
        faces = file_id_meas.(['geometry_' num2str(resolution)]).faces;
    end
    
%% Normalize the mesh coordinates to the range [0, 1]
% do it if the geometry is not in cm or if the unity is unknown

%     x_min = min(vertices(:, 1)); x_max = max(vertices(:, 1)); 
%     y_min = min(vertices(:, 2)); y_max = max(vertices(:, 2)); 
%     z_min = min(vertices(:, 3)); z_max = max(vertices(:, 3)); 
%     
%     % Normalize vertices
%     vertices(:, 1) = (vertices(:, 1) - x_min) / (x_max - x_min); % Normalize x
%     vertices(:, 2) = (vertices(:, 2) - y_min) / (y_max - y_min); % Normalize y
%     vertices(:, 3) = (vertices(:, 3) - z_min) / (z_max - z_min); % Normalize z
%     
%     % Scale the vertices to desired physical dimensions
%     vertices(:, 1) = vertices(:, 1) * scaleX; % Scale x to desired dimension
%     vertices(:, 2) = vertices(:, 2) * scaleY; % Scale y to desired dimension
%     vertices(:, 3) = vertices(:, 3) * scaleZ; % Scale z to desired dimension
%% obtain the corresponding vertices for the given electrodes

    [vertex1] = get_electrode_position(el1, file_id_meas, resolution);
    [vertex2] = get_electrode_position(el2, file_id_meas, resolution);

    % create adjacency matrix
    numVertices = size(vertices, 1);
    adjacencyMatrix = sparse(numVertices, numVertices);

    for i = 1:size(faces, 1)
        % get the indices of the vertices
        v = faces(i, :);
        for j = 1:3
            % Connect vertices in the face
            v1 = v(j);
            v2 = v(mod(j, 3) + 1); % connect v1 to v2, v2 to v3, v3 to v1
            dist = norm(vertices(v1, :) - vertices(v2, :)); % Euclidean distance
            adjacencyMatrix(v1, v2) = dist;
            adjacencyMatrix(v2, v1) = dist;
        end
    end

    % Create graph from adjacency matrix
    G = graph(adjacencyMatrix);

    % Calculate geodesic distance (shortest path)
    distance = distances(G, vertex1, vertex2);
end
