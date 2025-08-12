function distance = geodesicDistance(el1, el2, projections, resolution, scaleX, scaleY, scaleZ)
    %
    % Author: Angélica Quadros
    % Affiliation: HeartLab, UFABC
    % Year: 2024
    %
    % Description:
    % This function computes the geodesic distance between two electrodes on a 
    % 3D heart geometry mesh using the shortest path method. It supports different 
    % resolutions and optional scaling of mesh coordinates.
    %
    % Inputs:
    %   el1, el2       - Numbers of the electrodes. The following values are expected:
    %                    * MEA1 (RA): Electrodes 1 to 16
    %                    * MEA2 (V): Electrodes 17 to 32
    %                    * MEA3 (LA): Electrodes 65 to 80
    %   projections    - Struct containing the following:
    %                    * Heart geometries at different resolutions: 1200, 2500,
    %                      10000, 20000, or 'HR' (high resolution).
    %                    * Mapping of MEA electrode numbers to vertices on the geometry.
    %   resolution     - Resolution of the geometry ('HR' for high-resolution or one of 
    %                    the numeric values: 1200, 2500, 10000, 20000).
    %   scaleX, scaleY, scaleZ - (Optional) Scaling factors for x, y, and z dimensions.
    %                           Use these only if the geometry unit is unknown or not 
    %                           in cm, and you want to manually set the size.
    %
    % Outputs:
    %   distance       - Geodesic distance between the two electrodes
    %                    (defalut is in cm)
    %
    % Usage:
    %   To calculate the geodesic distance between two electrodes on a heart geometry:
    %
    %   Example 1: Without scaling
    %   distance = geodesicDistance(1, 17, projections, 'HR');
    %
    %   Example 2: With scaling factors for x, y, and z dimensions
    %   distance = geodesicDistance(1, 17, projections, 2500, 1.2, 1.2, 1.1);
    %
    %   Here, electrode 1 (RA) and electrode 17 (V) are used. The 'projections' 
    %   struct should include the geometry and mapping information, while 
    %   scaling factors adjust the unit dimensions if necessary.
    %   



    if ~isnumeric(resolution)
        vertices = projections.geometry_HR.vertices;
        faces = projections.geometry_HR.faces;
  
    else
        vertices = projections.(['geometry_' num2str(resolution)]).vertices;
        faces = projections.(['geometry_' num2str(resolution)]).faces;
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

    [vertex1] = get_electrode_position(el1, projections, resolution);
    [vertex2] = get_electrode_position(el2, projections, resolution);

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
