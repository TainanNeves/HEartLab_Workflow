function V_interpolated = electric_interp_2025(dataFiltered, caso)
% ELECTRIC_INTERP - Interpolate electrical data from MEAs or Tank electrodes.
%
% Syntax:
%   V_interpolated = electric_interp(dataFiltered, caso)
%
% Description:
%   This function performs Laplacian interpolation on the electrical data recorded
%   from Multi-Electrode Arrays (MEAs) or Tank electrodes. The function returns
%   the interpolated and filtered data based on the specified electrode configuration.
%
% Input:
%   - dataFiltered: Electrode signals (rows: electrodes, columns: time points).
%   - caso: Case identifier (1-4) specifying the electrode configuration:
%       - 1: MEA 1 (Left Atria)
%       - 2: MEA 2 (Right Atria)
%       - 3: MEA 3 (Ventricle)
%       - 4: Tank
%
% Output:
%   - V_interpolated: Interpolated and filtered data for the specified electrode configuration.
%
% Author:
%   Angélica Quadros, Tainan Neves
%   HEartLab - UFABC
%
% Example:
%   dataFiltered = Data(:, start_sample:end_sample);
%   caso = 1; % Using MEA 1 configuration
%   V_interpolated = electric_interp(dataFiltered, caso);
%

switch caso
    case 1
        % MEA1 - Left Atria
        
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA1
        MEA1_list_indx = [4 3 2 1, ...
                        8 7 6 5, ...
                        12 11 10 9, ...
                        16 15 14 13];
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        % Laplacian interpolation for MEA1
        [lap, ~] = mesh_laplacian(MEA.vertices, MEA.faces);
        MEA1_plane_indx = MEA_plane_indx;
        MEA1_plane_indx([]) = []; % Exclude bad electrodes
        [int1, ~, ~] = mesh_laplacian_interp(lap, MEA1_plane_indx);
        
        % Interpolate and filter data for MEA1
%         MEA1_list_indx([1 8]) = []; % Exclude bad electrodes
        V_interpolated = int1 * sgolayfilt(dataFiltered(MEA1_list_indx,:), 2, 201, [], 2);
        
        
    case 2
        % MEA2 - RA
        
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA2
        MEA2_list_indx = [20 19 18 17, ...
                        24 23 22 21, ...
                        28 27 26 25, ...
                        32 31 30 29];
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        % Laplacian interpolation for MEA2
        [lap, ~] = mesh_laplacian(MEA.vertices, MEA.faces);
        MEA2_plane_indx = MEA_plane_indx;
        [int2, ~, ~] = mesh_laplacian_interp(lap, MEA2_plane_indx);
        
        % Interpolate and filter data for MEA2
        V_interpolated = int2 * sgolayfilt(dataFiltered(MEA2_list_indx,:), 2, 201, [], 2);
        
        
    case 3
        % MEA3 - Ventricle
        
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA3
        MEA3_list_indx = [84 83 82 81, ...
                        88 87 86 85, ...
                        92 91 90 89, ...
                        96 95 94 93];
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        % Laplacian interpolation for MEA3
        [lap, ~] = mesh_laplacian(MEA.vertices, MEA.faces);
        MEA3_plane_indx = MEA_plane_indx;
        MEA3_plane_indx([]) = []; % Exclude bad electrodes
        [int3, ~, ~] = mesh_laplacian_interp(lap, MEA3_plane_indx);
        
        % Interpolate and filter data for MEA3
%         MEA3_list_indx([4]) = []; % Exclude bad electrodes
        V_interpolated = int3 * sgolayfilt(dataFiltered(MEA3_list_indx,:), 2, 201, [], 2);
        

    case 4
        % Tank
        
        % Load 3D tank geometry
        load('Tank_geometry.mat');
        
        % Electrode indices and positions for the tank
        tank_list_indx = [145 146 155 156 165 166 129 130 139 140 181 182, ...
            147 157 167 131 141 183, ...
            148 149 158 159 168 169 132 133 142 143 184 185, ...
            150 151 160 161 170 171 134 135 144 177 186 187, ...
            152 162 172 136 178 188, ...
            153 154 163 164 173 174 137 138 179 180 189 190];
        tank_plane_indx = [77:2:99, ...
            153:4:173, ...
            227:2:249, ...
            402:2:424, ...
            478:4:498, ...
            552:2:574];
        
        % Laplacian interpolation for the tank
        [lap, ~] = mesh_laplacian(Plane.vertices, Plane.faces);
        [intTANK, ~, ~] = mesh_laplacian_interp(lap, tank_plane_indx);
        
        % Interpolate and filter data for the tank
        V_interpolated = intTANK * sgolayfilt(dataFiltered(tank_list_indx,:), 2, 201, [], 2);
        
end

end
