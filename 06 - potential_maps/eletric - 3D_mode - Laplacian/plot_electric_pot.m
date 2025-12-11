function plot_electric_pot (Data, lim, sample, caso)
% PLOT_ELECTRIC_POT - Visualize electric potential data for different electrode configurations.
%
% Usage:
%   plot_electric_pot(Data, lim, sample, caso)
%
% Inputs:
%   - Data: Matrix containing electric potential data. Rows represent electrodes,
%           and columns represent samples.
%   - lim: Two-element vector specifying the y-axis limits for the plot.
%   - sample: Sample indices for plotting.
%   - caso: Integer indicating the electrode configuration (1-4).
%
% Output:
%   The function produces a plot visualizing electric potential data based on the
%   specified electrode configuration.
%
% Configuration Options (caso):
%   1: MEA1 - Right Atrium
%   2: MEA2 - Ventricle
%   3: MEA3 - Left Atrium
%   4: TANK

switch caso
    case 1 % MEA1 - Right Atrium
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA1
        MEA1_list_indx = [13 14 15 16, ...
                        9 10 11 12, ...
                        5 6 7 8, ...
                        1 2 3 4];
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        % Laplacian interpolation for MEA1
        [lap, ~] = mesh_laplacian(MEA.vertices, MEA.faces);
        MEA1_plane_indx = MEA_plane_indx;
        MEA1_plane_indx([1 8]) = []; % Exclude bad electrodes
        [int1, ~, ~] = mesh_laplacian_interp(lap, MEA1_plane_indx);
        
        % Interpolate and filter data for MEA1
        dataFiltered = Data;
        MEA1_list_indx([1 8]) = [];%bad electrodes
        V_MEA1 = int1 * sgolayfilt(dataFiltered(MEA1_list_indx,:),2,201,[],2);
        
        % Plot MEA1
        plotMEA(V_MEA1, sample, MEA, MEA_plane_indx, 'MEA 1 (Right Atrium)', lim);

    case 2 % MEA2 - Ventricle
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA2
        MEA2_list_indx = [17 21 25 29, ...
                        18 22 26 30, ...
                        19 23 27 31, ...
                        20 24 28 32];
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        % Laplacian interpolation for MEA2
        [lap, ~] = mesh_laplacian(MEA.vertices, MEA.faces);
        MEA2_plane_indx = MEA_plane_indx;
        [int2, ~, ~] = mesh_laplacian_interp(lap,MEA2_plane_indx);

        
        % Interpolate and filter data for MEA2
        dataFiltered = Data;
        V_MEA2 = int2 * sgolayfilt(dataFiltered(MEA2_list_indx,:),2,201,[],2);

        % Plot MEA2
        plotMEA(V_MEA2, sample, MEA, MEA_plane_indx, 'MEA 2 (Ventricle)', lim);

    case 3 % MEA3 - Left Atrium
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA3
        MEA3_list_indx = [77 78 79 80, ...
                        73 74 75 76, ...
                        69 70 71 72, ...
                        65 66 67 68];
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        % Laplacian interpolation for MEA3
        [lap, ~] = mesh_laplacian(MEA.vertices, MEA.faces);
        MEA3_plane_indx = MEA_plane_indx;
        MEA3_plane_indx([4]) = [];%bad electrodes
        [int3, ~, ~] = mesh_laplacian_interp(lap,MEA3_plane_indx);

        % Interpolate and filter data for MEA3
        dataFiltered = Data;
        MEA3_list_indx([4]) = [];%bad electrodes
        V_MEA3 = int3 * sgolayfilt(dataFiltered(MEA3_list_indx,:),2,201,[],2);

        % Plot MEA3
        plotMEA(V_MEA3, sample, MEA, MEA_plane_indx, 'MEA 3 (Left Atrium)', lim);

    case 4 % TANK
        %loading the 3d tank geometry
        load('Tank_geometry.mat');
        
        %choosing the number of the tank electrodes
        tank_list_indx =[145 146 155 156 165 166 129 130 139 140 181 182, ...
                    147 157 167 131 141 183, ...
                    148 149 158 159 168 169 132 133 142 143 184 185, ...
                    150 151 160 161 170 171 134 135 144 177 186 187, ...
                    152 162 172 136 178 188, ...
                    153 154 163 164 173 174 137 138 179 180 189 190];
        %choosing the positions that represents each electrode
        tank_plane_indx =[77:2:99, ...
                    153:4:173, ...
                    227:2:249, ...
                    402:2:424, ...
                    478:4:498, ...
                    552:2:574];

        % TANK - Laplacian interpolator
        [lap,~] = mesh_laplacian(Plane.vertices,Plane.faces);
        [int, ~, ~] = mesh_laplacian_interp(lap,tank_plane_indx);

        % Tank Interpolation
        dataFiltered = Data;
        V_TANK = int * sgolayfilt(dataFiltered(tank_list_indx,:),2,201,[],2);
        
        % Plot the tank
        plotTank(V_TANK, sample, Plane, tank_plane_indx, 'TANK' ,lim);

end
end