function [] = plot_electric_phasemap(phase_el, lim, frame, caso)
% PLOT_ELECTRIC_PHASEMAP - Create visual phase maps from MEAs or Tank electrodes.
%
% Syntax:
%   plot_electric_phasemap(V_interpolated, lim, frame, caso)
%
% Description:
%   This function generates visual phase maps for Multi-Electrode Arrays (MEAs) or Tank
%   electrodes based on the given phase data, electrode layout, and filtering.
%
% Input:
%   - V_interpolated: Phase data for electrodes over time (rows: electrodes, columns: phase).
%   - lim: Color axis limits for the phase map.
%   - frame: Frame index for visualization.
%   - caso: Case identifier (1-4) specifying the electrode configuration.
%       -1 - MEA 1
%       -2 - MEA 2
%       -3 - MEA 3
%       -4 - TANK
%
% Author:
%   Tainan Neves, HEartLab - UFABC
%
% Example:
%   plot_phasemap(data, [-pi pi], 50, 3);

switch caso
    case 1
        % MEA1 - Right Atrium
        
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA1
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        % Plot MEA1
        plotMEA(phase_el, MEA, MEA_plane_indx, 'MEA 1 (Right Atrium)', frame, lim);
        
    case 2
        % MEA2 - Ventricle
        
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA2
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        % Plot MEA2
        plotMEA(phase_el, MEA, MEA_plane_indx, 'MEA 2 (Ventricle)', frame, lim);

    case 3
        % MEA3 - Left Atrium
        
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA3
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
               
        % Plot MEA3
        plotMEA(phase_el, MEA, MEA_plane_indx, 'MEA 3 (Left Atrium)', frame, lim);

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
             
        % Plot the tank
        plotTank(phase_el, Plane, tank_plane_indx, frame, lim);

end

end