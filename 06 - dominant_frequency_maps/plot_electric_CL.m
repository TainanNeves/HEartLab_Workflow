function [] = plot_electric_CL(CL, lim, caso)

switch caso
    case 1
        % MEA1 - Right Atrium
        
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
        dataFiltered = CL';
        MEA1_list_indx([1 8]) = []; % Exclude bad electrodes
        CL_MEA1 = int1 * dataFiltered(MEA1_list_indx,:);
        
        % Plot MEA1
        plotMEA_CL(CL_MEA1, MEA, MEA_plane_indx, 'MEA 1 (Right Atrium)', lim);


    case 2
        % MEA2 - Ventricle
        
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA1
        MEA2_list_indx = [17 21 25 29, ...
                        18 22 26 30, ...
                        19 23 27 31, ...
                        20 24 28 32];
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        % Laplacian interpolation for MEA1
        [lap, ~] = mesh_laplacian(MEA.vertices, MEA.faces);
        MEA2_plane_indx = MEA_plane_indx;
        [int2, ~, ~] = mesh_laplacian_interp(lap, MEA2_plane_indx);
        
        % Interpolate and filter data for MEA2
        dataFiltered = CL';
        CL_MEA2 = int2 * dataFiltered(MEA2_list_indx,:);
        
        % Plot MEA1
        plotMEA_CL(CL_MEA2, MEA, MEA_plane_indx, 'MEA 2 (Ventricle)', lim);


    case 3
        % MEA3 - Left Atrium
        
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA1
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
        MEA3_plane_indx([4]) = []; % Exclude bad electrodes
        [int3, ~, ~] = mesh_laplacian_interp(lap, MEA3_plane_indx);
        
        % Interpolate and filter data for MEA1
        dataFiltered = CL';
        MEA3_list_indx([4]) = []; % Exclude bad electrodes
        CL_MEA3 = int3 * dataFiltered(MEA3_list_indx,:);
        
        % Plot MEA3
        plotMEA_CL(CL_MEA3, MEA, MEA_plane_indx, 'MEA 3 (Left Atrium)', lim);


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
        dataFiltered = CL';
        CL_TANK = intTANK * dataFiltered(tank_list_indx,:);
        
        % Plot the tank
        plotTank_CL(CL_TANK, Plane, tank_plane_indx, 'TANK' ,lim);


end
end