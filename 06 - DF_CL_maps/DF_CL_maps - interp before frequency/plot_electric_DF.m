function [] = plot_electric_DF(MFFTi, lim, caso)

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
        plotMEA_DF(MFFTi, MEA, MEA_plane_indx, 'MEA 1 (Right Atrium)', lim);


    case 2
        % MEA2 - Ventricle
        
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA1
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        
        % Plot MEA1
        plotMEA_DF(MFFTi, MEA, MEA_plane_indx, 'MEA 2 (Ventricle)', lim);


    case 3
        % MEA3 - Left Atrium
        
        % Load 3D MEAs geometry
        load('MEA.mat');
        
        % Electrode indices and positions for MEA1
        MEA_plane_indx = [25:2:31, ...
                        47:2:53, ...
                        69:2:75, ...
                        91:2:97];
        
        % Plot MEA3
        plotMEA_DF(MFFTi, MEA, MEA_plane_indx, 'MEA 3 (Left Atrium)', lim);


    case 4
        % Tank
        
        % Load 3D tank geometry
        load('Tank_geometry.mat');
        
        % Electrode indices and positions for the tank
        tank_plane_indx = [77:2:99, ...
            153:4:173, ...
            227:2:249, ...
            402:2:424, ...
            478:4:498, ...
            552:2:574];
       
        % Plot the tank
        plotTank_DF(MFFTi, Plane, tank_plane_indx, 'TANK' ,lim);


end
end