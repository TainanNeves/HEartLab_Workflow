function [mea_vertex] = get_electrode_position(electrode, projections, resolution)
%
% Author: Angélica Quadros
% Affiliation: HeartLab, UFABC
% Year: 2024
%
% Description:
% This function maps a specified MEA electrode number to its corresponding 
% vertex on the heart geometry mesh.
%
% Inputs:
%   electrode  - The number of the electrode (1 to 16 for MEA1, 17 to 32 for MEA2, 
%                or 65 to 80 for MEA3).
%   projections    - Struct containing the following:
%                    * Heart geometries at different resolutions: 1200, 2500,
%                      10000, 20000, or 'HR' (high resolution).
%                    * Mapping of MEA electrode numbers to vertices on the geometry.
%   resolution     - Resolution of the geometry ('HR' for high-resolution or one of 
%                    the numeric values: 1200, 2500, 10000, 20000).
%
% Outputs:
%   mea_vertex - The vertex number on the heart geometry mesh corresponding 
%                to the specified electrode.
%
% Usage:
%   vertex = get_electrode_position(5, projections, 'HR');
%   This example retrieves the vertex corresponding to electrode 5 in the 
%   high-resolution geometry.
%
%
    % Determine which index set to use based on resolution
    if ~isnumeric(resolution)
        id_mea_set = projections.MEAS_HR_IDX;
    else
        id_mea_set = projections.(['MEAS_IDX_' num2str(resolution)]);
    end

    % Determine the segment of meas_signal and the index
    if electrode < 17
        id_mea = id_mea_set.MEA1;
        electrode_mea = electrode;
        mea_vertex = id_mea(electrode_mea);
    elseif electrode < 33
        id_mea = id_mea_set.MEA2;
        electrode_mea = electrode - 16;
        mea_vertex = id_mea(electrode_mea);
    else
        id_mea = id_mea_set.MEA3;
        electrode_mea = electrode - 64;
        mea_vertex = id_mea(electrode_mea);
    end
end



