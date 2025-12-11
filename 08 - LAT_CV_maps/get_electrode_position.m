function [mea_vertex] = get_electrode_position(electrode, file_id_meas, resolution)
    % Determine which index set to use based on resolution
    if ~isnumeric(resolution)
        id_mea_set = file_id_meas.MEAS_HR_IDX;
    else
        id_mea_set = file_id_meas.(['MEAS_IDX_' num2str(resolution)]);
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

