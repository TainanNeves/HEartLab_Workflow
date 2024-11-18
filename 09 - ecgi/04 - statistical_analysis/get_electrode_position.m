function [mea, id_mea, electrode_mea] = get_electrode_position(electrode, meas_signal, file_id_meas, resolution)
    % Determine which index set to use based on resolution
    if ~isnumeric(resolution)
        id_mea_set = file_id_meas.MEAS_HR_IDX;
    else
        id_mea_set = file_id_meas.(['MEAS_IDX_' num2str(resolution)]);
    end

    % Determine the segment of meas_signal and the index
    if electrode < 17
        mea = meas_signal(1:16, :);
        id_mea = id_mea_set.MEA1;
        electrode_mea = electrode;
    elseif electrode < 33
        mea = meas_signal(17:32, :);
        id_mea = id_mea_set.MEA2;
        electrode_mea = electrode - 16;
    else
        mea = meas_signal(65:80, :);
        id_mea = id_mea_set.MEA3;
        electrode_mea = electrode - 64;
    end
end

