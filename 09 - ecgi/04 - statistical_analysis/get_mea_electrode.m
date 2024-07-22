function [mea, id_mea, electrode_mea] = get_mea_electrode(electrode, meas_signal, file_id_meas, resolution)
    if electrode < 17
        mea = meas_signal(1:16,:);
        id_mea = file_id_meas.(['MEAS_IDX_' num2str(resolution)]).MEA1;
        electrode_mea = electrode;
    elseif electrode < 33
        mea = meas_signal(17:32,:);
        id_mea = file_id_meas.(['MEAS_IDX_' num2str(resolution)]).MEA2;
        electrode_mea = electrode - 16;
    else
        mea = meas_signal(65:80,:);
        id_mea = file_id_meas.(['MEAS_IDX_' num2str(resolution)]).MEA3;
        electrode_mea = electrode - 64;
    end
end
