%% Interpolatin Signals
% This code create a struct containing the interpolated signals for
% MEAs and Tank. Using a laplacian interpolation methodology based
% in a 3D surface .

%% Loading Electric signals
load(""); % Load the filtered electrical signals


%% Laplacian Interpolation
data = D_EL.Data;
Fsampling = D_EL.Header.sample_rate;

interpolated = struct();
for i = 1:4
    interpolated_signals = electric_interp(data, i);
    field_name = sprintf('case%d', i);
    interpolated.(field_name) = fill_matrix(interpolated_signals, i);
end

clear i Fsampling interpolated_signals field_name


%% Electrode positions
% Function that finds electrodes position in the matrix based on its number
[row, col] = getElectrodePosition(5);


%% Export Struct
% Filename to save
FileName = 'EXX_FXX_RXX';

% Fill Struct to export
InterpSignal.TTL = D_EL.TTL;
InterpSignal.Timestamp = D_EL.Timestamps;
InterpSignal.Opticalin = D_EL.opticalin;
InterpSignal.Header = D_EL.Header;
InterpSignal.Data.MEA1 = interpolated.case1;
InterpSignal.Data.MEA2 = interpolated.case2;
InterpSignal.Data.MEA3 = interpolated.case3;
InterpSignal.Data.TANK = interpolated.case4;

%Exporting
save(['InterpolatedSignals', FileName, '_filtered'], 'InterpSignal', '-v7.3');













