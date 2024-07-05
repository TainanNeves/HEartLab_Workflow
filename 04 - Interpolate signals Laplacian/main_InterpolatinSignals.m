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

clear i Fsampling interpolated_signals field_name data


%% Sincrhonization
S = 8; % Multiplication factor
Fsampling = 4000; % Sampling frequency

% MEA 1
D_3D = interpolated.case1;
[numX, numY, ~] = size(D_3D);
D_SYNC.case1 = zeros(numX, numY, round(Fsampling * S)+1);
for x = 1:numX
    for y = 1:numY
        channel = squeeze(D_3D(x, y, :)); % Extract signal for electrode at position (x, y)
        [~, ido] = min(abs(D_EL.Timestamps - D_EL.TTL(1, 1)));
        ido = round(ido + (Fsampling * 1));
        idf = round(ido + (Fsampling * S));
        % Extract the signal segment from ido to idf
        signal_segment = channel(ido:idf)';
        % Store the extracted segment
        D_SYNC.case1(x, y, :) = signal_segment;
    end
end

% MEA 2
D_3D = interpolated.case2;
[numX, numY, ~] = size(D_3D);
D_SYNC.case2 = zeros(numX, numY, round(Fsampling * S)+1);
for x = 1:numX
    for y = 1:numY
        channel = squeeze(D_3D(x, y, :)); % Extract signal for electrode at position (x, y)
        [~, ido] = min(abs(D_EL.Timestamps - D_EL.TTL(1, 1)));
        ido = round(ido + (Fsampling * 1));
        idf = round(ido + (Fsampling * S));
        % Extract the signal segment from ido to idf
        signal_segment = channel(ido:idf)';
        % Store the extracted segment
        D_SYNC.case2(x, y, :) = signal_segment;
    end
end

% MEA 3
D_3D = interpolated.case3;
[numX, numY, ~] = size(D_3D);
D_SYNC.case3 = zeros(numX, numY, round(Fsampling * S)+1);
for x = 1:numX
    for y = 1:numY
        channel = squeeze(D_3D(x, y, :)); % Extract signal for electrode at position (x, y)
        [~, ido] = min(abs(D_EL.Timestamps - D_EL.TTL(1, 1)));
        ido = round(ido + (Fsampling * 1));
        idf = round(ido + (Fsampling * S));
        % Extract the signal segment from ido to idf
        signal_segment = channel(ido:idf)';
        % Store the extracted segment
        D_SYNC.case3(x, y, :) = signal_segment;
    end
end

% TANK
D_3D = interpolated.case4;
[numX, numY, ~] = size(D_3D);
D_SYNC.case4 = zeros(numX, numY, round(Fsampling * S)+1);
for x = 1:numX
    for y = 1:numY
        channel = squeeze(D_3D(x, y, :)); % Extract signal for electrode at position (x, y)
        [~, ido] = min(abs(D_EL.Timestamps - D_EL.TTL(1, 1)));
        ido = round(ido + (Fsampling * 1));
        idf = round(ido + (Fsampling * S));
        % Extract the signal segment from ido to idf
        signal_segment = channel(ido:idf)';
        % Store the extracted segment
        D_SYNC.case4(x, y, :) = signal_segment;
    end
end

% Cleaning (optional)
clear Fsampling S ans i idf ido S x y channel numX numY signal_segment D_3D


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
InterpSignal.Sync.MEA1 = D_SYNC.case1;
InterpSignal.Sync.MEA2 = D_SYNC.case2;
InterpSignal.Sync.MEA3 = D_SYNC.case3;
InterpSignal.Sync.TANK = D_SYNC.case4;

%Exporting
save(['InterpolatedSignals', FileName, '_filtered'], 'InterpSignal', '-v7.3');













