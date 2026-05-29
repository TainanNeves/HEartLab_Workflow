function Cn = estimate_Cn(y, baseline_indices, diagonalOnly)
% estimate_Cn - Estimate noise covariance matrix Cn from baseline data.
%
% INPUTS:
%   y               : Matrix of torso signals [N_electrodes x N_samples].
%   baseline_indices: Vector of time indices corresponding to a baseline segment.
%   diagonalOnly    : (Optional) If true, return only a diagonal Cn (white noise assumption).
%                     Default: false.
%
% OUTPUT:
%   Cn : Noise covariance matrix [N_electrodes x N_electrodes].

    if nargin < 3
        diagonalOnly = false;
    end

    % Extract baseline data
    y_baseline = y(:, baseline_indices);

    % Remove mean (DC offset) for each electrode
    y_zero_mean = y_baseline - mean(y_baseline, 2);

    % Compute full covariance matrix
    Cn_full = (1 / size(y_zero_mean, 2)) * (y_zero_mean * y_zero_mean.');

    if diagonalOnly
        % Use only variance of each electrode (diagonal approximation)
        Cn = diag(diag(Cn_full));
    else
        % Return full covariance matrix
        Cn = Cn_full;
    end
end
