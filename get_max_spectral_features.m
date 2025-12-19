% function [max_val, max_freq] = get_max_spectral_features(spectra_cube, freq)
% % spectra_cube: Nf x N x M array (frequency x spatial dimensions)
% % freq: Nf x 1 array of frequency values (in Hz or MHz)
% % max_val: N x M matrix of max values in spectra
% % max_freq: N x M matrix of frequency (from freq) corresponding to max_val
% 
%     % Find max value and index along the frequency dimension
%     [max_val, max_idx] = max(spectra_cube, [], 1);  % size: 1 x N x M
%     max_val = squeeze(max_val);                     % size: N x M
%     max_idx = squeeze(max_idx);                     % size: N x M
% 
%     % Map indices to frequencies
%     max_freq = freq(max_idx);                       % size: N x M
% end

function [max_val_map, max_freq_map] = get_max_spectral_features(spectra_cube, freq, dB_vas)
% spectra_cube: Nf x N x M array (frequency x spatial dimensions)
% freq: Nf x 1 array of frequency values (in Hz or MHz)
% dB_vas: scalar threshold in dB (e.g., 40)
% max_val_map: N x M matrix of max values in linear scale (set to 0 if < threshold)
% max_freq_map: N x M matrix of frequency corresponding to max_val (set to 0 if < threshold)

    % Get max value and index along frequency dimension
    [max_val, max_idx] = max(spectra_cube, [], 1);  % size: 1 x N x M
    max_val = squeeze(max_val);                     % size: N x M
    max_val = max_val / max(max_val(:));
    max_idx = squeeze(max_idx);                     % size: N x M

    % Convert max_val to dB
    max_val_dB = 20 * log10(max_val + eps);         % Add eps to avoid log(0)

    % Create masks for valid pixels
    mask = max_val_dB >= -dB_vas;

    % Prepare outputs
    max_val_map = zeros(size(max_val));
    max_freq_map = zeros(size(max_val));

    % Assign values only where valid
    max_val_map(mask) = max_val(mask);
    max_freq_map(mask) = freq(max_idx(mask));
end
