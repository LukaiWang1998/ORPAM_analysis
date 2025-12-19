function [raw_data1_filtered] = orpam_filter(raw_data1, fs, min_freq, max_freq)
% raw_data1: 3D array (time x N x M)
% fs: sampling frequency in Hz
% min_freq, max_freq: frequency bounds for bandpass filtering
% raw_data1_filtered: filtered signal, same size as raw_data1

[data_len, N, M] = size(raw_data1);
n_fft = 2^nextpow2(data_len);      % Zero-padding length
freq = fs * (0:n_fft-1) / n_fft;   % Full frequency axis

% Get indices for frequency passband
passband_mask = (freq >= min_freq) & (freq <= max_freq);
passband_mask = passband_mask';
% passband_mask(n_fft:-1:n_fft/2+2) = passband_mask(2:n_fft/2);  % Mirror for symmetry

raw_data1_filtered = zeros(size(raw_data1));

for j = 1:M
    for i = 1:N
        signal = double(reshape(raw_data1(:, i, j), [], 1));  % Ensure column vector
        signal_mean = mean(signal);     % Store original mean (DC component)
        signal = signal - signal_mean;  % Remove DC

        Y = fft(signal, n_fft);         % FFT
        Y_filtered = Y .* passband_mask; % Apply bandpass filter

        signal_filtered = real(ifft(Y_filtered));  % Inverse FFT

        % Restore original mean (DC component)
        signal_filtered = signal_filtered + signal_mean;

        raw_data1_filtered(:, i, j) = signal_filtered(1:data_len);  % Crop to original length
    end
end

end


% function [raw_data1_filtered] = orpam_filter(raw_data1, fs, min_freq, max_freq)
% % raw_data1: 3D array (time x N x M)
% % fs: sampling frequency in Hz
% % min_freq, max_freq: frequency bounds for bandpass filtering
% % raw_data1_filtered: filtered signal, same size as raw_data1
% 
% [data_len, N, M] = size(raw_data1);
% n_fft = 2^nextpow2(data_len);            % Zero-padding length
% freq = fs * (0:n_fft-1) / n_fft;         % Full frequency axis
% passband_mask = (freq >= min_freq) & (freq <= max_freq);
% passband_mask = passband_mask(:);        % Ensure column vector
% 
% raw_data1_filtered = zeros(size(raw_data1), 'like', raw_data1);
% 
% parfor j = 1:M
%     temp_slice = zeros(data_len, N);
%     for i = 1:N
%         signal = double(raw_data1(:, i, j));       % Extract A-line
%         signal = signal - mean(signal);            % Remove DC
%         Y = fft(signal, n_fft);                    % FFT
%         Y_filtered = Y .* passband_mask;           % Apply filter
%         signal_filtered = real(ifft(Y_filtered));  % Inverse FFT
%         signal_filtered = signal_filtered(1:data_len);  % Crop to original length
%         signal_filtered = signal_filtered + mean(signal);  % Restore DC
%         temp_slice(:, i) = signal_filtered;
%     end
%     raw_data1_filtered(:, :, j) = temp_slice;
% end
% end
