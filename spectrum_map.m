function [freq, spectra_cube] = spectrum_map(raw_data1, fs)
% raw_data1: 3D array of size (256 x N x M), where 256 is the time axis
% fs: sampling frequency in Hz (e.g., 20e6 for 20 MHz)
% freq: frequency axis
% spectra_cube: 3D array of one-sided PSDs, size (n_freq x N x M)

[data_len, N, M] = size(raw_data1);
n_fft = 2^nextpow2(data_len);  % Zero-padding for FFT efficiency
n_freq = n_fft/2 + 1;          % Number of one-sided frequency bins

spectra_cube = zeros(n_freq, N, M);

for j = 1:M
    for i = 1:N
        signal = double(squeeze(raw_data1(:, i, j)));
        signal = signal - mean(signal);  % Remove DC
        Y = fft(signal, n_fft);
        P2 = abs(Y / n_fft).^2;          % Two-sided power spectrum
        P1 = P2(1:n_freq);
        P1(2:end-1) = 2 * P1(2:end-1);   % Convert to one-sided
        spectra_cube(:, i, j) = P1;
    end
end

% Frequency axis
freq = fs * (0:(n_fft/2)) / n_fft;

end
