function log_DAS_data = dynamic_ranging_log_compression(RF_DAS_env, dB_PEdas)
% Apply log compression with dynamic range limitation to RF_DAS_env
% Inputs:
%   RF_DAS_env : normalized DAS envelope (3D array)
%   dB_PEdas   : dynamic range in dB (e.g., 75)
% Output:
%   log_DAS_data : log-compressed and filtered envelope volume

[~, ~, BLinesPerCScan] = size(RF_DAS_env);
log_DAS_data = zeros(size(RF_DAS_env));
min_dBdas = 10^(-dB_PEdas/20);

fwb = waitbar(0, 'Dynamic ranging ... ');

for i = 1:BLinesPerCScan
    das_tmp = squeeze(RF_DAS_env(:,:,i));
    idx_tmp = das_tmp < min_dBdas;

    das_tmp = (20/dB_PEdas) * log10(das_tmp) + 1;  % Log compression
    das_tmp(idx_tmp) = 0;                         % Clip low-intensity signals
    das_tmp = medfilt2(das_tmp, [3, 3]);          % Median filter

    log_DAS_data(:,:,i) = das_tmp;

    perc = i / BLinesPerCScan;
    msg = sprintf('Dynamic ranging ... %0.3f', perc);
    waitbar(perc, fwb, msg);
end

close(fwb);
disp('>>>>>>>> LOG COMPRESSION COMPLETE.');
end