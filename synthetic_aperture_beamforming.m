function [RF_raw_all, RF_Sum_all] = synthetic_aperture_beamforming(raw_data1, info, h, N, N_repeat)
% Perform synthetic aperture beamforming on a 3D raw data volume
% Inputs:
%   raw_data1   : 3D array (depth x A-lines x B-scans)
%   info        : struct containing parameters (see main script)
%   HF_cancel   : high-pass FIR filter coefficients
%   DC_cancel   : low-pass FIR filter coefficients
%   h           : deconvolution kernel
%   N           : size of Gaussian mask (scalar)
%   info.Nzeropad    : zero-padding before time zero
%   N_repeat    : number of beamforming iterations
% Outputs:
%   RF_raw_all  : raw beamformed signals
%   RF_Sum_all  : final processed envelope

[depth, ALinesPerBScan, BLinesPerCScan] = size(raw_data1);

RF_raw_all = zeros(size(raw_data1));
RF_Sum_all = zeros(size(raw_data1));

z_focus = info.t_focus * info.c;

fwd = waitbar(0, 'Synthetic aperture beamforming ...');

for idx_bscan = 1:BLinesPerCScan
    bscan = squeeze(raw_data1(:,:,idx_bscan));
    bscan = conv2(info.HF_cancel, 1, bscan, 'same');

    info.z_sample = (1:(info.Nzeropad + depth)) * info.pixel_z;
    info.z_sample = info.z_sample(info.Nzeropad+1:end);
    RF_Sum = bscan;

    for idx_order = 1:N_repeat
        RF_Sum_tmp = RF_Sum;
        bf_radius_tmp = round((info.N_ele_bf - 1) / 2 / idx_order);

        parfor i = 1:ALinesPerBScan
            scan_idx_left = max(1, i - bf_radius_tmp);
            scan_idx_right = min(ALinesPerBScan, i + bf_radius_tmp);
            L_lr = scan_idx_right - scan_idx_left + 1;
            das_raw_data = zeros(depth, L_lr);

            for j = scan_idx_left:scan_idx_right
                rf_fil = RF_Sum(:, j);
                x = info.xstep * abs(j - i);
                R_z = z_focus + sign(info.z_sample - z_focus) .* sqrt(x^2 + (z_focus - info.z_sample).^2);
                TOF = R_z / info.c;
                rxReadPtr = round(TOF * info.fs) - info.Nzeropad;
                rxReadPtr = max(min(rxReadPtr, depth), 1);
                rf_tmp = rf_fil(rxReadPtr);
                rf_tmp(rxReadPtr > depth | rxReadPtr < 1) = 0;
                das_raw_data(:, j - scan_idx_left + 1) = rf_tmp;
            end

            W_tmp = info.SIR(:, bf_radius_tmp - i + scan_idx_left + 1 : bf_radius_tmp - i + scan_idx_right + 1);
            das_raw_data = das_raw_data .* W_tmp;
            area_z = sum(sum(W_tmp, 2), 2);
            area_z = max(area_z, 1 + 1e-2);

            das_raw_data2d = reshape(das_raw_data, depth, []);
            rf_sum0 = sum(das_raw_data2d, 2);
            rf_sum_squared = (rf_sum0 ./ area_z).^2;
            rf_squared_sum = sum(das_raw_data2d.^2, 2) ./ area_z;
            CF_das = rf_sum_squared ./ rf_squared_sum;
            CF_var_das = sqrt(CF_das ./ (1 + 1e-3 - CF_das));
            CF_var_das(isnan(CF_var_das)) = 0;
            RF_Sum_tmp(:, i) = rf_sum0 .* CF_var_das;
        end
        RF_Sum = RF_Sum_tmp;
    end

    RF_Sum = conv2(info.DC_cancel, 1, RF_Sum, 'same');
    RF_raw = RF_Sum;
    RF_Sum = abs(hilbert(RF_Sum));

    I_pad = padarray(RF_Sum, [floor(N/2), floor(N/2)]);
    I_pad = deconvlucy(wiener2(I_pad, [3, 3]), h, 3);
    RF_Sum = I_pad(floor(N/2)+1:end-floor(N/2), floor(N/2)+1:end-floor(N/2));
    RF_Sum(RF_Sum < 0) = 0;

    RF_raw_all(:,:,idx_bscan) = RF_raw;
    RF_Sum_all(:,:,idx_bscan) = RF_Sum;

    perc = idx_bscan / BLinesPerCScan;
    msg = sprintf('SAFT for B scan # %d / %d ... %0.2f', idx_bscan, BLinesPerCScan, perc);
    waitbar(perc, fwd, msg);
end

close(fwd);
end