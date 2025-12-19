% Modified from old BSCAN_ANALYSIS_batch.m (initial commit in
% https://github.com/LukaiWang1998/ORPAM_analysis)
clear;
load('info_ORPAM.mat');
data_total = 256;
ALinesPerBScan = 1000;
% Dynamic ranging threshold
dB = 40;
% n_avg_alines: number of A-lines to average on each side (half-width)
n_avg_alines = 3;
% Define cropping parameters
odd_top_crop = 120;
odd_bottom_crop = 100;
even_crop_end = ALinesPerBScan - odd_top_crop;
even_crop_start = even_crop_end - odd_top_crop - odd_bottom_crop + 1;

% Define the filter
f = [0 0.03 0.03 1]; m = [0 0.0 1 1]; DC_cancel = fir2(64,f,m); % HIGH PASS FILTER
f = [0 0.7 0.7 1]; m = [1 0.9 0.1 0]; HF_cancel = fir2(64,f,m); % LOW PASS FILTER
info.DC_cancel = DC_cancel; info.HF_cancel = HF_cancel;  

DataFolder = '\\10.229.121.108\DataArchive\ORPAM_Data\EndometrialCancer';
TargetFolder = 'C:\Users\Admin\Box\ORPAM_EC_DL\MATLAB_test';
Cases = getFolderNames(DataFolder);

for i = 1:numel(Cases)
    Dir = fullfile(DataFolder, Cases{i});
    S_temp = dir(fullfile(Dir,'*.bin'));
    N_bins = numel(S_temp);
    fprintf('>>>>>>>> FOUND %d SCANS\n', N_bins)

    %% SELECT SCAN AND LOAD DATA
    for sample_idx = 1:N_bins
        F = fullfile(Dir,S_temp(sample_idx).name);
        tokens = regexp(S_temp(sample_idx).name, '^(.*)\.bin$', 'tokens');
        filename = tokens{1}{1};
        SaveFolder = fullfile(TargetFolder,Cases{i},filename);
        % SaveFolder = fullfile('D:\ORPAM Sample\processed data',aqdate,filename);
        if exist(SaveFolder)
            files = getFileNames(SaveFolder);
            if any(contains(files, 'max_freq_map_avgAline'))
                fprintf('>>>> Skipping %s (already processed)\n', filename);
                continue;
            end
        else
            mkdir(SaveFolder);
        end
        fileInfo = dir(F);
        fileSizeBytes = fileInfo.bytes; % Get the file size in bytes
        BLinesPerCScan = floor(fileSizeBytes / 8 / ALinesPerBScan / data_total);
        fid = fopen(F,'r','b'); raw_data = fread(fid,[data_total,ALinesPerBScan*BLinesPerCScan] ,'double'); fclose(fid);
        raw_data1 = reshape(raw_data,[data_total,ALinesPerBScan,BLinesPerCScan]);

        clearvars raw_data;
        Noffset = 0;
        % Align alternating B-scans in a 3D OR-PAM volume
        raw_data1 = align_bscans(raw_data1, Noffset);
        if(mod(size(raw_data1,3), 2) == 1)
            raw_data1 = raw_data1(:,:,1:(size(raw_data1,3)-1));
        end

        raw_data1_odd = raw_data1(:,:,1:2:end);
        raw_data1_even = raw_data1(:,:,2:2:end);

        % Step 1: Crop
        raw_data1_odd_crop = raw_data1_odd(:, odd_top_crop+1:end-odd_bottom_crop, :);
        raw_data1_even_crop = raw_data1_even(:, even_crop_start:even_crop_end, :);

        % Step 2: Concatenate along the second (A-line) dimension
        raw_data1_crop_combined = cat(2, raw_data1_odd_crop, raw_data1_even_crop);  % size: 256 × N_combined × 1750

        % Step 3: Resize the second dimension to ALinesPerBScan
        % raw_data1_resized = zeros(size(raw_data1_odd,1), ALinesPerBScan, BLinesPerCScan);
        raw_data1_resized = imresize3( ...
            raw_data1_crop_combined, ...
            [size(raw_data1_odd,1), ALinesPerBScan, BLinesPerCScan] ...
            );

        fs = 200e6; % Sampling frequency (Hz)
        [freq, spectra_cube] = spectrum_map(raw_data1_resized, fs);
        [max_val_map, max_freq_map] = get_max_spectral_features(spectra_cube, freq, 50);
        [freq_avgAline, spectra_cube_avgAline] = spectrum_map_avgAline(raw_data1_resized, fs, n_avg_alines);
        [max_val_map_avgAline, max_freq_map_avgAline] = get_max_spectral_features(spectra_cube_avgAline, freq_avgAline, 50);
        %%
        save(fullfile(SaveFolder, 'raw_data1.mat'), 'raw_data1', '-v7.3');
        save(fullfile(SaveFolder, 'raw_data1_resized.mat'), 'raw_data1_resized', '-v7.3');
        clear raw_data1_resized raw_data1_crop_combined raw_data1_even_crop raw_data1_even raw_data1_odd_crop raw_data1_odd;        
        [im_plot, RF] = ORPAM_recon(raw_data1, info, dB, BLinesPerCScan);
        % Save the data
        save(fullfile(SaveFolder, 'im_plot.mat'), 'im_plot');
        save(fullfile(SaveFolder, 'RF.mat'), 'RF');
        save(fullfile(SaveFolder, 'spectra_cube.mat'), 'spectra_cube', '-v7.3');
        % Plot by the aspect ratio of the data
        fig = figure;
        imshow(im_plot',[0.0,1],'border','tight');
        colormap('hot');
        set(gca,'xtick','','ytick','','fontweight','bold');
        saveas(fig, fullfile(SaveFolder, 'im_plot.fig'));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, 'im_plot.png'));  % Save as PNG image
        % Save the im_plot by the default sapect ratio
        fig = figure;
        imagesc(im_plot);
        colormap('hot');
        title('im plot','fontsize',18);
        saveas(fig, fullfile(SaveFolder, 'im_plot_hot.fig'));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, 'im_plot_hot.png'));  % Save as PNG image

        % % Save the central frequency map
        % fig = figure;
        % imagesc(max_freq_map);
        % colorbar;
        % clim([5e5 6e6]);
        % title('max freq map','fontsize',18);
        % saveas(fig, fullfile(SaveFolder, 'max_freq_map.fig'));  % Save as MATLAB figure
        % saveas(fig, fullfile(SaveFolder, 'max_freq_map.png'));  % Save as PNG image
        % % Save the central frequency value map
        % fig = figure;
        % imagesc(max_val_map);
        % title('max val map','fontsize',18);
        % saveas(fig, fullfile(SaveFolder, 'max_val_map.fig'));  % Save as MATLAB figure
        % saveas(fig, fullfile(SaveFolder, 'max_val_map.png'));  % Save as PNG image
        % 
        % 
        % 
        % % Save the central frequency map
        % fig = figure;
        % imagesc(max_freq_map_avgAline);
        % colorbar;
        % clim([5e5 6e6]);
        % title('max freq map avgAline','fontsize',18);
        % saveas(fig, fullfile(SaveFolder, ['max_freq_map_avgAline_',num2str(n_avg_alines),'.fig']));  % Save as MATLAB figure
        % saveas(fig, fullfile(SaveFolder, ['max_freq_map_avgAline_',num2str(n_avg_alines),'.png']));  % Save as PNG image
        % % Save the central frequency value map
        % fig = figure;
        % imagesc(max_val_map_avgAline);
        % title('max val map avgAline','fontsize',18);
        % saveas(fig, fullfile(SaveFolder, ['max_val_map_avgAline_',num2str(n_avg_alines),'.fig']));  % Save as MATLAB figure
        % saveas(fig, fullfile(SaveFolder, ['max_val_map_avgAline_',num2str(n_avg_alines),'.png']));  % Save as PNG image
        % close all;
        clear im_mip_odd im_mip_odd im_mip_even_correct im_mip_odd_correct raw_data1 raw_data1_even raw_data1_odd raw_data1_odd_crop raw_data1_even_crop RF_DAS_env RF_raw_all RF_Sum_all spectra_cube spectra_cube_avgAline log_DAS_data max_freq_map max_freq_map_avgAline max_val_map max_val_map_avgAline;
    end
end