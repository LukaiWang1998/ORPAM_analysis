% Based on \\10.229.121.108\workspace\Lukai\ORPAM
% ANALYSIS\BSCAN_ANALYSIS_batch.m: Trying to improve the frequency analysis
% clc;
clear; close all;
TargetFolder = '\\10.229.121.108\Workspace\Lukai\Endometrium_Jun6';
TablePath = '\\10.229.121.108\Workspace\Lukai\EndometrialCancer\imagingRecord.xlsx';
imagingRecordTable = readtable(TablePath, 'Sheet', 'Apr18_processing', 'VariableNamingRule', 'preserve');
data_total = 256;
% BLinesPerCScan = 4000;  % Can be calculated from file size
ALinesPerBScan = 1000;
% dB_PE = 45;
fs = 20e6; % Sampling frequency (Hz)
min_freq = 468750;
max_freq = 5312500;
load('freq.mat');
% n_avg_alines: number of A-lines to average on each side (half-width)
n_avg_alines = 3;
% Define cropping parameters
odd_top_crop = 120;
odd_bottom_crop = 100;
even_crop_end = ALinesPerBScan - odd_top_crop;
even_crop_start = even_crop_end - odd_top_crop - odd_bottom_crop + 1;


filepath = '\\10.229.121.108\DataArchive\ORPAM_Data\EndometrialCancer';
% filepath = 'C:\Users\Admin\Box\ORPAM Data';

Cases = getFolderNames(filepath);
% Cases = {'20250505_McCourt_normal_49'};

for i = 1:size(Cases, 2)

    aqdate = Cases{i};
    parts = split(Cases{i}, '_');
    ID = parts{4};

    % Dir = 'D:\20250320_EC';
    Dir = fullfile(filepath, aqdate);

    S_temp = dir(fullfile(Dir,'*.bin'));
    N_bins = numel(S_temp);
    fprintf('>>>>>>>> FOUND %d SCANS\n', N_bins)

    %% SELECT SCAN AND LOAD DATA
    for sample_idx = 1:N_bins
        F = fullfile(Dir,S_temp(sample_idx).name);
        tokens = regexp(S_temp(sample_idx).name, '^(.*)\.bin$', 'tokens');
        filename = tokens{1}{1};
        row_idx = find(imagingRecordTable.("Patient ID") == str2double(ID) & ...
            strcmp(imagingRecordTable.Position, filename));

        % Display result
        if isempty(row_idx)
            disp("No match found.");
            continue;
        elseif length(row_idx) > 1
            warning("Multiple matches found:");
            continue;
        end


        SaveFolder = fullfile(TargetFolder,aqdate,filename);
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
        % raw_data1 = raw_data1(:,:,1:100); BLinesPerCScan = 100;
        clearvars raw_data;
        Noffset = 0;
        % Align alternating B-scans in a 3D OR-PAM volume
        raw_data1 = align_bscans(raw_data1, Noffset);
        if(mod(size(raw_data1,3), 2) == 1)
            raw_data1 = raw_data1(:,:,1:(size(raw_data1,3)-1));
        end

        % raw_data1 = raw_data1(:,:,1:50); BLinesPerCScan = 50;

        raw_data1_resized = reshape_3d(raw_data1, ...
            odd_top_crop, odd_bottom_crop, BLinesPerCScan);


        raw_data1_filtered = orpam_filter(raw_data1_resized, fs, min_freq, max_freq);

        save(fullfile(SaveFolder, 'raw_data1_resized.mat'), 'raw_data1_resized', '-v7.3');
        save(fullfile(SaveFolder, 'raw_data1_filtered.mat'), 'raw_data1_filtered', '-v7.3');
        clear raw_data1_resized raw_data1;

        % fs = 20e6; % Sampling frequency (Hz)

        
        [~, spectra_cube_filtered] = spectrum_map(raw_data1_filtered, fs);
        % spectra_cube_filtered = reshape_3d( ...
        %     spectra_cube_filtered, ...
        %     odd_top_crop, odd_bottom_crop, BLinesPerCScan);
        [max_val_map_filtered, max_freq_map_filtered] = get_max_spectral_features(spectra_cube_filtered, freq, 50);
        freq_spectra_filtered = squeeze(sum(sum(spectra_cube_filtered, 3), 2));
        figure; plot(freq, freq_spectra_filtered);
        mean_freq_filtered = sum(freq(:) .* freq_spectra_filtered(:)) / sum(freq_spectra_filtered(:));
        [~, max_idx] = max(freq_spectra_filtered);
        peak_freq_filtered = freq(max_idx);

        save(fullfile(SaveFolder, 'spectra_cube_filtered.mat'), 'spectra_cube_filtered', '-v7.3');
        save(fullfile(SaveFolder, 'max_val_map_filtered.mat'), 'max_val_map_filtered', '-v7.3');
        save(fullfile(SaveFolder, 'max_freq_map_filtered.mat'), 'max_freq_map_filtered', '-v7.3');
        save(fullfile(SaveFolder, 'freq_spectra_filtered.mat'), 'freq_spectra_filtered');
        save(fullfile(SaveFolder, 'mean_freq_filtered.mat'), 'mean_freq_filtered');
        save(fullfile(SaveFolder, 'peak_freq_filtered.mat'), 'peak_freq_filtered');

        fig = figure; 
        imagesc(max_freq_map_filtered); 
        colorbar;
        clim([0 max_freq]);
        saveas(fig, fullfile(SaveFolder, 'max_freq_map_filtered.fig'));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, 'max_freq_map_filtered.png'));  % Save as PNG image
        
        fig = figure; 
        imagesc(max_val_map_filtered); 
        colorbar;
        saveas(fig, fullfile(SaveFolder, 'max_val_map_filtered.fig'));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, 'max_val_map_filtered.png'));  % Save as PNG image
        
        fig = figure; 
        plot(freq, freq_spectra_filtered);
        saveas(fig, fullfile(SaveFolder, 'freq_spectra_filtered.fig'));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, 'freq_spectra_filtered.png'));  % Save as PNG image

        clear spectra_cube_filtered max_val_map_filtered max_freq_map_filtered freq_spectra_filtered mean_freq_filtered peak_freq_filtered;
        
        [~, spectra_cube_avgAline_filtered] = spectrum_map_avgAline(raw_data1_filtered, fs, n_avg_alines);
        % spectra_cube_avgAline_filtered = reshape_3d( ...
        %     spectra_cube_avgAline_filtered, ...
        %     odd_top_crop, odd_bottom_crop, BLinesPerCScan);
        [max_val_map_avgAline_filtered, max_freq_map_avgAline_filtered] = get_max_spectral_features(spectra_cube_avgAline_filtered, freq, 50);
        freq_spectra_avgAline_filtered = squeeze(sum(sum(spectra_cube_avgAline_filtered, 3), 2));
        figure; plot(freq, freq_spectra_avgAline_filtered);
        mean_freq_avgAline_filtered = sum(freq(:) .* freq_spectra_avgAline_filtered(:)) / sum(freq_spectra_avgAline_filtered(:));
        [~, max_idx] = max(freq_spectra_avgAline_filtered);
        peak_freq_avgAline_filtered = freq(max_idx);

        save(fullfile(SaveFolder, 'spectra_cube_avgAline_filtered.mat'), 'spectra_cube_avgAline_filtered', '-v7.3');
        save(fullfile(SaveFolder, 'max_val_map_avgAline_filtered.mat'), 'max_val_map_avgAline_filtered', '-v7.3');
        save(fullfile(SaveFolder, 'max_freq_map_avgAline_filtered.mat'), 'max_freq_map_avgAline_filtered', '-v7.3');
        save(fullfile(SaveFolder, 'freq_spectra_avgAline_filtered.mat'), 'freq_spectra_avgAline_filtered');
        save(fullfile(SaveFolder, 'mean_freq_avgAline_filtered.mat'), 'mean_freq_avgAline_filtered');
        save(fullfile(SaveFolder, 'peak_freq_avgAline_filtered.mat'), 'peak_freq_avgAline_filtered');
        
        fig = figure; 
        imagesc(max_freq_map_avgAline_filtered); 
        colorbar;
        clim([0 max_freq]);
        saveas(fig, fullfile(SaveFolder, 'max_freq_map_avgAline_filtered.fig'));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, 'max_freq_map_avgAline_filtered.png'));  % Save as PNG image

        fig = figure; 
        imagesc(max_val_map_avgAline_filtered); 
        colorbar;
        saveas(fig, fullfile(SaveFolder, 'max_val_map_avgAline_filtered.fig'));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, 'max_val_map_avgAline_filtered.png'));  % Save as PNG image

        fig = figure; 
        plot(freq, freq_spectra_avgAline_filtered);
        saveas(fig, fullfile(SaveFolder, 'freq_spectra_avgAline_filtered.fig'));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, 'freq_spectra_avgAline_filtered.png'));  % Save as PNG image

        clear spectra_cube_avgAline_filtered max_val_map_avgAline_filtered max_freq_map_avgAline_filtered freq_spectra_avgAline_filtered mean_freq_avgAline_filtered peak_freq_avgAline_filtered;
        % clear freq_map_even freq_map_odd val_map_odd val_map_even;
        %%
        info.c = 1500;                                                              info.fs = 200e6;
        info.N_ele_bf = 15;      % This is a tunable parameter, we can try 11, 15, 21 to see which one yields the best image
        info.ut_NA = 10e-2;
        info.xstep = 5e-6;                                                          info.ystep = 5e-6;
        info.pixel_z = info.c/info.fs;                                              bf_radius = (info.N_ele_bf-1)/2;
        Nzeropad = 810;
        info.t_focus = 4e-6;                                                        z_focus = info.t_focus*info.c;
        % RF_raw_all = zeros(size(raw_data1));                                        RF_Sum_all = zeros(size(raw_data1));
        W = zeros(data_total,info.N_ele_bf);

        for zi = 1:data_total
            r = floor(1e-1 + (info.ut_NA*info.pixel_z/info.xstep)*abs(zi + Nzeropad - info.t_focus*info.fs));
            r = min(r,bf_radius);
            W(zi,bf_radius+1-r:bf_radius+1+r) = 1;
        end
        info.SIR = W;

        f = [0 0.03 0.03 1]; m = [0 0.0 1 1]; DC_cancel = fir2(64,f,m); % HIGH PASS FILTER
        f = [0 0.7 0.7 1]; m = [1 0.9 0.1 0];   HF_cancel = fir2(64,f,m); % LOW PASS FILTER

        N = 21; %// Define size of Gaussian mask
        sigma = 2; %// Define sigma here
        %// Generate Gaussian mask
        ind = -floor(N/2) : floor(N/2);
        [X, Y] = meshgrid(ind, ind);
        h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
        h = h([2,4,6,8,11,14,16,18,20],:);
        h = [zeros(6,N);h;zeros(6,N)];
        h = h/sum(h(:));
        figure; imagesc(W);
        clearvars X Y f m
        N_repeat=1; idx_bscan = BLinesPerCScan/2 + 1;
        
        
        % Perform synthetic aperture beamforming on a 3D raw data volume
        [~, RF_Sum_all] = synthetic_aperture_beamforming(raw_data1_filtered, info, HF_cancel, DC_cancel, h, N, Nzeropad, N_repeat);
        clear raw_data1_filtered;
        %
        RF_DAS_env = RF_Sum_all / max(RF_Sum_all(:));
        %
        log_DAS_data = dynamic_ranging_log_compression(RF_DAS_env, imagingRecordTable.db_freq(row_idx));    % Used 40 for most cases

        signal_offset = 11;
        z_range = signal_offset:(data_total);
        projection_method = 'MIP';

        if strcmp(projection_method,'MIP') == 1
            im_mip_odd = squeeze(max((log_DAS_data(z_range,:,1:2:end)),[],1));
            im_mip_even = squeeze(max((log_DAS_data(z_range,:,2:2:end)),[],1));
        elseif strcmp(projection_method,'PIP') == 1
            im_mip_odd = squeeze(sum((log_DAS_data(z_range,:,1:2:end)),1));
            im_mip_even = squeeze(sum((log_DAS_data(z_range,:,2:2:end)),1));
        end

        im_mip_odd = mat2gray((im_mip_odd));
        im_mip_even = mat2gray((im_mip_even));
        figure; imagesc(im_mip_odd); title('mip odd');
        figure; imagesc(im_mip_even); title('mip even');
        save(fullfile(SaveFolder, 'im_mip_even.mat'), 'im_mip_even');
        save(fullfile(SaveFolder, 'im_mip_odd.mat'), 'im_mip_odd');


        %% PLOT MIP
        % Lukai code:
        l_wiener = 5;
        l_sharpen = 7;
        im_mip_odd_correct = (deconvlucy(wiener2(im_mip_odd,[l_wiener,l_wiener]),fspecial('gaussian',l_sharpen,l_sharpen),3));
        im_mip_even_correct = (deconvlucy(wiener2(im_mip_even,[l_wiener,l_wiener]),fspecial('gaussian',l_sharpen,l_sharpen),3));
        % Perform cropping and concatenation
        im_plot = [ ...
            im_mip_odd_correct(odd_top_crop+1:end-odd_bottom_crop,:); ...
            im_mip_even_correct(even_crop_start:even_crop_end,:) ...
            ];

        % im_plot = [im_mip_odd_correct(121:end-100,:); im_mip_even_correct(661:880,:)];
        im_plot = wiener2(im_plot,[5,5]);
        im_plot = mat2gray(imresize((im_plot),[ALinesPerBScan,BLinesPerCScan]));

        % Save the data
        save(fullfile(SaveFolder, 'im_plot.mat'), 'im_plot');
        % save(fullfile(SaveFolder, 'spectra_cube.mat'), 'spectra_cube', '-v7.3');


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

        % Crop out the ROI
        lower = imagingRecordTable{row_idx, 9}; upper = imagingRecordTable{row_idx, 10};
        if(isnan(lower) || isnan(upper))
            lower = 0; upper = 1;
        end
        lb = imagingRecordTable{row_idx, 17}; ub = imagingRecordTable{row_idx, 18};
        if(isnan(lb) || isnan(ub))
            lb = 1; ub = size(im_plot, 2);
        else
            lb = lb + 1;
        end

        im_plot_ROI = im_plot(:,lb:ub);
        im_plot_ROI = im_plot_ROI / max(im_plot_ROI(:));
        fig1 = figure;
        imshow(im_plot',[lower,upper],'border','tight');
        colormap('hot');
        set(gca,'xtick','','ytick','','fontweight','bold');
        saveas(fig1, fullfile(SaveFolder, 'im_plot_s.fig'));  % Save as MATLAB figure
        saveas(fig1, fullfile(SaveFolder, 'im_plot_s.png'));  % Save as PNG image

        fig2 = figure;
        imshow(im_plot_ROI',[lower,upper],'border','tight');
        colormap('hot');
        set(gca,'xtick','','ytick','','fontweight','bold');
        saveas(fig2, fullfile(SaveFolder, 'im_plot_ROI.fig'));  % Save as MATLAB figure
        saveas(fig2, fullfile(SaveFolder, 'im_plot_ROI.png'));  % Save as PNG image

        save(fullfile(SaveFolder, 'im_plot_ROI.mat'), 'im_plot_ROI');

        
        close all;
        clear im_plot im_plot_ROI RF_DAS_env RF_Sum_all log_DAS_data SaveFolder;
        clear max_idx lower upper lb ub;
    end
end
