% Modified from \\10.229.121.108\workspace\Lukai\ORPAM
% ANALYSIS\BSCAN_ANALYSIS.m to perform B-scan analysis in batch
% clc;
clear;
data_total = 256;
% BLinesPerCScan = 4000;  % Can be calculated from file size
ALinesPerBScan = 1000;
dB_PE=45;
% n_avg_alines: number of A-lines to average on each side (half-width)
n_avg_alines = 3;
% Define cropping parameters
odd_top_crop = 120;
odd_bottom_crop = 100;
even_crop_end = ALinesPerBScan - odd_top_crop;
even_crop_start = even_crop_end - odd_top_crop - odd_bottom_crop + 1;


% filepath = 'D:\ORPAM Sample';
filepath = '\\10.229.121.108\DataArchive\ORPAM_Data\EndometrialCancer';

% Cases = getFolderNames(filepath);
Cases = {'20250505_McCourt_normal_49'};

for i = 1:size(Cases, 2)

    aqdate = Cases{i};


    % Dir = 'D:\20250320_EC';
    Dir = fullfile(filepath, aqdate);

    S_temp = dir(fullfile(Dir,'*.bin'));
    N_bins = numel(S_temp);
    fprintf('>>>>>>>> FOUND %d SCANS\n', N_bins)

    %% SELECT SCAN AND LOAD DATA
    for sample_idx = 2:N_bins
        F = fullfile(Dir,S_temp(sample_idx).name);
        tokens = regexp(S_temp(sample_idx).name, '^(.*)\.bin$', 'tokens');
        filename = tokens{1}{1};
        SaveFolder = fullfile('\\10.229.121.108\Workspace\Lukai\EndometrialCancer',aqdate,filename);

        % if exist(SaveFolder)
        %     files = getFileNames(SaveFolder);
        %     if any(contains(files, 'max_freq_map_avgAline'))
        %         fprintf('>>>> Skipping %s (already processed)\n', filename);
        %         continue;
        %     end
        % else
        %     mkdir(SaveFolder);
        % end
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


        % fs = 20e6; % Sampling frequency (Hz)
        % [freq, spectra_cube] = spectrum_map(raw_data1, fs);
        % [max_val, max_freq] = get_max_spectral_features(spectra_cube, freq);
        % [max_val_map, max_freq_map] = get_max_spectral_features(spectra_cube, freq, 20);

        fs = 200e6; % Sampling frequency (Hz)
        [freq, spectra_cube] = spectrum_map(raw_data1_resized, fs);
        [max_val_map, max_freq_map] = get_max_spectral_features(spectra_cube, freq, 50);
        [freq_avgAline, spectra_cube_avgAline] = spectrum_map_avgAline(raw_data1_resized, fs, n_avg_alines);
        [max_val_map_avgAline, max_freq_map_avgAline] = get_max_spectral_features(spectra_cube_avgAline, freq_avgAline, 50);
        %%
        info.c = 1500;                                                              info.fs = 200e6;
        info.N_ele_bf = 15;      % This is a tunable parameter, we can try 11, 15, 21 to see which one yields the best image
        info.ut_NA = 10e-2;
        info.xstep = 5e-6;                                                          info.ystep = 5e-6;
        info.pixel_z = info.c/info.fs;                                              bf_radius = (info.N_ele_bf-1)/2;
        info.Nzeropad = 810;
        info.t_focus = 4e-6;                                                        z_focus = info.t_focus*info.c;
        % RF_raw_all = zeros(size(raw_data1));                                        RF_Sum_all = zeros(size(raw_data1));
        W = zeros(data_total,info.N_ele_bf);
        for zi = 1:data_total
            r = floor(1e-1 + (info.ut_NA*info.pixel_z/info.xstep)*abs(zi + info.Nzeropad - info.t_focus*info.fs));
            r = min(r,bf_radius);
            W(zi,bf_radius+1-r:bf_radius+1+r) = 1;
        end
        info.SIR = W;

        f = [0 0.03 0.03 1]; m = [0 0.0 1 1]; DC_cancel = fir2(64,f,m); % HIGH PASS FILTER
        f = [0 0.7 0.7 1]; m = [1 0.9 0.1 0];   HF_cancel = fir2(64,f,m); % LOW PASS FILTER
        info.HF_cancel = HF_cancel; info.DC_cancel = DC_cancel; 

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
        save(fullfile(SaveFolder, 'raw_data1_resized.mat'), 'raw_data1_resized', '-v7.3');
        clear raw_data1_resized raw_data1_crop_combined raw_data1_even_crop raw_data1_even raw_data1_odd_crop raw_data1_odd;
        % Perform synthetic aperture beamforming on a 3D raw data volume
        [RF_raw_all, RF_Sum_all] = synthetic_aperture_beamforming(raw_data1, info, h, N, N_repeat);
        %
        RF_DAS_env = RF_Sum_all / max(RF_Sum_all(:));
        %
        log_DAS_data = dynamic_ranging_log_compression(RF_DAS_env, 40);    % Used 40 for most cases

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
        % Yixiao code:
        % l_wiener = 5;
        % l_sharpen = 7;
        % im_mip_odd_correct = (deconvlucy(wiener2(im_mip_odd,[l_wiener,l_wiener]),fspecial('gaussian',l_sharpen,l_sharpen),3));
        % im_mip_even_correct = (deconvlucy(wiener2(im_mip_even,[l_wiener,l_wiener]),fspecial('gaussian',l_sharpen,l_sharpen),3));
        %
        % im_plot = (mat2gray(imresize((im_mip_odd),[ALinesPerBScan,BLinesPerCScan])));
        % figure
        % imshow((im_plot(121:end,:)),[0.0,0.4],'border','tight')
        % colormap('gray')
        % set(gca,'xtick','','ytick','','fontweight','bold')
        %
        % im_plot = mat2gray(imresize((im_mip_even_correct),[ALinesPerBScan,BLinesPerCScan]));
        % figure
        % imshow((im_plot(1:880,:)),[0.0,0.4],'border','tight')
        % colormap('gray')
        % set(gca,'xtick','','ytick','','fontweight','bold')
        %
        % im_plot = [2.25*im_mip_odd(121:end-100,:); im_mip_even(661:880,:)]; % im_plot = [im_mip_odd(121:end-100,:); im_mip_even(661:880,:)];
        % im_plot = wiener2(im_plot,[5,5]);
        % im_plot = mat2gray(imresize((im_plot),[ALinesPerBScan,BLinesPerCScan]));
        % figure
        % imshow(im_plot(:,:),[0.0,0.15],'border','tight')
        % colormap('gray')
        % set(gca,'xtick','','ytick','','fontweight','bold')

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
        save(fullfile(SaveFolder, 'max_val_map.mat'), 'max_val_map');
        save(fullfile(SaveFolder, 'max_freq_map.mat'), 'max_freq_map');
        save(fullfile(SaveFolder, 'spectra_cube.mat'), 'spectra_cube', '-v7.3');
        save(fullfile(SaveFolder, ['spectra_cube_avgAline_', num2str(n_avg_alines), '.mat']), 'spectra_cube_avgAline', '-v7.3');
        save(fullfile(SaveFolder, 'freq.mat'), 'freq');
        save(fullfile(SaveFolder, ['max_val_map_avgAline_',num2str(n_avg_alines),'.mat']), 'max_val_map_avgAline');
        save(fullfile(SaveFolder, ['max_freq_map_avgAline_',num2str(n_avg_alines),'.mat']), 'max_freq_map_avgAline');

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

        % Save the central frequency map
        fig = figure;
        imagesc(max_freq_map);
        colorbar;
        clim([5e5 6e6]);
        title('max freq map','fontsize',18);
        saveas(fig, fullfile(SaveFolder, 'max_freq_map.fig'));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, 'max_freq_map.png'));  % Save as PNG image
        % Save the central frequency value map
        fig = figure;
        imagesc(max_val_map);
        title('max val map','fontsize',18);
        saveas(fig, fullfile(SaveFolder, 'max_val_map.fig'));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, 'max_val_map.png'));  % Save as PNG image



        % Save the central frequency map
        fig = figure;
        imagesc(max_freq_map_avgAline);
        colorbar;
        clim([5e5 6e6]);
        title('max freq map avgAline','fontsize',18);
        saveas(fig, fullfile(SaveFolder, ['max_freq_map_avgAline_',num2str(n_avg_alines),'.fig']));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, ['max_freq_map_avgAline_',num2str(n_avg_alines),'.png']));  % Save as PNG image
        % Save the central frequency value map
        fig = figure;
        imagesc(max_val_map_avgAline);
        title('max val map avgAline','fontsize',18);
        saveas(fig, fullfile(SaveFolder, ['max_val_map_avgAline_',num2str(n_avg_alines),'.fig']));  % Save as MATLAB figure
        saveas(fig, fullfile(SaveFolder, ['max_val_map_avgAline_',num2str(n_avg_alines),'.png']));  % Save as PNG image
        close all;
        clear im_mip_odd im_mip_odd im_mip_even_correct im_mip_odd_correct raw_data1 raw_data1_even raw_data1_odd raw_data1_odd_crop raw_data1_even_crop RF_DAS_env RF_raw_all RF_Sum_all spectra_cube spectra_cube_avgAline log_DAS_data max_freq_map max_freq_map_avgAline max_val_map max_val_map_avgAline;
    end
end