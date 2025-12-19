% Apr 4, 2025: analysis of endometrium ORPAM frequency
clear;
close all;
RootFolder = '\\10.229.121.108\Workspace\Lukai\EndometrialCancer';
Case = getFolderNames(RootFolder);
min_freq = 8e5;
max_freq = 5e6;
freq_itv = 1e5;
% Threshold: keep pixels with intensity > -40 dB
threshold_dB = 40;
min_region_size = 20;  % Minimum number of pixels to keep
for i = 15:size(Case, 2)
    ParentFolder = fullfile(RootFolder, Case{i});
    Position = getFolderNames(ParentFolder);
    for j = 1:size(Position, 2)
        DataFolder = fullfile(ParentFolder, Position{j});
        % Load data
        load(fullfile(DataFolder, 'max_freq_map_avgAline_2.mat'));
        load(fullfile(DataFolder, 'im_plot.mat'));
        [freq_values, freq_mean, im_plot_filtered] = extract_filtered_freq(max_freq_map_avgAline, im_plot, min_freq, max_freq, threshold_dB, min_region_size);

        % Parameters
        min_freq = 8e5;
        max_freq = 5e6;
        freq_itv = 1e4;

        % Define bin edges
        edges = min_freq:freq_itv:max_freq;

        % Histogram count
        [counts, edges] = histcounts(freq_values, edges);

        % Find the index of the max bin
        [~, max_idx] = max(counts);

        % Get frequency range of that bin
        peak_bin_start = edges(max_idx);
        peak_bin_end   = edges(max_idx + 1);

        % Display
        fprintf('Peak frequency bin: %.2e Hz to %.2e Hz\n', peak_bin_start, peak_bin_end);
        fprintf('Number of pixels in this bin: %d\n', counts(max_idx));

        % Plot histogram
        figure;
        bar(edges(1:end-1) + freq_itv/2, counts, 'FaceColor', [0.2 0.2 0.8]); % Center the bins

        % Highlight peak bin
        hold on;
        bar(edges(max_idx) + freq_itv/2, counts(max_idx), 'FaceColor', 'r');
        hold off;

        xlabel('Frequency (Hz)');
        ylabel('Count');
        title('Frequency Distribution of max\_freq\_map');
        grid on;

    end

end
