function [freq_values_1d, freq_mean, im_plot_filtered] = extract_filtered_freq( ...
    max_freq_map_avgAline, im_plot, min_freq, max_freq, threshold_dB, min_region_size)

    % Convert im_plot to dB (avoid log(0) by adding a small value)
    epsilon = 1e-12;
    im_plot_dB = 20 * log10(im_plot + epsilon);

    % Binary mask from intensity thresholding
    binary_map = im_plot_dB > -threshold_dB;

    % Frequency range mask
    freq_map = (max_freq_map_avgAline > min_freq) & (max_freq_map_avgAline < max_freq);

    % Filter out small connected regions in freq_map
    CC_freq = bwconncomp(freq_map);
    stats_freq = regionprops(CC_freq, 'Area');
    large_freq_idx = find([stats_freq.Area] >= min_region_size);
    filtered_freq_map = ismember(labelmatrix(CC_freq), large_freq_idx);

    % Filter out small connected regions in binary_map
    CC_bin = bwconncomp(binary_map);
    stats_bin = regionprops(CC_bin, 'Area');
    large_bin_idx = find([stats_bin.Area] >= min_region_size);
    filtered_binary_map = ismember(labelmatrix(CC_bin), large_bin_idx);

    % Combined map: only keep pixels that satisfy both masks
    filtered_combined_map = filtered_freq_map & filtered_binary_map;

    % Extract valid frequency and intensity values
    freq_values_1d = max_freq_map_avgAline(filtered_combined_map);
    freq_mean = mean(freq_values_1d);
    
    % Set invalid pixels in im_plot to 0
    im_plot_filtered = zeros(size(im_plot));
    im_plot_filtered(filtered_combined_map) = im_plot(filtered_combined_map);
end

