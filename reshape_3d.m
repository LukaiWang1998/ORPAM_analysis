function reshaped_cube = reshape_3d(original_cube, ...
    odd_top_crop, odd_bottom_crop, BLinesPerCScan)
    
    ALinesPerBScan = size(original_cube, 2);
    ALinesDimension = size(original_cube, 1);
    even_crop_end = ALinesPerBScan - odd_top_crop;
    even_crop_start = even_crop_end - odd_top_crop - odd_bottom_crop + 1;
    % Split into odd and even B-scans
    original_cube_odd = original_cube(:,:,1:2:end);
    original_cube_even = original_cube(:,:,2:2:end);

    % Crop the second dimension (A-lines)
    original_cube_odd_crop = original_cube_odd(:, odd_top_crop+1:end-odd_bottom_crop, :);
    original_cube_even_crop = original_cube_even(:, even_crop_start:even_crop_end, :);

    % Concatenate along the second (A-line) dimension
    reshaped_cube = cat(2, original_cube_odd_crop, original_cube_even_crop);

    % Resize to target size
    reshaped_cube = imresize3(reshaped_cube, [ALinesDimension, ALinesPerBScan, BLinesPerCScan]);

end
