function raw_data_aligned = align_bscans(raw_data1, Noffset)
% Align alternating B-scans in a 3D OR-PAM volume
% Inputs:
%   raw_data1: 3D array [depth × A-lines × B-scans]
%   Noffset: horizontal shift amount (in pixels)
% Output:
%   raw_data_aligned: aligned 3D array

[~, ~, num_bscans] = size(raw_data1);
raw_data_aligned = raw_data1;

f_align = waitbar(0, 'Aligning B scans ...');

for i_scan = 1:num_bscans
    bscan = squeeze(raw_data1(:,:,i_scan));
    if mod(i_scan, 2) == 1
        bscan = imtranslate(bscan, [-round(Noffset/2), 0]);  % Odd: shift left
    else
        bscan = imtranslate(fliplr(bscan), [round(Noffset/2), 0]);  % Even: flip and shift right
    end
    raw_data_aligned(:,:,i_scan) = bscan;

    % Progress bar update
    perc = i_scan / num_bscans;
    mesg = sprintf('Aligning B scans ... %0.3f', perc);
    waitbar(perc, f_align, mesg);
end

close(f_align);
disp('>>>>>>>> DATA ALIGNMENT COMPLETE');
end
