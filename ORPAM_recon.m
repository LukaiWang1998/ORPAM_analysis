function [im_plot, RF] = ORPAM_recon(raw_data1, info, dB, BLinesPerCScan)
N_repeat=1;
data_total = 256;
ALinesPerBScan = 1000;
% Define cropping parameters
odd_top_crop = 120;
odd_bottom_crop = 100;
even_crop_end = ALinesPerBScan - odd_top_crop;
even_crop_start = even_crop_end - odd_top_crop - odd_bottom_crop + 1;
% Define Gaussian deconvolution kernel
N = 21; %// Define size of Gaussian mask
sigma = 2; %// Define sigma here
%// Generate Gaussian mask
ind = -floor(N/2) : floor(N/2);
[X, Y] = meshgrid(ind, ind);
h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
h = h([2,4,6,8,11,14,16,18,20],:);
h = [zeros(6,N);h;zeros(6,N)];
h = h/sum(h(:));
% figure; imagesc(info.SIR);
clearvars X Y f m

% Perform synthetic aperture beamforming on a 3D raw data volume
[RF_raw_all, RF_Sum_all] = synthetic_aperture_beamforming(raw_data1, info, h, N, N_repeat);
%
RF_DAS_env = RF_Sum_all / max(RF_Sum_all(:));
%
log_DAS_data = dynamic_ranging_log_compression(RF_DAS_env, dB);    % Used 40 for most cases

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

%% CALCULATE MIP
l_wiener = 5;
l_sharpen = 7;
im_mip_odd_correct = (deconvlucy(wiener2(im_mip_odd,[l_wiener,l_wiener]),fspecial('gaussian',l_sharpen,l_sharpen),3));
im_mip_even_correct = (deconvlucy(wiener2(im_mip_even,[l_wiener,l_wiener]),fspecial('gaussian',l_sharpen,l_sharpen),3));

% Perform cropping and concatenation
im_plot = [ ...
    im_mip_odd_correct(odd_top_crop+1:end-odd_bottom_crop,:); ...
    im_mip_even_correct(even_crop_start:even_crop_end,:) ...
    ];

im_plot = wiener2(im_plot,[5,5]);
im_plot = mat2gray(imresize((im_plot),[ALinesPerBScan,BLinesPerCScan]));

RF.RF_raw_all = RF_raw_all; 
RF.RF_Sum_all = RF_Sum_all; 
RF.RF_DAS_env = RF_DAS_env;