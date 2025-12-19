function [a_dist, b_dist, c_dist] = fit_gaussian_cube(spectra_cube, freq)
    % FIT_GAUSSIAN_CUBE - Fits Gaussian to each pixel's spectrum in a 3D cube
    %
    % Inputs:
    %   spectra_cube - 3D array (frequency × y × x)
    %   freq - 1D frequency array (monotonically increasing)
    %
    % Outputs:
    %   a_dist - 2D array of amplitude parameters (y × x)
    %   b_dist - 2D array of center frequency parameters (y × x)
    %   c_dist - 2D array of width parameters (y × x)
    
    % Get dimensions
    [nFreq, nY, nX] = size(spectra_cube);
    
    % Preallocate output arrays
    a_dist = zeros(nY, nX);
    b_dist = zeros(nY, nX);
    c_dist = zeros(nY, nX);
    
    % Define Gaussian function (without offset)
    gaussian = @(x, a, b, c) a*exp(-((x-b)/c).^2);
    
    % Set options for fminsearch
    options = optimset('Display', 'off', 'MaxFunEvals', 1000, 'MaxIter', 1000);
    
    % Loop through each pixel
    for y = 1:nY
        for x = 1:nX
            % Extract spectrum for current pixel
            sp = squeeze(spectra_cube(:, y, x));
            
            % Get initial parameter guesses
            [maxVal, maxIdx] = max(sp);
            initialA = maxVal;             % Amplitude
            initialB = freq(maxIdx);       % Center frequency
            initialC = (max(freq)-min(freq))/6; % Width estimate
            
            % Perform the fit
            try
                params = fminsearch(@(p) sum((gaussian(freq, p(1), p(2), p(3)) - sp).^2), ...
                                  [initialA, initialB, initialC], options);
                
                % Store parameters
                a_dist(y, x) = params(1);
                b_dist(y, x) = params(2);
                c_dist(y, x) = params(3);
            catch
                % If fitting fails, use initial guesses
                % warning('Fit failed for pixel (%d, %d), using initial guesses', y, x);
                a_dist(y, x) = initialA;
                b_dist(y, x) = initialB;
                c_dist(y, x) = initialC;
            end
        end
    end
    
    % Display progress
    % fprintf('Gaussian fitting completed for %d x %d pixels\n', nY, nX);
end