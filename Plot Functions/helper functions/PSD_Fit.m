function [final_fit, aperiodic_params, Peak_freq, Peak_bw, peak_height, rmse] = PSD_Fit(data_vec, f, peak_distance_Hz, min_peak_height_STD, plotting)
% PSD_Fit.m Fits a power-law model to the power spectral density (PSD) data and performs data cleaning.
% Roni Naivelt | May 2024.
% INPUTS:
%   data_vec: 1D array, input data vector containing power spectral density values.
%   f: 1D array, frequency vector corresponding to data_vec.
%   freq_clean_range: 1x2 array, frequency range (in Hz) to focus on during peak cleaning. Default: [4 30].
%   clean_window: Scalar, window size (in Hz) around each peak for cleaning. Default: 2.5.
%   clean_thresh: Scalar, threshold multiplier for peak cleaning. Default: 1.25.
%   plotting: Logical, specifies whether to plot the results. Default: false.

% OUTPUTS:
  % C: Scalar, coefficient of the fitted power-law model.
  % A: Scalar, exponent of the fitted power-law model.
  % Curve_vec: 1D array, fitted power-law curve based on the model.
  % Delta_vec: 1D array, difference between original data and the fitted curve.

%% function handling
if nargin < 5
    plotting = false;  % Default value for plotting
end

if nargin < 4
    min_peak_height_STD = 1;
end

if nargin < 2
    peak_distance_Hz = 2.5;
end

%% Initializing
% Initialize arrays to store peak centers, widths, and Gaussian parameters
peak_centers = [];
peak_widths = [];
gaussian_params = {};
gaussian_region = {};

% Function for Gaussian fitting
gauss_eq = @(a, b, c, x) a * exp(-((x - b).^2) / (2 * c^2));
% function for log fit
% powerLawModel = fittype('(C / (f^A)) + B', 'independent', 'f', 'coefficients', {'C', 'A', 'B'});
powerLawModel = fittype('(C / (f^A))', 'independent', 'f', 'coefficients', {'C', 'A'});

% adjusting frequency and data vectors
f = double(f);
data_vec = double(data_vec);

if find(f==0)                % check for f=0
    zero_indx = find(f==0);
    f = f(zero_indx+1:end);
    data_vec = data_vec(zero_indx+1:end);
end

if size(f,1) ~= 1            % check frequency input size
    f = f';
end

if size(data_vec,1) ~= 1     % check data input size
    data_vec = data_vec';
end

%% Step 1: fit and remove aperiodic fit

% [fit_result,gof] = fit(f',data_vec',powerLawModel,...
%                        'StartPoint',[1, 1, 0],... %  StartPoint = Initial guess for [C, A, B]
%                        'Lower',[0,-inf, -inf], ...      % Lower bounds for [C, A, B]
%                        'Upper',[max(data_vec), inf, inf]);    % Upper bounds for [C, A, B]
[fit_result,gof] = fit(f',data_vec',powerLawModel,...
                       'StartPoint',[1, 1],... %  StartPoint = Initial guess for [C, A]
                       'Lower',[0,-inf], ...      % Lower bounds for [C, A]
                       'Upper',[max(data_vec), inf]);    % Upper bounds for [C, A]

% Extract fitted parameters
C = fit_result.C;
A = fit_result.A;
% B = fit_result.B;
% aperiodic_fit = (C ./ (f.^A)) + B;
aperiodic_fit = (C ./ (f.^A));
periodic_vec = data_vec - aperiodic_fit;

% % Trim start if needed (sometimes the start is lower):
% trim_thresh = mean(periodic_vec) - std(periodic_vec);
% trim_indx = 1;
% while periodic_vec(trim_indx) < trim_thresh
%     trim_indx = trim_indx + 1;
% end
% f = f(trim_indx:end);
% periodic_vec = periodic_vec(trim_indx:end);
% aperiodic_fit = aperiodic_fit(trim_indx:end);
% data_vec = data_vec(trim_indx:end);

% figure()
% plot(f,data_vec)
% hold on
% plot(f,aperiodic_fit)
% plot(f,periodic_vec)
% legend('original', 'aperiodic fit', 'periodic')

%% Step 2: fit a gaussian, find peaks, width and height

% Step 2.1: Find peaks higher than 1 standard deviation
std_dev = mean(periodic_vec) + min_peak_height_STD*std(periodic_vec);
[pks, locs] = findpeaks(periodic_vec, 'MinPeakHeight', std_dev, 'MinPeakDistance',peak_distance_Hz/(mean(diff(f))));

% Step 2.2: Fit Gaussian and subtract from original signal
final_vec = periodic_vec;
gauss_indx = 1; % Index to manage storage in gaussian_params and other arrays

for i = 1:length(pks)
    local_mean  = movmean(final_vec, 5/mean(diff(f)));
    local_std = movstd(final_vec,  5/mean(diff(f)));

    % Define dynamic region around the peak
    peak_center = locs(i);
    left_bound = peak_center;
    right_bound = peak_center;
    
    % Find left bound
    while left_bound > 1 &&...
        final_vec(left_bound) > local_mean(peak_center) % OPT: +0.2*local_std(peak_center)
        if gauss_indx > 1
            if left_bound == gaussian_region{gauss_indx-1}(end)+1
                break;
            end
        end
        left_bound = left_bound - 1;
    end
    
    % Find right bound
    while right_bound < length(final_vec) &&...
            final_vec(right_bound) > local_mean(peak_center) % OPT: +0.2*local_std(peak_center)
        right_bound = right_bound + 1;
    end
    
    % Define the region
    region = left_bound:right_bound;
    if length(region) <3
        continue
    end

    % Fit Gaussian to the region
    x_data = region';
    y_data = final_vec(region)';
    
    
    % Initial guess for parameters [amplitude, mean, stddev]
    init_guess = [pks(i), peak_center, 1];  
    % Fit the Gaussian with bounds
    gauss_fit = fit(x_data, y_data, gauss_eq,...
                    'StartPoint', init_guess,...
                    'Lower', [0, left_bound, 0],...  
                    'Upper', [max([0, max(final_vec(region))]), right_bound, length(final_vec)]); 
                     % Lower bounds: Amplitude >= 0, Mean >= left_bound, Stddev >= 0 
                     % Upper bounds: Amplitude <= max(data), Mean <= right_bound, Stddev <= length(data)

    % Save peak center, width (stddev), and Gaussian parameters
    
    peak_centers = [peak_centers; gauss_fit.b];
    peak_widths = [peak_widths; gauss_fit.c];
    gaussian_params{gauss_indx} = coeffvalues(gauss_fit);
    gaussian_region{gauss_indx} = region;
    gauss_indx = gauss_indx + 1; % Increment index for next valid Gaussian

    % Subtract Gaussian from the signal
    fitted_gaussian = gauss_eq(gauss_fit.a, gauss_fit.b, gauss_fit.c, x_data);
    final_vec(region) = final_vec(region) - fitted_gaussian';
end


% % Plot the original graph, the final graph, and all Gaussians
% figure;
% hold on;
% 
% % Plot original vector
% plot(periodic_vec, 'Color', [1 0.5 0], 'DisplayName', 'Original Signal');
% 
% % plot noise threshold
% yline(mean(periodic_vec) + std(periodic_vec), 'LineWidth', 1.5, 'DisplayName', 'Noise Threshold');
% 
% % Plot final vector
% plot(final_vec, 'Color', [0 0 1], 'LineWidth', 2, 'DisplayName', 'Final Signal after Subtraction');
% 
% % Plot each fitted Gaussian
% for i = 1:length(gaussian_params)
%     gauss_fit = gaussian_params{i};
%     x_data = gaussian_region{i};
%     fitted_gaussian = gauss_eq(gauss_fit(1), gauss_fit(2), gauss_fit(3), x_data);
%     plot(x_data, fitted_gaussian, '--', 'DisplayName', sprintf('Gaussian Fit %d', i));
% end
% 
% legend;
% title('Original Signal, Final Signal, and Gaussian Fits');
% xlabel('Sample Index');
% ylabel('Amplitude');
% grid on;
% hold off;


%% Step 3: Remove peaks from original vector and fit again
% Subtract Gaussian from the signal
data_vec_flattened = data_vec;
for i = 1:length(gaussian_params)
    gauss_fit = gaussian_params{i};
    region = gaussian_region{i};
    fitted_gaussian = gauss_eq(gauss_fit(1), gauss_fit(2), gauss_fit(3), region);
    data_vec_flattened(region) = data_vec_flattened(region) - fitted_gaussian;
end

% make the fit again on the flattened vec
% [fit_result,gof] = fit(f',data_vec_flattened',powerLawModel,...
%                        'StartPoint',[1, 1, 0],... %  StartPoint = Initial guess for [C, A, B]
%                        'Lower',[0,-inf, -inf], ...      % Lower bounds for [C, A, B]
%                        'Upper',[max(data_vec), inf, inf]);    % Upper bounds for [C, A, B]
[fit_result,gof] = fit(f',data_vec_flattened',powerLawModel,...
                       'StartPoint',[1, 1],... %  StartPoint = Initial guess for [C, A]
                       'Lower',[0,-inf], ...      % Lower bounds for [C, A]
                       'Upper',[max(data_vec), inf]);    % Upper bounds for [C, A]

% Extract fitted parameters
C = fit_result.C; % C = offset coeffitiant of the logarithmic fit (C/f^A)
A = fit_result.A; % A = exponent of the logirithmic fit (C/f^A)
% B = fit_result.B;



% flattened_aperiodic_fit = (C ./ (f.^A)) + B;
flattened_aperiodic_fit = (C ./ (f.^A));
flatenned_periodic_vec = data_vec_flattened - flattened_aperiodic_fit;

final_fit = flattened_aperiodic_fit;
for i = 1:length(gaussian_params)
    gauss_fit = gaussian_params{i};
    region = gaussian_region{i};
    fitted_gaussian = gauss_eq(gauss_fit(1), gauss_fit(2), gauss_fit(3), region);
    final_fit(region) = final_fit(region) + fitted_gaussian;
end

% goodness of fit test (rmse)
rmse = sqrt(mean((data_vec - final_fit).^2));


if plotting
    figure()
    plot(f, data_vec, 'DisplayName', 'original')
    hold on
    plot(f,data_vec_flattened, 'DisplayName', 'flattened')
    plot(f,flattened_aperiodic_fit, 'DisplayName', 'aperiodic')
    plot(f,final_fit, 'LineWidth', 2, 'DisplayName', 'final fit')
    % Plot each fitted Gaussian
    for i = 1:length(gaussian_params)
        gauss_fit = gaussian_params{i};
        x_data = gaussian_region{i};
        fitted_gaussian = flattened_aperiodic_fit(x_data) + gauss_eq(gauss_fit(1), gauss_fit(2), gauss_fit(3), x_data);
        plot(f(x_data), fitted_gaussian, '--', 'DisplayName', sprintf('Gaussian Fit %d', i));
    end
     xlabel('Frequency [Hz]')
    ylabel('PSD [uV^2/Hz]')
    title('Fit function result')
    legend()
    set(gcf,'Color','w')
end

%% step 4: save parameters

% save aperiodic parameters
aperiodic_params.A = A;
% aperiodic_params.B = B;
aperiodic_params.C = C;

% gaussians
Peak_freq = [];
for i = 1:length(peak_centers)
    Peak_freq = [Peak_freq f(round(peak_centers(i)))];
end

Peak_bw = [];
for i = 1:length(gaussian_region)
    region = gaussian_region{i};
    f_range = f(region(end))-f(region(1));
    Peak_bw = [Peak_bw f_range];
end

peak_height = [];
for i = 1:length(peak_centers)
    peak_height = [peak_height final_fit(round(peak_centers(i)))- flattened_aperiodic_fit(round(peak_centers(i)))];
end