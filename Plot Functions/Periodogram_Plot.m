function Periodogram_Plot(DBS_data, patient_name, hemi, periodogram_days, varargin)
% Periodogram_Plot - Plots PSD periodogram for LFP data with frequency/contact info.
%
% Inputs:
%   DBS_data            - Nested data structure containing DBS data
%   patient_name        - String patient identifier
%   hemi                - Hemisphere: 'Left' or 'Right'
%   periodogram_days  : - (Optional) array with [start_day, end_day]

% Name-Value Parameters:
%   'colors'             - RGB vector or color string for patient (default turbo colormap)
%

if nargin < 4
    periodogram_days = [];
end

% Parse optional parameters
p = inputParser;
addParameter(p, 'colors', []);
parse(p, varargin{:});
colors = p.Results.colors;

% Determine color
if isempty(colors)
    colors = turbo(1); % Default colorful
end
Pat_Color = colors(1,:);

data_out = OpenDataHelper(DBS_data, patient_name, hemi, periodogram_days);
key = [patient_name '_' hemi];
LFP_vec_trimmed = data_out.(key).LFP_vec_trimmed;
Datetime_vec_trimmed = data_out.(key).Datetime_vec_trimmed;
sensing_vec_trimmed = data_out.(key).sensing_vec_trimmed;
contacts_vec_trimmed = data_out.(key).contacts_vec_trimmed;


% Parse frequency and contact from sensing_vec_trimmed and contacts_vec_trimmed
label_str = Parse_Sensing_Contact_Info(sensing_vec_trimmed, contacts_vec_trimmed);



% Calculate sample frequency & time resolution
sample_interval = abs(hours(Datetime_vec_trimmed(2) - Datetime_vec_trimmed(1)));
sample_freq = 1 / sample_interval;
time_res = 0.16; % hours

% Define period and frequency axes for periodogram
min_period = max([2 * sample_interval, time_res]);
max_period = 80; % hours
time_ax_vec = min_period:time_res:max_period;
freq_ax_vec = 1 ./ time_ax_vec;

% Calculate recording duration in days
rec_days = round(days(Datetime_vec_trimmed(end) - Datetime_vec_trimmed(1)));

% Define window size and overlap for pwelch (in bins)
win_size_days = min(ceil(0.6 * rec_days - 1), 7);
win_size_bins = win_size_days * 6 * 24; % Assuming 10-min bins => 6 per hour, 24 hours/day
win_overlap_days = max(ceil(win_size_days / 2), 1);
win_overlap_bins = win_overlap_days * 6 * 24;

% Compute PSD estimate using pwelch
[psd_estimate, ~] = pwelch(LFP_vec_trimmed, win_size_bins, win_overlap_bins, freq_ax_vec, sample_freq);

% Plot
figure('Color', 'w');
plot(time_ax_vec, psd_estimate, 'LineWidth', 2.5, 'Color', Pat_Color, ...
    'DisplayName', sprintf('%s %s PSD at %s', label_str));
hold on;

grid off;
box off;
set(gca, 'LineWidth', 2, 'FontSize', 10);
xlabel('Cycle Duration [Hours]');
ylabel('PSD Estimate [μV^2/Hz]');
title(sprintf('Periodogram - %s (%s Hemisphere) - %s', patient_name, hemi, label_str), 'Interpreter', 'none');
xlim([0 40]);

legend('show', 'Interpreter', 'none');

end


function label_str = Parse_Sensing_Contact_Info(sensing_vec, contacts_vec)
% Parse_Sensing_Contact_Info - Creates a label string based on sensing frequencies and contacts.
%
% Inputs:
%   sensing_vec  - Vector of sensing frequencies (Hz)
%   contacts_vec - Vector of contact pairs (assumes pairs are stored as e.g. [1 3] or [2 2])
%
% Output:
%   label_str - Combined string describing sensing frequency and contacts

    sensing_unique = unique(sensing_vec);
    contacts_unique = unique(contacts_vec, 'rows');  % assumes Nx2 array

    if numel(sensing_unique) == 1 && size(contacts_unique, 1) == 1
        % Case 1: Single sensing freq and contact pair
        contact_pair = contacts_unique(1,:);
        label_str = sprintf('%2.2f Hz - Contacts %s', sensing_unique, contact_pair);

    elseif numel(sensing_unique) > 1 && size(contacts_unique, 1) == 1
        % Multiple freqs, single contact pair
        contact_pair = contacts_unique(1,:);
        label_str = sprintf('%2.2f-%.2f Hz - Contact %s', ...
            min(sensing_unique), max(sensing_unique), contact_pair);

    elseif numel(sensing_unique) == 1 && size(contacts_unique, 1) > 1
        % Single freq, multiple contact pairs
        label_str = sprintf('%2.2f Hz', sensing_unique);

    else
        % Multiple freqs and multiple contacts
        label_str = sprintf('%2.2f-%2.2f Hz', min(sensing_unique), max(sensing_unique));
    end

end
