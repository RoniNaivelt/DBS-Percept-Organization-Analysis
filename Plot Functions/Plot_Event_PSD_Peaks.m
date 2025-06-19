function Plot_Event_PSD_Peaks(DBS_data, patient_name, hemi, sensing_range, day_range)
% Plot_Event_PSD_Bubble
% This function generates a bubble chart of event-related PSD peak characteristics.
%
% Inputs:
%   DBS_data            - Nested data structure containing DBS data
%   patient_name        - String or cell array of strings (must match DBS_data entries)
%   hemi                - (Optional) 'Left', 'Right', or 'both' (default)
%   sensing_range       - (Optional) array indicating a band of sensing
%                         frequencies such as [Min_frequency Max_frequency]
%   'day_range'         - [start end], or cell array of ranges per patient



% Handle default inputs
if nargin < 3 || isempty(hemi)
    hemi = 'both';
end
if nargin < 4 || isempty(sensing_range)
    sensing_range = [0 inf];
end
if nargin < 5 || isempty(day_range)
    day_range = [];
end

% Convert to cell arrays if needed
if ~iscell(patient_name)
    patient_name = {patient_name};
end
n_patients = numel(patient_name);



for i_pat = 1:n_patients
    name = patient_name{i_pat};
    % Choose hemispheres
    if strcmpi(hemi, 'both')
        hemis = {'Left', 'Right'};
    else
        hemis = {hemi};
    end

    for i_hemi = 1:length(hemis)
        h_name = hemis{i_hemi};

        % Resolve day range
        if isempty(day_range)
            current_range = [1 numel(DBS_data.(name).TrendLogs.([h_name '_Hemi']).LFP_table.Properties.VariableNames)];
        elseif iscell(day_range)
            if i_pat > numel(day_range)
                error('Not enough entries in DayRange for patient %d', i_pat);
            end
            current_range = day_range{i_pat};
        elseif isnumeric(day_range) && numel(day_range) == 2
            current_range = day_range;
        else
            error('Invalid format for DayRange.');
        end

        % Use OpenDataHelper to extract band_LFP_mat
        try
            data_out = OpenDataHelper(DBS_data, name, h_name, current_range, sensing_range);
            key = [name '_' h_name];
            % For Events
            Events = data_out.(key).Events;
            Events_Days_Since_Surge = data_out.(key).Event_Summary.Event_dates_since_surg;
        catch
            warning('Failed to extract LFP data for %s (%s hemisphere). Skipping.', name, h_name);
            continue;
        end

        % Initialize storage
        Event_dates = datetime.empty;
        Event_names = {};
        Event_FFT = {};
        Event_fit = {};
        Event_Peak_freq = {};
        Event_Peak_bw = {};
        Event_peak_height = {};
        aperiodic_params = {};


        % Extract and fit PSD data from events
        for i = 1:size(Events, 1)
            snapshot = Events.([h_name '_Snapshot'])(i);
            if Events.Snapshot(i) == 1 && ~isempty(snapshot.FFT_Val)
                Event_dates(end+1) = Events.DateTime(i);
                Event_names{end+1} = Events.EventName{i};

                FFT = snapshot.FFT_Val;
                Freq = snapshot.FFT_Freq;

                % Fit PSD
                [fit_vals, aperiodic, peak_freqs, peak_bw, peak_heights, ~] = ...
                    PSD_Fit(FFT(6:end), Freq(6:end), 1.5, 1, false);

                Event_FFT{end+1} = fit_vals;
                Event_fit{end+1} = fit_vals;
                aperiodic_params{end+1} = aperiodic;
                Event_Peak_freq{end+1} = peak_freqs;
                Event_Peak_bw{end+1} = peak_bw;
                Event_peak_height{end+1} = peak_heights;
            end
        end

        % Prepare plotting data
        x_data = [];
        y_data = [];
        bubble_colors = [];
        bubble_sizes = [];
        aperiodic_exponent = [];
        aperiodic_offset = [];

        for i = 1:length(Events_Days_Since_Surge)
            num_peaks = length(Event_Peak_freq{i});
            x_data = [x_data; repmat(Events_Days_Since_Surge(i), num_peaks, 1)];
            y_data = [y_data; Event_Peak_freq{i}'];
            bubble_colors = [bubble_colors; Event_peak_height{i}'];
            bubble_sizes = [bubble_sizes; Event_Peak_bw{i}'];
            aperiodic_exponent = [aperiodic_exponent; repmat(aperiodic_params{i}.A, num_peaks, 1)];
            aperiodic_offset = [aperiodic_offset; repmat(aperiodic_params{i}.C, num_peaks, 1)];
        end

        % Normalize colors and scale bubble sizes
        bubble_colors = normalize(bubble_colors);
        min_bw = min(bubble_sizes);
        max_bw = max(bubble_sizes);
        scaled_bubble_sizes = 350 * (bubble_sizes - min_bw) / (max_bw - min_bw) + 1;

        % Plot
        figure()
        scatter(x_data, y_data, scaled_bubble_sizes, bubble_colors, 'filled');
        colormap(jet);
        c = colorbar;
        clim([0 1]);
        ylabel(c, 'Peak Amplitude', 'FontSize', 12);
        ylim([1 60]);

        xlabel('Days Since Surgery');
        ylabel('Frequency [Hz]');
        title({[h_name ' Hemisphere Peaks of ' name], ...
            'Color - Peak Amplitude ; Radius - Peak Bandwidth'}, ...
            'FontSize', 14, 'Interpreter', 'none');

        set(gca, 'FontSize', 12, 'XMinorGrid', 'on');
        set(gcf, 'Color', 'w');
        grid minor;
        box off;
    end
end
end