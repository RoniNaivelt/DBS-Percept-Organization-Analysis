function Diurnal_Plots = Diurnal_LFP_Plot(DBS_data, patient_name, hemi, sensing_range, day_range, varargin)
% Diurnal_LFP_Plot - Plots diurnal LFP data for one or more patients and hemispheres.
%
% Inputs:
%   DBS_data            - Nested data structure containing DBS data
%   patient_name        - String or cell array of strings (must match DBS_data entries)
%   hemi                - (Optional) 'Left', 'Right', or 'both' (default)
%   sensing_range       - (Optional) array indicating a band of sensing
%                         frequencies such as [Min_frequency Max_frequency]     
%   'day_range'         - [start end], or cell array of ranges per patient

% Name-Value Parameters:
%   'PatientIndex'      - Index for specific patient (if looping manually)
%   'PlotMean'          - true/false
%   'DiurnalPlots'      - Existing struct to continue adding plots
%   'colors'             - RGB vector or color string

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


% Parse optional name-value pairs
p = inputParser;
addParameter(p, 'DayRange', []);
addParameter(p, 'PatientIndex', 1);
addParameter(p, 'PlotMean', false);
addParameter(p, 'DiurnalPlots', struct());
addParameter(p, 'colors', []);
parse(p, varargin{:});

plot_mean = p.Results.PlotMean;
Diurnal_Plots = p.Results.DiurnalPlots;
colors = p.Results.colors;

% Convert to cell arrays if needed
if ~iscell(patient_name)
    patient_name = {patient_name};
end
n_patients = numel(patient_name);

% If no color input, generate a colorful grid of colors for all patients
if isempty(colors)
    colors = turbo(n_patients);
end

figure();

for i_pat = 1:n_patients
    name = patient_name{i_pat};
    base_color = colors(i_pat, :);  % base color for this patient
    % Choose hemispheres
    if strcmpi(hemi, 'both')
        hemis = {'Left', 'Right'};
    else
        hemis = {hemi};
    end

    for i_hemi = 1:length(hemis)
        h_name = hemis{i_hemi};
        if strcmpi(h_name, 'Left')
            % Slightly brighter or lighter for left
            color = min(base_color + 0.3, 1);
        else
            % Slightly darker for right
            color = max(base_color - 0.3, 0);
        end

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
            band_LFP_mat = data_out.([name '_' h_name]).band_LFP_mat;
        catch
            warning('Failed to extract LFP data for %s (%s hemisphere). Skipping.', name, h_name);
            continue;
        end

        if isempty(band_LFP_mat)
            warning('No LFP data found for %s (%s). Skipping.', name, h_name);
            continue;
        end



        % Compute mean LFP across selected days
        LFP_mean = mean(band_LFP_mat, 2, 'omitnan');
        key_name = [name '_' h_name];
        Diurnal_Plots.(key_name) = LFP_mean;

        % Create polar plot
        polarplot(deg2rad(linspace(1, 360, length(LFP_mean))), LFP_mean, ...
            'LineWidth', 2, 'Color', color, 'DisplayName', key_name);
        hold on;

        % Axis formatting
        ax = gca;
        ax.ThetaTick = 0:15:345;
        ax.ThetaTickLabel = arrayfun(@(h) sprintf('%02d:00', h), 0:23, 'UniformOutput', false);
        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir = 'clockwise';
        set(ax, 'FontSize', 10);
        title('Diurnal LFP Power', 'FontSize', 12, 'Interpreter', 'none');
        set(gcf, 'Color', 'w');
        rlim([0 1]);
    end
end

% Optionally plot mean curve across all patients for that hemisphere
if plot_mean
    hemi_fields = fieldnames(Diurnal_Plots);
    hemi_fields = hemi_fields(contains(hemi_fields, ['_' h_name]));

    all_curves = [];
    for f = 1:length(hemi_fields)
        all_curves = [all_curves, Diurnal_Plots.(hemi_fields{f})];
    end

    mean_curve = mean(all_curves, 2, 'omitnan');
    polarplot(deg2rad(linspace(1, 360, length(mean_curve))), mean_curve, ...
        '--', 'Color', 'k', 'LineWidth', 4, 'DisplayName', 'Mean');
end

legend('show', 'Interpreter', 'none');

end
