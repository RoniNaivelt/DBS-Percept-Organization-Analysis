function Histogram_Plot(DBS_data, patient_name, sensing_range, day_range, varargin)
% Histogram_Plot - Plots histograms comparing LFP activity:
% 1. Left vs Right
% 2. Left: Day vs Night
% 3. Right: Day vs Night

% Default arguments
if nargin < 3 || isempty(sensing_range)
    sensing_range = [0 inf];
end
if nargin < 4 || isempty(day_range)
    day_range = [];
end

% Prepare data containers
data_struct = struct();
hemis = {'Left', 'Right'};
figure('Color','w', 'Position', [100 100 1200 400]);

% Collect band_LFP data for each hemisphere
for i_hemi = 1:length(hemis)
    h_name = hemis{i_hemi};

    % Determine current_range based on day_range input format
    if isempty(day_range)
        all_vars = DBS_data.(patient_name).TrendLogs.([h_name '_Hemi']).LFP_table.Properties.VariableNames;
        current_range = [1, numel(all_vars)];
    elseif iscell(day_range)
        current_range = day_range{i_hemi};
    elseif isnumeric(day_range) && numel(day_range) == 2
        current_range = day_range;
    else
        error('Invalid format for day_range.');
    end

    % Try to load LFP data
    try
        data_out = OpenDataHelper(DBS_data, patient_name, h_name, current_range, sensing_range);
        key = [patient_name '_' h_name];

        % Extract and clean
        band_all = data_out.(key).band_LFP_mat;
        band_day = data_out.(key).band_LFP_mat_day;
        band_night = data_out.(key).band_LFP_mat_night;

        data_struct.(h_name).all = band_all(~isnan(band_all));
        data_struct.(h_name).day = band_day(~isnan(band_day));
        data_struct.(h_name).night = band_night(~isnan(band_night));

    catch ME
        warning('Failed to extract LFP data for %s (%s hemisphere): %s', patient_name, h_name, ME.message);
        continue; 
    end
end

% ----------- PLOTTING ------------

% --- 1. Left vs Right ---
subplot(2,2,[1 2]);
if isfield(data_struct, 'Left') && isfield(data_struct, 'Right')
    histogram(data_struct.Left.all, 'Normalization', 'probability', 'FaceAlpha', 0.6, 'DisplayName', 'Left');
    hold on;
    histogram(data_struct.Right.all, 'Normalization', 'probability', 'FaceAlpha', 0.6, 'DisplayName', 'Right');
    title('LFP: Left vs Right');
    xlabel('LFP amplitude');
    ylabel('Probability');
    legend;
else
    title('LFP: Left vs Right - Incomplete Data');
end

% --- 2. Left: Day vs Night ---
subplot(2,2,3);
if isfield(data_struct, 'Left')
    histogram(data_struct.Left.day, 'Normalization', 'probability', 'FaceAlpha', 0.6, 'DisplayName', 'Day');
    hold on;
    histogram(data_struct.Left.night, 'Normalization', 'probability', 'FaceAlpha', 0.6, 'DisplayName', 'Night');
    title('Left Hemi: Day vs Night', 'Interpreter','none');
    xlabel('LFP amplitude');
    ylabel('Probability');
    legend;
else
    title('Left Hemi: Day vs Night - No Data');
end

% --- 3. Right: Day vs Night ---
subplot(2,2,4);
if isfield(data_struct, 'Right')
    histogram(data_struct.Right.day, 'Normalization', 'probability', 'FaceAlpha', 0.6, 'DisplayName', 'Day');
    hold on;
    histogram(data_struct.Right.night, 'Normalization', 'probability', 'FaceAlpha', 0.6, 'DisplayName', 'Night');
    title('Right Hemi: Day vs Night', 'Interpreter','none');
    xlabel('LFP amplitude');
    ylabel('Probability');
    legend;
else
    title('Right Hemi: Day vs Night - No Data');
end

sgtitle(sprintf('LFP Band Histogram Comparison - %s', patient_name), 'Interpreter','none');

end
