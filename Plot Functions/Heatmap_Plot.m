function Heatmap_Plot(DBS_data, patient_name, hemi, day_range, varargin)
% Heatmap_Plot - Plots Heatmap of all LFP data of a given hemisphere.
%
% Inputs:
%   DBS_data            - Nested data structure containing DBS data
%   patient_name        - String patient identifier
%   hemi                - Hemisphere: 'Left' or 'Right'
%   day_range           - (Optional) array with [start_day, end_day]

% Name-Value Parameters:
%   none for this function

if nargin < 4
    day_range = [];
end

% Parse optional parameters
p = inputParser;
addParameter(p, 'colors', []);
parse(p, varargin{:});
colors = p.Results.colors;

% Determine color
if isempty(colors)
    colors = turbo(1);
end

% Get LFP table
LFP_table = DBS_data.(patient_name).TrendLogs.([hemi '_Hemi']).LFP_table;

if ~isempty(day_range)
    LFP_table = LFP_table(:,day_range(1):day_range(2));
end

% Get Surgery Date
surg_date = DBS_data.(patient_name).Info.deviceInfo.ImplantDate;
surg_date = datetime(regexprep(surg_date(1:end-1), 'T', ' '));

time_of_day = datetime(LFP_table.Properties.RowNames, 'Format', 'HH:mm');
date_strs = LFP_table.Properties.VariableNames;
date_nums = cellfun(@(x) datetime(x, 'InputFormat', 'dd-MMM-yy'), ...
    date_strs, 'UniformOutput', false);
days_since = round(days([date_nums{:}] - surg_date));
LFP_matrix = table2array(LFP_table);

% Plot
figure('Color', 'w');
h = heatmap(days_since, time_of_day, LFP_matrix);
h.Title = [hemi ' Hemisphere LFP Heatmap'];
h.XLabel = 'Days Since Surgery';
h.YLabel = 'Time of Day';
end