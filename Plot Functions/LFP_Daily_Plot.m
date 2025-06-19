function LFP_Daily_Plot(DBS_data, patient_name, hemi, sensing_range, day_range, varargin)
% LFP_Daily_Plot -
%
% Inputs:
%   DBS_data            - Nested data structure containing DBS data
%   patient_name        - String or cell array of strings (must match DBS_data entries)
%   hemi                - (Optional) 'Left', 'Right', or 'both' (default)
%   sensing_range       - (Optional) array indicating a band of sensing
%                         frequencies such as [Min_frequency Max_frequency]
%   'day_range'         - [start end], or cell array of ranges per patient

% Name-Value Parameters:
%   'Plot_TEED'          - If false plot LFP, If true plot TEED
%   'Show_Events'        - If false plot only daily LFP, if true, show
%                          event snapshots
%   'ShowEventTypes'     - Plot each event as a vertical line, indicating
%                          the time when it occurred
%   'ShowEventBands'     - Plot each event as several dots, each indicating
%                          different bandwidth


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
addParameter(p, 'Plot_TEED', false);
addParameter(p, 'Show_Events', true);
addParameter(p, 'ShowEventTypes', true);
addParameter(p, 'ShowEventBands', false);
parse(p, varargin{:});

Plot_TEED = p.Results.Plot_TEED;
Show_Events = p.Results.Show_Events;
ShowEventTypes = p.Results.ShowEventTypes;
ShowEventBands = p.Results.ShowEventBands;



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
            % For LFP / TEED Plotting
            LFP_vec = data_out.(key).LFP_vec_trimmed;
            sensing_vec = data_out.(key).sensing_vec_trimmed;
            contacts_vec = data_out.(key).contacts_vec_trimmed;
            TEED_vec = data_out.(key).TEED_vec_trimmed;
            days_since_surg = data_out.(key).days_since_surg;
            % day mean vecs
            LFP_day_mean = mean(data_out.(key).LFP_day,1);
            TEED_day_mean = mean(data_out.(key).TEED_day,1);
            days_since_surg_vec = data_out.(key).days_since_surg_vec;
            % For Events
            Events = data_out.(key).Events;
            Events_Days_Since_Surge = data_out.(key).Events_Days_Since_Surge;
            Event_Summary = data_out.(key).Event_Summary;
        catch
            warning('Failed to extract LFP data for %s (%s hemisphere). Skipping.', name, h_name);
            continue;
        end

        
        % === Plotting ===
        figure('Color', 'w');

        % === Plot TEED ===
        yyaxis right
        if Plot_TEED
            plot(days_since_surg_vec, TEED_day_mean, '--o', 'LineWidth',2, 'DisplayName', 'TEED');
            ylabel('TEED [E/sec]');
        else
            plot(days_since_surg_vec, LFP_day_mean, '--o', 'LineWidth',2, 'DisplayName', 'LFP');
            ylabel('LFP');

        end

        xlabel('Days Since Surgery');
        hold on

        % === Plot Events (optional) ===
        if Show_Events
            if ShowEventBands
                yyaxis left
                colors = jet(4);
                ylabel('Relative AUC');

                % Plot each frequency band's AUC as a scatter
                plot(Event_Summary.Event_dates_since_surg, Event_Summary.Beta_AUC,  '*', 'Color', colors(1,:), 'LineWidth',2, 'DisplayName', 'Beta');
                plot(Event_Summary.Event_dates_since_surg, Event_Summary.Theta_AUC, '*', 'Color', colors(2,:), 'LineWidth',2, 'DisplayName', 'Theta');
                plot(Event_Summary.Event_dates_since_surg, Event_Summary.Alpha_AUC, '*', 'Color', colors(3,:), 'LineWidth',2, 'DisplayName', 'Alpha');
                plot(Event_Summary.Event_dates_since_surg, Event_Summary.Gamma_AUC, '*', 'Color', colors(4,:), 'LineWidth',2, 'DisplayName', 'Gamma');
     
            elseif ShowEventTypes
                % Get unique event types
                unique_events = unique(Events.EventName);
                colors = turbo(numel(unique_events));

                % Assign each event type a color
                for j = 1:length(unique_events)
                    current_event = unique_events{j};
                    current_color = colors(j,:);
                    event_idx = strcmp(Events.EventName, current_event);
                    current_dates = Events_Days_Since_Surge(event_idx);
                    for k = 1:length(current_dates)
                        if k == 1
                            % Only add to legend for the first occurrence
                            xline(current_dates(k), '--', 'Color', current_color, 'LineWidth', 2, ...
                                'DisplayName', current_event);
                        else
                            xline(current_dates(k), '--', 'Color', current_color, 'LineWidth', 2, ...
                                'HandleVisibility', 'off');
                        end
                    end
                end
            end
        end

        grid on
        lgd = legend(unique_events, 'Location','northeastoutside');
        set(lgd, 'Color', 'none');
        title([name ' - ' h_name ' Hemisphere Events'], 'FontSize',14, 'Interpreter', 'none');


        % find unique periods
        combined_strings = strcat(num2str(sensing_vec(:)), '_', contacts_vec(:));
        % Filter out missing entries first
        valid_idx = ~ismissing(combined_strings);
        clean_strings = combined_strings(valid_idx);
        [unique_pairs, ~, pair_idx_clean] = unique(clean_strings, 'stable');

        % === Draw sensing frequency patches ===
        if numel(unique_pairs) > 1
            yLimits = ylim;
            clr = jet(length(unique_pairs));
            for i = 1:length(unique_pairs)
                current_pair = unique_pairs(i);
                parts = split(current_pair, '_');
                if numel(parts) < 4
                    continue;  % skip malformed entries
                end
                current_freq = parts(1);
                contact1 = parts(2);
                contact2 = parts(4);

                idx = (pair_idx_clean == i);
                sensing_times = days_since_surg(idx);

                if isempty(sensing_times)
                    continue
                end

                date_start = sensing_times(1);
                date_end   = sensing_times(end);

                patch([date_start date_start date_end date_end], ...
                    [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], ...
                    clr(i,:), 'FaceAlpha', 0.2, ...
                    'DisplayName', ['SF - ' char(current_freq) ' Hz; Contacts ' char(contact1) '-' char(contact2)]);
            end
        else
            title([name ' - ' h_name ' ; ' char(unique_pairs)], 'FontSize', 14);
        end

    end
end
end