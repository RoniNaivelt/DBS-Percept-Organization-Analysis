 function Event_table = Convert_Event_to_feature(DBS_data,Patient_name, fill_in_mode, bins)
%%%%% fill in modes: %%%%%%
% 'Min_dopa' - adds ones from 2hr to 1 hr before the postmed event (which is 
% 1hr after medication intake).
% 'Max_dopa' -  rising exponentialy untill the postmed event, than
% decreases exponentially for 5hrs.
% 'Tremor' - add '1' in the bun before the event, and after the event -
% total three bin marking 30 minutes arround
% 'Fog' - add '1' in the bun before the event, and the bin itself - total of 20 minutes arround

%%%%% Bins %%%%%
% input as [bins before bins after] - only is needed to chnge the default
% values.
% for example - to mark only the event itself - use [0 0]

% Patient_initials = 'RN';
% Event_name = 'Premed';
% bins_before = 9;
% bins_after = -3;
%fill_in_mode = 'increasing_gaussian';
% if nargin < 6
%     fill_in_mode = 'Binary';  % Set a default value
% end

% 1. Initiation of parameters:
LFP_right = DBS_data.(Patient_name).TrendLogs.Right_Hemi.LFP_table;
Events_data = DBS_data.(Patient_name).Events.Events_data;
cols_names = LFP_right.Properties.VariableNames;
rows_names = LFP_right.Properties.RowNames;
n_days = size(LFP_right,2); %  = no. of columns 


% 2. time bins vector:
timebins = {};
for i = 1:length(cols_names)
    for j = 1:length(rows_names)
        timebins{end+1} = [cols_names{i} ' ' rows_names{j}];
    end
end


% 3. Event feature vector:
Event_feature = zeros(size(timebins));

% fill in according to the fill in method:
if strcmp(fill_in_mode, 'Min_dopa') == 1
    for i = 1:size(Events_data,1)
        % for every event, check if it is the type of event we are looking for
        if strcmp(Events_data.EventName(i), 'Med') || strcmp(Events_data.EventName(i), 'Took Medication') ...
                || strcmp(Events_data.EventName(i), 'Medication')
            date_time = Events_data.DateTime(i);
            % column
            char_date = char(date_time);
            col_name = char_date(1:11);
            % row
            row_name = char_date(13:17);
            if row_name(5) ~= '0'
                row_name(5) = '0';
            end
            event_time = [col_name ' ' row_name];
            event_time_idx = find(strcmp(timebins, event_time));
            if nargin<4
                bins_before = 6;
                bins_after = 0;
            else
                bins_before = bins(1);
                bins_after = bins(2);
            end
            
            Event_feature(event_time_idx-bins_before:event_time_idx+bins_after) = 1;
        end
        if strcmp(Events_data.EventName(i), 'Postmed')
            date_time = Events_data.DateTime(i);
            % column
            char_date = char(date_time);
            col_name = char_date(1:11);
            % row
            row_name = char_date(13:17);
            if row_name(5) ~= '0'
                row_name(5) = '0';
            end
            event_time = [col_name ' ' row_name];
            event_time_idx = find(strcmp(timebins, event_time));
            if nargin<4
                bins_before = 12;
                bins_after = -6;
            else
                bins_before = bins(1);
                bins_after = bins(2);
            end
            
            Event_feature(event_time_idx-bins_before:event_time_idx+bins_after) = 1;
    
        end
    end
    
    
elseif strcmp(fill_in_mode, 'Max_dopa') == 1
    for i = 1:size(Events_data,1)
        % for every event, check if it is the type of event we are looking for
        if strcmp(Events_data.EventName(i), 'Postmed')
            date_time = Events_data.DateTime(i);
            % column
            char_date = char(date_time);
            col_name = char_date(1:11);
            % row
            row_name = char_date(13:17);
            if row_name(5) ~= '0'
                row_name(5) = '0';
            end
            event_time = [col_name ' ' row_name];
            event_time_idx = find(strcmp(timebins, event_time));

            if nargin<4
                bins_before = 6;
                bins_after = 30;
            else
                bins_before = bins(1);
                bins_after = bins(2);
            end
            
            a = linspace(0, 1,7);
            a_exp = (exp(a)-1)/max(exp(a)-1);
            b = linspace(1,0,31);
            b_exp = (exp(b)-1)/max(exp(b)-1);
            final_vector = [a_exp(1:end-1) b_exp];
            Event_feature(event_time_idx-bins_before:event_time_idx+bins_after) = final_vector;
        end
    end


elseif strcmp(fill_in_mode, 'Tremor') == 1
    for i = 1:size(Events_data,1)
        % for every event, check if it is the type of event we are looking for
        if strcmp(Events_data.EventName(i), 'Tremor')
            date_time = Events_data.DateTime(i);
            % column
            char_date = char(date_time);
            col_name = char_date(1:11);
            % row
            row_name = char_date(13:17);
            if row_name(5) ~= '0'
                row_name(5) = '0';
            end
            event_time = [col_name ' ' row_name];
            event_time_idx = find(strcmp(timebins, event_time));

            if nargin<4
                bins_before = 1;
                bins_after = 1;
            else
                bins_before = bins(1);
                bins_after = bins(2);
            end
            
            Event_feature(event_time_idx-bins_before:event_time_idx+bins_after) = 1;
        end
    end
elseif strcmp(fill_in_mode, 'Fog') == 1
    for i = 1:size(Events_data,1)
        % for every event, check if it is the type of event we are looking for
        if strcmp(Events_data.EventName(i), 'Tremor')
            date_time = Events_data.DateTime(i);
            % column
            char_date = char(date_time);
            col_name = char_date(1:11);
            % row
            row_name = char_date(13:17);
            if row_name(5) ~= '0'
                row_name(5) = '0';
            end
            event_time = [col_name ' ' row_name];
            event_time_idx = find(strcmp(timebins, event_time));

            if nargin<4
                bins_before = 1;
                bins_after = 0;
            else
                bins_before = bins(1);
                bins_after = bins(2);
            end
            
            Event_feature(event_time_idx-bins_before:event_time_idx+bins_after) = 1;
        end
    end
end

% 4. creating a table for the event:
vartypes = cell(1, n_days);
for i = 1:n_days
    vartypes{i} = 'double';
end

Event_table = table('Size',[144 n_days],...
    'VariableTypes',vartypes,...
    'VariableNames',cols_names, 'RowNames', cellstr(rows_names));

%%% Fill in the table:
for i = 1:length(Event_feature)
    Event_Time = timebins{i};
    col_name = Event_Time(1:11);
    row_name = Event_Time(end-4:end);
    col_idx = find(strcmp(cols_names,col_name));
    row_idx = find(strcmp(cellstr(rows_names),row_name));
    Event_table{row_idx,col_idx} = Event_feature(i);
end

end

