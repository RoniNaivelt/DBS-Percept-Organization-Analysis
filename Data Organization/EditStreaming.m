function Streaming_table = EditStreaming(data, file_name, DBS_data)

%% Open current Events
Patiant_Initials = [data.PatientInformation.Initial.PatientFirstName(1) data.PatientInformation.Initial.PatientLastName(1)];
Patient_name = ['Patient_' Patiant_Initials];
Streaming_old = DBS_data.(Patient_name).Streaming;


%% Create new table

if isfield(data, 'IndefiniteStreaming') && ~isempty(fieldnames(data.IndefiniteStreaming))
    data_streaming = data.IndefiniteStreaming;
    nstreams = size(data_streaming, 1);
    Streaming_new = table('Size',[nstreams 10],'VariableTypes',...
        {'datetime', 'double', 'string', 'string', 'string', 'double', 'double', 'cell','cell', 'cell'},...
        'VariableNames',{'StartTime','Duration_Sec', 'ActivityType', 'Hemisphere', 'Contacts', 'Gain', 'SampleRate', 'TimeData', 'FFTPower', 'FFTFreq'});



    for row_num = 1:nstreams
        row_data = data_streaming(row_num);
        
        % 1. Start Time
        Streaming_new.StartTime(row_num) = datetime(regexprep(row_data.FirstPacketDateTime(1:end-1),'T',' '));

        % 2. Duration (in seconds)
        Streaming_new.Duration_Sec(row_num) = size(row_data.TimeDomainData,1) / row_data.SampleRateInHz;

        % 3. Activity type
        % we extract the text inside the file name - after the date and
        % before the '.json' ending
        if contains(file_name, 'json')
            file_name = extractBefore(file_name, '.json');
        end
        [tokens, startIdx, endIdx] = regexp(file_name, '\d+', 'match', 'start', 'end');
        if ~isempty(endIdx)
            last_end_idx = endIdx(end); % Get the end index of the last number
            ActivityType = file_name(last_end_idx+2:end); % Extract substring after the last number
        else
            ActivityType = '';
        end

        Streaming_new.ActivityType(row_num) = ActivityType;

        % 4-5. Hemisphere and Contacts
        Channel_data = row_data.Channel;
        parsed_channel_data = strsplit(Channel_data, '_');
        
        Hemisphere = parsed_channel_data{3};
        Contacts = [parsed_channel_data{1} '_' parsed_channel_data{2}];

        Streaming_new.Hemisphere(row_num) = Hemisphere;
        Streaming_new.Contacts(row_num) = Contacts;

        % 6. Gain 
        Streaming_new.Gain(row_num) = row_data.Gain;

        % 7. Sample rate in Hz
        Streaming_new.SampleRate(row_num) = row_data.SampleRateInHz;

        % 8. Time Domain Data
        Streaming_new.TimeData{row_num} = row_data.TimeDomainData;

        % 9. PSD data (fft)
        [pxx, f] = ExtractPSD_Streaming(row_data.TimeDomainData, row_data.SampleRateInHz);

        Streaming_new.FFTPower{row_num} = pxx;
        Streaming_new.FFTFreq{row_num} = f;


    end



    %% Add new events
    full_table = [Streaming_old; Streaming_new];

    % Create a unique identifier combining StartTime, Hemisphere, and Contacts
    unique_combination = [datestr(full_table.StartTime), full_table.ActivityType, full_table.Hemisphere, full_table.Contacts];

    % Identify Duplicates based on the combination of StartTime, Hemisphere, and Contacts
    [unique_comb_dates, ~, idx] = unique(unique_combination, 'rows', 'stable');
    duplicates = histcounts(idx, 1:numel(unique_comb_dates)) > 1;
    duplicate_comb_dates = unique_comb_dates(duplicates, :);

    % Process Duplicates
    for i = 1:size(duplicate_comb_dates, 1)
        combination = duplicate_comb_dates(i, :);

        % Find the rows with the duplicate combination
        duplicate_streams = full_table(full_table.StartTime == combination(1) & ...
            full_table.ActivityType == combination(2) &...
            full_table.Hemisphere == combination(3) & ...
            full_table.Contacts == combination(4), :);

        % Find the first non-empty index (if any)
        non_empty_indx = find(~isempty(duplicate_streams.TimeData), 1);

        if ~isempty(non_empty_indx)
            stream_to_keep = duplicate_streams(non_empty_indx, :); % Keep the one with TimeData = 1
        else
            stream_to_keep = duplicate_streams(1, :); % If none have TimeData = 1, keep the first
        end

        % Remove duplicate rows and append the chosen one
        full_table(full_table.StartTime == combination(1) & ...
            full_table.ActivityType == combination(2) &...
            full_table.Hemisphere == combination(3) & ...
            full_table.Contacts == combination(4), :) = [];
        full_table = [full_table; stream_to_keep]; % Append the chosen event
    end

    % Remove duplicates and keep only unique rows, sorted by StartTime
    [~, idx_unique] = unique(full_table(:, {'StartTime', 'ActivityType', 'Hemisphere', 'Contacts'}), 'rows', 'stable');
    Streaming_table = sortrows(full_table(idx_unique, :), 'StartTime'); % Sort by StartTime



else
   Streaming_table = Streaming_old; 
end