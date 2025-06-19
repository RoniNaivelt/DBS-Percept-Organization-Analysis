function Streaming_table = CreateStreaming(data, file_name)

if isfield(data, 'IndefiniteStreaming') && ~isempty(fieldnames(data.IndefiniteStreaming))
    data_streaming = data.IndefiniteStreaming;
    nstreams = size(data_streaming, 1);
    Streaming_table = table('Size',[nstreams 10],'VariableTypes',...
        {'datetime', 'double', 'string', 'string', 'string', 'double', 'double', 'cell','cell', 'cell'},...
        'VariableNames',{'StartTime','Duration_Sec', 'ActivityType', 'Hemisphere', 'Contacts', 'Gain', 'SampleRate', 'TimeData', 'FFTPower', 'FFTFreq'});



    for row_num = 1:nstreams
        row_data = data_streaming(row_num);
        
        % 1. Start Time
        Streaming_table.StartTime(row_num) = datetime(regexprep(row_data.FirstPacketDateTime(1:end-1),'T',' '));

        % 2. Duration (in seconds)
        Streaming_table.Duration_Sec(row_num) = size(row_data.TimeDomainData,1) / row_data.SampleRateInHz;

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

        Streaming_table.ActivityType(row_num) = ActivityType;

        % 4-5. Hemisphere and Contacts
        Channel_data = row_data.Channel;
        parsed_channel_data = strsplit(Channel_data, '_');
        
        Hemisphere = parsed_channel_data{3};
        Contacts = [parsed_channel_data{1} '_' parsed_channel_data{2}];

        Streaming_table.Hemisphere(row_num) = Hemisphere;
        Streaming_table.Contacts(row_num) = Contacts;

        % 6. Gain 
        Streaming_table.Gain(row_num) = row_data.Gain;

        % 7. Sample rate in Hz
        Streaming_table.SampleRate(row_num) = row_data.SampleRateInHz;

        % 8. Time Domain Data
        Streaming_table.TimeData{row_num} = row_data.TimeDomainData;

        % 9. PSD data (fft)
        [pxx, f] = ExtractPSD_Streaming(row_data.TimeDomainData, row_data.SampleRateInHz);

        Streaming_table.FFTPower{row_num} = pxx;
        Streaming_table.FFTFreq{row_num} = f;


    end
else
   Streaming_table = struct(); 
end


end


