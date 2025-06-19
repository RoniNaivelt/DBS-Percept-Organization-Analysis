function TrendLogs = CreateTrendLogs(data)
TrendLogs = struct();
% create tables and save the following info:
    % 1. LFP table
    % 2. Amplitude table
    % 3. statistics table
    % 4. plots
% does not fill in the struct if:
    % 1. no 'DiagnosticData' in data file
    % 2. no 'LFPTrendLogs' in data file
    % 3. LFP vector is empty
    % 4. ignore files with only one field in LFPTrendLogs (this data alway already exist elsewhere)

    

hemi_names = {'Right', 'Left'};
for hemiIdx = 1:numel(hemi_names)
    %%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%
    hemi = hemi_names{hemiIdx};
    hemi_name = sprintf('HemisphereLocationDef_%s', hemi);
    LFP = [];
    Amp = [];
    DateTime = [];

    if isfield(data.DiagnosticData.LFPTrendLogs, hemi_name) &&... % Check if data for this side is not empty
        numel(fieldnames(data.DiagnosticData.LFPTrendLogs.(hemi_name))) > 1 % ignore files with only one line (this data alway already exist elsewhere)

        %%%%%%%%%%%%%%%% Extracting data from file %%%%%%%%%%%%%%%%%%%

        data_hemi = data.DiagnosticData.LFPTrendLogs.(hemi_name);
        field_names = fieldnames(data_hemi);
        for j = 1:numel(field_names)
            data_field = data_hemi.(field_names{j});

            LFP = [LFP data_field.LFP];
            Amp = [Amp data_field.AmplitudeInMilliAmps];
            Date_time_temp = arrayfun(@(x) datetime(regexprep(x.DateTime(1:end-1),'T',' ')), data_field);
            DateTime = [DateTime; Date_time_temp];
        end


        %%%%%%%%%%%%%%%% Preprocessing %%%%%%%%%%%%%%%%%%%

        % 1. Conversion from [LSB] to [uVP] (= uV^2)
        LFP = LFP / 100;

        % 2. root transformation (to handle right skewed data) [= uV]
        LFP = sqrt(LFP);

        % 3. Denoising
        LFP_filtered = filloutliers(LFP, "spline", "movmedian", 24); % 4 hours window

        %%%%%%%%%%%%%%%% Add an organized table %%%%%%%%%%%%%%%%

        % 1. find unique dates (to be used as columns)
        DMY = cell(1, numel(DateTime));
        for i = 1:length(DateTime)
            char_date = char(DateTime(i));
            DMY{i} = char_date(1:11);
        end
        DMY_unique_raw = unique(DMY, 'stable');
        % Sort based on Month and day order:
        dateNumbers = datenum(DMY_unique_raw, 'dd-mmm-yyyy');
        [~, idx] = sort(dateNumbers);
        DMY_unique = DMY_unique_raw(idx);


        % 2 - creating hhmm format vector (to be used as rows)
        round_HHMM = cell(1, numel(DateTime));
        for i = 1:length(DateTime)
            char_date = char(DateTime(i));
            tmp = char_date(13:17);
            if tmp(5) ~= '0'
                tmp(5) = '0';
            end
            round_HHMM{i} = tmp;
        end
        HHMM = datetime(0, 0, 0, 0, 0:10:1430, 0);
        row_headers = datestr(HHMM, 'HH:MM'); % = rows names

        % 3. create tables and add row headers
        vartypes = cell(1, numel(DMY_unique));
        for i = 1:numel(DMY_unique)
            vartypes{i} = 'double';
        end
        %
        TrendLogs.([hemi '_Hemi']).LFP_table = table('Size', [144 numel(DMY_unique)],...
            'VariableTypes', vartypes,...
            'VariableNames', DMY_unique, 'RowNames', cellstr(row_headers));

        % 4. fill table
        for i = 1:length(LFP_filtered)
            col_name = DMY{i};
            row_name = round_HHMM{i};
            col_idx = find(strcmp(TrendLogs.([hemi '_Hemi']).LFP_table.Properties.VariableNames, col_name));
            row_idx = find(strcmp(TrendLogs.([hemi '_Hemi']).LFP_table.Properties.RowNames, row_name));
            TrendLogs.([hemi '_Hemi']).LFP_table{row_idx, col_idx} = LFP_filtered(i);
        end

        % fill empty cells with NaN
        TrendLogs.([hemi '_Hemi']).LFP_table{:,:}(TrendLogs.([hemi '_Hemi']).LFP_table{:,:} == 0) = NaN;

        % repeat table creation for StimAmp:
        TrendLogs.([hemi '_Hemi']).Amp_table = table('Size', [144 numel(DMY_unique)],...
            'VariableTypes', vartypes,...
            'VariableNames', DMY_unique, 'RowNames', cellstr(row_headers));

        for i = 1:length(Amp)
            col_name = DMY{i};
            row_name = round_HHMM{i};
            col_idx = find(strcmp(TrendLogs.([hemi '_Hemi']).Amp_table.Properties.VariableNames, col_name));
            row_idx = find(strcmp(TrendLogs.([hemi '_Hemi']).Amp_table.Properties.RowNames, row_name));
            TrendLogs.([hemi '_Hemi']).Amp_table{row_idx, col_idx} = Amp(i);
        end

        % fill empty cells with NaN
        TrendLogs.([hemi '_Hemi']).Amp_table{:,:}(TrendLogs.([hemi '_Hemi']).LFP_table{:,:} == 0) = NaN;

        % 5. Statistics
        TrendLogs.([hemi '_Hemi']).LFP_Statistics = table('Size', [144 5],...
            'VariableTypes', {'double', 'double', 'double', 'double', 'double'},...
            'VariableNames', {'min', 'max', 'mean', 'median', 'STD'},...
            'RowNames', cellstr(row_headers));

        TrendLogs.([hemi '_Hemi']).LFP_Statistics{:, 1} = min(TrendLogs.([hemi '_Hemi']).LFP_table{:,:}, [], 2, 'omitnan');
        TrendLogs.([hemi '_Hemi']).LFP_Statistics{:, 2} = max(TrendLogs.([hemi '_Hemi']).LFP_table{:,:}, [], 2, 'omitnan');
        TrendLogs.([hemi '_Hemi']).LFP_Statistics{:, 3} = mean(TrendLogs.([hemi '_Hemi']).LFP_table{:,:}, 2, 'omitnan');
        TrendLogs.([hemi '_Hemi']).LFP_Statistics{:, 4} = median(TrendLogs.([hemi '_Hemi']).LFP_table{:,:}, 2, 'omitnan');
        TrendLogs.([hemi '_Hemi']).LFP_Statistics{:, 5} = std(TrendLogs.([hemi '_Hemi']).LFP_table{:,:}, 0, 2, 'omitnan');

        % 6. Add time vector
        TrendLogs.([hemi '_Hemi']).Time_vector = DateTime;
    end
end


end
