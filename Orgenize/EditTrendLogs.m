function TrendLogs = EditTrendLogs(data, DBS_data)
%% Open current TrendLogs
Patiant_Initials = [data.PatientInformation.Initial.PatientFirstName(1) data.PatientInformation.Initial.PatientLastName(1)];
Patient_name = ['Patient_' Patiant_Initials];
TrendLogs = DBS_data.(Patient_name).TrendLogs;
session_date = [ 'session_' data.SessionDate(9:10) '_' data.SessionDate(6:7)];

%% open current file  

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


        %%%%%%%%%%%%%%%% Preprocessing new data %%%%%%%%%%%%%%%%%%%

        % 1. Conversion from [LSB] to [uVP] (= uV^2)
        LFP = LFP / 100;

        % 2. root transformation (to handle right skewed data) [= uV]
        LFP = sqrt(LFP);

        % 3. Denoising
        LFP_filtered = filloutliers(LFP, "spline", "movmedian", 24); % 4 hours window

        figure(); plot(LFP);
hold on; plot(LFP_filtered)


        %% %%%%%%%%%%%%% Update LFP table %%%%%%%%%%%%%%%%


        % we need to update 2 things - the columns and the content.
        %%%%%%%%%%%%%%%% Update columns %%%%%%%%%%%%%%%%
        old_cols = TrendLogs.([hemi '_Hemi']).LFP_table.Properties.VariableNames;
        % now cols:
        DMY = cell(1, numel(DateTime));
        for i = 1:length(DateTime)
            char_date = char(DateTime(i));
            DMY{i} = char_date(1:11);
        end
        % unite dates:
        cols_to_update = unique(DMY, 'stable');
        new_cols_raw = [old_cols cols_to_update];
        new_cols_raw = unique(new_cols_raw, 'stable');
        % Sort based on serial date numbers:
        dateNumbers = datenum(new_cols_raw, 'dd-mmm-yyyy');
        [~, idx] = sort(dateNumbers);
        new_cols = new_cols_raw(idx);

        % create row headers (just a vector of times - from 00:00 to 23:50)
        HHMM = datetime(0,0,0,0,0:10:1430,0);
        row_headers = datestr(HHMM, 'HH:MM'); % = rows names

        %%%%%%%%%%%%%%%% Update data %%%%%%%%%%%%%%%%
        % updating data in 3 steps -
        % 1. create new empty tables with the new col_names.
        % 2. fill table with values from old table.
        % 3. fill in new data using fill table loop.


        % 1. create tables and add row headers
        % vartypes: % a list with a changing length is needed 
        vartypes = cell(1, numel(new_cols));
        for i = 1:numel(new_cols)
            vartypes{i} = 'double';
        end
        %
        new_LFP_table = table('Size',[144 numel(new_cols)],...
            'VariableTypes',vartypes,...
            'VariableNames',new_cols, 'RowNames', cellstr(row_headers));


        % 2. fill table with values from old table
        old_LFP_table = TrendLogs.([hemi '_Hemi']).LFP_table;
        for row = 1:size(row_headers,1) % rows
            for col = 1:size(old_cols,2) % cols
                col_name = old_LFP_table.Properties.VariableNames(col);
                row_name = old_LFP_table.Properties.RowNames(row);
                %
                col_idx = find(strcmp(new_cols,col_name));
                row_idx = find(strcmp(cellstr(row_headers),row_name));
                % 
                new_LFP_table{row_idx,col_idx} = old_LFP_table{row,col};
            end
        end

        % 3. fill in new data using a loop (find row for each LFP value)
        for i = 1:length(DateTime)
            char_date = char(DateTime(i));
            tmp = char_date(13:17);
            if tmp(5) ~= '0'
                tmp(5) = '0';
            end
            round_HHMM{i} = tmp;
        end

        % 4.   fill table with new data:
        % DateTime - vector containing all timepoints in orde
        for i = 1:length(LFP)
            col_name = DMY{i};
            row_name = round_HHMM{i};
            col_idx = find(strcmp(new_cols,col_name));
            row_idx = find(strcmp(cellstr(row_headers),row_name));
            new_LFP_table{row_idx,col_idx} = LFP_filtered(i);
        end

        % 5. fill empty cells with NaN
        new_LFP_table{:,:}(new_LFP_table{:,:} == 0) = NaN;

        % 6. update TrendLogs
        TrendLogs.([hemi '_Hemi']).LFP_table = new_LFP_table;


        %%%%%%%%%%%%%%%% repeat table creation for StimAmp %%%%%%%%%%%%%%%%

        % 1. create tables and add row headers
        new_Amp_table = table('Size',[144 numel(new_cols)],...
            'VariableTypes',vartypes,...
            'VariableNames',new_cols, 'RowNames', cellstr(row_headers));

        % 2. fill table with values from old table
        old_Amp_table = TrendLogs.([hemi '_Hemi']).Amp_table;
        for row = 1:size(row_headers,1) % rows
            for col = 1:size(old_cols,2) % cols
                col_name = old_Amp_table.Properties.VariableNames(col);
                row_name = old_Amp_table.Properties.RowNames(row);
                %
                col_idx = find(strcmp(new_cols,col_name));
                row_idx = find(strcmp(cellstr(row_headers),row_name));
                new_Amp_table{row_idx,col_idx} = old_Amp_table{row,col};
            end
        end

        % 3. fill in new data using fill table loop (find row for each Amp value)
        for i = 1:length(DateTime)
            char_date = char(DateTime(i));
            tmp = char_date(13:17);
            if tmp(5) ~= '0'
                tmp(5) = '0';
            end
            round_HHMM{i} = tmp;
        end

        % 4. fill table with new data:
        for i = 1:length(Amp)
            col_name = DMY{i};
            row_name = round_HHMM{i};
            col_idx = find(strcmp(new_cols,col_name));
            row_idx = find(strcmp(cellstr(row_headers),row_name));
            new_Amp_table{row_idx,col_idx} = Amp(i);

        end

        % 5. fill empty cells with NaN
        new_Amp_table{:,:}(new_LFP_table{:,:} == 0) = NaN;

        % 6. update TrendLogs
        TrendLogs.([hemi '_Hemi']).Amp_table = new_Amp_table;


        %% %%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%
        TrendLogs.([hemi '_Hemi']).LFP_Statistics{:,1} = min(TrendLogs.([hemi '_Hemi']).LFP_table{:,:},[],2, 'omitnan');
        TrendLogs.([hemi '_Hemi']).LFP_Statistics{:,2} = max(TrendLogs.([hemi '_Hemi']).LFP_table{:,:},[],2, 'omitnan');
        TrendLogs.([hemi '_Hemi']).LFP_Statistics{:,3} = mean(TrendLogs.([hemi '_Hemi']).LFP_table{:,:},2,'omitnan');
        TrendLogs.([hemi '_Hemi']).LFP_Statistics{:,4} = median(TrendLogs.([hemi '_Hemi']).LFP_table{:,:},2,'omitnan');
        TrendLogs.([hemi '_Hemi']).LFP_Statistics{:,5} = std(TrendLogs.([hemi '_Hemi']).LFP_table{:,:},0,2,'omitnan');

        %% Add time vector
        dates = TrendLogs.([hemi '_Hemi']).LFP_table.Properties.VariableNames;
        hours = TrendLogs.([hemi '_Hemi']).LFP_table.Properties.RowNames;
        [datesGrid, hoursGrid] = meshgrid(dates, hours);
        dateTimeStrings = datesGrid + " " + hoursGrid;
        Datetime_vec = datetime(dateTimeStrings(:), 'InputFormat', 'dd-MMM-yyyy HH:mm', 'Format', 'dd-MMM-yyyy HH:mm');
        TrendLogs.([hemi '_Hemi']).Time_vector = Datetime_vec;
    end
end  

end
