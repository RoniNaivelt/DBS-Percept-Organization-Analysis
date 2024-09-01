function Groups = EditGroups(data, DBS_data)
%% Open current TrendLogs
Patiant_Initials = [data.PatientInformation.Initial.PatientFirstName(1) data.PatientInformation.Initial.PatientLastName(1)];
Patient_name = ['Patient_' Patiant_Initials];
Groups = DBS_data.(Patient_name).Groups;


hemisphere_names = {'Left', 'Right'};
for hemiIdx = 1:numel(hemisphere_names)
hemi = hemisphere_names{hemiIdx};

    if isfield(data.DiagnosticData.LFPTrendLogs, ['HemisphereLocationDef_' hemi]) &&... % Check if data for this side is not empty
        numel(fieldnames(data.DiagnosticData.LFPTrendLogs.(['HemisphereLocationDef_' hemi]))) > 1 % ignore files with only one line (this data alway already exist elsewhere)

        %%%%%%%%%%%%%%%% open current file %%%%%%%%%%%%%%%%
        data_hemisphere = data.DiagnosticData.LFPTrendLogs.(['HemisphereLocationDef_' hemi]);
        Field_dates = fieldnames(data_hemisphere);
        n_fields = numel(Field_dates);
        DateTime = table;    

         % Step 1: find the initial group
        Group_indx = find([data.Groups.Initial.ActiveGroup]);  
        Group_params = data.Groups.Initial(Group_indx).ProgramSettings.SensingChannel;

        % Step 2: Create a DateTime vector
        for idx = 1:n_fields
            data_field = data_hemisphere.(Field_dates{idx});
            field = struct2table(data_hemisphere.(Field_dates{idx}));
            DateTime = [DateTime ; field.DateTime];
        end
        DateTime = sortrows(DateTime, 1);
        DateTime = table2array(DateTime);
        DateTime = arrayfun(@(Idx) datetime(regexprep(DateTime{Idx,:}(1:end-1),'T',' ')), 1:size(DateTime, 1), 'UniformOutput', false);


    %%%%%%%%%%%%%%%% Updated table Parameters %%%%%%%%%%%%%%%%

        % 1. Column headers

        old_cols = Groups.([hemi '_Hemi']).Pulse_width.Properties.VariableNames;
        % new cols:
        DMY = [];
        for i = 1:length(DateTime)
            char_date = char(DateTime{i});
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

        % 2. Row headers:

        HHMM = datetime(0,0,0,0,0:10:1430,0);
        row_headers = datestr(HHMM, 'HH:MM'); % = rows names

        % vartypes: % a list with a changing length is needed %
        vartypes = cell(1, numel(new_cols));
        for i = 1:numel(new_cols)
            vartypes{i} = 'double';
        end

    %%%%%%%%%%%%%%%% create the new data vectors %%%%%%%%%%%%%%%%
        numRows = 144;
        numCols = length(new_cols);
        if iscell(Group_params)
            Pulse_width_vec = repelem({Group_params{hemiIdx}.PulseWidthInMicroSecond}, numRows*numCols, 1);
            Stim_Rate_vec = repelem({Group_params{hemiIdx}.RateInHertz}, numRows*numCols, 1);
            Sensing_Freq_vec = repelem({Group_params{hemiIdx}.SensingSetup.FrequencyInHertz}, numRows*numCols, 1);
            contact = regexp(Group_params{hemiIdx}.Channel, '\.', 'split');
            Contacts_vec = repelem(contact(end), numRows*numCols, 1);
        else  
            Pulse_width_vec = repelem({Group_params(hemiIdx).PulseWidthInMicroSecond}, numRows*numCols, 1);
            Stim_Rate_vec = repelem({Group_params(hemiIdx).RateInHertz}, numRows*numCols, 1);
            Sensing_Freq_vec = repelem({Group_params(hemiIdx).SensingSetup.FrequencyInHertz}, numRows*numCols, 1);
            contact = regexp(Group_params(hemiIdx).Channel, '\.', 'split');
            Contacts_vec = repelem(contact(end), numRows*numCols, 1);
        end


        %% also extract impedance
        Impedance_struct = data.Impedance.Hemisphere;
        hemi_row = find(contains({Impedance_struct.Hemisphere}, hemi));
        Impedance_table = data.Impedance.Hemisphere(hemi_row).SessionImpedance.Monopolar;

        sensing_impedance = [];
        if iscell(Group_params)
            Sensing_electrodes = Group_params{hemiIdx}.ElectrodeState;
        else
            Sensing_electrodes = Group_params(hemiIdx).ElectrodeState;
        end
        for i = 1:length(Sensing_electrodes)-1
            electrode = strsplit(Sensing_electrodes{i}.Electrode, '_');
            electrode = electrode{end};
            electrode_row = find(contains({Impedance_table.Electrode2}, electrode));
            electrode_impedance = Impedance_table(electrode_row).ResultValue;
            sensing_impedance(end+1) = electrode_impedance;
        end
        Impedance = mean(sensing_impedance);
        Impedance_vec = repelem({Impedance}, numRows*numCols, 1);

    %%%%%%%%%%%%%%%% Update data %%%%%%%%%%%%%%%%
        % updating data in 3 steps -
        % 1. create new empty tables with the new col_names.
        % 2. fill table with values from old table.
        % 3. fill in new data using fill table loop.

        % 1. create tables and add row headers
        new_Pulse_width = table('Size',[numRows numCols],...
            'VariableTypes',vartypes,...
            'VariableNames',new_cols, 'RowNames', cellstr(row_headers));

        new_Stim_Rate = table('Size',[numRows numCols],...
            'VariableTypes',vartypes,...
            'VariableNames',new_cols, 'RowNames', cellstr(row_headers));

        new_Sensing_Freq = table('Size',[numRows numCols],...
            'VariableTypes',vartypes,...
            'VariableNames',new_cols, 'RowNames', cellstr(row_headers));

        new_Contacts = table('Size',[numRows numCols],...
            'VariableTypes',repmat({'string'},1,numel(new_cols)),...
            'VariableNames',new_cols, 'RowNames', cellstr(row_headers));

        new_Impedance = table('Size',[numRows numCols],...
            'VariableTypes',vartypes,...
            'VariableNames',new_cols, 'RowNames', cellstr(row_headers));


        % 2. fill table with values from old table
        if ~isempty(Groups.([hemi '_Hemi']).Pulse_width)

            old_Pulse_width = Groups.([hemi '_Hemi']).Pulse_width;
            old_Stim_Rate = Groups.([hemi '_Hemi']).Stim_Rate;
            old_Sensing_Freq = Groups.([hemi '_Hemi']).Sensing_Freq;
            old_Contacts = Groups.([hemi '_Hemi']).Contacts;
            old_Impedance = Groups.([hemi '_Hemi']).Impedance;

            for row = 1:size(row_headers,1) % rows
                for col = 1:size(old_cols,2) % cols
                    col_name = old_Pulse_width.Properties.VariableNames(col);
                    row_name = old_Pulse_width.Properties.RowNames(row);
                    %
                    col_idx = find(strcmp(new_cols,col_name));
                    row_idx = find(strcmp(cellstr(row_headers),row_name));
                    %
                    new_Pulse_width{row_idx,col_idx} = old_Pulse_width{row,col};
                    new_Stim_Rate{row_idx,col_idx} = old_Stim_Rate{row,col};
                    new_Sensing_Freq{row_idx,col_idx} = old_Sensing_Freq{row,col};
                    new_Contacts{row_idx,col_idx} = old_Contacts{row,col};
                    new_Impedance{row_idx,col_idx} = old_Impedance{row,col};
                end
            end
        end

        % 3. fill in new data using a loop
        %   find row for each LFP value:
        for i = 1:length(DateTime)
            char_date = char(DateTime{i});
            tmp = char_date(13:17);
            if tmp(5) ~= '0'
                tmp(5) = '0';
            end
            round_HHMM{i} = tmp;
        end

        %   fill table with new data:
        % DateTime - vector containing all timepoints in orde
        for i = 1:length(DateTime)
            col_name = DMY{i};
            row_name = round_HHMM{i};
            col_idx = find(strcmp(new_cols,col_name));
            row_idx = find(strcmp(cellstr(row_headers),row_name));

            new_Pulse_width{row_idx,col_idx} = Pulse_width_vec{i};
            new_Stim_Rate{row_idx,col_idx} = Stim_Rate_vec{i};
            new_Sensing_Freq{row_idx,col_idx} = Sensing_Freq_vec{i};
            new_Contacts{row_idx,col_idx} = Contacts_vec(i);
            new_Impedance{row_idx,col_idx} = Impedance_vec{i};
        end

        % fill empty cells with NaN
        new_Pulse_width{:,:}(new_Pulse_width{:,:} == 0) = NaN;
        new_Stim_Rate{:,:}(new_Stim_Rate{:,:} == 0) = NaN;
        new_Sensing_Freq{:,:}(new_Sensing_Freq{:,:} == 0) = NaN;

        % update Groups
        Groups.([hemi '_Hemi']).Pulse_width = new_Pulse_width;
        Groups.([hemi '_Hemi']).Stim_Rate = new_Stim_Rate;
        Groups.([hemi '_Hemi']).Sensing_Freq = new_Sensing_Freq;
        Groups.([hemi '_Hemi']).Contacts = new_Contacts;
        Groups.([hemi '_Hemi']).Impedance = new_Impedance;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Group Switching
%     for Events = data.DiagnosticData.EventLogs
%         if isfield(Events{i}.ParameterTrendId) && strcomp(Events{i}.ParameterTrendId, 'ParameterTrendIdDef.ActiveGroup') == 1
%             Change_timepoint = Events{i}.DateTime;
%             Groups = GroupChangeEvent(data, Groups);
%         end
%     end

        end
end

end
