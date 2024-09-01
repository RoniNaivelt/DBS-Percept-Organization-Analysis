function Groups = CreateGroups(data)
Groups = struct();

hemisphere_names = {'Left', 'Right'};
for hemiIdx = 1:numel(hemisphere_names)
hemi = hemisphere_names{hemiIdx};
Event_log = data.DiagnosticData.EventLogs;
switch_times = [];
new_groups = [];
new_group_row = [];

    if isfield(data.DiagnosticData.LFPTrendLogs, ['HemisphereLocationDef_' hemi]) &&... % Check if data for this side is not empty
        numel(fieldnames(data.DiagnosticData.LFPTrendLogs.(['HemisphereLocationDef_' hemi]))) > 1 % ignore files with only one line (this data alway already exist elsewhere)

        %%%%%%%%%%%%%%%% open current file %%%%%%%%%%%%%%%%
        data_hemisphere = data.DiagnosticData.LFPTrendLogs.(['HemisphereLocationDef_' hemi]);
        Field_dates = fieldnames(data_hemisphere);
        n_fields = numel(Field_dates);
        DateTime = table;    
        for idx = 1:n_fields
            data_field = data_hemisphere.(Field_dates{idx});
            field = struct2table(data_hemisphere.(Field_dates{idx}));
            DateTime = [DateTime ; field.DateTime];
        end
        DateTime = sortrows(DateTime, 1);
        DateTime = table2array(DateTime);
        DateTime = arrayfun(@(Idx) datetime(regexprep(DateTime{Idx,:}(1:end-1),'T',' ')), 1:size(DateTime, 1), 'UniformOutput', false);

         % Step 1: find the initial group and initial time
        Group_indx = find([data.Groups.Initial.ActiveGroup]);  
        Group_name = data.Groups.Initial(Group_indx).GroupId;
        split_string = split(Group_name, '_');
        Group_name = cell2mat(split_string(end));
        new_groups = [new_groups ; Group_name];
        new_group_row = [new_group_row Group_indx];
        switch_times = [switch_times ; DateTime{1}];

        % Step 2: find all group switching events
        for i = 1:length(Event_log)
            if isfield(Event_log{i}, 'ParameterTrendId') && contains(Event_log{i}.ParameterTrendId, 'ActiveGroup') == 1
                % time of change
                Change_timepoint = Event_log{i}.DateTime;
                Change_timepoint = datetime(regexprep(Change_timepoint(1:end-1),'T',' '));
                switch_times = [switch_times ; Change_timepoint];
                % new group
                new_group =  Event_log{i}.NewGroupId;
                split_string = split(new_group, '_');
                new_group = cell2mat(split_string(end));
                new_groups = [new_groups ; new_group];
                % new group row
                try
                    for i = 1:length(data.Groups.Initial)
                        Group_name = data.Groups.Initial(i).GroupId;
                        split_string = split(Group_name, '_');
                        Group_name = cell2mat(split_string(end));
                        if char(Group_name) == char(new_group)
                            Group_indx = i;
                        end
                    end
                catch
                    for i = 1:length(data.Groups.Final)
                        Group_name = data.Groups.Final(i).GroupId;
                        split_string = split(Group_name, '_');
                        Group_name = cell2mat(split_string(end));
                        if Group_name == new_group
                            Group_indx = i;
                        end
                    end
                end
                new_group_row = [new_group_row ; Group_indx];
            end
        end


        % Step 3: create empty tables
        % coulmns:
        DMY = {};
        for i = 1:length(DateTime)
            char_date = char(DateTime{i});
            DMY{i} = char_date(1:11);
        end
        DMY_unique_raw = unique(DMY, 'stable');
        % Sort based on Month and day order:
        dateNumbers = datenum(DMY_unique_raw, 'dd-mmm-yyyy');
        [~, idx] = sort(dateNumbers);
        DMY_unique = DMY_unique_raw(idx);

        % rows:
        HHMM = datetime(0,0,0,0,0:10:1430,0);
        row_headers = datestr(HHMM, 'HH:MM'); % = rows names

        % vartypes: 
        % a list with a changing length is needed %
        vartypes = cell(1, numel(DMY_unique));
        for i = 1:numel(DMY_unique)
            vartypes{i} = 'double';
        end
        %
        numRows = 144;
        numCols = length(DMY_unique);
        %
        Groups.([hemi '_Hemi']).Pulse_width = table('Size', [numRows numCols], ...
            'VariableTypes', vartypes, ...
            'VariableNames', DMY_unique, 'RowNames', cellstr(row_headers));

        Groups.([hemi '_Hemi']).Stim_Rate = table('Size', [numRows numCols], ...
            'VariableTypes', vartypes, ...
            'VariableNames', DMY_unique, 'RowNames', cellstr(row_headers));

        Groups.([hemi '_Hemi']).Sensing_Freq = table('Size', [numRows numCols], ...
            'VariableTypes', vartypes, ...
            'VariableNames', DMY_unique, 'RowNames', cellstr(row_headers));

        Groups.([hemi '_Hemi']).Contacts = table('Size', [numRows numCols], ...
            'VariableTypes', repmat({'string'}, 1, numel(DMY_unique)), ...
            'VariableNames', DMY_unique, 'RowNames', cellstr(row_headers));
        
        Groups.([hemi '_Hemi']).Impedance = table('Size', [numRows numCols], ...
            'VariableTypes', vartypes, ...
            'VariableNames', DMY_unique, 'RowNames', cellstr(row_headers));



        % Step 4: create constant data vectors
        % for i = 1:length(new_groups)
        %     try
        %         Group_params = data.Groups.Initial(new_group_row(i)).ProgramSettings.SensingChannel;
        %     catch
        %         Group_params = data.Groups.Final(new_group_row(i)).ProgramSettings.SensingChannel;
        %     end
       
        Group_params = data.Groups.Initial(new_group_row(1)).ProgramSettings.SensingChannel;

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

        %Step 5: filling the tables
        Groups.([hemi '_Hemi']).Pulse_width{:,:} = cell2mat(Pulse_width_vec(1));
        Groups.([hemi '_Hemi']).Stim_Rate{:,:} = cell2mat(Stim_Rate_vec(1));
        Groups.([hemi '_Hemi']).Sensing_Freq{:,:} = cell2mat(Sensing_Freq_vec(1));
        Groups.([hemi '_Hemi']).Contacts{:,:} = repmat(Contacts_vec(1), numRows, numel(DMY_unique));

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
        Impedance_vec = repelem(Impedance, numRows*numCols, 1);
        Groups.([hemi '_Hemi']).Impedance{:,:} = repmat(Impedance_vec(1), numRows, numel(DMY_unique));

    end
end

end
