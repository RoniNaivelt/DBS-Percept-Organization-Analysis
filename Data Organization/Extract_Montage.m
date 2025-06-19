function Events_data = Extract_Montage(data, DBS_data)

%%%%%%%%%%%%%%%% Initial Params %%%%%%%%%%%%%%%%
try
    Patiant_Initials = [data.PatientInformation.Initial.PatientFirstName(1) data.PatientInformation.Initial.PatientLastName(1)];
catch
    Patiant_Initials = strsplit(data.PatientInformation.Initial.PatientLastName, ' ');
    Patiant_Initials = [Patiant_Initials{1}(1) Patiant_Initials{2}(1)];
end
Patient_name = ['Patient_' Patiant_Initials];

try
    TrendLogs = DBS_data.(Patient_name).TrendLogs;
    Groups = DBS_data.(Patient_name).Groups;
    Events_data = DBS_data.(Patient_name).Events.Events_data;
catch
    TrendLogs = struct();
    Groups = struct();
    Events_data = table();
end

session_date = datetime([data.SessionDate(1:10) ' ' data.SessionDate(12:15) '0']);



if isfield(data, 'LFPMontage') || isfield(data, 'LfpMontageTimeDomain')

    Montage_table = table('Size',[1 8],'VariableTypes',...
        {'cell', 'double', 'cell', 'logical', 'logical', 'logical', 'struct','struct'},...
        'VariableNames',{'DateTime','EventID','EventName','LFP','Cycling', 'Snapshot', 'Left_Snapshot', 'Right_Snapshot'});

    Montage_table.DateTime = session_date;
    Montage_table.EventID = 42;
    Montage_table.EventName = {'Montage'};
    Montage_table.LFP = 0;
    Montage_table.Cycling = 0;
    Montage_table.Snapshot = 0;


    hemisphere_names = {'Left', 'Right'};
    for hemiIdx = 1:numel(hemisphere_names)
        hemi = hemisphere_names{hemiIdx};
        % make empty fields
        Montage_table.([hemi '_Snapshot']).FFT_Val = [];
        Montage_table.([hemi '_Snapshot']).FFT_Freq = [];
        Montage_table.([hemi '_Snapshot']).SensingElectrode = [];
        Montage_table.([hemi '_Snapshot']).Theta_data = [];
        Montage_table.([hemi '_Snapshot']).Alpha_data = [];
        Montage_table.([hemi '_Snapshot']).Beta_data = [];
        Montage_table.([hemi '_Snapshot']).Gamma_data = [];
        Montage_table.([hemi '_Snapshot']).TimeDomain = [];

        %% %%%%%%%%%%%%%% open current file %%%%%%%%%%%%%%%%
        % 1. Load vectors from DBS_data
        if ~isfield(TrendLogs, [hemi '_Hemi']) || ~isfield(Groups, [hemi '_Hemi'])
            try
                Group_indx = find([data.Groups.Initial.ActiveGroup]);
                Group_params = data.Groups.Initial(Group_indx).ProgramSettings.SensingChannel;
                if iscell(Group_params)
                    contact = regexp(Group_params{hemiIdx}.Channel, '\.', 'split');
                    contact = contact(end);
                else
                    contact = regexp(Group_params(hemiIdx).Channel, '\.', 'split');
                    contact = contact(end);
                end
            catch
                continue
            end
        else
            Time_vector = TrendLogs.([hemi '_Hemi']).Time_vector;
            Contacts = reshape(table2array(Groups.([hemi '_Hemi']).Contacts),1,[]);

            % 2. Find correct contacts based on session time
            time_diff = Time_vector - session_date;
            [min_diff,min_diff_loc] = min(abs(time_diff));
            if min_diff > hours(5)
                % stop the function if the session time does not correspond to
                % any recorded time in a window of 5 hours.
                continue
            end
            contact = Contacts(min_diff_loc);
            if ismissing(contact)
                nonNanIndices = find(~ismissing(Contacts));
                distances = abs(nonNanIndices - min_diff_loc);
                [~, minIdx] = min(distances);
                closest_idx = nonNanIndices(minIdx);
                contact = Contacts(closest_idx);
            end
        end

            %% %%%%%%%%%%% Find LFPmontage row %%%%%%%%%%%%%%%%%
            if isfield(data, 'LFPMontage')
                Montage_table.Snapshot = 1;

                % by contact
                if iscell(data.LFPMontage)
                    montage_electrodes_cell = {};
                    montage_hemi_cell = {};
                    for i = 1:length(data.LFPMontage)
                        montage_electrodes_cell{end+1} = data.LFPMontage{i}.SensingElectrodes;
                        montage_hemi_cell{end+1} = data.LFPMontage{i}.Hemisphere;

                    end
                else
                    montage_electrodes_cell = {data.LFPMontage.SensingElectrodes};
                    montage_hemi_cell = {data.LFPMontage.Hemisphere};
                end
                montage_electrodes = cellfun(@(x) split(x, '.'), montage_electrodes_cell, 'UniformOutput', false);
                montage_electrodes = cellfun(@(x) x{2}, montage_electrodes, 'UniformOutput', false);
                row_index_by_contact = find(strcmp(montage_electrodes, contact));

                % by hemisphere
                montage_hemi = cellfun(@(x) split(x, '.'), montage_hemi_cell, 'UniformOutput', false);
                montage_hemi = cellfun(@(x) x{2}, montage_hemi, 'UniformOutput', false);
                row_index_by_hemi = find(strcmp(montage_hemi, hemi));

                % combine the two
                row_index = intersect(row_index_by_hemi, row_index_by_contact);
                if ~isempty(row_index)
                    % extract data vectors from row
                    if iscell(data.LFPMontage)
                        Montage_Freq = data.LFPMontage{row_index}.LFPFrequency;
                        Montage_FFT = data.LFPMontage{row_index}.LFPMagnitude;
                    else
                        Montage_Freq = data.LFPMontage(row_index).LFPFrequency;
                        Montage_FFT = data.LFPMontage(row_index).LFPMagnitude;
                    end
                    % add fields to table
                    Montage_table.([hemi '_Snapshot']).FFT_Val = Montage_FFT;
                    Montage_table.([hemi '_Snapshot']).FFT_Freq = Montage_Freq;
                    Montage_table.([hemi '_Snapshot']).SensingElectrode = contact;

                    SnapPlot_data.FFT_Val = Montage_FFT;
                    SnapPlot_data.FFT_Freq = Montage_Freq;

                    [Theta_data, Alpha_data, Beta_data, Gamma_data] = Extract_SnapPlot(SnapPlot_data);
                    Montage_table.([hemi '_Snapshot']).Theta_data = Theta_data;
                    Montage_table.([hemi '_Snapshot']).Alpha_data = Alpha_data;
                    Montage_table.([hemi '_Snapshot']).Beta_data = Beta_data;
                    Montage_table.([hemi '_Snapshot']).Gamma_data = Gamma_data;
                end
            end


            %% %%%%%%%%%%% Find LFPmontage TimeDomain row %%%%%%%%%%%%%%%%%
            if isfield(data, 'LfpMontageTimeDomain')
                % by contact
                montage_electrodes_cell = {data.LfpMontageTimeDomain.Channel};
                montage_electrodes = cellfun(@(x) split(x, '_'), montage_electrodes_cell, 'UniformOutput', false);
                montage_electrodes = cellfun(@(x) [x{1} '_' x{2} '_' x{3}], montage_electrodes, 'UniformOutput', false);
                row_index_by_contact = find(strcmp(montage_electrodes, contact));

                % by hemisphere
                montage_hemi_cell = {data.LfpMontageTimeDomain.Channel};
                montage_hemi = cellfun(@(x) split(x, '_'), montage_hemi_cell, 'UniformOutput', false);
                montage_hemi = cellfun(@(x) [x{4}], montage_hemi, 'UniformOutput', false);
                row_index_by_hemi = find(strcmp(lower(montage_hemi), lower(hemi)));

                % combine the two
                row_index = intersect(row_index_by_hemi, row_index_by_contact);
                if ~isempty(row_index)
                    % extract data vectors from row
                    Montage_Time = data.LfpMontageTimeDomain(row_index).TimeDomainData;

                    % add fields to table
                    exact_datetime = data.LfpMontageTimeDomain(row_index).FirstPacketDateTime;
                    Montage_table.DateTime = datetime([exact_datetime(1:10) ' ' exact_datetime(12:15) '0']);
                    Montage_table.([hemi '_Snapshot']).TimeDomain = Montage_Time;
                end
            end
        end

        %% Add monate row to Events table
        if isfield(DBS_data.(Patient_name), 'Events') && istable(DBS_data.(Patient_name).Events.Events_data)
            Events_old = DBS_data.(Patient_name).Events.Events_data;
            full_events = [Events_old; Montage_table];
            Events_data = sortrows(full_events, 1);
        else
            Events_data = Montage_table;
        end

    end
end
