function Events_data = EditEvents(data, DBS_data)
%% Open current Events
Patiant_Initials = [data.PatientInformation.Initial.PatientFirstName(1) data.PatientInformation.Initial.PatientLastName(1)];
Patient_name = ['Patient_' Patiant_Initials];
Events_old = DBS_data.(Patient_name).Events.Events_data;
%% Open new Events
if isfield(data, 'DiagnosticData') && isfield(data.DiagnosticData, 'LFPTrendLogs')  &&...
        isfield(data.DiagnosticData, 'LfpFrequencySnapshotEvents')

    hemisphere_names = {'Left', 'Right'};
    %%%%%%%%%% Create Events table and add basic data from file %%%%%%%%%%%
    Data_events = data.DiagnosticData.LfpFrequencySnapshotEvents;
    nEvents = size(Data_events, 1);
    Events_new = table('Size',[nEvents 8],'VariableTypes',...
        {'cell', 'double', 'cell', 'logical', 'logical', 'logical', 'struct','struct'},...
        'VariableNames',{'DateTime','EventID','EventName','LFP','Cycling', 'Snapshot', 'Left_Snapshot', 'Right_Snapshot'});

    for event_num = 1:nEvents
        if iscell(Data_events) %depending on software version
            thisEvent = struct2table(Data_events{event_num}, 'AsArray', true);
        else
            thisEvent = struct2table(Data_events(event_num), 'AsArray', true);
        end
        Events_new(event_num, 1:5) = thisEvent(:, 1:5); %remove potential 'LfpFrequencySnapshotEvents'
        Events_new.Snapshot(event_num) = 0; % will change to 1 if it is an event with a snapshot


        for hemiIdx = 1:numel(hemisphere_names)
            hemi = hemisphere_names{hemiIdx};
            if any(ismember(thisEvent.Properties.VariableNames, 'LfpFrequencySnapshotEvents')) &&...
                isfield(thisEvent.LfpFrequencySnapshotEvents, ['HemisphereLocationDef_' hemi])

                Events_new.Snapshot(event_num) = 1;
                Events_new.([hemi '_Snapshot'])(event_num).FFT_Val = thisEvent.LfpFrequencySnapshotEvents.(['HemisphereLocationDef_' hemi]).FFTBinData;
                Events_new.([hemi '_Snapshot'])(event_num).FFT_Freq = thisEvent.LfpFrequencySnapshotEvents.(['HemisphereLocationDef_' hemi]).Frequency;
                Sensing = split(thisEvent.LfpFrequencySnapshotEvents.(['HemisphereLocationDef_' hemi]).SenseID, '.');
                Events_new.([hemi '_Snapshot'])(event_num).SensingElectrode = Sensing{2};
                Group = split(thisEvent.LfpFrequencySnapshotEvents.(['HemisphereLocationDef_' hemi]).GroupId, '.');
                Events_new.([hemi '_Snapshot'])(event_num).SensingElectrode = Group(2);

                [Theta_data, Alpha_data, Beta_data, Gamma_data] = Extract_SnapPlot(Events_new.([hemi '_Snapshot'])(event_num));
                Events_new.([hemi '_Snapshot'])(event_num).Theta_data = Theta_data;
                Events_new.([hemi '_Snapshot'])(event_num).Alpha_data = Alpha_data;
                Events_new.([hemi '_Snapshot'])(event_num).Beta_data = Beta_data;
                Events_new.([hemi '_Snapshot'])(event_num).Gamma_data = Gamma_data;
                Events_new.([hemi '_Snapshot'])(event_num).TimeDomain = [];
            else
                
                Events_new.([hemi '_Snapshot'])(event_num).FFT_Val = [];
                Events_new.([hemi '_Snapshot'])(event_num).FFT_Freq = [];
                Events_new.([hemi '_Snapshot'])(event_num).SensingElectrode = [];
                Events_new.([hemi '_Snapshot'])(event_num).SensingElectrode = [];
                Events_new.([hemi '_Snapshot'])(event_num).Theta_data = [];
                Events_new.([hemi '_Snapshot'])(event_num).Alpha_data = [];
                Events_new.([hemi '_Snapshot'])(event_num).Beta_data = [];
                Events_new.([hemi '_Snapshot'])(event_num).Gamma_data = [];
                Events_new.([hemi '_Snapshot'])(event_num).TimeDomain = [];
            end

        end
    end

    % adjust Datetime column to real datetime value
    Events_new.DateTime = cellfun(@(x) datetime(regexprep(x(1:end-1),'T',' ')), Events_new.DateTime);

    %% Add new events
    full_events = [Events_old; Events_new];
    % Identify Duplicate Dates:
    [unique_dates,~,idx] = unique(full_events.DateTime, 'stable');
    duplicates = histcounts(idx, 1:numel(unique_dates)) > 1;
    duplicate_dates = unique_dates(duplicates);
    % Process Duplicates and check LfpFrequencySnapshotEvents
    for i = 1:length(duplicate_dates)
        date = duplicate_dates(i);
        duplicate_events = full_events(full_events.DateTime == date, :);
        non_empty_indx = find(duplicate_events.Snapshot == 1, 1);
        if ~isempty(non_empty_indx)
            event_to_keep = duplicate_events(non_empty_indx,:);
        else
            event_to_keep = duplicate_events(1,:);
        end
        % Replace in full_events or create a new list to accumulate
        full_events(full_events.DateTime == date, :) = [];
        full_events = [full_events; event_to_keep]; % Append the chosen event
    end
    [~, idx_unique] = unique(full_events.DateTime, 'stable');
    Events_data = sortrows(full_events(idx_unique, :), 1);

else
    Events_data = Events_old;
end
end
