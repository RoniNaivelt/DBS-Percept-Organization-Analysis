function Events_table = CreateEvents(data)

if isfield(data, 'DiagnosticData') && isfield(data.DiagnosticData, 'LFPTrendLogs')  &&...
        isfield(data.DiagnosticData, 'LfpFrequencySnapshotEvents')
    
    hemisphere_names = {'Left', 'Right'};
    %%%%%%%%%% Create Events table and add basic data from file %%%%%%%%%%%
    data_events = data.DiagnosticData.LfpFrequencySnapshotEvents;
    nEvents = size(data_events, 1);
    Events_table = table('Size',[nEvents 8],'VariableTypes',...
        {'cell', 'double', 'cell', 'logical', 'logical', 'logical', 'struct','struct'},...
        'VariableNames',{'DateTime','EventID','EventName','LFP','Cycling', 'Snapshot', 'Left_Snapshot', 'Right_Snapshot'});
    
    for event_num = 1:nEvents
        if iscell(data_events) %depending on software version
            thisEvent = struct2table(data_events{event_num}, 'AsArray', true);
        else
            thisEvent = struct2table(data_events(event_num), 'AsArray', true);
        end
        Events_table(event_num, 1:5) = thisEvent(:, 1:5); %remove potential 'LfpFrequencySnapshotEvents'
        Events_table.Snapshot(event_num) = 0; % will change to 1 if it is an event with a snapshot


        
        for hemiIdx = 1:numel(hemisphere_names)
            hemi = hemisphere_names{hemiIdx};

            if any(ismember(thisEvent.Properties.VariableNames, 'LfpFrequencySnapshotEvents')) &&...
                isfield(thisEvent.LfpFrequencySnapshotEvents, ['HemisphereLocationDef_' hemi])
                Events_table.Snapshot(event_num) = 1;
                Events_table.([hemi '_Snapshot'])(event_num).FFT_Val = thisEvent.LfpFrequencySnapshotEvents.(['HemisphereLocationDef_' hemi]).FFTBinData;
                Events_table.([hemi '_Snapshot'])(event_num).FFT_Freq = thisEvent.LfpFrequencySnapshotEvents.(['HemisphereLocationDef_' hemi]).Frequency;
                Sensing = split(thisEvent.LfpFrequencySnapshotEvents.(['HemisphereLocationDef_' hemi]).SenseID, '.');
                Events_table.([hemi '_Snapshot'])(event_num).SensingElectrode = Sensing{2};
                Group = split(thisEvent.LfpFrequencySnapshotEvents.(['HemisphereLocationDef_' hemi]).GroupId, '.');
                Events_table.([hemi '_Snapshot'])(event_num).SensingElectrode = Group(2);
            
                [Theta_data, Alpha_data, Beta_data, Gamma_data] = Extract_SnapPlot(Events_table.([hemi '_Snapshot'])(event_num));
                Events_table.([hemi '_Snapshot'])(event_num).Theta_data = Theta_data;
                Events_table.([hemi '_Snapshot'])(event_num).Alpha_data = Alpha_data;
                Events_table.([hemi '_Snapshot'])(event_num).Beta_data = Beta_data;
                Events_table.([hemi '_Snapshot'])(event_num).Gamma_data = Gamma_data;
                Events_table.([hemi '_Snapshot'])(event_num).TimeDomain = [];
            
            end
        end
    end
    
    % adjust Datetime column to real datetime value
    Events_table.DateTime = cellfun(@(x) datetime(regexprep(x(1:end-1),'T',' ')), Events_table.DateTime);
else
   Events_table = struct(); 
end
end

