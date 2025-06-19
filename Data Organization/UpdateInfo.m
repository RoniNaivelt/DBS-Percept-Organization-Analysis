function DBS_data = UpdateInfo(data, DBS_data, file_name, file)
% get patient info
try
    Patient_Initials = [data.PatientInformation.Initial.PatientFirstName(1) data.PatientInformation.Initial.PatientLastName(1)];
catch
    Patient_Initials_split = strsplit(data.PatientInformation.Initial.PatientLastName, ' ');
    if numel(Patient_Initials_split) >= 2
        Patient_Initials = [Patient_Initials_split{1}(1) Patient_Initials_split{2}(1)];
    else
        Patient_Initials = [Patient_Initials_split{1}(1) 'X']; % default fallback
    end
end
Patient_name = ['Patient_' Patient_Initials];

%% Section 1: Update file info
dt = [data.SessionDate(1:10) ' ' data.SessionDate(12:end-1)];

if ~any(strcmp(DBS_data.(Patient_name).Info.Files{:, 1}, file_name))
    newRow = table(string(file_name), datetime(dt), 0, 0, 0, ...
        'VariableNames', {'Filename', 'Session Date', 'Trendlogs', 'Events', 'SnapshotEvents'});
    
    % 1. Check for Trendlogs (LFP Data)
    if isfield(data, 'DiagnosticData') &&...
       isfield(data.DiagnosticData,'LFPTrendLogs') &&...
       isfield(data.DiagnosticData.LFPTrendLogs, 'HemisphereLocationDef_Left') &&...
       isfield(data.DiagnosticData.LFPTrendLogs, 'HemisphereLocationDef_Right') &&...
       (size(fieldnames(data.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Left),1) > 1 ||...
       size(fieldnames(data.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Right),1) > 1)
        newRow.Trendlogs = true;
    end
    
    % 2. Check for Events
    if isfield(data, 'DiagnosticData') &&...
            isfield(data.DiagnosticData, 'LfpFrequencySnapshotEvents')
        newRow.Events = true;
        
        % 2.1. Check for SnapshotEvents
        if isfield(data.DiagnosticData.LfpFrequencySnapshotEvents, 'LfpFrequencySnapshotEvents')
            newRow.SnapshotEvents = true;
        end
    end
    % Append the new row to the table
    if isempty(DBS_data.(Patient_name).Info.Files)
        DBS_data.(Patient_name).Info.Files = newRow;
    else
        DBS_data.(Patient_name).Info.Files = [DBS_data.(Patient_name).Info.Files; newRow];
    end
end

%% Section 2: Update info about DBS_data
for hemi = {'Left', 'Right'}
    hemi_name = [cell2mat(hemi) '_Hemi'];
    Summary_name = ['Summary_' cell2mat(hemi)];
    if isfield(data, 'DiagnosticData') && isfield(data.DiagnosticData, 'LFPTrendLogs') &&...
       isfield(data.DiagnosticData.LFPTrendLogs, ['HemisphereLocationDef_' cell2mat(hemi)]) &&... % Check if data for this side is not empty
       numel(fieldnames(data.DiagnosticData.LFPTrendLogs.(['HemisphereLocationDef_' cell2mat(hemi)]))) > 1 % ignore files with only one line (this data alway already exist elsewhere)

        if ~isfield(DBS_data.(Patient_name), Summary_name)
            DBS_data.(Patient_name).Info.(Summary_name) = table('Size',[0 10],'VariableTypes',...
                {'double', 'cell', 'double', 'datetime', 'datetime', 'double','cell', 'cell','double','double'},...
                'VariableNames',{'Sensing_Freq', 'Contact_name', 'Sensing_Days','Date_Start','Date_End','Num_LFP_Bins','Amplitude_Range(min,max)','Event_types','Num_Events','Num_Events_with_SNPSHT'});
        end

        % find unique sensing frequencies
        try
            unq_frq = unique(reshape(table2array(DBS_data.(Patient_name).Groups.(hemi_name).Sensing_Freq),[],1),'stable');
            unq_frq = unq_frq(~isnan(unq_frq));
        catch
            unq_frq = [];

        for frq = 1:length(unq_frq)
            % sensing_days
            sensing_days = 0;
            date_start_indx = [];
            date_end_indx = [];
            for i = 1:size(DBS_data.(Patient_name).Groups.(hemi_name).Sensing_Freq,2)
                if any(table2array(DBS_data.(Patient_name).Groups.(hemi_name).Sensing_Freq(:,i))==unq_frq(frq))
                    sensing_days = sensing_days+1;
                    if isempty(date_start_indx) % this is how we find the first date - the first time the if statement is true
                        date_start_indx = i;
                    end
                    date_end_indx = i; % this is how we find the last date - the last time the is statement is true
                end
            end
            all_contacts = reshape(table2array(DBS_data.(Patient_name).Groups.(hemi_name).Contacts),[],1);
            Contacts = all_contacts(reshape(table2array(DBS_data.(Patient_name).Groups.(hemi_name).Sensing_Freq),[],1)  == unq_frq(frq));
            Contact_name = {unique(Contacts)};
            % date_start - easy to find after we saved the index
            date_start = DBS_data.(Patient_name).Groups.(hemi_name).Sensing_Freq.Properties.VariableNames(date_start_indx);
            % date_start - easy to find after we saved the index
            date_end = DBS_data.(Patient_name).Groups.(hemi_name).Sensing_Freq.Properties.VariableNames(date_end_indx);

            % Num_LFP_Bins
            Num_LFP_Bins = length(find(reshape(table2array(DBS_data.(Patient_name).Groups.(hemi_name).Sensing_Freq),[],1) == unq_frq(frq)));

            % Amplitude_Range
            all_amp = reshape(table2array(DBS_data.(Patient_name).TrendLogs.(hemi_name).Amp_table),[],1);
            amp = all_amp(table2array(DBS_data.(Patient_name).Groups.(hemi_name).Sensing_Freq) == unq_frq(frq));
            min_amp = min(min(amp));
            max_amp = max(max(amp));
            Amplitude_Range = {min_amp max_amp};

            if ~isstruct(DBS_data.(Patient_name).Events.Events_data)
                % Event_types
                events_indx = (DBS_data.(Patient_name).Events.Events_data.DateTime<date_end & DBS_data.(Patient_name).Events.Events_data.DateTime>date_start);
                events_sense = DBS_data.(Patient_name).Events.Events_data(events_indx,:);
                num_events = size(unique(events_sense.EventName),1);
                Event_types = cell(1, 5);
                if ischar(events_sense.EventName)
                    Event_types(1) = {events_sense.EventName};
                else
                    Event_types(1:num_events) = unique(events_sense.EventName);
                end

                %Num_Events
                Num_Events = size(events_sense,1);
                Num_Events_with_SNPSHT = sum(events_sense.Snapshot);
            else
                Event_types = cell(1, 4);
                Num_Events = 0;
                Num_Events_with_SNPSHT = 0;
            end

            newRow_data = table(unq_frq(frq),Contact_name ,sensing_days ,datetime(date_start), datetime(date_end), Num_LFP_Bins,...
                Amplitude_Range,Event_types,Num_Events,Num_Events_with_SNPSHT,...
                'VariableNames', {'Sensing_Freq', 'Contact_name', 'Sensing_Days','Date_Start','Date_End','Num_LFP_Bins',...
                'Amplitude_Range(min,max)','Event_types','Num_Events','Num_Events_with_SNPSHT'});

            % Append the new row to the table
            DBS_data.(Patient_name).Info.(Summary_name) = [DBS_data.(Patient_name).Info.(Summary_name); newRow_data];
        end
    end
end

%% Section 3: Device Info
path = data.DeviceInformation.Initial;    
% Initialize struct
deviceInfo = struct();

%
if isfield(path, 'DeviceName')
    deviceInfo.DeviceName = path.DeviceName;
end
%
if isfield(path, 'ImplantDate')
    deviceInfo.ImplantDate = path.ImplantDate;
end
%
if isfield(path, 'NeurostimulatorLocation')
    deviceInfo.NeurostimulatorLocation = path.NeurostimulatorLocation;
end
%
if isfield(data, 'AbnormalEnd')
    deviceInfo.AbnormalEnd = data.AbnormalEnd;
end
%
if isfield(data, 'SessionDate')
    deviceInfo.SessionDate = data.SessionDate;
end
%
if isfield(data, 'ProgrammerTimezone')
    deviceInfo.ProgrammerTimezone = data.ProgrammerTimezone;
end
%
if isfield(data, 'ProgrammerUtcOffset')
    deviceInfo.ProgrammerUtcOffset = data.ProgrammerUtcOffset;
end
%
if isfield(data.BatteryInformation, 'BatteryPercentage')
    deviceInfo.BatteryPercentage = data.BatteryInformation.BatteryPercentage;
end
%
if isfield(data.Stimulation, 'FinalStimStatus')
    deviceInfo.FinalStimStatus = data.Stimulation.FinalStimStatus;
end

DBS_data.(Patient_name).Info.deviceInfo = deviceInfo;


%% Section 4: Patient Info
% Initials
DBS_data.(Patient_name).Info.PatientInfo.Initials = Patient_Initials;

% Gender
if isfield(data.PatientInformation.Initial, 'PatientGender')
    DBS_data.(Patient_name).Info.PatientInfo.Gender = data.PatientInformation.Initial.PatientGender;
end

%%% To be Continued

end