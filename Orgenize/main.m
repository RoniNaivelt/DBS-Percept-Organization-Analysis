%% For opening single \ specific json 

  % [filenames, data_pathnames] = uigetfile('*.json', 'MultiSelect', 'on');
  % cd(data_pathnames)


 % data = jsondecode(fileread(filenames));
% data1 = jsondecode(fileread(filenames{1}));
% data2 = jsondecode(fileread(filenames{2}));
% data3 = jsondecode(fileread(filenames{3}));
% data4 = jsondecode(fileread(filenames{4}));
% data5 = jsondecode(fileread(filenames{5}));
% data6 = jsondecode(fileread(filenames{6}));
% data7 = jsondecode(fileread(filenames{7}));
% data8 = jsondecode(fileread(filenames{8}));
% data9 = jsondecode(fileread(filenames{9}));
% data10 = jsondecode(fileread(filenames{10}));



%% This is the main file - follow the steps to make sure all configuration are correct
%1. Initialize or load the main structure
% choose the folder in which the data file will be saved (or where it already is)

%DBS_data_Directory = uigetdir('Select a folder containing the MAT file'); % option 1 - choose file manually  
DBS_data_Directory = pwd; % option 2 - use current folder

if exist('DBS_data.mat', 'file')
    load('DBS_data.mat', 'DBS_data');
else
    DBS_data = struct();
end

% 2. choose action for existing files
read_again = 1;
    % this variable is for choosing what to do with files that already were opened and read:
    % read_again = 1   ->   process a file again even if it was already processed previously.
    % read_again = 0   ->   don't process the file if it was already processed.

% 3. Open JSON files to process
    % unmark one of the following options:
% %%% Option 1 - choose files manualy:
% [filenames, data_pathnames] = uigetfile('*.json', 'MultiSelect', 'on');
% if ~iscell(filenames)
%     filenames = {filenames};
% end

% %%% Option 2 (recomended) - choose a folder, all JSON files inside will be processed.
folder_path = uigetdir('Select a folder containing JSON files');
if folder_path == 1 % Check if the user canceled the operation
    disp('no files were chosen, operation canceled');
    return;
end
json_files = dir(fullfile(folder_path, '**', '*.json'));
filenames = {json_files.name};
data_pathnames = {json_files.folder};


% 4. The file loop - extracts data from each file:
    % note the 'fill in data' selction inside of the loop - mark/unmark the
    % subsection (1,2,3..) to choose what structure you want to extract.
    % you can start with one type of data and edd another later, the
    % DBS_data file will be updated, it will not over-write.
    
for file = 1:numel(filenames)
    %%%%%%%%%%%%%%% open file and directory %%%%%%%%%%%%%%%
    if iscell(data_pathnames)
        cd(data_pathnames{file});
    else
        cd(data_pathnames);
    end
    %
    if iscell(filenames)
        file_name = filenames{file};
    else
        file_name = filenames;
    end
    %
    data = jsondecode(fileread(file_name));
    
    
    %%%%%%%%%%%%%%% create patient folder %%%%%%%%%%%%%%%
    % get patient info
    Patiant_Initials = [data.PatientInformation.Initial.PatientFirstName(1) data.PatientInformation.Initial.PatientLastName(1)];
    Patient_name = ['Patient_' Patiant_Initials];
    % If the patient is new, create a patient subfolder
    if ~isfield(DBS_data, Patient_name)
        DBS_data.(Patient_name) = struct();
    end
    
    
    %%%%%%%%%%%%%%% Patient's Info %%%%%%%%%%%%%%%
    % if no info substruct -  create one:
    if ~isfield(DBS_data.(Patient_name), 'Info')
        DBS_data.(Patient_name).Info = struct();
    end
    % if no file history substruct -  create one:
    if ~isfield(DBS_data.(Patient_name).Info, 'Files')
    DBS_data.(Patient_name).Info.Files = table('Size',[0 5],'VariableTypes',...
        {'string', 'datetime', 'logical', 'logical', 'logical'},...
        'VariableNames',{'Filename','dDate','Trendlogs','Events','SnapshotEvents'});
    end

    
    %%%%%%%%%%%%%%% check file %%%%%%%%%%%%%%%
    % If file was already processed before AND read again is 0 - skip file.
    % If file was not processed before OR read again is 1 - start processing file.
    if (isempty(DBS_data.(Patient_name).Info) && isfield(DBS_data.(Patient_name).Info,'Files') &&...
            any(strcmp(DBS_data.(Patient_name).Info.Files.Filename, file_name))) && read_again == 0
        disp(['file: ' file_name ' aready processed']);
        continue
    elseif (~isempty(DBS_data.(Patient_name).Info) && isfield(DBS_data.(Patient_name).Info,'Files')...
            && ~any(strcmp(DBS_data.(Patient_name).Info.Files.Filename, file_name))) || read_again == 1
        disp(['Processing file ' num2str(file) '/' num2str(numel(filenames)) ':' file_name]);
        
        %%%%%%%%%%%%%%% Skip Flags %%%%%%%%%%%%%%%
        if ~isfield(data, 'Groups') ||...
                ~isfield(data, 'DiagnosticData') ||...
                ~isfield(data.DiagnosticData, 'LFPTrendLogs') ||...
                (isempty(data.Groups.Initial) && isempty(data.Groups.Final))
           
            cd(DBS_data_Directory)
            DBS_data.(Patient_name).Events.Events_data = Extract_Montage(data, DBS_data);
            disp('File format is not compatible')
            continue
        end
        
        try
            Group_indx = find([data.Groups.Initial.ActiveGroup]);
            if ~isfield(data.Groups.Initial(Group_indx).ProgramSettings, 'SensingChannel')
                disp('File format is not compatible (problem in Groups)')
                continue
            end
        catch
            Group_indx = find([data.Groups.Final.ActiveGroup]);
            if ~isfield(data.Groups.Final(Group_indx).ProgramSettings, 'SensingChannel')
                disp('File format is not compatible (problem in Groups)')
                continue
            end
        end


        
        %%%%%%%%%%%%%%% fill in data %%%%%%%%%%%%%%%
        cd(DBS_data_Directory)
        %% 1. TrendLogs
        % create tables and save the following info:
        % 1. LFP table
        % 2. Amplitude table
        % 3. statistics table
        % 4. Time vector
        if ~isfield(DBS_data.(Patient_name), 'TrendLogs') || isempty(fieldnames(DBS_data.(Patient_name).TrendLogs))
            DBS_data.(Patient_name).TrendLogs = CreateTrendLogs(data);
        else % Add relevat info to existing tables and update statistics
            DBS_data.(Patient_name).TrendLogs = EditTrendLogs(data, DBS_data);
        end
        
        
        %% 2. Events
        % organized table with plot in the last column
        if ~isfield(DBS_data.(Patient_name), 'Events') ||...
                isempty(DBS_data.(Patient_name).Events.Events_data)
            DBS_data.(Patient_name).Events = struct();
            DBS_data.(Patient_name).Events.Events_data = CreateEvents(data);
        else
            DBS_data.(Patient_name).Events.Events_data = EditEvents(data, DBS_data);
        end

        %% 3. Groups
        if ~isfield(DBS_data.(Patient_name), 'Groups') || isempty(fieldnames(DBS_data.(Patient_name).Groups))
            DBS_data.(Patient_name).Groups = CreateGroups(data);
        else % Add relevat info to existing tables and update statistics
            DBS_data.(Patient_name).Groups = EditGroups(data, DBS_data);
        end
        
        %% 4. Montages (as contol event)
        DBS_data.(Patient_name).Events.Events_data = Extract_Montage(data, DBS_data);
        
        %% 5. Update Info
        DBS_data = UpdateInfo(data, DBS_data, file_name, file);
    end
end

% Save the updated main struct
disp('Finished loading files. Saving DBS_Data...')
cd(DBS_data_Directory)
save('DBS_data.mat', 'DBS_data');




















