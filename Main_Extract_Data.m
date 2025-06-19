%% ========================================================================
% This is the main file that is used to convert json files 
%  and create a DBS Data Structure.
% follow the steps to make sure all configurations are correct.

%% Option 1: choose save path for the DBS_data Structure
use_current_path = 1;
% use_current_path = 1  -> search for the DBS_data.mat in the current
%                          path or create a new copy if it wasn't found.
% use_current_path = 0  -> choose the folder in which the data file will be 
%                          saved manually

%% Option 2: choose action for existing files
read_again = 1;
% When using an existing DBS_data Structure, this variable is for 
% choosing what to do with files that already were opened and read.
% read_again = 1   ->   process a file again even if it was already processed previously.
% read_again = 0   ->   don't process the file if it was already processed.

%% Option 3: Choose Processing Method
Process_all = 1;
% Process_all = 1   ->  Select a folder - the program will process all JSON
%                       files in the folder and the sub-folders inside it.
% Process_all = 0   ->   Choose One specific JSON file to process.

%% Option 4: Choose Types of data to process
% Control which types of data you would like to extract from the JSON files.
Process_TrendLogs =  1;
Process_Events =     1;
Process_Groups =     1;
Process_Streaming =  1;
Process_Montages =   1;

%% ========================================================================
% Initialization
%% ========================================================================
addpath("Data Orgenization\")    % Folder for organizing and processing data

%% Step 1: Choosing a data path
% Determine whether to use the current working directory or prompt the user to select a folder
if use_current_path
    DBS_data_Directory = pwd;    % Use the current directory as the data folder
else
    % Prompt the user to select a folder that contains the MAT data files
    DBS_data_Directory = uigetdir('', 'Select the folder that contains your MAT data files'); 
end

% Check if the user canceled the folder selection
if DBS_data_Directory == 0
    error('No folder selected. Operation cancelled.');  % Stop execution if no folder was selected
else
    fprintf('%s - was selected as main path: \n', DBS_data_Directory);  % Confirm selected path to the user
end

% Load previously saved data if it exists, otherwise initialize a new structure
dbs_data_file = fullfile(DBS_data_Directory, 'DBS_data.mat');
if exist(dbs_data_file, 'file')
    load(dbs_data_file, 'DBS_data');
else
    DBS_data = struct();
end

%% Step 2: choose action for existing files
if Process_all
    % Clear and informative folder selection
    folder_path = uigetdir('', 'Select the top-level folder that contains the JSON files (including subfolders)');
    
    if folder_path == 0  % Check if the user canceled the operation (uigetdir returns 0 on cancel)
        disp('No folder was selected. Operation canceled.');
        return;
    end

    % Find all JSON files recursively
    json_files = dir(fullfile(folder_path, '**', '*.json'));

    if isempty(json_files)
        warning('No JSON files were found in the selected folder or its subfolders.');
        return;
    end

    % Extract filenames and paths
    filenames = {json_files.name};
    data_pathnames = {json_files.folder};

    % Confirm selection to user
    fprintf('Selected folder: %s\nFound %d JSON files.\n', folder_path, numel(filenames));

else
    % File selection dialog
    [filenames, data_pathnames] = uigetfile('*.json', ...
        'Select one or more JSON files', ...
        'MultiSelect', 'on');
    
    % Handle cancellation
    if isequal(filenames, 0)
        disp('No files were selected. Operation canceled.');
        return;
    end

    % Normalize output
    if ~iscell(filenames)
        filenames = {filenames};
        data_pathnames = {data_pathnames};
    else
        % If multiple files are selected, data_pathnames is a single string
        data_pathnames = repmat({data_pathnames}, size(filenames));
    end

    % Confirm selection
    fprintf('Selected %d JSON file(s):\n', numel(filenames));
    for i = 1:numel(filenames)
        fprintf(' - %s\n', fullfile(data_pathnames{i}, filenames{i}));
    end
end

%% ========================================================================
% The file loop - extracts data from each file
%% ========================================================================

for file = 1:numel(filenames)
    %%%%%%%%%%%%%%% open file and directory %%%%%%%%%%%%%%%
    %
    if iscell(filenames)
        file_name = filenames{file};
    else
        file_name = filenames;
    end
    %
    if iscell(data_pathnames)
        file_path = fullfile(data_pathnames{file}, file_name);
    else
        file_path = fullfile(data_pathnames, file_name);
    end
    %
    try
        data = jsondecode(fileread(file_path));
    catch ME
        warning('Error decoding %s: %s', file_path, ME.message);
        continue;
    end    
    
    %%%%%%%%%%%%%%% create patient folder %%%%%%%%%%%%%%%
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
    file_already_processed = any(strcmp(DBS_data.(Patient_name).Info.Files.Filename, file_name));
    if file_already_processed && read_again == 0
        disp(['File: ' file_name ' already processed. Skipping.']);
        continue
    end
    

    %%%%%%%%%%%%%%% Skip Flags %%%%%%%%%%%%%%%
    skip_LFP = false;
    skip_streaming = false;

    % Check if necessary fields exist and contain data
    if ~isfield(data, 'Groups') || ...
            ~isfield(data, 'DiagnosticData') || ...
            ~isfield(data.DiagnosticData, 'LFPTrendLogs') || ...
            (isempty(data.Groups.Initial) && isempty(data.Groups.Final))

        disp('File format is not compatible for LFP extraction');
        skip_LFP = true;
    end

    % Try to find active group in Initial or Final and check for SensingChannel field
    try
        Group_indx = find([data.Groups.Initial.ActiveGroup]);
        if ~isfield(data.Groups.Initial(Group_indx).ProgramSettings, 'SensingChannel')
            disp('File format is not compatible for stimulation data extraction');
            skip_LFP = true;
        end
    catch
        try
            Group_indx = find([data.Groups.Final.ActiveGroup]);
            if ~isfield(data.Groups.Final(Group_indx).ProgramSettings, 'SensingChannel')
                disp('File format is not compatible for stimulation data extraction');
                skip_LFP = true;
            end
        catch
            disp('File format is not compatible for stimulation data extraction');
            skip_LFP = true;
        end
    end

    % =====================================================================
    %  fill in data 
    % =====================================================================

    fprintf('Processing file %d of %d: %s\n', file, numel(filenames), file_name);
    %% 1. TrendLogs
    if ~skip_LFP && Process_TrendLogs
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
    end

    %% 2. Events
    % organized table with plot in the last column
    if Process_Events
        if ~isfield(DBS_data.(Patient_name), 'Events') ||...
                isempty(DBS_data.(Patient_name).Events.Events_data)
            DBS_data.(Patient_name).Events = struct();
            DBS_data.(Patient_name).Events.Events_data = CreateEvents(data);
        else
            DBS_data.(Patient_name).Events.Events_data = EditEvents(data, DBS_data);
        end
    end

    %% 3. Groups
    if ~skip_LFP && Process_Groups
        if ~isfield(DBS_data.(Patient_name), 'Groups') || isempty(fieldnames(DBS_data.(Patient_name).Groups))
            DBS_data.(Patient_name).Groups = CreateGroups(data);
        else % Add relevat info to existing tables and update statistics
            DBS_data.(Patient_name).Groups = EditGroups(data, DBS_data);
        end
    end

    %% 3. Streaming
    if ~skip_streaming && Process_Streaming
        if ~isfield(DBS_data.(Patient_name), 'Streaming') || isempty(fieldnames(DBS_data.(Patient_name).Streaming))
            DBS_data.(Patient_name).Streaming = CreateStreaming(data, file_name);
        else % Add relevat info to existing tables and update statistics
            DBS_data.(Patient_name).Streaming = EditStreaming(data, file_name, DBS_data);
        end
    end

    %% 4. Montages (processed as event)
    if Process_Montages
        DBS_data.(Patient_name).Events.Events_data = Extract_Montage(data, DBS_data);
    end

    %% 5. Update Info
    DBS_data = UpdateInfo(data, DBS_data, file_name, file);

end

% Save the updated main struct
fprintf('Saving DBS_data\n');
save(dbs_data_file, 'DBS_data');
fprintf('Saved updated DBS_data to %s\n', dbs_data_file);




















