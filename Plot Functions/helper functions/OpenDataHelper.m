function data_out = OpenDataHelper(DBS_data, Patient_names, Hemi, day_range, sensing_band)
% Universal function to extract LFP and related DBS data from the DBS_data struct.
%
% INPUTS:
%   DBS_data          : Struct loaded from DBS .mat file
%   Patient_names     : String (one patient) or cell array of patient names
%   Hemi              : 'Left' or 'Right'
%   periodogram_days  : (Optional) Nx2 array with [start_day, end_day] per patient
%   sensing_band      : (Optional) array indicating a band of sensing
%                       frequencies such as [Min_frequency Max_frequency]
% OUTPUT:
%   data_out : struct with fields for each patient, containing:
%              - LFP_vec, Amp_vec, sensing_vec, etc.
%              - Daily/Nightly means if periodogram_days provided
%              - Trimmed vectors with NaNs removed
%              - Theta/Beta-filtered LFP matrices

if ischar(Patient_names)
    Patient_names = {Patient_names};
end

if nargin < 5
    sensing_band = [0 inf];
end
if nargin < 4
    day_range = [];
end

for i_pat = 1:length(Patient_names)
    Patient_name = Patient_names{i_pat};
    hemi_key = [Hemi '_Hemi'];
    key = [Patient_name '_' Hemi];

    % Surgery date
    surg_date = DBS_data.(Patient_name).Info.deviceInfo.ImplantDate;
    surg_date = datetime(regexprep(surg_date(1:end-1), 'T', ' '));

    % Extract core tables
    LFP_T = DBS_data.(Patient_name).TrendLogs.(hemi_key).LFP_table;
    Amp_T = DBS_data.(Patient_name).TrendLogs.(hemi_key).Amp_table;
    sensing_T = DBS_data.(Patient_name).Groups.(hemi_key).Sensing_Freq;
    contacts_T = DBS_data.(Patient_name).Groups.(hemi_key).Contacts;
    Stimrate_T = DBS_data.(Patient_name).Groups.(hemi_key).Stim_Rate;
    PulseWidth_T = DBS_data.(Patient_name).Groups.(hemi_key).Pulse_width;
    Events = DBS_data.(Patient_name).Events.Events_data;
    Impedance_T = DBS_data.(Patient_name).Groups.(hemi_key).Impedance;
    TEED_T = (Amp_T.^2 .* Impedance_T .* PulseWidth_T .* Stimrate_T).*10e-6;

    % Trim vectors
    if nargin > 3 && ~isempty(day_range)
        if iscell(day_range)
            day_range = day_range{i_pat};
        else
            day_range = day_range;
        end
        LFP_T = LFP_T(:, day_range(1):day_range(2));
        Amp_T = Amp_T(:, day_range(1):day_range(2));
        sensing_T = sensing_T(:, day_range(1):day_range(2));
        contacts_T = contacts_T(:, day_range(1):day_range(2));
        Stimrate_T = Stimrate_T(:, day_range(1):day_range(2));
        PulseWidth_T = PulseWidth_T(:, day_range(1):day_range(2));
        TEED_T = TEED_T(:, day_range(1):day_range(2));
        % Trim Events
        event_day_offsets = days(Events.DateTime - surg_date);  % Use same reference start time
        event_trimmed_idx = event_day_offsets >= day_range(1) & event_day_offsets <= day_range(2);
        Events = Events(event_trimmed_idx, :);

    end


    % Trimmed Datetime vec
    col_names = string(LFP_T.Properties.VariableNames);
    time_strings = string(LFP_T.Properties.RowNames);
    Datetime_vec = NaT(numel(col_names)*numel(time_strings),1);
    dates = datetime(col_names, 'InputFormat', 'dd-MMM-yyyy');
    idx = 1;
    for d = 1:length(dates)
        date_only = dates(d);
        for t = 1:length(time_strings)
            dt = datetime([char(date_only) ' ' char(time_strings(t))], ...
                'InputFormat', 'dd-MMM-yyyy HH:mm', ...
                'Format', 'dd-MMM-yyyy HH:mm');
            Datetime_vec(idx, 1) = dt;
            idx = idx + 1;
        end
    end

    % Raw Tables (only after trimming by days)
    data_out.(key).Tables.LFP = LFP_T;
    data_out.(key).Tables.Amp = Amp_T;
    data_out.(key).Tables.Sensing = sensing_T;
    data_out.(key).Tables.Stimrate = Stimrate_T;
    data_out.(key).Tables.PulseWidth = PulseWidth_T;
    data_out.(key).Tables.TEED = TEED_T;

    % Reshape and normalize
    LFP_T = normalize(LFP_T(:,:), 'range');
    LFP_vec = reshape(table2array(LFP_T), [], 1);
    Amp_vec = reshape(table2array(Amp_T), [], 1);
    sensing_vec = reshape(table2array(sensing_T), [], 1);
    contacts_vec = reshape(table2array(contacts_T), [], 1);
    Stimrate_vec = reshape(table2array(Stimrate_T), [], 1);
    PulseWidth_vec = reshape(table2array(PulseWidth_T), [], 1);
    TEED_vec = reshape(table2array(TEED_T), [], 1);
    days_since_surg = days(Datetime_vec - surg_date);

    % Remove NaNs (trimmed vectors)
    nanIdx = isnan(LFP_vec);

    % Trimmed versions
    LFP_vec_trimmed = LFP_vec(~nanIdx);
    Amp_vec_trimmed = Amp_vec(~nanIdx);
    sensing_vec_trimmed = sensing_vec(~nanIdx);
    contacts_vec_trimmed = contacts_vec(~nanIdx);
    Datetime_vec_trimmed = Datetime_vec(~nanIdx);
    days_since_surg_trimmed = days_since_surg(~nanIdx);
    Stimrate_vec_trimmed = Stimrate_vec(~nanIdx);
    PulseWidth_vec_trimmed = PulseWidth_vec(~nanIdx);
    TEED_vec_trimmed = TEED_vec(~nanIdx);

    % % 
    % % % Setup output struct
    % % data_out.(key).LFP_vec = LFP_vec;
    % % data_out.(key).Amp_vec = Amp_vec;
    % % data_out.(key).sensing_vec = sensing_vec;
    % % data_out.(key).contacts_vec = contacts_vec;
    % % data_out.(key).Datetime_vec = Datetime_vec;

    % Add trimmed vectors
    data_out.(key).LFP_vec_trimmed = LFP_vec_trimmed;
    data_out.(key).Amp_vec_trimmed = Amp_vec_trimmed;
    data_out.(key).sensing_vec_trimmed = sensing_vec_trimmed;
    data_out.(key).contacts_vec_trimmed = contacts_vec_trimmed;
    data_out.(key).Datetime_vec_trimmed = Datetime_vec_trimmed;
    data_out.(key).Stimrate_vec_trimmed = Stimrate_vec_trimmed;
    data_out.(key).PulseWidth_vec_trimmed = PulseWidth_vec_trimmed;
    data_out.(key).TEED_vec_trimmed = TEED_vec_trimmed;
    data_out.(key).days_since_surg = days_since_surg_trimmed;
    data_out.(key).Events = Events;
    data_out.(key).Events_Days_Since_Surge = days(Events.DateTime - surg_date);

    %% Event AUC extraction for Beta, Theta, Alpha, Gamma
    hemispheres = {'Left', 'Right'};
    i_hemi = find(strcmpi(Hemi, hemispheres));

    % Initialize
    Event_dates = datetime.empty;
    Event_names = {};
    Beta_AUC = [];
    Theta_AUC = [];
    Alpha_AUC = [];
    Gamma_AUC = [];

    for i = 1:size(Events,1)
        hemi_field = [hemispheres{i_hemi} '_Snapshot'];
        if (Events.Snapshot(i) == 1) && ...
            isfield(Events.(hemi_field)(i), 'Beta_data') && ...
            ~isempty(Events.(hemi_field)(i).Beta_data)

            Event_dates(end+1) = Events.DateTime(i);
            Event_names{end+1} = Events.EventName{i};

            Beta_AUC(end+1)  = Events.(hemi_field)(i).Beta_data.Beta_AUC_relative;
            Theta_AUC(end+1) = Events.(hemi_field)(i).Theta_data.Theta_AUC_relative;
            Alpha_AUC(end+1) = Events.(hemi_field)(i).Alpha_data.Alpha_AUC_relative;
            Gamma_AUC(end+1) = Events.(hemi_field)(i).Gamma_data.Gamma_AUC_relative;
        end
    end

    Event_dates_since_surg = days(Event_dates - surg_date);

    % Save into output struct
    data_out.(key).Event_Summary.Event_dates = Event_dates;
    data_out.(key).Event_Summary.Event_names = Event_names;
    data_out.(key).Event_Summary.Event_dates_since_surg = Event_dates_since_surg;
    data_out.(key).Event_Summary.Beta_AUC = Beta_AUC;
    data_out.(key).Event_Summary.Theta_AUC = Theta_AUC;
    data_out.(key).Event_Summary.Alpha_AUC = Alpha_AUC;
    data_out.(key).Event_Summary.Gamma_AUC = Gamma_AUC;

    %% Day/Night mean extraction

    % Time indices
    day_idx = 8*6+1 : 20*6;       % 08:00 to 20:00
    night_idx = [1 : 6*6, 20*6+1 : 24*6];  % 00:00–06:00 and 20:00–24:00

    % Day mean
    data_out.(key).LFP_day = table2array(LFP_T(day_idx, :));
    data_out.(key).Amp_day = table2array(Amp_T(day_idx, :));
    data_out.(key).sensing_day = table2array(sensing_T(day_idx, :));
    data_out.(key).contscts_day = table2array(contacts_T(day_idx, :));
    data_out.(key).Stimrate_day = table2array(Stimrate_T(day_idx, :));
    data_out.(key).PulseWidth_day = table2array(PulseWidth_T(day_idx, :));
    data_out.(key).TEED_day = table2array(TEED_T(day_idx, :));

    % Night mean
    data_out.(key).LFP_night = table2array(LFP_T(night_idx, :));
    data_out.(key).Amp_night = table2array(Amp_T(night_idx, :));
    data_out.(key).sensing_night = table2array(sensing_T(night_idx, :));
    data_out.(key).contscts_night = table2array(contacts_T(night_idx, :));
    data_out.(key).Stimrate_night = table2array(Stimrate_T(night_idx, :));
    data_out.(key).PulseWidth_night = table2array(PulseWidth_T(night_idx, :));
    data_out.(key).TEED_night = table2array(TEED_T(night_idx, :));

    data_out.(key).days_vec = datetime(LFP_T.Properties.VariableNames);
    data_out.(key).days_since_surg_vec = days(data_out.(key).days_vec - surg_date);


    %% Theta & Beta LFP Matrices
    band_LFP_mat = NaN(size(LFP_T));
    sensing_array = table2array(sensing_T);
    lfp_array = table2array(LFP_T);

    if ~isempty(sensing_band)
        band_mask = sensing_array > sensing_band(1) & sensing_array < sensing_band(2);
    else
        band_mask = sensing_array > 0 & sensing_array < inf;
    end

    band_LFP_mat(band_mask) = lfp_array(band_mask);
    data_out.(key).band_LFP_mat = band_LFP_mat;

    %% Day and Night band LFP
    band_LFP_mat_day = NaN(size(data_out.(key).LFP_day));
    band_LFP_mat_night = NaN(size(data_out.(key).LFP_night));

    % Compute day/night masks
    band_mask_day = data_out.(key).sensing_day > sensing_band(1) & ...
        data_out.(key).sensing_day < sensing_band(2);
    band_mask_night = data_out.(key).sensing_night > sensing_band(1) & ...
        data_out.(key).sensing_night < sensing_band(2);

    % Apply masks
    band_LFP_mat_day(band_mask_day) = data_out.(key).LFP_day(band_mask_day);
    band_LFP_mat_night(band_mask_night) = data_out.(key).LFP_night(band_mask_night);

    data_out.(key).band_LFP_mat_day = band_LFP_mat_day;
    data_out.(key).band_LFP_mat_night = band_LFP_mat_night;

end
end
