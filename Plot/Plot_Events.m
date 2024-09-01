close all
Patient_name = 'Patient_EC';
hemisphers = {'Left', 'Right'};

% save surgery date
surg_date = DBS_data.(Patient_name).Info.deviceInfo.ImplantDate;
surg_date = datetime(regexprep(surg_date(1:end-1),'T',' '));
            
for hemi_idx = 1:numel(hemisphers)
    LFP_T = DBS_data.(Patient_name).TrendLogs.([hemisphers{hemi_idx} '_Hemi']).LFP_table;
    LFP_vec = reshape(table2array(LFP_T),[],1);
    Amp_T = DBS_data.(Patient_name).TrendLogs.([hemisphers{hemi_idx} '_Hemi']).Amp_table;
    Amp_vec = reshape(table2array(Amp_T),[],1);
    sensing_T = DBS_data.(Patient_name).Groups.([hemisphers{hemi_idx} '_Hemi']).Sensing_Freq;
    sensing_vec = reshape(table2array(sensing_T),[],1);
    contacts_T = DBS_data.(Patient_name).Groups.([hemisphers{hemi_idx} '_Hemi']).Contacts;
    contacts_vec = reshape(table2array(contacts_T),[],1);
    Stimrate_T = DBS_data.(Patient_name).Groups.([hemisphers{hemi_idx} '_Hemi']).Stim_Rate;
    Stimrate_vec = reshape(table2array(Stimrate_T),[],1);
    PulseWidth_T = DBS_data.(Patient_name).Groups.([hemisphers{hemi_idx} '_Hemi']).Pulse_width;
    PulseWidth_vec = reshape(table2array(PulseWidth_T),[],1);
    Impedance_T = DBS_data.(Patient_name).Groups.([hemisphers{hemi_idx} '_Hemi']).Impedance;
    Impedance_vec = reshape(table2array(Impedance_T),[],1);
    Datetime_vec = DBS_data.(Patient_name).TrendLogs.([hemisphers{hemi_idx} '_Hemi']).Time_vector;
    days_since_surge_vec = days(Datetime_vec - surg_date);
    TEED_T = (Amp_T.^2 .* Impedance_T .* PulseWidth_T .* Stimrate_T).*10e-6;
    TEED_vec = (Amp_vec.^2 .* Impedance_vec .* PulseWidth_vec .* Stimrate_vec).*10e-6;
    %
    cols_names = LFP_T.Properties.VariableNames;
    rows_names = LFP_T.Properties.RowNames;
    % nans
    nanIdx = isnan(LFP_vec);
    LFP_vec = LFP_vec(~nanIdx);
    Amp_vec = Amp_vec(~nanIdx);
    sensing_vec = sensing_vec(~nanIdx);
    contacts_vec = contacts_vec(~nanIdx);
    if ~strcmp(Patient_name, 'Patient_AD')
        Datetime_vec = Datetime_vec(~nanIdx);
        days_since_surge_vec = days_since_surge_vec(~nanIdx);

    end

    % Events - Snapshots
    Events = DBS_data.(Patient_name).Events.Events_data;
    Event_dates = datetime.empty;
    Event_names = {};
    Beta_AUC = [];
    Theta_AUC = [];
    Alpha_AUC = [];
    Gamma_AUC = [];

    for i = 1:size(Events,1)
        if  (Events.Snapshot(i) == 1) && (~isempty(Events.([hemisphers{hemi_idx} '_Snapshot'])(i).Beta_data))  
            Event_dates(end+1) = Events.DateTime(i);
            Event_names{end+1} = Events.EventName{i};
            Beta_AUC(end+1) = Events.([hemisphers{hemi_idx} '_Snapshot'])(i).Beta_data.Beta_AUC_relative;
            Theta_AUC(end+1) = Events.([hemisphers{hemi_idx} '_Snapshot'])(i).Theta_data.Theta_AUC_relative;
            Alpha_AUC(end+1) = Events.([hemisphers{hemi_idx} '_Snapshot'])(i).Alpha_data.Alpha_AUC_relative;
            Gamma_AUC(end+1) = Events.([hemisphers{hemi_idx} '_Snapshot'])(i).Gamma_data.Gamma_AUC_relative;
        end
    end
    Event_dates_since_surg = days(Event_dates - surg_date);
    % for j = 1:size(Event_names,2)
    %     color_indx = find(strcmp(event_types,Event_names(j)));
    %     xline(Events.DateTime(j), '-','Color',colors{color_indx});
    % end

    % find unique periods
    combined_strings = strcat(num2str(sensing_vec(:)), '_', contacts_vec(:));
    [unique_pairs, ~, pair_idx] = unique(combined_strings,'stable');


    %% only Beta
    Amp_mean = table2array(mean(Amp_T,1));
    TEED_mean = table2array(mean(TEED_T,1));
    sensing_mean = table2array(mean(sensing_T,1));
    Stimrate_mean = table2array(mean(Stimrate_T,1));
    PulseWidth_mean = table2array(mean(PulseWidth_T,1));
    days_vec = datetime(LFP_T.Properties.VariableNames);
    days_since_surge = days(days_vec - surg_date);

    event_types = unique(Events.EventName, 'stable')';

    figure(1)
    subplot(2,1,hemi_idx)
    yyaxis right
    plot(days_since_surge,TEED_mean,'--o', 'LineWidth',2, 'DisplayName', 'TEED');
    ylabel('TEED [E/sec]')
    xlabel('Days Since Surgery')

    hold on
    yyaxis left
    colors = jet(numel(event_types));
    for i = 1:length(event_types)
        event_ind = strcmp(Event_names,event_types(i));
        plot(Event_dates_since_surg(event_ind), Beta_AUC(event_ind), '*', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', event_types{i});
    end
    ylabel('Relative AUC')
    grid on
    lgd = legend('Location','northeastoutside');
    set(lgd, 'Color', 'none');
    title([hemisphers{hemi_idx} ' Beta AUC'], 'FontSize',14, 'Interpreter', 'none');


    if numel(unique_pairs) > 1
        yLimits = ylim;
        clr = jet(length(unique_pairs));
        for i = 1:length(unique_pairs)
           current_pair = unique_pairs(i);
           pair_parts = split(unique_pairs(i), '_');
           current_freq = pair_parts(1);
           current_contact1 = pair_parts(2);
           current_contact2 = pair_parts(4);
           LFP_sense =  LFP_vec(pair_idx == i);
           sensing_datetime = days_since_surge_vec(pair_idx == i);
           date_start = sensing_datetime(1);
           date_end = sensing_datetime(end);
           patch([date_start date_start date_end date_end], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], clr(i,:), 'FaceAlpha', 0.2, 'DisplayName', ['SF - ' char(current_freq) ' [Hz] ; Contacts ' char(current_contact1) ' and ' char(current_contact2)])
        end
    else
        title([hemisphers{hemi_idx} ' Beta AUC ; ' char(unique_pairs)])
    end
    set(gcf, 'Color', 'w')




    % subplot(6,1,sub_idx_2)
    % yyaxis left
    % plot(days_vec,Stimrate_mean,'-', 'LineWidth',4)
    % yyaxis right
    % plot(days_vec,PulseWidth_mean,'-', 'LineWidth',4)
    % lgd = legend('Stim rate', 'Pulse Width','Location','northeastoutside');
    % set(lgd, 'Color', 'none');
    % grid on

   %% only Theta
    Amp_mean = table2array(mean(Amp_T,1));
    TEED_mean = table2array(mean(TEED_T,1));
    sensing_mean = table2array(mean(sensing_T,1));
    Stimrate_mean = table2array(mean(Stimrate_T,1));
    PulseWidth_mean = table2array(mean(PulseWidth_T,1));
    days_vec = datetime(LFP_T.Properties.VariableNames);
    days_since_surge = days(days_vec - surg_date);


    event_types = unique(Events.EventName, 'stable')';

    figure(2)
    subplot(2,1,hemi_idx)

    yyaxis right
    plot(days_since_surge,TEED_mean,'--o', 'LineWidth',2, 'DisplayName', 'TEED');
    ylabel('TEED [E/sec]')
    xlabel('Days Since Surgery')
    hold on

    yyaxis left
    colors = jet(numel(event_types));
    for i = 1:length(event_types)
        event_ind = strcmp(Event_names,event_types(i));
        plot(Event_dates_since_surg(event_ind), Theta_AUC(event_ind), '*', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', event_types{i});
    end
    ylabel('Relative AUC')

    grid on
    lgd = legend('Location','northeastoutside');
    set(lgd, 'Color', 'none');
    title([hemisphers{hemi_idx} ' Theta AUC'], 'FontSize',14, 'Interpreter', 'none');


    if numel(unique_pairs) > 1
        yLimits = ylim;
        clr = jet(length(unique_pairs));
        for i = 1:length(unique_pairs)
           current_pair = unique_pairs(i);
           pair_parts = split(unique_pairs(i), '_');
           current_freq = pair_parts(1);
           current_contact1 = pair_parts(2);
           current_contact2 = pair_parts(4);
           LFP_sense =  LFP_vec(pair_idx == i);
           sensing_datetime = days_since_surge_vec(pair_idx == i);
           date_start = sensing_datetime(1);
           date_end = sensing_datetime(end);
           patch([date_start date_start date_end date_end], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], clr(i,:), 'FaceAlpha', 0.2, 'DisplayName', ['SF - ' char(current_freq) ' [Hz] ; Contacts ' char(current_contact1) ' and ' char(current_contact2)])
        end
    else
        title([hemisphers{hemi_idx} ' Theta AUC ; ' char(unique_pairs)])
    end
    set(gcf, 'Color', 'w')


    % subplot(6,1,sub_idx_2)
    % yyaxis left
    % plot(days_vec,Stimrate_mean,'-o', 'LineWidth',4)
    % yyaxis right
    % plot(days_vec,PulseWidth_mean,'-', 'LineWidth',4)
    % lgd = legend('Stim rate', 'Pulse Width','Location','northeastoutside');
    % set(lgd, 'Color', 'none');
    % grid on

    %% All bands
    Amp_mean = table2array(mean(Amp_T,1));
    TEED_mean = table2array(mean(TEED_T,1));
    sensing_mean = table2array(mean(sensing_T,1));
    Stimrate_mean = table2array(mean(Stimrate_T,1));
    PulseWidth_mean = table2array(mean(PulseWidth_T,1));

    %
    days_vec = datetime(LFP_T.Properties.VariableNames);
    days_since_surge = days(days_vec - surg_date);

    figure(3)
    subplot(2,1,hemi_idx)
    yyaxis right
    plot(days_since_surge,TEED_mean,'--o', 'LineWidth',2, 'DisplayName', 'TEED')
    ylabel('TEED [E/sec]')
    xlabel('Days Since Surgery')

    hold on
    yyaxis left
    colors = jet(4);
    plot(Event_dates_since_surg,Beta_AUC, '*', 'Color', colors(1,:), 'LineWidth',2, 'DisplayName', 'Beta')
    plot(Event_dates_since_surg,Theta_AUC,'*', 'Color', colors(2,:), 'LineWidth',2, 'DisplayName', 'Theta')
    plot(Event_dates_since_surg,Alpha_AUC,'*', 'Color', colors(3,:), 'LineWidth',2, 'DisplayName', 'Alpha')
    plot(Event_dates_since_surg,Gamma_AUC,'*', 'Color', colors(4,:), 'LineWidth',2, 'DisplayName', 'Gamma')
    ylabel('Relative AUC')

    grid on
    lgd = legend('Location','northeastoutside');
    set(lgd, 'Color', 'none');
    title([hemisphers{hemi_idx} ' Events'], 'FontSize',14, 'Interpreter', 'none');
    
    % draw patches:
    if numel(unique_pairs) > 1
        yLimits = ylim;
        clr = jet(length(unique_pairs));
        for i = 1:length(unique_pairs)
           current_pair = unique_pairs(i);
           pair_parts = split(unique_pairs(i), '_');
           current_freq = pair_parts(1);
           current_contact1 = pair_parts(2);
           current_contact2 = pair_parts(4);
           LFP_sense =  LFP_vec(pair_idx == i);
           sensing_datetime = days_since_surge_vec(pair_idx == i);
           date_start = sensing_datetime(1);
           date_end = sensing_datetime(end);
           patch([date_start date_start date_end date_end], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], clr(i,:), 'FaceAlpha', 0.2, 'DisplayName', ['SF - ' char(current_freq) ' [Hz] ; Contacts ' char(current_contact1) ' and ' char(current_contact2)])
        end
    else
        title([hemisphers{hemi_idx} ' Events ; ' char(unique_pairs)])
    end

    % subplot(6,1,sub_idx_2)
    % yyaxis left
    % plot(days_vec,Stimrate_mean,'-', 'LineWidth',4)
    % yyaxis right
    % plot(days_vec,PulseWidth_mean,'-', 'LineWidth',4)
    % lgd = legend('Stim rate', 'Pulse Width','Location','northeastoutside');
    % set(lgd, 'Color', 'none');
    % grid on
    set(gcf, 'Color', 'w')


    %% PSD FIT + Plot

    Event_dates = datetime.empty;
    Event_names = {};
    Event_FFT = {};
    %
    Event_fit = {};
    Event_Peak_freq = {};
    Event_Peak_bw = {};
    Event_peak_height = {};
    aperiodic_params = {};

    for i = 1:size(Events,1)
        if  (Events.Snapshot(i) == 1) && (~isempty(Events.([hemisphers{hemi_idx} '_Snapshot'])(i).FFT_Val))  
            Event_dates(end+1) = Events.DateTime(i);
            Event_names{end+1} = Events.EventName{i};

            FFT = Events.([hemisphers{hemi_idx} '_Snapshot'])(i).FFT_Val;
            Event_FFT{end+1} = FFT;
            Freq = Events.([hemisphers{hemi_idx} '_Snapshot'])(i).FFT_Freq;
            [Event_FFT{end+1}, aperiodic_params{end+1}, Event_Peak_freq{end+1}, Event_Peak_bw{end+1}, Event_peak_height{end+1}, rmse] = PSD_Fit(FFT(6:end), Freq(6:end), 1.5, 1, false);                
        end
    end

    Event_dates_since_surg = days(Event_dates - surg_date);
    % Prepare data for the bubble chart
    x_data = [];
    y_data = [];
    bubble_colors = [];
    bubble_sizes = [];
    aperiodic_exponent = [];
    aperiodic_offset = [];

    for i = 1:length(Event_dates_since_surg)
        num_peaks = length(Event_Peak_freq{i});
        x_data = [x_data; repmat(Event_dates_since_surg(i), num_peaks, 1)];
        y_data = [y_data; Event_Peak_freq{i}'];
        bubble_colors = [bubble_colors; Event_peak_height{i}'];
        bubble_sizes = [bubble_sizes; Event_Peak_bw{i}'];
        aperiodic_exponent = [aperiodic_exponent ; repmat(aperiodic_params{i}.A, num_peaks, 1)];
        aperiodic_offset = [aperiodic_offset ; repmat(aperiodic_params{i}.B, num_peaks, 1)];
    end

    bubble_colors = normalize(bubble_colors);
    % Scale the bubble sizes for better visualization
    min_bw = min(bubble_sizes);
    max_bw = max(bubble_sizes);
    scaled_bubble_sizes = 350 * (bubble_sizes - min_bw) / (max_bw - min_bw) + 1; % Adjust scaling factors as needed

    % Plot the bubble chart with sizes and colors
    figure(4)
    subplot(1,2,3-hemi_idx)
    scatter(x_data, y_data, scaled_bubble_sizes, bubble_colors, 'filled');
    cmap = colormap(jet); % Choose a colormap, 'jet' is used here
    c = colorbar; % Show the colorbar to indicate bandwidth values
    clim([0 0.9])
    ylim([4 30])
    ylabel(c ,'Peak Amplitude', 'FontSize', 12)
    xlabel('Days Since Surgery')
    ylabel('Frequency [Hz]');
    title({[hemisphers{hemi_idx} ' Hemisphere Peaks of ' Patient_name], 'Color - Peak Amplitude ; Radius - Peak Bandwith'}, 'FontSize',14, 'Interpreter', 'none');
    set(gca, 'FontSize', 12)

    set(gcf,'color', 'w')
    grid minor
    set(gca, 'XMinorGrid', 'on')
    box off

    savefig(['Peaks Plot - ' Patient_name '.fig'])


    %% Calc error and plot The fit aperiodic components
 
    % Function to calculate R^2 error
    calc_r_squared = @(y, y_fit) 1 - (sum((y - y_fit).^2) / sum((y - mean(y)).^2));

    figure(5)
    % Trend line for aperiodic exponent
    p_exponent = polyfit(x_data, aperiodic_exponent, 1);
    exponent_fit = polyval(p_exponent, x_data);
    r_squared_exponent = calc_r_squared(aperiodic_exponent, exponent_fit);

    % Plot aperiodic exponent
    subplot(2, 2, 3 - hemi_idx);
    plot(x_data, aperiodic_exponent, '.', 'MarkerSize', 15);
    hold on;
    plot(x_data, exponent_fit, 'r', 'LineWidth', 1.5);
    ylabel('Aperiodic Exponent', 'FontSize', 12);
    xlabel('Days Since Surgery', 'FontSize', 12);
    title({[hemisphers{hemi_idx} ' Hemisphere'], ['Fit Exponent'], ...
        ['R squared: ' num2str(r_squared_exponent, '%.4f')]}, 'FontSize', 14);
    grid on;

    % Trend line for aperiodic offset
    p_offset = polyfit(x_data, aperiodic_offset, 1);
    offset_fit = polyval(p_offset, x_data);
    r_squared_offset = calc_r_squared(aperiodic_offset, offset_fit);
    
    % Plot aperiodic offset
    subplot(2, 2, 5 - hemi_idx);
    plot(x_data, aperiodic_offset, '.', 'MarkerSize', 15);
    hold on;
    plot(x_data, offset_fit, 'r', 'LineWidth', 1.5);
    ylabel('Aperiodic Offset', 'FontSize', 12);
    xlabel('Days Since Surgery', 'FontSize', 12);
    title({[hemisphers{hemi_idx} ' Hemisphere'], ['Fit Offset'], ...
        ['R squared: ' num2str(r_squared_offset, '%.4f')]}, 'FontSize', 14);
    grid on;

    set(gcf, 'Color', 'w');

    savefig(['Peaks Statistics - fit params - ' Patient_name '.fig'])


    %% Fit parameters
    
    % Extract only beta band peaks
    beta_freq_indx = cellfun(@(x) x> 13 & x< 30, Event_Peak_freq, 'UniformOutput', false);
    theta_freq_indx = cellfun(@(x) x> 4 & x< 8, Event_Peak_freq, 'UniformOutput', false);
    %
    Event_Peak_freq_beta = cell(size(beta_freq_indx));
    Event_Peak_freq_theta = cell(size(theta_freq_indx));
    %
    Event_peak_height_beta = cell(size(beta_freq_indx));
    Event_peak_height_theta = cell(size(theta_freq_indx));
    %
    for i = 1:length(beta_freq_indx)
        Event_peak_height_beta{i} = Event_peak_height{i}(logical(beta_freq_indx{i}));
        Event_Peak_freq_beta{i} = Event_Peak_freq{i}(logical(beta_freq_indx{i}));
        %
        Event_peak_height_theta{i} = Event_peak_height{i}(logical(theta_freq_indx{i}));
        Event_Peak_freq_theta{i} = Event_Peak_freq{i}(logical(theta_freq_indx{i}));
    end


    % Extract Events Stats from theta and beta
    % beta freq
    peak_freq_num_beta = cell2mat(cellfun(@(x) numel(x), Event_Peak_freq_beta, 'UniformOutput', false));
    peak_freq_var_beta = cell2mat(cellfun(@(x) var(x), Event_Peak_freq_beta, 'UniformOutput', false));
    % beta height
    peak_height_var_beta = cell2mat(cellfun(@(x) var(x), Event_peak_height_beta, 'UniformOutput', false));
    peak_mean_height_beta = cell2mat(cellfun(@(x) mean(x), Event_peak_height_beta, 'UniformOutput', false));
    % theta freq
    peak_freq_num_theta = cell2mat(cellfun(@(x) numel(x), Event_Peak_freq_theta, 'UniformOutput', false));
    peak_freq_var_theta = cell2mat(cellfun(@(x) var(x), Event_Peak_freq_theta, 'UniformOutput', false));
    % theta height
    peak_height_var_theta = cell2mat(cellfun(@(x) var(x), Event_peak_height_theta, 'UniformOutput', false));
    peak_mean_height_theta = cell2mat(cellfun(@(x) mean(x), Event_peak_height_theta, 'UniformOutput', false));


    % Function to calculate RMS error
    calc_r_squared = @(y, y_fit) 1 - (sum((y - y_fit).^2) / sum((y - mean(y)).^2));

    % Plots for beta
    figure(6)
    
    % Remove NaNs for beta peaks
    valid_idx_beta_freq = ~isnan(peak_freq_num_beta);
    valid_idx_beta_height = ~isnan(peak_mean_height_beta);

    % Trend line for number of beta peaks
    p_freq = polyfit(Event_dates_since_surg(valid_idx_beta_freq), peak_freq_num_beta(valid_idx_beta_freq), 1);
    beta_freq_fit = polyval(p_freq, Event_dates_since_surg(valid_idx_beta_freq));
    r2_beta_freq = calc_r_squared(peak_freq_num_beta(valid_idx_beta_freq), beta_freq_fit);
    
    % Trend line for mean beta peak height
    p_height = polyfit(Event_dates_since_surg(valid_idx_beta_height), peak_mean_height_beta(valid_idx_beta_height), 1);
    beta_height_fit = polyval(p_height, Event_dates_since_surg(valid_idx_beta_height));
    r2_beta_height = calc_r_squared(peak_mean_height_beta(valid_idx_beta_height), beta_height_fit);

    % Remove NaNs for beta peak frequency variance and height variance
    valid_idx_beta_freq_var = ~isnan(peak_freq_var_beta);
    valid_idx_beta_height_var = ~isnan(peak_height_var_beta);

    % Trend line for beta peak frequency variance
    p_freq_var = polyfit(Event_dates_since_surg(valid_idx_beta_freq_var), peak_freq_var_beta(valid_idx_beta_freq_var), 1);
    beta_freq_var_fit = polyval(p_freq_var, Event_dates_since_surg(valid_idx_beta_freq_var));
    r2_beta_freq_var = calc_r_squared(peak_freq_var_beta(valid_idx_beta_freq_var), beta_freq_var_fit);

    % Trend line for beta peak height variance
    p_height_var = polyfit(Event_dates_since_surg(valid_idx_beta_height_var), peak_height_var_beta(valid_idx_beta_height_var), 1);
    beta_height_var_fit = polyval(p_height_var, Event_dates_since_surg(valid_idx_beta_height_var));
    r2_beta_height_var = calc_r_squared(peak_height_var_beta(valid_idx_beta_height_var), beta_height_var_fit);

    % Plot number of beta peaks
    subplot(4, 2, hemi_idx);
    plot(Event_dates_since_surg, peak_freq_num_beta, '.', 'MarkerSize', 15);
    hold on;
    plot(Event_dates_since_surg(valid_idx_beta_freq), beta_freq_fit, 'r', 'LineWidth', 1.5);
    ylabel('Number of Peaks', 'FontSize', 12);
    xlabel('Days Since Surgery', 'FontSize', 12);
    title({[hemisphers{hemi_idx} ' Hemisphere'], ['Number of Beta Peaks'], ...
        ['R^2 Error: ' num2str(r2_beta_freq, '%.4f')]}, 'FontSize', 14);
    grid on;

    % Plot mean beta peak height
    subplot(4, 2, hemi_idx + 2);
    plot(Event_dates_since_surg, peak_mean_height_beta, '.', 'MarkerSize', 15);
    hold on;
    plot(Event_dates_since_surg(valid_idx_beta_height), beta_height_fit, 'r', 'LineWidth', 1.5);
    ylabel('Mean Peak Height [μV]', 'FontSize', 12);
    xlabel('Days Since Surgery', 'FontSize', 12);
    title({[hemisphers{hemi_idx} ' Hemisphere'], ['Beta Mean Peak Height'], ...
        ['R^2 Error: ' num2str(r2_beta_height, '%.4f')]}, 'FontSize', 14);
    grid on;

    % Plot beta peak frequency variance
    subplot(4, 2, hemi_idx + 4);
    plot(Event_dates_since_surg, peak_freq_var_beta, '.', 'MarkerSize', 15);
    hold on;
    plot(Event_dates_since_surg(valid_idx_beta_freq_var), beta_freq_var_fit, 'r', 'LineWidth', 1.5);
    ylabel('Peak Frequency Variance', 'FontSize', 12);
    xlabel('Days Since Surgery', 'FontSize', 12);
    title({[hemisphers{hemi_idx} ' Hemisphere'], ['Beta Peak Frequency Variance'], ...
        ['R^2 Error: ' num2str(r2_beta_freq_var, '%.4f')]}, 'FontSize', 14);
    grid on;

    % Plot beta peak height variance
    subplot(4, 2, hemi_idx + 6);
    plot(Event_dates_since_surg, peak_height_var_beta, '.', 'MarkerSize', 15);
    hold on;
    plot(Event_dates_since_surg(valid_idx_beta_height_var), beta_height_var_fit, 'r', 'LineWidth', 1.5);
    ylabel('Peak Height Variance', 'FontSize', 12);
    xlabel('Days Since Surgery', 'FontSize', 12);
    title({[hemisphers{hemi_idx} ' Hemisphere'], ['Beta Peak Height Variance'], ...
        ['R^2 Error: ' num2str(r2_beta_height_var, '%.4f')]}, 'FontSize', 14);
    grid on;

    set(gcf, 'Color', 'w');
    savefig(['Peaks Statistics - Beta - ' Patient_name '.fig'])

    % Plots for theta
    figure(7)
    
    % Remove NaNs for theta peaks
    valid_idx_theta_freq = ~isnan(peak_freq_num_theta);
    valid_idx_theta_height = ~isnan(peak_mean_height_theta);

    % Trend line for number of theta peaks
    p_freq = polyfit(Event_dates_since_surg(valid_idx_theta_freq), peak_freq_num_theta(valid_idx_theta_freq), 1);
    theta_freq_fit = polyval(p_freq, Event_dates_since_surg(valid_idx_theta_freq));
    r2_theta_freq = calc_r_squared(peak_freq_num_theta(valid_idx_theta_freq), theta_freq_fit);
    
    % Trend line for mean theta peak height
    p_height = polyfit(Event_dates_since_surg(valid_idx_theta_height), peak_mean_height_theta(valid_idx_theta_height), 1);
    theta_height_fit = polyval(p_height, Event_dates_since_surg(valid_idx_theta_height));
    r2_theta_height = calc_r_squared(peak_mean_height_theta(valid_idx_theta_height), theta_height_fit);

    % Remove NaNs for theta peak frequency variance and height variance
    valid_idx_theta_freq_var = ~isnan(peak_freq_var_theta);
    valid_idx_theta_height_var = ~isnan(peak_height_var_theta);

    % Trend line for theta peak frequency variance
    p_freq_var = polyfit(Event_dates_since_surg(valid_idx_theta_freq_var), peak_freq_var_theta(valid_idx_theta_freq_var), 1);
    theta_freq_var_fit = polyval(p_freq_var, Event_dates_since_surg(valid_idx_theta_freq_var));
    r2_theta_freq_var = calc_r_squared(peak_freq_var_theta(valid_idx_theta_freq_var), theta_freq_var_fit);

    % Trend line for theta peak height variance
    p_height_var = polyfit(Event_dates_since_surg(valid_idx_theta_height_var), peak_height_var_theta(valid_idx_theta_height_var), 1);
    theta_height_var_fit = polyval(p_height_var, Event_dates_since_surg(valid_idx_theta_height_var));
    r2_theta_height_var = calc_r_squared(peak_height_var_theta(valid_idx_theta_height_var), theta_height_var_fit);

    % Plot number of theta peaks
    subplot(4, 2, hemi_idx);
    plot(Event_dates_since_surg, peak_freq_num_theta, '.', 'MarkerSize', 15);
    hold on;
    plot(Event_dates_since_surg(valid_idx_theta_freq), theta_freq_fit, 'r', 'LineWidth', 1.5);
    ylabel('Number of Peaks', 'FontSize', 12);
    xlabel('Days Since Surgery', 'FontSize', 12);
    title({[hemisphers{hemi_idx} ' Hemisphere'], ['Number of theta Peaks'], ...
        ['R^2 Error: ' num2str(r2_theta_freq, '%.4f')]}, 'FontSize', 14);
    grid on;

    % Plot mean theta peak height
    subplot(4, 2, hemi_idx + 2);
    plot(Event_dates_since_surg, peak_mean_height_theta, '.', 'MarkerSize', 15);
    hold on;
    plot(Event_dates_since_surg(valid_idx_theta_height), theta_height_fit, 'r', 'LineWidth', 1.5);
    ylabel('Mean Peak Height [μV]', 'FontSize', 12);
    xlabel('Days Since Surgery', 'FontSize', 12);
    title({[hemisphers{hemi_idx} ' Hemisphere'], ['theta Mean Peak Height'], ...
        ['R^2 Error: ' num2str(r2_theta_height, '%.4f')]}, 'FontSize', 14);
    grid on;

    % Plot theta peak frequency variance
    subplot(4, 2, hemi_idx + 4);
    plot(Event_dates_since_surg, peak_freq_var_theta, '.', 'MarkerSize', 15);
    hold on;
    plot(Event_dates_since_surg(valid_idx_theta_freq_var), theta_freq_var_fit, 'r', 'LineWidth', 1.5);
    ylabel('Peak Frequency Variance', 'FontSize', 12);
    xlabel('Days Since Surgery', 'FontSize', 12);
    title({[hemisphers{hemi_idx} ' Hemisphere'], ['theta Peak Frequency Variance'], ...
        ['R^2 Error: ' num2str(r2_theta_freq_var, '%.4f')]}, 'FontSize', 14);
    grid on;

    % Plot theta peak height variance
    subplot(4, 2, hemi_idx + 6);
    plot(Event_dates_since_surg, peak_height_var_theta, '.', 'MarkerSize', 15);
    hold on;
    plot(Event_dates_since_surg(valid_idx_theta_height_var), theta_height_var_fit, 'r', 'LineWidth', 1.5);
    ylabel('Peak Height Variance', 'FontSize', 12);
    xlabel('Days Since Surgery', 'FontSize', 12);
    title({[hemisphers{hemi_idx} ' Hemisphere'], ['theta Peak Height Variance'], ...
        ['R^2 Error: ' num2str(r2_theta_height_var, '%.4f')]}, 'FontSize', 14);
    grid on;

    set(gcf, 'Color', 'w');
    savefig(['Peaks Statistics - theta - ' Patient_name '.fig'])


    %% some more plots (for EC)
    After_surgery_indx = 1:87;
    Current_indx = 89:numel(peak_freq_num_beta);

    % Dividing each vector
    peak_freq_num_beta_after = peak_freq_num_beta(After_surgery_indx);
    peak_freq_num_beta_current = peak_freq_num_beta(Current_indx);
    [p_peaks_num, ~] = ranksum(peak_freq_num_beta_after,peak_freq_num_beta_current);
    
    figure()
    bar([mean(peak_freq_num_beta_after), mean(peak_freq_num_beta_current)], 'FaceColor', [124 169 165]/255)
    xticklabels({'After surgery', 'Currently'})
    ylabel('Number of Peaks')
    title('Number of Beta Peaks')
    box off
    set(gca, 'LineWidth',2)
    set(gca, 'FontSize', 14)
    set(gcf, 'Color', 'w')
    
    peak_freq_var_beta_after = peak_freq_var_beta(After_surgery_indx);
    peak_freq_var_beta_current = peak_freq_var_beta(Current_indx);
    
    peak_height_var_beta_after = peak_height_var_beta(After_surgery_indx);
    peak_height_var_beta_current = peak_height_var_beta(Current_indx);
    
    peak_mean_height_beta_after = peak_mean_height_beta(After_surgery_indx);
    peak_mean_height_beta_current = peak_mean_height_beta(Current_indx);
    

    %%
    aperiodic_exponent = unique(aperiodic_exponent, 'stable');
    Current_indx = 89:numel(aperiodic_exponent);
    aperiodic_exponent_after = aperiodic_exponent(After_surgery_indx);
    aperiodic_exponent_current = aperiodic_exponent(Current_indx);
    
    aperiodic_offset_after = aperiodic_offset(After_surgery_indx);
    aperiodic_offset_current = aperiodic_offset(Current_indx);

    figure()
    histogram(aperiodic_exponent_after, 'normalization', 'probability');
    hold on
    histogram(aperiodic_exponent_current, 'normalization', 'probability')
    [p_exp, ~] = ranksum(aperiodic_exponent_after,aperiodic_exponent_current);

    figure()
    histogram(aperiodic_offset_after, 'normalization', 'probability');
    hold on
    histogram(aperiodic_offset_current, 'normalization', 'probability')
    [p_offs, real_perry] = ttest2(aperiodic_offset_after,aperiodic_offset_current);


%%

    % Prepare data for bar plot
    means = [   mean(aperiodic_exponent_after,"omitnan"), mean(aperiodic_exponent_current,"omitnan");
        mean(aperiodic_offset_after,"omitnan"), mean(aperiodic_offset_current,"omitnan");];

    
    %means = normalize(means, 2, "scale");
    % Create a figure
    figure;
    
    % Bar plot
    b = bar(means, 'grouped');
    
    % Set x-axis labels
    set(gca, 'XTickLabel', {'Peak Freq Num', 'Peak Freq Var', 'Peak Height Var', 'Peak Mean Height'});
    
    % Add legend
    legend({'After Surgery', 'Current'});
    
    % Add titles and labels
    title('Comparison of Metrics: After Surgery vs. Current');
    ylabel('Mean Value');
    xlabel('Parameters');
    
    % Improve layout
    grid on;
end