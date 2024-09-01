%% VARIABLES
close all
Patient_name = 'Patient_EC';
hemisphers = {'Right', 'Left'};

% save surgery date
surg_date = DBS_data.(Patient_name).Info.deviceInfo.ImplantDate;
surg_date = datetime(regexprep(surg_date(1:end-1),'T',' '));

for hemi_idx = 1:numel(hemisphers)
    LFP_T = DBS_data.(Patient_name).TrendLogs.([hemisphers{hemi_idx} '_Hemi']).LFP_table;
    LFP_T = normalize(LFP_T(:,:),'range');
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
    Datetime_vec = DBS_data.(Patient_name).TrendLogs.([hemisphers{hemi_idx} '_Hemi']).Time_vector;
    days_since_surge_vec = days(Datetime_vec - surg_date);

    %
    cols_names = LFP_T.Properties.VariableNames;
    rows_names = LFP_T.Properties.RowNames;
    % nans
    nanIdx = isnan(LFP_vec);
    LFP_vec = LFP_vec(~nanIdx);
    Amp_vec = Amp_vec(~nanIdx);
    sensing_vec = sensing_vec(~nanIdx);
    contacts_vec = contacts_vec(~nanIdx);
    %Datetime_vec = Datetime_vec(~nanIdx);
    % Events
    Events = DBS_data.(Patient_name).Events.Events_data;
    
    % find unique periods
    combined_strings = strcat(num2str(sensing_vec(:)), '_', contacts_vec(:));
    [unique_pairs, ~, pair_idx] = unique(combined_strings,'stable');

    % finding only unique sensing freqs:
    %uniq_sensing_freq = unique(sensing_vec,'stable');
    %uniq_sensing_freq = uniq_sensing_freq(~isnan(uniq_sensing_freq));
   
   %% Plot LFP + Events 
   for i = 1:length(unique_pairs)
       current_pair = unique_pairs(i);
       pair_parts = split(unique_pairs(i), '_');
       current_freq = pair_parts(1);
       current_contact1 = pair_parts(2);
       current_contact2 = pair_parts(4);
       LFP_sense =  LFP_vec(pair_idx == i);
       sensing_datetime = Datetime_vec(pair_idx == i);
        
       % remove outilying dates
       sensing_datetime_diff = diff(sensing_datetime);
       mean_datetime_diff = mean(sensing_datetime_diff);
       leap_indx = find(sensing_datetime_diff>mean_datetime_diff);
       if ~isempty(leap_indx)
           LFP_sense = LFP_sense(leap_indx(1)+1:end);
           sensing_datetime = sensing_datetime(leap_indx(1)+1:end);
       end

       
       % events
       events_sense_indx = (Events.DateTime<sensing_datetime(end) & Events.DateTime>sensing_datetime(1) & Events.LFP == 1);
       events_sense = Events(events_sense_indx,:);
       events_sense = events_sense(strcmp(events_sense.EventName,{'Med'}) == 1 | strcmp(events_sense.EventName, {'Postmed'}) == 1 |...
           strcmp(events_sense.EventName, {'Took Medication'}) == 1 | strcmp(events_sense.EventName,{'Medication'}) == 1,:);
       event_types = unique(events_sense.EventName, 'stable')';
       
       figure(i)
       colors = {'magenta', 'green', 'red', 'blue'};
       subplot(2,1,hemi_idx)
       plot(sensing_datetime, LFP_sense);
       xlabel('Time of Day');
       ylabel('LFP Amplitude [\muV]');
       for j = 1:size(events_sense,1)
           color_indx = find(strcmp(event_types,events_sense.EventName(j)));
           xline(events_sense.DateTime(j), '-','Color',colors{color_indx});
       end
       lgd = legend(['data', event_types]);   
       title({[hemisphers{hemi_idx} ' LFP'], ['Sensing frequency - ' char(current_freq) ' [Hz] '], ['Between contacts ' char(current_contact1) ' and ' char(current_contact2)]});
       sgtitle('LFP Plots');
       set(gcf,'Color','w')
   end
   %% PLOT LFP CLUSTERS -
   
   % printing first vs last 2 days of LFP
   
   % figure(800)
   % subplot(2,1,hemi_idx)
   % plot(LFP_vec(1:288), '*');
   % xlabel('Bins');
   % ylabel('LFP Amplitude [\muV]');
   % hold on
   % plot(LFP_vec(end-288:end), '*');
   % xlabel('Bins');
   % ylabel('LFP Amplitude [\muV]');
   % title([hemisphers{hemi_idx} ' LFP Right After Surgery V.S LFP Now' ]);
   % 
   % set(gcf,'Color','w')
   % 
   % % printing 
   % LFP_mat = table2array(LFP_T);
   % LFP_Start = LFP_mat(:,1:4);
   % LFP_End = LFP_mat(:,end-3:end);
   % 
   % LFP_Start_mean = mean(LFP_Start, 2, 'omitnan');
   % LFP_End_mean = mean(LFP_End, 2, 'omitnan');
   % 
   % [h,p_value] = ttest(LFP_Start_mean,LFP_End_mean);
   % time_ax = datetime(0, 0, 0, 0, 0:10:1430, 0);
   % 
   %    figure(801)
   %     subplot(2,1,hemi_idx)
   %     plot(time_ax, LFP_Start_mean, '*');
   %     xlabel('Bins');
   %     ylabel('LFP Amplitude [\muV]');
   %     hold on
   %     plot(time_ax, LFP_End_mean, '*');
   %     xlabel('Time Of Day');
   %     ylabel('LFP Amplitude [\muV]');
   %     title({[hemisphers{hemi_idx} ' LFP Right After Surgery V.S LFP Now'], ['P-Value - ' num2str(p_value)]});
   %     ax = gca;
   %     ax.XAxis.TickLabelFormat = 'HH:mm';
   %     legend('Right After Surgery', 'Now');
   %     set(gcf,'Color','w')
       
   %% Plot Heatmap
    % Extracting row names and column names
    time_of_day = LFP_T.Properties.RowNames;
    LFP_table_dates = LFP_T.Properties.VariableNames;
    time_of_day_dt = datetime(time_of_day, 'Format', 'HH:mm');
    dates_dt = cellfun(@(x) datetime(x, 'InputFormat', 'dd-MMM-yy'), LFP_table_dates, 'UniformOutput', false);
    dates_dt_vector = [dates_dt{:}];
    days_since_suge = round(days(dates_dt_vector - surg_date));

    LFP_matrix = table2array(LFP_T);

    figure;
    h = heatmap(days_since_suge, time_of_day_dt, LFP_matrix);
    h.XLabel = 'Date';
    h.YLabel = 'Time of Day';
    h.Title = 'LFP Heatmap';
    set(gcf, 'Color', 'w');
    
    % Formatting x-axis and y-axis for better readability
    %ax = gca;
    %ax.YDisplayLabels = datestr(time_of_day_dt(round(linspace(1, numel(time_of_day_dt), 10))), 'HH:MM');
    %ax.XAxis.FontSize = 12;  % Adjust the font size as needed

    %% Histogram & lillies
    theta_LFP_mat = zeros(size(LFP_T));
    beta_LFP_mat = zeros(size(LFP_T));
    %
    theta_indx = sensing_T{:,:} > 5 & sensing_T{:,:} <13;
    beta_indx = sensing_T{:,:} > 13 & sensing_T{:,:} < 30;

    for i=1:length(rows_names)
        for j=1:(length(cols_names)-1)
            %THETA
            if theta_indx(i,j) == 1 
                theta_LFP_mat(i,j) = LFP_T{i,j};
            else
                theta_LFP_mat(i,j) = NaN;
            end
            % BETA
            if beta_indx(i,j) == 1
                beta_LFP_mat(i,j) = LFP_T{i,j};
            else
                beta_LFP_mat(i,j) = NaN;
            end
        end
    end

    % 0.5 reshape to vectors
    theta_LFP_vec = reshape(theta_LFP_mat, [], 1);
    beta_LFP_vec = reshape(beta_LFP_mat, [], 1);

    % 1. checking normal distribution of the full LFP vector
    if ~isempty(theta_LFP_vec) 
        lilli_theta.(hemisphers{hemi_idx}) = lillietest(theta_LFP_vec);
    end
    if ~isempty(beta_LFP_vec)
        lilli_beta.(hemisphers{hemi_idx}) = lillietest(beta_LFP_vec);
    end
    % both are not normal -> right skewed dist.


    % 2. checking normal distribution of the day/night LFP vectors seperately

    LFP_day_theta = theta_LFP_mat(37:127,:);
    LFP_night_theta = [theta_LFP_mat(1:37,:); theta_LFP_mat(127:144,:)];
    %
    LFP_day_beta = beta_LFP_mat(37:127,:);
    LFP_night_beta = [beta_LFP_mat(1:37,:); beta_LFP_mat(127:144,:)];

    %
    lilli_day_theta.(hemisphers{hemi_idx}) = lillietest(reshape(LFP_day_theta, [], 1));
    lilli_night_theta.(hemisphers{hemi_idx}) = lillietest(reshape(LFP_night_theta, [], 1));
    %
    lilli_day_beta.(hemisphers{hemi_idx}) = lillietest(reshape(LFP_day_beta, [], 1));
    lilli_night_beta.(hemisphers{hemi_idx}) = lillietest(reshape(LFP_night_beta, [], 1));

    %%% Everything isn't normal -> right skewed dist.

    % Skweness measures
    skew_day_theta.(hemisphers{hemi_idx}) = skewness(reshape(LFP_day_theta, [], 1));
    skew_night_theta.(hemisphers{hemi_idx}) = skewness(reshape(LFP_night_theta, [], 1));
    %
    skew_day_beta.(hemisphers{hemi_idx}) = skewness(reshape(LFP_day_beta, [], 1));
    skew_night_beta.(hemisphers{hemi_idx}) = skewness(reshape(LFP_night_beta, [], 1));


    % for plotting multiple patients:
    % LFP_theta_day_struct.([Patient_name '_' hemisphers{hemi_idx}]) = reshape(LFP_day_theta,[],1);
    % LFP_theta_night_struct.([Patient_name '_' hemisphers{hemi_idx}]) = reshape(LFP_night_theta,[],1);
    % LFP_beta_day_struct.([Patient_name '_' hemisphers{hemi_idx}]) = reshape(LFP_day_beta,[],1);
    % LFP_beta_night_struct.([Patient_name '_' hemisphers{hemi_idx}]) = reshape(LFP_night_beta,[],1);
    % for plotting multiple patients:
    LFP_theta_day_struct.([Patient_name '_' hemisphers{hemi_idx}]) = mean(LFP_day_theta, 'omitnan');
    LFP_theta_night_struct.([Patient_name '_' hemisphers{hemi_idx}]) = mean(LFP_night_theta, 'omitnan');
    LFP_beta_day_struct.([Patient_name '_' hemisphers{hemi_idx}]) = mean(LFP_day_beta, 'omitnan');
    LFP_beta_night_struct.([Patient_name '_' hemisphers{hemi_idx}]) = mean(LFP_night_beta, 'omitnan');

    total_LFP_theta_day = [];
    total_LFP_theta_night = [];
    total_LFP_beta_day = [];
    total_LFP_beta_night = [];
    patients = fieldnames(LFP_beta_night_struct);
    for i = 1:length(fieldnames(LFP_beta_night_struct))
        total_LFP_theta_day = [total_LFP_theta_day , LFP_theta_day_struct.(patients{i})];
        total_LFP_theta_night = [total_LFP_theta_night , LFP_theta_night_struct.(patients{i})];
        total_LFP_beta_day = [total_LFP_beta_day , LFP_beta_day_struct.(patients{i})];
        total_LFP_beta_night = [total_LFP_beta_night , LFP_beta_night_struct.(patients{i})];
    end

    total_LFP_theta_day = total_LFP_theta_day(~isnan(total_LFP_theta_day));
    total_LFP_theta_night = total_LFP_theta_night(~isnan(total_LFP_theta_night));
    total_LFP_beta_day = total_LFP_beta_day(~isnan(total_LFP_beta_day));
    total_LFP_beta_night = total_LFP_beta_night(~isnan(total_LFP_beta_night));




    % Plotting
    figure(8)
    set(gcf, 'Color','w')
    
    day_color = [124 169 165]/255;  % Light blue
    night_color = [0.85 0.33 0.1];  % Orange

    subplot(2,2,hemi_idx*2-1)
    histogram(reshape(total_LFP_theta_day, [], 1),30, 'Normalization', 'probability', 'FaceColor', day_color, 'BinLimits',[0 1])
    hold on
    histogram(reshape(total_LFP_theta_night, [], 1),30, 'Normalization', 'probability', 'FaceColor', night_color, 'BinLimits',[0 1])
    title([hemisphers{hemi_idx} ' Hemisphere LFP distribution - Theta'], 'FontSize', 18)
    legend('day - theta', 'night - theta', 'FontSize', 14)

    subplot(2,2,hemi_idx*2)
    histogram(reshape(total_LFP_beta_day, [], 1),30, 'Normalization', 'probability', 'FaceColor', day_color, 'BinLimits',[0 1])
    hold on
    histogram(reshape(total_LFP_beta_night, [], 1),30, 'Normalization', 'probability', 'FaceColor', night_color, 'BinLimits',[0 1])
    title([hemisphers{hemi_idx} ' Hemisphere LFP distribution - Beta'], 'FontSize', 18)
    legend('day - beta', 'night - beta', 'FontSize', 14)
    xlabel('LFP Amplitude [µV]', 'FontSize', 14)
    ylabel('probability', 'FontSize', 14)
    set(gca,'FontSize',14)
    set(gca, 'LineWidth',2)
    box off
    % savefig('Hemisphere LFP distribution - AL.fig')


    %% Diurnal plot
    %
    LFP_theta_mean = mean(theta_LFP_mat,2, 'omitnan');
    LFP_theta_std = std(theta_LFP_mat,0,2, 'omitnan');

    LFP_beta_mean = mean(beta_LFP_mat,2, 'omitnan');
    LFP_beta_std = std(beta_LFP_mat,0,2, 'omitnan');
    
    % for plotting multiple patients:
    % LFP_theta_mean_struct.([Patient_name '_' hemisphers{hemi_idx}]) = mean(theta_LFP_mat,2, 'omitnan');
    % LFP_theta_std_struct.([Patient_name '_' hemisphers{hemi_idx}]) = std(theta_LFP_mat,0,2, 'omitnan');
    % 
    % LFP_beta_mean_struct.([Patient_name '_' hemisphers{hemi_idx}]) = mean(beta_LFP_mat,2, 'omitnan');
    % LFP_beta_std_struct.([Patient_name '_' hemisphers{hemi_idx}]) = std(beta_LFP_mat,0,2, 'omitnan');
    % 
    % total_LFP_theta_means = [];
    % total_LFP_theta_stds = [];
    % total_LFP_beta_means = [];
    % total_LFP_bets_stds = [];
    % patients = fieldnames(LFP_theta_mean_struct);
    % for i = 1:length(fieldnames(LFP_theta_mean_struct))
    %     total_LFP_theta_means = [total_LFP_theta_means ; LFP_theta_mean_struct.(patients{i})'];
    %     total_LFP_theta_stds = [total_LFP_theta_stds ; LFP_theta_std_struct.(patients{i})'];
    %     total_LFP_beta_means = [total_LFP_beta_means ; LFP_beta_mean_struct.(patients{i})'];
    %     total_LFP_bets_stds = [total_LFP_bets_stds ; LFP_beta_std_struct.(patients{i})'];
    % end
    % final_LFP_theta_mean = mean(total_LFP_theta_means,1);
    % final_LFP_theta_mean = mean(total_LFP_theta_stds,1);
    % final_LFP_theta_mean = mean(total_LFP_beta_means,1);
    % final_LFP_theta_mean = mean(total_LFP_bets_stds,1);

    % rename instead of changing the names in the plot:
    % LFP_theta_mean = final_LFP_theta_mean;
    % LFP_theta_std = final_LFP_theta_mean;
    % LFP_beta_mean = final_LFP_theta_mean;
    % LFP_beta_std = final_LFP_theta_mean;



    % figure(9)
    % 
    % subplot(2,2,hemi_idx*2-1)
    % polarplot(deg2rad(linspace(1,360,144)),LFP_theta_mean, 'LineWidth', 2);
    % hold on
    % polarplot(deg2rad(linspace(1,360,144)),LFP_theta_mean-LFP_theta_std,'--', 'Color','c', 'LineWidth', 1);
    % polarplot(deg2rad(linspace(1,360,144)),LFP_theta_mean+LFP_theta_std,'--', 'Color','c', 'LineWidth', 1);
    % ax = gca;
    % ax.ThetaTick = (0:15:345); % Customize theta ticks for each hour
    % ax.ThetaTickLabel = arrayfun(@(h) sprintf('%02d:00', h), 0:23, 'UniformOutput', false);
    % set(ax, 'FontSize', 12);
    % %set(ax, 'FontName', 'Times New Roman');
    % ax.ThetaZeroLocation = 'top'; % Adjust the 0-degree location to the top
    % ax.ThetaDir = 'clockwise'; % Adjust the direction of rotation
    % title(['Theta  Mean Amplitude - ' hemisphers{hemi_idx} ' Hemisphere']);
    % rlim_1 = rlim; 
    % 
    % subplot(2,2,hemi_idx*2)
    % polarplot(deg2rad(linspace(1,360,144)),LFP_beta_mean, 'LineWidth', 2);
    % hold on
    % polarplot(deg2rad(linspace(1,360,144)),LFP_beta_mean-LFP_beta_std,'--', 'Color','c', 'LineWidth', 1);
    % polarplot(deg2rad(linspace(1,360,144)),LFP_beta_mean+LFP_beta_std,'--', 'Color','c', 'LineWidth', 1);
    % ax = gca;
    % ax.ThetaTick = (0:15:345); % Customize theta ticks for each hour
    % ax.ThetaTickLabel = arrayfun(@(h) sprintf('%02d:00', h), 0:23, 'UniformOutput', false);
    % set(ax, 'FontSize', 14);
    % % set(ax, 'FontName', 'Times New Roman');
    % ax.ThetaZeroLocation = 'top'; % Adjust the 0-degree location to the top
    % ax.ThetaDir = 'clockwise'; % Adjust the direction of rotation
    % title(['Beta Mean Amplitude - ' hemisphers{hemi_idx} ' Hemisphere'], 'FontSize', 18);
    % rlim_2 = rlim; 
    % %
    % sgtitle('24 hours Diurnal Cycle', 'FontSize', 18, 'FontName', 'Times New Roman');
    % set(gcf,'Color','w')
    % legend('Mean Amplitude', '+- 1 STD','FontSize', 16);
    % % Set the same ylim for both subplots
    % common_rlim = [min(rlim_1(1), rlim_2(1)), max(rlim_1(2), rlim_2(2))];
    % subplot(2,2,hemi_idx*2-1)
    % rlim(common_rlim);  
    % subplot(2,2,hemi_idx*2)
    % rlim(common_rlim);


    figure(10)
    subplot(1,2,3-hemi_idx)
    polarplot(deg2rad(linspace(1,360,144)),LFP_theta_mean, 'Color','#77AC30', 'LineWidth', 2);
    hold on
    polarplot(deg2rad(linspace(1,360,144)),LFP_theta_mean-LFP_theta_std,'--', 'Color','g', 'LineWidth', 1);
    polarplot(deg2rad(linspace(1,360,144)),LFP_theta_mean+LFP_theta_std,'--', 'Color','g', 'LineWidth', 1);
    polarplot(deg2rad(linspace(1,360,144)),LFP_beta_mean, 'Color','c', 'LineWidth', 2);
    polarplot(deg2rad(linspace(1,360,144)),LFP_beta_mean-LFP_beta_std,'--', 'Color',[124 169 165]/255, 'LineWidth', 1);
    polarplot(deg2rad(linspace(1,360,144)),LFP_beta_mean+LFP_beta_std,'--', 'Color',[124 169 165]/255, 'LineWidth', 1);
    ax = gca;
    ax.ThetaTick = (0:15:345); % Customize theta ticks for each hour
    ax.ThetaTickLabel = arrayfun(@(h) sprintf('%02d:00', h), 0:23, 'UniformOutput', false);
    set(ax, 'FontSize', 12);
    ax.ThetaZeroLocation = 'top'; % Adjust the 0-degree location to the top
    ax.ThetaDir = 'clockwise'; % Adjust the direction of rotation
    title(['Band Power - ' hemisphers{hemi_idx} ' Hemisphere'],'FontSize', 18);
    legend('Theta', '+- 1 STD Theta', '','Beta', '+- 1 STD Beta','FontSize', 16);
    sgtitle('24 hours Diurnal Cycle', 'FontSize', 18, 'FontName', 'Times New Roman');
    set(gcf,'Color','w')

    % Set the same ylim for both subplots
    if hemi_idx == 2
    subplot(1,2,1)
    rlim_1 = rlim;
    subplot(1,2,2)
    rlim_2 = rlim;

    common_rlim = [min(rlim_1(1), rlim_2(1)), max(rlim_1(2), rlim_2(2))];
    subplot(1,2,1)
    rlim(common_rlim);  
    subplot(1,2,2)
    rlim(common_rlim);
    end

    savefig(['Diurnal Clock - ' Patient_name '.fig'])

end

%% Avarage day plot
% % Assuming your LFP data is in the LFP_table
% lfpDataRight = table2array(TrendLogs.Right_Hemi.LFP_table);
% lfpDataLeft = table2array(TrendLogs.Left_Hemi.LFP_table);
%
% % Replace NaN values with a placeholder (e.g., -1)
% lfpDataRight(isnan(lfpDataRight)) = -1;
% lfpDataLeft(isnan(lfpDataLeft)) = -1;
%
% fig3 = figure('Visible', 'off');
% sgtitle('Diurnal change of LFP power')
% subplot(2,1,1)
% boxplot(transpose(lfpDataRight),'MedianStyle','line','PlotStyle', 'compact','DataLim', [0 2000])
% title('Right hemi')
% ylim([0 2000])
% xlabel('10-min bins')
% ylabel('LFP Power')
% hold on
% plot(TrendLogs.Right_Hemi.LFP_Statistics{:,3},'LineWidth',2,'Color','r')
% subplot(2,1,2)
% boxplot(transpose(lfpDataLeft),'MedianStyle','line','PlotStyle', 'compact','DataLim', [0 2000])
% title('Left hemi')
% xlabel('10-min bins')
% ylabel('LFP Power')
% ylim([0 2000])
% hold on
% plot(TrendLogs.Left_Hemi.LFP_Statistics{:,3},'LineWidth',2,'Color','r')
% TrendLogs.plots.Daily_Boxplot.(session_date) = fig3;
