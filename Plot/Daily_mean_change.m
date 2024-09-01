close all
Patient_name = 'Patient_EC';
hemisphers = {'Left', 'Right'};

% save surgery date
surg_date = DBS_data.(Patient_name).Info.deviceInfo.ImplantDate;
surg_date = datetime(regexprep(surg_date(1:end-1),'T',' '));

for hemi_idx = 1:numel(hemisphers)
    hemi = hemisphers{hemi_idx};

    Events = DBS_data.(Patient_name).Events.Events_data;
    LFP_T = DBS_data.(Patient_name).TrendLogs.([hemisphers{hemi_idx} '_Hemi']).LFP_table;
    % For multiple patients analysis - normalization to [0 1] is required:
    % Change those 2 lines in the future for single patient analysis:
    LFP_T = normalize(LFP_T(:,:),'range');
    LFP_vec = normalize(reshape(table2array(LFP_T),[],1),'range');
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
    TEED_vec = TEED_vec(~nanIdx);
    sensing_vec = sensing_vec(~nanIdx);
    contacts_vec = contacts_vec(~nanIdx);
    try
        Datetime_vec = Datetime_vec(~nanIdx);
        days_since_surge_vec = days_since_surge_vec(~nanIdx);
    catch

    end
    
    % Only Beta
    beta_LFP_mat = LFP_T;
    beta_indx = sensing_T{:,:} > 12 & sensing_T{:,:} < 30;
    for i=1:length(rows_names)
        for j=1:(length(cols_names)-1)
            if beta_indx(i,j) == 1
                beta_LFP_mat{i,j} = LFP_T{i,j};
            else
                beta_LFP_mat{i,j} = NaN;
            end
        end
    end

    
    % find unique periods
    combined_strings = strcat(num2str(sensing_vec(:)), '_', contacts_vec(:));
    [unique_pairs, ~, pair_idx] = unique(combined_strings,'stable');

    % finding only unique sensing freqs:
    %uniq_sensing_freq = unique(sensing_vec,'stable');
    %uniq_sensing_freq = uniq_sensing_freq(~isnan(uniq_sensing_freq));
    savefig(['Daily LFP mean - ALL ' Patient_name '.fig'])


    day_start_indx = 6*6+1; % day start at 6 in the morning
    day_end_indx = 21*6; %  day end at 21 in the evening
    night_start_indx = day_end_indx+1;
    night_end_indx = day_start_indx-1;
    %% Daily-mean ( from 5:00 to 21:00)
    LFP_day = LFP_T(day_start_indx:day_end_indx,:);
    Amp_day = Amp_T(day_start_indx:day_end_indx,:);
    TEED_day = TEED_T(day_start_indx:day_end_indx,:);
    sensing_day = sensing_T(day_start_indx:day_end_indx,:);
    Stimrate_day = Stimrate_T(day_start_indx:day_end_indx,:);
    PulseWidth_day = PulseWidth_T(day_start_indx:day_end_indx,:);

    LFP_mean_day.([Patient_name '_' hemi]) = table2array(mean(LFP_day,1));
    Amp_mean = table2array(mean(Amp_day,1));
    TEED_mean = table2array(mean(TEED_day,1),"omitnan");
    sensing_mean = table2array(mean(sensing_day,1));
    Stimrate_mean = table2array(mean(Stimrate_day,1));
    PulseWidth_mean = table2array(mean(PulseWidth_day,1));
    %
    days_vec = datetime(LFP_T.Properties.VariableNames);
    days_since_surg_vec = days(days_vec - surg_date);

    % take only part of the days
    % LFP_mean_day = LFP_mean_day([1:23 89:97]);
    % days_since_surg_vec = days_since_surg_vec([1:23 89:97]);
    % TEED_mean = TEED_mean([1:23 89:97]);

    figure(hemi_idx)
    subplot(2,1,1)
    plot(days_since_surg_vec,LFP_mean_day.([Patient_name '_' hemi]),'-o', 'LineWidth',2)
    
    % xticks(1:32)
    %xticklabels(nights_since_surg_vec)

    ylim_l_1 = ylim;   
    ylabel('LFP Amplitude [\muV]')
    xlabel('Days Since Surgery')
    hold on

    yyaxis right
    plot(days_since_surg_vec,TEED_mean,'o', 'LineWidth',2)
    ylim_r_1 = ylim;   
    ylabel('TEED [E\sec]')
    grid on
    lgd = legend('LFP', 'TEED', 'Location','northeastoutside');
    set(lgd, 'Color', 'none');
    titletime = [num2str((day_start_indx-1)/6) ':00 - ' num2str(day_end_indx/6) ':00 '];
    title([hemisphers{hemi_idx} ' Daily LFP (' titletime ')'], 'FontSize',14, 'Interpreter', 'none');
    
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
        titletime = [num2str((day_start_indx-1)/6) ':00 - ' num2str(day_end_indx/6) ':00 '];
        title([hemisphers{hemi_idx} ' Daily LFP (' titletime ') - '  char(unique_pairs)], 'Interpreter', 'none')
    end

    % subplot(7,1,4)
    % yyaxis left
    % plot(days_vec,Stimrate_mean,'o', 'LineWidth',4)
    % yyaxis right
    % plot(days_vec,PulseWidth_mean,'o', 'LineWidth',4)
    % lgd = legend('Stim rate', 'Pulse Width','Location','northeast');
    % set(lgd, 'Color', 'none');
    % grid on
     set(gcf, 'Color', 'w')
    savefig(['Nightly LFP mean - ALL ' Patient_name '.fig'])


    %% Night-mean ( from 21:00 to 5:00 next day)
    LFP_night = [LFP_T(1:night_end_indx,:) ; LFP_T(night_start_indx:end,:)];
    Amp_night = [Amp_T(1:night_end_indx,:) ; Amp_T(night_start_indx:end,:)];
    TEED_night = [TEED_T(1:night_end_indx,:) ; TEED_T(night_start_indx:end,:)];
    sensing_night = [sensing_T(1:night_end_indx,:) ; sensing_T(night_start_indx:end,:)];
    Stimrate_night = [Stimrate_T(1:night_end_indx,:) ; Stimrate_T(night_start_indx:end,:)];
    PulseWidth_night = [PulseWidth_T(1:night_end_indx,:) ; PulseWidth_T(night_start_indx:end,:)];

    LFP_mean_night.([Patient_name '_' hemi]) = table2array(mean(LFP_night,1));
    Amp_mean = table2array(mean(Amp_night,1));
    TEED_mean = table2array(mean(TEED_night,1, "omitnan"));
    sensing_mean = table2array(mean(sensing_night,1));
    Stimrate_mean = table2array(mean(Stimrate_night,1));
    PulseWidth_mean = table2array(mean(PulseWidth_night,1));

    %
    nights_vec = datetime(LFP_T.Properties.VariableNames);
    nights_since_surg_vec = days(nights_vec - surg_date);

    % take only part of the days
    % LFP_mean_night = LFP_mean_night([1:23 89:97]);
    % nights_since_surg_vec = nights_since_surg_vec([1:23 89:97]);
    % TEED_mean = TEED_mean([1:23 89:97]);

    subplot(2,1,2)
    plot(nights_since_surg_vec,LFP_mean_night.([Patient_name '_' hemi]),'-o', 'LineWidth',2)
    
    % xticks(1:32)
    %xticklabels(nights_since_surg_vec)
    
    ylabel('LFP Amplitude [\muV]')
    xlabel('Days Since Surgery')
    ylim_l_2 = ylim;   
    hold on

    yyaxis right
    plot(nights_since_surg_vec,TEED_mean,'o', 'LineWidth',2)
    ylim_r_2 = ylim;   
    ylabel('TEED [E\sec]')

    grid on
    lgd = legend('LFP', 'TEED', 'Location','northeastoutside');
    set(lgd, 'Color', 'none');
    title([hemisphers{hemi_idx} ' Night LFP (' titletime ')'], 'FontSize',14, 'Interpreter', 'none');

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
        titletime = [num2str(day_end_indx/6) ':00 - ' num2str((day_start_indx-1)/6)  ':00 '];
        title([hemisphers{hemi_idx} ' Night LFP (' titletime ') - '  char(unique_pairs)], 'Interpreter', 'none')
    end


    % Set the same ylim for both subplots
    common_ylim_r = [min(ylim_r_1(1), ylim_r_2(1)), max(ylim_r_1(2), ylim_r_2(2))];
    common_ylim_l = [min(ylim_l_1(1), ylim_l_2(1)), max(ylim_l_1(2), ylim_l_2(2))];

    subplot(2,1,1)
    yyaxis left
    ylim(common_ylim_l);  
    yyaxis right
    ylim(common_ylim_r);  

    subplot(2,1,2)
    yyaxis left
    ylim(common_ylim_l);  
    yyaxis right
    ylim(common_ylim_r);

    % subplot(6,1,6)
    % yyaxis left
    % plot(days_vec,Stimrate_mean,'o', 'LineWidth',4)
    % yyaxis right
    % plot(days_vec,PulseWidth_mean,'o', 'LineWidth',4)
    % lgd = legend('Stim rate', 'Pulse Width','Location','northeast');
    % set(lgd, 'Color', 'none');
    % grid on
    set(gcf, 'Color', 'w')


    %% connected dot plot for day vs night
    figure(3)
    subplot(1,2,hemi_idx)
    total_LFP_mean_night = [];
    patients = fieldnames(LFP_mean_night);
    for i = 1:length(fieldnames(LFP_mean_night))
        total_LFP_mean_night = [total_LFP_mean_night LFP_mean_night.([patients{i}])];
    end

    total_LFP_mean_day = [];
    for i = 1:length(fieldnames(LFP_mean_day))
        total_LFP_mean_day = [total_LFP_mean_day LFP_mean_day.([patients{i}])];
    end
    % For specific patient:
    % daboxplot([LFP_mean_day', LFP_mean_night.(Patient_name)'], 'xtlabels', {'Day', 'Night'}, 'withinlines', 1, 'scatter', 2, 'fill', 0);
    % For all patients:
    daboxplot([total_LFP_mean_day', total_LFP_mean_night'], 'xtlabels', {'Day', 'Night'}, 'withinlines', 1, 'scatter', 2, 'fill', 0);
    title([hemisphers{hemi_idx}]);
    xlabel('Time of Day');
    ylabel('LFP Mean Amplitude [µV]');
    box off
    set(gca, 'LineWidth', 3)
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 14)
    sgtitle('LFP means at day and night', 'FontSize', 18)
  
    % set(gcf, 'color', [124 169 165]/255);    
    % set(gca, 'color', [124 169 165]/255);
    % savefig(['Hemisphere LFP BoxPlot - ' Patient_name '.fig'])

    % export_fig('LFP Mean Amplitude', '-png')
    % exportgraphics(gcf,'LFP Mean Amplitude.png')
    
    %% save statistics
    
    % For one patient:
    p_name = [Patient_name '_' hemisphers{hemi_idx}];
    [pval_unpaired.(p_name), ~] = ranksum(LFP_mean_day.([Patient_name '_' hemi]), LFP_mean_night.([Patient_name '_' hemi]));
    [pval_paired.(p_name), ~] = signrank(LFP_mean_day.([Patient_name '_' hemi]), LFP_mean_night.([Patient_name '_' hemi]));
    
    % For all patients:
    [pval_unpaired_all.(hemisphers{hemi_idx}), ~] = ranksum(total_LFP_mean_day, total_LFP_mean_night);
    [pval_paired_all.(hemisphers{hemi_idx}), ~] = signrank(total_LFP_mean_day, total_LFP_mean_night);

    
    % save('pval_daily_mean_unpaired.mat', 'pval_unpaired')
    % save('pval_daily_mean_paired.mat', 'pval_paired')
end