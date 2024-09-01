Patient_name = 'Patient_YB';
hemisphers = {'Right', 'Left'};
for hemi_idx = 1:numel(hemisphers)
    LFP_T = DBS_data.(Patient_name).TrendLogs.([hemisphers{hemi_idx} '_Hemi']).LFP_table;
    LFP_vec = reshape(table2array(LFP_T),[],1);
    Amp_vec = reshape(table2array(DBS_data.(Patient_name).TrendLogs.([hemisphers{hemi_idx} '_Hemi']).Amp_table),[],1);
    sensing_T = DBS_data.(Patient_name).Groups.([hemisphers{hemi_idx} '_Hemi']).Sensing_Freq;
    sensing_vec = reshape(table2array(sensing_T),[],1);
    contacts_T = DBS_data.(Patient_name).Groups.([hemisphers{hemi_idx} '_Hemi']).Contacts;
    contacts_vec = reshape(table2array(contacts_T),[],1);
    Datetime_vec = DBS_data.(Patient_name).TrendLogs.([hemisphers{hemi_idx} '_Hemi']).Time_vector;

    %
    cols_names = LFP_T.Properties.VariableNames;
    rows_names = LFP_T.Properties.RowNames;

    % Take part of the data
    % LFP_vec = LFP_vec(57*144:67*144);
    % Amp_vec = Amp_vec(57*144:67*144);
    % sensing_vec = sensing_vec(57*144:67*144);
    % contacts_vec = contacts_vec(57*144:67*144);
    % Datetime_vec = Datetime_vec(57*144:67*144);

    % nans
    % nanIdx = isnan(LFP_vec);
    % LFP_vec = LFP_vec(~nanIdx);
    % Amp_vec = Amp_vec(~nanIdx);
    % sensing_vec = sensing_vec(~nanIdx);
    % contacts_vec = contacts_vec(~nanIdx);
    % Datetime_vec = Datetime_vec(~nanIdx);

    % Take part of the data (8 days)
    idx_days_first = 1*144;
    idx_days_last = 23*144;
    LFP_vec_9d = LFP_vec(idx_days_first:idx_days_last);
    Amp_vec_9d = Amp_vec(idx_days_first:idx_days_last);
    sensing_vec_9d = sensing_vec(idx_days_first:idx_days_last);
    contacts_vec_9d = contacts_vec(idx_days_first:idx_days_last);
    Datetime_vec_9d = Datetime_vec(idx_days_first:idx_days_last);

    % Remove NaNs
    nanIdx = isnan(LFP_vec_9d);
    LFP_vec = LFP_vec_9d(~nanIdx);
    Amp_vec = Amp_vec_9d(~nanIdx);
    sensing_vec = sensing_vec_9d(~nanIdx);
    contacts_vec = contacts_vec_9d(~nanIdx);
    Datetime_vec = Datetime_vec_9d(~nanIdx);
    LFP_rand = LFP_vec(randperm(length(LFP_vec)));
    % Replicate the data 4 times as data augmentation technique (total 24
    % days) - help improve the frequency resolution of the spectral analysis
    % LFP_vec = repmat(LFP_vec_9d, 3, 1);
    % Amp_vec = repmat(Amp_vec_9d, 3, 1);
    % sensing_vec = repmat(sensing_vec_9d, 3, 1);
    % contacts_vec = repmat(contacts_vec_9d, 3, 1);
    % Datetime_vec = repmat(Datetime_vec_9d, 3, 1);


    % Events
    Events = DBS_data.(Patient_name).Events.Events_data;
      
    % find unique periods
    combined_strings = strcat(num2str(sensing_vec(:)), '_', contacts_vec(:));
    [unique_pairs, ~, pair_idx] = unique(combined_strings,'stable');
   
%% Plot LFP + Events 
for i = 1:length(unique_pairs)
    current_pair = unique_pairs(i);
    pair_parts = split(unique_pairs(i), '_');
    current_freq = pair_parts(1);
    current_contact1 = pair_parts(2);
    current_contact2 = pair_parts(4);
    LFP_sense = LFP_vec(pair_idx == i);
    sensing_datetime = Datetime_vec(pair_idx == i);

    % Remove outlying dates
    % sensing_datetime_diff = diff(sensing_datetime);
    % mean_datetime_diff = mean(sensing_datetime_diff);
    % leap_indx = find(sensing_datetime_diff > mean_datetime_diff);
    % if ~isempty(leap_indx)
    %     LFP_sense = LFP_sense(leap_indx(1)+1:end);
    %     LFP_sense_shuffled = LFP_sense_shuffled(leap_indx(1)+1:end);
    %     sensing_datetime = sensing_datetime(leap_indx(1)+1:end);
    % end

    sample_interval = abs(hours(sensing_datetime(2) - sensing_datetime(1)));
    sample_freq = 1 / sample_interval;
    time_res = 0.16; % resolution of the time axis in hours
    min_period = max([2 * sample_interval, time_res]);
    max_period = 80;
    time_ax_vec = min_period:time_res:max_period;
    freq_ax_vec = 1 ./ time_ax_vec;

    rec_days = round(days(sensing_datetime(end) - sensing_datetime(1)));
    win_size_days = min(ceil(0.6 * rec_days - 1), 7);
    win_size_bins = win_size_days * 6 * 24;
    win_overlap_days = max(ceil(win_size_days / 2), 1);
    win_overlap_bins = win_overlap_days * 6 * 24;

    % Calculate periodogram for unshuffled data
    [psd_estimate, f_welch] = pwelch(LFP_sense, win_size_bins, win_overlap_bins, freq_ax_vec, sample_freq);

    % Create shuffled data
    LFP_sense_shuffled = LFP_rand(pair_idx == i);

    % Calculate periodogram for shuffled data
    [psd_estimate_shuffled, f_welch_shuffled] = pwelch(LFP_sense_shuffled, win_size_bins, win_overlap_bins, freq_ax_vec, sample_freq);

    % Plot 2x2 figure for each pair
    figure();
    set(gcf, 'Color', 'w')

    % Plot 1: Time plot of unshuffled data
    subplot(2, 2, 1);
    plot(sensing_datetime, LFP_sense, 'LineWidth', 1, 'color', [124 169 165]/255);
    grid off;
    box off;
    set(gca, 'LineWidth', 1.5);
    set(gca, 'FontSize', 10);
    xlabel('Time');
    ylabel('LFP [μV]');
    title('Time Plot - Unshuffled Data');

    % Plot 2: Time plot of shuffled data
    subplot(2, 2, 2);
    plot(sensing_datetime, LFP_sense_shuffled, 'LineWidth', 1, 'color', [124 169 165]/255);
    grid off;
    box off;
    set(gca, 'LineWidth', 1.5);
    set(gca, 'FontSize', 10);
    xlabel('Time');
    ylabel('LFP [μV]');
    title('Time Plot - Shuffled Data');

    % Plot 2: Periodogram of unshuffled data
    subplot(2, 2, [3 4]);
    plot(time_ax_vec, psd_estimate, 'LineWidth', 2.5, 'color', [124 169 165]/255, 'DisplayName', 'PSD of Original LFP');
    hold on
    plot(time_ax_vec, psd_estimate_shuffled, 'LineWidth', 2.5, 'DisplayName', 'PSD of Shuffled LFP');
    grid off;
    box off;
    legend();
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 12);
    xlabel('Cycle Duration [Hours]');
    ylabel('PSD Estimate [μV^2/Hz]');
    title('Periodograms');
    xlim([0 40])
    % seperate plots:
    % Plot 2: Periodogram of unshuffled data
    % subplot(2, 2, 3);
    % plot(time_ax_vec, psd_estimate, 'LineWidth', 2.5, 'color', [124 169 165]/255);
    % grid off;
    % box off;
    % set(gca, 'LineWidth', 2);
    % set(gca, 'FontSize', 12);
    % xlabel('Cycle Duration [Hours]');
    % ylabel('PSD Estimate [μV^2/Hz]');
    % title('Periodogram - Unshuffled Data');
    % 
    % % Plot 4: Periodogram of shuffled data
    % subplot(2, 2, 4);
    % plot(time_ax_vec, psd_estimate_shuffled, 'LineWidth', 2.5, 'color', [124 169 165]/255);
    % grid off;
    % box off;
    % set(gca, 'LineWidth', 2);
    % set(gca, 'FontSize', 12);
    % xlabel('Cycle Duration [Hours]');
    % ylabel('PSD Estimate [μV^2/Hz]');
    % title('Periodogram - Shuffled Data');

    sgtitle({[hemisphers{hemi_idx} ' Hemisphere '], ...
        ['Frequency - ' char(current_freq)] ...
        ['Between contacts ' char(current_contact1) ' and ' char(current_contact2)]});

end


end

% savefig(['Periodogram - Feb-Mar 24 theta (RECENT) - ' Patient_name '.fig'])
