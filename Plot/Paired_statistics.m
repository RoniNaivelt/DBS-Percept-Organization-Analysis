
%%  Paired Events 

Patient_name = 'Patient_EC';
hemisphers = {'Right', 'Left'};
Events = DBS_data.(Patient_name).Events.Events_data;

Postmed_indx = find(strcmp(Events.EventName,'Postmed'));
Premed_indx = Postmed_indx - 1;
num_events = length(Premed_indx);
%
Premed_events = Events(Premed_indx,:);
Postmed_events = Events(Postmed_indx,:);
% for counting:
freq = Postmed_events.Left_Snapshot(1,1).FFT_Freq;

for hemi_idx = 1:numel(hemisphers)
    % Comperison Plot 
    Theta_peak_amp_post = [];
    Theta_peak_amp_pre = [];
    Theta_AUC_post = [];
    Theta_AUC_pre = [];
    Beta_peak_amp_post = [];
    Beta_peak_amp_pre = [];
    Beta_AUC_post = [];
    Beta_AUC_pre = [];
    
    for i = 1:size(Premed_events,1)
        fft_pre = Premed_events.([hemisphers{hemi_idx} '_Snapshot'])(i,1).FFT_Val;
        fft_post = Postmed_events.([hemisphers{hemi_idx} '_Snapshot'])(i,1).FFT_Val;    
        %
        Theta_peak_amp_post(end+1) = Postmed_events.([hemisphers{hemi_idx} '_Snapshot'])(i,1).Theta_data.Theta_peak;
        Theta_peak_amp_pre(end+1) = Premed_events.([hemisphers{hemi_idx} '_Snapshot'])(i,1).Theta_data.Theta_peak;
    
        %
        Theta_AUC_post(end+1) = Postmed_events.([hemisphers{hemi_idx} '_Snapshot'])(i,1).Theta_data.Theta_AUC_relative;
        Theta_AUC_pre(end+1) = Premed_events.([hemisphers{hemi_idx} '_Snapshot'])(i,1).Theta_data.Theta_AUC_relative;
        %
        Beta_peak_amp_post(end+1) = Postmed_events.([hemisphers{hemi_idx} '_Snapshot'])(i,1).Beta_data.Beta_peak;
        Beta_peak_amp_pre(end+1) = Premed_events.([hemisphers{hemi_idx} '_Snapshot'])(i,1).Beta_data.Beta_peak;
        %
        Beta_AUC_post(end+1) = Postmed_events.([hemisphers{hemi_idx} '_Snapshot'])(i,1).Beta_data.Beta_AUC_relative;
        Beta_AUC_pre(end+1) = Premed_events.([hemisphers{hemi_idx} '_Snapshot'])(i,1).Beta_data.Beta_AUC_relative;
        %
        Premed.(hemisphers{hemi_idx}).Theta_peak_amp_pre = Theta_peak_amp_pre;
        Premed.(hemisphers{hemi_idx}).Theta_AUC_pre = Theta_AUC_pre;
        Premed.(hemisphers{hemi_idx}).Beta_peak_amp_pre = Beta_peak_amp_pre;
        Premed.(hemisphers{hemi_idx}).Beta_AUC_pre = Beta_AUC_pre;
        Postmed.(hemisphers{hemi_idx}).Theta_peak_amp_post = Theta_peak_amp_post;
        Postmed.(hemisphers{hemi_idx}).Theta_AUC_post = Theta_AUC_post;
        Postmed.(hemisphers{hemi_idx}).Beta_peak_amp_post = Beta_peak_amp_post;
        Postmed.(hemisphers{hemi_idx}).Beta_AUC_post = Beta_AUC_post;

        % figure()
        % sgtitle('PSD Plot of Pre/Post-Medication Events')
        % subplot(2,4,[1 4])
        % plot(freq,fft_pre,'linewidth',2)
        % hold on
        % plot(freq,fft_post,'linewidth',2)
        % title('FFT Plot')
        % xlabel('Frequency [Hz]')
        % ylabel('Power spectral density [uV^2/Hz]')
        % legend('Pre Medication','Post Medication')
        % box off
        % set(gca, 'FontSize', 12)
        % %
        % subplot(2,4,5)
        % bar([Theta_peak_amp_pre(i) Theta_peak_amp_post(i)])
        % xticklabels({'pre', 'post'})
        % title('Theta peak amplitude')
        % ylim([0 2])
        % box off
        % set(gca, 'FontSize', 12)
        % %
        % subplot(2,4,6)
        % bar([Theta_AUC_pre(i) Theta_AUC_post(i)])
        % xticklabels({'pre', 'post'})
        % title('Theta AUC')
        % ylim([0 0.3])
        % box off
        % set(gca, 'FontSize', 12)
        % %
        % subplot(2,4,7)
        % bar([Beta_peak_amp_pre(i) Beta_peak_amp_post(i)])
        % xticklabels({'pre', 'post'})
        % title('Beta peak amplitude')
        % ylim([0 2])
        % box off
        % set(gca, 'FontSize', 12)
        % %
        % subplot(2,4,8)
        % bar([Beta_AUC_pre(i) Beta_AUC_post(i)])
        % xticklabels({'pre', 'post'})
        % title('Beta AUC')
        % ylim([0 0.3])
        % box off
        % set(gca, 'FontSize', 12)
        % 
        % set(gcf, 'Color', 'w')
    end
end
%%
figure()
bar([mean(Premed.Left.Theta_peak_amp_pre) mean(Postmed.Left.Theta_peak_amp_post);...
    mean(Premed.Left.Theta_AUC_pre) mean(Postmed.Left.Theta_AUC_post);... 
    mean(Premed.Left.Beta_peak_amp_pre) mean(Postmed.Left.Beta_peak_amp_post);
    mean(Premed.Left.Beta_AUC_pre) mean(Postmed.Left.Beta_AUC_post)])
xticklabels({'Theta peak amplitude', 'Theta AUC', 'Beta peak amplitude', 'Beta AUC'})
ylabel('count post<pre')
title('fraction of post < pre events (n=16)')
legend('right', 'left')

 %%%%% t tests %%%%%%%

% paired pre - post (EC)

[h_beta_peak_amp_left, p_beta_peak_amp_left] = ttest(Premed.Left.Theta_peak_amp_pre,Postmed.Left.Theta_peak_amp_post);
[h_beta_peak_amp_right, p_beta_peak_amp_right] = ttest(Premed.Left.Theta_AUC_pre,Postmed.Left.Theta_AUC_post);

[h_beta_AUC_left, p_beta_AUC_left] = ttest(Premed.Left.Beta_peak_amp_pre,Postmed.Left.Beta_peak_amp_post);
[h_beta_AUC_right, p_beta_AUC_right] = ttest(Premed.Left.Beta_AUC_pre,Postmed.Left.Beta_AUC_post);


%% Paired LFP - Medications
% Step 1 - combine all 'med' or 'premed' events
hemi_labels = {'Right_Hemi', 'Left_Hemi'};
hemisphers = {'Right', 'Left'};

patients = {'AL'}; %{'EC', 'AD', 'YB', 'AL'};
right_lfp_mean_before = [];
right_lfp_mean_after = [];
left_lfp_mean_before = [];
left_lfp_mean_after = [];
left_lfp_around = [];
right_lfp_around = [];
right_lfp_mean_long_before = [];
right_lfp_mean_long_after = [];
left_lfp_mean_long_before = [];
left_lfp_mean_long_after = [];
right_lfp_rand =[];
left_lfp_rand = [];


for i = 1 % :length(patients)    %- for multiple patients
    Patient_name = ['Patient_' patients{i}];
    Event_table = Convert_Event_to_feature(DBS_data, Patient_name, 'Min_dopa', [0 0]); % for med - Min_dopa
    event_vec = reshape(table2array(Event_table), [], 1);
    event_idx = find(event_vec == 1);
    for hemi_idx = 1:numel(hemisphers)
        hemi_label = hemi_labels{hemi_idx};
        hemi = hemisphers{hemi_idx};
        
        LFP_vec = normalize(reshape(table2array(DBS_data.(Patient_name).TrendLogs.(hemi_label).LFP_table), [], 1),'range');
        Sensing_Freq = reshape(table2array(DBS_data.(Patient_name).Groups.(hemi_label).Sensing_Freq), [], 1);
        LFP_vec = LFP_vec(Sensing_Freq > 13 & Sensing_Freq < 30); % ONLY BETA SENSING
        event_vec_hemi = event_vec(Sensing_Freq > 13 & Sensing_Freq < 30);
        event_idx_hemi = find(event_vec_hemi == 1);
        event_idx_hemi = event_idx_hemi(2:end-1);
        
        % Create random vector indices
        all_rand_idx = randperm(length(LFP_vec)-13);
        rand_idx_hemi = [];       
        for i = 1:length(event_idx_hemi)
            exclusion_range = (event_idx_hemi(i)-12):(event_idx_hemi(i)+12);
            all_rand_idx(ismember(all_rand_idx, exclusion_range)) = [];
            all_rand_idx(all_rand_idx<13) = [];
        end
        rand_idx_hemi = all_rand_idx(1:length(event_idx_hemi));


        for idx = 1:length(event_idx_hemi)
            if strcmp(hemi, 'Right')
                %right_lfp_mean_before(end+1) = mean(LFP_vec(event_idx_hemi(idx)-3:event_idx_hemi(idx)-1));
                %right_lfp_mean_after(end+1) = mean(LFP_vec(event_idx_hemi(idx)+6:event_idx_hemi(idx)+8));

                right_lfp_mean_before(end+1) = mean(LFP_vec(event_idx_hemi(idx)-3:event_idx_hemi(idx)-1));
                %right_lfp_mean_long_before(end+1) = mean(LFP_vec(event_idx_hemi(idx)-23:event_idx_hemi(idx)-4));
                right_lfp_mean_after(end+1) = mean(LFP_vec(event_idx_hemi(idx)+6:event_idx_hemi(idx)+8));
                right_lfp_mean_long_after(end+1) = mean(LFP_vec(event_idx_hemi(idx)+9:event_idx_hemi(idx)+23));
                right_lfp_around = [right_lfp_around ; LFP_vec(event_idx_hemi(idx)-12:event_idx_hemi(idx)+12)'];
                right_lfp_rand = [right_lfp_rand ; LFP_vec(rand_idx_hemi(idx)-12:rand_idx_hemi(idx)+12)'];
            else
                %left_lfp_mean_before(end+1) = mean(LFP_vec(event_idx_hemi(idx)-3:event_idx_hemi(idx)-1));
                %left_lfp_mean_after(end+1) = mean(LFP_vec(event_idx_hemi(idx)+6:event_idx_hemi(idx)+8));
                
                left_lfp_mean_before(end+1) = mean(LFP_vec(event_idx_hemi(idx)-3:event_idx_hemi(idx)-1));
                %left_lfp_mean_long_before(end+1) = mean(LFP_vec(event_idx_hemi(idx)-23:event_idx_hemi(idx)-4));
                left_lfp_mean_after(end+1) = mean(LFP_vec(event_idx_hemi(idx)+6:event_idx_hemi(idx)+8));
                left_lfp_mean_long_after(end+1) = mean(LFP_vec(event_idx_hemi(idx)+9:event_idx_hemi(idx)+23));
                left_lfp_around = [left_lfp_around ; LFP_vec(event_idx_hemi(idx)-12:event_idx_hemi(idx)+12)'];
                left_lfp_rand = [left_lfp_rand ; LFP_vec(rand_idx_hemi(idx)-12:rand_idx_hemi(idx)+12)'];

            end
        end
        
        if strcmp(hemi, 'Right')
            [p_right_lfp_s, h_right_lfp_s] = signrank(right_lfp_mean_before, right_lfp_mean_after);
            results.(Patient_name).RightLFP_ttest = h_right_lfp_s;
            results.(Patient_name).RightLFP_pvalue = p_right_lfp_s;
            results.(Patient_name).RightEvents = length(event_idx_hemi);
        else
            [p_left_lfp_s, h_left_lfp_s] = signrank(left_lfp_mean_before, left_lfp_mean_after);
            results.(Patient_name).LeftLFP_ttest = h_left_lfp_s;
            results.(Patient_name).LeftLFP_pvalue = p_left_lfp_s;
            results.(Patient_name).LeftEvents = length(event_idx_hemi);
        end
    end
end

save('medication_paired_pvals.mat', 'results')
%%%%% final t-test %%%%
[p_right_lfp, h_right_lfp] = signrank(right_lfp_mean_before, right_lfp_mean_after);
[p_left_lfp, h_left_lfp] = signrank(left_lfp_mean_before, left_lfp_mean_after);

% to plot several patients:
lfp_mean_before_struct.([Patient_name '_right']) = right_lfp_mean_before;
lfp_mean_after_struct.([Patient_name '_right']) = right_lfp_mean_after;
lfp_mean_long_after_struct.([Patient_name '_right']) = right_lfp_mean_long_after;
lfp_around_struct.([Patient_name '_right']) = right_lfp_around;
lfp_rand_struct.([Patient_name '_right']) = right_lfp_rand;

lfp_mean_before_struct.([Patient_name '_left']) = left_lfp_mean_before;
lfp_mean_after_struct.([Patient_name '_left']) = left_lfp_mean_after;
lfp_mean_long_after_struct.([Patient_name '_left']) = left_lfp_mean_long_after;
lfp_around_struct.([Patient_name '_left']) = left_lfp_around;
lfp_rand_struct.([Patient_name '_left']) = left_lfp_rand;


%% plot LFP around medication events

lfp_mean_around = [];
lfp_mean_rand = [];
fields = fieldnames(lfp_around_struct);
for i = 1:length(fields)
    lfp_mean_around = [lfp_mean_around ; lfp_around_struct.(fields{i})];
    lfp_mean_rand = [lfp_mean_rand ; lfp_rand_struct.(fields{i})];
end


% boxchart plot:
figure()
set(gcf, 'Color', 'w')
b = boxchart(lfp_mean_around,'Notch','on');
b.JitterOutliers = 'on';
b.MarkerStyle = '.';

hold on
plot(mean(left_lfp_around, 1), '-o', 'Color', 'r', 'LineWidth', 2.5)
xticklabels(-120:10:120)
ylabel('LFP Amplitude [\muV]', 'FontSize',14)
xlabel('Time around medication intake [minutes]', 'FontSize',14)
title('LFP around medication events', 'FontSize', 18)


% lines with mean plot:
figure()
set(gcf, 'Color', 'w')
colors = summer(size(lfp_mean_around, 1)); % Generate colormap
% for i = 1:size(lfp_mean_around, 1)
%     plot(1:25, lfp_mean_around(i,:), 'LineWidth', 1, 'Color', [colors(i,:), 0.2]); % 0.8 for opacity
%     hold on
% end
ylabel('LFP Amplitude [\muV]', 'FontSize',14)
xlabel('Time around medication intake [minutes]', 'FontSize',14)
title('LFP around medication events', 'FontSize', 18)
plot(mean(lfp_mean_around, 1, 'omitnan'), '-ok', 'LineWidth', 2.5)
colormap(summer(size(lfp_mean_around, 1)))
set(gca, 'LineWidth',2)
box off
xticks(1:25)
xticklabels(-120:10:120)
xlim([1 25])


% plot mse + correlation to mean
mean_med_vec = mean(lfp_mean_around, 1, 'omitnan');
corr_vec_med = [];
corr_vec_rand = [];
mse_vec_med = [];
mse_vec_rand = [];
% % % MI_vec_med = [];
% % % MI_vec_rand = [];
for i = 1:size(lfp_mean_around, 1)
    LFP_mean_curr = lfp_mean_around(i,:);
    LFP_rand_curr = lfp_mean_rand(i,:);
    % Calculate the Spearman correlation
    corr_vec_med(end+1) = corr(mean_med_vec(:), LFP_mean_curr(:), 'Type', 'Spearman');
    corr_vec_rand(end+1) = corr(mean_med_vec(:), LFP_rand_curr(:), 'Type', 'Spearman');

    % Calculate the Mean Squared Error (MSE)
    diff = lfp_mean_around(i,:) - mean_med_vec;
    mse_vec_med(i) = mean(diff.^2);
   
    diff = lfp_mean_rand(i,:) - mean_med_vec;
    mse_vec_rand(i) = mean(diff.^2);

    % % % % % MI
    % % % % lfp_mean_around_discrt = discretize(lfp_mean_around, 5);
    % % % % lfp_mean_rand_discrt = discretize(lfp_mean_rand, 5);
    % % % % mean_med_discrt = discretize(mean_med_vec, 5);
    % % % % LFP_mean_curr_discrt = lfp_mean_around_discrt(i,:);
    % % % % lfp_rand_discrt = lfp_mean_rand_discrt(i,:);
    % % % % MI_med = mutualinfo([mean_med_discrt; LFP_mean_curr_discrt], (1:2)');
    % % % % MI_rand = mutualinfo([mean_med_discrt; lfp_rand_discrt], (1:2)');
    % % % % MI_vec_med(end+1) = MI_med;
    % % % % MI_vec_rand(end+1) = MI_rand;
end
mse_vec_rand(isnan(mse_vec_rand)) = 1;

figure()
hold on; % Hold on to plot multiple box charts
b1 = boxchart(ones(size(MI_vec_med)), MI_vec_med);
b2 = boxchart(2*ones(size(MI_vec_rand)), MI_vec_rand);
% daboxplot([corr_vec_med', corr_vec_rand'], 'xtlabels', {'Med', 'Random'}, 'scatter', 2, 'fill', 0);
% daboxplot([mse_vec_med', mse_vec_rand'], 'xtlabels', {'Med', 'Random'}, 'scatter', 2, 'fill', 0);

b1.BoxFaceColor = [124 169 165]/255; % Color for the first box
set(gca, 'XTick', [1 2], 'XTickLabel', {'Med', 'Rand'}); % Label X-axis
ylabel('MSE', 'FontSize', 14);
title('MSE of medication intake LFPs to the mean', 'FontSize', 14);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 14);
box off;set(gca, 'LineWidth',2)
box off
set(gcf, 'Color', 'w')

[p_corr, h_corr] = signrank(corr_vec_med, corr_vec_rand);
[p_mse, h_mse] = signrank(mse_vec_med, mse_vec_rand);

%% periodogram
% shuffle rows
row_ind_shuffle = randperm(size(lfp_mean_around,1));
lfp_mean_around = lfp_mean_around(row_ind_shuffle,:);

lfp_med_combined = [];
for i = 1:length(lfp_mean_around)
    lfp_vec = lfp_mean_around(i,:);
    lfp_shuffpled = lfp_vec(randperm(25));
    lfp_med_combined = [lfp_med_combined , lfp_vec];
end

Datetime_vec = 1:length(lfp_med_combined);

sample_freq = 25;
min_period = 1;
max_period = 50;
time_ax_vec = 1:1:50;
freq_ax_vec = 1 ./ time_ax_vec;

rec_days = Datetime_vec(end)/25;
win_size_cicles = min(ceil(0.6 * rec_days - 1), 7);
win_size_bins = win_size_cicles * 25;
win_overlap_cicles = max(ceil(win_size_cicles / 2), 1);
win_overlap_bins = win_overlap_cicles * 25;

% Calculate periodogram
[psd_estimate, f_welch] = pwelch(lfp_med_combined, win_size_bins, win_overlap_bins, freq_ax_vec, sample_freq);

% Plot Periodogram 
figure()
plot(time_ax_vec, psd_estimate, 'LineWidth', 2.5, 'color', [124 169 165]/255);
grid off;
box off;
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 12);
xlabel('Cycle Duration [Hours]');
ylabel('PSD Estimate [μV^2/Hz]');
title('Periodogram');
set(gcf, 'Color', 'w')


%% Recreating a medication model
% Exponential rise parameters
[max_value, peak_index] = max(mean_med_vec);
min_value = min(mean_med_vec);
mean_value = mean(mean_med_vec);
line_length = length(mean_med_vec);

a = max_value;
b = log(min_value / a) / (1 - peak_index); % Adjust this based on mean value

% Logarithmic decrease parameters
c = max_value;
d = (45 - peak_index) / log(c);  % Adjust this to fit the decay

% Generate the rising and falling parts
x_rise = 1:peak_index;
y_rise = a * exp(b * (x_rise - peak_index));

x_fall = peak_index:line_length;
y_fall = mean_value + (c - mean_value) * (1 - log(d * (x_fall - peak_index) + 1) / log(d * (line_length - peak_index) + 1));

% Combine both parts
x_full = [x_rise, x_fall(2:end)]; % Avoid duplicating the peak point
y_full = [y_rise, y_fall(2:end)];

% Plot the result
figure;
plot(x_full, y_full, 'LineWidth', 2);
hold on
plot(mean_med_vec)
xlabel('Index');
ylabel('Value');
title('Manual Curve Fitting: Exponential Rise and Logarithmic Decrease');
set(gca, 'FontSize', 12);
grid on;



%% per-hemispher individuals plot
figure()
set(gcf, 'Color', 'w')
subplot(2,1,1)
colors = turbo(size(left_lfp_around, 1)); % Generate colormap
for i = 1:size(left_lfp_around, 1)
    plot(left_lfp_around(i,:), 'LineWidth', 1, 'Color', colors(i,:));
    hold on
end
ylabel('LFP Amplitude [\muV]', 'FontSize',14)
xlabel('Time around medication intake [minutes]', 'FontSize',14)
title('Left Hemisphere', 'FontSize', 16)
xticks(1:21)
xticklabels(-100:10:100)
xlim([1 21])
plot(mean(left_lfp_around, 1), 'k', 'LineWidth', 2.5)
colormap(turbo(size(left_lfp_around, 1)))
colorbar('Ticks', [0, 1], 'TickLabels', {'Oldest', 'Newest'}, 'FontSize',14)

subplot(2,1,2)
colors = turbo(size(right_lfp_around, 1)); % Generate colormap
for i = 1:size(right_lfp_around, 1)
    plot(right_lfp_around(i,:), 'LineWidth', 1, 'Color', colors(i,:));
    hold on
end
ylabel('LFP Amplitude [\muV]', 'FontSize',14)
xlabel('Time around medication intake [minutes]', 'FontSize',14)
title('Right Hemisphere', 'FontSize', 16)
xticks(1:21)
xticklabels(-100:10:100)
xlim([1 21])
plot(mean(right_lfp_around, 1), 'k', 'LineWidth', 2.5)
colormap(turbo(size(left_lfp_around, 1)))
colorbar('Ticks', [0, 1], 'TickLabels', {'Oldest', 'Newest'}, 'FontSize',14)
sgtitle('LFP around medication events', 'FontSize', 18)

% check var
var_before_left = mean(peak2peak(left_lfp_around(:,1:10),2)');
var_after_left = mean(peak2peak(left_lfp_around(:,12:21),2)');

var_before_right = mean(peak2peak(right_lfp_around(:,1:10),2)');
var_after_right = mean(peak2peak(right_lfp_around(:,12:21),2)');

%% connected dot plot around medication
figure()
subplot(1,2,1)
daboxplot([right_lfp_mean_before', right_lfp_mean_after', right_lfp_mean_long_after'], 'xtlabels', {'Before', '1hr After', '1-4 Hr After'}, 'withinlines', 1, 'scatter', 2, 'fill', 0);
[p_right_lfp_s, h_right_lfp_s] = signrank(right_lfp_mean_before, right_lfp_mean_after);
title('Right Hemisphere');
xlabel('Condition');
ylabel('LFP Mean Amplitude [µV]');
box off
set(gca, 'LineWidth', 3)
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)
% set(gcf, 'color', [124 169 165]/255);    
% set(gca, 'color', [124 169 165]/255);

subplot(1,2,2)
daboxplot([left_lfp_mean_before', left_lfp_mean_after', left_lfp_mean_long_after'], 'xtlabels', {'Before', '1hr After', '1-4 Hr After'}, 'withinlines', 1, 'scatter', 2, 'fill', 0);
% 'xtlabels', {'Before', 'After'}
[p_left_lfp_s, h_left_lfp_s] = signrank(left_lfp_mean_before, left_lfp_mean_after);
title('Left Hemisphere');
xlabel('Condition');
ylabel('LFP Mean Amplitude [µV]');
box off
set(gca, 'LineWidth', 3)
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)
sgtitle('LFP Means arround Medication intake', 'FontSize', 18)

% set(gcf, 'color', [124 169 165]/255);    
% set(gca, 'color', [124 169 165]/255);
savefig(['LFP Mean Amplitude' num2str(string(patients)) '.fig'])
% export_fig(['LFP Mean Amplitude' num2str(string(patients))], '-png')
% exportgraphics(gcf,'LFP Mean Amplitude.png')




lfp_mean_before = [];
lfp_mean_after = [];
lfp_mean_long_after = [];
fields = fieldnames(lfp_mean_before_struct);
for i = 1:length(fields)
    lfp_mean_before = [lfp_mean_before lfp_mean_before_struct.(fields{i})];
    lfp_mean_after = [lfp_mean_after lfp_mean_after_struct.(fields{i})];
    lfp_mean_long_after = [lfp_mean_long_after lfp_mean_long_after_struct.(fields{i})];
end

figure()
daboxplot([lfp_mean_before', lfp_mean_after', lfp_mean_long_after'], 'xtlabels', {'Before', '1hr After', '2-4 Hr After'}, 'withinlines', 1, 'scatter', 2, 'fill', 0);
[p_right_lfp_s, h_right_lfp_s] = signrank(lfp_mean_before, lfp_mean_after);
title('LFP Plot Around Medication Intake');
xlabel('Condition');
ylabel('LFP Mean Amplitude [µV]');
box off
ylim([0 1])
set(gca, 'LineWidth', 3)
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)


%% final Wilcoxon test %%%

p_w_right_lfp = signrank(right_lfp_mean_before,right_lfp_mean_after);
p_w_left_lfp = signrank(left_lfp_mean_before,left_lfp_mean_after);

%% Paired LFP - Tremor

right_lfp_mean_tremor = [];
right_lfp_mean_normal = [];
left_lfp_mean_tremor = [];
left_lfp_mean_normal = [];

Patient_name = 'Patient_EC';
Event_table = Convert_Event_to_feature(DBS_data,Patient_name, 'Tremor', [0 0]);
event_vec = reshape(table2array(Event_table),[],1);
event_idx = find(event_vec == 1);

% RIGHT
LFP_right_vec = reshape(table2array(DBS_data.(Patient_name).TrendLogs.Right_Hemi.LFP_table),[],1);
Sensing_Freq_right = reshape(table2array(DBS_data.(Patient_name).Groups.Right_Hemi.Sensing_Freq),[],1);
LFP_right_vec = LFP_right_vec(Sensing_Freq_right > 5 & Sensing_Freq_right < 30);
event_vec_right = event_vec(Sensing_Freq_right> 5 & Sensing_Freq_right < 30);
event_idx_r = find(event_vec_right == 1);
%event_idx_r = event_idx_r(2:end-1);

for idx = 1:length(event_idx_r)
    right_lfp_mean_tremor(end+1) = mean(LFP_right_vec(event_idx_r(idx)-3:event_idx_r(idx)));
    right_lfp_mean_normal(end+1) = mean(LFP_right_vec(event_idx_r(idx)+6:event_idx_r(idx)+9));
end

% LEFT
LFP_left_vec = reshape(table2array(DBS_data.(Patient_name).TrendLogs.Left_Hemi.LFP_table),[],1);
Sensing_Freq_left = reshape(table2array(DBS_data.(Patient_name).Groups.Left_Hemi.Sensing_Freq),[],1);
LFP_left_vec = LFP_left_vec(Sensing_Freq_left>5 & Sensing_Freq_left < 30);
event_vec_left = event_vec(Sensing_Freq_left> 5 & Sensing_Freq_left < 30);
event_idx_l = find(event_vec_left == 1);
%event_idx_l = event_idx_l(2:end-1);

for idx = 1:length(event_idx_l)
    left_lfp_mean_tremor(end+1) = mean(LFP_left_vec(event_idx_l(idx)-3:event_idx_l(idx)));
    left_lfp_mean_normal(end+1) = mean(LFP_left_vec(event_idx_l(idx)+6:event_idx_l(idx)+9));
end

%%%%% final t-tests %%%%%%%

[h_right_trem,p_right_trem] = ttest(right_lfp_mean_tremor,right_lfp_mean_normal);
[h_left_trem,p_left_trem] = ttest(left_lfp_mean_tremor,left_lfp_mean_normal);




