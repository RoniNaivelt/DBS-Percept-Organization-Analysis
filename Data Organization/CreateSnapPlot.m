function [Theta_data,Beta_data, SnapPlot] = CreateSnapPlot(data, event_num)
%%%%  OLD FUNCTION - USED TO BE INSIDE Create_Events


subject_name = [data.PatientInformation.Initial.PatientFirstName(1) '.' data.PatientInformation.Initial.PatientLastName(1)];
session_date = data.SessionDate(1:10);
%
events = data.DiagnosticData.LfpFrequencySnapshotEvents;
try 
   event_data = events{event_num};
catch
    event_data = events(event_num);   
end

freq = event_data.LfpFrequencySnapshotEvents.HemisphereLocationDef_Left.Frequency;

beta_freq_idx = find(freq > 14 & freq < 30);
theta_freq_idx = find(freq > 5 & freq < 9);

        
event_date = event_data.DateTime;
event_type = event_data.EventName;

% for each hemisphere we find the beta peak location (frequency) and the
% peak's amplitude and save them in a list.

%% RIGHT
right_hemi = event_data.LfpFrequencySnapshotEvents.HemisphereLocationDef_Right;
FFT_right = right_hemi.FFTBinData;
beta_fft_right = FFT_right(beta_freq_idx);
theta_fft_right = FFT_right(theta_freq_idx);

[peak_right_b, peak_idx_right_b] = max(beta_fft_right);
[peak_right_t, peak_idx_right_t] = max(theta_fft_right);

peak_freq_right_b = freq(beta_freq_idx(peak_idx_right_b));
peak_freq_right_t = freq(theta_freq_idx(peak_idx_right_t));

AUC_right_b = trapz(freq(beta_freq_idx),FFT_right(beta_freq_idx));
AUC_right_t = trapz(freq(theta_freq_idx),FFT_right(theta_freq_idx));
AUC_tot = trapz(FFT_right);

AUC_right_b_relative = AUC_right_b/AUC_tot;
AUC_right_t_relative = AUC_right_t/AUC_tot;

% saving
Theta_data.right_hemi.Peak_Amplitude = peak_right_t;
Theta_data.right_hemi.Peak_Frequency = peak_freq_right_t;
Theta_data.right_hemi.AUC = AUC_right_t;
Theta_data.right_hemi.relative_AUC = AUC_right_t_relative;

Beta_data.right_hemi.Peak_Amplitude = peak_right_b;
Beta_data.right_hemi.Peak_Frequency = peak_freq_right_b;
Beta_data.right_hemi.AUC = AUC_right_b;
Beta_data.right_hemi.relative_AUC = AUC_right_b_relative;


%% LEFT
left_hemi = event_data.LfpFrequencySnapshotEvents.HemisphereLocationDef_Left;
FFT_left = left_hemi.FFTBinData;
beta_fft_left = FFT_left(beta_freq_idx);
theta_fft_left = FFT_left(theta_freq_idx);

[peak_left_b, peak_idx_left_b] = max(beta_fft_left);
[peak_left_t, peak_idx_left_t] = max(theta_fft_left);

peak_freq_left_b = freq(beta_freq_idx(peak_idx_left_b));
peak_freq_left_t = freq(theta_freq_idx(peak_idx_left_t));

AUC_left_b = trapz(freq(beta_freq_idx),FFT_left(beta_freq_idx));
AUC_left_t = trapz(freq(theta_freq_idx),FFT_left(theta_freq_idx));
AUC_tot = trapz(FFT_left);

AUC_left_b_relative = AUC_left_b/AUC_tot;
AUC_left_t_relative = AUC_left_t/AUC_tot;

% saving
Theta_data.left_hemi.Peak_Amplitude = peak_left_t;
Theta_data.left_hemi.Peak_Frequency = peak_freq_left_t;
Theta_data.left_hemi.AUC = AUC_left_t;
Theta_data.left_hemi.relative_AUC = AUC_left_t_relative;

Beta_data.left_hemi.Peak_Amplitude = peak_left_b;
Beta_data.left_hemi.Peak_Frequency = peak_freq_left_b;
Beta_data.left_hemi.AUC = AUC_left_b;
Beta_data.left_hemi.relative_AUC = AUC_left_b_relative;


% PLOT
SnapPlot = figure('Visible', 'off');
sgtitle({['Subject - ' subject_name],['Session date: ' session_date],['Event: ' event_type ' ; Time: ' event_date(1:10) ' , ' event_date(end-8:end-4)]});

subplot(2,1,1)
plot(freq,FFT_right,'linewidth',2)
title('Right Hemisphere')
xlabel('Frequency [Hz]')
ylabel('Power spectral density [uVP/Hz]')
xline(peak_freq_right_b, '-', {['Beta Peak at ' num2str(peak_freq_right_b)]}, 'linewidth',2)
xline(peak_freq_right_t, '-', {['Theta Peak at ' num2str(peak_freq_right_t)]}, 'linewidth',2)

subplot(2,1,2)
plot(freq,FFT_left,'linewidth',2)
title('Left Hemisphere')
xlabel('Frequency [Hz]')
ylabel('Power spectral density [uVP/Hz]')
xline(peak_freq_left_b, '-', {['Beta Peak at  ' num2str(peak_freq_left_b)]}, 'linewidth',2);
xline(peak_freq_left_t, '-', {['Theta Peak at ' num2str(peak_freq_left_t)]}, 'linewidth',2);



end

