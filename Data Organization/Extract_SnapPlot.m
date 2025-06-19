function [Theta_data, Alpha_data, Beta_data, Gamma_data] = Extract_SnapPlot(SnapPlot_data)

% vectors:
FFT = SnapPlot_data.FFT_Val;
Freq = SnapPlot_data.FFT_Freq;
AUC_tot = trapz(FFT);

% bands:
Theta_freq_idx = find(Freq > 4 & Freq < 8);
Alpha_freq_idx = find(Freq > 8 & Freq < 13);
Beta_freq_idx = find(Freq > 14 & Freq < 30);
Gamma_freq_idx = find(Freq > 30 & Freq < 60);
        
%%%% extract band data:
% Theta
Theta_fft = FFT(Theta_freq_idx);
[Theta_data.Theta_peak, Theta_data.Theta_peak_indx] = max(Theta_fft);
Theta_data.Theta_peak_freq = Freq(Theta_freq_idx(Theta_data.Theta_peak_indx));
Theta_AUC = trapz(Freq(Theta_freq_idx),FFT(Theta_freq_idx));
Theta_data.Theta_AUC_relative = Theta_AUC/AUC_tot;

% Alpha
Alpha_fft = FFT(Alpha_freq_idx);
[Alpha_data.Alpha_peak, Alpha_data.Alpha_peak_indx] = max(Alpha_fft);
Alpha_data.Alpha_peak_freq = Freq(Alpha_freq_idx(Alpha_data.Alpha_peak_indx));
Alpha_AUC = trapz(Freq(Alpha_freq_idx),FFT(Alpha_freq_idx));
Alpha_data.Alpha_AUC_relative = Alpha_AUC/AUC_tot;

% Beta
Beta_fft = FFT(Beta_freq_idx);
[Beta_data.Beta_peak, Beta_data.Beta_peak_indx] = max(Beta_fft);
Beta_data.Beta_peak_freq = Freq(Beta_freq_idx(Beta_data.Beta_peak_indx));
Beta_AUC = trapz(Freq(Beta_freq_idx),FFT(Beta_freq_idx));
Beta_data.Beta_AUC_relative = Beta_AUC/AUC_tot;

% Gamma
Gamma_fft = FFT(Gamma_freq_idx);
[Gamma_data.Gamma_peak, Gamma_data.Gamma_peak_indx] = max(Gamma_fft);
Gamma_data.Gamma_peak_freq = Freq(Gamma_freq_idx(Gamma_data.Gamma_peak_indx));
Gamma_AUC = trapz(Freq(Gamma_freq_idx),FFT(Gamma_freq_idx));
Gamma_data.Gamma_AUC_relative = Gamma_AUC/AUC_tot;


end

