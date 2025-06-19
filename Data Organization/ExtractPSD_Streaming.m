function [pxx_normalized, f] = ExtractPSD_Streaming(TimeDomainData, SampleRateInHz)


% 1. Bandpass Filter (2-130 Hz)
bpFilt = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', 2, 'HalfPowerFrequency2', 125, ...
    'SampleRate', SampleRateInHz);
filtered_data = filtfilt(bpFilt, TimeDomainData);
filtered_data = filtered_data(551:end-550);


% 2. Outlier/Artifact Detection
% 2.1 Compute z-score
% 2.2 Flag samples beyond 3 standard deviations
% 2.3 Remove artifacts by interpolation or NaN substitution
zScores = (filtered_data - mean(filtered_data)) / std(filtered_data);
artifact_idx = abs(zScores) > 3; 
filtered_data(artifact_idx) = interp1(find(~artifact_idx), filtered_data(~artifact_idx), find(artifact_idx), 'linear', 'extrap');


% 3. Power Spectral Density Estimation using pwelch
freq_range = 0:0.5:(SampleRateInHz/2); % 0.5 Hz resolution
[pxx, f] = pwelch(filtered_data, [], [], freq_range, SampleRateInHz);


% 4. Normalization of PSD
totalPower = sum(pxx);
if totalPower > 0
    pxx_normalized = pxx / totalPower;
else
    error('Total power is zero, check preprocessing steps.');
end


end
