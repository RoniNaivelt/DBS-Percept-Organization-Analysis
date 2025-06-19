function LFP_All_Plots(DBS_data, patient_name, hemi)
% LFP_Daily_Plot -
%
% Inputs:
%   DBS_data            - Nested data structure containing DBS data
%   patient_name        - String or cell array of strings (must match DBS_data entries)
%   hemi                - (Optional) 'Left', 'Right', or 'both' (default)


% Handle default inputs
if nargin < 3 || isempty(hemi)
    hemi = 'both';
end



% Convert to cell arrays if needed
if ~iscell(patient_name)
    patient_name = {patient_name};
end
n_patients = numel(patient_name);



for i_pat = 1:n_patients
    name = patient_name{i_pat};
    % Choose hemispheres
    if strcmpi(hemi, 'both')
        hemis = {'Left', 'Right'};
    else
        hemis = {hemi};
    end

    for i_hemi = 1:length(hemis)
        h_name = hemis{i_hemi};

        % Use OpenDataHelper to extract band_LFP_mat
        try
            data_out = OpenDataHelper(DBS_data, name, h_name);
            key = [name '_' h_name];
            % For LFP / TEED Plotting
            LFP_vec = data_out.(key).LFP_vec_trimmed;
            AMP_vec = data_out.(key).Amp_vec_trimmed;
            sensing_vec = data_out.(key).sensing_vec_trimmed;
            contacts_vec = data_out.(key).contacts_vec_trimmed;
            Stimrate_vec = data_out.(key).Stimrate_vec_trimmed;
            PulseWidth_vec = data_out.(key).PulseWidth_vec_trimmed;
            TEED_vec = data_out.(key).TEED_vec_trimmed;
            days_since_surg = data_out.(key).days_since_surg;
        catch
            warning('Failed to extract LFP data for %s (%s hemisphere). Skipping.', name, h_name);
            continue;
        end


        % find unique periods
        combined_strings = strcat(num2str(sensing_vec(:)), '_', contacts_vec(:));
        % Filter out missing entries first
        valid_idx = ~ismissing(combined_strings);
        clean_strings = combined_strings(valid_idx);
        [unique_pairs, ~, pair_idx_clean] = unique(clean_strings, 'stable');

        % === Draw sensing frequency patches ===
        if numel(unique_pairs) > 1
            clr = turbo(length(unique_pairs));
            for i = 1:length(unique_pairs)
                current_pair = unique_pairs(i);
                parts = split(current_pair, '_');
                if numel(parts) < 4
                    continue;  % skip malformed entries
                end
                current_freq = parts(1);
                contact1 = parts(2);
                contact2 = parts(4);

                idx = (pair_idx_clean == i);
                sensing_times = days_since_surg(idx);
                Current_LFP_vec = LFP_vec(idx);
                Current_AMP_vec = AMP_vec(idx);
                Current_Stimrate_vec = Stimrate_vec(idx);
                Current_PulseWidth_vec = PulseWidth_vec(idx);

                if isempty(sensing_times)
                    continue
                end

                figure('Color', 'w')
                subplot(3,1,[1 2])
                plot(sensing_times, Current_LFP_vec)
                xlabel('Days Since Surgery');
                title({[name ' - ' h_name ' Hemisphere'],...
                    ['Sensing Frequency - ' char(current_freq) ' Hz; Contacts ' char(contact1) '-' char(contact2)]}, 'FontSize',14, 'Interpreter', 'none');
          
                subplot(3,1,3)
                plot(sensing_times, Current_AMP_vec, 'LineWidth', 1.5, 'DisplayName', 'Amplitude')
                hold on
                plot(sensing_times, Current_Stimrate_vec, 'LineWidth', 1.5, 'DisplayName', 'Stimulation Rate')
                plot(sensing_times, Current_PulseWidth_vec, 'LineWidth', 1.5, 'DisplayName', 'Pulse Width')
                xlabel('Days Since Surgery');
                legend()
            
            end
        end
    end
end