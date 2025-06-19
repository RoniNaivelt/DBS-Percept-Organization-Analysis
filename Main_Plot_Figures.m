%% ========================================================================
% This file contains scripts and helper documentation for plotting various
% figures derived from the DBS_data struct.

% Each plotting function is modular and accepts key inputs:
%   - the dataset (DBS_data)
%   - Patient names (given as their initials, like 'Patient_RN')
%   - Optional hemisphere or frequency range specifications
% some functions support additional name-value parameters for customizing
% plotting behavior or filtering the data.



% Initialization
load('DBS_data.mat');
addpath('Plot Functions/');
addpath('Plot Functions/helper functions/');


%% ========================================================================
% Function 1: Diurnal_LFP_Plot
%% ========================================================================
%   This function generates diurnal polar plots of LFP (Local Field Potential)
%   power across time for one or more patients and hemispheres. Each plot represents
%   a 24-hour LFP cycle, optionally averaged across selected days or patients.
%
% Inputs:
%   - DBS_data: structured input with LFP data for each patient and hemisphere.
%   - patient_names: string or cell array of patient name(s).
%   - hemi (optional): 'Left', 'Right', or 'both' (default).
%   - sensing_range (optional): e.g., [13 30] for beta band.
%   - day_range (optional): [start end] or cell array with ranges per patient.
%
% Name-Value Pairs:
%   - 'PlotMean': true/false (plot an average curve across all patients).
%   - 'DiurnalPlots': struct to continue adding plots.
%   - 'colors': color spec for each patient.
%
% Usage Example:
patient_names = {'Patient_RN', 'Patient_NR'};
hemisphere = 'Left';
sensing_range = [13 30];
day_range = {[1 10], [30 50]};
Diurnal_LFP_Plot(DBS_data, patient_names, hemisphere, sensing_range, day_range, 'PlotMean', true);

% Notes:
%   - This function accepts inputting several patients as a cell array, and
%   accordingly day_range can also be a cell array.
%   - Outputs a polar plot for each patient/hemisphere.
%   - Handles optional color assignment and overlay of mean curve.


%% ========================================================================
% Function 2: Heatmap_Plot
%% ========================================================================
%   This function creates a heatmap of LFP power over time for a specific
%   patient and hemisphere. Each column represents a different day, and
%   each row corresponds to a specific time of day, providing a clear
%   overview of temporal LFP dynamics.
%
% Inputs:
%   - DBS_data: structured input containing LFP data for the patient.
%   - patient_name: string specifying the patient's name.
%   - hemi: 'Left' or 'Right' hemisphere.
%   - day_range (optional): [start_day, end_day] to restrict displayed data.
%
% Name-Value Pairs:
%   - 'colors': custom color map (optional, defaults to 'turbo').
%
% Usage Example:
patient_name = 'Patient_RN';
hemi = 'Right';
day_range = [1 20];
Heatmap_Plot(DBS_data, patient_name, hemi, day_range);

% Notes:
%   - The X-axis shows days since surgery, based on the implant date.
%   - The Y-axis shows time of day using the LFP sample timestamps.
%   - No band selection is applied here; the full spectrum values are shown.


%% =========================================================================
% Function 3: Histogram_Plot
%% =========================================================================
%   This function creates histograms comparing LFP amplitude distributions:
%   - Left vs Right hemisphere
%   - Day vs Night for each hemisphere
%
% Inputs:
%   - DBS_data: structured input with LFP data.
%   - patient_name: string, patient ID.
%   - sensing_range (optional): frequency band, e.g., [13 30].
%   - day_range (optional): [start end] or cell array for each hemi.
%
% Usage:
patient_name = 'Patient_RN';
sensing_range = [13 30];
day_range = [10 40];
Histogram_Plot(DBS_data, patient_name, sensing_range, day_range);

%% =========================================================================
% Function 4: LFP_All_Plots
%% =========================================================================
% Plots longitudinal Local Field Potential (LFP) amplitude trends for one
% or more patients across different sensing configurations
% (e.g., frequencies and contact pairs), alongside corresponding
% stimulation parameters (amplitude, stimulation rate, and pulse width).

% Inputs
% DBS_data (struct):
% patient_name (char or cell array of chars):
% Patient name(s) corresponding to fields in DBS_data.
% hemi (char, optional): 'Left', 'Right', or 'both' (default is 'both').

patient_name = 'Patient_RN';
hemi = 'Right';
LFP_All_Plots(DBS_data, patient_name, hemi)


%% =========================================================================
% Function 5: LFP_Daily_Plot
%% =========================================================================
%   This function plots daily LFP or TEED values over time for a specific
%   patient (or group of patients) and hemisphere. It also overlays event
%   markers and sensing frequency bands to contextualize trends in neural
%   data across the postoperative timeline.
%
% Inputs:
%   - DBS_data: structured input containing DBS LFP/TEED data for patients.
%   - patient_name: string or cell array of strings specifying patient(s).
%   - hemi (optional): 'Left', 'Right', or 'both' (default = 'both').
%   - sensing_range (optional): [min_freq max_freq] range for LFP bands.
%   - day_range (optional): [start_day, end_day], or cell array of ranges.
%
% Name-Value Pairs:
%   - 'Plot_TEED': logical (default = false). If true, plots TEED instead of LFP.
%   - 'Show_Events': logical (default = true). If true, overlays event data.
%   - 'ShowEventTypes': logical. Show event markers by type using xlines.
%   - 'ShowEventBands': logical. Show event AUC per band using scatter.
%
% Usage Example:
patient_name = 'Patient_RN';
hemi = 'Left';
sensing_range = [13 30]; % e.g., beta band
day_range = [1 60];
LFP_Daily_Plot(DBS_data, patient_name, hemi, sensing_range, day_range, ...
    'Show_Events', true, 'Plot_TEED', false);

% Notes:
%   - The function uses OpenDataHelper to extract trimmed and averaged LFP or TEED data.
%   - If 'Show_Events' is true:
%       - With 'ShowEventTypes': Events appear as labeled vertical lines.
%       - With 'ShowEventBands': Events appear as dots showing AUC in different bands.
%   - A background color patch is used to indicate changes in sensing frequency or contact pairs.
%   - If only one sensing configuration is found, it is shown in the figure title.


%% =========================================================================
% Function 6: LFP_Daily_Plot
%% =========================================================================
%   This function plots the Power Spectral Density (PSD) periodogram of LFP
%   data for a given patient and hemisphere across a specified postoperative
%   day range. It is useful for identifying dominant neural oscillation
%   cycles (e.g., circadian or ultradian rhythms) in long-term recordings.
%
% Inputs:
%   - DBS_data: structured input containing DBS LFP data for patients.
%   - patient_name: string specifying the patient ID.
%   - hemi: 'Left' or 'Right' hemisphere to analyze.
%   - periodogram_days (optional): [start_day, end_day] defining the range
%                                  of days post-operation to include in analysis.

% Usage Example:
patient_name = 'Patient_RN';
hemi = 'Left';
sensing_range = [13 30]; % e.g., beta band
day_range = [1 60];
Periodogram_Plot(DBS_data, patient_name, hemi, day_range);


% Notes:
%   - It computes the PSD using Welch's method with overlapping windows.
%   - X-axis shows oscillation cycles in hours; Y-axis shows power.
%   - Contact pair and sensing frequency are auto-labeled using
%     `Parse_Sensing_Contact_Info`.


%% =========================================================================
% Function 7: Plot_Event_PSD_Peaks
%% =========================================================================
%   This function generates a bubble chart visualizing event-related PSD peak
%   characteristics for specified patients and hemispheres over a given day range.
%   Each bubble's position corresponds to days since surgery (X-axis) and peak
%   frequency (Y-axis), with bubble size representing peak bandwidth and color
%   representing peak amplitude.
%
% Inputs:
%   - DBS_data: nested structured data containing DBS LFP and event information.
%   - patient_name: string or cell array of strings specifying patient(s) to analyze.
%   - hemi (optional): 'Left', 'Right', or 'both' (default: 'both').
%   - sensing_range (optional): [Min_frequency Max_frequency] array defining frequency band of interest.
%   - day_range (optional): [start_day end_day] range or cell array of ranges for each patient.
%
% Usage Example:
patient_name = 'Patient_RN';
hemi = 'Left';
sensing_range = [13 30]; % beta band
day_range = [1 60];
Plot_Event_PSD_Peaks(DBS_data, patient_name, hemi, sensing_range, day_range);

% Notes:
%   - Uses PSD_Fit to fit PSD and extract peak parameters.
%   - Colors represent peak amplitude (normalized).
%   - Bubble sizes represent peak bandwidth (scaled).
%   - Handles multiple patients and hemispheres.
%   - Plots frequency peaks as a function of days since surgery.


%% =========================================================================