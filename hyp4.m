%% ERP
clear all; clc; tic
%%
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if isOctave
    addpath('/media/cygnuseco/ext4_files/research/src/')
    addpath('/media/cygnuseco/ext4_files/research/eeglab-2022.0')
    addpath('/media/cygnuseco/ext4_files/research/bva-io')
    pkg load signal
end

% eeglab;

%% Choose experiment conditions
hypothesis_data.event_info = struct();
hypothesis_data.experiment_name = 'Hypothesis_4_highpass_0p5_v2'
hypothesis_data.event_types = {'manmade_old_hit_subs_rem';% subsequently remembered
               'manmade_old_hit_subs_for';% subsequently forgotten
               'manmade_old_miss_subs_rem';
               'manmade_old_miss_subs_for';
               'natural_old_hit_subs_rem';
               'natural_old_hit_subs_for';
               'natural_old_miss_subs_rem';
               'natural_old_miss_subs_for'}

hypothesis_data.condition_values = [1111;
                    1110;
                    1121
                    1120
                    2111
                    2110
                    2121
                    2110];
hypothesis_data.num_trials = 1200;
hypothesis_data.sample_rate_hz = 512;
hypothesis_data.low_pass_upper_limit_hz = 30;
hypothesis_data.high_pass_lower_limit_hz = 0.5;
hypothesis_data.pre_stimulus_ms = 200;
hypothesis_data.post_stimulus_ms = 800;
hypothesis_data.baseline_correction_time_ms = 200;








%% Folder and file paths
use_reref = true;
if use_reref; reref = 'reref_'; else; reref = ''; end;

if isOctave
    folder = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision';
    folder_subject_match = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/*.vhdr';
    folder_generated_data = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_generated_matv7';
    folder_analysed_data = ['/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_analysed_' reref hypothesis_data.experiment_name]; 
else
    folder = 'w:\EMP_data\EMP_data\eeg_brainvision\';
    folder_subject_match = 'w:\EMP_data\EMP_data\eeg_brainvision\*.vhdr';
    folder_generated_data = 'w:\EMP_data\EMP_data\eeg_brainvision\_generated_matv7';
    folder_analysed_data = ['w:\EMP_data\EMP_data\eeg_brainvision\_analysed_' reref hypothesis_data.experiment_name];
end

if (exist(folder_analysed_data, 'file') == 0)
    mkdir(folder_analysed_data); 
end

folder_subject_root = fileparts(folder_subject_match); % not used/necessary in linux (!)
subject_files = ls(folder_subject_match);
subject_total = size(subject_files, 1);


%% ---------------------------------------------
subjects_to_use = 1:subject_total;
electrodes_to_use = 1:64% fronto_central_channels;
%% ---------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Processing I - Extract all ERPs based on events
process(hypothesis_data, folder_analysed_data, folder_subject_root, subject_files, subjects_to_use, folder_generated_data, electrodes_to_use, use_reref)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Post processing - calculate mean of all trials; per subject, electrode and event

% subjects_to_use = subjects_to_use;
% selected_electrodes = electrodes_to_use;
% subject_summary_filename = postprocess(folder_analysed_data, subjects_to_use, selected_electrodes, hypothesis_data.event_types, subject_files);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Post processing - spectral analysis, high resolution
subjects_to_use = 1:33;
selected_electrodes = 1:64% fronto_central_channels;
time_window_ms = [200 1000]; % unreferenced
summary_spectral = plot_spectral_analysis_high_res(hypothesis_data.sample_rate_hz, time_window_ms, folder_analysed_data, subjects_to_use, selected_electrodes, hypothesis_data.event_types, subject_files, hypothesis_data.pre_stimulus_ms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post processing II - FFT power spectral analysis

% selected_subjects = subjects_to_use;
% selected_electrodes = [58 59 60 62];
% time_window_ms = [200 800]; % unreferenced
% plot_spectral_analysis(subject_summary_filename, hypothesis_data.event_types, selected_subjects, selected_electrodes, time_window_ms, hypothesis_data.pre_stimulus_ms)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post processing III - Calculate mean

% selected_subjects = subjects_to_use;
% selected_electrodes = [58 59 60 62];
% plot_merged_subjects(subject_summary_filename, hypothesis_data.event_types, selected_subjects, selected_electrodes, hypothesis_data.pre_stimulus_ms);