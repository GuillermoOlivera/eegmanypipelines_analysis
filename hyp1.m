%% ERP
clear all; clc; tic;
%%
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if isOctave
    addpath('/media/cygnuseco/ext4_files/research/src/')
    addpath('/media/cygnuseco/ext4_files/research/eeglab-2022.0')
    addpath('/media/cygnuseco/ext4_files/research/bva-io')
    pkg load signal
end

% eeglab;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experiment conditions:

%%
% If we want to have accurate/testable results for the frequency analysis, its
% best to not filter the data and not average out all the ERPs. Of couse the
% disadvantage of this approach is a very big file that will be saved.
%
% Lets do this for one subject with limited amount electrodes to speed up and
% after confirming plausible results run this for all subjects.
%
% Note: Only Hypothesis I asks for a specific time-interval from where to extract
%       the theta, delta power. For the other hypothesis then time-range is free
%       to choose from which means we have to do a spectogram + topoplot to choose
%       correct time-range and electrodes.
DO_SPECTRAL_ANALYSIS = true; if DO_SPECTRAL_ANALYSIS; dospectralanalysis = '_dospectralanalysis'; else; dospectralanalysis = ''; end;

%%
% We do a another pass to re-reference all the electrado voltage signal.
%
% Warning: If this electrode is "invalid" it will introduce noise to all our data
USE_REREF = false; if USE_REREF; reref = '_reref'; else; reref = ''; end;

%%
hypothesis_data = struct();
hypothesis_data.use_reref = USE_REREF;
hypothesis_data.use_spectral_analysis = DO_SPECTRAL_ANALYSIS;
hypothesis_data.experiment_name = ['Hypothesis_1_v08' reref dospectralanalysis];
hypothesis_data.event_types = {'manmade' ;
                               'natural'};
% Use the digit '6' to ignore code position
hypothesis_data.condition_values = [1666;
                                    2666];
hypothesis_data.num_trials = 1200; % only used for memory allocation
hypothesis_data.sample_rate_hz = 512;
hypothesis_data.low_pass_upper_limit_hz = 30;
hypothesis_data.high_pass_lower_limit_hz = 1.8; % We extract 1 second segments so it makes sense to filter above ~2Hz
hypothesis_data.pre_stimulus_ms = 200;
hypothesis_data.post_stimulus_ms = 800;
hypothesis_data.baseline_correction_time_ms = 200;

%% Folder and file paths
if isOctave
    folder = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision';
    folder_subject_match = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/*.vhdr';
    folder_generated_data = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_generated_matv7';
    folder_analysed_data = ['/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_analysed_' hypothesis_data.experiment_name]; 
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
subjects_to_use = 1:subject_total;
electrodes_to_use = 1:64;% fronto_central_channels;

hypothesis_data.electrodes_to_use = electrodes_to_use;
hypothesis_data.subjects_to_use = subjects_to_use;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Processing I - Extract all ERPs based on event types
process(hypothesis_data, folder_analysed_data, folder_subject_root, subject_files, subjects_to_use, folder_generated_data, electrodes_to_use, USE_REREF, DO_SPECTRAL_ANALYSIS)

%%%%%  Processing II - Extract all ERPs based on event types
eegdata_to_struct(hypothesis_data, folder_analysed_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Post processing - calculate mean of all trials; per subject, electrode and event

% subjects_to_use = subjects_to_use;
% selected_electrodes = electrodes_to_use;
% subject_summary_filename = postprocess(folder_analysed_data, subjects_to_use, selected_electrodes, hypothesis_data.event_types, subject_files);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Post processing - spectral analysis, high resolution
% subjects_to_use = 1:33;
% selected_electrodes = 1:64% fronto_central_channels;
% time_window_ms = [500 700]; % unreferenced
% summary_spectral = plot_spectral_analysis_high_res(hypothesis_data.sample_rate_hz, time_window_ms, folder_analysed_data, subjects_to_use, selected_electrodes, hypothesis_data.event_types, subject_files, hypothesis_data.pre_stimulus_ms);

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