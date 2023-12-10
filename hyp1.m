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
% best to not filter the data and not average out all the ERPs.
%
% Lets do this for one subject with limited amount electrodes to speed up and
% after confirming plausible results run this for all subjects.
%
% Note: Only Hypothesis I asks for a specific time-interval from where to extract
%       the theta, delta power. For the other hypothesis then time-range is free
%       to choose which means we have to do a spectogram + topoplot to choose
%       correct time-range and electrodes.
DO_SPECTRAL_ANALYSIS = false; if DO_SPECTRAL_ANALYSIS; dospectralanalysis = '_dospectralanalysis'; else; dospectralanalysis = ''; end;

%%
% We do a another pass to re-reference all the electrado voltage signal.
%
% Warning: If this electrode is "invalid" it will introduce noise to all our data
USE_REREF = true; if USE_REREF; reref = '_reref'; else; reref = ''; end;

%%
info = struct();
info.use_reref = USE_REREF;
info.use_spectral_analysis = DO_SPECTRAL_ANALYSIS;
info.experiment_name = ['Hypothesis_1' reref dospectralanalysis];
info.event_types = {'manmade' ;
                    'natural'};
% Use the digit '6' to ignore code position
info.condition_values = [1666;
                        2666];
info.num_trials = 1200; % only used for memory allocation
info.sample_rate_hz = 512;
info.low_pass_upper_limit_hz = 30;
info.high_pass_lower_limit_hz = 4;
info.pre_stimulus_ms = 200;
info.post_stimulus_ms = 1100;
info.baseline_correction_time_ms = 200;

%% Folder and file paths
if isOctave
    folder = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision';
    folder_subject_match = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/*.vhdr';
    folder_generated_data = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_generated_matv7';
    folder_analysed_data = ['/media/cygnuseco/ext4_files/research/EMP_data/_final/' info.experiment_name]; 
else
    folder = 'w:\EMP_data\EMP_data\eeg_brainvision\';
    folder_subject_match = 'w:\EMP_data\EMP_data\eeg_brainvision\*.vhdr';
    folder_generated_data = 'w:\EMP_data\EMP_data\eeg_brainvision\_generated_matv7';
    folder_analysed_data = ['w:\EMP_data\EMP_data\eeg_brainvision\_analysed_' reref info.experiment_name];
end

if (exist(folder_analysed_data, 'file') == 0)
    mkdir(folder_analysed_data); 
end

folder_subject_root = fileparts(folder_subject_match); % not necessary in linux
subject_files = ls(folder_subject_match);
subject_total = size(subject_files, 1);
subjects_to_use = 1:subject_total;

electrodes_to_use = 1:64;
info.electrodes_to_use = electrodes_to_use;
info.subjects_to_use = subjects_to_use;

% Referenced to pre_stimulus_ms. For min frequency reconstruction at least
% two full cycles needed ie, 200 ms -> 100 ms -> 10 Hz
info.spectra_time_window_ms = [100 500];

info.spectra_selected_electrodes = info.electrodes_to_use;
info.spectra_subjects_to_use = info.subjects_to_use;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Processing I - Extract all ERPs based on event types
process(info, folder_analysed_data, folder_subject_root, subject_files, folder_generated_data, USE_REREF, DO_SPECTRAL_ANALYSIS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Processing II - Extract all ERPs based on event types
no_run = false;
only_log = false;  % skip saving just log
eegdata_to_struct(info, folder_analysed_data, only_log, no_run)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Post processing - spectral analysis
summary_spectral = spectral_analysis(info, folder_analysed_data, subject_files);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Post processing - voltage analysis
% plot_merged_subjects(subject_summary_filename, info.event_types, selected_subjects, selected_electrodes, info.pre_stimulus_ms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Post processing - Merge results across electrodes and subjects
% subject_summary_filename = postprocess(folder_analysed_data, subjects_to_use, selected_electrodes, info.event_types, subject_files);