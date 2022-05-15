%% ERPs N1
clear all
clc
tic

%%
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if isOctave
    addpath('/media/cygnuseco/ext4_files/research/src/')
    addpath('/media/cygnuseco/ext4_files/research/eeglab-2022.0')
    addpath('/media/cygnuseco/ext4_files/research/bva-io')
    pkg load signal
end

% eeglab;

%% Hardcoded folder and file paths
use_reref = false;
if use_reref; reref = 'reref_'; else; reref = ''; end;

if isOctave
    folder = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision';
    folder_subject_match = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/*.vhdr';
    folder_generated_data = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_generated_matv7';
    folder_analysed_data = ['/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_analysed_' reref 'hyp3']; 
else
    folder = 'w:\EMP_data\EMP_data\eeg_brainvision\';
    folder_subject_match = 'w:\EMP_data\EMP_data\eeg_brainvision\*.vhdr';
    folder_generated_data = 'w:\EMP_data\EMP_data\eeg_brainvision\_generated_matv7';
    folder_analysed_data = 'w:\EMP_data\EMP_data\eeg_brainvision\_analysed_hyp3';
end

folder_subject_root = fileparts(folder_subject_match); % not used/necessary in linux (!)
session_filename = 'all_session';
post_session_filename = 'all_post_session';
subject_files = ls(folder_subject_match);
subject_total = size(subject_files, 1);

if (exist(folder_analysed_data, 'file') == 0)
    mkdir(folder_analysed_data); 
end

%% Choose subjects and electrodes (channels)
subjects_to_use = 1:subject_total;
electrodes_to_use = 1:72;
    
event_info = struct();
event_types = {'manmade_old_hit' ;
               'manmade_old_miss' ;
               'natural_old_hit' ;
               'natural_old_miss'}

% Use the digit '6' to ignore code position
condition_values = [1116;
                    1126
                    2116
                    2126];

for id = 1:size(event_types, 1)
    event_type = event_types{id};
    event_info.(event_type).condition_to_use = condition_values(id);
    event_info.(event_type).condition_to_use_in_digits = fdec2base(condition_values(id), 10) - '0';
end

%% Channel position for hypothesis 2 - NOT USED YET
% Frontal central channels
%     9->FC5
%     10->FC3
%     11->FC1
%     44->FC6
%     45->FC4
%     46->FC2
%     47->FCz
fronto_central_channels = [9, 10, 11, 44, 45, 46, 47];

%%
num_trials = 1200;
sample_rate_hz = 512;
low_pass_upper_limit_hz = 30;
pre_stimulus_ms = 200;
post_stimulus_ms = 800;
baseline_correction_time_ms = 250;

baseline_correction_samples = floor(baseline_correction_time_ms / 1000 * 512);
pre_stimulus = floor(pre_stimulus_ms / 1000 * sample_rate_hz);
post_stimulus = floor(post_stimulus_ms / 1000 * sample_rate_hz);
length_segment = pre_stimulus + post_stimulus + 1;

info = struct();
info.num_trials = num_trials;
info.pre_stimulus_ms = pre_stimulus_ms;
info.post_stimulus_ms = post_stimulus_ms;
info.length_segment_sanity_check = length_segment;
info.sample_rate_hz = sample_rate_hz;
info.low_pass_upper_limit_hz = low_pass_upper_limit_hz;

%% Memory pre-allocation
all_segments_erp = struct();

for subject = subjects_to_use
    subject_root_vhdr_name = subject_files(subject, :);

    subject_root_name = erase(subject_root_vhdr_name, '.vhdr');

    if isOctave
       [dir, name, ext] = fileparts(subject_root_vhdr_name);
       subject_root_vhdr_name = [name ext];
       [dir, subject_root_name, ext] = fileparts(subject_root_name);
    end

    %% Check if file already exists and skip analysis if true
    subject_analysis_filename = fullfile(folder_analysed_data, [subject_root_name '_analysed.hdf']);
    isAnalysed = exist(subject_analysis_filename , 'file') == 2;

    if isAnalysed
        continue;
    end
    
    %% Save EEG data as .mat file
    subject_vhdr_filepath = fullfile(folder_subject_root, subject_root_vhdr_name);
    vmrk_to_mat(subject_vhdr_filepath, folder_generated_data);
    %% Get trigger times from .vmrk file
    subject_root_vmrk_name = [subject_root_name '.vmrk'];  % the name of the .vmrk file
    subject_vmrk_filepath = fullfile(folder_subject_root, subject_root_vmrk_name);
    all_triggers = zeros(num_trials,2);
    all_triggers = read_triggers_from_vmrk(subject_vmrk_filepath);
    number_of_trials = size(all_triggers, 1);

    %% Load previously saved EEG data
    file_name_data = fullfile(folder_generated_data, [subject_root_name '.out.mat']);
    tic; disp(['Loading mat file...' file_name_data])
    load(file_name_data);
    toc; disp('Loaded.')

    tic; disp(''); disp('Performing analysis...')
    
    %% EOG correction
    vEOG_data = apply_filters(loaded_raw_data_from_eeglab.data(71,:), sample_rate_hz, low_pass_upper_limit_hz);
    hEOG_data = apply_filters(loaded_raw_data_from_eeglab.data(72,:), sample_rate_hz, low_pass_upper_limit_hz);

    %% create ERP for each condition per electrode
    data = [];
    all_triggers_tmp = [];
    electrodes_erp = struct();
    for electrode = electrodes_to_use
        electrode_id_str = ['ch_' num2str(electrode)];
        for id = 1:size(event_types, 1)
            event_type = event_types{id};
            trial_str = [event_type '_trials_ok'];
            error_str = [event_type '_error_code'];
            electrodes_erp.(event_type).(electrode_id_str) = [];
            info.(trial_str) = zeros(max(electrodes_to_use), 1)';
            info.(error_str) = zeros(max(electrodes_to_use), 1)';
        end
    end

    for electrode = electrodes_to_use
        electrode
        electrode_id_str = ['ch_' num2str(electrode)];
        
        % ANG: re-ref
        if use_reref
            refM1 = loaded_raw_data_from_eeglab.data(69,:);
            refM2 = loaded_raw_data_from_eeglab.data(70,:);
            input_data = (loaded_raw_data_from_eeglab.data(electrode,:) - (refM1 + refM2 / 2));
        else
            input_data = loaded_raw_data_from_eeglab.data(electrode,:); % dont re-ref, just copy original
        end

        % Apply filter
        data = apply_filters(input_data, sample_rate_hz, low_pass_upper_limit_hz);
        all_triggers_tmp = all_triggers;

        % toc; tic;
        %% Important part starts here. Here we go through each defined event type (eg. manmade new, natural new, ...)
        for trigger_id = 1:number_of_trials
            trigger = all_triggers_tmp(trigger_id, 1);
            trigger_in_digits = fdec2base(trigger, 10) - '0';

            for id = 1:size(event_types, 1)
                event_type = event_types{id};
                trial_str = [event_type '_trials_ok'];
                error_str = [event_type '_error_code'];
                condition_in_digits = event_info.(event_type).condition_to_use_in_digits;

                if trigger_in_digits(condition_in_digits~=6) == condition_in_digits(condition_in_digits~=6)
                    stimulus_datapoint = all_triggers_tmp(trigger_id, 2); % vmrk file gives datapoint back! NOT timepoint
                    selected_datapoints = stimulus_datapoint - pre_stimulus:stimulus_datapoint + post_stimulus;
                    channel_segment = data(selected_datapoints); 

                    % disp(['condition=' num2str(condition_to_use) ' trigger_in_digits=' num2str(trigger_in_digits) ' trigger_time=' num2str(stimulus_datapoint)]);

                    %% check there is no eye blink or artifact on any of the channels
                    % peak-to-peak voltage vs threshold voltage
                    % peak-to-peak measure less distorted by slow changes in baseline voltage, reducing the impact of source of
                    % misses and false alarms, and increasing the sensitivity of the artifact rejection process (see Luck'book)
                    segment_ok = 1;
                    ch_not_ok = max(channel_segment) - min(channel_segment) > 200;
                    veog_not_ok = (max(vEOG_data(selected_datapoints)) - min(vEOG_data(selected_datapoints)) > 160);
                    heog_not_ok = (max(hEOG_data(selected_datapoints)) - min(hEOG_data(selected_datapoints)) > 160); % paul, 160(80),120(60),60(30)/luck,200(100),160(80),160(80)
                    if ch_not_ok || veog_not_ok || heog_not_ok
                        info.(error_str)(electrode) = 1*ch_not_ok + 10*veog_not_ok + 100*heog_not_ok;
                        segment_ok = 0;
                    end

                    % only takes good segments
                    if segment_ok == 1
                        info.(trial_str)(electrode) = info.(trial_str)(electrode) + 1; % count number of segments
                        bias = mean(channel_segment(1:baseline_correction_samples));
                        electrodes_erp.(event_type).(electrode_id_str) = [electrodes_erp.(event_type).(electrode_id_str), (channel_segment - bias)'];
                    end

                end % end if is_match
                
            end % for trigger_id
            
        end % for id = 1:size(event_types, 1)

    end % for electrodes_to_use

    electrodes_erp.info = info;
    save(subject_analysis_filename, 'electrodes_erp', '-mat7-binary') % in Octave can also use hdf5 output format '-hdf5'
    toc; disp('Analysis done.')

end % for subjects_to_use

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Post processing I - calculate mean of all trials per subject, electrode and event
subjects_to_use = 1:33;
selected_electrodes = 1:72;
subject_summary_filename = postprocess_step1(folder_analysed_data, subjects_to_use, selected_electrodes, event_types, subject_files);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse - FFT power spectral analysis
selected_subjects = 1:5;
selected_electrodes = [58 59 60 62];
time_window_ms = [0 801];
plot_spectral_analysis(subject_summary_filename, event_types, selected_subjects, selected_electrodes, time_window_ms, pre_stimulus_ms)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyse - Mean across subjects for all events (and selected electrodes)
selected_subjects = 1:33;
selected_electrodes = [58 59 60 62];
plot_merged_subjects(subject_summary_filename, event_types, selected_subjects, selected_electrodes, pre_stimulus_ms);