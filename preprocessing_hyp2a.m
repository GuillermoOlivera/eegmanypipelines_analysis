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
% folder = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\';
% folder_generated_data = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\_generated';
% folder_analysed_data = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\_analysis_hy2a_6';
% folder_subject_match = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\*.vhdr';

folder = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision';
folder_subject_match = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/*.vhdr';
folder_generated_data = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_generated_matv7';
folder_analysed_data = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_analysed_full_10'; % 5 is bad!

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
event_types = ['manmade_new' ;
              'manmade_old' ;
               'natural_new' ;
               'natural_old']

% Use the digit "6" to ignore code position
condition_values = [1066;
                    1166;
                    2066
                    2166];

%% For Hypothesis 3
% condition_values = [1066;
%                     2166];

for id = 1:size(event_types, 1)
    event_type = event_types(id, :);
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
    tic; disp(["Loading mat file..." file_name_data])
    load(file_name_data);
    toc; disp("Loaded.")

    tic; disp(""); disp("Performing analysis...")
    %% EOG correction
    % Data were already re-referenced to channel 30 (POz)
    % Consider if we want to imprement other re-reference: all-channels
    % average? Can we use symmetrical channels (eg left/right mastoid
    % given that it seems an extensive cap (72 chn!), if this is not possible,
    % it would probably make sense to go for the average too.. think about it)
    vEOG_data = apply_filters(loaded_raw_data_from_eeglab.data(71,:), sample_rate_hz, low_pass_upper_limit_hz);
    hEOG_data = apply_filters(loaded_raw_data_from_eeglab.data(72,:), sample_rate_hz, low_pass_upper_limit_hz);

    %% create ERP for each condition per electrode
    data = [];
    all_triggers_tmp = [];
    electrodes_erp = struct();
    for electrode = electrodes_to_use
        electrode_id_str = ['ch_' num2str(electrode)];
        for id = 1:size(event_types, 1)
            event_type = event_types(id, :);
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

        % Apply filter
        data = apply_filters(loaded_raw_data_from_eeglab.data(electrode, :), sample_rate_hz, low_pass_upper_limit_hz);
        all_triggers_tmp = all_triggers;

        % toc; tic;
        %% Important part starts here. Here we go through each defined event type (eg. manmade new, natural new, ...)
        for trigger_id = 1:number_of_trials
            trigger = all_triggers_tmp(trigger_id, 1);
            trigger_in_digits = fdec2base(trigger, 10) - '0';

            for id = 1:size(event_types, 1)
                event_type = event_types(id, :);
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
    save(subject_analysis_filename, 'electrodes_erp', '-hdf5') % in Octave can also use hdf5 output format '-hdf5'
    toc; disp("Analysis done.")

end % for subjects_to_use



%%%%%  Post processing I - mean across subjects
summary = struct();
subjects_to_use = 1:33;
selected_electrodes = 1:72;

subject_summary_filename = fullfile(folder_analysed_data, 'summary.hdf');
isSummarized = exist(subject_summary_filename, 'file') == 2;

% Memory allocation
for subject = subjects_to_use
    if isSummarized continue end

    subject_id_str = ['subj_' num2str(subject)];
    for electrode = selected_electrodes
        electrode_id_str = ['ch_' num2str(electrode)];
        for id = 1:size(event_types, 1)
            event_type = event_types(id, :);
            summary.(subject_id_str).(event_type).(electrode_id_str) = [];
        end
    end
end

for subject = subjects_to_use
    if isSummarized continue end
    disp(subject)
    subject_id_str = ['subj_' num2str(subject)];
    subject_root_vhdr_name = subject_files(subject, :);
    subject_root_name = erase(subject_root_vhdr_name, '.vhdr');
    if isOctave
        [dir, subject_root_name, ext] = fileparts(subject_root_name);
    end

    subject_analysis_filename = fullfile(folder_analysed_data, [subject_root_name '_analysed.hdf']);
    load(subject_analysis_filename); % variable name: electrodes_erp

    total_electrodes = 1;
    for electrode = selected_electrodes
        electrode_id_str = ['ch_' num2str(electrode)];
        channel_summary = struct();
        for id = 1:size(event_types, 1)
            event_type = event_types(id, :);
            erp_matrix = mean(electrodes_erp.(event_type).(electrode_id_str), 2);

            if isempty(summary.(subject_id_str).(event_type).(electrode_id_str))
                summary.(subject_id_str).(event_type).(electrode_id_str) = erp_matrix;
            end

            summary.(subject_id_str).(event_type).(electrode_id_str) = summary.(subject_id_str).(event_type).(electrode_id_str) + erp_matrix;
            total_electrodes = total_electrodes + 1;
        end
    end

    for electrode = selected_electrodes
        electrode_id_str = ['ch_' num2str(electrode)];
        for id = 1:size(event_types, 1)
            event_type = event_types(id, :);
            summary.(subject_id_str).(event_type).(electrode_id_str) / total_electrodes;
        end
    end

end

if ~isSummarized
    save(subject_summary_filename, 'summary', '-hdf5') % in Octave can also use hdf5 output format '-hdf5'
end

%%
load(subject_summary_filename);
for id = 1:size(event_types, 1)
    event_type = event_types(id, :);
    subject_erp_mean.(event_type) = [];
end

erp_matrix = [];
leg = {};
subjects_to_use = 1:33
total_subjects = 1;
for subject = subjects_to_use
    subject_id_str = ['subj_' num2str(subject)]
    for id = 1:size(event_types, 1)
        event_type = event_types(id, :);
        erp_cell = struct2cell(summary.(subject_id_str).(event_type));
        erp_matrix = cell2mat(erp_cell');

        if isempty(subject_erp_mean.(event_type))
            subject_erp_mean.(event_type) = erp_matrix;
        end

        subject_erp_mean.(event_type) = (subject_erp_mean.(event_type) + erp_matrix);
        total_subjects = total_subjects + 1;
    end
    % plot(mean(subject_erp_mean.(event_type)(:, [58 59]), 2), 'LineWidth', 4); hold on;
    % leg{end+1} = num2str(subject)
    % legend(leg, 'FontSize', 14);
end

for id = 1:size(event_types, 1)
    event_type = event_types(id, :);
    subject_erp_mean.(event_type) / total_subjects;
end

%% Hypothesis checking
selected_electrodes = [58 59 60 62];
hyp_1 = mean(subject_erp_mean.manmade_new(:, selected_electrodes), 2);
hyp_2 = mean(subject_erp_mean.manmade_old(:, selected_electrodes), 2);
hyp_3 = mean(subject_erp_mean.natural_new(:, selected_electrodes), 2);
hyp_4 = mean(subject_erp_mean.natural_old(:, selected_electrodes), 2);
figure; plot(1/512*(1:length(hyp_1)), hyp_1, 'LineWidth', 4)
hold on
plot(1/512*(1:length(hyp_1)), hyp_2, 'LineWidth', 4)
plot(1/512*(1:length(hyp_1)), hyp_3, 'LineWidth', 4)
plot(1/512*(1:length(hyp_1)), hyp_4, 'LineWidth', 4)
legend('manmade new', 'manmade old', 'natural new', 'natural old');
 
 
 
%%%%%  Post processing II - mean for all subjects
%%%%  all_segments_erp_summary actually not needed because we load it either way
%%% tic
%%% for id = 1:size(event_types, 1)
%%%     event_type = event_types(id, :);
%%%     add_to_summary(event_type, session_filename, post_session_filename, folder_analysed_data, selected_electrodes);
%%% end
%%% toc



%% Post processing III - do welch power spectral analysis

%%% tic
%%% disp("Start spectral analysis...")
%%% 
%%% Legend=cell(2,1)%  two positions 
%%% iter = 1;
%%% erp_segment_test = [];
%%% erp_segment_test2 = [];
%%% for subject = subjects_to_use
%%%     subject_root_vhdr_name = subject_files(subject, :);
%%%     subject_root_name = erase(subject_root_vhdr_name, '.vhdr');
%%%     if isOctave
%%%         [dir, subject_root_name, ext] = fileparts(subject_root_name);
%%%     end
%%% 
%%%     subject_analysis_filename = fullfile(folder_analysed_data, [subject_root_name '_analysed.hdf']);
%%%     load(subject_analysis_filename) % variable name: all_segments_erp
%%%     disp(subject)
%%%     for electrode = selected_electrodes
%%%         channel_spectral_power_density = struct();
%%%         for id = 1:size(event_types, 1)
%%%             event_type = event_types(id, :);
%%%             electrode_id_str = ['ch_' num2str(electrode)];
%%%             erp_segments_matrix = all_segments_erp.all_electrodes.(electrode_id_str).(event_type);
%%% 
%%%             num_trials_in_event = size(erp_segments_matrix, 1);
%%% 
%%%             %%-----------------
%%%             if id == 2 && electrode == 11
%%%                 % window = 0.6 * sample_rate_hz; % 600ms
%%%                 data = mean(erp_segments_matrix);
%%%                 [spectra, Pxx_ci, freq] = pwelch(data, [], [], [], sample_rate_hz);
%%%                 keyboard()
%%%                 % erp_segment = mean(erp_segments_matrix);
%%%                 erp_segment_test = [erp_segment_test; Pxx_ci];
%%%                 erp_segment_test2 = [erp_segment_test2; freq];
%%%             end
%%%             % low, high = 0.5, 4 % alpha limits
%%%             %%-----------------
%%% 
%%%             % channel_spectral_power_density.(event_type).alpha_power = pwelch_alpha;
%%%             % channel_spectral_power_density.(event_type).beta_power = pwelch_beta;
%%% 
%%%         end
%%% 
%%%         % all_segments_erp_summary.session.(['subject_' num2str(subject)]).(['channel_' num2str(electrode)]) = channel_summary;
%%%     end
%%% end
%%% plot(mean(erp_segment_test)', mean(erp_segment_test2)');
% Legend{iter}=strcat('Electrode', num2str(electrode));
% iter = iter + 1;
% legend(Legend);
% Define delta lower and upper limits

toc