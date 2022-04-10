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

eeglab;

%% Hardcoded folder and file paths
% folder = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\';
% folder_generated_data = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\_generated';
% folder_analysed_data = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\_analysis_hy2a_6';
% folder_subject_match = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\*.vhdr';

folder = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision';
folder_subject_match = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/*.vhdr';
folder_generated_data = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_generated_matv7';
folder_analysed_data = '/media/cygnuseco/ext4_files/research/EMP_data/EMP_data/eeg_brainvision/_analysed_3';

folder_subject_root = fileparts(folder_subject_match); % not used/necessary in linux (!)
session_filename = 'all_session.mat';
post_session_filename = 'all_post_session';
subject_files = ls(folder_subject_match);
subject_total = size(subject_files, 1);

if (exist(folder_analysed_data, 'file') == 0)
    mkdir(folder_analysed_data); 
end

%% Choose subjects and electrodes (channels)
subjects_to_use = 1:subject_total;
electrodes_to_use = 1:72;
    
% % num_subjects = length(subjects_to_use);
% % num_electrodes = 72;
% % num_conditions = 4; % new vs old in man-made/natural-enviro
number_of_trials = 1200;
sample_rate_hz = 512;
low_pass_upper_limit_hz = 30;
pre_stimulus_ms = 100;
post_stimulus_ms = 600;
pre_stimulus = floor(pre_stimulus_ms / 1000 * sample_rate_hz);
post_stimulus = floor(post_stimulus_ms / 1000 * sample_rate_hz);
length_segment = pre_stimulus + post_stimulus + 1;

all_segments_erp_summary = struct();
info = struct();
info.number_of_trials = number_of_trials;
info.pre_stimulus_ms = pre_stimulus_ms;
info.post_stimulus_ms = post_stimulus_ms;
info.length_segment_sanity_check = length_segment;
info.sample_rate_hz = sample_rate_hz;
info.low_pass_upper_limit_hz = low_pass_upper_limit_hz;

%% Memory pre-allocation
all_segments_erp = struct();
% trial_data_manmade_new = [];
% trial_data_manmade_old = [];
% trial_data_natural_new = [];
% trial_data_natural_old = [];

for subject = subjects_to_use
    subject_root_vhdr_name = subject_files(subject, :);
    subject_root_name = erase(subject_root_vhdr_name, '.vhdr');
    if isOctave
       [dir, name, ext] = fileparts(subject_root_vhdr_name);
       subject_root_vhdr_name = [name ext];
       [dir, subject_root_name, ext] = fileparts(subject_root_name);
    end
    
    %% Save EEG data as .mat file
    subject_vhdr_filepath = fullfile(folder_subject_root, subject_root_vhdr_name);
    vmrk_to_mat(subject_vhdr_filepath, folder_generated_data);
    %% Get trigger times from .vmrk file
    subject_root_vmrk_name = [subject_root_name '.vmrk'];  % the name of the .vmrk file
    subject_vmrk_filepath = fullfile(folder_subject_root, subject_root_vmrk_name);
    all_triggers = zeros(number_of_trials,2);
    all_triggers = read_triggers_from_vmrk(subject_vmrk_filepath);
    %% Load previously saved EEG data
    file_name_data = fullfile(folder_generated_data, [subject_root_name '.out.mat']);
    tic
    disp(["Loading mat file..." file_name_data])
    load(file_name_data);
    toc
    disp("Loaded.")

    %% EOG correction
    % Data were already re-referenced to channel 30 (POz)
    % Consider if we want to imprement other re-reference: all-channels
    % average? Can we use symmetrical channels (eg left/right mastoid
    % given that it seems an extensive cap (72 chn!), if this is not possible,
    % it would probably make sense to go for the average too.. think about it)
    vEOG_data = apply_filters(loaded_raw_data_from_eeglab.data(71,:), sample_rate_hz, low_pass_upper_limit_hz);
    hEOG_data = apply_filters(loaded_raw_data_from_eeglab.data(72,:), sample_rate_hz, low_pass_upper_limit_hz);

    %% create ERP for each condition per electrode
    for electrode = electrodes_to_use
        % Apply filter
        data = apply_filters(loaded_raw_data_from_eeglab.data(electrode,:), sample_rate_hz, low_pass_upper_limit_hz);

        BB = struct();
        BB.('manmade_new') = zeros(number_of_trials, length_segment);
        BB.('manmade_old') = zeros(number_of_trials, length_segment);
        BB.('natural_new') = zeros(number_of_trials, length_segment);
        BB.('natural_old') = zeros(number_of_trials, length_segment);

        num_events_ok = 0; % count variable to record number of segments
        
        for condition = (1:4)
            % new (0 in 2 location) VS old (1 in 2 location)
            if condition == 1
                conditions_to_use = [1030 1031 1039 1040 1041 1049]; % new man-made
            elseif condition == 2
                conditions_to_use = [2030 2031 2039 2040 2041 2049 2091]; % new natutal env.
            elseif condition == 3
                conditions_to_use = [1110 1111 1119 1120 1121 1129]; % old man-made
            elseif condition == 4
                conditions_to_use = [2110 2111 2119 2120 2121 2129]; % old man-made
            end
            
            for trial = (1:number_of_trials) % loop through all trials
                
                if (any(all_triggers(trial,1) == conditions_to_use))
                    stimulus_timepoint = all_triggers(trial,2);
                    selected_datapoints = stimulus_timepoint - pre_stimulus:stimulus_timepoint + post_stimulus;
                    channel_segment = data(selected_datapoints); % original_ not electrode but channel!!
                    vEOG_segment = vEOG_data(selected_datapoints); % get EOG segment
                    hEOG_segment = hEOG_data(selected_datapoints);
                    
                    %% check there is no eye blink or artifact on any of the channels
                    % peak-to-peak voltage vs threshold voltage
                    % peak-to-peak measure less distorted by slow changes in baseline voltage, reducing the impact of source of
                    % misses and false alarms, and increasing the sensitivity of the artifact rejection process (see Luck'book)
                    segment_ok = 1;
                    
                    if (max(channel_segment(1: length_segment)) - min(channel_segment(1: length_segment))) > 200 ||...
                            ((max(vEOG_segment(1: length_segment)) - min(vEOG_segment(1: length_segment))) > 160) ||...
                            ((max(hEOG_segment(1: length_segment)) - min(hEOG_segment(1: length_segment))) > 160) % Paul, 160(80),120(60),60(30)/Luck,200(100),160(80),160(80)
                        
                        segment_ok = 0; % not OK
                    end
                    
                    % only take the good segments
                    if  segment_ok == 1

                        num_events_ok = num_events_ok + 1; % count number of segments
                        
                        if condition == 1
                            BB.('manmade_new')(num_events_ok, 1:length_segment) = channel_segment;
                        end
                        
                        if condition == 2
                            BB.('natural_new')(num_events_ok, 1:length_segment) = channel_segment;
                        end
                        
                        if condition == 3
                            BB.('manmade_old')(num_events_ok, 1:length_segment) = channel_segment;
                        end
                        
                        if condition == 4
                            BB.('natural_old')(num_events_ok, 1:length_segment) = channel_segment;
                        end
                    end 
                end 
            end            
        end
        
        info.num_events_ok = num_events_ok;

        BB.('manmade_new') = BB.('manmade_new')(~all(BB.('manmade_new') == 0, 2),:);
        BB.('manmade_old') = BB.('manmade_old')(~all(BB.('manmade_old') == 0, 2),:);
        BB.('natural_new') = BB.('natural_new')(~all(BB.('natural_new') == 0, 2),:);
        BB.('natural_old') = BB.('natural_old')(~all(BB.('natural_old') == 0, 2),:);

        % Save result of each electrode into cell
        % trial_data_manmade_new[electrode,:] = BB.('manmade_new');
        % trial_data_manmade_old[electrode,:] = BB.('manmade_old');
        % trial_data_natural_new[electrode,:] = BB.('natural_new');
        % trial_data_natural_old[electrode,:] = BB.('natural_old');

        erp_manmade_new_mean = mean(BB.('manmade_new'));
        erp_manmade_old_mean = mean(BB.('manmade_old'));
        erp_natural_new_mean = mean(BB.('natural_new'));
        erp_natural_old_mean = mean(BB.('natural_old'));
        
        channel_summary = struct();
        channel_summary.info = info;
        try
            channel_summary.('manmade_new') = erp_manmade_new_mean - mean(erp_manmade_new_mean(1 : pre_stimulus));
            channel_summary.('manmade_old') = erp_manmade_old_mean - mean(erp_manmade_old_mean(1 : pre_stimulus));
            channel_summary.('natural_new') = erp_natural_new_mean - mean(erp_natural_new_mean(1 : pre_stimulus));
            channel_summary.('natural_old') = erp_natural_old_mean - mean(erp_natural_old_mean(1 : pre_stimulus));
        catch err
            % warning(err.identifier, err.message);
            disp(["Error at subject (" num2str(subject) ") and  electrode ( " num2str(electrode) ")"])
        end
        
        all_segments_erp_summary.session.(['subject_' num2str(subject)]).(['channel_' num2str(electrode)]) = channel_summary;
    end

    % Save individual results
    all_segments_erp_manmade_new.info = info;
    all_segments_erp_natural_new.info = info;
    all_segments_erp_manmade_old.info = info;
    all_segments_erp_natural_new.info = info;
    all_segments_erp_manmade_new.all_electrodes = BB.('manmade_new');
    all_segments_erp_manmade_old.all_electrodes = BB.('manmade_old');
    all_segments_erp_natural_new.all_electrodes = BB.('natural_new');
    all_segments_erp_natural_old.all_electrodes = BB.('natural_old');
    save(fullfile(folder_analysed_data, [subject_root_name '_manmade_new']),'all_segments_erp_manmade_new', '-hdf5')
    save(fullfile(folder_analysed_data, [subject_root_name '_manmade_old']),'all_segments_erp_manmade_old', '-hdf5')
    save(fullfile(folder_analysed_data, [subject_root_name '_natural_new']),'all_segments_erp_natural_new', '-hdf5')
    save(fullfile(folder_analysed_data, [subject_root_name '_natural_old']),'all_segments_erp_natural_old', '-hdf5')

end

save(fullfile(folder_analysed_data, session_filename), 'all_segments_erp_summary', '-hdf5')
toc

%% Post processing
% all_segments_erp_summary actually not needed because we load it either
% way
tic
add_to_summary('manmade_new', session_filename, post_session_filename, folder_analysed_data, electrodes_to_use);
add_to_summary('manmade_old', session_filename, post_session_filename, folder_analysed_data, electrodes_to_use);
add_to_summary('natural_new', session_filename, post_session_filename, folder_analysed_data, electrodes_to_use);
add_to_summary('natural_old', session_filename, post_session_filename, folder_analysed_data, electrodes_to_use);
toc
