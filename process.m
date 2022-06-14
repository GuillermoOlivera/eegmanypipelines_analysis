function process(hypothesis_data, folder_analysed_data, folder_subject_root, subject_files, subjects_to_use, folder_generated_data, electrodes_to_use, use_reref, skip_filtering)
% 1. Convert to mat file from brainvision file format with help of eeglab/plugin
% 2. Load mat file and save ERPs for each condition in matrix

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

event_types = hypothesis_data.event_types;
condition_values = hypothesis_data.condition_values;
num_trials = hypothesis_data.num_trials;
sample_rate_hz = hypothesis_data.sample_rate_hz;
low_pass_upper_limit_hz = hypothesis_data.low_pass_upper_limit_hz;
high_pass_lower_limit_hz = hypothesis_data.high_pass_lower_limit_hz;
pre_stimulus_ms = hypothesis_data.pre_stimulus_ms;
post_stimulus_ms = hypothesis_data.post_stimulus_ms;
baseline_correction_time_ms = hypothesis_data.baseline_correction_time_ms;

for id = 1:size(event_types, 1)
    event_type = event_types{id};
    event_info.(event_type).condition_to_use = condition_values(id);
    event_info.(event_type).condition_to_use_in_digits = fdec2base(condition_values(id), 10) - '0';
end

baseline_correction_samples = floor(baseline_correction_time_ms / 1000 * 512);
pre_stimulus = floor(pre_stimulus_ms / 1000 * sample_rate_hz);
post_stimulus = floor(post_stimulus_ms / 1000 * sample_rate_hz);
length_segment = pre_stimulus + post_stimulus + 1;

%hypothesis_data = struct();
%hypothesis_data.num_trials = num_trials;
%hypothesis_data.pre_stimulus_ms = pre_stimulus_ms;
%hypothesis_data.post_stimulus_ms = post_stimulus_ms;
%hypothesis_data.length_segment_sanity_check = length_segment;
%hypothesis_data.sample_rate_hz = sample_rate_hz;
%hypothesis_data.low_pass_upper_limit_hz = low_pass_upper_limit_hz;

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
    vEOG_data = apply_filters(loaded_raw_data_from_eeglab.data(71,:), sample_rate_hz, low_pass_upper_limit_hz, high_pass_lower_limit_hz);
    hEOG_data = apply_filters(loaded_raw_data_from_eeglab.data(72,:), sample_rate_hz, low_pass_upper_limit_hz, high_pass_lower_limit_hz);

    %% create ERP for each condition per electrode
    data = [];
    all_triggers_tmp = [];
    electrodes_erp = struct();
    for electrode= electrodes_to_use
        % electrode = electrodes_to_use(electrode_idx);
        electrode_id_str = ['ch_' num2str(electrode)];
        for id = 1:size(event_types, 1)
            event_type = event_types{id};
            trial_str = [event_type '_trials_ok'];
            error_str = [event_type '_error_code'];
            hypothesis_data.(trial_str) = zeros(length(electrodes_to_use), 1)';
            hypothesis_data.(error_str) = zeros(length(electrodes_to_use), 1)';
        end
    end
    electrodes_erp = zeros(size(event_types,1), length(electrodes_to_use), number_of_trials, ceil(sample_rate_hz*((pre_stimulus_ms + post_stimulus_ms)/1000)));
    counter = zeros(size(event_types, 1), length(electrodes_to_use));
    counter_electrode = 0;
    for electrode = electrodes_to_use
        fprintf('%d ', electrode);
        electrode_id_str = ['ch_' num2str(electrode)];

        counter_electrode = counter_electrode + 1;

        % ANG: re-ref
        if use_reref
            refM1 = loaded_raw_data_from_eeglab.data(69,:);
            refM2 = loaded_raw_data_from_eeglab.data(70,:);
            input_data = (loaded_raw_data_from_eeglab.data(electrode,:) - (refM1 + refM2 / 2));
        else
            input_data = loaded_raw_data_from_eeglab.data(electrode,:);
        end

        % Apply filter or not
        if skip_filtering
            data = input_data;
        else
            data = apply_filters(input_data, sample_rate_hz, low_pass_upper_limit_hz, high_pass_lower_limit_hz);
        end
        
        all_triggers_tmp = all_triggers;

        %% Important part starts here. Here we go through each defined event type (eg. manmade new, natural new, ...)
        %counter_trial_ok = 0;
        for trigger_id = 1:number_of_trials
            trigger = all_triggers_tmp(trigger_id, 1);
            trigger_in_digits = fdec2base(trigger, 10) - '0';

            for id = 1:size(event_types, 1)
                event_type = event_types{id};
                counter_ok_segments = 1;
                trial_str = [event_type '_trials_ok'];
                error_str = [event_type '_error_code'];
                condition_in_digits = event_info.(event_type).condition_to_use_in_digits;

                if trigger_in_digits(condition_in_digits~=6) == condition_in_digits(condition_in_digits~=6)
                    stimulus_datapoint = all_triggers_tmp(trigger_id, 2); % vmrk file gives datapoint back! NOT timepoint
                    selected_datapoints = stimulus_datapoint - pre_stimulus:stimulus_datapoint + post_stimulus;
                    channel_segment = data(selected_datapoints); 

                    %% check there is no eye blink or artifact on any of the channels
                    % peak-to-peak voltage vs threshold voltage
                    % peak-to-peak measure less distorted by slow changes in baseline voltage, reducing the impact of source of
                    % misses and false alarms, and increasing the sensitivity of the artifact rejection process (see Luck'book)
                    segment_ok = 1;
                    ch_not_ok = max(channel_segment) - min(channel_segment) > 200;
                    veog_not_ok = (max(vEOG_data(selected_datapoints)) - min(vEOG_data(selected_datapoints)) > 160);
                    heog_not_ok = (max(hEOG_data(selected_datapoints)) - min(hEOG_data(selected_datapoints)) > 160); % paul, 160(80),120(60),60(30)/luck,200(100),160(80),160(80)
                    if ch_not_ok || veog_not_ok || heog_not_ok
                        hypothesis_data.(error_str)(electrode) = 1*ch_not_ok + 10*veog_not_ok + 100*heog_not_ok;
                        segment_ok = 0;
                    end

                    % only takes good segments
                    if segment_ok == 1
                        %counter_trial_ok = counter_trial_ok + 1;
                        counter(id, counter_electrode) = counter(id, counter_electrode) + 1;
                        hypothesis_data.(trial_str)(counter_electrode) = hypothesis_data.(trial_str)(counter_electrode) + 1; % count number of segments
                        bias = mean(channel_segment(1:baseline_correction_samples));
                        electrodes_erp(id,  counter_electrode, counter(id, counter_electrode), :) = channel_segment - bias;
                    end
                end % end if is_match
            end % for trigger_id
        end % for id = 1:size(event_types, 1)
    end % for electrodes_to_use

    electrodes_erp_info = hypothesis_data;
    electrodes_erp_counter = counter;
    save(subject_analysis_filename, 'electrodes_erp*', '-mat7-binary') % in Octave can also use hdf5 output format '-hdf5'
    toc; disp('Analysis done.')

end % for subjects_to_use

end