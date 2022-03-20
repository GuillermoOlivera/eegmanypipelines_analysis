%% ERPs N1
clear all
clc

%% Hardcoded folder and file paths
folder = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\';
folder_generated_data = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\_generated';
folder_analysed_data = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\_analysis_hy2a_2';
folder_subject_match = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\*.vhdr';
folder_subject_root = fileparts(folder_subject_match);
subject_files = ls(folder_subject_match);
subject_total = size(subject_files, 1);

if (exist(folder_analysed_data, 'file') == 0)
    mkdir(folder_analysed_data); 
end

%% Choose subjects and electrodes (channels)
subjects_to_use = 1:3 %1:subject_total;
electrodes_to_use = [2,10,22,33,44] % 1:72;
    
% % num_subjects = length(subjects_to_use);
% % num_electrodes = 72;
% % num_conditions = 4; % new vs old in man-made/natural-enviro
sample_rate_hz = 512;
low_pass_upper_limit_hz = 30;
pre_stimulus_ms = 100;
post_stimulus_ms = 600;
pre_stimulus = floor(pre_stimulus_ms / 1000 * sample_rate_hz);
post_stimulus = floor(post_stimulus_ms / 1000 * sample_rate_hz);
length_segment = pre_stimulus + post_stimulus + 1;

info = struct();
info.pre_stimulus_ms = pre_stimulus_ms;
info.post_stimulus_ms = post_stimulus_ms;
info.length_segment_sanity_check = length_segment;
info.sample_rate_hz = sample_rate_hz;
info.low_pass_upper_limit_hz = low_pass_upper_limit_hz;

%% Memory pre-allocation
all_segments_erp = struct();

%%
for subject = subjects_to_use
    subject_root_vhdr_name = subject_files(subject, :);
    subject_root_name = erase(subject_root_vhdr_name, '.vhdr');
    
    %% Save EEG data as .mat file
    subject_vhdr_filepath = fullfile(folder_subject_root, subject_root_vhdr_name);
    vmrk_to_mat(subject_vhdr_filepath, folder_generated_data);
    
    %% Trials
    number_of_trials = 1200;
    total_trigger_count = 0;
    
    %% Get trigger times from .vmrk file
    all_triggers = zeros(number_of_trials,2);
    subject_root_vmrk_name = [subject_root_name '.vmrk'];  % the name of the .vmrk file
    subject_vmrk_filepath = fullfile(folder_subject_root, subject_root_vmrk_name);
    all_triggers = read_triggers_from_vmrk(subject_vmrk_filepath);
    
    %% Load previously saved EEG data
    file_name_data = fullfile(folder_generated_data, [subject_root_name '.out.mat']);
    load(file_name_data);

    %% Apply filter
    data = apply_filters(loaded_raw_data_from_eeglab.data, sample_rate_hz, low_pass_upper_limit_hz);
    
    % Data were already re-referenced to channel 30 (POz)
    % Consider if we want to imprement other re-reference: all-channels
    % average? Can we use symmetrical channels (eg left/right mastoid
    % given that it seems an extensive cap (72 chn!), if this is not possible,
    % it would probably make sense to go for the average too.. think about it)
    vEOG_data = data(71,:);
    hEOG_data = data(72,:);
    

    
    %% create ERP for each condition
    for electrode = electrodes_to_use
        
        ERP_matrix_manmadenew = zeros(number_of_trials, length_segment); % empty matrix for each segment, each segment is 400ms
        ERP_matrix_naturalold = zeros(number_of_trials, length_segment); % empty matrix for each segment, each segment is 400ms
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
%                     stimulus_timepoint = floor(stimulus_time_ms / 1000 * sample_rate_hz);
                    selected_datapoints = stimulus_timepoint - pre_stimulus:stimulus_timepoint + post_stimulus;
                    channel_segment = data(electrode, selected_datapoints); % original_ not electrode but channel!!
                    channel_segment = channel_segment - mean(channel_segment(1 : pre_stimulus));
                    
                    vEOG_segment = vEOG_data(selected_datapoints); % get EOG segment
                    vEOG_segment = vEOG_segment - mean(vEOG_segment(1 : pre_stimulus)); % baseline correct
                    
                    hEOG_segment = hEOG_data(selected_datapoints);
                    hEOG_segment = hEOG_segment - mean(hEOG_segment(1 : pre_stimulus));
                    
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
                            ERP_matrix_manmadenew(num_events_ok, 1:length_segment) = channel_segment; % put the good segment into matrix
                        end
                        
                        if condition == 2
                            ERP_matrix_naturalnew(num_events_ok, 1:length_segment) = channel_segment; % put the good segment into matrix
                        end
                        
                        if condition == 3
                            ERP_matrix_manmadeold(num_events_ok, 1:length_segment) = channel_segment; % put the good segment into matrix
                        end
                        
                        if condition == 4
                            ERP_matrix_naturalold(num_events_ok, 1:length_segment) = channel_segment; % put the good segment into matrix
                        end
                    end
                    
                end
                
            end
            
        end
        
        ERP_matrix_manmadenew = ERP_matrix_manmadenew(~all(ERP_matrix_manmadenew == 0, 2),:);
        ERP_matrix_naturalnew = ERP_matrix_naturalnew(~all(ERP_matrix_naturalnew == 0, 2),:);
        ERP_matrix_manmadeold = ERP_matrix_manmadeold(~all(ERP_matrix_manmadeold == 0, 2),:);
        ERP_matrix_naturalold = ERP_matrix_naturalold(~all(ERP_matrix_naturalold == 0, 2),:);
        
        % Save result of each electrode into cell
        trial_data_manmadenew{electrode,:} = ERP_matrix_manmadenew;
        trial_data_naturalnew{electrode,:} = ERP_matrix_naturalnew;
        trial_data_manmadeold{electrode,:} = ERP_matrix_manmadeold;
        trial_data_naturalold{electrode,:} = ERP_matrix_naturalold;

        erp_manmadenew_mean = mean(ERP_matrix_manmadenew);
        erp_naturalnew_mean = mean(ERP_matrix_naturalnew);
        erp_manmadeold_mean = mean(ERP_matrix_manmadeold);
        erp_naturalold_mean = mean(ERP_matrix_naturalold);
        trial_data_summary = struct();
        trial_data_summary.erp_manmadenew_mean = erp_manmadenew_mean - mean(erp_manmadenew_mean(1: pre_stimulus));
        trial_data_summary.erp_naturalnew_mean = erp_naturalnew_mean - mean(erp_naturalnew_mean(1: pre_stimulus));
        trial_data_summary.erp_manmadeold_mean = erp_manmadeold_mean - mean(erp_manmadeold_mean(1: pre_stimulus));
        trial_data_summary.erp_naturalold_mean = erp_naturalold_mean - mean(erp_naturalold_mean(1: pre_stimulus));
        %trial_data_summary.erp_manmade_mean_without_correction = erp_manmade_mean;
        %trial_data_summary.erp_natural_mean_without_correction = erp_natural_mean;
        
        all_segments_erp_summary.data.(['electrode' num2str(electrode)]) = trial_data_summary;
    end
    
    % Save results into struct
    all_segments_erp_summary.info = info;
    all_segments_erp_manmadenew.info = info;
    all_segments_erp_naturalnew.info = info;
    all_segments_erp_manmadeold.info = info;
    all_segments_erp_naturalnew.info = info;
    
    all_segments_erp_manmadenew.all_electrodes = trial_data_manmadenew;
    all_segments_erp_naturalnew.all_electrodes = trial_data_naturalnew;
    all_segments_erp_manmadeold.all_electrodes = trial_data_manmadeold;
    all_segments_erp_naturalnew.all_electrodes = trial_data_naturalold;
        
    save(fullfile(folder_analysed_data, [subject_root_name '_manmadenew']),'all_segments_erp_manmadenew')
    save(fullfile(folder_analysed_data, [subject_root_name '_naturalnew']),'all_segments_erp_naturalnew')
    save(fullfile(folder_analysed_data, [subject_root_name '_manmadeold']),'all_segments_erp_manmadeold')
    save(fullfile(folder_analysed_data, [subject_root_name '_naturalnew']),'all_segments_erp_naturalnew')
    save(fullfile(folder_analysed_data, [subject_root_name '_summary']),'all_segments_erp_summary')
    
end

match_trial_type = ls(strcat(folder_analysed_data, '\*_manmadenew*'));

clear erp;
count = 1;
mean_all_subjects = struct();
for electrode_selected = electrodes_to_use
    erp{count} = mean_erp_for_all_subjects_from_specific_electrode(folder_analysed_data, match_trial_type, electrode_selected);
    count = count + 1;
end

all_segments_erp_summary.manmade_new_mean_metadata = electrodes_to_use;
all_segments_erp_summary.manmade_new_mean_electrodes = erp;

save(fullfile(folder_analysed_data, 'test_all'),'all_segments_erp_summary')



