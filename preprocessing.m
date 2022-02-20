%% ERPs N1
clear all
clc

folder = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\';
folder_generated_data = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\_generated';
folder_subject_match = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\*.vhdr';
folder_subject_root = fileparts(folder_subject_match);
subject_files = ls(folder_subject_match);
subject_total = size(subject_files, 1);

%% Choose subjects and electrodes (channels)
subjects_to_use = 1:subject_total;
num_subjects = length(subjects_to_use);
num_electrodes = 72;
num_conditions = 2; % man-made/natural-enviro

%% Memory pre-allocation
ERPs_valid= zeros(num_electrodes, num_conditions, num_subjects);
all_ERPs = zeros(400, num_subjects, num_conditions, num_electrodes);

all_segments_erp = struct();

%%
for subject = subjects_to_use
    subject_root_vhdr_name = subject_files(subject, :);
    subject_root_name = erase(subject_root_vhdr_name, '.vhdr');
    
    %% Save EEG data as .mat file
    subject_vhdr_filepath = fullfile(folder_subject_root, subject_root_name);
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
    data = apply_filters(loaded_raw_data_from_eeglab);

    electrodes_to_use = (1:72);
    
    % Data were already re-referenced to channel 30 (POz)
    % Consider if we want to imprement other re-reference: all-channels
    % average? Can we use symmetrical channels (eg left/right mastoid
    % given that it seems an extensive cap (72 chn!), if this is not possible,
    % it would probably make sense to go for the average too.. think about it)
    vEOG_data = data(71,:);
    hEOG_data = data(72,:);
    
    %% create ERP for each condition
    for electrode = electrodes_to_use
        
        ERP_matrix_manmade = zeros(number_of_trials, 400); % empty matrix for each segment, each segment is 400ms
        ERP_matrix_natural = zeros(number_of_trials, 400); % empty matrix for each segment, each segment is 400ms
        num_events_ok = 0; % count variable to record number of segments
            
        for condition = (1:2)
            
            if condition == 1
                conditions_to_use = [1030 1031 1039 1040 1041 1049 1110 1111 1119 1120 1121 1129];
            elseif condition == 2
                conditions_to_use = [2030 2031 2039 2040 2041 2049 2091 2110 2111 2119 2120 2121 2129];
            end

            for trial = (1:number_of_trials) % loop through all trials
                
                if (any(all_triggers(trial,1) == conditions_to_use))
                    stimulus_time = all_triggers(trial,2);
                    channel_segment = data(electrode,stimulus_time-100:stimulus_time+299); % original_ not electrode but channel!!
                    channel_segment = channel_segment - mean(channel_segment(1:100));
                    
                    vEOG_segment = vEOG_data(stimulus_time-100:stimulus_time+299); % get EOG segment
                    vEOG_segment = vEOG_segment - mean(vEOG_segment(1:100)); % baseline correct
                    
                    hEOG_segment = hEOG_data(stimulus_time-100:stimulus_time+299);
                    hEOG_segment = hEOG_segment - mean(hEOG_segment(1:100));

                    %% check there is no eye blink or artifact on any of the channels
                    % peak-to-peak voltage vs threshold voltage
                    % peak-to-peak measure less distorted by slow changes in baseline voltage, reducing the impact of source of
                    % misses and false alarms, and increasing the sensitivity of the artifact rejection process (see Luck'book)
                    segment_ok = 1;
                    
                    channel_segment = data(electrode,stimulus_time-100:stimulus_time+299); % original_ not electrode but channel!!
                    if (max(channel_segment(1:400)) - min(channel_segment(1:400))) > 200 ||...
                            ((max(vEOG_segment(1:400)) - min(vEOG_segment(1:400))) > 160) ||...
                            ((max(hEOG_segment(1:400)) - min(hEOG_segment(1:400))) > 160) % Paul, 160(80),120(60),60(30)/Luck,200(100),160(80),160(80)
                        
                        segment_ok = 0; % not OK                       
                    end

                    % only take the good segments
                    if  segment_ok == 1
                        num_events_ok = num_events_ok + 1; % count number of segments
                        
                        if condition == 1
                            ERP_matrix_manmade(num_events_ok,:) = channel_segment; % put the good segment into matrix
                        end
                        
                        if condition == 2
                            ERP_matrix_natural(num_events_ok,:) = channel_segment; % put the good segment into matrix
                        end
                    end
                    
                end
                        
            end

        end

        ERP_matrix_manmade = ERP_matrix_manmade(~all(ERP_matrix_manmade == 0, 2),:);
        ERP_matrix_natural = ERP_matrix_natural(~all(ERP_matrix_natural == 0, 2),:);

        trial_data_manmade = struct();
        trial_data_manmade.subject = subject;
        trial_data_manmade.erp = ERP_matrix_manmade;
        all_segments_erp_manmade.(['electrode' num2str(electrode)]) = trial_data_manmade;

        trial_data_natural = struct();
        trial_data_natural.subject = subject;
        trial_data_natural.erp = ERP_matrix_natural;
        all_segments_erp_natural.(['electrode' num2str(electrode)]) = trial_data_natural;
        
        erp_manmade_mean = mean(ERP_matrix_manmade);
        erp_natural_mean = mean(ERP_matrix_natural);
        trial_data_summary = struct();
        trial_data_summary.erp_manmade_mean = erp_manmade_mean - mean(erp_manmade_mean(1:100));    
        trial_data_summary.erp_natural_mean = erp_natural_mean - mean(erp_natural_mean(1:100));
        trial_data_summary.erp_manmade_mean_without_correction = erp_manmade_mean;    
        trial_data_summary.erp_natural_mean_without_correction = erp_natural_mean;    

        all_segments_erp_summary.(['electrode' num2str(electrode)]) = trial_data_summary;
    end
   
save(fullfile(folder_generated_data, [subject_root_name '_manmade']),'all_segments_erp_manmade')
save(fullfile(folder_generated_data, [subject_root_name '_natural']),'all_segments_erp_natural') 

end
