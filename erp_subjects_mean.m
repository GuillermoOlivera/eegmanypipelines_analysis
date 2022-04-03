function [length_preprocessed_subjects, erp] = erp_subjects_mean(folder_analysed_data, event_name_str, match_trial_type, specified_electrode)

length_preprocessed_subjects = size(match_trial_type, 1);

a_electrodes = [];
for current_subject = 1:length_preprocessed_subjects
    subject_data = load(fullfile(folder_analysed_data, match_trial_type(current_subject,:)));
    datafield_name = fieldnames(subject_data);
    electrode_data = subject_data.(datafield_name{1}).data.(['channel_' num2str(specified_electrode)]).(event_name_str);
    a_electrodes = [a_electrodes; electrode_data];
    mean_erp = mean(a_electrodes);
end

erp = mean_erp;

end