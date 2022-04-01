function [length_preprocessed_subjects erp] = mean_erp_for_all_subjects_from_specific_electrode(folder_analysed_data, match_trial_type, specified_electrode)

length_preprocessed_subjects = size(match_trial_type, 1);

a_electrodes = [];
for current_subject = 1:length_preprocessed_subjects
    subject_data = load(fullfile(folder_analysed_data, match_trial_type(current_subject,:)));
    datafield_name = fieldnames(subject_data);
    electrodes_all = subject_data.(datafield_name{1}).all_electrodes;
    electrode_selected = cell2mat(electrodes_all(specified_electrode));
    a_electrodes = [a_electrodes; electrode_selected];
    mean_erp = mean(a_electrodes);
end

erp = mean_erp;

end