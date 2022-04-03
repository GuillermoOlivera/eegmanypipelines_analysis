function [total_subjects, erp] = erp_subjects_mean(session_data, event_name_str, specified_electrode)

a_electrodes = [];
subject_names = fieldnames(session_data.all_segments_erp_summary.session);
total_subjects = size(subject_names, 1);
for subject = 1:total_subjects
    electrode_data = session_data.all_segments_erp_summary.session.(subject_names{subject}).(['channel_' num2str(specified_electrode)]).(event_name_str);
    a_electrodes = [a_electrodes; electrode_data];
    mean_erp = mean(a_electrodes);
end

erp = mean_erp;

end