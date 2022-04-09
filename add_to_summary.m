function result = add_to_summary( event_name_str, input_filename, output_filename, folder_analysed_data, electrodes_to_use)

filepath_input = fullfile(folder_analysed_data, [input_filename '.mat']);
filepath_output = fullfile(folder_analysed_data, [output_filename '.mat']);

disp(filepath_input)
session_data = load(filepath_input, 'all_segments_erp_summary');

if (exist(filepath_output, 'file') == 2)
    load(filepath_output, 'all_subject_summary');
else
    all_subject_summary = struct();
end

count = 1;
mean_all_subjects = struct();
erp_matrix = [];
for electrode_selected = electrodes_to_use
    [total_subjects, erp] = erp_subjects_mean(session_data, event_name_str, electrode_selected);
    erp_matrix = [erp_matrix; erp];
    count = count + 1;
end

all_subject_summary.(event_name_str).total_subjects = total_subjects;
all_subject_summary.(event_name_str).electrodes_used = electrodes_to_use;
all_subject_summary.(event_name_str).data = erp_matrix;

save(filepath_output, 'all_subject_summary')

end