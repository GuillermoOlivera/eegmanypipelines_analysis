function subject_erp_mean = mean_subjects(data_wrapper, selected_subjects, event_types)

summary = data_wrapper.summary;

for id = 1:size(event_types, 1)
    event_type = event_types{id};
    subject_erp_mean.(event_type) = [];
end

erp_matrix = [];
total_subjects = 0;
for subject = selected_subjects 
    subject_id_str = ['subj_' num2str(subject)];
    for id = 1:size(event_types, 1)
        event_type = event_types{id};
        erp_cell = struct2cell(summary.(subject_id_str).(event_type));
        erp_matrix = cell2mat(erp_cell');

	if isempty(erp_matrix)
	   disp(['Warning (@mean_subjects): Skipping empty matrix. Subject=' num2str(subject) ' type=' event_type]);
	   continue;
	end

        if isempty(subject_erp_mean.(event_type))
            subject_erp_mean.(event_type) = erp_matrix;
        end

        subject_erp_mean.(event_type) = (subject_erp_mean.(event_type) + erp_matrix);

    end
    
    total_subjects = total_subjects + 1;
end

for id = 1:size(event_types, 1)
    event_type = event_types{id};
    subject_erp_mean.(event_type) = subject_erp_mean.(event_type) / total_subjects;
end

end
