function subject_summary_filename = postprocess(folder_analysed_data, subjects_to_use, selected_electrodes, event_types, subject_files)
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

summary = struct();

subject_summary_filename = fullfile(folder_analysed_data, 'summary.hdf');
isSummarized = exist(subject_summary_filename, 'file') == 2;

if isSummarized
    disp('Skipping processing -> calculate mean of all trials per subject, electrode and event');
end

% Memory allocation
for subject = subjects_to_use
    % skip if we have already the summary file
    if isSummarized
        continue
    end

    subject_id_str = ['subj_' num2str(subject)];
    for electrode = selected_electrodes
        electrode_id_str = ['ch_' num2str(electrode)];
        for id = 1:size(event_types, 1)
            event_type = event_types{id};
            summary.(subject_id_str).(event_type).(electrode_id_str) = [];
        end
    end
end

for subject = subjects_to_use
    % skip if we have already the summary file
    if isSummarized
        continue
    end

    disp(['Summarizing subject: ' num2str(subject)]);
    subject_id_str = ['subj_' num2str(subject)];
    subject_root_vhdr_name = subject_files(subject, :);
    subject_root_name = erase(subject_root_vhdr_name, '.vhdr');
    if isOctave
        [dir, subject_root_name, ext] = fileparts(subject_root_name);
    end

    subject_analysis_filename = fullfile(folder_analysed_data, [subject_root_name '_analysed.hdf']);
    load(subject_analysis_filename); % variable name: electrodes_erp

    for id = 1:size(event_types, 1)
        event_type = event_types{id};
        erp_matrix = [];
        for electrode = selected_electrodes
            electrode_id_str = ['ch_' num2str(electrode)];
            erp_matrix = mean(electrodes_erp.(event_type).(electrode_id_str), 2);

            if ~isempty(erp_matrix)
                summary.(subject_id_str).(event_type).(electrode_id_str) = erp_matrix;
            else
                disp(["Warning: no valid ERPs found. subject=" num2str(subject) ", event=" event_type]);
            end
        end
    end
end

if ~isSummarized
    save(subject_summary_filename, 'summary', '-mat7-binary') % in Octave can also use hdf5 output format '-hdf5'
end




end