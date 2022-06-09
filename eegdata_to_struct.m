function datastruct = eegdata_to_struct(hypothesis_data, folder, only_log, no_run)

if no_run == true
    return;
end

% Remove all empty vectors and make a struct with data

subject_files = ls([folder '/*.hdf']); % processed
subject_total = size(subject_files, 1);

for subj=1:subject_total
    selected_subject_file = subject_files(subj, :)
    savefile = [selected_subject_file(1:end-4) '_struct.mat'];
    load(selected_subject_file);

    % Convert to struct
    erp = struct();
    ev = electrodes_erp_info.event_types;
    for event_id=1:length(ev)
        event_name = ev{event_id};
        erp_cell = {};
        erp_electrode_struct = struct();

        el = electrodes_erp_info.electrodes_to_use;
        for electrode_id=1:length(el)
            electrode= el(electrode_id);

            counter = 1;
            trials = electrodes_erp_info.num_trials;
            for t = 1:trials
                erp_vector = squeeze(electrodes_erp(event_id, electrode_id, t, :));
                if mean(erp_vector) ~= 0
                    if only_log == false
                        erp_cell{electrode_id, counter} = erp_vector;
                        erp_electrode_struct.(['electrode_' num2str(electrode_id)])(counter,:) = erp_vector;
                    end
                    counter = counter + 1;
                end
            end
            % printf('Subject: %d, Event: %s, Electrode: %d, Trial: %d\n', subj, event_name, electrode_id, counter);
            printf('LOG: %d, %s, %d, %d\n', subj, event_name, electrode_id, counter);
        end

        if only_log == false
            erp.(event_name) = erp_electrode_struct;
        end
    end

    if only_log == false
        save(savefile, 'erp', '-mat7-binary');
    end

    disp('-------------------------------')
end

end