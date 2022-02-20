function result = vmrk_to_mat(vmrk_filepath, result_filepath)
%     vmrk_filepath = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\EMP01.vhdr';
%     result_filepath = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\_generated\';
    
    [filepath, name, extension] = fileparts(vmrk_filepath);
    savefile = fullfile(result_filepath, [name '.out.mat']);
    
    if (exist(savefile, 'file') == 2)
        warning('File already exists, skipping: %s', savefile);
        return;
    end
    
    mkdir(result_filepath); 

    EEG = pop_loadbv(filepath, [name extension], [], []); % open file
    loaded_raw_data_from_eeglab = EEG.data;

    save(savefile, 'loaded_raw_data_from_eeglab', '-v7.3');
    clear EEG;
    clear loaded_raw_data_from_eeglab;
end