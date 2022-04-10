function result = vmrk_to_mat(vhdr_filepath, result_filepath)
%     vmrk_filepath = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\EMP01.vhdr';
%     result_filepath = 'G:\_EEGManyPipelines\EMP_data\eeg_brainvision\_generated\';
    
    [filepath, name, extension] = fileparts(vhdr_filepath);
    savefile = fullfile(result_filepath, [name '.out.mat']);
    
    if (exist(savefile, 'file') == 2)
        % warning('File already exists, skipping: %s', savefile);
        return;
    end
    
    mkdir(result_filepath); 

    loaded_raw_data_from_eeglab = struct();
    EEG = pop_loadbv(filepath, [name extension], [], []); % open file
    loaded_raw_data_from_eeglab.data = EEG.data;
    sample_rate_hz = 512;
    datapoints = size(loaded_raw_data_from_eeglab.data,2);    
    loaded_raw_data_from_eeglab.times = [1:datapoints];
    loaded_raw_data_from_eeglab.xmin = 0;
    loaded_raw_data_from_eeglab.xmax = datapoints/sample_rate_hz;
    
    save(savefile, 'loaded_raw_data_from_eeglab', '-mat7-binary');
    clear EEG;
    clear loaded_raw_data_from_eeglab;
end