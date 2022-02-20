%% save eegdata for each subject

clear all
clc

function result = vmrk_to_mar()

folder = 'C:\Users\eeg\Desktop\eeg';

subject_names = {'EMP01'};

subjects_to_use = 1;

num_subjects = length(subjects_to_use);

for subject = subjects_to_use
    
    subject_name =  char(subject_names(subject));
    
    folder_name = [folder '\subj_' num2str(subject)];
    
    file_name = fullfile(folder_name,subject_name);
    
    %% load EEG data
    
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % EEG lab
    EEG = pop_loadbv(folder_name, [subject_name '.vhdr'], [], []); % open file
    [ALLEEG EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
    
    data = EEG.data; 
   
    save([num2str(subject_name) '_data'],'data') 
   
end

end