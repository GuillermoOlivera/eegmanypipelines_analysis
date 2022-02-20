function filtered_data = apply_filters(input_data)

data_d = double(input_data); 
        
%% Filtering
sample_rate = 512; % acquired at 1000 Hz and downsampled to 512 Hz

%  only low pass at 30Hz
high_cutoff_freq = 30;
filter_order = 3;
[b, a] = butter(filter_order, high_cutoff_freq / (sample_rate / 2)); % low pass digital filter design
dataOut_low_pass = filtfilt(b, a, data_d); % zero-phase filtering     % if only low-pass
dataOut_low_pass = dataOut_low_pass - mean(dataOut_low_pass);
data_lp =  dataOut_low_pass; % replace all_data with filtered data

%  notch filter
notch_freq = 50;
[b, a] = iirnotch(notch_freq / (sample_rate / 2), 0.01);
notch_filtered_data = filtfilt(b, a, data_lp); % zero-phase filtering
notch_filtered_data = notch_filtered_data - mean(notch_filtered_data); % baseline correct

filtered_data = notch_filtered_data;

% % low pass
% high_cutoff_freq = 40;
% filter_order = 2;
% [b,a]=butter(filter_order,high_cutoff_freq/(sample_rate/2)); % low pass digital filter design
% dataIn_vEOG = notch_filtered_data_vEOG;
% dataOut_zero_phase_vEOG = filtfilt(b,a,dataIn_vEOG); % zero-phase filtering
% dataOut_zero_phase_vEOG = dataOut_zero_phase_vEOG - mean(dataOut_zero_phase_vEOG);  % baseline correct
         
 %         data(electrode,:) = all_data; % never comment CHECK THIS!!!
end
