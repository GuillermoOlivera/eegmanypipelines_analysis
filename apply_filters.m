function filtered_data = apply_filters(input_data, sample_rate_hz, high_cutoff_freq)

input_data_as_double = double(input_data);

% high pass
low_cutoff_freq = 0.1; % use 1hz for ICA pre-proces (good results in terms of SNR, better decomposition)
% 1hz too high for ERP analysis? check Luck - ICA / use 0.1 according to Luck and Zakeri
% ("Influence of signal pre-proc on ICA-based EEG decomposition)insted of 0.01 (original version) and 1.
filter_order = 2; % use 2nd order (= 12 dB/octave, Luck)
[b, a] = butter(filter_order,low_cutoff_freq /(sample_rate_hz / 2),'high'); % high pass digital filter design
dataOut_high_pass = filtfilt(b,a,input_data_as_double); % zero-phase filtering
dataOut_high_pass = dataOut_high_pass - mean(dataOut_high_pass);

%%  notch filter
notch_freq = 50;
wo = notch_freq/(sample_rate_hz/2);
qfactor = 35;
bw = wo/qfactor;
[b, a] = iirnotch(wo, bw);
notch_filtered_data = filtfilt(b, a, dataOut_high_pass); % zero-phase filtering
notch_filtered_data = notch_filtered_data - mean(notch_filtered_data); % baseline correct

%%  low pass ideally at 30Hz
filter_order = 3;
[b, a] = butter(filter_order, high_cutoff_freq / (sample_rate_hz / 2)); % low pass digital filter design
dataOut_low_pass = filtfilt(b, a, notch_filtered_data); % zero-phase filtering     % if only low-pass
dataOut_low_pass = dataOut_low_pass - mean(dataOut_low_pass);
data_lp =  dataOut_low_pass; % replace all_data with filtered data

% %%  notch filter
% notch_freq = 50;
% wo = notch_freq/(sample_rate_hz/2);
% qfactor = 35;
% bw = wo/qfactor;
% [b, a] = iirnotch(wo, bw);
% notch_filtered_data = filtfilt(b, a, data_lp); % zero-phase filtering
% notch_filtered_data = notch_filtered_data - mean(notch_filtered_data); % baseline correct

%%
filtered_data = data_lp;

% % low pass
% high_cutoff_freq = 40;
% filter_order = 2;
% [b,a]=butter(filter_order,high_cutoff_freq/(sample_rate/2)); % low pass digital filter design
% dataIn_vEOG = notch_filtered_data_vEOG;
% dataOut_zero_phase_vEOG = filtfilt(b,a,dataIn_vEOG); % zero-phase filtering
% dataOut_zero_phase_vEOG = dataOut_zero_phase_vEOG - mean(dataOut_zero_phase_vEOG);  % baseline correct

%         data(electrode,:) = all_data; % never comment CHECK THIS!!!
end
