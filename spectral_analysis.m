function spectra = spectral_analysis(info, folder, subject_files)

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

use_periodogram = true;

time_window_ms = info.spectra_time_window_ms;
subjects_to_use = info.spectra_subjects_to_use;
selected_electrodes = info.spectra_selected_electrodes;

sample_rate = info.sample_rate_hz;
event_types = info.event_types;
pre_stimulus_ms = info.pre_stimulus_ms;

theta = struct(); theta.min = 4; theta.max = 8;
alpha = struct(); alpha.min = 8; alpha.max = 12;
betha = struct(); betha.min = 12; betha.max = 30;

% Here we reload the data and apply a short time fft or pwelch
% on EACH trial, which we consequently sum to make an average

% In this step it would also make sense to already extract the
% theta and gamma power spectral density

disp('Start spectral analysis high resolution...')

subject_files = ls([folder '/*struct.mat']); % processed
subject_total = size(subject_files, 1);


spectra = struct();
for subj=1:subject_total
    selected_subject_file = subject_files(subj, :)
    load(selected_subject_file);

    % Prepare time window
    time_window_s = [(time_window_ms(1) + erp_data.pre_stimulus_ms) / 1000 (time_window_ms(2) + erp_data.pre_stimulus_ms) / 1000]; % erp_data from loaded file (!)
    if time_window_ms(1) < 0
        savefile = [folder '/spectra.mat']
        first_sample = 0;
        last_sample = 512; % FIXME: hardcoded
    else    
        savefile = [folder '/spectra_' num2str(time_window_ms(1)) '_' num2str(time_window_ms(2)) '.mat']
        first_sample = floor(time_window_s(1) * sample_rate);
        last_sample = ceil(time_window_s(2) * sample_rate);
    end

    fields_event = fieldnames(erp);
    length_event = size(fields_event, 1);

    for event_id=1:length_event
        event_name = fields_event{event_id, :}

        fields_electrode = fieldnames(erp.(event_name));
        length_electrode = size(fields_electrode, 1);
        for electrode_id=1:length_electrode
            electrode_name = fields_electrode{electrode_id, :};
            
            freq = 1:50; % FIXME: check low frequencies ok when interpolating?
            psd_trial_average = zeros(length(freq),1);
            psd_dB_trial_average = zeros(length(freq),1);
            vec_power_total = []; vec_power_band_theta = []; vec_power_band_alpha = []; vec_power_band_betha = []; vec_per_power_theta = []; vec_per_power_alpha = []; vec_per_power_betha = [];

            data_mat = erp.(event_name).(electrode_name);
            length_trials = size(data_mat, 1);
            for trial_id = 1:length_trials
                data = data_mat(trial_id,:);
                data_segment = data(first_sample:last_sample);
                data_len = length(data_segment);

                if use_periodogram
                    [psd, freq_] = periodogram(data_segment, hanning(length(data_segment)), length(data_segment), sample_rate);
                    psd = interp1(freq_, psd, freq);
                else % will currently not work because of incorrectly-sized container vector
                    s_c = fft(data_segment);
                    psd = s_c.*conj(s_c) * 1/(sample_rate * length(data));
                    freq = 0 : sample_rate/data_len : sample_rate/2;
                end

                psd_trial_average = psd_trial_average + psd';
                psd_dB_trial_average = psd_dB_trial_average + 10*log10(psd)';

                power_total = bandpower(psd, freq, 'psd');

                power_band_theta = bandpower(psd, freq, [theta.min theta.max], 'psd');
                power_band_alpha = bandpower(psd, freq, [alpha.min alpha.max], 'psd');
                power_band_betha = bandpower(psd, freq, [betha.min betha.max], 'psd');

                per_power_theta = 100*(power_band_theta/power_total);
                per_power_alpha = 100*(power_band_alpha/power_total);
                per_power_betha = 100*(power_band_betha/power_total);

                vec_power_total = [vec_power_total power_total];
                vec_power_band_theta = [vec_power_band_theta power_band_theta];
                vec_power_band_alpha = [vec_power_band_alpha power_band_alpha];
                vec_power_band_betha = [vec_power_band_betha power_band_betha];
                vec_per_power_theta = [vec_per_power_theta per_power_theta];
                vec_per_power_alpha = [vec_per_power_alpha per_power_alpha];
                vec_per_power_betha = [vec_per_power_betha per_power_betha];

            end
            
            psd_trial_average = psd_trial_average / length_trials;
            spectra.freq = freq; % is the same throughout this analysis
            spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).num_trials = length_trials;

            spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).psd = psd_trial_average;
            spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).psd_dB = psd_dB_trial_average;

            spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).power_total = vec_power_total;
            spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).power_total = vec_power_band_theta;
            spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).power_total = vec_power_band_alpha;
            spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).power_total = vec_power_band_betha;

            %std_per_power_theta = std(vec_per_power_theta);
            %mean_per_power_theta = mean(vec_per_power_theta);
            %std_per_power_alpha = std(vec_per_power_alpha);
            %mean_per_power_alpha = mean(vec_per_power_alpha);
            %std_per_power_betha = std(vec_per_power_betha);
            %mean_per_power_betha = mean(vec_per_power_betha);

            %spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).theta_per_power.std   = std_per_power_theta;
            %spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).theta_per_power.mean  = mean_per_power_theta;
            %spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).alpha_per_power.std   = std_per_power_alpha;
            %spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).alpha_per_power.mean  = mean_per_power_alpha;
            %spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).betha_per_power.std   = std_per_power_betha;
            %spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).betha_per_power.mean  = mean_per_power_betha;

            % Some extra handy variables for plotting
            spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).plot_friendly.f = freq(freq > 5 && freq < 25);
            spectra.(['subject_' num2str(subj)]).(event_name).(electrode_name).plot_friendly.psd = psd(freq > 5 && freq < 25);

            % printf('I: %d, %s, %d, %d, theta_mean = %d, theta_std = %d, alpha_mean = %d alpha_std = %d, betha_mean = %d, betha_std = %d\n\n', ...
            %        subj, event_name, electrode_id, ...
            %        mean_per_power_theta, std_per_power_theta, mean_per_power_alpha, std_per_power_alpha, mean_per_power_betha, std_per_power_betha);
            printf('I: %d, %d, %s, theta_mean = %d, theta_std = %d, alpha_mean = %d alpha_std = %d, betha_mean = %d, betha_std = %d\n', ...
                   subj, electrode_id, event_name, ...
                   mean_per_power_theta, std_per_power_theta,...
                   mean_per_power_alpha, std_per_power_alpha, mean_per_power_betha, std_per_power_betha);
        end

    end

    disp('-------------------------------')
end

spectra_info = info;
printf('Saving spectra to file: %s ... ', savefile)
save(savefile, 'spectra', 'spectra_info', '-mat7-binary');
printf('Done.\n')
end
% psd_mean = mean(psd, 2) / counter;
% data_length = length(psd_mean);

% 
% f = (0 : data_length - 1) * (sample_rate / data_length);     % frequency range
% plot(f, psd_mean); hold on;
