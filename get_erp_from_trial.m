function [erp_manmade, erp_natural] = function get_erp_from_trial(all_triggers, data, vEOG_data)

stimulus_time = all_triggers(trial,2);
channel_segment = data(electrode,stimulus_time-100:stimulus_time+299); % original_ not electrode but channel!!
channel_segment = channel_segment - mean(channel_segment(1:100));

vEOG_segment = vEOG_data(stimulus_time-100:stimulus_time+299); % get EOG segment
vEOG_segment = vEOG_segment - mean(vEOG_segment(1:100)); % baseline correct

hEOG_segment = hEOG_data(stimulus_time-100:stimulus_time+299);
hEOG_segment = hEOG_segment - mean(hEOG_segment(1:100));

%% check there is no eye blink or artifact on any of the channels
% peak-to-peak voltage vs threshold voltage
% peak-to-peak measure less distorted by slow changes in baseline voltage, reducing the impact of source of
% misses and false alarms, and increasing the sensitivity of the artifact rejection process (see Luck'book)
segment_ok = 1;

channel_segment = data(electrode,stimulus_time-100:stimulus_time+299); % original_ not electrode but channel!!
if (max(channel_segment(1:400)) - min(channel_segment(1:400))) > 200 ||...
        ((max(vEOG_segment(1:400)) - min(vEOG_segment(1:400))) > 160) ||...
        ((max(hEOG_segment(1:400)) - min(hEOG_segment(1:400))) > 160) % Paul, 160(80),120(60),60(30)/Luck,200(100),160(80),160(80)
    
    segment_ok = 0; % not OK                       
end

% only take the good segments
if  segment_ok == 1
    num_events_ok = num_events_ok + 1; % count number of segments
    
    if condition == 1
        ERP_matrix_manmade(num_events_ok,:) = channel_segment; % put the good segment into matrix
    end
    
    if condition == 2
        ERP_matrix_natural(num_events_ok,:) = channel_segment; % put the good segment into matrix
    end
end
        
end
