function is_match = check_if_condition_match()

event_type = struct();
event_type.scene_category = 6
event_type.old = 6
event_type.behaviour = 6
event_type.subsequent_memory = 6

test_trigger_code = 1029;


end


%
% 1. get code
% 2. check if code equals desired conditions
% 
% - how to define what condition I want... just use a struct and functions "knows" everything to check if match 
%         