function is_match = check_if_condition_match(condition, trial_code)

is_match = false;

trial = struct();
condition_check = struct();
is_match_scene_category = true;
is_match_old = true;
is_match_behaviour = true;
is_match_subsequent_memory = true;

code = trial_code;
trial_code_digits = dec2base(code, 10) - '0';
trial.scene_category = trial_code_digits(1);
trial.old = trial_code_digits(2);
trial.behaviour = trial_code_digits(3);
trial.subsequent_memory = trial_code_digits(4);

condition_digits = dec2base(condition, 10) - '0';
condition_check.scene_category = condition_digits(1);
condition_check.old = condition_digits(2);
condition_check.behaviour = condition_digits(3);
condition_check.subsequent_memory = condition_digits(4);

if (condition_check.scene_category ~= trial.scene_category) && (condition_check.scene_category ~= 6)
    is_match_scene_category = false;
end

if (condition_check.old ~= trial.old) && (condition_check.old ~= 6)
    is_match_old = false;
end

if (condition_check.behaviour ~= trial.behaviour) && (condition_check.behaviour ~= 6)
    is_match_behaviour = false;
end

if (condition_check.subsequent_memory ~= trial.subsequent_memory) && (condition_check.subsequent_memory ~= 6)
    is_match_subsequent_memory = false;
end

is_match = is_match_scene_category && is_match_old && is_match_behaviour && is_match_subsequent_memory;

if (is_match == true)
    return
end

% TODO: return some sanity check debug info to verify logic has no edge cases
end