function all_triggers = read_triggers_from_vmrk(marker_filename)
    fidInFile = fopen(marker_filename, 'r');    % Open input file for reading
    
    trigger_count = 0;
    
    nextLine = fgets(fidInFile);  % Get the next line of input
    while nextLine >= 0           % Loop until getting -1 (end of file)
        
        if strfind(nextLine, 'Stimulus') > 0  %check if trigger is a stimulus
            
            trigger_count = trigger_count +1;
            comma_locations = strfind(nextLine, ',');  % get the location of the commas
            
            % Description: number_of_digits_for_trigger = 4; -> comma_locations(2)- 4
            %all_triggers(trigger_count,1) =  str2num(nextLine(comma_locations(2)-4:comma_locations(2)-1));  % get trigger name
            %all_triggers(trigger_count,2) =  str2num(nextLine(comma_locations(2)+1:comma_locations(3)-1));  % get trigger time
            
            % for hyp2a
            all_triggers(trigger_count,1) =  str2num(nextLine(comma_locations(2)-4:comma_locations(2)-1));  % get trigger name
            all_triggers(trigger_count,2) =  str2num(nextLine(comma_locations(2)+1:comma_locations(3)-1));  % get trigger time
            
        end
        
        nextLine = fgets(fidInFile);              % Get the next line of input
    end
    
    fclose(fidInFile);
end