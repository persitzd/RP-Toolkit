function is_valid = HPZ_check_file_format_correctness (data)

% this function makes sure that the data matrix read from the file is in
% the required format.
% if it is not in the required format, it may lead either to crashing, or
% erroneous results, or both, therefore this function checks this ahead and
% crashes immediately, with an informative error to the user as what caused
% the problem and how to fix it.

% data - a matrix that its 1st column is subject number and the 2nd is the
% observation number, the remaining columns are quantities and prices.

% initialization
is_valid = true;

% a list of subjects. we could have get it using the command: 
subjects = unique(data(:,1), 'stable');
num_subjects = length(subjects);

% total number of observations in all subjects
[total_obs , ~] = size(data);

% initializations
current_subject = 1;
current_observation = 1;

% checking validity is in this loop
for i=1:total_obs
    
    if data(i,2) == current_observation
        % if it is the same observation number, then we are still working
        % on the same subject
        if data(i,1) == subjects(current_subject)
            % then it's ok, all clear
            current_observation = current_observation + 1;
        else
            is_valid = false;
            error_str = char(strcat('Data Format Error - the data file is not formatted correctly. In row', ...
                                    {' '}, num2str(i), {' '}, ...
                                    'in the file, the number of observation continues to grow, while the subject number (', ...
                                    num2str(subjects(current_subject)), ...
                                    ') has changed. If it is a new subject, observation number should equal 1.'));
            %error(error_str);
            msgbox(error_str);
            return
        end
    elseif data(i,2) == 1
        % then we moved to the next subject - let's make sure this is what
        % actually happened
        if data(i,1) == subjects(current_subject)
            % the observation turned to 1, but the subject number still hadn't changed
            is_valid = false;
            error_str = char(strcat('Data Format Error - the data file is not formatted correctly. In row', ...
                                    {' '}, num2str(i), {' '}, ...
                                    'in the file, the number of observation turned to 1, while the subject number (', ...
                                    num2str(subjects(current_subject)), ...
                                    ') has not changed. If it is still the same subject,', ...
                                    {' '}, 'observation number should be the consecutive of the previous one.'));
            %error(error_str);
            msgbox(error_str);
            return
        elseif current_subject == num_subjects || data(i,1) ~= subjects(current_subject+1)
            % we either already went threw all the subjects, or the subject
            % next to the previous subject is not what's supposed to be.
            % one of these things is due to happen if a subject appears
            % twice in the data in 2 locations, e.g. in rows 120-155 and
            % again in rows 320-355 (with the exact same subject number).
            is_valid = false;
            error_str = char(strcat('Data Format Error - the data file is not formatted correctly. In row', ...
                                    {' '}, num2str(i), {' '}, ...
                                    'in the file, there is a new subject number (', ...
                                    num2str(data(i,1)), ...
                                    '), but it seems that this subject number already appeared earlier in the file,', {' '}, ...
                                    'meaning: this subject was duplicated, or there are two subjects with the same ID number.'));
            %error(error_str);
            msgbox(error_str);
            return
        else
            % then it's ok, all clear, we prepare for the next iteration
            current_observation = 2;
            current_subject = current_subject + 1;
        end
    else
        % observation number is not the consecutive number of the previous
        % and not 1 - this can't be
        if i==1
            is_valid = false;
            error_str = char(strcat('Data Format Error - the data file is not formatted correctly. In row 1', ...
                                    {' '}, 'in the file, the number of observation is not 1 as required.'));
        else
            is_valid = false;
            error_str = char(strcat('Data Format Error - the data file is not formatted correctly. In row', ...
                                    {' '}, num2str(i), {' '}, ...
                                    'in the file, the number of observation (', num2str(data(i,2)), ...
                                    ') is neither the consecutive of the previous (', num2str(data(i-1,2)), ...
                                    ') nor is it 1 (for if it was a new subject).'));
        end
        %error(error_str);
        msgbox(error_str);
        return
    end
    
end

end