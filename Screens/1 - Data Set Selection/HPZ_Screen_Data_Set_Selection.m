function [data_matrix, data_pref_class, data_choice_set_type, ok] = HPZ_Screen_Data_Set_Selection(main_folder, runs_counter)

% this function promotes a user-interface screen to the user, allowing her
% to choose a data on which to perform consistency checks or parameters 
% estimation.
% it also allows the user to choose whether to automatically fix the
% endowments to be precisely 1, to avoid miscalculations.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

% it is recommended to use the numeric variables in the beginning of the 
% function when making changes to the screen.
% it is also recommended to use numeric variables in the same manner when
% adding new elements to the screen.



% initialization to prevent erros
data_matrix = 0;
data_pref_class = 0;
data_choice_set_type = 0;
ok = 0;

%% get the data settings
[data_list_str, data_list_path, data_list_prefs, data_list_choice_set_types, data_list_subject, data_list_obs, data_list_quantity1, data_list_quantity2, data_list_maxquantity1, data_list_maxquantity2, data_set, fix_endowments] = HPZ_Data_Set_Selection_Settings_Read(main_folder);

% We check if all these files actually exist (in the defined paths).
% For example, it could be that the file's location was changed, or that
% the program package was transferred to another computer
% We delete the files that we don't manage to find from the list.
list_size = length(data_list_path);
file_exists = zeros(1, list_size);
for i=1:list_size
    if exist(data_list_path{i}, 'file')
        file_exists(i) = 1;
    end
end

% creating the new list (without the deleted items),
% for the list of names and the list of paths and the list of locations
helper_vector = 1:list_size;
helper_file_exists = helper_vector .* file_exists;
data_list_str = data_list_str(helper_vector == helper_file_exists);
data_list_path = data_list_path(helper_vector == helper_file_exists);
data_list_prefs = data_list_prefs(helper_vector == helper_file_exists);
data_list_choice_set_types = data_list_choice_set_types(helper_vector == helper_file_exists);
data_list_subject = data_list_subject(helper_vector == helper_file_exists);
data_list_obs = data_list_obs(helper_vector == helper_file_exists);
data_list_quantity1 = data_list_quantity1(helper_vector == helper_file_exists);
data_list_quantity2 = data_list_quantity2(helper_vector == helper_file_exists);
data_list_maxquantity1 = data_list_maxquantity1(helper_vector == helper_file_exists);
data_list_maxquantity2 = data_list_maxquantity2(helper_vector == helper_file_exists);
if helper_file_exists(data_set)
    helper_vector = helper_vector(helper_vector == helper_file_exists);
    data_set = find(helper_vector == data_set, 1);
elseif sum(helper_file_exists) == 0
    data_set = 0;
else
    data_set = 1;
end


% this little cell helps us to convert 0 and 1 to 'off' and 'on', respectively 
enable = {'off','on'};



%% location and size parameters (use these to easily make changes to the screen)

% the height of the list inside the screen
list_height = min(max(56 , 13.5*(max(size(data_list_str))+1)), 250); 

% screen size
figure_width = 450;
figure_height = 320 + list_height;

% limit of height as percentage of computer screen height
max_height_percent = HPZ_Constants.max_height_percent;

% height of the label in the head and of the more options label
label_height = 20;
% distance of highest element from top
top_dist = 15;

% height of almost-bottom space designated for OK and Cancel buttons
buttons_space_height = 50;
% buttons height as percentage of the space they are in
buttons_height = 0.6;
buttons_height_more_options = 0.8;
% the distant between a button to the edge of the screen,
% and is also half the distance between buttons.
% the size of the buttons is designed to fit this and the buttons_height 
% and the width of the screen
buttons_dists = 10;
% number of buttons
buttons_num = 2;

% height of bottom space designated for buttons that lead to advanced options 
advanced_options_space_height = 40;
% total height of "more options" section
more_options_total_height = 2*(top_dist/2) + advanced_options_space_height;

% width of each element
element_width = 320;
% normal_width for a sub element
sub_element_width = 80;

% relative location of elements:
% how much distance between parts (different indexes) that are one below each other 
move_down = 15;

% each radio options will be with height 20
radio_height = 20;
% distance of  radio option from the bottom of the radio frame
radio_bottom = 5;
% offsets to the right of yes/no horizontal radio options
yes_no_offsets = [25 , 150];

% font size of label
label_font_size = 12;
% general font size
font_size = 8;
% general but bigger font size
big_font_size = 10;





%% create the figure, with a slider if needed
% scroll bar width (if scroll is needed)
scroll_width = 20;
bottom_space_height = buttons_space_height + more_options_total_height;
top_space_height = 0;
panel_height = figure_height - bottom_space_height - top_space_height;
figure_title = char(strcat('Dataset Selection (', HPZ_Constants.current_run_screen, {' '}, num2str(runs_counter), ')'));
[fh , panel] = ui_scroll_screen(figure_width, figure_height, scroll_width, max_height_percent, top_space_height, bottom_space_height, figure_title);

more_options_panel = uipanel('Parent',fh,...
                'backgroundc',get(fh,'color'),...
                'units','pix',...
                'position',[0 , 0 , figure_width , more_options_total_height]);


% width including scroll bar if there is one
pos = get(fh,'position');
full_width = pos(3);





%% Head Label
% current bottom coordinate
current_bottom = panel_height - top_dist - label_height;
label_DA = uicontrol('Parent',panel, ...
    'style','text',...
    'units','pix',...
    'position',[0 , current_bottom , figure_width , label_height],...
    'backgroundc',get(fh,'color'),...
    'fontsize',label_font_size,'fontweight','bold',...
    'string','Dataset Selection'); %#ok<NASGU>



%% Select Dataset
% current height of element
current_height = 70 + list_height;
% current bottom coordinate
current_bottom = current_bottom - move_down - current_height;
% a panel to contain the list and the buttons
dataset_panel = uipanel('Parent',panel,...
                'backgroundc',get(fh,'color'),...
                'units','pix',...
                'position',[0 , current_bottom , figure_width+2 , current_height]);

% List of Dataset Files to select from
lb_DS = uicontrol('Parent',dataset_panel,...
    'value', data_set,...
    'style','listbox',...
    'enable','on',...
    'unit','pix',...
    'position',[figure_width/10 , buttons_space_height , figure_width*8/10 , list_height],...
    'fontsize',font_size,...
    'string',data_list_str);

% Button to Remove file from list Button
remove_file = uicontrol('Parent',dataset_panel, ...
    'style','push',...
    'unit','pix',...
    'position',[buttons_dists , (1-buttons_height)/2*buttons_space_height , (figure_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_height*buttons_space_height],...
    'backgroundc',[255,184,184]/255,...
    'fontsize',font_size,...
    'string','Remove File',...
    'callback',{@remove_file_call});

% button to Browse for a new file Button
browse_file = uicontrol('Parent',dataset_panel, ...
    'style','push',...
    'unit','pix',...
    'position',[(figure_width/buttons_num)+buttons_dists , (1-buttons_height)/2*buttons_space_height , (figure_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_height*buttons_space_height],...
    'backgroundc',[184,255,184]/255,...
    'fontsize',font_size,...
    'string','Add File...',...
    'callback',{@browse_file_call}); %#ok<NASGU>





%% Fix Endowments (or not)
% current height of element
current_height = 50;
% current bottom coordinate
current_bottom = current_bottom - move_down - current_height;
bg_ff = uibuttongroup(panel,'units','pix',...
    'title', 'Fix Quantities so Endowments will be Exactly 1', ...
    'fontsize',font_size,...
    'pos',[(figure_width-element_width)/2 , current_bottom , element_width , current_height]);

fix_endowments_rd(1) = uicontrol(bg_ff,...
    'value', (fix_endowments == 0),...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_height],...
    'fontsize',font_size,...
    'string',' No');
fix_endowments_rd(2) = uicontrol(bg_ff,...
    'value', (fix_endowments ~= 0),...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_height],...
    'fontsize',font_size,...
    'string',' Yes');

% help button
fix_endowments_help_button = uicontrol('Parent',panel, 'style','push', ...
    'unit','pix', ...
    'position',[(figure_width-element_width)/2+element_width-50 , current_bottom+10 , 20 , 20], ...
    'string','?', ...
    'fontsize',big_font_size, ...
    'callback',{@fix_endowments_help_button_call}); %#ok<NASGU>





%% OK Button
ok_button = uicontrol('Parent',fh, 'style','push',...
    'enable',enable{(max(size(data_list_str))>0)+1},...
    'unit','pix',...
    'position',[buttons_dists , (1-buttons_height)/2*buttons_space_height + more_options_total_height , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_height*buttons_space_height],...
    'string','OK',...
    'fontsize',big_font_size,...
    'callback',{@ok_button_call});

%% Cancel Button
cancel_button = uicontrol('Parent',fh, 'style','push',...
    'enable','on',...
    'unit','pix',...
    'position',[(full_width/2)+buttons_dists , (1-buttons_height)/2*buttons_space_height + more_options_total_height , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_height*buttons_space_height],...
    'string','Cancel',...
    'fontsize',big_font_size,...
    'callback',{@cancel_button_call}); %#ok<NASGU>



%% Button for entering Advanced Settings
advanced_settings_button = uicontrol('Parent',more_options_panel, 'style','push',...
    'enable','on',...
    'unit','pix',...
    'position',[buttons_dists , (1-buttons_height_more_options)/2*advanced_options_space_height + top_dist/2 , (full_width-2*buttons_dists) , buttons_height_more_options*advanced_options_space_height],...
    'string','Advanced Settings',...
    'fontsize',font_size,...
    'callback',{@advanced_settings_button_call}); %#ok<NASGU>





uiwait(fh)  % Prevent all other processes from starting until closed.





%% for fix endowments help button
function [] = fix_endowments_help_button_call(varargin)
    msgbox({'Whether to perform a fixing procedure to make sure that the chosen bundles are exactly on the budget line.', ...
            'It is crucial, as this program always performs its calculations under this assumption.', ...
            'In some experimental software there is inaccuracy in the reported quantities of the bundles chosen by the subjects, and as a result, the bundles reported in the file data are slightly above or under the budget line.', ...
            'The fixing procedure changes both quantities (x1,x2) by multiplying both by the same number (by: 1 / (x1*p1 + x2*p2)).', ...
            'This fixing procedure ensures that corner bundles will remain this way, that bundles with equal quantities (on the 45 degrees line) will stay this way, and in general, that the ratio between the goods before and after fixing will be the same.', ...
            'All the explanations above were made for the case of 2 goods for simplicity, but the procedure takes place with any number of goods (by multiplying each quantity by: 1 / (x1*p1 + ... + xn*pn)).', ...
            '', 'If you wish the procedure not to be used, select the "NO" option.', ...
            'If you wish to fix the quantities in another manner, do it in the data file itself, and only then read it and use this program on it.', ...
            '', 'Note that if the deviation from the budget line is of more than 5% (in even one observation in the dataset), the program will refuse to use the data file, as in this case the data file is assumed to contain significant errors.'});
end



%% Program for remove file button
function [] = remove_file_call(varargin)
    % Callback for pushbutton.
    
    % current list
    %current_list = get(lb_DS,'string');
    % current chosen index in list (will be deleted)
    data_set = get(lb_DS,'value');
    % current size of list
    list_size = max(size(data_list_str));
    if list_size == 1
        % then the list will be empty after the delete - 
        % removing file and pressing ok should be disabled
        set(remove_file,'enable', 'off');
        set(ok_button,'enable', 'off');
    end
    if list_size == data_set
        % then the current value will be outside the bounds of the list
        set(lb_DS,'value',data_set-1);
    end
    
    % creating the new list (without the deleted item),
    % for the list of names and the list of paths and the list of locations
    
    helper_vector = 1:list_size;
    
    data_list_str = data_list_str(helper_vector ~= helper_vector(data_set));
    data_list_path = data_list_path(helper_vector ~= helper_vector(data_set));
    data_list_prefs = data_list_prefs(helper_vector ~= helper_vector(data_set));
    data_list_choice_set_types = data_list_choice_set_types(helper_vector ~= helper_vector(data_set));
    data_list_subject = data_list_subject(helper_vector ~= helper_vector(data_set));
    data_list_obs = data_list_obs(helper_vector ~= helper_vector(data_set));
    data_list_quantity1 = data_list_quantity1(helper_vector ~= helper_vector(data_set));
    data_list_quantity2 = data_list_quantity2(helper_vector ~= helper_vector(data_set));
    data_list_maxquantity1 = data_list_maxquantity1(helper_vector ~= helper_vector(data_set));
    data_list_maxquantity2 = data_list_maxquantity2(helper_vector ~= helper_vector(data_set));

    % updating the list
    set(lb_DS,'string', data_list_str);
end



%% Program for browse for file button
function [] = browse_file_call(varargin)
    % Callback for pushbutton.
    
    % browse for a data file
    [FileName, PathName, ~] = uigetfile('../*.csv', 'Select a CSV Data File', strcat(main_folder, '/', HPZ_Constants.data_files_dir));
    
    % if the user chose a file (didn't cancel)
    if (FileName ~= 0)
         
        [data_name, choice_set_type, locations, pref_class, ok] = HPZ_Screen_Data_Set_Adding(FileName);
        if (ok == 1)
            % then the list will not be empty after the addition - 
            % removing file and pressing ok should be disabled
            set(remove_file,'enable', 'on');
            set(ok_button,'enable', 'on');

            % current list
            %current_list = get(lb_DS,'string');
            % current size of list
            list_size = max(size(data_list_str));
            % if the list is empty, we will get 1 while we want to get 0 
            if isempty(data_list_str)
                list_size = 0;
            end
            
            
            % creating a new list (with the new item),
            % for the list of names and the list of paths and the list of locations 
            data_list_str{list_size+1} = data_name;
            data_list_path{list_size+1} = strcat(PathName, FileName);
            data_list_prefs{list_size+1} = pref_class;
            data_list_choice_set_types{list_size+1} = choice_set_type;
            
            if choice_set_type == HPZ_Constants.choice_set_finite_set
               locations = 1:6;   % just to prevent an error 
            end
            data_list_subject{list_size+1} = locations(1);
            data_list_obs{list_size+1} = locations(2);
            data_list_quantity1{list_size+1} = locations(3);
            data_list_quantity2{list_size+1} = locations(4);
            data_list_maxquantity1{list_size+1} = locations(5);
            data_list_maxquantity2{list_size+1} = locations(6);
            
            % updating the list
            set(lb_DS,'string', data_list_str);
            
            % automatically have the new data file be the selected one
            set(lb_DS,'value',list_size+1);
        end
    end
end





%% Program for OK button
function [] = ok_button_call(varargin)
    % Callback for OK pushbutton.

    % Setting the dataset selection
    data_set = get(lb_DS, 'value');

    % Setting whether to fix endwoments
    switch findobj(get(bg_ff,'selectedobject'))
        case fix_endowments_rd(1)
            % No
            fix_endowments = 0;

        case fix_endowments_rd(2)
            % Yes
            fix_endowments = 1;
    end

    
    file_path = data_list_path{data_set};
    locations = [data_list_subject{data_set},...
                data_list_obs{data_set},...
                data_list_quantity1{data_set},...
                data_list_quantity2{data_set},...
                data_list_maxquantity1{data_set},...
                data_list_maxquantity2{data_set}];
    data_choice_set_type = data_list_choice_set_types{data_set};
    
    if data_choice_set_type == HPZ_Constants.choice_set_budget_line
        [data_matrix, success, is_valid] = HPZ_Data_Format (file_path, locations);
    elseif data_choice_set_type == HPZ_Constants.choice_set_finite_set
        [data_matrix, success, is_valid] = HPZ_Data_Format_Finite_Set (file_path);
    end
    
    if (~success)
        % the file does not exist or is not available - we end this run
        msgbox(char(strcat(HPZ_Constants.could_not_read_file_1, {' '}, file_path, HPZ_Constants.could_not_read_file_2)));
        % remove the file from the list
        remove_file_call();
    elseif (~is_valid)
        % remove the file from the list (on second thought - leave it, but don't allow the user to continue with this file until issues in the file will be solved)    
        %remove_file_call();
    else
        if data_choice_set_type == HPZ_Constants.choice_set_budget_line
            
            % check if endowments are (approximately) equal to 1, otherwise print a warning 
            % the endowments (may not equal to 1)
            %endowments = data_matrix(:,3) .* data_matrix(:,5) + data_matrix(:,4) .* data_matrix(:,6);
            [~, num_of_columns] = size(data_matrix);
            num_of_goods = (num_of_columns - 2) / 2;
            endowments = sum( data_matrix(:, (2+1):(2+num_of_goods)) .* data_matrix(:, (2+num_of_goods+1):(2+num_of_goods+num_of_goods)) , 2 );
            % number of endowment significantly different from 1
            num_of_errors = sum ( (endowments > 1+HPZ_Constants.fix_endowments_error) | (endowments < 1/(1+HPZ_Constants.fix_endowments_error)) );
            if (num_of_errors)
                % the file does not exist or is not available - we end this run
                msgbox(char(strcat({'There were '}, num2str(num_of_errors), {' cases that the chosen bundle was significantly not on the budget line (the expenditure of the chosen bundle is more than '}, num2str(1+HPZ_Constants.fix_endowments_error), {' or less than '}, num2str(1/(1+HPZ_Constants.fix_endowments_error)), {' times the budget constraint). Please fix your data such that all observations will be (at least approximately) on the budget line.'})));
                % remove the file from the list (on second thought - leave it, but don't allow the user to continue with this file until issues in the file will be solved)   
                %remove_file_call();
            else
                % fixing the bundles so their endowments will be equal exactly to 1.
                % this fix is needed because of rounding that makes the endowment a little
                % different from 1, which may damage the estimations to some extent
                if fix_endowments && data_choice_set_type == HPZ_Constants.choice_set_budget_line
                    % if the user chose to fix all bundles so their endwoments will be exactly 1, 
                    % if the difference (in percentage) between 1 and the endowment is more
                    % than this threshold, a warning will be given
                    data_matrix = HPZ_Fix_Endowments_To_One(data_matrix, 1);
                end

                data_pref_class = data_list_prefs{data_set};
                ok = 1;

                HPZ_Data_Set_Selection_Settings_Write(data_list_str, data_list_path, data_list_prefs, data_list_choice_set_types, data_list_subject, data_list_obs, data_list_quantity1, data_list_quantity2, data_list_maxquantity1, data_list_maxquantity2, data_set, fix_endowments, main_folder);
                % close the window
                close(fh);
            end
            
        elseif data_choice_set_type == HPZ_Constants.choice_set_finite_set
            
            data_pref_class = data_list_prefs{data_set};
            ok = 1;

            HPZ_Data_Set_Selection_Settings_Write(data_list_str, data_list_path, data_list_prefs, data_list_choice_set_types, data_list_subject, data_list_obs, data_list_quantity1, data_list_quantity2, data_list_maxquantity1, data_list_maxquantity2, data_set, fix_endowments, main_folder);
            % close the window
            close(fh);
            
        end

    end

end



%% Program for Cancel button
function [] = cancel_button_call(varargin)
    % Callback for Cancel pushbutton.

    ok = 0;

    % close the window
    close(fh);

end




%% Program for Advanced Settings button
function [] = advanced_settings_button_call(varargin)

    HPZ_Screen_Advanced_Options(main_folder);

end



end

