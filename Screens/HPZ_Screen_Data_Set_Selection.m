function [data_matrix, data_pref_class, fix_endowments, ok] = HPZ_Screen_Data_Set_Selection(fix_endowments, main_folder, runs_counter)

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
ok = 0;

%% get the data settings
[data_list_str, data_list_path, data_list_prefs, data_list_subject, data_list_obs, data_list_quantity1, data_list_quantity2, data_list_maxquantity1, data_list_maxquantity2, data_set] = HPZ_Data_Settings_Read(main_folder);

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
data_list_subject = data_list_subject(helper_vector == helper_file_exists);
data_list_obs = data_list_obs(helper_vector == helper_file_exists);
data_list_quantity1 = data_list_quantity1(helper_vector == helper_file_exists);
data_list_quantity2 = data_list_quantity2(helper_vector == helper_file_exists);
data_list_maxquantity1 = data_list_maxquantity1(helper_vector == helper_file_exists);
data_list_maxquantity2 = data_list_maxquantity2(helper_vector == helper_file_exists);



% this little cell helps us to convert 0 and 1 to 'off' and 'on', respectively 
enable = {'off','on'};



%% location and size parameters (use these to easily make changes to the screen)

% the height of the list inside the screen
list_height = min(max(56 , 14*(max(size(data_list_str))+1)), 210); 

% screen size
figure_width = 450;
figure_height = 260 + list_height;

% limit of height as percentage of computer screen height
max_height_percent = HPZ_Constants.max_height_percent;

% height of bottom space designated for OK and Cancel buttons
buttons_space_height = 50;
% buttons height as percentage of the space they are in
buttons_height = 0.6;
% the distant between a button to the edge of the screen,
% and is also half the distance between buttons.
% the size of the buttons is designed to fit this and the buttons_height 
% and the width of the screen
buttons_dists = 10;
% number of buttons
buttons_num = 2;

% width of each element
element_width = 320;
% normal_width for a sub element
sub_element_width = 80;

% relative location of elements:
% how much distance between parts (different indexes) that are one below each other 
move_down = 15;

% distance of highest element from top
top_dist = 15;

% height of the label in the head
label_height = 20;

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
%big_font_size = 10;





%% create the figure, with a slider if needed
% scroll bar width (if scroll is needed)
scroll_width = 20;
bottom_space_height = buttons_space_height;
top_space_height = 0;
panel_height = figure_height - bottom_space_height - top_space_height;
figure_title = strcat('Dataset Selection (', HPZ_Constants.current_run_screen, num2str(runs_counter), ')');
[S.fh , S.panel] = ui_scroll_screen(figure_width, figure_height, scroll_width, max_height_percent, top_space_height, bottom_space_height, figure_title);

% width including scroll bar if there is one
pos = get(S.fh,'position');
full_width = pos(3);





%% Head Label
% current bottom coordinate
current_bottom = panel_height-top_dist - label_height;
S.label_DA = uicontrol('Parent',S.panel, ...
    'style','text',...
    'units','pix',...
    'position',[0 , current_bottom , figure_width , label_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',label_font_size,'fontweight','bold',...
    'string','Dataset Selection');



%% Select Dataset
% current height of element
current_height = 70 + list_height;
% current bottom coordinate
current_bottom = current_bottom - move_down - current_height;
% a panel to contain the list and the buttons
S.dataset_panel = uipanel('Parent',S.panel,...
                'backgroundc',get(S.fh,'color'),...
                'units','pix',...
                'position',[0 , current_bottom , figure_width+2 , current_height]);

% List of Dataset Files to select from
S.lb_DS = uicontrol('Parent',S.dataset_panel,...
    'value', data_set,...
    'style','listbox',...
    'enable','on',...
    'unit','pix',...
    'position',[figure_width/10 , buttons_space_height , figure_width*8/10 , list_height],...
    'fontsize',font_size,...
    'string',data_list_str);

% Button to Remove file from list Button
S.remove_file = uicontrol('Parent',S.dataset_panel, ...
    'style','push',...
    'unit','pix',...
    'position',[buttons_dists , (1-buttons_height)/2*buttons_space_height , (figure_width-(buttons_num+1)*buttons_dists)/2 , buttons_height*buttons_space_height],...
    'backgroundc',[255,184,184]/255,...
    'fontsize',font_size,...
    'string','Remove File',...
    'callback',{@remove_file_call,S});

% button to Browse for a new file Button
S.browse_file = uicontrol('Parent',S.dataset_panel, ...
    'style','push',...
    'unit','pix',...
    'position',[(figure_width/2)+buttons_dists , (1-buttons_height)/2*buttons_space_height , (figure_width-(buttons_num+1)*buttons_dists)/2 , buttons_height*buttons_space_height],...
    'backgroundc',[184,255,184]/255,...
    'fontsize',font_size,...
    'string','Add File...',...
    'callback',{@browse_file_call,S});





%% Fix Endowments (or not)
% current height of element
current_height = 60;
% current bottom coordinate
current_bottom = current_bottom - move_down - current_height;
S.bg_ff = uibuttongroup(S.panel,'units','pix',...
    'title', 'Fix Quantities so Endowments will be Exactly 1', ...
    'fontsize',font_size,...
    'pos',[(figure_width-element_width)/2 , current_bottom , element_width , current_height]);

S.fix_endowments_rd(1) = uicontrol(S.bg_ff,...
    'value',1-fix_endowments,...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_height],...
    'fontsize',font_size,...
    'string',' No');
S.fix_endowments_rd(2) = uicontrol(S.bg_ff,...
    'value',fix_endowments,...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_height],...
    'fontsize',font_size,...
    'string',' Yes');





%% OK Button
S.ok_button = uicontrol('Parent',S.fh, 'style','push',...
    'enable',enable{(max(size(data_list_str))>0)+1},...
    'unit','pix',...
    'position',[buttons_dists , (1-buttons_height)/2*buttons_space_height , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_height*buttons_space_height],...
    'string','OK',...
    'fontsize',font_size,...
    'callback',{@ok_button_call,S});

%% Cancel Button
S.cancel_button = uicontrol('Parent',S.fh, 'style','push',...
    'enable','on',...
    'unit','pix',...
    'position',[(full_width/2)+buttons_dists , (1-buttons_height)/2*buttons_space_height , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_height*buttons_space_height],...
    'string','Cancel',...
    'fontsize',font_size,...
    'callback',{@cancel_button_call,S});





uiwait(S.fh)  % Prevent all other processes from starting until closed.





%% Program for remove file button
function [] = remove_file_call(varargin)
    % Callback for pushbutton.
    
    % current list
    %current_list = get(S.lb_DS,'string');
    % current chosen index in list (will be deleted)
    data_set = get(S.lb_DS,'value');
    % current size of list
    list_size = max(size(data_list_str));
    if list_size == 1
        % then the list will be empty after the delete - 
        % removing file and pressing ok should be disabled
        set(S.remove_file,'enable', 'off');
        set(S.ok_button,'enable', 'off');
    end
    if list_size == data_set
        % then the current value will be outside the bounds of the list
        set(S.lb_DS,'value',data_set-1);
    end
    
    % creating the new list (without the deleted item),
    % for the list of names and the list of paths and the list of locations
    
    helper_vector = 1:list_size;
    
    data_list_str = data_list_str(helper_vector ~= helper_vector(data_set));
    data_list_path = data_list_path(helper_vector ~= helper_vector(data_set));
    data_list_prefs = data_list_prefs(helper_vector ~= helper_vector(data_set));
    data_list_subject = data_list_subject(helper_vector ~= helper_vector(data_set));
    data_list_obs = data_list_obs(helper_vector ~= helper_vector(data_set));
    data_list_quantity1 = data_list_quantity1(helper_vector ~= helper_vector(data_set));
    data_list_quantity2 = data_list_quantity2(helper_vector ~= helper_vector(data_set));
    data_list_maxquantity1 = data_list_maxquantity1(helper_vector ~= helper_vector(data_set));
    data_list_maxquantity2 = data_list_maxquantity2(helper_vector ~= helper_vector(data_set));

    % updating the list
    set(S.lb_DS,'string', data_list_str);
end



%% Program for browse for file button
function [] = browse_file_call(varargin)
    % Callback for pushbutton.
    
    % browse for a data file
    [FileName, PathName, ~] = uigetfile('../*.csv', 'Select a CSV Data File', strcat(main_folder, '/', HPZ_Constants.data_files_dir));
    
    % if the user chose a file (didn't cancel)
    if (FileName ~= 0)
        
        [data_name, locations, pref_class, ok] = HPZ_Screen_Data_Set_Adding(FileName);
        if (ok == 1)
            % then the list will not be empty after the addition - 
            % removing file and pressing ok should be disabled
            set(S.remove_file,'enable', 'on');
            set(S.ok_button,'enable', 'on');

            % current list
            %current_list = get(S.lb_DS,'string');
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
            
            data_list_subject{list_size+1} = locations(1);
            data_list_obs{list_size+1} = locations(2);
            data_list_quantity1{list_size+1} = locations(3);
            data_list_quantity2{list_size+1} = locations(4);
            data_list_maxquantity1{list_size+1} = locations(5);
            data_list_maxquantity2{list_size+1} = locations(6);
            
            % updating the list
            set(S.lb_DS,'string', data_list_str);
            
            % automatically have the new data file be the selected one
            set(S.lb_DS,'value',list_size+1);
        end
    end
end





%% Program for OK button
function [] = ok_button_call(varargin)
    % Callback for OK pushbutton.

    % Setting the dataset selection
    data_set = get(S.lb_DS, 'value');

    % Setting whether to fix endwoments
    switch findobj(get(S.bg_ff,'selectedobject'))
        case S.fix_endowments_rd(1)
            % No
            fix_endowments = false;

        case S.fix_endowments_rd(2)
            % Yes
            fix_endowments = true;
    end

    
    file_path = data_list_path{data_set};
    locations = [data_list_subject{data_set},...
                data_list_obs{data_set},...
                data_list_quantity1{data_set},...
                data_list_quantity2{data_set},...
                data_list_maxquantity1{data_set},...
                data_list_maxquantity2{data_set}];
    
    [data_matrix, success] = HPZ_Data_Format (file_path, locations);
    if (~success)
        % the file does not exist or is not available - we end this run
        msgbox(char(strcat(HPZ_Constants.could_not_read_file_1, {' '}, file_path, HPZ_Constants.could_not_read_file_2)));
        % remove the file from the list
        remove_file_call();
    else
        % check if endowments are (approximately) equal to 1, otherwise print a warning 
        % the endowments (may not equal to 1)
        endowments = data_matrix(:,3) .* data_matrix(:,5) + data_matrix(:,4) .* data_matrix(:,6);
        % number of endowment significantly different from 1
        num_of_errors = sum ( (endowments > 1+HPZ_Constants.fix_endowments_error) | (endowments < 1/(1+HPZ_Constants.fix_endowments_error)) );
        if (num_of_errors)
            warning('There were %d cases that the endowment was significantly (more than %0.5g) different from 1. If your data endowments are not meant to be set to 1, the estimation may give erroneous results. If you chose the "fix endowments" option, the endowments were fixed to 1 in some manner, but you may still want to check your data.', num_of_errors, HPZ_Constants.fix_endowments_error);
        end

        % fixing the bundles so their endowments will be equal exactly to 1.
        % this fix is needed because of rounding that makes the endowment a little
        % different from 1, which may damage the estimations to some extent
        if (fix_endowments == true)
            % if the user chose to fix all bundles so their endwoments will be exactly 1, 
            % if the difference (in percentage) between 1 and the endowment is more
            % than this threshold, a warning will be given
            data_matrix = HPZ_Fix_Endowments_To_One(data_matrix, 1);
        end

        data_pref_class = data_list_prefs{data_set};
        ok = 1;
        HPZ_Data_Settings_Write(data_list_str, data_list_path, data_list_prefs, data_list_subject, data_list_obs, data_list_quantity1, data_list_quantity2, data_list_maxquantity1, data_list_maxquantity2, data_set, main_folder);
        % close the window
        close(S.fh);
    end

end



%% Program for Cancel button
function [] = cancel_button_call(varargin)
    % Callback for Cancel pushbutton.

    ok = 0;

    % close the window
    close(S.fh);

end



end

