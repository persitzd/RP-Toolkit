function [numeric_flag, function_flag, param1_restrictions, param2_restrictions, ok] = HPZ_Screen_Functional_Form_Settings_OR(numeric_flag, function_flag, param1_restrictions, param2_restrictions, runs_counter)

% this function promotes a user-interface screen to the user, allowing her
% to make choices regarding the functional form and the solution approach
% (numeric or analytic)

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

% it is recommended to use the numeric variables in the beginning of the 
% function when making changes to the screen.
% it is also recommended to use numeric variables in the same manner when
% adding new elements to the screen.



% this little cell helps us to convert 0 and 1 to 'off' and 'on', respectively 
enable = {'off','on'}; %#ok<NASGU>



% screen size
figure_width = 450;
figure_height = 270;
% limit of height as percentage of computer screen height
max_height_percent = HPZ_Constants.max_height_percent;

% height of bottom space designated for OK and Cancel buttons
buttons_space_height = 50;
% buttons height as percentage of the space they are in
buttons_height = 0.6;
% the distant between a button to the edge of the screen,
% and is also the distance between buttons.
% the size of the buttons is designed to fit this and the buttons_height 
% and the width of the screen
buttons_dists = 20;
% number of buttons
buttons_num = 2;

% distance of highest element from top
top_dist = 10;

% distance of the labels from the left
% left_label = 15;
% distance of other elements from the left
left_other = 45;

% each radio options will be with height 25
radio_height = 25;
% distance of bottom radio option from the bottom
radio_bottom = 10;
% offset to the right of radio options that are vertical
% radio_offset = 15;
% offsets to the right of yes/no horizontal radio options
yes_no_offsets = [25 , 180];

% width of each element
element_width = 360;
% normal_width for a sub element
sub_element_width = 80;

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
figure_title = strcat('Functional Form Settings (', HPZ_Constants.current_run_screen, num2str(runs_counter), ')');
[S.fh , S.panel] = ui_scroll_screen(figure_width, figure_height, scroll_width, max_height_percent, top_space_height, bottom_space_height, figure_title);

% width including scroll bar if there is one
pos = get(S.fh,'position');
full_width = pos(3);





%% Other Regarding
% current height of element
current_height = 30;
% current bottom coordinate
current_bottom = panel_height-top_dist - current_height;
S.label_OR = uicontrol('Parent',S.panel, 'style','text',...
    'units','pix',...
    'position',[0 , current_bottom , figure_width , current_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',label_font_size, 'fontweight','bold',...
    'string','Other Regarding');



%% Functional Forms
% current height of element
current_height = 60;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.bg_ff = uibuttongroup('Parent',S.panel, ...
    'units','pix',...
    'title', 'Functional Forms', ...
    'fontsize',font_size,...
    'pos',[left_other , current_bottom , element_width , current_height]);

S.functional_form_rd(1) = uicontrol(S.bg_ff,...
    'value',(function_flag == HPZ_Constants.CES_func),...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'fontsize',font_size,...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_height],...
    'string',' CES');

%% Numeric optimization or Analytical solution
% current height of element
current_height = 60;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.bg_na = uibuttongroup('Parent',S.panel, ...
    'units','pix',...
    'title','Solution Options', ...
    'fontsize',font_size,...
    'pos',[left_other , current_bottom , element_width , current_height]);

S.solution_option_rd(1) = uicontrol(S.bg_na,...
    'value', numeric_flag,...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'fontsize',font_size,...
    'position',[yes_no_offsets(1) , radio_bottom , element_width , radio_height],...
    'string',' Numerical Approach');

S.solution_option_rd(2) = uicontrol(S.bg_na,...
    'value', 1-numeric_flag,...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'fontsize',font_size,...
    'position',[yes_no_offsets(2) , radio_bottom , element_width , radio_height],...
    'string',' Analytical Approach');

% set(S.functional_form_rd(:),'callback',{@rd_call,S})  % Set callback.

%% end of settings





%% OK Button
S.ok_button = uicontrol('Parent',S.fh, 'style','push',...
    'unit','pix',...
    'position',[buttons_dists , (1-buttons_height)/2*buttons_space_height , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_height*buttons_space_height],...
    'string','OK',...
    'fontsize',font_size,...
    'callback',{@ok_button_call,S});

%% Cancel Button
S.cancel_button = uicontrol('Parent',S.fh, 'style','push',...
    'unit','pix',...
    'position',[(full_width/2)+buttons_dists , (1-buttons_height)/2*buttons_space_height , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_height*buttons_space_height],...
    'string','Cancel',...
    'fontsize',font_size,...
    'callback',{@cancel_button_call,S});





uiwait(S.fh)  % Prevent all other processes from starting until closed.





%% CES settings
% function [] = rd_call(varargin)
%     % Callback for pushbutton.
%     S = varargin{3};  % Get the structure.
% 
%     if (get(S.functional_form_rd(1),'value') == 1)
%         set(S.solution_option_rd(:), 'enable', 'on');
%         set(S.solution_option_rd(2), 'value', 1);
%     end
% end





%% OK button
function [] = ok_button_call(varargin)
    
    ok = 1;
    
    % Callback for pushbutton.
    S = varargin{3};  % Get the structure.
    % Instead of switch, we could use num2str on:
    % find(get(S.bg,'selectedobject')==S.rd)      (or similar)
    % Note the use of findobj.  This is because of a BUG in MATLAB, whereby if
    % the user selects the same button twice, the selectedobject property will
    % not work correctly.

    %% Setting Functional Form
    switch findobj(get(S.bg_ff,'selectedobject'))
        case S.functional_form_rd(1)
            % CES
            function_flag = 1;
    end

    %% Solution Options
    if get(S.solution_option_rd(1), 'value') == 1.0
        numeric_flag = true;
    else
        numeric_flag = false;
    end

    % close the window
    close(S.fh);
end



%% Cancel button
function [] = cancel_button_call(varargin)
    % Callback for Cancel button.

    ok = 0;

    % close the window
    close(S.fh);

end



end   % end of main function