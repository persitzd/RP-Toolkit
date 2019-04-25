function [numeric_flag, function_flag, param1_restrictions, param2_restrictions, fix_corners, ok] = HPZ_Screen_Functional_Form_Settings_Risk(action_flag, numeric_flag, function_flag, param1_restrictions, param2_restrictions, fix_corners, runs_counter)

% this function promotes a user-interface screen to the user, allowing her
% to make choices regarding the functional form and the solution approach
% (numeric or analytic), and also some other choices and restrictions.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

% it is recommended to use the numeric variables in the beginning of the 
% function when making changes to the screen.
% it is also recommended to use numeric variables in the same manner when
% adding new elements to the screen.



% this little cell helps us to convert 0 and 1 to 'off' and 'on', respectively 
enable = {'off','on'};



% screen size
figure_width = 500;
figure_height = 420;
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

% distance of highest element from top
top_dist = 10;

% distance of the labels from the left
%left_label = 15;
% distance of other elements from the left
left_other = 45;

% each radio options will be with height 25
radio_height = 25;
% distance of bottom radio option from the bottom
radio_bottom = 10;
% offset to the right of radio options that are vertical
radio_offset = 15;
% offsets to the right of yes/no horizontal radio options
yes_no_offsets = [25 , 180];
% offsets to the right of one/two/three horizontal radio options
three_offsets = [10 , 130, 250];

% width of each element
element_width = 410;
% normal_width for a sub element
sub_element_width = 80;

% font size of label
label_font_size = 12;
% general font size
font_size = 8;
% general but bigger font size
big_font_size = 10;



%% create the figure, with a slider if needed
% scroll bar width (if scroll is needed)
scroll_width = 20;
bottom_space_height = buttons_space_height;
top_space_height = 0;
panel_height = figure_height - bottom_space_height - top_space_height;
figure_title = char(strcat('Functional Form Settings (', HPZ_Constants.current_run_screen, {' '}, num2str(runs_counter), ')'));
[S.fh , S.panel] = ui_scroll_screen(figure_width, figure_height, scroll_width, max_height_percent, top_space_height, bottom_space_height, figure_title);

% width including scroll bar if there is one
pos = get(S.fh,'position');
full_width = pos(3);





%% Disappointment Aversion
% current height of element
current_height = 30;
% current bottom coordinate
current_bottom = panel_height-top_dist - current_height;
S.label_DA = uicontrol('Parent',S.panel, 'style','text',...
    'units','pix',...
    'position',[0 , current_bottom , figure_width , current_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',label_font_size, 'fontweight','bold',...
    'string','Disappointment Aversion');

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
    'value',(function_flag == HPZ_Constants.CRRA_func)*1,...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'fontsize',font_size,...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_height],...
    'string',' CRRA');

S.functional_form_rd(2) = uicontrol(S.bg_ff,...
    'value',(function_flag == HPZ_Constants.CARA_func)*1,...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'fontsize',font_size,...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_height],...
    'string',' CARA');

%% Numeric optimization or Analytical solution
% current height of element
current_height = 60;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.bg_na = uibuttongroup('Parent',S.panel, ...
    'units','pix',...
    'title', 'Solution Options', ...
    'fontsize',font_size,...
    'pos',[left_other , current_bottom , element_width , current_height]);

S.solution_option_rd(1) = uicontrol(S.bg_na,...
    'value', (numeric_flag == HPZ_Constants.numeric)*1,...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'fontsize',font_size,...
    'position',[three_offsets(1) , radio_bottom , element_width , radio_height],...
    'string',' Numerical Approach');

S.solution_option_rd(2) = uicontrol(S.bg_na,...
    'value', (numeric_flag == HPZ_Constants.analytic || (numeric_flag == HPZ_Constants.semi_numeric && ~(action_flag == HPZ_Constants.MMI_action || action_flag == HPZ_Constants.BI_action)))*1,...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'fontsize',font_size,...
    'position',[three_offsets(2) , radio_bottom , element_width , radio_height],...
    'string',' Analytical Approach');

S.solution_option_rd(3) = uicontrol(S.bg_na,...
    'value', (numeric_flag == HPZ_Constants.semi_numeric && (action_flag == HPZ_Constants.MMI_action || action_flag == HPZ_Constants.BI_action))*1,...
    'enable', enable{1 + 1*(action_flag == HPZ_Constants.MMI_action || action_flag == HPZ_Constants.BI_action)},...
    'style','rad',...
    'unit','pix',...
    'fontsize',font_size,...
    'position',[three_offsets(3) , radio_bottom , element_width , radio_height],...
    'string',' Semi-Numerical Approach');

%% Parameter Settings (Zeros, and negative parameter setting)
% current height of element
current_height = 90;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.beta_choice_group = uibuttongroup('Parent',S.panel, ...
    'units','pix',...
    'title', 'Parameter Setting', ...
    'fontsize',font_size,...
    'pos',[left_other , current_bottom , element_width , current_height]);

negative_beta_value = ((param1_restrictions(1) == -1) && param1_restrictions(2) == inf);
param_zero_value = ((param1_restrictions(1) == 0) && param1_restrictions(2) == 0);

S.negative_beta_ch = uicontrol(S.beta_choice_group, ...
    'value', negative_beta_value,... 
    'enable', enable{2 - param_zero_value},... 
    'style','check',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+radio_height , element_width , radio_height],...
    'string',' Allow DA coefficient to be negative (bounded to -1)',...
    'fontsize',big_font_size);

S.param_zero_chk = uicontrol(S.beta_choice_group,...
    'value', param_zero_value,... 
    'enable', enable{2 - negative_beta_value},... 
    'style','check',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom , element_width , radio_height],...
    'string',' DA Coefficient = 0',...
    'fontsize',big_font_size);

%% Corners
% current height of element
current_height = 90;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.corners_rd_group = uibuttongroup('Parent',S.panel, ...
    'units','pix',...
    'title', 'Corners', ...
    'fontsize',font_size,...
    'pos',[left_other , current_bottom , element_width , current_height]);
S.corners_rd(1) = uicontrol(S.corners_rd_group, ...
    'value',1-fix_corners,...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+radio_height , element_width , radio_height],...
    'string',' Allow boundary choices',...
    'fontsize',big_font_size);
S.corners_rd(2) = uicontrol(S.corners_rd_group, ...
    'value',fix_corners,...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom , element_width , radio_height],...
    'string',' Adjust boundary choices [Choi et al. (2007)]',...
    'fontsize',big_font_size);

% this function makes sure that if CARA is chosen, the fix corner option is
% disabled
rd_call();

set(S.negative_beta_ch,'callback',{@ch_call_1,S})  % Set callback.
set(S.param_zero_chk,'callback',{@ch_call_2,S})  % Set callback.
set(S.functional_form_rd(:),'callback',{@rd_call,S})  % Set callback.

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





%% Negative Beta and Beta = 0
function [] = ch_call_1(varargin)
    % Callback for pushbutton.
    %S = varargin{3};  % Get the structure.

    if get(S.negative_beta_ch,'value') == 1.0
        set(S.param_zero_chk, 'enable', 'off');
        set(S.param_zero_chk, 'value', 0.0);              
    else
        set(S.param_zero_chk, 'enable', 'on');
    end

end

function [] = ch_call_2(varargin)
    % Callback for pushbutton.
    %S = varargin{3};  % Get the structure.

    if get(S.param_zero_chk,'value') == 1.0
        set(S.negative_beta_ch, 'enable', 'off');
        set(S.negative_beta_ch, 'value', 0.0);            
    else
        set(S.negative_beta_ch, 'enable', 'on');
    end

end

%% CRRA and Corners and parameter settings
function [] = rd_call(varargin)
    % Callback for pushbutton.
    %S = varargin{3};  % Get the structure.

    if (get(S.functional_form_rd(2),'value') == 1.0)
        % CARA 
        % Disable the corner correction
        set(S.corners_rd(:), 'enable', 'off');
        set(S.corners_rd(1), 'value', 1);
        set(S.corners_rd(2), 'value', 0);
        %set(S.negative_beta_ch, 'enable', 'on');
        %set(S.param_zero_chk, 'enable', 'on');
        %set(S.solution_option_rd(1), 'enable', 'on');
        %set(S.solution_option_rd(1), 'value', 0);
        %set(S.solution_option_rd(2), 'enable', 'on');
        %set(S.solution_option_rd(2), 'value', 1);

    else % in case of CRRA
        set(S.corners_rd(:), 'enable', 'on');
        %set(S.negative_beta_ch, 'enable', 'on');
        %set(S.param_zero_chk, 'enable', 'on');
        %set(S.solution_option_rd(:), 'enable', 'on');
        %set(S.solution_option_rd(2), 'value', 1);
    end
end





%% OK button
function [] = ok_button_call(varargin)
    
    ok = 1;
    
    % Callback for pushbutton.
    %S = varargin{3};  % Get the structure.
    % Instead of switch, we could use num2str on:
    % find(get(S.bg,'selectedobject')==S.rd)      (or similar)
    % Note the use of findobj.  This is because of a BUG in MATLAB, whereby if
    % the user selects the same button twice, the selectedobject property will
    % not work correctly.

    %% Setting Functional Form
    switch findobj(get(S.bg_ff,'selectedobject'))
        case S.functional_form_rd(1)
            % CRRA
            function_flag = HPZ_Constants.CRRA_func;

        case S.functional_form_rd(2)
            % CARA
            function_flag = HPZ_Constants.CARA_func;
            % Zeros_Flag
            fix_corners = false;   % corners are allowed, no need for fix
    end

    %% Solution Options
    if get(S.solution_option_rd(1), 'value') == 1.0
        numeric_flag = HPZ_Constants.numeric;
    elseif get(S.solution_option_rd(2), 'value') == 1.0
        numeric_flag = HPZ_Constants.analytic;
    elseif get(S.solution_option_rd(3), 'value') == 1.0
        numeric_flag = HPZ_Constants.semi_numeric;
    end

    %% Setting parameters for the seleceted functional form
    if get(S.negative_beta_ch, 'value') == 1.0
        % Beta in [-1, Inf)
        % (the widest possible range for beta)
        param1_restrictions = [-1 , inf];

    elseif get(S.negative_beta_ch, 'value') == 0.0 && get(S.param_zero_chk, 'value') == 1.0
        % Beta = 0
        param1_restrictions = [0 , 0];

    elseif get(S.negative_beta_ch, 'value') == 0.0
        % Beta in [0 , Inf)
        param1_restrictions = [0 , inf];

    end

    %% Corners
    if ~(get(S.functional_form_rd(2),'value') == 1.0)
        if get(S.corners_rd(1), 'value') == 1.0
            % corners are allowed, no need to fix corners
            fix_corners = false;
        else
            % Choi et al. (2007) for corner correction is required
            fix_corners = true;
        end
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



end
