function [residual_flag, in_sample_flag, out_sample_flag, ok] = HPZ_Screen_Residual_Calculation_Settings(residual_flag, in_sample_flag, out_sample_flag, runs_counter)

% this function promotes a user-interface screen to the user, allowing her
% to choose whether to calculate residuals, and if so - whether to
% calculate in-sample or out-of-sample residuals (or both).

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
figure_width = 420;
figure_height = 250;
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

% relative location of elements:
% how much distance between parts (different indexes) that are one below each other 
move_down = 15;

% distance of highest element from top
top_dist = 10;

% each radio options will be with height 25
radio_height = 25;
% distance of radio option from the bottom of the radio frame
radio_bottom = 10;
% offset to the right of radio options that are vertical
radio_offset = 15;
% offsets to the right of yes/no horizontal radio options
yes_no_offsets = [25 , 180];

% width of each element
element_width = 360;
% normal_width for a sub element
sub_element_width = 80;

% font size of label
%label_font_size = 12;
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
figure_title = strcat('Residuals Calculation Settings (', HPZ_Constants.current_run_screen, num2str(runs_counter), ')');
[S.fh , S.panel] = ui_scroll_screen(figure_width, figure_height, scroll_width, max_height_percent, top_space_height, bottom_space_height, figure_title);

% width including scroll bar if there is one
pos = get(S.fh,'position');
full_width = pos(3);





%% Calculate residuals or not
% current height of element
current_height = 60;
% current bottom coordinate
current_bottom = panel_height-top_dist - current_height;
S.residuals_rd_group = uibuttongroup('Parent',S.panel, ...
    'units','pix',...
    'title', 'Residuals', ...
    'fontsize',font_size,...
    'pos',[(figure_width-element_width)/2 , current_bottom , element_width , current_height]);

S.residuals_rd(1) = uicontrol(S.residuals_rd_group,...
    'value',residual_flag, ...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_height],...
    'fontsize',font_size,...
    'string',' YES');
S.residuals_rd(2) = uicontrol(S.residuals_rd_group,...
    'value',1-residual_flag, ...
    'enable','on', ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_height],...
    'fontsize',font_size,...
    'string',' NO');



%% Which residuals to calculate (in-sample, out-of-sample, or both)
% current height of element
current_height = 90;
% current bottom coordinate
current_bottom = current_bottom - move_down - current_height;
S.in_out_sample_rd_group = uibuttongroup('Parent',S.panel, ...
    'units','pix',...
    'title', 'Options', ...
    'fontsize',font_size,...
    'pos',[(figure_width-element_width)/2 , current_bottom , element_width , current_height]);
S.in_out_sample_rd(1) = uicontrol(S.in_out_sample_rd_group, ...
    'value',(in_sample_flag && ~out_sample_flag),...
    'enable',enable{residual_flag+1},...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_height*2 , element_width , radio_height],...
    'string',' In Sample Calculation',...
    'fontsize',big_font_size);
S.in_out_sample_rd(2) = uicontrol(S.in_out_sample_rd_group, ...
    'value',(~in_sample_flag && out_sample_flag),...
    'enable',enable{residual_flag+1},...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_height*1 , element_width , radio_height],...
    'string',' Out of Sample Calculation',...
    'fontsize',big_font_size);
S.in_out_sample_rd(3) = uicontrol(S.in_out_sample_rd_group, ...
    'value',(in_sample_flag && out_sample_flag),...
    'enable',enable{residual_flag+1},...
    'style','rad',...
    'unit','pix',...         
    'position',[radio_offset , radio_height*0 , element_width , radio_height],...
    'string',' Both',...
    'fontsize',big_font_size);

set(S.residuals_rd(:),'callback',{@residuals_rd_call,S})  % Set callback.

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





%% CRRA and Corners and parameter settings
function [] = residuals_rd_call(varargin)
    % Callback for pushbutton.
    %S = varargin{3};  % Get the structure.

    if (get(S.residuals_rd(2),'value') == 1.0)
        % Not Interested in Residuals 
        set(S.in_out_sample_rd(:), 'enable', 'off');
    else 
        set(S.in_out_sample_rd(:), 'enable', 'on');
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

    %% Residuals
    switch findobj(get(S.residuals_rd_group,'selectedobject'))
        case S.residuals_rd(1)
            residual_flag = true;
        case S.residuals_rd(2)
            residual_flag = false;
    end

    %% in-sample and/or out-of-sampe
    if get(S.in_out_sample_rd(1), 'value') == 1.0
        in_sample_flag = true;
        out_sample_flag = false;
    elseif get(S.in_out_sample_rd(2), 'value') == 1.0
        in_sample_flag = false;
        out_sample_flag = true;
    elseif get(S.in_out_sample_rd(3), 'value') == 1.0
        in_sample_flag = true;
        out_sample_flag = true;
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