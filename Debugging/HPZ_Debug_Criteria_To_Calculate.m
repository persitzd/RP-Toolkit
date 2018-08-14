function [output_file_config, ok] = HPZ_Debug_Criteria_To_Calculate()

% this function is based on "HPZ_Screen_Output_File_Format" with many
% deletions, and some minor modifications
%
% it is meant to be used in the debugging modul, to allow the user to
% calculate criteria without performing estimaation



file_val_str = {'Criterion', 'Criterion', 'Criterion', 'Criterion', 'Criterion', 'Criterion', 'Criterion'};
output_file_config = [1 , 1 , 1 , 1 , 1 , 1 , 1];
file_setting_op_enable = {'on','on','on','on','on','on','on'}; 



         

% screen size
figure_width = 450;
figure_height = 425;
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
%move_down = 10;

% distance of highest element from top
top_dist = 10;

% height of label
label_height = 30;

% distance of other elements from the left
left_other = 10;

% height of each option
option_height = 25;
% number of options
options_num = 7;
% height of the whole options frame
option_frame_height = option_height * options_num + 50;

% height of yes_no_radio_group
%radio_group_height = 50;
% each radio options will be with height 25
%radio_height = 20;
% distance of radio option from the bottom of the radio frame
radio_bottom = 10;
% offset to the right of radio options that are vertical
radio_offset = 30;
% offsets to the right of yes/no horizontal radio options
%yes_no_offsets = [25 , 180];

% width of each element
element_width = 400;
% normal_width for a sub element
%sub_element_width = 80;

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
figure_title = 'Which Criteria to Calculate?';
[S.fh , S.panel] = ui_scroll_screen(figure_width, figure_height, scroll_width, max_height_percent, top_space_height, bottom_space_height, figure_title);

% width including scroll bar if there is one
pos = get(S.fh,'position');
full_width = pos(3);


                


%% Options
% current bottom coordinate
current_bottom = panel_height - top_dist - option_frame_height;
S.file_setting_group = uibuttongroup('Parent',S.panel, ...
    'units','pix',...
    'title', 'Options', ...
    'pos',[(figure_width-element_width)/2 , current_bottom , element_width , option_frame_height]);

S.label_file_setting(1) = uicontrol(S.file_setting_group,...
    'style','text',...
    'units','pix',...
    'position',[left_other , radio_bottom+7*option_height , element_width , label_height],...
    'fontsize',label_font_size,'fontweight','bold',...
    'string','Select the values needed for the output file');
                 
S.file_setting_op(1) = uicontrol(S.file_setting_group,...
    'style','check',...
    'unit','pix',...
    'enable', file_setting_op_enable{1}, ...
    'value', output_file_config(1), ...
    'position',[radio_offset , radio_bottom+6*option_height , element_width , option_height],...
    'fontsize',font_size,...
    'string',strcat(' NLLS ',{' '}, file_val_str{1},' based on Euclidean Metric'));

S.file_setting_op(2) = uicontrol(S.file_setting_group,...
    'style','check',...
    'enable', file_setting_op_enable{2}, ...
    'value', output_file_config(2), ...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+5*option_height , element_width , option_height],...
    'fontsize',font_size,...
    'string',strcat(' NLLS ',{' '},  file_val_str{2},' based on CFGK Metric'));

S.file_setting_op(3) = uicontrol(S.file_setting_group,...
    'style','check',...
    'unit','pix',...
    'enable', file_setting_op_enable{3}, ...
    'value', output_file_config(3), ...
    'position',[radio_offset , radio_bottom+4*option_height , element_width , option_height],...
    'fontsize',font_size,...
    'string',strcat(' NLLS ',{' '}, file_val_str{3},' based on normalized-Euclidean Metric'));

S.file_setting_op(4) = uicontrol(S.file_setting_group,...
    'style','check',...
    'unit','pix',...
    'enable', file_setting_op_enable{4}, ...
    'value', output_file_config(4), ...
    'position',[radio_offset , radio_bottom+3*option_height , element_width , option_height],...
    'fontsize',font_size,...
    'string',strcat(' MMI ',{' '}, file_val_str{4},' based on MAX Waste'));

S.file_setting_op(5) = uicontrol(S.file_setting_group,...
    'style','check',...
    'unit','pix',...
    'enable', file_setting_op_enable{5}, ...
    'value', output_file_config(5), ...
    'position',[radio_offset , radio_bottom+2*option_height , element_width , option_height],...
    'fontsize',font_size,...
    'string',strcat(' MMI ',{' '}, file_val_str{5},' based on MEAN Waste'));

S.file_setting_op(6) = uicontrol(S.file_setting_group,...
    'style','check',...
    'unit','pix',...
    'enable', file_setting_op_enable{6}, ...
    'value', output_file_config(6), ...
    'position',[radio_offset , radio_bottom+1*option_height , element_width , option_height],...
    'fontsize',font_size,...
    'string', strcat(' MMI ',{' '}, file_val_str{6},' based on AVG(SSQ(Wastes))'));

S.file_setting_op(7) = uicontrol(S.file_setting_group,...
    'style','check',...
    'unit','pix',...
    'enable', file_setting_op_enable{7}, ...
    'value', output_file_config(7), ...
    'position',[radio_offset , radio_bottom+0*option_height , element_width , option_height],...
    'fontsize',font_size,...
    'string',strcat(' BI ',{' '}, file_val_str{7}));
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





%% OK button
function [] = ok_button_call(varargin)
    
    ok = 1;   
    
    output_flag_temp = cell2mat(get(S.file_setting_op(:), 'value'));
    output_file_config = output_flag_temp';

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