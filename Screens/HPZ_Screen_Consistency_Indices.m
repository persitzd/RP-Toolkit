function [GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, ok] = HPZ_Screen_Consistency_Indices(GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, runs_counter)  

% this function creates a user interface screen from which the user chooses
% which of the consistency and inconsistency measures to be calculated, and
% also if to calculate the residuals for each observation, and also if to
% calculate the residuals using an "in sample" approach or an "out of sample"
% approach.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% this little cell helps us to convert 0 and 1 to 'off' and 'on', respectively 
enable = {'off','on'};



%% location and size parameters (use these to easily make changes to the screen)

% screen size
figure_width = 780;
figure_height = 690;
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
% how much the right part should be offset to the right in comparison to the left part 
move_right = element_width + 60;
% how much distance between parts (different indices) that are one below each other 
move_down = 20;

% distance of highest element from top
top_dist = 10;

% distance of the labels from the left
left_label = 15;
% distance of other elements from the left
left_other = 45;

% sizes of labels used in the head of each section
label_width = 300;
label_height = 30;

% each radio options will be with height 25
% an element with X radio options one below another will be of height 25X + 15  
% a "yes-no" element (except for the special ones in VARIAN and HM) with X radio options one right to another may be of height 35 + 15 
radio_height = 20;
radio_label = 25;
radio_yes_no = 35;
% distance of bottom radio option from the bottom
radio_bottom = 5;
% offset to the right of radio options that are vertical
radio_offset = 15;
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
bottom_space_height = buttons_space_height;
top_space_height = 0;
panel_height = figure_height - bottom_space_height - top_space_height;
figure_title = char(strcat('Inconsistency Indices (', HPZ_Constants.current_run_screen, {' '}, num2str(runs_counter), ')'));
[S.fh , S.panel] = ui_scroll_screen(figure_width, figure_height, scroll_width, max_height_percent, top_space_height, bottom_space_height, figure_title);

% width including scroll bar if there is one
pos = get(S.fh,'position');
full_width = pos(3);



%% GARP WARP and SARP settings
% current bottom coordinate
current_bottom = panel_height-top_dist-label_height;
S.label_GARP = uicontrol('Parent',S.panel, ...
    'style','text',...
    'units','pix',...
    'position',[left_label , current_bottom , label_width , label_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',label_font_size,'fontweight','bold',...
    'HorizontalAlignment','left',...
    'string','GARP, WARP and SARP violations');

%% calculate GARP WARP and SARP or not
% current bottom coordinate
current_bottom = current_bottom - (radio_yes_no+radio_label);
S.calculate_GARP = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Calculate', ...
    'fontsize',font_size,...
    'pos',[left_other , current_bottom , element_width , radio_yes_no+radio_label]);

S.calculate_GARP_choice(1) = uicontrol(S.calculate_GARP,...
    'value',GARP_flags(1), ...
    'enable','on', ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' YES');
S.calculate_GARP_choice(2) = uicontrol(S.calculate_GARP,...
    'value',1-GARP_flags(1), ...
    'enable','on', ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' NO');

%% With or Without residuals
enable_temp = enable{GARP_flags(1)+1};
% current height of element
current_height = radio_yes_no+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.residual_GARP = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Residuals', ...
    'fontsize',font_size,...
    'pos',[left_other , current_bottom , element_width , current_height]);

S.residual_GARP_choice(1) = uicontrol(S.residual_GARP,...
    'value',GARP_flags(2), ...
    'enable',enable_temp, ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' YES');
S.residual_GARP_choice(2) = uicontrol(S.residual_GARP,...
    'value',1-GARP_flags(2), ...
    'enable',enable_temp, ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' NO');

%% In Sample and/or Out Of Sample
enable_temp = enable{(GARP_flags(1)&&GARP_flags(2))+1};
% current height of element
current_height = 3*radio_height+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.residual_GARP_options_group = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Options', ...
    'fontsize',font_size,...
    'pos',[left_other , current_bottom , element_width , current_height]);
S.residual_GARP_options(1) = uicontrol(S.residual_GARP_options_group, ...
    'value',(GARP_flags(3) && ~GARP_flags(4)),...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+2*radio_height , element_width , radio_height],...
    'string',' In Sample Calculation',...
    'fontsize',big_font_size);
S.residual_GARP_options(2) = uicontrol(S.residual_GARP_options_group, ...
    'value',(~GARP_flags(3) && GARP_flags(4)),...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+1*radio_height , element_width , radio_height],...
    'string',' Out of Sample Calculation',...
    'fontsize',big_font_size);
S.residual_GARP_options(3) = uicontrol(S.residual_GARP_options_group, ...
    'value',(GARP_flags(3) && GARP_flags(4)),...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+0*radio_height , element_width , radio_height],...
    'string',' Both',...
    'fontsize',big_font_size);

%% Which residuals to print - WARP / GARP / SARP / All
enable_temp = enable{(GARP_flags(1)&&GARP_flags(2))+1};
% current height of element
current_height = 4*radio_height+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.which_residual_GARP_options_group = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Which residuals to print', ...
    'fontsize',font_size,...
    'pos',[left_other , current_bottom , element_width , current_height]);
S.which_residual_GARP_options(1) = uicontrol(S.which_residual_GARP_options_group, ...
    'value',(GARP_flags(5) && ~GARP_flags(6) && ~GARP_flags(7)),...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+3*radio_height , element_width , radio_height],...
    'string',' WARP',...
    'fontsize',big_font_size);
S.which_residual_GARP_options(2) = uicontrol(S.which_residual_GARP_options_group, ...
    'value',(~GARP_flags(5) && GARP_flags(6) && ~GARP_flags(7)),...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+2*radio_height , element_width , radio_height],...
    'string',' GARP',...
    'fontsize',big_font_size);
S.which_residual_GARP_options(3) = uicontrol(S.which_residual_GARP_options_group, ...
    'value',(~GARP_flags(5) && ~GARP_flags(6) && GARP_flags(7)),...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+1*radio_height , element_width , radio_height],...
    'string',' SARP',...
    'fontsize',big_font_size);
S.which_residual_GARP_options(4) = uicontrol(S.which_residual_GARP_options_group, ...
    'value',(GARP_flags(5) && GARP_flags(6) && GARP_flags(7)),...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+0*radio_height , element_width , radio_height],...
    'string',' All',...
    'fontsize',big_font_size);


set(S.calculate_GARP_choice(:),'callback',{@GARP_calc_call,S})  % Set callback.
set(S.residual_GARP_choice(:),'callback',{@GARP_rd_call,S})     % Set callback.

%% End of GARP WARP and SARP settings





%% MPI settings
% MPI will be ready for the next version (version 3), currently not implemented 
% % current bottom coordinate (MPI is just below GARP)
% current_bottom = current_bottom - label_height - move_down;
% 
% S.label_MPI = uicontrol('Parent',S.panel, 'style','text',...
%                  'units','pix',...
%                  'position',[left_label , current_bottom , label_width , label_height],...
%                  'backgroundc',get(S.fh,'color'),...
%                  'fontsize',label_font_size,'fontweight','bold',... 
%                  'HorizontalAlignment','left',...
%                  'string','MPI inconsistency index');
% 
% %% calculate MPI inconsistency index or not
% 
% % current height of element
% current_height = radio_yes_no+radio_label;
% % current bottom coordinate
% current_bottom = current_bottom - current_height;
% S.calculate_MPI = uibuttongroup('Parent',S.panel, 'units','pix',...
%     'title', 'Calculate', ...
%     'fontsize',font_size,...
%     'pos',[left_other , current_bottom , element_width , current_height]);
% 
% S.calculate_MPI_choice(1) = uicontrol(S.calculate_MPI,...
%     'value',MPI_flags(1), ...
%     'enable','on',...
%     'style','rad',...
%     'unit','pix',...
%     'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_yes_no],...
%     'fontsize',font_size,...
%     'string',' YES');
% S.calculate_MPI_choice(2) = uicontrol(S.calculate_MPI,...
%     'value',1-MPI_flags(1), ...
%     'enable','on', ...
%     'style','rad',...
%     'unit','pix',...
%     'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_yes_no],...
%     'fontsize',font_size,...
%     'string',' NO');
% 
% %% With or Without residuals
% enable_temp = enable{MPI_flags(1)+1};
% % current height of element
% current_height = radio_yes_no+radio_label;
% % current bottom coordinate
% current_bottom = current_bottom - current_height;
% S.residual_MPI = uibuttongroup('Parent',S.panel, 'units','pix',...
%     'title', 'Residuals', ...
%     'fontsize',font_size,...
%     'pos',[left_other , current_bottom , element_width , current_height]);
% 
% S.residual_MPI_choice(1) = uicontrol(S.residual_MPI,...
%     'value',MPI_flags(2), ...
%     'enable',enable_temp,...
%     'style','rad',...
%     'unit','pix',...
%     'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_yes_no],...
%     'fontsize',font_size,...
%     'string',' YES');
% S.residual_MPI_choice(2) = uicontrol(S.residual_MPI,...
%     'value',1-MPI_flags(2), ...
%     'enable',enable_temp,...
%     'style','rad',...
%     'unit','pix',...
%     'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_yes_no],...
%     'fontsize',font_size,...
%     'string',' NO');
% 
% %% In Sample and/or Out Of Sample
% enable_temp = enable{(MPI_flags(1)&&MPI_flags(2))+1};
% % current height of element
% current_height = radio_height+radio_label;
% % current bottom coordinate
% current_bottom = current_bottom - current_height;
% S.residual_MPI_options_group = uibuttongroup('Parent',S.panel, 'units','pix',...
%     'title', 'Options', ...
%     'fontsize',font_size,...
%     'pos',[left_other , current_bottom , element_width , current_height]);
% S.residual_MPI_options(1) = uicontrol(S.residual_MPI_options_group, ...
%     'value',MPI_flags(3), ...
%     'enable',enable_temp,...
%     'style','rad',...
%     'unit','pix',...
%     'position',[radio_offset , radio_bottom , element_width , radio_height],...
%     'string',' Out of Sample Calculation',...
%     'fontsize',big_font_size);
% 
% 
% set(S.calculate_MPI_choice(:),'callback',{@MPI_calc_call,S})  % Set callback.
% set(S.residual_MPI_choice(:),'callback',{@MPI_rd_call,S})     % Set callback.

%% End of MPI settings





%% AFRIAT settings

% current bottom coordinate (AFRIAT is at the bottom, to the right of GARP)
current_bottom = panel_height-label_height-top_dist;

S.label_AFRIAT = uicontrol('Parent',S.panel, 'style','text',...
                 'units','pix',...
                 'position',[left_label+move_right , current_bottom , label_width , label_height],...
                 'backgroundc',get(S.fh,'color'),...
                 'fontsize',label_font_size,'fontweight','bold',... 
                 'HorizontalAlignment','left',...
                 'string','AFRIAT inconsistency index');

%% calculate AFRIAT inconsistency index or not
% current height of element
current_height = radio_yes_no+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.calculate_AFRIAT = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Calculate', ...
    'fontsize',font_size,...
    'pos',[left_other+move_right , current_bottom , element_width , current_height]);

S.calculate_AFRIAT_choice(1) = uicontrol(S.calculate_AFRIAT,...
    'value',AFRIAT_flags(1), ...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' YES');
S.calculate_AFRIAT_choice(2) = uicontrol(S.calculate_AFRIAT,...
    'value',1-AFRIAT_flags(1), ...
    'enable','on', ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' NO');

%% With or Without residuals
enable_temp = enable{AFRIAT_flags(1)+1};
% current height of element
current_height = radio_yes_no+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.residual_AFRIAT = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Residuals', ...
    'fontsize',font_size,...
    'pos',[left_other+move_right , current_bottom , element_width , current_height]);

S.residual_AFRIAT_choice(1) = uicontrol(S.residual_AFRIAT,...
    'value',AFRIAT_flags(2), ...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' YES');
S.residual_AFRIAT_choice(2) = uicontrol(S.residual_AFRIAT,...
    'value',1-AFRIAT_flags(2), ...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' NO');

%% In Sample and/or Out Of Sample
enable_temp = enable{(AFRIAT_flags(1)&&AFRIAT_flags(2))+1};
% current height of element
current_height = radio_height+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.residual_AFRIAT_options_group = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Options', ...
    'fontsize',font_size,...
    'pos',[left_other+move_right , current_bottom , element_width , current_height]);
S.residual_AFRIAT_options(1) = uicontrol(S.residual_AFRIAT_options_group, ...
    'value',AFRIAT_flags(3), ...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom , element_width , radio_height],...
    'string',' Out of Sample Calculation',...
    'fontsize',big_font_size);


set(S.calculate_AFRIAT_choice(:),'callback',{@AFRIAT_calc_call,S})  % Set callback.
set(S.residual_AFRIAT_choice(:),'callback',{@AFRIAT_rd_call,S})     % Set callback.

%% End of AFRIAT settings





%% VARIAN settings

% current bottom coordinate (VARIAN is at the bottom of AFRIAT)
current_bottom = current_bottom - move_down - label_height;

S.label_VARIAN = uicontrol('Parent',S.panel, 'style','text',...
                 'units','pix',...
                 'position',[left_label+move_right , current_bottom , label_width , label_height],...
                 'backgroundc',get(S.fh,'color'),...
                 'fontsize',label_font_size,'fontweight','bold',... 
                 'HorizontalAlignment','left',...
                 'string','VARIAN inconsistency index');

%% calculate VARIAN inconsistency index or not
% current height of element
current_height = radio_yes_no+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.calculate_VARIAN = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Calculate', ...
    'fontsize',font_size,...
    'pos',[left_other+move_right , current_bottom , element_width , current_height]);

S.calculate_VARIAN_choice(1) = uicontrol(S.calculate_VARIAN,...
    'value',VARIAN_flags(1), ...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' YES');
S.calculate_VARIAN_choice(2) = uicontrol(S.calculate_VARIAN,...
    'value',1-VARIAN_flags(1), ...
    'enable','on',...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' NO');

%% With or Without residuals
enable_temp = enable{VARIAN_flags(1)+1};
% current height of element
current_height = radio_yes_no+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.residual_VARIAN = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Residuals', ...
    'fontsize',font_size,...
    'pos',[left_other+move_right , current_bottom , element_width , current_height]);

S.residual_VARIAN_choice(1) = uicontrol(S.residual_VARIAN,...
    'value',VARIAN_flags(2), ...
    'enable',enable_temp, ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' YES');
S.residual_VARIAN_choice(2) = uicontrol(S.residual_VARIAN,...
    'value',1-VARIAN_flags(2), ...
    'enable',enable_temp, ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' NO');

%% In Sample and/or Out Of Sample
enable_temp = enable{(VARIAN_flags(1)&&VARIAN_flags(2))+1};
% current height of element
current_height = 3*radio_height+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.residual_VARIAN_options_group = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Options', ...
    'fontsize',font_size,...
    'pos',[left_other+move_right , current_bottom , element_width , current_height]);
S.residual_VARIAN_options(1) = uicontrol(S.residual_VARIAN_options_group, ...
    'value',(VARIAN_flags(3) && ~VARIAN_flags(4)), ...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+2*radio_height , element_width , radio_height],...
    'string',' In Sample Calculation',...
    'fontsize',big_font_size);
S.residual_VARIAN_options(2) = uicontrol(S.residual_VARIAN_options_group, ...
    'value',(~VARIAN_flags(3) && VARIAN_flags(4)), ...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'position',[radio_offset , radio_bottom+1*radio_height , element_width , radio_height],...
    'string',' Out of Sample Calculation',...
    'fontsize',big_font_size);
S.residual_VARIAN_options(3) = uicontrol(S.residual_VARIAN_options_group, ...
    'value',(VARIAN_flags(3) && VARIAN_flags(4)), ...
    'enable',enable_temp,...
    'style','rad',...
    'unit','pix',...
    'value', 1.0,...
    'position',[radio_offset , radio_bottom+0*radio_height , element_width , radio_height],...
    'string',' Both',...
    'fontsize',big_font_size);

set(S.calculate_VARIAN_choice(:),'callback',{@VARIAN_calc_call,S})  % Set callback.
set(S.residual_VARIAN_choice(:),'callback',{@VARIAN_rd_call,S})     % Set callback.

%% End of VARIAN settings





%% HOUTMAN-MAKS settings

% current bottom coordinate (HOUTMAN-MAKS is just below VARIAN)
current_bottom = current_bottom - move_down;

% current bottom coordinate
current_bottom = current_bottom - label_height;
S.label_HOUTMAN = uicontrol('Parent',S.panel, 'style','text',...
                 'units','pix',...
                 'position',[left_label+move_right , current_bottom , label_width , label_height],...
                 'backgroundc',get(S.fh,'color'),...
                 'fontsize',label_font_size,'fontweight','bold',... 
                 'HorizontalAlignment','left',...
                 'string','HOUTMAN-MAKS inconsistency index');

%% calculate HOUTMAN-MAKS inconsistency index or not
% current height of element
current_height = radio_yes_no+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.calculate_HOUTMAN = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Calculate', ...
    'fontsize',font_size,...
    'pos',[left_other+move_right , current_bottom , element_width , current_height]);

S.calculate_HOUTMAN_choice(1) = uicontrol(S.calculate_HOUTMAN,...
    'value',HOUTMAN_flags(1), ...
    'enable','on', ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' YES');
S.calculate_HOUTMAN_choice(2) = uicontrol(S.calculate_HOUTMAN,...
    'value',1-HOUTMAN_flags(1), ...
    'enable','on', ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' NO');

%% With or Without residuals
enable_temp = enable{HOUTMAN_flags(1)+1};
% current height of element
current_height = radio_yes_no+radio_label;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.residual_HOUTMAN = uibuttongroup('Parent',S.panel, 'units','pix',...
    'title', 'Residuals', ...
    'fontsize',font_size,...
    'pos',[left_other+move_right , current_bottom , element_width , current_height]);

S.residual_HOUTMAN_choice(1) = uicontrol(S.residual_HOUTMAN,...
    'value',HOUTMAN_flags(2), ...
    'enable',enable_temp, ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(1) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' YES');
S.residual_HOUTMAN_choice(2) = uicontrol(S.residual_HOUTMAN,...
    'value',1-HOUTMAN_flags(2), ...
    'enable',enable_temp, ...
    'style','rad',...
    'unit','pix',...
    'position',[yes_no_offsets(2) , radio_bottom , sub_element_width , radio_yes_no],...
    'fontsize',font_size,...
    'string',' NO');

set(S.calculate_HOUTMAN_choice(:),'callback',{@HOUTMAN_calc_call,S})    % Set callback.

%% End of HOUTMAN settings





%% OK Button
S.pb = uicontrol('Parent',S.fh, 'style','push',...
    'unit','pix',...
    'position',[buttons_dists , (1-buttons_height)/2*buttons_space_height , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_height*buttons_space_height],...
    'string','OK',...
    'fontsize',font_size,...
    'callback',{@pb_call,S});

%% Cancel Button
S.pc = uicontrol('Parent',S.fh, 'style','push',...
    'unit','pix',...
    'position',[(full_width/2)+buttons_dists , (1-buttons_height)/2*buttons_space_height , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_height*buttons_space_height],...
    'string','Cancel',...
    'fontsize',font_size,...
    'callback',{@pc_call,S});





uiwait(S.fh)   % Prevent all other processes from starting until closed.





%% managing enable and disable when choosing to calculate GARP or not
function [] = GARP_calc_call(varargin)
    % Callback for pushbutton.
    S = varargin{3};  % Get the structure.

    if (get(S.calculate_GARP_choice(2),'value') == 1.0)
        % Not interested in calcualtion
        set(S.residual_GARP_choice(:), 'enable', 'off');
        set(S.residual_GARP_options(:), 'enable', 'off');
        set(S.which_residual_GARP_options(:), 'enable', 'off');
    else 
        set(S.residual_GARP_choice(:), 'enable', 'on');
        if (get(S.residual_GARP_choice(2),'value') == 0.0)
            set(S.residual_GARP_options(:), 'enable', 'on');
            set(S.which_residual_GARP_options(:), 'enable', 'on');
        end
    end
end

%% managing enable and disable when choosing with and without residuals in GARP
function [] = GARP_rd_call(varargin)
    % Callback for pushbutton.
    S = varargin{3};  % Get the structure.

    if (get(S.residual_GARP_choice(2),'value') == 1.0)
        % Not Interested in Residuals 
        set(S.residual_GARP_options(:), 'enable', 'off');
        set(S.which_residual_GARP_options(:), 'enable', 'off');
    else 
        set(S.residual_GARP_options(:), 'enable', 'on');
        set(S.which_residual_GARP_options(:), 'enable', 'on');
    end
end



%% MPI will be ready for the next version (version 3), currently not implemented 
% %% managing enable and disable when choosing to calculate MPI or not
% function [] = MPI_calc_call(varargin)
%     % Callback for pushbutton.
%     S = varargin{3};  % Get the structure.
% 
%     if (get(S.calculate_MPI_choice(2),'value') == 1.0)
%         % Not interested in calcualtion
%         set(S.residual_MPI_choice(:), 'enable', 'off');
%         set(S.residual_MPI_options(:), 'enable', 'off');
%     else 
%         set(S.residual_MPI_choice(:), 'enable', 'on');
%         if (get(S.residual_MPI_choice(2),'value') == 0.0)
%             set(S.residual_MPI_options(:), 'enable', 'on');
%         end
%     end
% end
% 
% %% managing enable and disable when choosing with and without residuals in MPI
% function [] = MPI_rd_call(varargin)
%     % Callback for pushbutton.
%     S = varargin{3};  % Get the structure.
% 
%     if (get(S.residual_MPI_choice(2),'value') == 1.0)
%         % Not Interested in Residuals 
%         set(S.residual_MPI_options(:), 'enable', 'off');
%     else 
%         set(S.residual_MPI_options(:), 'enable', 'on');
%     end
% end



%% managing enable and disable when choosing to calculate AFRIAT or not
function [] = AFRIAT_calc_call(varargin)
    % Callback for pushbutton.
    S = varargin{3};  % Get the structure.

    if (get(S.calculate_AFRIAT_choice(2),'value') == 1.0)
        % Not interested in calcualtion
        set(S.residual_AFRIAT_choice(:), 'enable', 'off');
        set(S.residual_AFRIAT_options(:), 'enable', 'off');
    else 
        set(S.residual_AFRIAT_choice(:), 'enable', 'on');
        if (get(S.residual_AFRIAT_choice(2),'value') == 0.0)
            set(S.residual_AFRIAT_options(:), 'enable', 'on');
        end
    end
end

%% managing enable and disable when choosing with and without residuals in AFRIAT
function [] = AFRIAT_rd_call(varargin)
    % Callback for pushbutton.
    S = varargin{3};  % Get the structure.

    if (get(S.residual_AFRIAT_choice(2),'value') == 1.0)
        % Not Interested in Residuals 
        set(S.residual_AFRIAT_options(:), 'enable', 'off');
    else 
        set(S.residual_AFRIAT_options(:), 'enable', 'on');
    end
end



%% managing enable and disable when choosing to calculate VARIAN or not
function [] = VARIAN_calc_call(varargin)
    % Callback for pushbutton.
    S = varargin{3};  % Get the structure.

    if (get(S.calculate_VARIAN_choice(2),'value') == 1.0)
        % Not interested in calcualtion
        set(S.residual_VARIAN_choice(:), 'enable', 'off');
        set(S.residual_VARIAN_options(:), 'enable', 'off');
    else 
        set(S.residual_VARIAN_choice(:), 'enable', 'on');
        if (get(S.residual_VARIAN_choice(2),'value') == 0.0)
            set(S.residual_VARIAN_options(:), 'enable', 'on');
        end
    end
end

%% managing enable and disable when choosing with and without residuals in VARIAN
function [] = VARIAN_rd_call(varargin)
    % Callback for pushbutton.
    S = varargin{3};  % Get the structure.

    if (get(S.residual_VARIAN_choice(2),'value') == 1.0)
        % Not Interested in Residuals 
        set(S.residual_VARIAN_options(:), 'enable', 'off');
    else 
        set(S.residual_VARIAN_options(:), 'enable', 'on');
    end
end



%% managing enable and disable when choosing to calculate HOUTMAN-MAKS or not
function [] = HOUTMAN_calc_call(varargin)
    % Callback for pushbutton.
    S = varargin{3};  % Get the structure.

    if (get(S.calculate_HOUTMAN_choice(2),'value') == 1.0)
        % Not interested in calcualtion
        set(S.residual_HOUTMAN_choice(:), 'enable', 'off');
    else 
        set(S.residual_HOUTMAN_choice(:), 'enable', 'on');
    end
end





%% OK button
function [] = pb_call(varargin)
    % Callback for ok button.

    ok = 1;

    S = varargin{3};  % Get the structure.
    % Instead of switch, we could use num2str on:
    % find(get(S.bg,'selectedobject')==S.rd)      (or similar)
    % Note the use of findobj.  This is because of a BUG in MATLAB, whereby if
    % the user selects the same button twice, the selectedobject property will
    % not work correctly.


    %% choices regarding GARP WARP and SARP

    GARP_calculate_flag = 0;
    GARP_residual_flag = 0;
    GARP_in_sample_flag = 0;
    GARP_out_sample_flag = 0;
    WARP_res = 0;
    GARP_res = 0;
    SARP_res = 0;

    switch findobj(get(S.calculate_GARP,'selectedobject'))
        case S.calculate_GARP_choice(1)
            % with residuals
            GARP_calculate_flag = 1;

        case S.calculate_GARP_choice(2)
            % without residuals
            GARP_calculate_flag = 0;
    end

    switch findobj(get(S.residual_GARP,'selectedobject'))
        case S.residual_GARP_choice(1)
            % with residuals
            GARP_residual_flag = 1;

        case S.residual_GARP_choice(2)
            % without residuals
            GARP_residual_flag = 0;
    end

    % in and out of sample
    %if ~(get(S.residual_GARP_choice(2),'value') == 1.0)
    if get(S.residual_GARP_options(1), 'value') == 1.0
        GARP_in_sample_flag = 1;
    elseif get(S.residual_GARP_options(2), 'value') == 1.0
        GARP_out_sample_flag = 1;
    elseif get(S.residual_GARP_options(3), 'value') == 1.0
        GARP_in_sample_flag = 1;
        GARP_out_sample_flag = 1;
    end
    %end

    % which residuals (WARP / GARP / SARP)
    %if ~(get(S.residual_GARP_choice(2),'value') == 1.0)
    if get(S.which_residual_GARP_options(1), 'value') == 1.0
        WARP_res = 1;
    elseif get(S.which_residual_GARP_options(2), 'value') == 1.0
        GARP_res = 1;
    elseif get(S.which_residual_GARP_options(3), 'value') == 1.0
        SARP_res = 1;
    elseif get(S.which_residual_GARP_options(4), 'value') == 1.0
        WARP_res = 1;
        GARP_res = 1;
        SARP_res = 1;
    end
    %end

    GARP_flags = [GARP_calculate_flag , GARP_residual_flag , GARP_in_sample_flag , GARP_out_sample_flag, WARP_res, GARP_res, SARP_res];


    %% choices regarding MPI consistency index

    MPI_calculate_flag = 0;
    MPI_residual_flag = 0;
    MPI_in_sample_flag = 0;
    MPI_out_sample_flag = 0;

    % MPI will be ready for the next version (version 3), currently not implemented 
%     switch findobj(get(S.calculate_MPI,'selectedobject'))
%         case S.calculate_MPI_choice(1)
%             % with residuals
%             MPI_calculate_flag = 1;
% 
%         case S.calculate_MPI_choice(2)
%             % without residuals
%             MPI_calculate_flag = 0;
%     end
% 
%     switch findobj(get(S.residual_MPI,'selectedobject'))
%         case S.residual_MPI_choice(1)
%             % with residuals
%             MPI_residual_flag = 1;
% 
%         case S.residual_MPI_choice(2)
%             % without residuals
%             MPI_residual_flag = 0;
%     end
% 
%     % in and out of sample (only out of sample is relevant to MPI) 
%     %if ~(get(S.residual_MPI_choice(2),'value') == 1.0)
%     if get(S.residual_MPI_options(1), 'value') == 1.0
%         MPI_out_sample_flag = 1;
%     end
%     %end

    MPI_flags = [MPI_calculate_flag , MPI_residual_flag , MPI_in_sample_flag , MPI_out_sample_flag];

    
    %% choices regarding AFRIAT consistency index

    AFRIAT_calculate_flag = 0;
    AFRIAT_residual_flag = 0;
    AFRIAT_in_sample_flag = 0;
    AFRIAT_out_sample_flag = 0;

    switch findobj(get(S.calculate_AFRIAT,'selectedobject'))
        case S.calculate_AFRIAT_choice(1)
            % with residuals
            AFRIAT_calculate_flag = 1;

        case S.calculate_AFRIAT_choice(2)
            % without residuals
            AFRIAT_calculate_flag = 0;
    end

    switch findobj(get(S.residual_AFRIAT,'selectedobject'))
        case S.residual_AFRIAT_choice(1)
            % with residuals
            AFRIAT_residual_flag = 1;

        case S.residual_AFRIAT_choice(2)
            % without residuals
            AFRIAT_residual_flag = 0;
    end

    % in and out of sample (only out of sample is relevant to AFRIAT) 
    if ~(get(S.residual_AFRIAT_choice(2),'value') == 1.0)
        if get(S.residual_AFRIAT_options(1), 'value') == 1.0
            AFRIAT_out_sample_flag = 1;
        end
    end

    AFRIAT_flags = [AFRIAT_calculate_flag , AFRIAT_residual_flag , AFRIAT_in_sample_flag , AFRIAT_out_sample_flag];


    %% choices regarding VARIAN consistency index

    VARIAN_calculate_flag = 0;
    VARIAN_residual_flag = 0;
    VARIAN_in_sample_flag = 0;
    VARIAN_out_sample_flag = 0;

    switch findobj(get(S.calculate_VARIAN,'selectedobject'))
        case S.calculate_VARIAN_choice(1)
            % with residuals
            VARIAN_calculate_flag = 1;

        case S.calculate_VARIAN_choice(2)
            % without residuals
            VARIAN_calculate_flag = 0;
    end

    switch findobj(get(S.residual_VARIAN,'selectedobject'))
        case S.residual_VARIAN_choice(1)
            % with residuals
            VARIAN_residual_flag = 1;

        case S.residual_VARIAN_choice(2)
            % without residuals
            VARIAN_residual_flag = 0;
    end

    %if ~(get(S.residual_VARIAN_choice(2),'value') == 1.0)
    % in and out of sample
    if get(S.residual_VARIAN_options(1), 'value') == 1.0
        VARIAN_in_sample_flag = 1;
    elseif get(S.residual_VARIAN_options(2), 'value') == 1.0
        VARIAN_out_sample_flag = 1;
    elseif get(S.residual_VARIAN_options(3), 'value') == 1.0
        VARIAN_in_sample_flag = 1;
        VARIAN_out_sample_flag = 1;
    end
    %end

    VARIAN_flags = [VARIAN_calculate_flag , VARIAN_residual_flag , VARIAN_in_sample_flag , VARIAN_out_sample_flag];


    %% choices regarding HOUTMAN-MAKS inconsistency index

    HOUTMAN_calculate_flag = 0;
    HOUTMAN_residual_flag = 0;

    switch findobj(get(S.calculate_HOUTMAN,'selectedobject'))
        case S.calculate_HOUTMAN_choice(1)
            % with residuals
            HOUTMAN_calculate_flag = 1;

        case S.calculate_HOUTMAN_choice(2)
            % without residuals
            HOUTMAN_calculate_flag = 0;
    end

    switch findobj(get(S.residual_HOUTMAN,'selectedobject'))
        case S.residual_HOUTMAN_choice(1)
            % with residuals
            HOUTMAN_residual_flag = 1;

        case S.residual_HOUTMAN_choice(2)
            % without residuals
            HOUTMAN_residual_flag = 0;
    end

    % these 2 flags are only for consistency with the other indices
    HOUTMAN_in_sample_flag = 0;
    HOUTMAN_out_sample_flag = HOUTMAN_residual_flag;
    
    HOUTMAN_flags = [HOUTMAN_calculate_flag , HOUTMAN_residual_flag , HOUTMAN_in_sample_flag , HOUTMAN_out_sample_flag];


    %% close the window
    close(S.fh);

end



%% Cancel button
function [] = pc_call(varargin)
    % Callback for Cancel button.

    ok = 0;

    % close the window
    close(S.fh);

end


    
end % end of function