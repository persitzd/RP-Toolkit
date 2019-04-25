function [aggregation_flag, metric_flag, max_time_estimation, min_counter, parallel_flag, max_starting_points, ok] = HPZ_Screen_Optimization_Settings(aggregation_flag, metric_flag, max_time_estimation, min_counter, parallel_flag, action_flag, numeric_flag, runs_counter)   

% this function promotes a user-interface screen to the user, allowing her
% to make choices regarding the time limit and/or others limits on the
% optimization process, as well as choosing the aggregation or metric to be
% used, depending on the chosen method.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

% it is recommended to use the numeric variables in the beginning of the 
% function when making changes to the screen.
% it is also recommended to use numeric variables in the same manner when
% adding new elements to the screen.



% this little cell helps us to convert 0 and 1 to 'off' and 'on', respectively 
enable = {'off','on'};



% maximal time of estimation per subject
if isinf(max_time_estimation) || max_time_estimation == 0
    max_time_str = '';
else
    max_time_str = num2str(max_time_estimation);
end



% NLLS
if action_flag == HPZ_Constants.NLLS_action
    enable_string_NLLS = 'on';
    enable_string_BI = 'off';
    enable_string_MMI = 'off';

% MMI
elseif action_flag == HPZ_Constants.MMI_action
    enable_string_NLLS = 'off';
    enable_string_BI = 'off';
    enable_string_MMI = 'on';

% BI
elseif action_flag == HPZ_Constants.BI_action
    enable_string_NLLS = 'off';
    enable_string_BI = 'on';
    enable_string_MMI = 'off';
end



% max_starting_points is initialized here, because it is related to 
%   min_counter that is determined by the user in this screen.
% also, if in the future one would want to let the user choose the
%   max_starting_points, this screen is the one it should be in.
if (numeric_flag == HPZ_Constants.numeric)
    max_starting_points = HPZ_Constants.max_starting_points_numeric;
elseif (numeric_flag == HPZ_Constants.analytic)
    max_starting_points = HPZ_Constants.max_starting_points_analytic;
elseif (numeric_flag == HPZ_Constants.semi_numeric)
    max_starting_points = HPZ_Constants.max_starting_points_semi_numeric;
end

% this str will appear to the user as an option for min_counter, meaning
% that min_counter (number of convergence points) should not limit the
% estimation progress
inf_str = 'inf';

% cell array of all possible values for min_counter
min_counter_values = HPZ_Constants.min_counter_values;
min_counter_size = max(size(min_counter_values));
if min_counter_values{min_counter_size} < max_starting_points
    min_counter_size = min_counter_size + 1;
    min_counter_values{min_counter_size} = num2str(max_starting_points);
end
%min_counter_size = min_counter_size + 1;
%min_counter_values{min_counter_size} = inf_str;

% finding the initial index for the screen
min_counter_vector = cell2mat(cellfun(@str2num, min_counter_values(1:end), 'un', 0).');
index_vector = 1:min_counter_size;
if min_counter == inf
    min_counter_index = min_counter_size;
else
    min_counter_index = index_vector(min_counter_vector == min_counter);
    % this if isempty is for possible problems and bugs
    if isempty(min_counter_index)
        min_counter_index = max(1, min_counter_size-1);
    end
end





% screen size
figure_width = 470;
figure_height = 520;
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

% how much distance between parts (different methods) that are one below each other 
move_down = 20;

% distance of the labels from the left
%left_label = 15;
% distance of other elements from the left
%left_other = 45;

% each radio options will be with height 25
radio_height = 25;
% distance of bottom radio option from the bottom
radio_bottom = 10;
% offset to the right of radio options that are vertical
%radio_offset = 15;
% offsets to the right of yes/no horizontal radio options
%yes_no_offsets = [25 , 180];
% offsets to the right of one/two/three horizontal radio options
three_offsets = [10 , 129, 230];
three_offsets2 = [10 , 117, 250];

% width of each element
element_width = 400;
% normal_width for a sub element
sub_element_width = 120;
sub_element_width2 = 180;
% label width for aggregation / metric
label_width = 170;

% height of elements that allows to used to type in a number
text_height = 25;

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
figure_title = char(strcat('Optimization Settings (', HPZ_Constants.current_run_screen, {' '}, num2str(runs_counter), ')'));
[S.fh , S.panel] = ui_scroll_screen(figure_width, figure_height, scroll_width, max_height_percent, top_space_height, bottom_space_height, figure_title);

% width including scroll bar if there is one
pos = get(S.fh,'position');
full_width = pos(3);



            
            
%% Optimization settings MMI
% current height of element
current_height = 30;
% current bottom coordinate
current_bottom = panel_height-top_dist - current_height;
S.label_MMI = uicontrol('Parent',S.panel, ...
    'style','text',...
    'units','pix',...
    'position',[15 , current_bottom , 250 , current_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',label_font_size,'fontweight','bold',...
    'HorizontalAlignment','left',...
    'string','Money Metric Index Method');

% current height of element
current_height = 40;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.initial_points_MMI = uicontrol('Parent',S.panel, ...
    'style','text',...
    'enable', enable_string_MMI, ...
    'units','pix',...
    'position',[15 , current_bottom , 150 , current_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',big_font_size,...
    'string','Number of convergence points:');

S.initial_points_MMI_pop = uicontrol('Parent',S.panel, ...
    'style','pop',...
    'value',min_counter_index,...
    'enable', enable_string_MMI, ...
    'units','pixels',...
    'position',[170 , current_bottom , 40 , current_height],...
    'backgroundc',get(S.fh,'color'),...
    'string',min_counter_values,...
    'fontsize',font_size,...
    'callback',{@on_convergence_change_check_time,S});


S.estimated_time_MMI_label = uicontrol('Parent',S.panel, ...
    'style','text',...
    'units','pix',...
    'enable', enable_string_MMI, ...
    'position',[220 , current_bottom , 150 , current_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',big_font_size,...
    'string','Allocated time (minutes):');

S.estimated_time_MMI_text = uicontrol('Parent',S.panel, ...
    'style','edit',...
    'string',max_time_str,...
    'enable', enable_string_MMI, ...
    'units','pix',...
    'position',[375 , current_bottom+current_height/2 , 80 , text_height],...
    'backgroundc','white',...
    'fontsize',big_font_size,...
    'callback',{@on_time_change_check_convergence,S});


%% Aggregation MMI
% current height of element
current_height = 20;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.ed_MMI = uicontrol('Parent',S.panel, ...
    'style','edit',...
    'unit','pix',...
    'enable', enable_string_MMI, ...
    'position',[20 , current_bottom , label_width , current_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',font_size,...
    'string','Aggregation  Method');

% current height of element
current_height = 40;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.bg_MMI = uibuttongroup('Parent',S.panel, ...
    'units','pix',...
    'pos',[20 , current_bottom , element_width , current_height]);

S.Aggregation_MMI_rd(1) = uicontrol(S.bg_MMI,...
    'value',(aggregation_flag == HPZ_Constants.MMI_Max),...
    'enable', enable_string_MMI, ...
    'style','rad',...
    'unit','pix',...
    'position',[three_offsets(1) , radio_bottom , sub_element_width , radio_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',font_size,...
    'string',' Maximum Waste');

S.Aggregation_MMI_rd(2) = uicontrol(S.bg_MMI,...
    'value',(aggregation_flag == HPZ_Constants.MMI_Mean),...
    'enable', enable_string_MMI, ...
    'style','rad',...
    'unit','pix',...
    'position',[three_offsets(2) , radio_bottom , sub_element_width , radio_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',font_size,...
    'string',' Mean Waste');

S.Aggregation_MMI_rd(3) = uicontrol(S.bg_MMI,...
    'value',(aggregation_flag == HPZ_Constants.MMI_AVGSSQ),...
    'enable', enable_string_MMI, ...
    'style','rad',...
    'unit','pix',...
    'position',[three_offsets(3) , radio_bottom , sub_element_width , radio_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',font_size,...
    'string',' AVG(SSQ(Wastes))');
%% End of MMI settings



%% Optimization settings BI

% current bottom coordinate (BI is just below MMI)
current_bottom = current_bottom - move_down;

% current height of element
current_height = 30;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.label_BI = uicontrol('Parent',S.panel, ...
                'style','text',...
                'units','pix',...
                'position',[10 , current_bottom , 175 , current_height],...
                'backgroundc',get(S.fh,'color'),...
                'fontsize',label_font_size,'fontweight','bold',... 
                'string','Binary Index Method');

% current height of element
current_height = 40;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.initial_points_BI = uicontrol('Parent',S.panel, ...
                'style','text',...
                'enable', enable_string_BI, ...
                'units','pix',...
                'position',[15 current_bottom 150 current_height],...
                'backgroundc',get(S.fh,'color'),...
                'fontsize',big_font_size,...
                'string','Number of convergence points:');

S.initial_points_BI_pop = uicontrol('Parent',S.panel, ...
                'style','pop',...
                'value',min_counter_index,...
                'units','pixels',...
                'enable', enable_string_BI, ...
                'position',[170 , current_bottom , 40 , current_height],...
                'fontsize',font_size,...
                'string',min_counter_values,...
                'callback',{@on_convergence_change_check_time,S});

S.estimated_time_BI_label = uicontrol('Parent',S.panel, ...
                'style','text',...
                'units','pix',...
                'enable', enable_string_BI, ...
                'position',[220 , current_bottom , 150 , current_height],...
                'backgroundc',get(S.fh,'color'),...
                'fontsize',big_font_size,...
                'string','Allocated time (minutes):');

S.estimated_time_BI_text = uicontrol('Parent',S.panel, ...
                'style','edit',...
                'string',max_time_str,...
                'enable', enable_string_BI, ...
                'units','pix',...
                'position',[375 , current_bottom+current_height/2 , 80 , text_height],...
                'backgroundc','white',...
                'fontsize',big_font_size,...
                'callback',{@on_time_change_check_convergence,S});
%% End of BI settings



%% Optimization settings NLLS

% current bottom coordinate (NLLS is just below BI)
current_bottom = current_bottom - move_down;

% current height of element
current_height = 30;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.label_NLLS = uicontrol('Parent',S.panel, ...
                'style','text',...
                'units','pix',...
                'position',[15 , current_bottom , 190 , current_height],...
                'backgroundc',get(S.fh,'color'),...
                'fontsize',label_font_size,'fontweight','bold',... 
                'string','Nonlinear Least Squares');

% current height of element
current_height = 40;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.initial_points_NLLS = uicontrol('Parent',S.panel, ...
                'style','text',...
                'enable', enable_string_NLLS, ...
                'units','pix',...
                'position',[15 , current_bottom , 150 , current_height],...
                'backgroundc',get(S.fh,'color'),...
                'fontsize',big_font_size,...
                'string','Number of convergence points:');

S.initial_points_NLLS_pop = uicontrol('Parent',S.panel, ...
                'style','pop',...
                'value',min_counter_index,...
                'enable', enable_string_NLLS, ...
                'units','pixels',...
                'position',[170 , current_bottom , 40 , current_height],...
                'string',min_counter_values,...
                'fontsize',font_size,...
                'callback',{@on_convergence_change_check_time,S});

S.estimated_time_NLLS_label = uicontrol('Parent',S.panel, ...
                'style','text',...
                'enable', enable_string_NLLS, ...
                'units','pix',...
                'position',[220 , current_bottom , 150 , current_height],...
                'backgroundc',get(S.fh,'color'),...
                'fontsize',big_font_size,...
                'string','Allocated time (minutes):');

S.estimated_time_NLLS_text = uicontrol('Parent',S.panel, ...
                'style','edit',...
                'string',max_time_str,...
                'units','pix',...
                'enable', enable_string_NLLS, ...
                'position',[375 , current_bottom+current_height/2 , 80 , text_height],...
                'backgroundc','white',...
                'fontsize',big_font_size,...
                'callback',{@on_time_change_check_convergence,S});

%% Metric Method NLLS
% current height of element
current_height = 20;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.ed_NLLS = uicontrol('Parent',S.panel, ...
                'style','edit',...
                'enable', enable_string_NLLS, ...
                'unit','pix',...
                'position',[20 , current_bottom , label_width , current_height],...
                'backgroundc',get(S.fh,'color'),...
                'fontsize',font_size,...
                'string','Metric  Selection');

% current height of element
current_height = 40;
% current bottom coordinate
current_bottom = current_bottom - current_height;
S.bg_NLSS = uibuttongroup('Parent',S.panel, ...
                'units','pix',...
                'pos',[20 , current_bottom , 430 , current_height]);

S.metric_NLLS_rd(1) = uicontrol(S.bg_NLSS,...
                'value',(metric_flag == HPZ_Constants.euclidean_metric),...
                'enable', enable_string_NLLS, ...
                'style','rad',...
                'unit','pix',...
                'position',[three_offsets2(1) , radio_bottom , sub_element_width2 , radio_height],...
                'fontsize',font_size,...
                'string',' Euclidean metric');

S.metric_NLLS_rd(2) = uicontrol(S.bg_NLSS,...
                'value',(metric_flag == HPZ_Constants.CFGK_metric),...
                'enable', enable_string_NLLS, ...
                'style','rad',...
                'unit','pix',...
                'position',[three_offsets2(2) , radio_bottom , sub_element_width2 , radio_height],...
                'fontsize',font_size,...
                'string',' Choi et al. (2007) metric');
            
S.metric_NLLS_rd(3) = uicontrol(S.bg_NLSS,...
                'value',(metric_flag == HPZ_Constants.normalized_euclidean_metric),...
                'enable', enable_string_NLLS, ...
                'style','rad',...
                'unit','pix',...
                'position',[three_offsets2(3) , radio_bottom , sub_element_width2 , radio_height],...
                'fontsize',font_size,...
                'string',' normalized-Euclidean metric');
%% End of NLLS settings

                

            
            
%% parallel computing                    
S.parallel_label = uicontrol('Parent',S.panel, ...
                'style','text',...
                'units','pix',...
                'position',[15 , 40 , 150 , 30],...
                'backgroundc',get(S.fh,'color'),...
                'fontsize',label_font_size,'fontweight','bold',... 
                'string','Parallel Processing');

% if analytic was chosen (numeric_flag = false), we never use parallel computing 
S.parallel_label_ch = uicontrol('Parent',S.panel, ...
                'style','check',...
                'unit','pix',...
                'position',[175 , 47 , 275 , 25],...
                'string',' Use matlab parallel computing package.',...
                'fontsize',font_size,...
                'value', ((numeric_flag == HPZ_Constants.numeric) && parallel_flag), ...
                'enable', enable{(numeric_flag == HPZ_Constants.numeric)+1},...
                'fontsize',10);

            
            



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


function [] = on_time_change_check_convergence(varargin)
    
    NLLS_time = str2num(get(S.estimated_time_NLLS_text,'string'));
    NLLS_convergence = str2num(min_counter_values{get(S.initial_points_NLLS_pop,'value')});
    if (isempty(NLLS_time) || isinf(NLLS_time)) && isinf(NLLS_convergence)
        % we changed time to inf, but max convergence points is also inf
        set(S.initial_points_NLLS_pop,'value',max(size(min_counter_values))-1);
    end
    MMI_time = str2num(get(S.estimated_time_MMI_text,'string'));
    MMI_convergence = str2num(min_counter_values{get(S.initial_points_MMI_pop,'value')});
    if (isempty(MMI_time) || isinf(MMI_time)) && isinf(MMI_convergence)
        % we changed time to inf, but max convergence points is also inf
        set(S.initial_points_MMI_pop,'value',max(size(min_counter_values))-1);
    end
    BI_time = str2num(get(S.estimated_time_BI_text,'string'));
    BI_convergence = str2num(min_counter_values{get(S.initial_points_BI_pop,'value')});
    if (isempty(BI_time) || isinf(BI_time)) && isinf(BI_convergence)
        % we changed time to inf, but max convergence points is also inf
        set(S.initial_points_BI_pop,'value',max(size(min_counter_values))-1);
    end
end

function [] = on_convergence_change_check_time(varargin)
    
    arbitrary_time = 30;
    
    NLLS_time = str2num(get(S.estimated_time_NLLS_text,'string'));
    NLLS_convergence = str2num(min_counter_values{get(S.initial_points_NLLS_pop,'value')});
    if (isempty(NLLS_time) || isinf(NLLS_time)) && isinf(NLLS_convergence)
        % we changed time to inf, but max convergence points is also inf
        set(S.estimated_time_NLLS_text,'string',arbitrary_time);
    end
    MMI_time = str2num(get(S.estimated_time_MMI_text,'string'));
    MMI_convergence = str2num(min_counter_values{get(S.initial_points_MMI_pop,'value')});
    if (isempty(MMI_time) || isinf(MMI_time)) && isinf(MMI_convergence)
        % we changed time to inf, but max convergence points is also inf
        set(S.estimated_time_MMI_text,'string',arbitrary_time);
    end
    BI_time = str2num(get(S.estimated_time_BI_text,'string'));
    BI_convergence = str2num(min_counter_values{get(S.initial_points_BI_pop,'value')});
    if (isempty(BI_time) || isinf(BI_time)) && isinf(BI_convergence)
        % we changed time to inf, but max convergence points is also inf
        set(S.estimated_time_BI_text,'string',arbitrary_time);
    end
end





%% OK button
function [] = ok_button_call(varargin)
    
    ok = 1;

    % Callback for the popup.
    S = varargin{3};  % Get the structure.

    if strcmpi(get(S.initial_points_NLLS, 'enable'), 'on') % NLLS
        rd_vals_NLLS = get (S.metric_NLLS_rd(:), 'value');
        if rd_vals_NLLS{1} == 1
            % Euclidean metric
            metric_flag = HPZ_Constants.euclidean_metric;
        elseif rd_vals_NLLS{2} == 1
            % Choi (2007) metric
            metric_flag = HPZ_Constants.CFGK_metric;
        elseif rd_vals_NLLS{3} == 1
            % normalized-Euclidean metric
            metric_flag = HPZ_Constants.normalized_euclidean_metric;
        end % end switch

        % set thte number of convergence points for NLLS
        P_NLLS = get(S.initial_points_NLLS_pop,{'string','val'});  % Get the users choice.
        if strcmp(P_NLLS{1}{P_NLLS{2}}, inf_str)
            min_counter = inf;
        else
            min_counter = str2num(P_NLLS{1}{P_NLLS{2}}); %#ok<*ST2NM>
        end
        
        % get the time limit. if not specified - set it to inf (no time limit) 
        if isempty(get(S.estimated_time_NLLS_text,'string')) || strcmp(get(S.estimated_time_NLLS_text,'string'), '') || strcmp(get(S.estimated_time_NLLS_text,'string'), '0')
            max_time_estimation = Inf;
        else
            max_time_estimation = str2num(get(S.estimated_time_NLLS_text,'string'));
        end

    elseif strcmpi(get(S.initial_points_MMI, 'enable'), 'on') % MMI
        rd_vals_MMI = get (S.Aggregation_MMI_rd(:), 'value');
        if rd_vals_MMI{1} == 1
                % Max Waste
                aggregation_flag = 1;

        elseif rd_vals_MMI{2} == 1
                % Mean Waste
                aggregation_flag = 2;

        else
                % Sum of Squares Wastes
                aggregation_flag = 3;
        end % end of switch

        % set the number of convergence points MMI
        P_MMI = get(S.initial_points_MMI_pop,{'string','val'});  % Get the users choice.
        if strcmp(P_MMI{1}{P_MMI{2}}, inf_str)
            min_counter = inf;
        else
            min_counter = str2num(P_MMI{1}{P_MMI{2}});
        end

        if isempty(get(S.estimated_time_MMI_text,'string')) || strcmp(get(S.estimated_time_MMI_text,'string'), '') || strcmp(get(S.estimated_time_MMI_text,'string'), '0')
            max_time_estimation = inf;
        else
            max_time_estimation = str2num(get(S.estimated_time_MMI_text,'string'));
        end

    elseif strcmpi(get(S.initial_points_BI, 'enable'), 'on') % BI

        % set the number of convergence points MMI
        P_BI = get(S.initial_points_BI_pop,{'string','val'});  % Get the users choice.
        if strcmp(P_BI{1}{P_BI{2}}, inf_str)
            min_counter = inf;
        else
            min_counter = str2num(P_BI{1}{P_BI{2}});
        end

        if isempty(get(S.estimated_time_BI_text,'string')) || strcmp(get(S.estimated_time_BI_text,'string'), '') || strcmp(get(S.estimated_time_BI_text,'string'), '0')
            max_time_estimation = inf;
        else
            max_time_estimation = str2num(get(S.estimated_time_BI_text,'string'));
        end

    end % end of if
    
    if get(S.parallel_label_ch, 'value') == 1.0   % Parallel Computing
        % Use matlabpool
        parallel_flag = true;
    else
        % Use serial version
        parallel_flag = false;
    end
    
    % close the window
    close(S.fh);

end % end of call back function



%% Cancel button
function [] = cancel_button_call(varargin)
    % Callback for Cancel button.

    ok = 0;

    % close the window
    close(S.fh);

end



end % end of function
