function [ok, another_run] = HPZ_Screen_Another_Run (runs_counter)

% this function promotes a user-interface screen to the user in the end of
% all the rest of the screens, asking the user whether she desires to
% define another run, or should the program start running all the runs
% defined so far.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 

% it is recommended to use the numeric variables in the beginning of the 
% function when making changes to the screen.
% it is also recommended to use numeric variables in the same manner when
% adding new elements to the screen.



% dimensions of the screen
full_width = 460;
full_height = 130;

% distance of highest element from top
top_dist = 25;
% height of header label
label_height = 25;

% buttons absolute height
buttons_abs_height = 30;
% the distant between a button to the edge of the screen,
% and is also half the distance between buttons.
% the size of the buttons is designed to fit this and the buttons_height 
% and the width of the screen
buttons_dists = 10;
% number of buttons
buttons_num = 3;

% font size of label
label_font_size = 12;
% general font size
font_size = 9;



%% create the gui figure
sz = [full_width , full_height]; % figure size
screensize = get(0,'ScreenSize'); % screen size
xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the center
ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the center
S.fh = figure('units','pixels',...
    'position',[xpos, ypos, sz(1), sz(2)],...
    'menubar','none',...
    'name',char(strcat('Running Plan (there are currently', {' '}, num2str(runs_counter), ' runs)')),...
    'numbertitle','off',...
    'resize','off');



% current bottom value
current_bottom = full_height-label_height-top_dist;
%% Head Label
S.label_DA = uicontrol('style','text',...
    'units','pix',...
    'position',[0 , current_bottom , full_width , label_height],...
    'backgroundc',get(S.fh,'color'),...
    'fontsize',label_font_size,'fontweight','bold',...
    'string','Define Another Run, or Start Running?');



% current bottom value (so that buttons will be vertically in the middle) 
current_bottom = (current_bottom - buttons_abs_height) / 2;

%% Run Now Button
S.run_now_button = uicontrol('style','push',...
    'unit','pix',...
    'position',[buttons_dists , current_bottom , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_abs_height],...
    'string','Run Now',...
    'fontsize',font_size,...
    'callback',{@run_now_call,S});

%% Another Run Button
S.another_run_button = uicontrol('style','push',...
    'unit','pix',...
    'position',[(full_width/3)+buttons_dists , current_bottom , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_abs_height],...
    'string','Define Another Run',...
    'fontsize',font_size,...
    'callback',{@another_run_call,S});

%% Cancel All Button
S.cancel_all_button = uicontrol('style','push',...
    'unit','pix',...
    'position',[(full_width*2/3)+buttons_dists , current_bottom , (full_width-(buttons_num+1)*buttons_dists)/buttons_num , buttons_abs_height],...
    'string','Cancel All',...
    'fontsize',font_size,...
    'callback',{@cancel_all_call,S});





uiwait(S.fh)  % Prevent all other processes from starting until closed.





%% Callback for Run Now button
function [] = run_now_call(varargin)

    ok = 1;
    another_run = false;
    
    % close the window
    close(S.fh);

end

%% Callback for Run Now button
function [] = another_run_call(varargin)

    ok = 1;
    another_run = true;
    
    % close the window
    close(S.fh);

end

%% Callback for Cancel All button
function [] = cancel_all_call(varargin)

    ok = 0;
    another_run = false;

    % close the window
    close(S.fh);

end



end