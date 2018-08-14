function [base_figure, base_panel] = ui_scroll_screen(figure_width, figure_height, scroll_width, max_height_percent, top_space_height, bottom_space_height, title)

% this function creates a screen with a built-in scroll procedure.
% the function is given a height for the screen, and a percent representing
% the maximum percent of the pc screen that the actual height of the screen
% should not pass.
% if the percentage limitation is smaller than the desired height - the
% function will create a screen that has the maximum height that can be, 
% and will create a slide so that the rest of the desired height will be
% reachable to the user.
% if the limit is not a problem, no scroll will be created

% the function returns both the screen itself, and the panel that can be
% scrolled (or can't if scroll was not needed).

% input:
%   figure_width - width of the screen (NOT including the width for the scroll bar) 
%   figure_height - required height for the screen (the screen itself might be smaller, but the scroll will allow this height) 
%   scroll_width - width of the scroll bar
%   max_height_percent - maximum height of the screen in percent of the computer's screen height 
%   top_space_height - height of area on the top of the screen that should not be scrollabe (e.g. used for header) 
%   bottom_space_height - height of area on the bottom of the screen that should not be scrollabe (e.g. used for ok/cancel buttons)  
%   title - title of the screen



% default number of steps of scroll UI
scroll_steps = 1000;

% we add a little epsilon to the panels, so that the right edge of the
% panel won't stand out
panel_epsilon = 2;

% dimensions of the screen of the computer
pc_screen_size = get(0,'ScreenSize');
pc_screen_width = pc_screen_size(3);
pc_screen_height = pc_screen_size(4);

if (figure_height <= pc_screen_height * max_height_percent)
    % there's no need in a scroll
    adjusted_width = figure_width;
    adjusted_height = figure_height;
    create_scroll = false;
else
    % we need a scroll
    adjusted_width = figure_width + scroll_width;
    adjusted_height = pc_screen_height * max_height_percent;
    create_scroll = true;
end

width_pos = ceil((pc_screen_width - adjusted_width)/2);     % center the figure on the screen 
height_pos = ceil((pc_screen_height - adjusted_height)/2);  % center the figure on the screen 


%% create the gui figure
base_figure = figure('units','pixels',...
              'position',[width_pos , height_pos , adjusted_width , adjusted_height],...
              'menubar','none',...
              'name',title,...
              'numbertitle','off',...
              'resize','off');
    
%% panel to contain the panel that will be "slided"
parent_panel = uipanel('Parent',base_figure,...
                'backgroundc',get(base_figure,'color'),...
                'units','pix',...
                'position',[adjusted_width-figure_width , bottom_space_height , figure_width+panel_epsilon , adjusted_height-bottom_space_height-top_space_height]);

if (~create_scroll)    
    base_panel = parent_panel;  
else         
    %% slider object
    slider_panel = uicontrol('Parent',base_figure, 'style','slider',...
                    'backgroundc',get(base_figure,'color'),...
                    'Min',0,'Max',scroll_steps,'Value',scroll_steps,...
                    'units','pix',...
                    'position',[0 , bottom_space_height , scroll_width , adjusted_height-bottom_space_height-top_space_height],...
                    'Callback', @slide_panel); %#ok<NASGU>
                
    %% panel to contain all elements, to be move by "sliding" 
    sliding_panel = uipanel('Parent',parent_panel,...
                    'backgroundc',get(base_figure,'color'),...
                    'units','pix',...
                    'position',[0 , adjusted_height-figure_height , figure_width+panel_epsilon , figure_height-bottom_space_height-top_space_height]);
                
    base_panel = sliding_panel;
end
  


%% slider function
function [] = slide_panel(source, event) %#ok<INUSD>
    val = scroll_steps - source.Value;
    step_size = (figure_height-adjusted_height) / scroll_steps;

    set(sliding_panel, 'position', [0 , adjusted_height-figure_height+val*step_size , figure_width , figure_height-bottom_space_height-top_space_height]);
end



end