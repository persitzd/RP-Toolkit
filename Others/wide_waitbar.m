function wb = wide_waitbar(bar_value, message, name, width_multiplier, varargin)
% this function creates a wider waitbar, that is: 
% a waitbar that its width is: (default_width * width_multiplier)

% shifts to the right (negative if to the left) and to the bottom (negative
% if to the top), in percentage of the size screen (e.g. 0.1 = 10% of screen size) 
if ~isempty(varargin)
    shifts = varargin{1};
else
    shifts = [0,0];
end

% Create the waitbar
wb = waitbar(bar_value, message, 'name', name);

if str2double(strip(strip(version('-release'), 'a'), 'b')) >= 2019
    % do nothing; starting 2019b it is not possible to change the size of
    % the waitbar after it is created (causes an error and a crash)
else
    % Change the Units Poreperty of the figure and all the children
    set(findall(wb),'Units', 'normalized');
    % Change the size of the figure
    wb_position = get(wb, 'Position');
    set(wb,'Position', [wb_position(1)-wb_position(3)*(width_multiplier-1)/2+shifts(1) , wb_position(2)-shifts(2) , wb_position(3)*width_multiplier , wb_position(4)]);
end

% this command is needed in order to refresh the waitbar graphics immediately 
waitbar(bar_value, wb);

end