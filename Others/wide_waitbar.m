function wb = wide_waitbar(bar_value, message, name, width_multiplier)
% this function creates a wider waitbar, that is: 
% a waitbar that its width is: (default_width * width_multiplier)

% Create the waitbar
wb = waitbar(bar_value, message, 'name', name);

% Change the Units Poreperty of the figure and all the children
set(findall(wb),'Units', 'normalized');

% Change the size of the figure
wb_position = get(wb, 'Position');
set(wb,'Position', [wb_position(1)-wb_position(3)*(width_multiplier-1)/2 , wb_position(2) , wb_position(3)*width_multiplier , wb_position(4)]);

% this command is needed in order to refresh the waitbar graphics immediately 
waitbar(bar_value, wb);

end