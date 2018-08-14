function [Xs,NumberXs,MaxX]=Collect(DataColumn)

% function [Xs,NumberXs,MaxX]=Collect(DataColumn)
% Extracts the distinct (unique) values from a column of input data,
% counts the number of distinct values, and finds the maximum value.

% Version: 3f
% Date: July 8, 2009

Xs=sort(unique(DataColumn)); % Collects distinct rows (and performs sort)
NumberXs=length(Xs); % Counts number of distinct rows
MaxX=max(Xs); % Finds max value among distinct rows

clear DataColumn

end