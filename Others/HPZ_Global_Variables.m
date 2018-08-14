function HPZ_Global_Variables(chosen_subjects_num)

% we try to avoid as much as possible from using global variables,
% but in some cases (such as wanting to get values from functions that are
% called through fminsearchbnd) we have no choice but to use them.

% in order to keep things under control, the global variables will all be
% initialized through here



% this cell array keeps the warnings for a single subject
% global warnings_cell_one_subject
% warnings_cell_one_subject = {};

% numbers (counters) of warnings per subject due to:
% Criterion = NaN
global warnings_nan
warnings_nan = zeros(chosen_subjects_num,1);
% Criterion = -Inf
global warnings_minus_inf
warnings_minus_inf = zeros(chosen_subjects_num,1);
% Criterion = +Inf
global warnings_plus_inf
warnings_plus_inf = zeros(chosen_subjects_num,1);
% Criterion > 1 + threshold (MMI only)
global warnings_bigger_than_1
warnings_bigger_than_1 = zeros(chosen_subjects_num,1);
% Criterion < 0 - threshold (MMI only)
global warnings_smaller_than_0
warnings_smaller_than_0 = zeros(chosen_subjects_num,1);
% (index of current subject)
global current_subject
current_subject = 1;



end