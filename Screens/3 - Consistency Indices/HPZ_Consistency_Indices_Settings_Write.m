function HPZ_Consistency_Indices_Settings_Write(Graph_flags, power_test_settings, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags, main_folder)

% This function is called when the user wishes (and it is also the
% default...) to save the choices he made now for the next time.

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% make sure the Settings Files directory exists, if not - create it
dir_exists = exist(strcat(main_folder, '/', HPZ_Constants.settings_files_dir) , 'dir');
if ~dir_exists
    mkdir(HPZ_Constants.settings_files_dir);
end

% creating/opening the settings file
settings_file = fopen(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.consistency_indices_settings_file_name, '.csv'), 'w');

% close file when this function end
cleanup_close_file = onCleanup(@() fclose(settings_file));



%% printing the new settings and choices

% Revealed Preference Graphs (per subject) - whether to create them and as .png or as .fig or both 
for i=1:max(size(Graph_flags))
    fprintf(settings_file, '%s,', num2str(Graph_flags(i)));   % line (1)
end
fprintf(settings_file, '\n');
% power test - perform or not, and number of simulations
for i=1:max(size(power_test_settings))
    fprintf(settings_file, '%s,', num2str(power_test_settings(i)));   % line (2)
end
fprintf(settings_file, '\n');
% GARP - perform GARP WARP and SARP, with/without residuals
for i=1:max(size(GARP_flags))
    fprintf(settings_file, '%s,', num2str(GARP_flags(i)));   % line (3)
end
fprintf(settings_file, '\n');
% AFRIAT - perform AFRIAT, with/without residuals
for i=1:max(size(AFRIAT_flags))
    fprintf(settings_file, '%s,', num2str(AFRIAT_flags(i)));   % line (4)
end
fprintf(settings_file, '\n');
% VARIAN - perform VARIAN, with/without residuals
for i=1:max(size(VARIAN_flags))
    fprintf(settings_file, '%s,', num2str(VARIAN_flags(i)));   % line (5)
end
fprintf(settings_file, '\n');
% HOUTMAN-MAKS - perform HOUTMAN-MAKS, with/without residuals
for i=1:max(size(HOUTMAN_flags))
    fprintf(settings_file, '%s,', num2str(HOUTMAN_flags(i)));   % line (6)
end
fprintf(settings_file, '\n');
% MPI - perform MPI, with/without residuals
for i=1:max(size(MPI_flags))
    fprintf(settings_file, '%s,', num2str(MPI_flags(i)));   % line (7)
end


end