function [Graph_flags, GARP_flags, AFRIAT_flags, VARIAN_flags, HOUTMAN_flags, MPI_flags] = HPZ_Consistency_Indices_Settings_Read(main_folder)

% this function reads the user-interface screens saved settings

% for detailed explanations about input/output variables that possess
% the same name and meaning in multiple functions (e.g. data, action_flag, 
% pref_class, etc.) see: HPZ_Variables_Documentation in the "Others" sub-folder 



% read settings from file 
try
    % reading from the file
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.consistency_indices_settings_file_name, '.csv'));
    
    % in the next line, we call the6th row and 7th column in the settings
    % matrix, because we have 6 rows and the max number of columns is 7
    % (it is 7 for GARP flags).
    % if it doesn't exists, it means that the settings file is corrupted,
    % cause it doesn't have all required data.
    % if in the future you will add more rows (>6) or use longer rows (num
    % of columns > 7), you should increase these numbers respectively.
    check_settings = settings(6 , 7); %#ok<NASGU>
catch
    % A) if we failed to read the file,
    % we assume it might be because it was accidently deleted, or because
    % it is not formatted in a way that allows it to be read.
    % B) if the file's data was corrupted, we assume someone might've
    % played with it and made changes to it by mistake.
    % A+B) either way, we reset the file, then read it again.
    HPZ_Consistency_Indices_Settings_Reset(main_folder);
    [settings] = csvread(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.consistency_indices_settings_file_name, '.csv'));
end

% assign initial values to variables in accordance with the settings
Graph_flags = settings(1, 1:4);
GARP_flags = settings(2, 1:7);
AFRIAT_flags = settings(3, 1:4);
VARIAN_flags = settings(4, 1:4);
HOUTMAN_flags = settings(5, 1:4);
MPI_flags = settings(6, 1:4);



%% for (almost) every variable, we check if it has a valid value, otherwise we reset it
% consistency flags
for i=1:length(Graph_flags)
    if (Graph_flags(i) ~= 0 && Graph_flags(i) ~= 1)
        Graph_flags(i) = 1;
    end
end
for i=1:length(GARP_flags)
    if (GARP_flags(i) ~= 0 && GARP_flags(i) ~= 1)
        GARP_flags(i) = 1;
    end
end
for i=1:length(AFRIAT_flags)
    if (AFRIAT_flags(i) ~= 0 && AFRIAT_flags(i) ~= 1)
        AFRIAT_flags(i) = 1;
    end
end
for i=1:length(VARIAN_flags)
    if (VARIAN_flags(i) ~= 0 && VARIAN_flags(i) ~= 1)
        VARIAN_flags(i) = 1;
    end
end
for i=1:length(HOUTMAN_flags)
    if (HOUTMAN_flags(i) ~= 0 && HOUTMAN_flags(i) ~= 1)
        HOUTMAN_flags(i) = 1;
    end
end


end