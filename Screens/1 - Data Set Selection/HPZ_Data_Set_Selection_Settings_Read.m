function [data_list_str, data_list_path, data_list_prefs, data_list_choice_set_types, data_list_subject, data_list_obs, data_list_quantity1, data_list_quantity2, data_list_maxquantity1, data_list_maxquantity2, data_set, fix_endowments] = HPZ_Data_Set_Selection_Settings_Read(main_folder)

% this file reads the dataset settings.
% for each dataset there is its name, path, columns numbers of the required 
% data, and preferences class. 
% also it saves the index of the last dataset that was chosen.



try
    % read the data settings from file
    data_settings = readtable(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.data_settings_file_name, '.csv'));
    
    % extracting the data as seperate cell arrays
    data_list_str = data_settings{:,HPZ_Constants.data_name};
    data_list_path = data_settings{:,HPZ_Constants.file_name};
    data_list_prefs = data_settings{:,HPZ_Constants.pref_class};
    data_list_choice_set_types = data_settings{:,HPZ_Constants.choice_set_type}; 
    data_list_subject = data_settings{:,HPZ_Constants.subject_index};
    data_list_obs = data_settings{:,HPZ_Constants.obs_index};
    data_list_quantity1 = data_settings{:,HPZ_Constants.quantity1_index};
    data_list_quantity2 = data_settings{:,HPZ_Constants.quantity2_index};
    data_list_maxquantity1 = data_settings{:,HPZ_Constants.maxquantity1_index};
    data_list_maxquantity2 = data_settings{:,HPZ_Constants.maxquantity2_index};
    % this is just for checking
    data_set = data_settings{:,HPZ_Constants.data_set}; %#ok<NASGU>
    fix_endowments = data_settings{:,HPZ_Constants.fix_endowments}; %#ok<NASGU>
catch
    % if failed - we reset the file then read it again
    HPZ_Data_Set_Selection_Settings_Reset(main_folder);
    
    % read the data settings from file (again)
    data_settings = readtable(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.data_settings_file_name, '.csv'));

    % extracting the data as separate cell arrays (again)
    data_list_str = data_settings{:,HPZ_Constants.data_name};
    data_list_path = data_settings{:,HPZ_Constants.file_name};
    data_list_prefs = data_settings{:,HPZ_Constants.pref_class};
    data_list_choice_set_types = data_settings{:,HPZ_Constants.choice_set_type}; 
    data_list_subject = data_settings{:,HPZ_Constants.subject_index};
    data_list_obs = data_settings{:,HPZ_Constants.obs_index};
    data_list_quantity1 = data_settings{:,HPZ_Constants.quantity1_index};
    data_list_quantity2 = data_settings{:,HPZ_Constants.quantity2_index};
    data_list_maxquantity1 = data_settings{:,HPZ_Constants.maxquantity1_index};
    data_list_maxquantity2 = data_settings{:,HPZ_Constants.maxquantity2_index};
end



if isempty(data_list_str)
    % we need this to handle the case that the file is empty
    data_set = data_settings{:,HPZ_Constants.data_set};
    fix_endowments = 1;
else
    data_set = data_settings{1,HPZ_Constants.data_set};
    fix_endowments = data_settings{1,HPZ_Constants.fix_endowments};
    
    length = max(size(data_list_str));
    
    % just in case, to prevent bugs
    if data_set > length
        data_set = length;
    elseif data_set < 1 || mod(data_set,1) ~= 0
        data_set = 1;
    end
    if (fix_endowments ~= 0 && fix_endowments ~= 1)
        fix_endowments = 1;
    end
    
    % we need to change 7 of the above from vectors to cell arrays
    
    data_list_prefs_temp = data_list_prefs;
    data_list_choice_set_types_temp = data_list_choice_set_types;
    data_list_subject_temp = data_list_subject;
    data_list_obs_temp = data_list_obs;
    data_list_quantity1_temp = data_list_quantity1;
    data_list_quantity2_temp = data_list_quantity2;
    data_list_maxquantity1_temp = data_list_maxquantity1;
    data_list_maxquantity2_temp = data_list_maxquantity2;
    
    data_list_prefs = cell(1,length);
    data_list_choice_set_types = cell(1,length);
    data_list_subject = cell(1,length);
    data_list_obs = cell(1,length);
    data_list_quantity1 = cell(1,length);
    data_list_quantity2 = cell(1,length);
    data_list_maxquantity1 = cell(1,length);
    data_list_maxquantity2 = cell(1,length);
    
    for i=1:length
        data_list_prefs{i} = data_list_prefs_temp(i);
        data_list_choice_set_types{i} = data_list_choice_set_types_temp(i);
        data_list_subject{i} = data_list_subject_temp(i);
        data_list_obs{i} = data_list_obs_temp(i);
        data_list_quantity1{i} = data_list_quantity1_temp(i);
        data_list_quantity2{i} = data_list_quantity2_temp(i);
        data_list_maxquantity1{i} = data_list_maxquantity1_temp(i);
        data_list_maxquantity2{i} = data_list_maxquantity2_temp(i);
    end
end

end