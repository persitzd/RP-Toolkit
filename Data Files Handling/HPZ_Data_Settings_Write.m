function HPZ_Data_Settings_Write(data_list_str, data_list_path, data_list_prefs, data_list_subject, data_list_obs, data_list_quantity1, data_list_quantity2, data_list_maxquantity1, data_list_maxquantity2, data_set, main_folder)

% this file write the dataset settings to the file, after the user may have 
% changed them (in the screen presented by HPZ_Data_Set_Selection).
% for each dataset there is its name, path, columns numbers of the required 
% data, and preferences class. 
% also it saves the index of the last dataset that was chosen.



% make sure the Settings Files directory exists, if not - create it
dir_exists = exist(strcat(main_folder, '/', HPZ_Constants.settings_files_dir) , 'dir');
if ~dir_exists
    mkdir(strcat(main_folder, '/', HPZ_Constants.settings_files_dir));
end



% we need to make data_set also a cell of the same length as the others
data_set = {data_set};
if (length(data_list_str) > 1)
    data_set{length(data_list_str)} = '';
end

% for each of the array cells - if it is a row cell, we turn it into a column cell 
sz = size(data_list_str);
if (sz(1) < sz(2))
    data_list_str = data_list_str';
end
sz = size(data_list_path);
if (sz(1) < sz(2))
    data_list_path = data_list_path';
end
sz = size(data_list_prefs);
if (sz(1) < sz(2))
    data_list_prefs = data_list_prefs';
end
sz = size(data_list_subject);
if (sz(1) < sz(2))
    data_list_subject = data_list_subject';
end
sz = size(data_list_obs);
if (sz(1) < sz(2))
    data_list_obs = data_list_obs';
end
sz = size(data_list_quantity1);
if (sz(1) < sz(2))
    data_list_quantity1 = data_list_quantity1';
end
sz = size(data_list_quantity2);
if (sz(1) < sz(2))
    data_list_quantity2 = data_list_quantity2';
end
sz = size(data_list_maxquantity1);
if (sz(1) < sz(2))
    data_list_maxquantity1 = data_list_maxquantity1';
end
sz = size(data_list_maxquantity2);
if (sz(1) < sz(2))
    data_list_maxquantity2 = data_list_maxquantity2';
end
sz = size(data_set);
if (sz(1) < sz(2))
    data_set = data_set';
end

% create the data settings table
data_settings = table(data_list_str, data_list_path, data_list_prefs, data_list_subject, data_list_obs, data_list_quantity1, data_list_quantity2, data_list_maxquantity1, data_list_maxquantity2, data_set);

% set the headers of the table
for i=1:max(size(HPZ_Constants.data_settings_headers))
    data_settings.Properties.VariableNames{i} = HPZ_Constants.data_settings_headers{i};
end

% print the table to the data settings file
print_table_to_file(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.data_settings_file_name), data_settings);

end