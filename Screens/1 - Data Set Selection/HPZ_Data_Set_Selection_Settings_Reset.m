function HPZ_Data_Set_Selection_Settings_Reset(main_folder)

% this file resets the dataset settings.
% it is automatically called when reading the dataset settings fails.



% make sure the Settings Files directory exists, if not - create it
dir_exists = exist(strcat(main_folder, '/', HPZ_Constants.settings_files_dir) , 'dir');
if ~dir_exists
    mkdir(strcat(main_folder, '/', HPZ_Constants.settings_files_dir));
end



% initialize the data list to be empty
data_list_str = {'Choi et al. (2007)','Halevy et al. (2016)'};
data_list_path = {strcat(HPZ_Constants.data_files_dir, '/Data_CFGK_2007.csv'),...
                  strcat(HPZ_Constants.data_files_dir, '/Data_HPZ_2016.csv')};
% data_list_str = {'Choi et al. (2007)','Kurtz et al. (2016)','Halevy et al. (2016)'};
% data_list_path = {strcat(main_folder, '/', HPZ_Constants.data_files_dir, '/Data_CFGK_2007.csv'),...
%                 strcat(main_folder, '/', HPZ_Constants.data_files_dir, '/Data_KLP_2016.csv'),...
%                 strcat(main_folder, '/', HPZ_Constants.data_files_dir, '/Data_HPZ_2016.csv')};
data_list_prefs = {'1','1'};
data_list_subject = {'1','1'};
data_list_obs = {'2','2'};
data_list_quantity1 = {'3','3'};
data_list_quantity2 = {'4','4'};
data_list_maxquantity1 = {'5','5'};
data_list_maxquantity2 = {'6','6'};
data_set = {'1',''};
fix_endowments = {'1',''};

% create the data settings table
data_settings = table(data_list_str', data_list_path', data_list_prefs', data_list_subject', data_list_obs', data_list_quantity1', data_list_quantity2', data_list_maxquantity1', data_list_maxquantity2', data_set', fix_endowments');

% set the headers of the table
for i=1:max(size(HPZ_Constants.data_settings_headers))
    data_settings.Properties.VariableNames{i} = HPZ_Constants.data_settings_headers{i};
end

% print the table to the data settings file
print_table_to_file(strcat(main_folder, '/', HPZ_Constants.settings_files_dir, '/', HPZ_Constants.data_settings_file_name), data_settings);

end