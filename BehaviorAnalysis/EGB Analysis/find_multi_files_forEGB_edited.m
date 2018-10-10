

cd('\\crash.dhe.duke.edu\data\home\andrew\Behavior\Data')

mouse_name = '*602*';%should be mouse name, must be starred for dir function to work 
first_pass_filenames = 'i602-18'; %should be mouse name and year
start_date = '180920';
end_date = '180928';

dir_all = dir(mouse_name);
mouse_idx = char(dir_all.name);
mouse_idx = cellstr(mouse_idx);
date_idx = strfind(mouse_idx, 'i602-18');
target_mouse_idx = find(~cellfun(@isempty,date_idx));
mouse_dates = mouse_idx(target_mouse_idx);
start_date_idx = strfind(mouse_dates,start_date);
start_date = find(~cellfun(@isempty,start_date_idx));
end_date_idx = strfind(mouse_dates, end_date);
end_date = find(~cellfun(@isempty,end_date_idx));
start_date = start_date(end);
end_date = end_date(end);

%celleqel2mat_padded(days_to_analyze);
%%
for date = start_date:end_date
load(mouse_dates{date});
days_to_analyze{date-(start_date-1)} = input;
end

days_to_analyze = days_to_analyze';
input = days_to_analyze;
clearvars -except input



    

    
    
    
    