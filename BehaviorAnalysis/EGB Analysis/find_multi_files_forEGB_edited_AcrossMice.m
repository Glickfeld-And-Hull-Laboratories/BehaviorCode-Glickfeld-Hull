% example of indexing name in cell array: mouse_names{1}{1}(1,:)
cd('Z:\Behavior\Data')
mouse_number = 2;
mouse_names{:,:} = {'1401'; '1402'};
mouse_start_dates{:,:} = {'190509';'190613'};
mouse_end_dates{:,:} = {'190617';'190628'};

for mouse = 1:mouse_number

mouse_name = ['*', mouse_names{1}{mouse}(1:end), '*'];%should be mouse name, must be starred for dir function to work 
first_pass_filenames = ['i',mouse_names{1}{mouse}(1:end),'-1']; %should be mouse name and year
start_date = mouse_start_dates{1}{mouse}(1:end);
end_date = mouse_end_dates{1}{mouse}(1:end);

%below- finds portions of strings like dates/names, finds file names containing these strings 
dir_all = dir(mouse_name); 
mouse_idx = char(dir_all.name);
mouse_idx = cellstr(mouse_idx);
date_idx = strfind(mouse_idx, ['i', mouse_names{1}{mouse}(1:end), '-1']); %Input used to be 'i1400-18' for the year, but can just get rid of last # to include multiple years
target_mouse_idx = find(~cellfun(@isempty,date_idx));
mouse_dates = mouse_idx(target_mouse_idx); %mouse_dates is cell array containing file names with dates of interest 
start_date_idx = strfind(mouse_dates,start_date);
start_date = find(~cellfun(@isempty,start_date_idx));
end_date_idx = strfind(mouse_dates, end_date);
end_date = find(~cellfun(@isempty,end_date_idx));
start_date = start_date(end);
end_date = end_date(end);

%celleqel2mat_padded(days_to_analyze);

for date = start_date:end_date
load(mouse_dates{date});
days_to_analyze{mouse}{date-(start_date-1)} = input; %days_to_analyze is a cell array full of structures
end

%days_to_analyze = days_to_analyze';
%input = days_to_analyze; %just a name change
%clearvars -except input
end

input = horzcat(days_to_analyze{1}, days_to_analyze{2})';
clearvars -except input
    

    
    
    
    