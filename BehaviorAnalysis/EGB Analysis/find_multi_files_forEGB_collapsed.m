cd('Z:\Behavior\Data')

mouse_name = '*1402*';%should be mouse name, must be starred for dir function to work 
first_pass_filenames = 'i1402-19'; %should be mouse name and year
start_date = '190614';
end_date = '190711';


%below- finds portions of strings like dates/names, finds file names containing these strings 
dir_all = dir(mouse_name);
mouse_idx = char(dir_all.name);
mouse_idx = cellstr(mouse_idx);
date_idx = strfind(mouse_idx, 'i1402-1'); %Input used to be 'i1400-18' for the year, but can just get rid of last # to include multiple years
target_mouse_idx = find(~cellfun(@isempty,date_idx));
mouse_dates = mouse_idx(target_mouse_idx); %mouse_dates is cell array containing file names with dates of interest 
start_date_idx = strfind(mouse_dates,start_date);
start_date = find(~cellfun(@isempty,start_date_idx));
end_date_idx = strfind(mouse_dates, end_date);
end_date = find(~cellfun(@isempty,end_date_idx));
start_date = start_date(end);
end_date = end_date(end);

onlyIn1_all =[];
onlyIn2_all =[];
diff_dates = [];
start = 1;
clear days_to_analyze
for date = start_date:end_date
    
load(mouse_dates{date});
input2 = input;
if date>start_date
    input1 = days_to_analyze(date-(start_date-1)-1);
    
    [isDiff onlyIn1 onlyIn2] = compareStructures(input1,input2);
    if isDiff
        diff_dates = [diff_dates date];
        onlyIn1_all{start} = onlyIn1;
        onlyIn2_all{start} = onlyIn2;
        start = start+1;
    end
        
    if isDiff & length(onlyIn1)>0
        for i = 1:length(onlyIn1)
            aField = char(onlyIn1(i));
            if length(input1.(aField)) == 1
                input2.(aField) = NaN;
            elseif length(input1.(aField)) == length(input1.trialOutcomeCell)
                input2.(aField) = cell(size(input2.trialOutcomeCell));
            end
        end
    elseif isDiff & length(onlyIn2)>0
        for id = start_date:date-1
            for i = 1:length(onlyIn2)
                aField = char(onlyIn2(i));
                if length(input2.(aField)) == 1
                    input1.(aField) = NaN;
                elseif length(input2.(aField)) == length(input2.trialOutcomeCell)
                    input1.(aField) = cell(size(input1.trialOutcomeCell));
                end
            end
            days_to_analyze(id-(start_date-1)) = input1;
        end
    end
end
days_to_analyze(date-(start_date-1)) = input2;
end
all_days = concatenateStructuresLG(days_to_analyze);            

input = all_days;
clearvars -except input


