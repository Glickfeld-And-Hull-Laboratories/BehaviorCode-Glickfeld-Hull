cd('Z:\Behavior\Data')

mouse_name = '*1402*';%should be mouse name, must be starred for dir function to work 
first_pass_filenames = 'i1402-19'; %should be mouse name and year
start_date = '190617';
end_date = '190828';


%below- finds portions of strings like dates/names, finds file names containing these strings 
dir_all = dir(mouse_name); %dir fxn loads all files from directory with given mouse name
mouse_idx = char(dir_all.name); %grabs just the name of the files from the directory
mouse_idx = cellstr(mouse_idx); %converts from character array to cell array
date_idx = strfind(mouse_idx, 'i1402-1'); %goes into array of file names & finds only 1402 exactly
%(Input here used to be 'i1400-18';got rid of last # for multiple years)
target_mouse_idx = find(~cellfun(@isempty,date_idx)); %list of positions with only 1402 exactly by finding target_mouse_idx cells that are not empty

mouse_dates = mouse_idx(target_mouse_idx); %goes back to directory and selects only relevant positions-list of all file names 
start_date_idx = strfind(mouse_dates,start_date); %find cell pos of start date in mouse_dates 
start_date = find(~cellfun(@isempty,start_date_idx)); %lists position of start date by finding cell that is not empty
end_date_idx = strfind(mouse_dates, end_date); %repeat for end date
end_date = find(~cellfun(@isempty,end_date_idx));

start_date = start_date(end); %??? does not appear to do anything
end_date = end_date(end);

%% Next part will deal with variables that differ on diff days
% Will also call upon the files from the directory, from the start pos to the end pos

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

% clearvars -except days_to_analyze

%Output days_to_analyze is a structure containing arrays of each day's data
%% Here is where we need to crop days
%Input lower and upper bounds for each day as b and c

% c = [384 150 0 175 409 390 466 476 400 465 400 384 275 0 300 446 250 200];
% b = [1 1 1 1 1 1 1 1 1 1 200 1 1 1 200 175 1 1];

%Flash Adapt 1402 only 6/17 - 8/28. Delete #2 - 6/18 6 trials. #34 - 8/05
%378 trials
x = [617 618 619 620 621 624 625 626 627 628 ...
    702 703 708 709 710 711 712 715 716 717 718 719 ...
    722 723 724 725 726 729 730 731 801 802 ...
    806 807 808 809 812 813 814 815 816 817 819 ...
    820 821 822 823 826 827 828];

c = [150 175 409 390 466 476 400 465 400 384 ...  
    275 0 300 446 250 200 0 0 0 0 0 0 ...
    470 509 403 441 591 332 352 300 457 431 ...
    377 200 300 383 0 0 0 0 0 0 0 ...
    324 412 0 350 378 367 159];

b = [21 21 21 21 21 21 21 21 200 21 ...
    21 21 200 175 21 21 1 1 1 1 1 1 ...
    21 21 21 21 21 21 21 21 21 21 ...
    21 21 21 21 1 1 1 1 1 1 1 ...
    21 21 1 21 21 21 21];

for i = 1:length(days_to_analyze) 
    days_to_analyze(i).trialOutcomeCell = days_to_analyze(i).trialOutcomeCell(b(i):c(i));
    days_to_analyze(i).tLeftTrial = days_to_analyze(i).tLeftTrial(b(i):c(i));
    days_to_analyze(i).tDecisionTimeMs = days_to_analyze(i).tDecisionTimeMs(b(i):c(i));
    days_to_analyze(i).tConsecCorrects = days_to_analyze(i).tConsecCorrects(b(i):c(i));
    days_to_analyze(i).tConsecErrors = days_to_analyze(i).tConsecErrors(b(i):c(i));
    days_to_analyze(i).tNTrialsCompleted = days_to_analyze(i).tNTrialsCompleted(b(i):c(i));
    days_to_analyze(i).stimTimestampMs = days_to_analyze(i).stimTimestampMs(b(i):c(i));
    days_to_analyze(i).qStimOn = days_to_analyze(i).qStimOn(b(i):c(i));
    days_to_analyze(i).qTrialStart = days_to_analyze(i).qTrialStart(b(i):c(i));
    days_to_analyze(i).qStartReact = days_to_analyze(i).qStartReact(b(i):c(i));
    days_to_analyze(i).tLeftResponse = days_to_analyze(i).tLeftResponse(b(i):c(i));
    days_to_analyze(i).tRightResponse = days_to_analyze(i).tRightResponse(b(i):c(i));
    days_to_analyze(i).tGratingDirectionStart = days_to_analyze(i).tGratingDirectionStart(b(i):c(i));
    days_to_analyze(i).adapterTimestampMs = days_to_analyze(i).adapterTimestampMs(b(i):c(i));
    days_to_analyze(i).aGratingContrast = days_to_analyze(i).aGratingContrast(b(i):c(i));
    days_to_analyze(i).aGratingDirectionDeg = days_to_analyze(i).aGratingDirectionDeg(b(i):c(i));
    days_to_analyze(i).quadratureValues = days_to_analyze(i).quadratureValues(b(i):c(i))
end

%% This part uses Lindsey's function to collapse into a single data structure 


all_days = concatenateStructuresLG(days_to_analyze);            

input = all_days;
clearvars -except input


