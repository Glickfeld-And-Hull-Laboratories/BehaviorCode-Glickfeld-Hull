clear all
close all

%% Get the file names we want
behav_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\';

data_path = 'Z:\home\LindseyW\MATLAB\Behavior_analysis';
dataset = 'CurrentMouseList';
eval(dataset)


mouselist = string({allmice(:).name});

behav_files_struct = dir(behav_path);
behav_file_names = string({behav_files_struct(:).name});
mouseIx = zeros(1,length(behav_file_names));
mouseIxMat = zeros(length(mouselist),length(behav_file_names));

for i = 1:length(mouselist);
    for file = 1:length(behav_file_names);
        if contains(behav_file_names(file),mouselist(i));
            mouseIx(file) = 1;
        else
            mouseIx(file) = 0;
        end
    mouseIxMat(i,:) = mouseIx;
    end
end

validIx = sum(mouseIxMat, 1);
filenames = behav_file_names(validIx == 1)'

%% Generate/Save BehaviorIndices table

metadata = strings(length(filenames),1);
for file = 1:length(filenames);
    metadata(file) = extractBetween(filenames(file), "data-", ".mat");
end

metadata = split(metadata, "-");
mouseNames = metadata(:,1);
dates = metadata(:,2);
times = metadata(:,3);

BehaviorIndices = table(filenames,mouseNames,dates,times)
writetable(BehaviorIndices, 'BehaviorIndices.xlsx', 'WriteRowNames',true)



clear behav_file_names behav_files_struct mouseIx mouseIxMat


%Figure out how to add good trials programmatically