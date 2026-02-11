clear all
close all
clc

%%
%mouse = "i2174";
%dart = "pre"
%filenames = ["data-i2174-250225-1539.mat", "data-i2174-250226-1037.mat"]
%dart = "post"
%filenames = ["data-i2174-250228-1047.mat", "data-i2174-250301-1312.mat"]

%mouse = "i2182";
%dart = "pre"
%filenames = ["data-i2182-250304-1456.mat", "data-i2182-250305-1055.mat"]
%dart = "post"
%filenames = ["data-i2182-250306-1140.mat", "data-i2182-250307-0944"]

%mouse = "i2179";
%dart = "pre"
%filenames = ["data-i2179-250304-1459.mat", "data-i2179-250305-1059.mat"]
%dart = "post"
%filenames = ["data-i2179-250306-1144.mat", "data-i2179-250307-0948.mat"]';

%mouse = "i2140"
%dart = "pre"
%filenames = ["data-i2140-240904-1152.mat", "data-i2140-240905-1151.mat"]
%dart = "post"
%filenames = ["data-i2140-240906-1055.mat", "data-i2140-240907-1130.mat"]

mouse = "i2167"
%dart = "pre"
%filenames = ["data-i2167-241211-1203.mat", "data-i2167-241213-1051.mat"] 
dart = "post"
filenames = ["data-i2167-241214-1130.mat", "data-i2167-241215-1035.mat"]

%%

dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
n_days = length(filenames);

for i = 1:length(filenames)
    load(fullfile(dir,filenames(i)));
    temp(i) = input;
end

clear input
input = concatenateDataBlocks(temp)

%%
 savedir = "Z:\All_Staff\Behavior\Data";
 if dart == "pre"
     datestr = "000000"
     timestr = "0000";
 elseif dart == "post"
     datestr = "999999"
     timestr = "9999";
 end

 filename = join(["data", mouse, datestr, timestr], "-")
 fullname = join([savedir, filename], "\")

 save(fullname,"input")